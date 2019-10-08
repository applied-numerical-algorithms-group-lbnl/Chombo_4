#include <cmath>
#include <cstdio>
#include <iostream>


#include "EBProto.H"
#include "EBLevelBoxData.H"
#include "LevelData.H"
#include "BaseFab.H"

#include "ParmParse.H"
#include "LoadBalance.H"
#include "ProtoInterface.H"
#include "BRMeshRefine.H"
#include "GeometryService.H"
#include "EBDictionary.H"
#include "EBMultigrid.H"
#include "EBChombo.H"
#include <iomanip>


using std::cout;
using std::endl;
using std::shared_ptr;

typedef Proto::Box Bx;
using   Proto::Point;
using   Proto::BoxData;
using   Proto::Stencil;
using   ProtoCh::getPoint;
using   ProtoCh::getProtoBox;
using   ProtoCh::getIntVect;
using   ProtoCh::getBox;
using     std::cout;
using     std::endl;
using     std::shared_ptr;
using   Proto::BaseIF;
using   Proto::SimpleEllipsoidIF;
using   Proto::CENTERING;
using   Proto::CELL;
using Proto::PointSet;
using Proto::PointSetIterator;
void 
dumpOrigin(EBBoxData<CELL, Real, 1>* dataPtr)
{
  if(dataPtr != NULL)
  {
    cout    << setprecision(6)
            << setiosflags(ios::showpoint)
            << setiosflags(ios::scientific);
    typedef EBIndex<CELL> VolIndex;
    IntVect iv(1, 1);
    BoxData<Real, 1> & data = dataPtr->getRegData();
    Box region = grow(Box(iv, iv), 2);
    IntVect lo = region.smallEnd();
    IntVect hi = region.bigEnd();
    cout << "data region contains:" << endl;
    for(int j = lo[1]; j <= hi[1]; j++)
    {
      for(int i = lo[0]; i <= hi[0]; i++)
      {
        Point pt(i,j);
        cout << pt << ":" << data(pt, 0) << "  ";
      }
      cout << endl;
    }
  }
}


void 
dumpArea(EBBoxData<CELL, Real, 1>* dataPtr)
{
  if(dataPtr != NULL)
  {
    cout    << setprecision(6)
            << setiosflags(ios::showpoint)
            << setiosflags(ios::scientific);
    typedef EBIndex<CELL> VolIndex;
    IntVect iv(26, 17);
    EBBoxData<CELL, Real, 1> & data = *dataPtr;
    EBGraph  ebgraph = data.ebgraph();
    Box region = grow(Box(iv, iv), 2);
    IntVect lo = region.smallEnd();
    IntVect hi = region.bigEnd();
    cout << "data region contains:" << endl;
    for(int j = lo[1]; j <= hi[1]; j++)
    {
      for(int i = lo[0]; i <= hi[0]; i++)
      {
        Point pt(i,j);
        vector<VolIndex> vofs = ebgraph.getVoFs(pt);
        for(int ivof = 0; ivof < vofs.size(); ivof++)
        {
          VolIndex vof = vofs[ivof];
          cout << pt << ":" << data(vof, 0) << "  ";
        }
      }
      cout << endl;
    }
  }
}
void
dumpPPS(const PointSet* a_ivs)
{
  for(PointSetIterator ivsit(*a_ivs);  ivsit.ok(); ++ivsit)
  {
    std::cout << ivsit() << " " ;
  }
  std::cout << std::endl;
 }

typedef Proto::Var<Real, 1> Sca;
PROTO_KERNEL_START 
void  addCorToPhiF(Sca     a_phi,
                   Sca     a_cor)
{
  a_phi(0) = a_phi(0) + a_cor(0);
}
PROTO_KERNEL_END(addCorToPhiF, addCorToPhi)

int
runTest(int a_argc, char* a_argv[])
{
  dumpArea(NULL);

  Real coveredval = -1;
  int nx      = 32;
  int maxGrid = 32;
  Real x0 = 0.5;
  Real y0 = 0.5;
  Real z0 = 0.5;
  Real A = 1.0;
  Real B = 1.0;
  Real C = 1.0;
  Real R = 0.25;
  Real tol = 0.00001;
  Real alpha = 1.0;
  Real beta = -0.001;
    
  int  maxIter = 10;
  int nStream    = 8;
  int numSmooth  = 2;
  int dombc = 1;
  int ebbc  = 1;
  ParmParse pp;
    
  pp.get("domainBC"  , dombc);
  pp.get("EBBC"      , ebbc);
  pp.get("nx"        , nx);
  pp.get("maxIter"   , maxIter);
  pp.get("nStream"   , nStream);
  pp.get("maxGrid"   , maxGrid);
  pp.get("alpha"     , alpha);
  pp.get("beta"      , beta);
  pp.get("tolerance" , tol);
  pp.get("x0"        , x0);
  pp.get("y0"        , y0);
  pp.get("z0"        , z0);
  pp.get("A"         , A);
  pp.get("B"         , B);
  pp.get("C"         , C);
  pp.get("R"         , R);         
  pp.get("coveredval", coveredval);         
  pp.get("numSmooth" , numSmooth);         

  pout() << "nx        = " << nx       << endl;
  pout() << "maxGrid   = " << maxGrid  << endl;
  pout() << "x0        = " << x0       << endl;
  pout() << "y0        = " << y0       << endl;
  pout() << "z0        = " << z0       << endl;
  pout() << "A         = " << A        << endl;
  pout() << "B         = " << B        << endl;
  pout() << "C         = " << C        << endl;
  pout() << "R         = " << R        << endl;
  pout() << "alpha     = " << alpha    << endl;
  pout() << "beta      = " << beta    << endl;
  pout() << "tolerance = " << tol    << endl;
  pout() << "numSmooth = " << numSmooth    << endl;

  pout() << "maxIter   = " << maxIter    << endl;
  pout() << "nstream = " << nStream  << endl;

#ifdef PROTO_CUDA
  Proto::DisjointBoxLayout::setNumStreams(nStream);
#endif

  RealVect ABC, X0;
  ABC[0] = A;
  ABC[1] = B;
  X0[0] = x0;
  X0[1] = y0;
#if DIM==3
  ABC[2] = C;
  X0[2] = z0;
#endif
  IntVect domLo = IntVect::Zero;
  IntVect domHi  = (nx - 1)*IntVect::Unit;

// EB and periodic do not mix
  ProblemDomain domain(domLo, domHi);

  Vector<DisjointBoxLayout> vecgrids;
  pout() << "making grids" << endl;
  GeometryService<2>::generateGrids(vecgrids, domain.domainBox(), maxGrid);

  DisjointBoxLayout grids = vecgrids[0];
  grids.printBalance();

  IntVect dataGhostIV =   IntVect::Unit;
  Point   dataGhostPt = getPoint(dataGhostIV); 
  int geomGhost = 4;
  RealVect origin = RealVect::Zero();
  Real dx = 1.0/nx;
//  Real dx = 1.0;
  shared_ptr<BaseIF>                       impfunc(new SimpleEllipsoidIF(ABC, X0, R, false));
  Bx domainpr = getProtoBox(domain.domainBox());

  pout() << "defining geometry" << endl;
  GeometryService<2>* geomptr = new GeometryService<2>(impfunc, origin, dx, domain.domainBox(), vecgrids, geomGhost);
//  GeometryService<2>* geomptr = new GeometryService<2>(impfunc, origin, dx, domain.domainBox(), vecgrids[0], geomGhost);
  shared_ptr< GeometryService<2> >  geoserv(geomptr);

  pout() << "making dictionary" << endl;
  shared_ptr<EBDictionary<2, Real, CELL, CELL> > 
    dictionary(new EBDictionary<2, Real, CELL, CELL>(geoserv, grids, domain.domainBox(), dataGhostPt, dataGhostPt, dx));
  
//  typedef EBStencil<2, Real, CELL, CELL> ebstencil_t;
  string stenname("Second_Order_Poisson");
  string dombcname, ebbcname;
  if(dombc == 0)
  {
    dombcname = string("Neumann");
    pout() << "using Neumann BCs at domain" << endl;
  }
  else
  {
    dombcname = string("Dirichlet");
    pout() << "using Dirichlet BCs at domain" << endl;
  }

  if(ebbc == 0)
  {
    ebbcname = string("Neumann");
    pout() << "using Neumann BCs at EB" << endl;
  }
  else
  {
    ebbcname = string("Dirichlet");
    pout() << "using Dirichlet BCs at EB" << endl;
  }

  shared_ptr<LevelData<EBGraph> > graphs = geoserv->getGraphs(domain.domainBox());

  pout() << "making data" << endl;
  EBLevelBoxData<CELL,   1>  phi(grids, dataGhostIV, graphs);
  EBLevelBoxData<CELL,   1>  rhs(grids, dataGhostIV, graphs);
  EBLevelBoxData<CELL,   1>  res(grids, dataGhostIV, graphs);
  EBLevelBoxData<CELL,   1>  cor(grids, dataGhostIV, graphs);

  EBMultigrid solver(dictionary, geoserv, alpha, beta, dx, grids, stenname, dombcname, ebbcname, domain.domainBox(), dataGhostIV, dataGhostIV);
  EBMultigrid::s_numSmoothUp   = numSmooth;
  EBMultigrid::s_numSmoothDown = numSmooth;
  DataIterator dit = grids.dataIterator();
  pout() << "setting values" << endl;
  for(int ibox = 0; ibox < dit.size(); ibox++)
  {
    EBBoxData<CELL, Real, 1>& phibd = phi[dit[ibox]];
    EBBoxData<CELL, Real, 1>& rhsbd = rhs[dit[ibox]];
    EBBoxData<CELL, Real, 1>& corbd = cor[dit[ibox]];
    phibd.setVal(1.0);
    rhsbd.setVal(0.0);
    corbd.setVal(0.0);
  }

  solver.solve(phi, rhs, tol, maxIter);

  pout() << "writing to file " << endl;
  
  phi.writeToFileHDF5("phi.hdf5",  -1.0);
  rhs.writeToFileHDF5("rhs.hdf5",   0.0);
  res.writeToFileHDF5("res.hdf5",   0.0);
  
  pout() << "exiting " << endl;
  return 0;
}


int main(int a_argc, char* a_argv[])
{
#ifdef CH_MPI
  MPI_Init(&a_argc, &a_argv);
  pout() << "MPI INIT called" << std::endl;
#endif
  //needs to be called after MPI_Init
  CH_TIMER_SETFILE("ebapply.time.table");
  {
    if (a_argc < 2)
    {
      cerr<< " usage " << a_argv[0] << " <input_file_name> " << endl;
      exit(0);
    }
    char* in_file = a_argv[1];
    ParmParse  pp(a_argc-2,a_argv+2,NULL,in_file);
    runTest(a_argc, a_argv);
  }

  pout() << "printing time table " << endl;
  CH_TIMER_REPORT();
#ifdef CH_MPI
  pout() << "about to call MPI Finalize" << std::endl;
  MPI_Finalize();
#endif
  return 0;
}
