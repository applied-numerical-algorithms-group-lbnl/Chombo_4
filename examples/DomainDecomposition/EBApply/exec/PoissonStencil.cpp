#include <cmath>
#include <cstdio>
#include <iostream>


#include "EBProto.H"
#include "Chombo_EBLevelBoxData.H"
#include "Chombo_LevelData.H"
#include "Chombo_BaseFab.H"

#include "Chombo_ParmParse.H"
#include "Chombo_LoadBalance.H"
#include "Chombo_ProtoInterface.H"
#include "Chombo_BRMeshRefine.H"
#include "Chombo_GeometryService.H"
#include "Chombo_GeometryService.H"
#include "Chombo_EBDictionary.H"
#include "Chombo_EBChombo.H"
#include "Proto_SimpleImplicitFunctions.H"
#include <iomanip>

#include "Chombo_NamespaceHeader.H"

#define MAX_ORDER 2

using std::cout;
using std::endl;
using std::shared_ptr;


void
dumpPPS(const PointSet* a_ivs)
{
  for(PointSetIterator ivsit(*a_ivs);  ivsit.ok(); ++ivsit)
  {
    std::cout << ivsit() << " " ;
  }
  std::cout << std::endl;
 }

int
runTest(int a_argc, char* a_argv[])
{
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
  int nIter       = 10;
  ParmParse pp;
    
  pp.get("nx"     , nx);
  pp.get("niter"  , nIter);
  pp.get("maxGrid", maxGrid);
  pp.get("x0"     , x0);
  pp.get("y0"     , y0);
  pp.get("z0"     , z0);
  pp.get("A"      , A);
  pp.get("B"      , B);
  pp.get("C"      , C);
  pp.get("R"      , R);         
  pp.get("coveredval"      , coveredval);         

  pout() << "nx      = " << nx       << endl;
  pout() << "maxGrid = " << maxGrid  << endl;
  pout() << "x0      = " << x0       << endl;
  pout() << "y0      = " << y0       << endl;
  pout() << "z0      = " << z0       << endl;
  pout() << "A       = " << A        << endl;
  pout() << "B       = " << B        << endl;
  pout() << "C       = " << C        << endl;
  pout() << "R       = " << R        << endl;

  pout() << "nIter   = " << nIter    << endl;

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

  Vector<Box> boxes;
  unsigned int blockfactor = 8;
  domainSplit(domain, boxes, maxGrid, blockfactor);
  
  Vector<int> procs;
  pout() << "making grids" << endl;
  LoadBalance(procs, boxes);
  DisjointBoxLayout grids(boxes, procs, domain);
  grids.printBalance();

  IntVect dataGhostIV =   4*IntVect::Unit;
  Point   dataGhostPt = ProtoCh::getPoint(dataGhostIV); 
  int geomGhost = 4;
  RealVect origin = RealVect::Zero();
  Real dx = 1.0/nx;
  shared_ptr<BaseIF>                       impfunc(new Proto::SimpleEllipsoidIF(ABC, X0, R, false));
//  Bx domainpr = getProtoBox(domain.domainBox());
  pout() << "defining geometry" << endl;
  shared_ptr<GeometryService<MAX_ORDER> >  geoserv(new GeometryService<MAX_ORDER>(impfunc, origin, dx, domain.domainBox(), grids, geomGhost));

  pout() << "making dictionary" << endl;
  EBDictionary<2, Real, CELL, CELL> dictionary(geoserv, grids, domain.domainBox(), dx, dataGhostPt);
  typedef EBStencil<2, Real, CELL, CELL> ebstencil_t;
  string stenname("Second_Order_Poisson");
  string dombcname("Dirichlet");
  string  ebbcname("Dirichlet");

  pout() << "registering stencil" << endl;
  dictionary.registerStencil(stenname, dombcname, ebbcname, domain.domainBox(), domain.domainBox());

  shared_ptr<LevelData<EBGraph> > graphs = geoserv->getGraphs(domain.domainBox());

  pout() << "making data" << endl;
  EBLevelBoxData<CELL,   1>  srcData(grids, dataGhostIV, graphs);
  EBLevelBoxData<CELL,   1>  dstData(grids, dataGhostIV, graphs);

  DataIterator dit = grids.dataIterator();
  pout() << "setting values" << endl;
//#pragma omp parallel for
  for(int ibox = 0; ibox < dit.size(); ibox++)
  {
    EBBoxData<CELL, Real, 1>& srcebbd = srcData[dit[ibox]];
    EBBoxData<CELL, Real, 1>& dstebbd = dstData[dit[ibox]];
    srcebbd.setVal(0.0);
    dstebbd.setVal(0.0);
  }
  //not really necessary but I wanted to see if it would work
//  Copier exchangeCopier;
//  exchangeCopier.exchangeDefine(grids, dataGhostIV);
//  srcData.exchange(exchangeCopier);

  pout() << "applying stencil" << endl;
  for(int iiter = 0; iiter < nIter; iiter++)
  {    
//#pragma omp parallel for
    for(int ibox = 0; ibox < dit.size(); ibox++)
    {
      shared_ptr<ebstencil_t> stencil = dictionary.getEBStencil(stenname, ebbcname, domain.domainBox(), domain.domainBox(), ibox);
      stencil->apply(dstData[dit[ibox]], srcData[dit[ibox]]);
    }
  }
//  srcData.writeToFileHDF5("srcData.hdf5", -1.0);
//  dstData.writeToFileHDF5("dstData.hdf5", -1.0);
  pout() << "exiting runtest" << endl;
  return 0;
}

#include "Chombo_NamespaceFooter.H"

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
    Chombo4::runTest(a_argc, a_argv);
  }

  pout() << "printing time table " << endl;
  CH_TIMER_REPORT();
#ifdef CH_MPI
  pout() << "about to call MPI Finalize" << std::endl;
  MPI_Finalize();
#endif
  return 0;
}
