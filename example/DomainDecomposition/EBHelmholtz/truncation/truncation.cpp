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
#include "Chombo_EBDictionary.H"
#include "Chombo_EBChombo.H"
#include "EBMultigrid.H"
#include "Proto_DebugHooks.H"
#include "DebugFunctions.H"
#include <iomanip>


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
  maxGrid = nx;
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
  Chombo4::ProblemDomain domain(domLo, domHi);

  Vector<Chombo4::DisjointBoxLayout> vecgrids;
  pout() << "making grids" << endl;
  GeometryService<2>::generateGrids(vecgrids, domain.domainBox(), maxGrid);

  Chombo4::DisjointBoxLayout grids = vecgrids[0];
  grids.printBalance();

  IntVect dataGhostIV =   2*IntVect::Unit;
  Point   dataGhostPt = ProtoCh::getPoint(dataGhostIV); 
  int geomGhost = 4;
  RealVect origin = RealVect::Zero();
  Real dx = 1.0/nx;
//  Real dx = 1.0;
  shared_ptr<BaseIF>    impfunc(new Proto::SimpleEllipsoidIF(ABC, X0, R, true));
//  Bx domainpr = getProtoBox(domain.domainBox());

  pout() << "defining geometry" << endl;
  GeometryService<2>* geomptr = new GeometryService<2>(impfunc, origin, dx, domain.domainBox(), vecgrids, geomGhost);
//  GeometryService<2>* geomptr = new GeometryService<2>(impfunc, origin, dx, domain.domainBox(), vecgrids[0], geomGhost);
  shared_ptr< GeometryService<2> >  geoserv(geomptr);

  pout() << "making dictionary" << endl;

  vector<Chombo4::Box>    vecdomain(vecgrids.size(), domain.domainBox());
  vector<Real>   vecdx    (vecgrids.size(), dx);
  for(int ilev = 1; ilev < vecgrids.size(); ilev++)
  {
    vecdomain[ilev] = coarsen(vecdomain[ilev-1], 2);
    vecdx    [ilev] =           2*vecdx[ilev-1];
  }
  shared_ptr<EBDictionary<2, Real, CELL, CELL> > 
    dictionary(new EBDictionary<2, Real, CELL, CELL>(geoserv, vecgrids, vecdomain, vecdx, dataGhostPt));
  
//  typedef EBStencil<2, Real, CELL, CELL> ebstencil_t;
  string stenname = StencilNames::Poisson2;
  string dombcname, ebbcname;
  if(dombc == 0)
  {
    dombcname = StencilNames::Neumann;
    pout() << "using Neumann BCs at domain" << endl;
  }
  else
  {
    dombcname = StencilNames::Dirichlet;
    pout() << "using Dirichlet BCs at domain" << endl;
  }

  if(ebbc == 0)
  {
    ebbcname = StencilNames::Neumann;
    pout() << "using Neumann BCs at EB" << endl;
  }
  else
  {
    ebbcname = StencilNames::Dirichlet;
    pout() << "using Dirichlet BCs at EB" << endl;
  }
  Chombo4::Box dombox = domain.domainBox();
  shared_ptr<LevelData<EBGraph> > graphs = geoserv->getGraphs(dombox);

  pout() << "making data" << endl;
  EBLevelBoxData<CELL,   1>  phi(grids, dataGhostIV, graphs);
  EBLevelBoxData<CELL,   1>  rhs(grids, dataGhostIV, graphs);
  EBLevelBoxData<CELL,   1>  res(grids, dataGhostIV, graphs);
  EBLevelBoxData<CELL,   1>  cor(grids, dataGhostIV, graphs);

  EBMultigrid solver(dictionary, geoserv, alpha, beta, dx, grids, stenname, dombcname, ebbcname, dombox, dataGhostIV);
  EBMultigrid::s_numSmoothUp   = numSmooth;
  EBMultigrid::s_numSmoothDown = numSmooth;
  Chombo4::DataIterator dit = grids.dataIterator();
  pout() << "setting values" << endl;
  for(int ibox = 0; ibox < dit.size(); ibox++)
  {
    EBBoxData<CELL, Real, 1>& phibd = phi[dit[ibox]];
    EBBoxData<CELL, Real, 1>& rhsbd = rhs[dit[ibox]];
    EBBoxData<CELL, Real, 1>& corbd = cor[dit[ibox]];
    phibd.setVal(0.0);
    rhsbd.setVal(1.0);
    corbd.setVal(0.0);
  }

  solver.solve(phi, rhs, tol, maxIter, false);

  pout() << "writing to file " << endl;
  
//  auto& kappa = solver.getKappa();
//  writeEBLevelHDF5<1>(string("phi.hdf5"), phi, kappa, dombox, graphs, coveredval, dx, 1.0, 0.0);
//  writeEBLevelHDF5<1>(string("rhs.hdf5"), rhs, kappa, dombox, graphs, coveredval, dx, 1.0, 0.0);
//  writeEBLevelHDF5<1>(string("res.hdf5"), res, kappa, dombox, graphs, coveredval, dx, 1.0, 0.0);
  
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
  CH_TIMER_SETFILE("trunc.time.table");
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
