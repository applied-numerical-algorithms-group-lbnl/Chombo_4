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

#include "DebugFunctions.H"
#include <iomanip>


int
runTest(int a_argc, char* a_argv[])
{
  using Chombo4::pout;
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

  std::vector<Chombo4::DisjointBoxLayout> vecgrids;
  Chombo4::pout() << "making grids" << endl;
  GeometryService<2>::generateGrids(vecgrids, domain.domainBox(), maxGrid);

  Chombo4::DisjointBoxLayout grids = vecgrids[0];
  grids.printBalance();

  IntVect dataGhostIV =   2*IntVect::Unit;
  Point   dataGhostPt = ProtoCh::getPoint(dataGhostIV); 
  int geomGhost = 4;
  RealVect origin = RealVect::Zero();
  Real dx = 1.0/nx;
//  Real dx = 1.0;
  bool insideRegular = false;
  pp.get("inside_regular", insideRegular);
                          
  shared_ptr<BaseIF>    impfunc(new Proto::SimpleEllipsoidIF(ABC, X0, R, insideRegular));
//  Bx domainpr = getProtoBox(domain.domainBox());

  Chombo4::pout() << "defining geometry" << endl;
  GeometryService<2>* geomptr = new GeometryService<2>(impfunc, origin, dx, domain.domainBox(), vecgrids, geomGhost);
//  GeometryService<2>* geomptr = new GeometryService<2>(impfunc, origin, dx, domain.domainBox(), vecgrids[0], geomGhost);
  shared_ptr< GeometryService<2> >  geoserv(geomptr);

  Chombo4::pout() << "making dictionary" << endl;

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
    Chombo4::pout() << "using Neumann BCs at domain" << endl;
  }
  else
  {
    dombcname = StencilNames::Dirichlet;
    Chombo4::pout() << "using Dirichlet BCs at domain" << endl;
  }

  if(ebbc == 0)
  {
    ebbcname = StencilNames::Neumann;
    Chombo4::pout() << "using Neumann BCs at EB" << endl;
  }
  else
  {
    ebbcname = StencilNames::Dirichlet;
    Chombo4::pout() << "using Dirichlet BCs at EB" << endl;
  }
  Chombo4::Box dombox = domain.domainBox();
  auto graphs = geoserv->getGraphs(dombox);

  Chombo4::pout() << "making data" << endl;
  EBLevelBoxData<CELL,   1>  phi(grids, dataGhostIV, graphs);
  EBLevelBoxData<CELL,   1>  rhs(grids, dataGhostIV, graphs);
  EBLevelBoxData<CELL,   1>  res(grids, dataGhostIV, graphs);
  EBLevelBoxData<CELL,   1>  cor(grids, dataGhostIV, graphs);

  bool direct_to_bottom = false;
  bool useWCycle = false;
  bool printStuff = true;
  string prefix;
  string bottom_solver;
  pp.get("direct_to_bottom", direct_to_bottom);
  pp.get("bottom_solver", bottom_solver);
  pp.get("useWCycle", useWCycle);
  EBMultigrid solver(dictionary, geoserv, alpha, beta, dx, grids,
                     stenname, dombcname, ebbcname, dombox, dataGhostIV,
                     bottom_solver, direct_to_bottom, prefix, useWCycle, numSmooth,
                     printStuff);
  

  Chombo4::DataIterator dit = grids.dataIterator();
  Chombo4::pout() << "setting values" << endl;
  int numsolves = 1;
  pp.query("num_solves", numsolves);
  Chombo4::pout() << "number of solves = " <<  numsolves << endl;
  for(int isolve =0; isolve < numsolves; isolve++)
  {
    Chombo4::pout() << "going into solve number " << isolve << endl;
    for(int ibox = 0; ibox < dit.size(); ibox++)
    {
      EBBoxData<CELL, Real, 1>& phibd = phi[dit[ibox]];
      EBBoxData<CELL, Real, 1>& rhsbd = rhs[dit[ibox]];
      EBBoxData<CELL, Real, 1>& corbd = cor[dit[ibox]];
      phibd.setVal(0.0);
      rhsbd.setVal(1.0);
      corbd.setVal(0.0);
    }

    bool initToZero = true;
    bool printStuff = false;//turn on for debugging
    solver.solve(phi, rhs, tol, maxIter, initToZero, printStuff);
  }
  Chombo4::pout() << "done with solves " << endl;
  
  
  CH_TIMER_REPORT();
  Chombo4::pout() << "exiting " << endl;
  return 0;
}


int main(int a_argc, char* a_argv[])
{
#ifdef CH_USE_PETSC  
  //because of some kind of solipsistic madness, PetscInitialize calls MPI_INIT
   PetscInt ierr = PetscInitialize(&a_argc, &a_argv, "./.petscrc",PETSC_NULL); CHKERRQ(ierr);
#else  
#ifdef CH_MPI
  MPI_Init(&a_argc, &a_argv);
  Chombo4::pout() << "MPI INIT called" << std::endl;
#endif
#endif
  
  //needs to be called after MPI_Init
  CH_TIMER_SETFILE("helmholtz.time.table");
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
  using Chombo4::pout;
  Chombo4::pout() << "printing time table " << endl;
  CH_TIMER_REPORT();
#ifdef CH_MPI
#ifdef CH_USE_PETSC
  Chombo4::pout() << "about to call petsc Finalize" << std::endl;
  PetscFinalize();
#else  
  Chombo4::pout() << "about to call MPI Finalize" << std::endl;
  MPI_Finalize();
#endif
#endif
  return 0;
}
