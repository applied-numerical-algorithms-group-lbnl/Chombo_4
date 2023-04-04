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
#include "EBPetscSolver.H"
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
    
  int dombc = 1;
  int ebbc  = 1;
  ParmParse pp;
    
  pp.get("domainBC"  , dombc);
  pp.get("EBBC"      , ebbc);
  pp.get("nx"   , nx);
  pp.get("maxGrid"   , maxGrid);
  pp.get("x0"        , x0);
  pp.get("y0"        , y0);
  pp.get("z0"        , z0);
  pp.get("A"         , A);
  pp.get("B"         , B);
  pp.get("C"         , C);
  pp.get("R"         , R);         
  pp.get("coveredval", coveredval);         

  Chombo4::pout() << "nx        = " << nx       << endl;
  Chombo4::pout() << "maxGrid   = " << maxGrid  << endl;
  Chombo4::pout() << "x0        = " << x0       << endl;
  Chombo4::pout() << "y0        = " << y0       << endl;
  Chombo4::pout() << "z0        = " << z0       << endl;
  Chombo4::pout() << "A         = " << A        << endl;
  Chombo4::pout() << "B         = " << B        << endl;
  Chombo4::pout() << "C         = " << C        << endl;
  Chombo4::pout() << "R         = " << R        << endl;



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
  string alldombc[2*DIM];
  for(int iface = 0; iface < 2*DIM; iface++)
  {
    alldombc[iface] = dombcname;
  }
  Chombo4::Box dombox = domain.domainBox();
  auto graphs = geoserv->getGraphs(dombox);

  Real alpha, beta;
  pp.get("alpha", alpha);
  pp.get("beta",  beta);
  EBPetscSolver<2> solver(geoserv, dictionary, graphs, grids, domain.domainBox(),
                          stenname, alldombc, ebbcname,
                          dx, alpha, beta, dataGhostPt);
  Chombo4::pout() << "making data" << endl;
  EBLevelBoxData<CELL,   1>  phi(grids, dataGhostIV, graphs);
  EBLevelBoxData<CELL,   1>  rhs(grids, dataGhostIV, graphs);
  EBLevelBoxData<CELL,   1>  res(grids, dataGhostIV, graphs);

  Chombo4::pout() << "setting values" << endl;
  Chombo4::DataIterator dit = grids.dataIterator();
  for(int ibox = 0; ibox < dit.size(); ibox++)
  {
    EBBoxData<CELL, Real, 1>& phibd = phi[dit[ibox]];
    EBBoxData<CELL, Real, 1>& rhsbd = rhs[dit[ibox]];
    phibd.setVal(0.0);
    rhsbd.setVal(1.0);
  }

  solver.solve(phi, rhs);

  EBMultigrid ebmg(dictionary, geoserv, alpha, beta, dx, grids,
                   stenname, dombcname, ebbcname, domain.domainBox(),
                   dataGhostIV, string("testmg"), false);

  ebmg.residual(res, phi, rhs);
  Real errnorm = res.maxNorm(0);
  
  Chombo4::pout() << "norm of computed residual =  " << errnorm << endl;
  Chombo4::pout() << "writing to file " << endl;

  phi.writeToFileHDF5("phi.hdf5", 0.0);
  rhs.writeToFileHDF5("rhs.hdf5", 0.0);
  res.writeToFileHDF5("res.hdf5", 0.0);
  Chombo4::pout() << "exiting " << endl;
  return 0;
}


int main(int a_argc, char* a_argv[])
{
  PetscInt ierr;
  //because of some kind of solipsistic madness, PetscInitialize calls MPI_INIT
  ierr = PetscInitialize(&a_argc, &a_argv, "./.petscrc",PETSC_NULL); CHKERRQ(ierr); 
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

  Chombo4::pout() << "printing time table " << endl;
  ierr = PetscFinalize(); CHKERRQ(ierr); 
  CH_TIMER_REPORT();
  return 0;
}
