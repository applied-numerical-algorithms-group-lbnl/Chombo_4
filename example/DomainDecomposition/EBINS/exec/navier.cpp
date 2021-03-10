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
#include "Chombo_EBEncyclopedia.H"
#include "Chombo_EBDictionary.H"
#include "Chombo_EBChombo.H"
#include "EBINS.H"
#include "EBIBC.H"
#include "SetupFunctions.H"

#include <iomanip>

#define MAX_ORDER 2
/***/
EBIBC getIBCs()
{
  string veloIC, scalIC;
  string loDomBC[DIM];
  string hiDomBC[DIM];
  ParmParse pp;
  pp.get("initial_velo", veloIC);
  pp.get("initial_scal", scalIC);
  using std::to_string;
  for(unsigned int idir = 0; idir < DIM; idir++)
  {
    string lostr = "domain_bc_lo_" + to_string(idir);
    string histr = "domain_bc_hi_" + to_string(idir);
    pp.get(lostr.c_str(), loDomBC[idir]);
    pp.get(histr.c_str(), hiDomBC[idir]);
  }
  string ebbc("NoSlipWall");
  EBIBC retval(veloIC, scalIC, loDomBC, hiDomBC, ebbc);
  return retval;
}
/***/
int
runNavierStokes()
{

  Real coveredval = -1;
  Real nu         = -1.0;
  int nx          = 32;
  int  max_step   = 10;
  Real max_time   = 1.0;
  int numSmooth  = 4;
  int nStream    = 8;
  int outputInterval = -1;
  bool useWCycle = false;
  ParmParse pp;

  pp.get("nstream", nStream);

  pp.get("viscosity" , nu);
  pp.get("max_step"  , max_step);
  pp.get("max_time"  , max_time);
  pp.get("output_interval", outputInterval);
  pp.get("covered_value", coveredval);
  pp.get("num_smooth", numSmooth);
  pp.get("use_w_cycle", useWCycle);
  EBMultigridLevel::s_numSmoothUp   = numSmooth;
  EBMultigridLevel::s_numSmoothDown = numSmooth;
  EBMultigridLevel::s_useWCycle     = useWCycle;
  pout() << "nStream         = " << nStream         << endl;
  pout() << "max_step        = " << max_step        << endl;
  pout() << "max_time        = " << max_time        << endl;
  pout() << "output interval = " << outputInterval  << endl;

  pout() << "defining geometry" << endl;
  shared_ptr<GeometryService<MAX_ORDER> >  geoserv;

  pp.get("nx"        , nx);

  pout() << "nx       = " << nx     << endl;
  Real dx = 1.0/Real(nx);

  Vector<Chombo4::DisjointBoxLayout> vecgrids;
  Vector<Chombo4::Box>               vecdomains;
  Vector<Real> vecdx;
  int whichGeom;

  Real geomCen, geomRad;
  defineGeometry(vecgrids, vecdomains, vecdx, geoserv, geomCen, geomRad, whichGeom, dx, nx);

  IntVect dataGhostIV =   4*IntVect::Unit;
  Point   dataGhostPt = ProtoCh::getPoint(dataGhostIV); 
  
  
  pout() << "making dictionary" << endl;
  shared_ptr<EBEncyclopedia<2, Real> > 
    brit(new EBEncyclopedia<2, Real>(geoserv, vecgrids, vecdomains, vecdx, dataGhostPt));


  pout() << "inititializing data"   << endl;
  
  Chombo4::Box domain              = vecdomains[0];
  Chombo4::DisjointBoxLayout grids =   vecgrids[0];
  Real tol = 0.00001;
  int  maxIter = 10;

  Real blobCen, blobRad, maxVelMag, maxVelRad, viscosity;
  Real cfl            = 0.5;
  int pIters = 1;

  pp.get("maxIter"   , maxIter);
  pp.get("tolerance" , tol);
  pp.get("covered_value", coveredval);
  pp.get("blob_cen", blobCen);
  pp.get("blob_rad", blobRad);
  pp.get("viscosity", viscosity);
  pp.get("max_vel_mag", maxVelMag);
  pp.get("max_vel_rad", maxVelRad);
  pp.get("max_step"  , max_step);
  pp.get("max_time"  , max_time);
  pp.get("output_interval", outputInterval);
  pp.get("pressure_iterations", pIters);
  pp.get("cfl"  , cfl);
  int whichSolver;
  pp.get("parabolic_solver", whichSolver);
  EBINS::ParabolicSolverType paraSolver;
  if(whichSolver == 0)
  {
    paraSolver = EBINS::BackwardEuler;
    pout() << "using backward Euler for parabolic solver"  << endl;
  }
  else if(whichSolver == 1)
  {
    paraSolver = EBINS::CrankNicolson;
    pout() << "using Crank Nicolson for parabolic solver"  << endl;
  }
  else if(whichSolver == 2)
  {
    paraSolver = EBINS::TGA;
    pout() << "using TGA for parabolic solver"  << endl;
  }
  else
  {
    Chombo4::MayDay::Error("unrecognized solver type input");
  }

  pout() << "=============================================="  << endl;

  pout() << "tolerance       = " << tol        << endl;
  pout() << "maxIter         = " << maxIter    << endl;
  pout() << "blob cen        = " << blobCen    << endl;
  pout() << "geom cen        = " << geomCen    << endl;
  pout() << "max vel mag     = " << maxVelMag  << endl;
  pout() << "max vel rad     = " << maxVelRad  << endl;
  pout() << "num_streams     = " << nStream    << endl;
  pout() << "max_step        = " << max_step   << endl;
  pout() << "max_time        = " << max_time   << endl;
  pout() << "pressure iter   = " << pIters     << endl;
  pout() << "viscosity       = " << viscosity  << endl;
  pout() << "=============================================="  << endl;

  
  EBIBC ibc = getIBCs();
  pout() << "initializing solver " << endl;
  int num_species = 0;
  pp.query("num_species", num_species);
  vector<Real> diffusion_coeffs(num_species);
  for(int ispec = 0; ispec < num_species; ispec++)
  {
    string diffname = string("diffusion_coeff_") + to_string(ispec);
    Real thisco;
    pp.get(diffname.c_str(), thisco);
    diffusion_coeffs[ispec] = thisco;
  }
  
  EBINS solver(brit, geoserv, grids, domain,  dx, viscosity, dataGhostIV, 
               paraSolver, ibc, num_species, diffusion_coeffs);


   pout() << "initializing data " << endl;
   initializeData(solver, grids, dx, geomCen, geomRad, blobCen, blobRad, maxVelMag, maxVelRad, ibc);

  Real fixedDt = -1.0;//signals varaible dt

  solver.run(max_step, max_time, cfl, fixedDt, tol, pIters,  maxIter, outputInterval, coveredval);
  pout() << "exiting run" << endl;
  return 0;
}

/**********/

int main(int a_argc, char* a_argv[])
{
#ifdef CH_MPI
  MPI_Init(&a_argc, &a_argv);
  pout() << "MPI INIT called" << std::endl;
#endif
  //needs to be called after MPI_Init
  CH_TIMER_SETFILE("navier.time.table");
  {
    if (a_argc < 2)
    {
      cerr<< " usage " << a_argv[0] << " <input_file_name> " << endl;
      exit(0);
    }
    char* in_file = a_argv[1];
    ParmParse  pp(a_argc-2,a_argv+2,NULL,in_file);
    runNavierStokes();
  }

  pout() << "printing time table " << endl;
  CH_TIMER_REPORT();
#ifdef CH_MPI
  pout() << "about to call MPI Finalize" << std::endl;
  MPI_Finalize();
#endif
  return 0;
}


