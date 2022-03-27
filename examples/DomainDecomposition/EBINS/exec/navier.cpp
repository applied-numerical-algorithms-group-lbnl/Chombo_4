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
  int checkpointInterval = -1;
  int   plotfileInterval = -1;
  bool useWCycle = false;
  ParmParse pp;
  using Chombo4::pout;

  pp.get("viscosity" , nu);
  pp.get("max_step"  , max_step);
  pp.get("max_time"  , max_time);
  pp.get("checkpoint_interval", checkpointInterval);
  pp.get("plotfile_interval"  ,   plotfileInterval);
  pp.get("covered_value", coveredval);
  pp.get("num_smooth", numSmooth);
  pp.get("use_w_cycle", useWCycle);
  EBMultigridLevel::s_numSmoothUp   = numSmooth;
  EBMultigridLevel::s_numSmoothDown = numSmooth;
  EBMultigridLevel::s_useWCycle     = useWCycle;
  pout() << "max_step        = " << max_step        << endl;
  pout() << "max_time        = " << max_time        << endl;
  pout() << "checkpoint interval = " << checkpointInterval  << endl;
  pout() << "plotfile   interval = " <<   plotfileInterval  << endl;

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
//begin debug  
//  pout() << "leaving after geometry" << endl;
//  return 0;
//end debug    
  
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
  bool stokesFlowInitialization;
  pp.get("stokes_flow_initialization", stokesFlowInitialization);
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
  pp.get("checkpoint_interval", checkpointInterval);
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
  pout() << "max_step        = " << max_step   << endl;
  pout() << "max_time        = " << max_time   << endl;
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

  bool printStuff = true;
  EBINS solver(brit, geoserv, grids, domain,  dx, viscosity, dataGhostIV, 
               paraSolver, ibc, num_species, diffusion_coeffs, printStuff);


//begin debug  
  pout() << "leaving after EBINS constructor" << endl;
  return 0;
//end debug    
  unsigned int starting_step = 0;
  Real         starting_time = 0;
  string checkpointFile;
  //If the input file specifies a checkpoint restart, do that.
  if(pp.query("checkpoint_restart", checkpointFile))
  {
    //the current step and time are also stashed in the checkpoint file
    pout() << "Reading all data from a checkpoint file " << checkpointFile << endl;
    solver.readDataFromCheckpoint(starting_time, starting_step, checkpointFile);
  }
  else
  {
    //step and time already initialized to zero
    pout() << "going into initialize data " << endl;
    initializeData(solver, grids, dx, geomCen, geomRad, blobCen, blobRad, maxVelMag, maxVelRad, ibc);
  }
  if(starting_step == 0)
  {
    if(stokesFlowInitialization)
    {
      pout() << "Initializing pressure with gph = nu lapl(v)  (stokes flow initialization)" << endl;
    }
    else
    {
      if(pIters > 0)
      {
        pout() << "Using fixed point interation for initial pressure with "
               << pIters << "iterations."  << endl;
      }
      else
      {
        pout() << "Standard Treb pressure initializtion:" << endl;
        pout() << "initializing pressure with (I-P)(v*)." << endl;
        pout() << "(gphi out of initial projection).    " << endl;
      }
    }
  }
  /**
     For convergence tests and other things, fixed time steps can be useful.
     The fact that fixedDt is a negative number signals that we are using varaible dt in this case.
     This is our weird interface that says sending something non-sensical as an argument
     turns off that bit of functionality.    This quixotic interface has been standard
     for at least thirty years.   At the very least, it is not over-specified 
     and that can save a lot of code. --dtg 3-12-2021
  **/
  Real fixedDt = -1.0;

  pout() << "startiing run" << endl;
  solver.run(max_step, max_time, starting_step, starting_time,
             cfl, fixedDt, tol, pIters,  maxIter,
             plotfileInterval, checkpointInterval,
             stokesFlowInitialization, coveredval);
  pout() << "finished run" << endl;

  return 0;
}

/**********/

int main(int a_argc, char* a_argv[])
{
#ifdef CH_USE_PETSC  
  //In what must be some kind of solipsistic madness, PetscInitialize calls MPI_INIT
   PetscInt ierr = PetscInitialize(&a_argc, &a_argv, "./.petscrc",PETSC_NULL); CHKERRQ(ierr);
#else  
#ifdef CH_MPI
  MPI_Init(&a_argc, &a_argv);
  Chombo4::pout() << "MPI INIT called" << std::endl;
#endif
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
  using Chombo4::pout;
  pout() << "printing time table " << endl;
  CH_TIMER_REPORT();
#ifdef CH_MPI
#ifdef CH_USE_PETSC
  pout() << "about to call petsc Finalize" << std::endl;
  PetscFinalize();
#else  
  pout() << "about to call MPI Finalize" << std::endl;
  MPI_Finalize();
#endif
#endif
  return 0;
}


