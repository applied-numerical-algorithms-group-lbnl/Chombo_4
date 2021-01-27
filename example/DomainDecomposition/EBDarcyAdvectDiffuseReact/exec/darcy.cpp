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
#include "EBDarcy.H"
#include "SetupFunctions.H"

#include <iomanip>

#define MAX_ORDER 2
/***/
int
runNavierStokes()
{

  Real coveredval = -1;
  Real permeability       = -1.0;
  Real diffusivity        = -1.0;
  int nx          = 32;
  int  max_step   = 10;
  Real max_time   = 1.0;
  int numSmooth  = 4;
  int nStream    = 8;
  int outputInterval = -1;
  bool useWCycle = false;
  ParmParse pp;

  pp.get("nstream", nStream);

  pp.get("permeability" , permeability);
  pp.get("diffusivity" ,  diffusivity);
  pp.get("max_step"  , max_step);
  pp.get("max_time"  , max_time);
  pp.get("output_interval", outputInterval);
  pp.get("covered_value", coveredval);
  pp.get("num_smooth", numSmooth);
  pp.get("use_w_cycle", useWCycle);
  EBMultigrid::s_numSmoothUp   = numSmooth;
  EBMultigrid::s_numSmoothDown = numSmooth;
  EBMultigrid::s_useWCycle     = useWCycle;
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
  pp.get("cfl"  , cfl);
  int whichSolver;
  pp.get("parabolic_solver", whichSolver);
  EBDarcy::ParabolicSolverType paraSolver;
  if(whichSolver == 0)
  {
    paraSolver = EBDarcy::BackwardEuler;
    pout() << "using backward Euler for parabolic solver"  << endl;
  }
  else if(whichSolver == 1)
  {
    paraSolver = EBDarcy::CrankNicolson;
    pout() << "using Crank Nicolson for parabolic solver"  << endl;
  }
  else if(whichSolver == 2)
  {
    paraSolver = EBDarcy::TGA;
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

  
  pout() << "initializing solver " << endl;
  pout() << "inflow outflow xdirection, noflow all other directions" << endl;
  EBDarcy solver(brit, geoserv, grids, domain,  dx,
                 viscosity, permerability, diffusivity,
                 dataGhostIV, paraSolver);


 auto &  velo = *(solver.m_velo);
 auto &  scal = *(solver.m_scal);


  pout() << "initializing data " << endl;
  initializeData(scal, velo, grids, dx, geomCen, geomRad, blobCen, blobRad, maxVelMag, maxVelRad);

  Real fixedDt = -1.0;//signals varaible dt

  solver.run(max_step, max_time, cfl, fixedDt, tol,   maxIter, outputInterval, coveredval);
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


