#include <cmath>
#include <cstdio>
#include <iostream>


#include "EBProto.H"
#include "Chombo_EBLevelBoxData.H"

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
  int outputInterval = -1;
  bool useWCycle = false;
  ParmParse pp;

  pp.get("permeability" , permeability);
  pp.get("diffusivity" ,  diffusivity);
  pp.get("max_step"  , max_step);
  pp.get("max_time"  , max_time);
  pp.get("output_interval", outputInterval);
  pp.get("covered_value", coveredval);
  pp.get("num_smooth", numSmooth);
  pp.get("use_w_cycle", useWCycle);
  EBMultigridLevel::s_numSmoothUp   = numSmooth;
  EBMultigridLevel::s_numSmoothDown = numSmooth;
  EBMultigridLevel::s_useWCycle     = useWCycle;
  Chombo4::pout() << "max_step        = " << max_step        << endl;
  Chombo4::pout() << "max_time        = " << max_time        << endl;
  Chombo4::pout() << "output interval = " << outputInterval  << endl;

  Chombo4::pout() << "defining geometry" << endl;
  shared_ptr<GeometryService<MAX_ORDER> >  geoserv;

  pp.get("nx"        , nx);

  Chombo4::pout() << "nx       = " << nx     << endl;
  Real dx = 1.0/Real(nx);

  std::vector<Chombo4::DisjointBoxLayout> vecgrids;
  std::vector<Chombo4::Box>               vecdomains;
  std::vector<Real> vecdx;
  int whichGeom;

  Real geomCen, geomRad;
  defineGeometry(vecgrids, vecdomains, vecdx, geoserv, geomCen, geomRad, whichGeom, dx, nx);

  IntVect dataGhostIV =   4*IntVect::Unit;
  Point   dataGhostPt = ProtoCh::getPoint(dataGhostIV); 
  
  
  Chombo4::pout() << "making dictionary" << endl;
  shared_ptr<EBEncyclopedia<2, Real> > 
    brit(new EBEncyclopedia<2, Real>(geoserv, vecgrids, vecdomains, vecdx, dataGhostPt));


  Chombo4::pout() << "inititializing data"   << endl;
  
  Chombo4::Box domain              = vecdomains[0];
  Chombo4::DisjointBoxLayout grids =   vecgrids[0];
  Real tol = 0.00001;
  int  maxIter = 10;

  Real blobCen, blobRad, viscosity;
  Real cfl            = 0.5;

  pp.get("maxIter"   , maxIter);
  pp.get("tolerance" , tol);
  pp.get("covered_value", coveredval);
  pp.get("blob_cen", blobCen);
  pp.get("blob_rad", blobRad);
  pp.get("viscosity", viscosity);
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
    Chombo4::pout() << "using backward Euler for parabolic solver"  << endl;
  }
  else if(whichSolver == 1)
  {
    paraSolver = EBDarcy::CrankNicolson;
    Chombo4::pout() << "using Crank Nicolson for parabolic solver"  << endl;
  }
  else if(whichSolver == 2)
  {
    paraSolver = EBDarcy::TGA;
    Chombo4::pout() << "using TGA for parabolic solver"  << endl;
  }
  else
  {
    Chombo4::MayDay::Error("unrecognized solver type input");
  }

  Chombo4::pout() << "=============================================="  << endl;

  Chombo4::pout() << "tolerance       = " << tol          << endl;
  Chombo4::pout() << "maxIter         = " << maxIter      << endl;
  Chombo4::pout() << "blob cen        = " << blobCen      << endl;
  Chombo4::pout() << "geom cen        = " << geomCen      << endl;
  Chombo4::pout() << "max_step        = " << max_step     << endl;
  Chombo4::pout() << "max_time        = " << max_time     << endl;
  Chombo4::pout() << "viscosity       = " << viscosity    << endl;
  Chombo4::pout() << "diffusivity     = " << diffusivity  << endl;
  Chombo4::pout() << "pereability     = " << permeability << endl;
  Chombo4::pout() << "=============================================="  << endl;

  
  Chombo4::pout() << "initializing solver " << endl;
  Chombo4::pout() << "inflow outflow xdirection, noflow all other directions" << endl;
  //Here I am creating an EBIBC to show the solvers what the boundary conditions are.
  //This application forces this particular bit of the problem so this should not be
  //part of the public interface.   It does allow for a lot of code reuse, however,
  //so I remain unapologetic.  --dtg
  string  loDomainBC[DIM];
  string  hiDomainBC[DIM];
  loDomainBC[0] = string("inflow");
  hiDomainBC[0] = string("inflow");
  loDomainBC[1] = string("slip_wall");
  hiDomainBC[1] = string("slip_wall");
#if DIM==3  
  loDomainBC[2] = string("slip_wall");
  hiDomainBC[2] = string("slip_wall");
#endif
  string ebbc("no_slip_wall");       //just to put something in there (unneeded)
  EBIBC ibcs(string("no_velo_ic"),   //just to put something in there (unneeded)
             string("some_scal_ic"), //just to put something in there (hardwired for now)
             loDomainBC, hiDomainBC);
                    
  EBDarcy solver(brit, geoserv, grids, domain,  dx,
                 viscosity, permeability, diffusivity,
                 dataGhostIV, paraSolver, ibcs);


// auto &  velo = *(solver.m_velo);
 auto &  scal = *(solver.m_scal);

 
  Chombo4::pout() << "initializing data " << endl;
  initializeData(scal, grids, dx, geomCen, geomRad, blobCen, blobRad);

  Real fixedDt = -1.0;//signals varaible dt

  Chombo4::pout() << "starting      EBDarcy::run" << endl; 
  solver.run(max_step, max_time, cfl, fixedDt, tol,   maxIter, outputInterval, coveredval);
  Chombo4::pout() << "finished with EBDarcy::run" << endl; 
  return 0;
}

/**********/

int main(int a_argc, char* a_argv[])
{
#ifdef CH_MPI
  MPI_Init(&a_argc, &a_argv);
  Chombo4::pout() << "MPI INIT called" << std::endl;
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

  Chombo4::pout() << "printing time table " << endl;
  CH_TIMER_REPORT();
#ifdef CH_MPI
  Chombo4::pout() << "about to call MPI Finalize" << std::endl;
  MPI_Finalize();
#endif
  return 0;
}


