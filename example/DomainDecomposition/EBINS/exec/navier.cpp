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
#include "EBAdvection.H"
#include "SetupFunctions.H"

#include <iomanip>

#define MAX_ORDER 2


int
runNavierStokes(int a_argc, char* a_argv[])
{

  Real coveredval = -1;
  Real nu         = -1.0;
  int nx          = 32;
  int  max_step   = 10;
  Real max_time   = 1.0;

  int nStream    = 8;
  int outputInterval = -1;
  ParmParse pp;

  pp.get("nStream", nStream);

  pp.get("viscosity" , nu);
  pp.get("max_step"  , max_step);
  pp.get("max_time"  , max_time);
  pp.get("output_interval", outputInterval);
  pp.get("covered_value", coveredval);

  pout() << "nStream         = " << nStream         << endl;
  pout() << "max_step        = " << max_step        << endl;
  pout() << "max_time        = " << max_time        << endl;
  pout() << "output interval = " << outputInterval  << endl;

#ifdef PROTO_CUDA
  Proto::DisjointBoxLayout::setNumStreams(nStream);
#endif



  pout() << "defining geometry" << endl;
  shared_ptr<GeometryService<MAX_ORDER> >  geoserv;

  Real diffCoef;
  pp.get("nx"        , nx);

  pout() << "nx       = " << nx     << endl;
  Real dx = 1.0/Real(nx);

  Vector<DisjointBoxLayout> vecgrids;
  Vector<Box>               vecdomains;
  Vector<Real> vecdx;
  defineGeometry(vecgrids, vecdomains, vecdx, geoserv, dx, nx);

  IntVect dataGhostIV =   4*IntVect::Unit;
  Point   dataGhostPt = ProtoCh::getPoint(dataGhostIV); 
  
  
  pout() << "making dictionary" << endl;
  shared_ptr<EBEncyclopedia<2, Real> > 
    brit(new EBEncyclopedia<2, Real>(geoserv, vecgrids, vecdomains, vecdx, dataGhostPt));


  pout() << "inititializing data"   << endl;
  
  Box domain              = vecdomains[0];
  DisjointBoxLayout grids =   vecgrids[0];
  Real tol = 0.00001;
  int  maxIter = 10;
  Real coveredval = -0.0001;
  Real blobCen, blobRad, maxVelMag, maxVelRad;
  int  max_step       = 10;
  Real max_time       = 1.0;
  int  outputInterval = -1;
  Real cfl            = 0.5;

  pp.get("maxIter"   , maxIter);
  pp.get("tolerance" , tol);
  pp.get("covered_value", coveredval);
  pp.get("blob_cen", blobCen);
  pp.get("geom_cen", geomCen);
  pp.get("max_vel_mag", maxVelMag);
  pp.get("max_vel_rad", maxVelRad);
  pp.get("max_step"  , max_step);
  pp.get("max_time"  , max_time);
  pp.get("output_interval", outputInterval);
  pp.get("cfl"  , cfl);

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
  pout() << "=============================================="  << endl;

  pout() << "initializing solver " << endl;
  EBINS solver(brit, geoserv, grids, domain,  dx, tol, 
               maxIter, coveredval,dataGhostIV);

  shared_ptr<EBLevelBoxData<CELL, DIM> velo = solver.m_velo;
  shared_ptr<EBLevelBoxData<CELL,   1> scal = solver.m_scal;

  
  pout() << "initializing data " << endl;
  initializeData(scal, velo, grids, dx, geomCen, geomRad, blobCen, blobRad, maxVelMag, maxVelRad);

  solver.run(max_step, max_time, cfl, outputInterval);
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
  CH_TIMER_SETFILE("ebadvect.time.table");
  {
    if (a_argc < 2)
    {
      cerr<< " usage " << a_argv[0] << " <input_file_name> " << endl;
      exit(0);
    }
    char* in_file = a_argv[1];
    ParmParse  pp(a_argc-2,a_argv+2,NULL,in_file);
    runNavierStokes(a_argc, a_argv);
  }

  pout() << "printing time table " << endl;
  CH_TIMER_REPORT();
#ifdef CH_MPI
  pout() << "about to call MPI Finalize" << std::endl;
  MPI_Finalize();
#endif
  return 0;
}


