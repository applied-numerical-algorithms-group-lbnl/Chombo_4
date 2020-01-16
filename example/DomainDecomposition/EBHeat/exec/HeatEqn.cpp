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

using std::cout;
using std::endl;
using std::shared_ptr;
using Proto::Var;
using Proto::SimpleEllipsoidIF;

typedef Var<Real,DIM> Vec;
typedef Var<Real,  1> Sca;

//=================================================
void 
initializeSource(EBLevelBoxData<CELL,   1>   &  a_source,
                 const DisjointBoxLayout     &  a_grids,
                 const Real                  &  a_dx,
                 const Real                  &  a_geomCen,
                 const Real                  &  a_geomRad,
                 const Real                  &  a_blobCen,
                 const Real                  &  a_blobRad)

{
  DataIterator dit = a_grids.dataIterator();
  int ideb = 0;
  pout() << "calling initializespot for scalar rhs " << endl;

  for(int ibox = 0; ibox < dit.size(); ibox++)
  {
    auto& scalfab = a_source[dit[ibox]];
    unsigned long long int numflopsscal = 5*DIM +3;
    Bx scalbox = scalfab.box();

    ebforallInPlace_i(numflopsscal, "IntializeSpot", InitializeSpot,  scalbox, 
                      scalfab, a_blobCen, a_blobRad, a_dx);
  }
}

int
runHeatEqn(int a_argc, char* a_argv[])
{
  Real coveredval = -1;
  Real cfl    = 0.5;
  int nx      = 32;
  int  max_step   = 10;
  Real max_time   = 1.0;

  int nStream    = 8;
  int outputInterval = -1;
  ParmParse pp;

  pp.get("nstream", nStream);

    
  pp.get("max_step"  , max_step);
  pp.get("max_time"  , max_time);
  pp.get("output_interval", outputInterval);
  pp.get("covered_value", coveredval);

  pout() << "num_streams     = " << nStream         << endl;
  pout() << "max_step        = " << max_step        << endl;
  pout() << "max_time        = " << max_time        << endl;
  pout() << "output interval = " << outputInterval  << endl;

#ifdef PROTO_CUDA
  Proto::DisjointBoxLayout::setNumStreams(nStream);
#endif

  Real dx;
  DisjointBoxLayout grids;

  pout() << "defining geometry" << endl;
  shared_ptr<GeometryService<MAX_ORDER> >  geoserv;

  Real geomCen;
  Real geomRad;
  Real blobCen;
  Real blobRad;
  int whichGeom;
  defineGeometry(grids, dx, geomCen, geomRad, blobCen, blobRad, whichGeom, nx,  geoserv);

  IntVect dataGhostIV =   4*IntVect::Unit;
  Point   dataGhostPt = ProtoCh::getPoint(dataGhostIV); 

  pout() << "making dictionary" << endl;
  Box domain = grids.physDomain().domainBox();
  shared_ptr<EBEncyclopedia<2, Real> > 
    brit(new EBEncyclopedia<2, Real>(geoserv, grids, domain, dx, dataGhostPt));


  pout() << "inititializing data"   << endl;
  
  shared_ptr<LevelData<EBGraph> > graphs = geoserv->getGraphs(domain);
  EBLevelBoxData<CELL,   1>  scalcell(grids, dataGhostIV, graphs);
  initializeData(scalcell,  grids, dx, blobCen, blobRad);


  int step = 0; Real time = 0;
  Real dt = dx;
  pout() << "setting the dt = dx" << endl;

  EBAdvection advectOp(brit, geoserv, velocell, grids, domain,  dx, dataGhostIV, dataGhostIV);
  const EBLevelBoxData<CELL, 1> & kappa = advectOp.getKappa();

  if(outputInterval > 0)
  {
    string filep("scal.0.hdf5");
    writeEBLevelHDF5<1>(  filep,  scalcell, kappa, domain, graphs, coveredval, dx, dt, time);
  }

  pout() << "running advection operator " << endl;

  while((step < max_step) && (time < max_time))
  {
    advectOp.advance(scalcell, dt);

    pout() <<" step = " << step << " time = " << time << " time step = " << dt << endl;
    step++;
    time += dt;

    if((outputInterval > 0) && ( (step%outputInterval == 0) || step == (max_step-1)))
    {
      string filep = string("scal.") + std::to_string(step) + string(".hdf5");
      writeEBLevelHDF5<1>(  filep,  scalcell, kappa, domain, graphs, coveredval, dx, dt, time);
    }
  }
  pout() << "exiting runAdvection" << endl;
  return 0;
}


int main(int a_argc, char* a_argv[])
{
#ifdef CH_MPI
  MPI_Init(&a_argc, &a_argv);
  pout() << "MPI INIT called" << std::endl;
#endif
  //needs to be called after MPI_Init
  CH_TIMER_SETFILE("heateqn.time.table");
  {
    if (a_argc < 2)
    {
      cerr<< " usage " << a_argv[0] << " <input_file_name> " << endl;
      exit(0);
    }
    char* in_file = a_argv[1];
    ParmParse  pp(a_argc-2,a_argv+2,NULL,in_file);
    runHeatEqn(a_argc, a_argv);
  }

  pout() << "printing time table " << endl;
  CH_TIMER_REPORT();
#ifdef CH_MPI
  pout() << "about to call MPI Finalize" << std::endl;
  MPI_Finalize();
#endif
  return 0;
}


