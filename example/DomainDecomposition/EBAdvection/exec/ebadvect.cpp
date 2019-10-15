#include <cmath>
#include <cstdio>
#include <iostream>


#include "EBProto.H"
#include "EBLevelBoxData.H"
#include "LevelData.H"
#include "BaseFab.H"

#include "ParmParse.H"
#include "LoadBalance.H"
#include "ProtoInterface.H"
#include "BRMeshRefine.H"
#include "GeometryService.H"
#include "EBDictionary.H"
#include "EBChombo.H"
#include <iomanip>

#define MAX_ORDER 2

using std::cout;
using std::endl;
using std::shared_ptr;
using Proto::Var;
using Proto::SimpleEllipsoidIF;
#define PI 3.141592653589793
typedef Var<Real,DIM> Vec;
typedef Var<Real,  1> Sca;

//=================================================
PROTO_KERNEL_START 
void InitializeSpotF(int       a_p[DIM],
                     Sca       a_phi,
                     Real      a_X0,
                     Real      a_rad,
                     Real      a_dx)
{
  Real rlocsq = 0;
  for(int idir = 0; idir < DIM; idir++)
  {
    Real xrel = (a_p[idir] + 0.5)*a_dx - a_X0;
    rlocsq += xrel*xrel;
  }
  Real val = 0;
  if(rlocsq < (a_rad*a_rad))
  {
    val = a_rad*a_rad - rlocsq;
  }
  a_phi(0) = val;
}
PROTO_KERNEL_END(InitializeSpotF, InitializeSpot)

//=================================================
PROTO_KERNEL_START 
void InitializeVCellF(int       a_p[DIM],
                      Vec       a_vel,
                      Real      a_cen,
                      Real      a_rad,
                      Real      a_maxr,
                      Real      a_mag,
                      Real      a_dx)
{
  Real rlocsq = 0;
  Real xrel[DIM];
  for(int idir = 0; idir < DIM; idir++)
  {
    xrel[idir] = (a_p[idir] + 0.5)*a_dx - a_cen;
    rlocsq += xrel[idir]*xrel[idir];
  }
  Real raddiff = rlocsq-a_maxr*a_maxr;
  Real velmag = a_mag*raddiff*raddiff;
#if DIM==2
  a_vel(0) =  velmag*xrel[1];
  a_vel(1) = -velmag*xrel[0];
#else
  a_vel(0) =  velmag*( xrel[1] + xrel[2]);
  a_vel(1) =  velmag*(-xrel[0] - xrel[2]);
  a_vel(2) = -velmag*( xrel[0] - xrel[1]);
#endif

}
PROTO_KERNEL_END(InitializeVCellF, InitializeVCell)

//=================================================

void makeGrids(DisjointBoxLayout& a_grids,
               Real             & a_dx,
               const int        & a_nx,
               const int        & a_maxGrid)
{
  pout() << "making grids" << endl;

  IntVect domLo = IntVect::Zero;
  IntVect domHi  = (a_nx - 1)*IntVect::Unit;

  // EB and periodic do not mix
  ProblemDomain domain(domLo, domHi);

  Vector<Box> boxes;
  unsigned int blockfactor = 8;
  domainSplit(domain, boxes, a_maxGrid, blockfactor);
  
  Vector<int> procs;

  a_dx = 1.0/a_nx;
  LoadBalance(procs, boxes);
  a_grids = DisjointBoxLayout(boxes, procs, domain);
  a_grids.printBalance();

}

//=================================================
void computeDt(Real                        &  a_dt,
               EBLevelBoxData<CELL, DIM>   &  a_velocell,
               const   Real                &  a_dx)
{
//HERE  
}

//=================================================
void initializeData(EBLevelBoxData<CELL,   1>   &  a_scalcell,
                    EBLevelBoxData<CELL, DIM>   &  a_velocell,
                    const DisjointBoxLayout     &  a_grids,
                    const Real                  &  a_dx,
                    const Real                  &  a_geomCen,
                    const Real                  &  a_geomRad,
                    const Real                  &  a_blobCen,
                    const Real                  &  a_blobRad)
{
  MayDay::Error("not implemented");
}
//=================================================
void defineGeometry(DisjointBoxLayout& a_grids,
                    Real             & a_dx,
                    Real             & a_geomCen,
                    Real             & a_geomRad,
                    Real             & a_blobCen,
                    Real             & a_blobRad,
                    int              & a_nx,
                    shared_ptr<GeometryService<MAX_ORDER> >&  a_geoserv)
{
  pout() << "defining geometry" << endl;

  ParmParse pp;
  int maxGrid = 32;
    
  pp.get("nx"        , a_nx);
  pp.get("maxGrid", maxGrid);
  pp.get("geom_cen", a_geomCen);
  pp.get("blob_cen", a_blobCen);
  pp.get("geom_rad", a_geomRad);
  pp.get("blob_rad", a_blobRad);

  pout() << "nx       = " << a_nx     << endl;
  pout() << "maxGrid  = " << maxGrid  << endl;
  pout() << "geom_cen = " << a_geomCen       << endl;
  pout() << "geom_rad = " << a_geomRad       << endl;
  pout() << "blob_cen = " << a_blobCen       << endl;
  pout() << "blob_rad = " << a_blobRad       << endl;



  RealVect ABC = RealVect::Unit();
  RealVect  X0 = RealVect::Unit();
  X0 *= a_geomCen;


  makeGrids(a_grids, a_dx, a_nx, maxGrid);
  Box domain = a_grids.physDomain().domainBox();
  int geomGhost = 4;
  RealVect origin = RealVect::Zero();

  pout() << "creating implicit function" << endl;
  shared_ptr<BaseIF>  impfunc(new SimpleEllipsoidIF(ABC, X0, a_geomRad, true));//true is for inside regular

  pout() << "creating geometry service" << endl;
  a_geoserv  = shared_ptr<GeometryService<MAX_ORDER> >(new GeometryService<MAX_ORDER>(impfunc, origin, a_dx, domain, a_grids, geomGhost));
}

int
runAdvection(int a_argc, char* a_argv[])
{
  Real coveredval = -1;
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

  pout() << "num_streams     = " << nStream   << endl;
  pout() << "max_step        = " << max_step  << endl;
  pout() << "max_time        = " << max_time  << endl;
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
  defineGeometry(grids, dx, geomCen, geomRad, blobCen, blobRad, nx,  geoserv);

  IntVect dataGhostIV =   IntVect::Unit;
  Point   dataGhostPt = ProtoCh::getPoint(dataGhostIV); 

  pout() << "making dictionary" << endl;
  Box domain = grids.physDomain().domainBox();
  EBDictionary<2, Real, CELL, CELL> dictionary(geoserv, grids, domain, dataGhostPt, dataGhostPt, dx);


  pout() << "inititializing data"   << endl;
  
  shared_ptr<LevelData<EBGraph> > graphs = geoserv->getGraphs(domain);
  EBLevelBoxData<CELL,   1>  scalcell(grids, dataGhostIV, graphs);
  EBLevelBoxData<CELL, DIM>  velocell(grids, dataGhostIV, graphs);
  initializeData(scalcell, velocell, grids, dx, geomCen, geomRad, blobCen, blobRad);


  int step = 0; Real time = 0;
  Real dt = 0;
  pout() << "computing the time step"  << endl;
  computeDt(dt, velocell, dx);

  if(outputInterval > 0)
  {
    string filev("velo.hdf5");
    velocell.writeToFileHDF5(filev, coveredval);
    string filep("scal.0.hdf5");
    scalcell.writeToFileHDF5(filep, coveredval);
  }

  pout() << "running advection operator " << endl;
//  EBAdvectionOp advectOp(dictionary, geoserv, velocell, grids, domain,  dx, dataGhostIV, dataGhostIV);

  while((step < max_step) && (time < max_time))
  {
//    advectOp.advance(dt, scalcell);

    pout() <<" step = " << step << " time = " << time << " time step = " << dt << endl;
    if((outputInterval > 0) && (step%outputInterval == 0))
    {
      string filep = string("scal.") + std::to_string(step) + string(".hdf5");
      scalcell.writeToFileHDF5(filep, coveredval);
    }

    step++;
    time += dt;
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
  CH_TIMER_SETFILE("ebadvect.time.table");
  {
    if (a_argc < 2)
    {
      cerr<< " usage " << a_argv[0] << " <input_file_name> " << endl;
      exit(0);
    }
    char* in_file = a_argv[1];
    ParmParse  pp(a_argc-2,a_argv+2,NULL,in_file);
    runAdvection(a_argc, a_argv);
  }

  pout() << "printing time table " << endl;
  CH_TIMER_REPORT();
#ifdef CH_MPI
  pout() << "about to call MPI Finalize" << std::endl;
  MPI_Finalize();
#endif
  return 0;
}


