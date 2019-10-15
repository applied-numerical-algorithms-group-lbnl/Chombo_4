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


#define PI 3.141592653589793
typedef Var<Real,DIM> Vec;
typedef Var<Real,  1> Sca;

//=================================================
PROTO_KERNEL_START 
unsigned int InitializeSpotF(int       a_p[DIM],
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
unsigned int InitializeVCellF(int       a_p[DIM],
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
    rlocsq += xrel*xrel;
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
               const int        & a_nx)
{
  pout() << "making grids" << endl;

  IntVect domLo = IntVect::Zero;
  IntVect domHi  = (a_nx - 1)*IntVect::Unit;

  // EB and periodic do not mix
  ProblemDomain domain(domLo, domHi);

  Vector<Box> boxes;
  unsigned int blockfactor = 8;
  domainSplit(domain, boxes, maxGrid, blockfactor);
  
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
                    Real 
                    Real             & a_geomCen,
                    Real             & a_geomRad,
                    Real             & a_blobCen,
                    Real             & a_blobRad)
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
                    int              & a_max_step,
                    int              & a_max_time,
                    shared_ptr<GeometryService<MAX_ORDER> >&  a_geoserv)
{
  pout() << "defining geometry" << endl;

  Real A = 1.0;
  Real B = 1.0;
  Real C = 1.0;
  int nStream    = 8;
  ParmParse pp;
  int maxGrid = 32;
    
  pp.get("nx"        , a_nx);
  pp.get("max_step"  , a_max_step);
  pp.get("max_time"  , a_max_time);
  pp.get("nstream", nStream);
  pp.get("maxGrid", maxGrid);
  pp.get("geom_cen", a_geomCen);
  pp.get("blob_cen", a_blobCen);
  pp.get("geom_rad", a_geomRad);
  pp.get("blob_rad", a_blobRad);

  pp.get("A"      , A);
  pp.get("B"      , B);
  pp.get("C"      , C);


  pout() << "nx       = " << a_nx     << endl;
  pout() << "maxGrid  = " << maxGrid  << endl;
  pout() << "geom_cen = " << a_geomCen       << endl;
  pout() << "geom_rad = " << a_geomRad       << endl;
  pout() << "blob_cen = " << a_blobCen       << endl;
  pout() << "blob_rad = " << a_blobRad       << endl;

  pout() << "A        = " << A        << endl;
  pout() << "B        = " << B        << endl;
  pout() << "C        = " << C        << endl;

  pout() << "max_step= " << a_max_step  << endl;
  pout() << "max_time= " << a_max_time  << endl;

#ifdef PROTO_CUDA
  Proto::DisjointBoxLayout::setNumStreams(nStream);
#endif

  RealVect ABC, X0;
  ABC[0] = A;
  ABC[1] = B;
#if DIM==3
  ABC[2] = C;
#endif
  Real  X0 = a_geomCen*RealVect::Unit();


  makeGrids(a_grids, a_dx, a_nx);
  Box domain = grids.physDomain().domainBox();
  int geomGhost = 4;
  RealVect origin = RealVect::Zero();

  pout() << "creating implicit function"
    shared_ptr<BaseIF>  impfunc(new SimpleEllipsoidIF(ABC, X0, a_geomRad, false));

  pout() << "creating geometry service"
    a_geoserv  = shared_ptr<GeometryService<MAX_ORDER> >(new GeometryService<MAX_ORDER>(impfunc, origin, a_dx, domain, a_grids, geomGhost));
}

int
runAdvection(int a_argc, char* a_argv[])
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
  int  max_step   = 10;
  Real max_time   = 1.0;
  int nStream    = 8;
  ParmParse pp;
    
  pp.get("nx"        , nx);
  pp.get("max_step"  , max_step);
  pp.get("max_time"  , max_time);
  pp.get("nstream", nStream);
  pp.get("maxGrid", maxGrid);
  pp.get("x0"     , x0);
  pp.get("y0"     , y0);
  pp.get("z0"     , z0);
  pp.get("A"      , A);
  pp.get("B"      , B);
  pp.get("C"      , C);
  pp.get("R"      , R);         
  pp.get("coveredval"      , coveredval);         

  pout() << "nx      = " << nx       << endl;
  pout() << "maxGrid = " << maxGrid  << endl;
  pout() << "x0      = " << x0       << endl;
  pout() << "y0      = " << y0       << endl;
  pout() << "z0      = " << z0       << endl;
  pout() << "A       = " << A        << endl;
  pout() << "B       = " << B        << endl;
  pout() << "C       = " << C        << endl;
  pout() << "R       = " << R        << endl;

  pout() << "max_step= " << max_step  << endl;
  pout() << "max_time= " << max_time  << endl;

#ifdef PROTO_CUDA
  Proto::DisjointBoxLayout::setNumStreams(nStream);
#endif

  RealVect ABC, X0;
  ABC[0] = A;
  ABC[1] = B;
  X0[0] = x0;
  X0[1] = y0;
#if DIM==3
  ABC[2] = C;
  X0[2] = z0;
#endif
  Real dx;
  DisjointBoxLayout& grids;

  pout() << "defining geometry" << endl;
  shared_ptr<GeometryService<MAX_ORDER> >  geoserv;

  Real geomCen;
  Real geomRad;
  Real blobCen;
  Real blobRad;
  defineGeometry(grids, dx, geomCen, geomRad, blobCen, blobRad, nx, max_step, max_time, geoserv);

  IntVect dataGhostIV =   IntVect::Unit;
  Point   dataGhostPt = getPoint(dataGhostIV); 

  pout() << "making dictionary" << endl;
  Box domain = grids.physDomain().domainBox();
  EBDictionary<2, Real, CELL, CELL> dictionary(geoserv, grids, domain, dataGhostPt, dataGhostPt, dx);


  pout() << "inititializing data"   << endl;
  
  EBLevelBoxData<CELL,   1>  scalcell(grids, dataGhostIV, graphs);
  EBLevelBoxData<CELL, DIM>  velocell(grids, dataGhostIV, graphs);
  initializeData(scalcell, velocell, grids, dx, , geomCen, geomRad, blobCen, blobRad);

  ParmParse pp;
  int outputInterval = -1;
  pp.get("output_interval", outputInterval);

  int step = 0; Real time = 0;
  Real dt = 0;
  pout() << "computing the time step"  << endl;
  computeDt(dt, velo, dx);

  if(outputInterval > 0)
  {
    string filev("velo.hdf5");
    velocell.writeToFileHDF5(filev, -1.);
    string filep("scal.0.hdf5");
    scalcell.writeToFileHDF5(filep, -1.);
  }

  pout() << "running advection operator " << endl;
  EBAdvectionOp advectOp(dictionary, geoserv, velocell, grids, domain,  dx, dataGhostIV, dataGhostIV);

  while((step < max_step) && (time < max_time))
  {
    advectOp.advance(dt, scalcell);

    pout() <<" step = " << step << " time = " << time << " time step = " << dt << endl;
    if((outputInterval > 0) && (k%outputInterval == 0))
    {
      string filep = string("scal.") + std::to_string(step) + string(".hdf5");
      scalcell.writeToFileHDF5(filep, -1.);
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
#include <cstdio>
#include <cstring>
#include <cassert>
#include <cmath>

#include <vector>
#include <memory>

#include <iostream>
#include <fstream>
#include <sstream>

#include "Proto.H"
#include "EulerRK4.H"
#include "Proto_WriteBoxData.H"
#include "Proto_Timer.H"


using namespace std;
using namespace Proto;




/***/
int main(int argc, char* argv[])
{
  //have to do this to get a time table
  PR_TIMER_SETFILE("proto.time.table");
  {
    viewDataNC(NULL);
    printDataNC(NULL, 0);

    PR_TIME("main");
    double tstop;
    int size1D, maxStep, outputInterval;
    parseCommandLine(tstop, size1D, maxStep, outputInterval, argc, argv);


    int nGhost = NGHOST;
    EulerOp::s_gamma = 1.4;
    EulerRK4Op::s_count = 0;
    Point lo = Point::Zeros();
    Point hi = Point::Ones(size1D - 1);
    Box dbx0(lo,hi);
    EulerOp::s_dx = 1./size1D;
    EulerState state(dbx0);
    RK4<EulerState,EulerRK4Op,EulerDX> rk4;
    Box dbx = dbx0.grow(nGhost);
    Box dbx1 = dbx.grow(1);
    BoxData<double,NUMCOMPS> UBig(dbx1);
    BoxData<double,DIM> x(dbx1);
    forallInPlace_p(iotaFunc, dbx1, x, EulerOp::s_dx);

    BoxData<double,NUMCOMPS>& U = state.m_U;
    //iota(x,EulerOp::s_dx);
    double dt = .25/size1D;
    Stencil<double> Lap2nd = Stencil<double>::Laplacian();
    cout << "before initializestate"<< endl;
    forallInPlace(InitializeState,dbx1,UBig,x);
    cout << "after initializestate"<< endl;

    U |= Lap2nd(UBig,dbx,1.0/24.0); 
    U += UBig;
    double time = 0.;
    string resStr = "_"+std::to_string(size1D);
    string fileRoot = "outfile";
    cout << "starting time loop"<< endl;
  }    
  PR_TIMER_REPORT();

}
