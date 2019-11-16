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
#if DIM==2
void 
dumpBlob(BoxData<Real, 1>* dataPtr)
{
  if(dataPtr != NULL)
  {
    cout    << setprecision(6)
            << setiosflags(ios::showpoint)
            << setiosflags(ios::scientific);
    Point lo(6,  7);
    Point hi(10, 9);
    Bx area(lo, hi);

    BoxData<Real, 1> & data = *dataPtr;
    Bx databox = dataPtr->box();
    cout << "data region contains:" << endl;
    for(int j = hi[1]; j >= lo[1]; j--)
    {
      for(int i = lo[0]; i <= hi[0]; i++)
      {
        Point pt(i,j);
        if(databox.contains(pt))
        {
          cout << pt << ":" << data(pt, 0) << "  ";
        }
      }
      cout << endl;
    }
  }
}
#endif

//=================================================
void initializeData(EBLevelBoxData<CELL,   1>   &  a_scalcell,
                    EBLevelBoxData<CELL, DIM>   &  a_velocell,
                    const DisjointBoxLayout     &  a_grids,
                    const Real                  &  a_dx,
                    const Real                  &  a_geomCen,
                    const Real                  &  a_geomRad,
                    const Real                  &  a_blobCen,
                    const Real                  &  a_blobRad,
                    const Real                  &  a_maxVelMag,
                    const Real                  &  a_maxVelRad)
{
  DataIterator dit = a_grids.dataIterator();
  int ideb = 0;
  for(int ibox = 0; ibox < dit.size(); ibox++)
  {
    {
      auto& scalfab = a_scalcell[dit[ibox]];
      unsigned long long int numflopsscal = 5*DIM +3;
      Bx scalbox = scalfab.box();

      ebforallInPlace_i(numflopsscal, "IntializeSpot", InitializeSpot,  scalbox, 
                        scalfab, a_blobCen, a_blobRad, a_dx);
      ideb++;
    }

    {
      auto& velofab = a_velocell[dit[ibox]];
      unsigned long long int numflopsvelo = (DIM+5)*DIM +4;
      Bx velobox = velofab.box();
      ebforallInPlace_i(numflopsvelo, "IntializeVCell", InitializeVCell,  velobox, 
                        velofab, a_geomCen, a_geomRad, a_maxVelMag, a_maxVelRad, a_dx);

    
      ideb++;
    }
  }
}
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
               const   Real                &  a_dx,
               const   Real                &  a_cfl)
{
  Real maxvel = 0;
  for(int idir = 0; idir < DIM; idir++)
  {
    maxvel = std::max(maxvel, a_velocell.maxNorm(idir));
  }
  if(maxvel > 1.0e-16)
  {
    pout() << "maxvel = " << maxvel << ", dx = " << a_dx << ", dt = " << a_dt << endl;
    a_dt = a_cfl*a_dx/maxvel;
  }    
  else
  {
    pout() << "velocity seems to be zero--setting dt to dx" << endl;
    a_dt = a_dx;
  }
}
///
shared_ptr<BaseIF>  getImplicitFunction(Real  & a_geomCen,
                                        Real  & a_geomRad,
                                        int   & a_whichGeom)

{
  using Proto::BaseIF;
  shared_ptr<BaseIF>  retval;
  ParmParse pp;
  
  a_geomCen = 0;
  a_geomRad = 1;
  pp.get("which_geom", a_whichGeom);
  if(a_whichGeom == -1)
  {
    using Proto::AllRegularIF;
    pout() << "all regular geometry" << endl;
    retval = shared_ptr<BaseIF>(new AllRegularIF());
  }
  else if(a_whichGeom == 0)
  {
    using Proto::SimpleEllipsoidIF;
    pout() << "sphere" << endl;

    pp.get("geom_cen", a_geomCen);
    pp.get("geom_rad", a_geomRad);
    pout() << "geom_cen = " << a_geomCen       << endl;
    pout() << "geom_rad = " << a_geomRad       << endl;

    RealVect ABC = RealVect::Unit(); //this is what it makes it a sphere instead of an ellipse
    RealVect  X0 = RealVect::Unit();
    X0 *= a_geomCen;

    retval = shared_ptr<BaseIF>(new SimpleEllipsoidIF(ABC, X0, a_geomRad, true));//true is for inside regular
  }
  else if(a_whichGeom ==  1)
  {
    using Proto::PlaneIF;
    pout() << "plane" << endl;
    RealVect normal, startPt;
    vector<double> v_norm, v_start;
    pp.getarr("geom_normal", v_norm, 0, DIM);
    pp.getarr("geom_start_pt", v_start, 0, DIM);
    for(int idir = 0; idir < DIM; idir++)
    {
      normal[ idir] = v_norm[ idir];
      startPt[idir] = v_start[idir];
      pout() << "normal ["<< idir << "] = " << normal [idir]  << endl;
      pout() << "startPt["<< idir << "] = " << startPt[idir]  << endl;
    }
    retval = shared_ptr<BaseIF>(new PlaneIF(startPt, normal));
  }
  else
  {
    MayDay::Error("bogus geometry");
  }
  return retval;
}
//=================================================
void defineGeometry(DisjointBoxLayout& a_grids,
                    Real             & a_dx,
                    Real             & a_geomCen,
                    Real             & a_geomRad,
                    Real             & a_blobCen,
                    Real             & a_blobRad,
                    int              & a_whichGeom,
                    int              & a_nx,
                    shared_ptr<GeometryService<MAX_ORDER> >&  a_geoserv)
{
  pout() << "defining geometry" << endl;

  ParmParse pp;
  int maxGrid = 32;
    
  pp.get("nx"        , a_nx);
  pp.get("maxGrid", maxGrid);
  pp.get("blob_cen", a_blobCen);
  pp.get("blob_rad", a_blobRad);

  pout() << "nx       = " << a_nx     << endl;
  pout() << "maxGrid  = " << maxGrid  << endl;
  pout() << "blob_cen = " << a_blobCen       << endl;
  pout() << "blob_rad = " << a_blobRad       << endl;

  makeGrids(a_grids, a_dx, a_nx, maxGrid);
  Box domain = a_grids.physDomain().domainBox();
  int geomGhost = 4;
  RealVect origin = RealVect::Zero();

  pout() << "creating implicit function" << endl;
  shared_ptr<BaseIF>  impfunc = getImplicitFunction(a_geomCen, a_geomRad, a_whichGeom);

  pout() << "creating geometry service" << endl;
  a_geoserv  = shared_ptr<GeometryService<MAX_ORDER> >(new GeometryService<MAX_ORDER>(impfunc, origin, a_dx, domain, a_grids, geomGhost));
}

int
runAdvection(int a_argc, char* a_argv[])
{
#if DIM==2
  dumpBlob(NULL);
#endif

  Real coveredval = -1;
  Real cfl    = 0.5;
  int nx      = 32;
  int  max_step   = 10;
  Real max_time   = 1.0;
  Real max_vel_mag = 1.0;
  Real max_vel_rad = 0.25;
  int nStream    = 8;
  int outputInterval = -1;
  ParmParse pp;

  pp.get("nstream", nStream);

    
  pp.get("max_step"  , max_step);
  pp.get("max_time"  , max_time);
  pp.get("output_interval", outputInterval);
  pp.get("covered_value", coveredval);
  pp.get("cfl"  , cfl);
  pp.get("max_vel_mag"  , max_vel_mag);
  pp.get("max_vel_rad"  , max_vel_rad);

  pout() << "num_streams     = " << nStream         << endl;
  pout() << "max_step        = " << max_step        << endl;
  pout() << "max_time        = " << max_time        << endl;
  pout() << "output interval = " << outputInterval  << endl;
  pout() << "cfl             = " << cfl             << endl;
  pout() << "max_vel_mag     = " << max_vel_mag     << endl;
  pout() << "max_vel_rad     = " << max_vel_rad     << endl;

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
    brit(new EBEncyclopedia<2, Real>(geoserv, grids, domain, dx, dataGhostPt, dataGhostPt));


  pout() << "inititializing data"   << endl;
  
  shared_ptr<LevelData<EBGraph> > graphs = geoserv->getGraphs(domain);
  EBLevelBoxData<CELL,   1>  scalcell(grids, dataGhostIV, graphs);

  shared_ptr<EBLevelBoxData<CELL, DIM> >  velocell(new EBLevelBoxData<CELL, DIM>(grids, dataGhostIV, graphs));
  initializeData(scalcell, *velocell, grids, dx, geomCen, geomRad, blobCen, blobRad, max_vel_mag, max_vel_rad);


  int step = 0; Real time = 0;
  Real dt = 0;
  pout() << "computing the time step"  << endl;
  computeDt(dt, *velocell, dx, cfl);

  if(outputInterval > 0)
  {
    string filev("velo.hdf5");
    velocell->writeToFileHDF5(filev, coveredval);
    string filep("scal.0.hdf5");
    scalcell.writeToFileHDF5(filep, coveredval);

//    DataIterator dit = grids.dataIterator();
//    for(int ibox = 0; ibox < dit.size(); ibox++)
//    {
//      auto& scalfab = scalcell[dit[ibox]];
//      dumpBlob(&scalfab.getRegData());
//    }
  }

  pout() << "running advection operator " << endl;
  EBAdvection advectOp(brit, geoserv, velocell, grids, domain,  dx, dataGhostIV, dataGhostIV);

  while((step < max_step) && (time < max_time))
  {
    advectOp.advance(scalcell, dt);

    pout() <<" step = " << step << " time = " << time << " time step = " << dt << endl;
    step++;
    time += dt;

    if((outputInterval > 0) && ( (step%outputInterval == 0) || step == (max_step-1)))
    {
      string filep = string("scal.") + std::to_string(step) + string(".hdf5");
      scalcell.writeToFileHDF5(filep, coveredval);
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


