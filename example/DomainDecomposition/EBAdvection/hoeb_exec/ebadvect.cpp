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
#include "EBAdvection.H"
#include "SetupFunctions.H"

#include <iomanip>

#define MAX_ORDER 4
namespace Chombo4
{
  EBIBC getIBCs()
  {
    string veloIC("does_not_matter");
    string scalIC("does_not_matter");
    string loDomBC[DIM];
    string hiDomBC[DIM];
    ParmParse pp;
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


///
  void initializeData(EBLevelBoxData<CELL, 1>     &  a_scalcell,
                      EBLevelFluxData<1>          &  a_veloface,
                      const Hoeb_MAC_Projector    &  a_projectionOperator)
  {
    a_scalcell.setVal(1.);
    a_veloface.setVal(1.);
    a_projectionOperator.project(a_veloface);
  }
///

  void makeGrids(Chombo4::DisjointBoxLayout& a_grids,
                 Real                      & a_dx,
                 const int                 & a_nx,
                 const int                 & a_maxGrid)
  {
    Chombo4::pout() << "making grids" << endl;

    IntVect domLo = IntVect::Zero;
    IntVect domHi  = (a_nx - 1)*IntVect::Unit;

    // EB and periodic do not mix
    ProblemDomain domain(domLo, domHi);

    std::vector<Box> boxes(1, domain.domainBox());
    std::vector<int> procs(1, 0);

    a_dx = 1.0/a_nx;
    a_grids = DisjointBoxLayout(boxes, procs, domain);
    a_grids.printBalance();
  }

///
  void computeDt(Real                        &  a_dt,
                 EBLevelFluxData<1>          &  a_veloface,
                 const   Real                &  a_dx,
                 const   Real                &  a_cfl)
  {
    Real maxvel = 1.;
    a_dt = a_cfl*a_dx/maxvel;
    Chombo4::pout() << "maxvel = " << maxvel << ", dx = " << a_dx << ", dt = " << a_dt << endl;
  }
///
  shared_ptr<BaseIF>  getImplicitFunction()
  {
    shared_ptr<BaseIF>  retval;
    ParmParse pp;
  
    using Proto::AllRegularIF;
    Chombo4::pout() << "all regular geometry" << endl;
    retval = shared_ptr<BaseIF>(new AllRegularIF());
    return retval;
  }
///
  void defineGeometry(DisjointBoxLayout                      & a_grids,
                      Real                                   & a_dx,
                      int                                    & a_nx,
                      shared_ptr<GeometryService<MAX_ORDER> >&  a_geoserv)
  {
    Chombo4::pout() << "defining geometry" << endl;
    ParmParse pp;
    int maxGrid = 32;
    pp.get("nx"        , a_nx);
    pp.get("maxGrid", maxGrid);
    Chombo4::pout() << "nx       = " << a_nx     << endl;
    Chombo4::pout() << "maxGrid  = " << maxGrid  << endl;
    makeGrids(a_grids, a_dx, a_nx, maxGrid);
    Chombo4::Box domain = a_grids.physDomain().domainBox();

    int geomGhost = 4;
    RealVect origin = RealVect::Zero();

    Chombo4::pout() << "creating implicit function" << endl;
    shared_ptr<BaseIF>  impfunc = getImplicitFunction();

    Chombo4::pout() << "creating geometry service" << endl;
    a_geoserv  = shared_ptr<GeometryService<MAX_ORDER> >(new GeometryService<MAX_ORDER>(impfunc, origin, a_dx, domain, a_grids, geomGhost));
  }
} //namespace Chombo4
int
runAdvection(int a_argc, char* a_argv[])
{

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

  Chombo4::pout() << "num_streams     = " << nStream         << endl;
  Chombo4::pout() << "max_step        = " << max_step        << endl;
  Chombo4::pout() << "max_time        = " << max_time        << endl;
  Chombo4::pout() << "output interval = " << outputInterval  << endl;
  Chombo4::pout() << "cfl             = " << cfl             << endl;
  Chombo4::pout() << "max_vel_mag     = " << max_vel_mag     << endl;
  Chombo4::pout() << "max_vel_rad     = " << max_vel_rad     << endl;

  Real dx;
  Chombo4::DisjointBoxLayout grids;

  Chombo4::pout() << "defining geometry" << endl;
  shared_ptr<GeometryService<MAX_ORDER> >  geoserv;

  Real geomCen;
  Real geomRad;
  Real blobCen;
  Real blobRad;
  int whichGeom;
  defineGeometry(grids, dx, geomCen, geomRad, blobCen, blobRad, whichGeom, nx,  geoserv);

  IntVect dataGhostIV =   4*IntVect::Unit;
  Point   dataGhostPt = ProtoCh::getPoint(dataGhostIV); 

  Chombo4::pout() << "making dictionary" << endl;
  Chombo4::Box domain = grids.physDomain().domainBox();
  shared_ptr<EBEncyclopedia<2, Real> > 
    brit(new EBEncyclopedia<2, Real>(geoserv, grids, domain, dx, dataGhostPt));

  Chombo4::pout() << "inititializing data"   << endl;
  
  auto graphs = geoserv->getGraphs(domain);
  shared_ptr<EBLevelBoxData<CELL,   1>  > scalcell( new EBLevelBoxData<CELL, 1>(grids, dataGhostIV, graphs));
 
  shared_ptr<EBLevelFluxData<1> >  advectiveVel(new EBLevelFluxData<1>(grids, dataGhostIV, graphs));
  initializeData(*scalcell, *velocell);


  int step = 0; Real time = 0;
  Real dt = 0;
  Chombo4::pout() << "computing the time step"  << endl;
  computeDt(dt, *velocell, dx, cfl);
  EBIBC bc = getIBCs();
  Hoeb_Advection advectOp(brit, geoserv, grids, domain,  dx, bc, dataGhostIV);
  const EBLevelBoxData<CELL, 1> & kappa = *advectOp.m_kappa;

  if(outputInterval > 0)
  {
    string filev("velo.hdf5");
    string filep("scal.0.hdf5");
    writeEBLevelHDF5<DIM>(filev, *velocell, kappa, domain, graphs, coveredval, dx, dt, time);
    writeEBLevelHDF5<1>(  filep,  scalcell, kappa, domain, graphs, coveredval, dx, dt, time);
  }

  Chombo4::pout() << "running advection operator " << endl;

  while((step < max_step) && (time < max_time))
  {
    Real fluxval;
    ParmParse pp;
    pp.get("scalar_inflow_value",   fluxval);
    advectOp.advance(scalcell,  dt, fluxval);

    Chombo4::pout() <<" step = " << step << " time = " << time << " time step = " << dt << endl;
    step++;
    time += dt;

    if((outputInterval > 0) && ( (step%outputInterval == 0) || step == (max_step-1)))
    {
      string filep = string("scal.") + std::to_string(step) + string(".hdf5");
      writeEBLevelHDF5<1>(  filep,  scalcell, kappa, domain, graphs, coveredval, dx, dt, time);
    }
  }
  Chombo4::pout() << "exiting runAdvection" << endl;
  return 0;
}


int main(int a_argc, char* a_argv[])
{
#ifdef CH_MPI
  MPI_Init(&a_argc, &a_argv);
  Chombo4::pout() << "MPI INIT called" << std::endl;
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

  Chombo4::pout() << "printing time table " << endl;
  CH_TIMER_REPORT();
#ifdef CH_MPI
  Chombo4::pout()  << "about to call MPI Finalize" << std::endl;
  MPI_Finalize();
#endif
  return 0;
}


