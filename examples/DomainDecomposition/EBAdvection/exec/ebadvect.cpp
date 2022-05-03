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

#define MAX_ORDER 2

using std::cout;
using std::endl;
using std::shared_ptr;
using Proto::Var;
using Proto::SimpleEllipsoidIF;

typedef Var<Real,DIM> Vec;
typedef Var<Real,  1> Sca;
/***/
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


//=================================================
void initializeData(EBLevelBoxData<CELL,   1>   &  a_scalcell,
                    EBLevelBoxData<CELL, DIM>   &  a_velocell,
                    const Chombo4::DisjointBoxLayout     &  a_grids,
                    const Real                  &  a_dx,
                    const Real                  &  a_geomCen,
                    const Real                  &  a_geomRad,
                    const Real                  &  a_blobCen,
                    const Real                  &  a_blobRad,
                    const Real                  &  a_maxVelMag,
                    const Real                  &  a_maxVelRad)
{
  using Chombo4::ProblemDomain;
  using Chombo4::DisjointBoxLayout;
  using Chombo4::LevelBoxData;
  using Chombo4::Copier;
  using Chombo4::DataIterator;
  using Chombo4::MayDay;

  DataIterator dit = a_grids.dataIterator();
  int ideb = 0;
  int whichvelo = 0;
  int whichscal = 0;
  ParmParse pp;
  pp.query("which_velo", whichvelo);
  if(whichvelo == 0)
  {
    Chombo4::pout() << "calling initializevel for velocity" << endl;
  }
  else if(whichvelo == -1)
  {
    Chombo4::pout() << "calling initializevelconst for velocity" << endl;
  }
  else
  {
    Chombo4::MayDay::Error("bogus velo");
  }

  pp.query("which_scal", whichscal);
  if(whichscal == 0)
  {
    Chombo4::pout() << "calling initializespot for scalar " << endl;
  }
  else if(whichscal == -1)
  {
    Chombo4::pout() << "calling initializeline for scalar" << endl;
  }
  else
  {
    Chombo4::MayDay::Error("bogus scal");
  }

  for(int ibox = 0; ibox < dit.size(); ibox++)
  {
    {
      auto& scalfab = a_scalcell[dit[ibox]];
      unsigned long long int numflopsscal = 5*DIM +3;
      Bx scalbox = scalfab.box();

      if(whichscal == 0)
      {
        ebforallInPlace_i(numflopsscal, "IntializeSpot", InitializeSpot,  scalbox, 
                          scalfab, a_blobCen, a_blobRad, a_dx);
        ideb++;
      }
      else
      {
        vector<double> v_norm;
        pp.getarr("geom_normal", v_norm, 0, DIM);
        ebforallInPlace_i(numflopsscal, "IntializeLine", InitializeLine,  scalbox, 
                          scalfab, a_blobCen, a_blobRad, a_dx);
        ideb++;
      }
    }

    {
    
      auto& velofab = a_velocell[dit[ibox]];
      unsigned long long int numflopsvelo = (DIM+5)*DIM +4;
      Bx velobox = velofab.box();
      if(whichvelo == 0)
      {
        bool solidBody = false;
        pp.query("solid_body_rot", solidBody);
        if(solidBody)
        {
          Chombo4::pout() << "solid body rotation ON" << endl;
        }
        else
        {
          Chombo4::pout() << "solid body rotation OFF" << endl;
        }
        ebforallInPlace_i(numflopsvelo, "IntializeVCell", InitializeVCell,  velobox, 
                          velofab, a_geomCen, a_geomRad, a_maxVelMag, a_maxVelRad, 
                          solidBody, a_dx);
      }
      else if(whichvelo == -1)
      {
        vector<double> v_norm;
        pp.getarr("geom_normal", v_norm, 0, DIM);
        Real xvel =  v_norm[1];
        Real yvel = -v_norm[0];
        ebforallInPlace_i(numflopsvelo, "IntializeVCellConst", InitializeVCellConst,  velobox, 
                          velofab, xvel, yvel);
        ideb++;
      }
      
      ideb++;
    }
  }
}
//=================================================

void makeGrids(Chombo4::DisjointBoxLayout& a_grids,
               Real             & a_dx,
               const int        & a_nx,
               const int        & a_maxGrid)
{
  Chombo4::pout() << "making grids" << endl;

  using Chombo4::ProblemDomain;
  using Chombo4::DisjointBoxLayout;
  using Chombo4::LevelBoxData;
  using Chombo4::Copier;
  using Chombo4::DataIterator;
  using Chombo4::MayDay;
  
  IntVect domLo = IntVect::Zero;
  IntVect domHi  = (a_nx - 1)*IntVect::Unit;

  // EB and periodic do not mix
  ProblemDomain domain(domLo, domHi);

  std::vector<Chombo4::Box> boxes;
  unsigned int blockfactor = 8;
  domainSplit(domain.domainBox(), boxes, a_maxGrid, blockfactor);
  
  std::vector<int> procs;

  a_dx = 1.0/a_nx;
  LoadBalance(procs, boxes);
  a_grids = Chombo4::DisjointBoxLayout(boxes, procs, domain);
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
    a_dt = a_cfl*a_dx/maxvel;
    Chombo4::pout() << "maxvel = " << maxvel << ", dx = " << a_dx << ", dt = " << a_dt << endl;
  }    
  else
  {
    Chombo4::pout() << "velocity seems to be zero--setting dt to dx" << endl;
    a_dt = a_dx;
  }
}
///
shared_ptr<BaseIF>  getImplicitFunction(Real  & a_geomCen,
                                        Real  & a_geomRad,
                                        int   & a_whichGeom)

{
  using Chombo4::ProblemDomain;
  using Chombo4::DisjointBoxLayout;
  using Chombo4::LevelBoxData;
  using Chombo4::Copier;
  using Chombo4::DataIterator;
  using Chombo4::MayDay;
  using Proto::BaseIF;
  
  shared_ptr<BaseIF>  retval;
  ParmParse pp;
  
  a_geomCen = 0;
  a_geomRad = 1;
  pp.get("which_geom", a_whichGeom);
  if(a_whichGeom == -1)
  {
    using Proto::AllRegularIF;
    Chombo4::pout() << "all regular geometry" << endl;
    retval = shared_ptr<BaseIF>(new AllRegularIF());
  }
  else if(a_whichGeom == 0)
  {
    using Proto::SimpleEllipsoidIF;
    Chombo4::pout() << "sphere" << endl;

    pp.get("geom_cen", a_geomCen);
    pp.get("geom_rad", a_geomRad);
    Chombo4::pout() << "geom_cen = " << a_geomCen       << endl;
    Chombo4::pout() << "geom_rad = " << a_geomRad       << endl;

    RealVect ABC = RealVect::Unit(); //this is what it makes it a sphere instead of an ellipse
    RealVect  X0 = RealVect::Unit();
    X0 *= a_geomCen;

    retval = shared_ptr<BaseIF>(new SimpleEllipsoidIF(ABC, X0, a_geomRad, true));//true is for inside regular
  }
  else if(a_whichGeom ==  1)
  {
    using Proto::PlaneIF;
    Chombo4::pout() << "plane" << endl;
    RealVect normal, startPt;
    vector<double> v_norm, v_start;
    pp.getarr("geom_normal", v_norm, 0, DIM);
    pp.getarr("geom_start_pt", v_start, 0, DIM);
    for(int idir = 0; idir < DIM; idir++)
    {
      normal[ idir] = v_norm[ idir];
      startPt[idir] = v_start[idir];
      Chombo4::pout() << "normal ["<< idir << "] = " << normal [idir]  << endl;
      Chombo4::pout() << "startPt["<< idir << "] = " << startPt[idir]  << endl;
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
void defineGeometry(Chombo4::DisjointBoxLayout& a_grids,
                    Real             & a_dx,
                    Real             & a_geomCen,
                    Real             & a_geomRad,
                    Real             & a_blobCen,
                    Real             & a_blobRad,
                    int              & a_whichGeom,
                    int              & a_nx,
                    shared_ptr<GeometryService<MAX_ORDER> >&  a_geoserv)
{
  Chombo4::pout() << "defining geometry" << endl;

  ParmParse pp;
  int maxGrid = 32;
    
  pp.get("nx"        , a_nx);
  pp.get("maxGrid", maxGrid);
  pp.get("blob_cen", a_blobCen);
  pp.get("blob_rad", a_blobRad);

  Chombo4::pout() << "nx       = " << a_nx     << endl;
  Chombo4::pout() << "maxGrid  = " << maxGrid  << endl;
  Chombo4::pout() << "blob_cen = " << a_blobCen       << endl;
  Chombo4::pout() << "blob_rad = " << a_blobRad       << endl;

  makeGrids(a_grids, a_dx, a_nx, maxGrid);
  Chombo4::Box domain = a_grids.physDomain().domainBox();
  int geomGhost = 4;
  RealVect origin = RealVect::Zero();

  Chombo4::pout() << "creating implicit function" << endl;
  shared_ptr<BaseIF>  impfunc = getImplicitFunction(a_geomCen, a_geomRad, a_whichGeom);

  Chombo4::pout() << "creating geometry service" << endl;
  a_geoserv  = shared_ptr<GeometryService<MAX_ORDER> >(new GeometryService<MAX_ORDER>(impfunc, origin, a_dx, domain, a_grids, geomGhost));
}

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
  EBLevelBoxData<CELL,   1>  scalcell(grids, dataGhostIV, graphs);

  shared_ptr<EBLevelBoxData<CELL, DIM> >  velocell(new EBLevelBoxData<CELL, DIM>(grids, dataGhostIV, graphs));
  velocell->setVal(0.);
  initializeData(scalcell, *velocell, grids, dx, geomCen, geomRad, blobCen, blobRad, max_vel_mag, max_vel_rad);


  int step = 0; Real time = 0;
  Real dt = 0;
  Chombo4::pout() << "computing the time step"  << endl;
  computeDt(dt, *velocell, dx, cfl);
  EBIBC bc = getIBCs();
  EBAdvection advectOp(brit, geoserv, velocell, grids, domain,  dx, bc, dataGhostIV);
  const EBLevelBoxData<CELL, 1> & kappa = advectOp.getKappa();

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


