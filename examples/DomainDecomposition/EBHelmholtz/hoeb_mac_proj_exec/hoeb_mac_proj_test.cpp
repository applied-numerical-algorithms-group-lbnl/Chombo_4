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
#include "Chombo_EBLevelFluxData.H"
#include "Hoeb_MAC_Projector.H"
#include "SetupFunctions.H"

#include <iomanip>

#define MAX_ORDER 4

using std::cout;
using std::endl;
using std::shared_ptr;
using Proto::Var;
using Proto::SimpleEllipsoidIF;

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


//=================================================
void initializeData(EBLevelFluxData<1>          &  a_velo,
                    const Chombo4::DisjointBoxLayout     &  a_grids,
                    const Real                  &  a_dx,
                    const Real                  &  a_geomCen,
                    const Real                  &  a_geomRad,
                    const Real                  &  a_blobCen,
                    const Real                  &  a_blobRad,
                    const Real                  &  a_maxVelMag,
                    const Real                  &  a_maxVelRad)
{
  Chombo4::DataIterator dit = a_grids.dataIterator();
  int ideb = 0;

  for(int ibox = 0; ibox < dit.size(); ibox++)
  {
    auto& velofab = a_velo[dit[ibox]];
    unsigned long long int numflopsvelo = (DIM+5)*DIM +4;

    {
      auto& xvel = *velofab.m_xflux;
      unsigned int facedir = 0;
      ebforallInPlace_i(numflopsvelo, "InitializeVFace", InitializeVFace,  xvel.box(), xvel, 
                        a_geomCen, a_geomRad, a_blobCen, a_blobRad, a_maxVelMag, a_maxVelRad, a_dx, facedir);
    }
    {
      auto& yvel = *velofab.m_yflux;
      unsigned int facedir = 1;
      ebforallInPlace_i(numflopsvelo, "InitializeVFace", InitializeVFace,  yvel.box(), yvel, 
                        a_geomCen, a_geomRad, a_blobCen, a_blobRad, a_maxVelMag, a_maxVelRad, a_dx, facedir);
    }

#if DIM==3
    {
      auto& zvel = *velofab.m_zflux;
      unsigned int facedir = 2;
      ebforallInPlace_i(numflopsvelo, "InitializeVFace", InitializeVFace,  zvel.box(), zvel, 
                        a_geomCen, a_geomRad, a_blobCen, a_blobRad, a_maxVelMag, a_maxVelRad, a_dx, facedir);
    }
#endif
    ideb++;
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
    Chombo4::MayDay::Error("bogus geometry");
  }
  return retval;
}
//=================================================
void defineGeometry(std::vector<Chombo4::DisjointBoxLayout>& a_grids,
                    const Chombo4::Box        & a_finestDomain,
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
    
  pp.get("blob_cen", a_blobCen);
  pp.get("blob_rad", a_blobRad);

  Chombo4::pout() << "blob_cen = " << a_blobCen       << endl;
  Chombo4::pout() << "blob_rad = " << a_blobRad       << endl;

  int geomGhost = 4;
  RealVect origin = RealVect::Zero();

  Chombo4::pout() << "creating implicit function" << endl;
  shared_ptr<BaseIF>  impfunc = getImplicitFunction(a_geomCen, a_geomRad, a_whichGeom);

  Chombo4::pout() << "creating geometry service" << endl;
  Chombo4::Box domain = a_finestDomain;
  a_geoserv  = shared_ptr<GeometryService<MAX_ORDER> >(new GeometryService<MAX_ORDER>(impfunc, origin, a_dx, domain, a_grids, geomGhost));
}

int
runProjection(int a_argc, char* a_argv[])
{
  int nx      = 32;
  Real max_vel_mag = 1.0;
  Real max_vel_rad = 0.25;
  int numSmooth;
  ParmParse pp;
  int maxGrid = 32;
  pp.get("maxGrid", maxGrid); 

  pp.get("numSmooth"  , numSmooth);
  pp.get("max_vel_mag"  , max_vel_mag);
  pp.get("max_vel_rad"  , max_vel_rad);
  pp.get("nx"           , nx);

  Chombo4::pout() << "nx       = " << nx     << endl;
  Chombo4::pout() << "maxGrid  = " << maxGrid  << endl;
  Chombo4::pout() << "max_vel_mag     = " << max_vel_mag     << endl;
  Chombo4::pout() << "max_vel_rad     = " << max_vel_rad     << endl;

  Real dx = 1.0/nx;
  std::vector<Chombo4::DisjointBoxLayout> vecgrids;

  Chombo4::pout() << "making grids" << endl;
  IntVect domLo = IntVect::Zero;
  IntVect domHi  = (nx - 1)*IntVect::Unit;

// EB and periodic do not mix
  Chombo4::ProblemDomain domain(domLo, domHi);
  GeometryService<2>::generateGrids(vecgrids, domain.domainBox(), maxGrid);

  Chombo4::pout() << "defining geometry" << endl;
  shared_ptr<GeometryService<MAX_ORDER> >  geoserv;

  Real geomCen;
  Real geomRad;
  Real blobCen;
  Real blobRad;
  int whichGeom;

  Chombo4::Box domainb = domain.domainBox();
  defineGeometry(vecgrids, domainb, dx, geomCen, geomRad, blobCen, blobRad, whichGeom, nx,  geoserv);

  IntVect dataGhostIV =   4*IntVect::Unit;
  Point   dataGhostPt = ProtoCh::getPoint(dataGhostIV); 

  Chombo4::pout() << "making dictionary" << endl;
  vector<Chombo4::Box>    vecdomain(vecgrids.size(), domain.domainBox());
  vector<Real>   vecdx    (vecgrids.size(), dx);
  for(int ilev = 1; ilev < vecgrids.size(); ilev++)
  {
    vecdomain[ilev] = coarsen(vecdomain[ilev-1], 2);
    vecdx    [ilev] =           2*vecdx[ilev-1];
  }

  shared_ptr<EBEncyclopedia<2, Real> > 
    brit(new EBEncyclopedia<2, Real>(geoserv, vecgrids, vecdomain, vecdx, dataGhostPt));


  Chombo4::pout() << "inititializing data"   << endl;
  
  auto graphs = geoserv->getGraphs(domain.domainBox());
  Chombo4::DisjointBoxLayout& grids = vecgrids[0];
  shared_ptr< EBLevelFluxData<1>  > velo( new EBLevelFluxData<1>(grids, dataGhostIV, graphs));
  shared_ptr< EBLevelFluxData<1>  > gphi( new EBLevelFluxData<1>(grids, dataGhostIV, graphs));

  initializeData(velo, grids, dx, geomCen, geomRad, blobCen, blobRad, max_vel_mag, max_vel_rad);


  EBIBC bc = getIBCs();
  Hoeb_MAC_Projector<MAX_ORDER> proj(brit, geoserv, grids, domain.domainBox(), dx, dataGhostIV, bc);
  Real tol = 1.0e-8;
  unsigned int maxiter = 27;
  proj.project(velo, gphi, tol, maxiter, true);

  EBLevelBoxData<CELL,   1>&  divergence = proj.getRHSHolder();
  proj.kappaDivU(divergence, velo, true);
  
  Real divnorm = divergence.maxNorm(0);
  
  Chombo4::pout() << "max norm of post projection divergence(vel) = " << divnorm << endl;
   
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
    runProjection(a_argc, a_argv);
  }

  Chombo4::pout() << "printing time table " << endl;
  CH_TIMER_REPORT();
#ifdef CH_MPI
  Chombo4::pout() << "about to call MPI Finalize" << std::endl;
  MPI_Finalize();
#endif
  return 0;
}


