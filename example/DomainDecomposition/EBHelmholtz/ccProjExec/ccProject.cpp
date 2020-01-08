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
#include "EBCCProjector.H"
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
void initializeData(EBLevelBoxData<CELL, DIM>   &  a_velo,
                    const DisjointBoxLayout     &  a_grids,
                    const Real                  &  a_dx,
                    const Real                  &  a_geomCen,
                    const Real                  &  a_geomRad,
                    const Real                  &  a_maxVelMag,
                    const Real                  &  a_maxVelRad)
{
  DataIterator dit = a_grids.dataIterator();
  int ideb = 0;

  for(int ibox = 0; ibox < dit.size(); ibox++)
  {
    auto& velofab = a_velo[dit[ibox]];
    unsigned long long int numflopsvelo = (DIM+5)*DIM +4;

    ebforallInPlace_i(numflopsvelo, "InitializeVCell", InitializeVCell,  velofab.box(), velofab, 
                      a_geomCen, a_geomRad, a_maxVelMag, a_maxVelRad, a_dx);
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
void defineGeometry(Vector<DisjointBoxLayout>& a_grids,
                    const Box        & a_finestDomain,
                    Real             & a_dx,
                    Real             & a_geomCen,
                    Real             & a_geomRad,
                    int              & a_whichGeom,
                    int              & a_nx,
                    shared_ptr<GeometryService<MAX_ORDER> >&  a_geoserv)
{
  pout() << "defining geometry" << endl;

  ParmParse pp;
    
  int geomGhost = 4;
  RealVect origin = RealVect::Zero();

  pout() << "creating implicit function" << endl;
  shared_ptr<BaseIF>  impfunc = getImplicitFunction(a_geomCen, a_geomRad, a_whichGeom);

  pout() << "creating geometry service" << endl;
  Box domain = a_finestDomain;
  a_geoserv  = shared_ptr<GeometryService<MAX_ORDER> >(new GeometryService<MAX_ORDER>(impfunc, origin, a_dx, domain, a_grids, geomGhost));
}

int
runProjection(int a_argc, char* a_argv[])
{
#if DIM==2
  //debugging hooks to get these functions into the symbol table
  dumpBlob(NULL);
  dumpEB1(NULL);
  dumpEBDIM(NULL);
  dumpXFace(NULL);
  dumpYFace(NULL);
#endif

  int nx      = 32;
  Real max_vel_mag = 1.0;
  Real max_vel_rad = 0.25;
  int nStream    = 8;
  int numSmooth;
  ParmParse pp;
  int maxGrid = 32;
  pp.get("maxGrid", maxGrid); 

  pp.get("nstream", nStream);

    
  pp.get("numSmooth"  , numSmooth);
  pp.get("max_vel_mag"  , max_vel_mag);
  pp.get("max_vel_rad"  , max_vel_rad);
  pp.get("nx"           , nx);

  pout() << "nx       = " << nx     << endl;
  pout() << "maxGrid  = " << maxGrid  << endl;
  pout() << "max_vel_mag     = " << max_vel_mag     << endl;
  pout() << "max_vel_rad     = " << max_vel_rad     << endl;
  pout() << "num_streams     = " << nStream         << endl;

  EBMultigrid::s_numSmoothUp   = numSmooth;
  EBMultigrid::s_numSmoothDown = numSmooth;
#ifdef PROTO_CUDA
  Proto::DisjointBoxLayout::setNumStreams(nStream);
#endif

  Real dx = 1.0/nx;
  Vector<DisjointBoxLayout> vecgrids;

  pout() << "making grids" << endl;
  IntVect domLo = IntVect::Zero;
  IntVect domHi  = (nx - 1)*IntVect::Unit;

// EB and periodic do not mix
  ProblemDomain domain(domLo, domHi);
  GeometryService<2>::generateGrids(vecgrids, domain.domainBox(), maxGrid);

  pout() << "defining geometry" << endl;
  shared_ptr<GeometryService<MAX_ORDER> >  geoserv;

  Real geomCen;
  Real geomRad;

  int whichGeom;

  Box domainb = domain.domainBox();
  defineGeometry(vecgrids, domainb, dx, geomCen, geomRad, whichGeom, nx,  geoserv);

  IntVect dataGhostIV =   4*IntVect::Unit;
  Point   dataGhostPt = ProtoCh::getPoint(dataGhostIV); 

  pout() << "making dictionary" << endl;
  vector<Box>    vecdomain(vecgrids.size(), domain.domainBox());
  vector<Real>   vecdx    (vecgrids.size(), dx);
  for(int ilev = 1; ilev < vecgrids.size(); ilev++)
  {
    vecdomain[ilev] = coarsen(vecdomain[ilev-1], 2);
    vecdx    [ilev] =           2*vecdx[ilev-1];
  }

  shared_ptr<EBEncyclopedia<2, Real> > 
    brit(new EBEncyclopedia<2, Real>(geoserv, vecgrids, vecdomain, vecdx, dataGhostPt));


  pout() << "inititializing data"   << endl;
  
  shared_ptr<LevelData<EBGraph> > graphs = geoserv->getGraphs(domain.domainBox());
  DisjointBoxLayout& grids = vecgrids[0];
  EBLevelBoxData<CELL, DIM>  velo(grids, dataGhostIV, graphs);
  EBLevelBoxData<CELL, DIM>  gphi(grids, dataGhostIV, graphs);

  initializeData(velo, grids, dx, geomCen, geomRad, max_vel_mag, max_vel_rad);

  EBCCProjector proj(brit, geoserv, grids, domain.domainBox(), dx, dataGhostIV);
  Real tol = 1.0e-8;
  unsigned int maxiter = 27;
  proj.project(velo, gphi, tol, maxiter);

  EBLevelBoxData<CELL,   1>&  divergence = proj.getRHSHolder();
  proj.divergence(divergence, velo);
  
  Real divnorm = divergence.maxNorm(0);
  
  pout() << "max norm of post projection divergence(vel) = " << divnorm << endl;
   
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
    runProjection(a_argc, a_argv);
  }

  pout() << "printing time table " << endl;
  CH_TIMER_REPORT();
#ifdef CH_MPI
  pout() << "about to call MPI Finalize" << std::endl;
  MPI_Finalize();
#endif
  return 0;
}


