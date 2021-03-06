#ifndef __SetupFunctionss__
#define __SetupFunctionss__
#define PI 3.141592653589793
#include "DebugFunctions.H"
#include "Chombo_NamespaceHeader.H"
using Proto::Var;
using std::cout;
using std::endl;
using std::shared_ptr;
using Proto::Var;
using Proto::SimpleEllipsoidIF;

typedef Var<Real,DIM> Vec;
typedef Var<Real,  1> Sca;


//=================================================
PROTO_KERNEL_START 
void InitializeSpotF(int       a_p[DIM],
                     Var<Real, 1>       a_phi,
                     Real      a_X0,
                     Real      a_rad,
                     Real      a_dx)
{
  Real val = 0;
  Real rlocsq = 0;
  for(int idir = 0; idir < DIM; idir++)
  {
    Real xrel = (a_p[idir] + 0.5)*a_dx - a_X0;
    rlocsq += xrel*xrel;
  }
  Real rloc = sqrt((Real) rlocsq);
  
  if(rloc < a_rad)
  {
    Real carg = 0.5*PI*(rloc)/a_rad;
    Real cosval = cos(carg);
    val = cosval*cosval;
  }
  a_phi(0) = val;
}
PROTO_KERNEL_END(InitializeSpotF, InitializeSpot)


//=================================================
PROTO_KERNEL_START 
void InitializeLineF(int          a_p[DIM],
                     Var<Real, 1> a_phi,
                     Real      a_px,
                     Real      a_rad,
                     Real      a_dx)
{
  Real val = 0;
  Real xloc[DIM];
  for(int idir = 0; idir < DIM; idir++)
  {
    xloc[idir] = (a_p[idir] + 0.5)*a_dx;
    
  }
  
  Real funcval = std::abs(a_px - xloc[0]);
  Real rloc = funcval;
  
  if(rloc < a_rad && (xloc[1] < 0.75) )
  {
    Real carg = 0.5*PI*(rloc)/a_rad;
    Real cosval = cos(carg);
    val = std::abs(cosval);
//    val = 1.;
  }
  a_phi(0) = val;
}
PROTO_KERNEL_END(InitializeLineF, InitializeLine)

//=================================================
PROTO_KERNEL_START 
void InitializeVCellF(int       a_p[DIM],
                      Var<Real, DIM>       a_vel,
                      Real      a_cen,  //geom center
                      Real      a_rad,  //geom rad
                      Real      a_mag,  //max vel
                      Real      a_maxr, //radius for max vel
                      Real      a_dx)
{

  Real rlocsq = 0;
  Real xrel[DIM];
  for(int idir = 0; idir < DIM; idir++)
  {
    xrel[idir] = (a_p[idir] + 0.5)*a_dx - a_cen;
    rlocsq += xrel[idir]*xrel[idir];
  }

  Real rloc = sqrt((Real) rlocsq);
  Real carg = 0.5*PI*(rloc - a_maxr)/a_maxr;
  Real cosval = cos(carg);
  Real velmag = a_mag*cosval*cosval;
//debug
//  velmag = 1;
//end debug
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
PROTO_KERNEL_START 
void InitializeVCellConstF(int       a_p[DIM],
                           Var<Real, DIM>       a_vel,
                           Real      a_xvel,
                           Real      a_yvel)
{
  a_vel(0) =  a_xvel;
  a_vel(1) =  a_yvel;
#if DIM==3  
  a_vel(2) =  0;
#endif
//end debug
}
PROTO_KERNEL_END(InitializeVCellConstF, InitializeVCellConst)
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
    bool inside_regular;
    pp.get("inside_regular",inside_regular);
    retval = shared_ptr<BaseIF>(new SimpleEllipsoidIF(ABC, X0, a_geomRad, inside_regular));
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
void defineGeometry(Vector<DisjointBoxLayout >      & a_grids,
                    Vector<Box>                     & a_vecdom,
                    Vector<Real>                    & a_vecdx,
                    shared_ptr<GeometryService<2> > & a_geoserv,
                    Real  & a_geomCen,
                    Real  & a_geomRad,
                    int   & a_whichGeom,
                    const Real                      & a_dx,
                    const int                       & a_nx)
{
  CH_TIME("defineGeometry");
  pout() << "defining geometry" << endl;

  ParmParse pp;
  int maxGrid = 32;
  pp.get("maxGrid", maxGrid);
  pout() << "maxGrid  = " << maxGrid  << endl;

  IntVect domLo = IntVect::Zero;
  IntVect domHi  = (a_nx - 1)*IntVect::Unit;

// EB and periodic do not mix
  Box domain(domLo, domHi);


  pout() << "creating grid hierarchy" << endl;
  GeometryService<2>::generateGrids(a_grids, domain, maxGrid);
  a_vecdom.resize(a_grids.size());
  a_vecdx.resize( a_grids.size());
  a_vecdom[0] = a_grids[0].physDomain().domainBox();
  PR_assert(a_vecdom[0]==domain)
  a_vecdx[0]  = a_dx;
  for(int ilev = 1; ilev < a_grids.size(); ilev++)
  {

    a_vecdx[ ilev] = 2.*a_vecdx[ilev-1];
    a_vecdom[ilev] = a_grids[ilev].physDomain().domainBox();
  }
    


  int geomGhost = 4;
  RealVect origin = RealVect::Zero();

  pout() << "creating implicit function" << endl;
  
  shared_ptr<BaseIF>  impfunc = getImplicitFunction(a_geomCen,
                                                    a_geomRad,
                                                    a_whichGeom);

  pout() << "creating geometry service" << endl;
  a_geoserv  = shared_ptr<GeometryService<2> >(new GeometryService<2>(impfunc, origin, a_dx, domain, a_grids, geomGhost));
}

//=================================================
void initializeData(EBINS                       &  a_ins,
                    const DisjointBoxLayout     &  a_grids,
                    const Real                  &  a_dx,
                    const Real                  &  a_geomCen,
                    const Real                  &  a_geomRad,
                    const Real                  &  a_blobCen,
                    const Real                  &  a_blobRad,
                    const Real                  &  a_maxVelMag,
                    const Real                  &  a_maxVelRad,
                    EBIBC  a_ibc)
{
  CH_TIME("InitialData");
  auto& scalcell = *a_ins.m_scal;
  auto& velocell = *a_ins.m_velo;
  auto& species  = a_ins.m_species;
  DataIterator dit = a_grids.dataIterator();

  string whichvelo = a_ibc.m_veloIC;
  string whichscal = a_ibc.m_scalIC;

  if(whichvelo == string("zero"))
  {
    pout() << "zero velocity in the interior (potential flow ic)" << endl;
  }
  else if(whichvelo == string("vortex_0"))
  {
    pout() << "vortex ic " << endl;
  }
  else
  {
    MayDay::Error("bogus velo ic");
  }

  if(whichscal == string("blob_scal_0"))
  {
    pout() << "calling initializespot for scalar " << endl;
  }
  else if(whichscal == string("line_scal_0"))
  {
    pout() << "calling initializeline for scalar" << endl;
  }
  else
  {
    MayDay::Error("bogus scal ic");
  }

  pout() << "just making all the species spots" << endl;
  for(int ibox = 0; ibox < dit.size(); ibox++)
  {
    {
      auto& scalfab = scalcell[dit[ibox]];
      unsigned long long int numflopsscal = 5*DIM +3;
      Bx scalbox = scalfab.box();

      if(whichscal == string("blob_scal_0"))
      {
        ebforallInPlace_i(numflopsscal, "IntializeSpot", InitializeSpot,  scalbox, 
                          scalfab, a_blobCen, a_blobRad, a_dx);
      }
      else if(whichscal == string("line_scal_0"))
      {
        vector<double> v_norm;
        ParmParse pp;
        pp.getarr("geom_normal", v_norm, 0, DIM);
        ebforallInPlace_i(numflopsscal, "IntializeLine", InitializeLine,  scalbox, 
                          scalfab, a_blobCen, a_blobRad, a_dx);
      }
      else
      {
        MayDay::Error("unknown scalar ic");
      }
      for(int ispec = 0; ispec < species.size(); ispec++)
      {
        auto& spec = *species[ispec];
        auto& specfab = spec[dit[ibox]];

        ebforallInPlace_i(numflopsscal, "IntializeSpot", InitializeSpot,  scalbox, 
                          specfab, a_blobCen, a_blobRad, a_dx);
      }
    
    }

    
    auto& velofab = velocell[dit[ibox]];
    Bx velobox = velofab.box();
    if(whichvelo == string("vortex_0"))
    {
      unsigned long long int numflopsvelo = (DIM+5)*DIM +4;
      ebforallInPlace_i(numflopsvelo, "IntializeVCell", InitializeVCell,  velobox, 
                        velofab, a_geomCen, a_geomRad, a_maxVelMag, a_maxVelRad, a_dx); // 
    }
    else if(whichvelo == string("zero"))
    {
      Real yvel =  0;
      Real xvel =  0;
      unsigned long long int numflopsvelo = 0;
      //setting to zero allows the potential flow solution as an initial condition.
      ebforallInPlace_i(numflopsvelo, "IntializeVCellConst", InitializeVCellConst,  velobox, 
                        velofab, xvel, yvel);
    }
    else
    {
      MayDay::Error("unknown velo ic");
    }
  }
}

#include "Chombo_NamespaceFooter.H"

#endif
