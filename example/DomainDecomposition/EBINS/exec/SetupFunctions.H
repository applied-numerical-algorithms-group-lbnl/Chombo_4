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
#if DIM==2
  int ideb = 0;
  IntVect ivdeb(1, 12);
  IntVect curiv(a_p[0], a_p[1]);
  if(curiv == ivdeb)
  {
    ideb = 1;
  }
#endif
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
  int whichvelo = 0;
  int whichscal = 0;
  ParmParse pp;
  pp.query("which_velo", whichvelo);
  if(whichvelo == 0)
  {
    pout() << "calling initializevel for velocity" << endl;
  }
  else if(whichvelo == -1)
  {
    pout() << "calling initializevelconst for velocity" << endl;
  }
  else
  {
    MayDay::Error("bogus velo");
  }

  pp.query("which_scal", whichscal);
  if(whichscal == 0)
  {
    pout() << "calling initializespot for scalar " << endl;
  }
  else if(whichscal == -1)
  {
    pout() << "calling initializeline for scalar" << endl;
  }
  else
  {
    MayDay::Error("bogus scal");
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
        ebforallInPlace_i(numflopsvelo, "IntializeVCell", InitializeVCell,  velobox, 
                          velofab, a_geomCen, a_geomRad, a_maxVelMag, a_maxVelRad, a_dx); // 
        ideb++;
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

#include "Chombo_NamespaceFooter.H"

#endif
