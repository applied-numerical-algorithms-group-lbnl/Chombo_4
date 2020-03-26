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
inline Real getSpotVal(int       a_p[DIM],
                       Real      a_X0,
                       Real      a_rad,
                       Real      a_dx)
{
  Real val =  0;
  Real rlocsq = 0;
  for(int idir = 0; idir < DIM; idir++)
  {
    Real xrel = (a_p[idir] + 0.5)*a_dx - a_X0;
    rlocsq += xrel*xrel;
  }
  Real rloc = sqrt((Real) rlocsq);
  

  Real carg = 0.5*PI*(rloc)/a_rad;
  Real cosval = cos(carg);
  val = cosval*cosval;
  return val;
}


//=================================================
PROTO_KERNEL_START 
void InitializeSpotF(int       a_p[DIM],
                     Var<Real, 1>       a_phi,
                     Real      a_X0,
                     Real      a_rad,
                     Real      a_dx)
{
  Real val =  getSpotVal(a_p, a_X0, a_rad, a_dx);
  a_phi(0) = val;
}
PROTO_KERNEL_END(InitializeSpotF, InitializeSpot)



//=================================================
PROTO_KERNEL_START 
void checkSpotF(int       a_p[DIM],
                     Var<Real, 1>       a_phi,
                     Real      a_X0,
                     Real      a_rad,
                     Real      a_dx)
{
  Real val =  getSpotVal(a_p, a_X0, a_rad, a_dx);
  Real phi = a_phi(0);
  Real diff = phi - val;
  if(std::abs(diff) > 1.0e-6)
  {
    using std::cout;
    using std::endl;
    cout << "bogus value at point: (" << a_p[0] << "," <<  a_p[1] << ")" << endl;;
  }
}
PROTO_KERNEL_END(checkSpotF, checkSpot)


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


#include "Chombo_NamespaceFooter.H"

#endif
