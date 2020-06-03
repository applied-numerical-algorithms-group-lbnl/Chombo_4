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
void InitializeSpotDIMF(int       a_p[DIM],
                        Var<Real, DIM>       a_phi,
                        Real      a_X0,
                        Real      a_rad,
                        Real      a_dx)
{
  Real val =  getSpotVal(a_p, a_X0, a_rad, a_dx);
  for(unsigned int idir = 0; idir < DIM; idir++)
  {
    a_phi(idir) = val*Real(idir+1);
  }
}
PROTO_KERNEL_END(InitializeSpotDIMF, InitializeSpotDIM)



//=================================================
PROTO_KERNEL_START 
void checkSpotDIMF(int       a_p[DIM],
                   Var<Real, DIM>       a_phi,
                   Real      a_X0,
                   Real      a_rad,
                   Real      a_dx)
{
  Real val0 =  getSpotVal(a_p, a_X0, a_rad, a_dx);
  for(unsigned int idir = 0; idir < DIM; idir++)
  {
    Real val = val0*Real(idir+1);
    Real phi = a_phi(idir);
    Real diff = phi - val;
    if(std::abs(diff) > 1.0e-6)
    {
      using std::cout;
      using std::endl;
      cout << "bogus value at variable " << idir << ", point: (" << a_p[0] << "," <<  a_p[1] << ")" << endl;;
    }
  }
}
PROTO_KERNEL_END(checkSpotDIMF, checkSpotDIM)
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




#include "Chombo_NamespaceFooter.H"

#endif
