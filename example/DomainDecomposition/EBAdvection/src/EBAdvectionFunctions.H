#ifndef __EBAdvectionFunctions__
#define __EBAdvectionFunctions__
#define PI 3.141592653589793
#include "Chombo_NamespaceHeader.H"
using Proto::Var;
#ifdef PROTO_CUDA
__device__
#endif
Real minVal(Real vala, Real valb)
{
  Real minval;
  if(vala < valb)
  {
    minval = vala;
  }
  else
  {
    minval = valb;
  }
  return minval;
}

#ifdef PROTO_CUDA
__device__
#endif
Real absVal(Real vala)
{
  Real absval;
  if(vala < 0.0)
  {
    absval = -vala;
  }
  else
  {
    absval =  vala;
  }
  return absval;
}


#ifdef PROTO_CUDA
__device__
#endif
Real cpySign(Real vala, Real signval)
{
  Real retval;
  if(signval < 0.0)
  {
    retval = -vala;
  }
  else
  {
    retval =  vala;
  }
  return retval;
}

////
PROTO_KERNEL_START 
void HybridDivergenceF(Var<Real, 1>    a_hybridDiv,
                       Var<Real, 1>    a_kappaDiv,
                       Var<Real, 1>    a_nonConsDivF,
                       Var<Real, 1>    a_deltaM,
                       Var<Real, 1>    a_kappa)
{
  a_hybridDiv(0) = a_kappaDiv(0) +(1. - a_kappa(0))*a_nonConsDivF(0);
  a_deltaM(0)    =
    a_kappaDiv(0)*   (1. - a_kappa(0)) +
    a_nonConsDivF(0)*(1. - a_kappa(0))*a_kappa(0);
//debug
//  a_hybridDiv(0) = kappaConsDiv;
//  a_deltaM(0)    = 0;
//end debug
}
PROTO_KERNEL_END(HybridDivergenceF, HybridDivergence)

//// pick scalup = upwind scalar and define flux = scalup*faceVel
PROTO_KERNEL_START 
void GetUpwindFluxF(Var<Real, 1>    a_flux,
                    Var<Real, 1>    a_scalLo,
                    Var<Real, 1>    a_scalHi,
                    Var<Real, 1>    a_faceVel)
{
  Real upwindScal;
  if(a_faceVel(0) > 0.)
  {
    upwindScal = a_scalLo(0);
  }
  else
  {
    upwindScal = a_scalHi(0);
  }
  a_flux(0) = upwindScal*(a_faceVel(0));
}
PROTO_KERNEL_END(GetUpwindFluxF, GetUpwindFlux)
////
PROTO_KERNEL_START 
void AdvanceScalarF(Var<Real, 1>    a_scal,
                    Var<Real, 1>    a_divF,
                    Real            a_dt)
{
  a_scal(0) = a_scal(0) - a_dt*a_divF(0);
}
PROTO_KERNEL_END(AdvanceScalarF, AdvanceScalar)
///
PROTO_KERNEL_START 
void ExtrapolateScalF(Var<Real, 1>       a_scal_imh_nph,
                      Var<Real, 1>       a_scal_iph_nph,
                      Var<Real, 1>       a_scal_ijk_n,
                      Var<Real, DIM>     a_slopeLo,
                      Var<Real, DIM>     a_slopeHi,
                      Var<Real, DIM>     a_ccVel,
                      unsigned int       a_dir,
                      Real               a_dt)
{
  Real mono = (a_slopeLo(a_dir) * a_slopeHi(a_dir));
  Real vanLeer = 0;
  if(mono >= 0)
  {
    Real cent = 0.5*(a_slopeLo(0) + a_slopeHi(0));
    vanLeer = minVal(absVal(cent)   , absVal(2.*a_slopeHi(a_dir)));
    vanLeer = minVal(absVal(vanLeer), absVal(2.*a_slopeLo(a_dir)));
    vanLeer = cpySign(vanLeer, cent);
  }
  //standard BCG Godunov dance to extrapolate in space and time
  a_scal_iph_nph(0) = a_scal_ijk_n(0) + 0.5*vanLeer*( 1.0 + a_dt*a_ccVel(a_dir));
  a_scal_imh_nph(0) = a_scal_ijk_n(0) + 0.5*vanLeer*(-1.0 + a_dt*a_ccVel(a_dir));
  for(int idir = 0; idir < DIM; idir++)
  {
    if(idir != a_dir)
    {
      Real upwind = a_slopeLo(idir);
      if(a_ccVel(idir) < 0)
      {
        upwind = a_slopeHi(idir);
      }
      a_scal_imh_nph(0) = a_scal_imh_nph(0) + 0.5*a_dt*upwind;
      a_scal_iph_nph(0) = a_scal_iph_nph(0) + 0.5*a_dt*upwind;
    }
  }
}
PROTO_KERNEL_END(ExtrapolateScalF, ExtrapolateScal)
#include "Chombo_NamespaceFooter.H"

#endif
