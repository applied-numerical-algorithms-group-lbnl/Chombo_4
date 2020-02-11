#ifndef __EBAdvectionFunctions__
#define __EBAdvectionFunctions__
#define PI 3.141592653589793
#include "Chombo_NamespaceHeader.H"
using Proto::Var;
#ifdef PROTO_CUDA
__device__
#endif
inline Real minVal(Real vala, Real valb)
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
inline Real absVal(Real vala)
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
inline Real cpySign(Real vala, Real signval)
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
                       Var<Real, 1>    a_nonConsDiv,
                       Var<Real, 1>    a_deltaM,
                       Var<Real, 1>    a_kappa)
{
  Real kapdiv = a_kappaDiv(0);
  Real  ncdiv = a_nonConsDiv(0);
  Real kappa  = a_kappa(0);
  Real hybdiv = kapdiv +(1. - kappa)*ncdiv;
  Real deltam = (1. - kappa)*kapdiv  - kappa*(1. - kappa)*ncdiv;

  a_hybridDiv(0) = hybdiv;
  a_deltaM(0) = deltam;
}
PROTO_KERNEL_END(HybridDivergenceF, HybridDivergence)


////
PROTO_KERNEL_START 
void HybridDivergencePtF(int a_pt[DIM],
                         Var<Real, 1>    a_hybridDiv,
                         Var<Real, 1>    a_kappaDiv,
                         Var<Real, 1>    a_nonConsDiv,
                         Var<Real, 1>    a_deltaM,
                         Var<Real, 1>    a_kappa)     
{

//begin debug
#if DIM==2
  Point debpt (17, 19);
  int ideb = 0;
  if((a_pt[0] == debpt[0]) && (a_pt[1] == debpt[1]))
  {
    ideb = 1;
  }
#endif  
//end   debug
  HybridDivergence(a_hybridDiv, 
                   a_kappaDiv,  
                   a_nonConsDiv,
                   a_deltaM,    
                   a_kappa);     

}
PROTO_KERNEL_END(HybridDivergencePtF, HybridDivergencePt)

//// pick scalup = upwind scalar and define flux = scalup*faceVel
PROTO_KERNEL_START 
void GetUpwindFluxF(Var<Real, 1>    a_flux,
                    Var<Real, 1>    a_scalLo,
                    Var<Real, 1>    a_scalHi,
                    Var<Real, 1>    a_faceVel)
{
  Real upwindScal;
  Real advVel  = a_faceVel(0);
  Real scaLo = a_scalLo(0);
  Real scaHi = a_scalHi(0);
  if(advVel > 0.)
  {
    upwindScal = scaLo;
  }
  else
  {
    upwindScal = scaHi;
  }
  a_flux(0) = upwindScal*advVel;
}
PROTO_KERNEL_END(GetUpwindFluxF, GetUpwindFlux)
//// pick scalup = upwind scalar and define flux = scalup*faceVel
PROTO_KERNEL_START 
void FluxMultiplyF(Var<Real, 1>    a_flux,
                   Var<Real, 1>    a_scal,
                   Var<Real, 1>    a_faceVel)
{
  Real upwindScal = a_scal(0);
  Real advVel  = a_faceVel(0);

  a_flux(0) = upwindScal*advVel;
}
PROTO_KERNEL_END(FluxMultiplyF, FluxMultiply)
//// pick scalup = upwind scalar and define flux = scalup*faceVel
PROTO_KERNEL_START 
void UpwindedF(Var<Real, 1>    a_upwind,
               Var<Real, 1>    a_scalLo,
               Var<Real, 1>    a_scalHi,
               Var<Real, 1>    a_faceVel)
{
  Real upwindScal;
  Real advVel  = a_faceVel(0);
  Real scaLo = a_scalLo(0);
  Real scaHi = a_scalHi(0);
  if(advVel > 0.)
  {
    upwindScal = scaLo;
  }
  else
  {
    upwindScal = scaHi;
  }
  a_upwind(0) = upwindScal;
}
PROTO_KERNEL_END(UpwindedF, Upwinded)
PROTO_KERNEL_START 
void CopiedF(Var<Real, 1>    a_dst,
             Var<Real, 1>    a_src)
{
  a_dst(0) = a_src(0);
}
PROTO_KERNEL_END(CopiedF, Copied)
//// pick scalup = upwind scalar and define flux = scalup*faceVel
PROTO_KERNEL_START 
void GetUpwindFluxPtF(int             a_pt[DIM],
                      Var<Real, 1>    a_flux,
                      Var<Real, 1>    a_scalLo,
                      Var<Real, 1>    a_scalHi,
                      Var<Real, 1>    a_faceVel)
{
//begin debug
#if DIM==2
  Point debpt (17, 19);
  int ideb = 0;
  if((a_pt[0] == debpt[0]) && (a_pt[1] == debpt[1]))
  {
    ideb = 1;
  }
#endif  
//end   debug
  GetUpwindFlux(a_flux,
                a_scalLo,
                a_scalHi,
                a_faceVel);
}
PROTO_KERNEL_END(GetUpwindFluxPtF, GetUpwindFluxPt)
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
                      Var<Real, 1  >     a_source,
                      unsigned int       a_dir,
                      Real               a_dt,
                      Real               a_dx)
{
  //van leer slopes
  Real vanLeer = 0;
  Real sloplo = a_slopeLo(a_dir);
  Real slophi = a_slopeHi(a_dir);
  Real cent = 0.5*(sloplo + slophi);

  Real mono   = sloplo*slophi;
  if(mono >= 0)
  {
    vanLeer = minVal(absVal(cent)   , absVal(2.*a_slopeHi(a_dir)));
    vanLeer = minVal(absVal(vanLeer), absVal(2.*a_slopeLo(a_dir)));
    vanLeer = cpySign(vanLeer, cent);
  }

  //normal predictor
  Real dtbydx    = a_dt/a_dx;
  Real ccveldir  = a_ccVel(a_dir);
  //standard BCG Godunov dance to extrapolate in space and time
  //normal predictor
  Real ccval = a_scal_ijk_n(0) + 0.5*a_dt*a_source(0);
  a_scal_iph_nph(0)  = ccval + 0.5*vanLeer*minVal((1. - ccveldir*dtbydx), 1.0);
  a_scal_imh_nph(0)  = ccval - 0.5*vanLeer*minVal((1. + ccveldir*dtbydx), 1.0);

  //transverse predictor
  for(int idir = 0; idir < DIM; idir++)
  {
    if(idir != a_dir)
    {
      Real upwind = a_slopeLo(idir);
      if(a_ccVel(idir) < 0)
      {
        upwind = a_slopeHi(idir);
      }
      Real tanvel = a_ccVel(idir);
      a_scal_imh_nph(0) = a_scal_imh_nph(0) - 0.5*dtbydx*upwind*tanvel;
      a_scal_iph_nph(0) = a_scal_iph_nph(0) - 0.5*dtbydx*upwind*tanvel;
    }
  }
}
PROTO_KERNEL_END(ExtrapolateScalF, ExtrapolateScal)


PROTO_KERNEL_START 
void ExtrapolateScalPtF(int                a_pt[DIM],
                        Var<Real, 1>       a_scal_imh_nph,
                        Var<Real, 1>       a_scal_iph_nph,
                        Var<Real, 1>       a_scal_ijk_n,
                        Var<Real, DIM>     a_slopeLo,
                        Var<Real, DIM>     a_slopeHi,
                        Var<Real, DIM>     a_ccVel,
                        Var<Real, 1>       a_source,
                        unsigned int       a_dir,
                        Real               a_dt,
                        Real               a_dx)
{

//begin debug
#if DIM==2
  Point debpt (17, 18);
  int ideb = 0;
  if((a_dir==1) && (a_pt[0] == debpt[0]) && (a_pt[1] == debpt[1]))
  {
    ideb = 1;
  }
#endif  
//end   debug
  ExtrapolateScal(a_scal_imh_nph,
                  a_scal_iph_nph,
                  a_scal_ijk_n,
                  a_slopeLo,
                  a_slopeHi,
                  a_ccVel,
                  a_source,
                  a_dir,
                  a_dt,
                  a_dx);
}
PROTO_KERNEL_END(ExtrapolateScalPtF, ExtrapolateScalPt)
#include "Chombo_NamespaceFooter.H"

#endif
