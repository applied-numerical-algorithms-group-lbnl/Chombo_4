#ifndef __EBAdvectionFunctions__
#define __EBAdvectionFunctions__
#define PI 3.141592653589793
//#include "Chombo_NamespaceHeader.H"
using Proto::Var;
#ifdef PROTO_CUDA
CUDA_DECORATION
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
CUDA_DECORATION
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
CUDA_DECORATION
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
unsigned int setFluxValF(Var<Real, 1>    a_flux,
                         Real            a_value)
{
  a_flux(0) = a_value;
  return 0;
}
PROTO_KERNEL_END(setFluxValF, setFluxVal)
////
PROTO_KERNEL_START 
unsigned int copyExtrapF(Var<Real, 1>    a_fluxval,
                         Var<Real, 1>    a_scalarlo,
                         Var<Real, 1>    a_scalarhi,
                         int             a_sideSign)
{
  if(a_sideSign > 0)
  {
    a_fluxval(0) = a_scalarlo(0);
  }
  else
  {
    a_fluxval(0) = a_scalarhi(0);
  }
  return 0;
}
PROTO_KERNEL_END(copyExtrapF, copyExtrap)

////
PROTO_KERNEL_START 
unsigned int HybridDivergenceF(Var<Real, 1>    a_hybridDiv,
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
  return 0;
}
PROTO_KERNEL_END(HybridDivergenceF, HybridDivergence)



//// pick scalup = upwind scalar and define flux = scalup*faceVel
PROTO_KERNEL_START 
unsigned int FluxMultiplyF(Var<Real, 1>    a_flux,
                           Var<Real, 1>    a_scal,
                           Var<Real, 1>    a_faceVel)
{
  Real upwindScal = a_scal(0);
  Real advVel  = a_faceVel(0);

  a_flux(0) = upwindScal*advVel;
  return 0;
}
PROTO_KERNEL_END(FluxMultiplyF, FluxMultiply)


PROTO_KERNEL_START 
unsigned int UpwindedF(Var<Real, 1>    a_upwind,
                       Var<Real, 1>    a_scalLo,
                       Var<Real, 1>    a_scalHi,
                       Var<Real, 1>    a_faceVel,
                       unsigned int a_facedir,
                       unsigned int a_curcomp,
                       unsigned int a_doingvel)
{
  Real upwindScal;
  Real advVel  = a_faceVel(0);
  Real scaLo = a_scalLo(0);
  Real scaHi = a_scalHi(0);
  Real tol  = 1.0e-12;
  bool doingvel   = (a_doingvel == 1);
  bool normalcomp = (a_curcomp == a_facedir);
  if(advVel > tol)
  {
    upwindScal = scaLo;
  }
  else if (advVel < (-tol))
  {
    upwindScal = scaHi;
  }
  else if(doingvel && normalcomp)
  {
    upwindScal = 0.0;
  }
  else
  {
    upwindScal = 0.5*(scaLo + scaHi);
  }
  a_upwind(0) = upwindScal;
  return 0;
}
PROTO_KERNEL_END(UpwindedF, Upwinded)

PROTO_KERNEL_START 
unsigned int CopiedF(Var<Real, 1>    a_dst,
                     Var<Real, 1>    a_src)
{
  a_dst(0) = a_src(0);
  return 0;
}
PROTO_KERNEL_END(CopiedF, Copied)


PROTO_KERNEL_START 
unsigned int MinionedF(Var<Real, 1>    a_minion,
                       Var<Real, 1>    a_scalar,
                       Var<Real, 1>    a_source,
                       Real a_dt)
{
  a_minion(0) = a_scalar(0) + 0.5*a_dt*a_source(0);
  return 0;
}
PROTO_KERNEL_END(MinionedF, Minioned)
//// pick scalup = upwind scalar and define flux = scalup*faceVel
////
PROTO_KERNEL_START 
unsigned int AdvanceScalarF(Var<Real, 1>    a_scal,
                            Var<Real, 1>    a_divF,
                            Real            a_dt)
{
  a_scal(0) = a_scal(0) - a_dt*a_divF(0);
  return 0;
}
PROTO_KERNEL_END(AdvanceScalarF, AdvanceScalar)
///
PROTO_KERNEL_START 
unsigned int ExtrapolateScalF(Var<Real, 1>       a_scal_imh_nph,
                              Var<Real, 1>       a_scal_iph_nph,
                              Var<Real, 1>       a_scal_ijk_n,
                              Var<Real, DIM>     a_slopeLoNor,
                              Var<Real, DIM>     a_slopeHiNor,
                              Var<Real, DIM>     a_slopeLoTan,
                              Var<Real, DIM>     a_slopeHiTan,
                              Var<Real, DIM>     a_ccVel,
                              Var<Real, 1  >     a_source,
                              unsigned int       a_dir,
                              Real               a_dt,
                              Real               a_dx)
{
  //normal slopes get straight scal.
  //tangential slopes get the minion fix
  //van leer slopes
  Real vanLeer = 0;
  Real sloplo = a_slopeLoNor(a_dir);
  Real slophi = a_slopeHiNor(a_dir);
  Real cent = 0.5*(sloplo + slophi);

  Real mono   = sloplo*slophi;
  if(mono >= 0)
  {
    vanLeer = minVal(absVal(cent)   , absVal(2.*a_slopeHiNor(a_dir)));
    vanLeer = minVal(absVal(vanLeer), absVal(2.*a_slopeLoNor(a_dir)));
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
      Real upwind = a_slopeLoTan(idir);
      if(a_ccVel(idir) < 0)
      {
        upwind = a_slopeHiTan(idir);
      }
      Real tanvel = a_ccVel(idir);
      a_scal_imh_nph(0) = a_scal_imh_nph(0) - 0.5*dtbydx*upwind*tanvel;
      a_scal_iph_nph(0) = a_scal_iph_nph(0) - 0.5*dtbydx*upwind*tanvel;
    }
  }
  return 0;
}
PROTO_KERNEL_END(ExtrapolateScalF, ExtrapolateScal)


//#include "Chombo_NamespaceFooter.H"

#endif
