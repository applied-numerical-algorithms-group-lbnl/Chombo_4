#ifndef __EBMultiGridFunctions__
#define __EBMultiGridFunctions__
#include "Chombo_NamespaceHeader.H"
using Proto::Var;

/****/
//lph comes in holding beta*div(F)--leaves holding alpha phi + beta div(F)
PROTO_KERNEL_START 
unsigned int  addAlphaPhiPtF(int a_pt[DIM],
                             Var<Real, 1>     a_lph,
                           Var<Real, 1>     a_phi,
                           Var<Real, 1>     a_kap,
                           Real    a_alpha,
                           Real    a_beta)
{
  //kappa and beta are already in lph
  //kappa because we did not divide by kappa
  //beta was sent in to ebstencil::apply
  

  
  Real kappdiv = a_lph(0);
  Real bkdivF  = a_beta*kappdiv;
  Real phival  = a_phi(0);
  Real kapval  = a_kap(0);
  a_lph(0) = a_alpha*phival*kapval + bkdivF;
  return 0;
}
PROTO_KERNEL_END(addAlphaPhiPtF, addAlphaPhiPt)
PROTO_KERNEL_START 
unsigned int  addAlphaPhiF(Var<Real, 1>     a_lph,
                           Var<Real, 1>     a_phi,
                           Var<Real, 1>     a_kap,
                           Real    a_alpha,
                           Real    a_beta)
{
  //kappa and beta are already in lph
  //kappa because we did not divide by kappa
  //beta was sent in to ebstencil::apply
  Real kappdiv = a_lph(0);
  Real bkdivF  = a_beta*kappdiv;
  Real phival  = a_phi(0);
  Real kapval  = a_kap(0);
  a_lph(0) = a_alpha*phival*kapval + bkdivF;
  return 0;
}
PROTO_KERNEL_END(addAlphaPhiF, addAlphaPhi)


/****/
//res comes in holding lphi.   leaves holding res= rhs-lphi
PROTO_KERNEL_START 
unsigned int  subtractRHSF(Var<Real, 1>     a_res,
                           Var<Real, 1>     a_rhs)
{
  a_res(0) = a_rhs(0) - a_res(0);
  return 0;
}
PROTO_KERNEL_END(subtractRHSF, subtractRHS)
//res comes in holding lphi.   leaves holding res= rhs-lphi
PROTO_KERNEL_START 
unsigned int  copyDiagF(Var<Real, 1>     a_dst,
                        Var<Real, 1>     a_src)
{
  a_dst(0) = a_src(0);
  return 0;
}
PROTO_KERNEL_END(copyDiagF, copyDiag)
/****/
/****/
PROTO_KERNEL_START 
unsigned int  gsrbResidF(int     a_pt[DIM],
                         Var<Real, 1>     a_phi,
                         Var<Real, 1>     a_res,
                         Var<Real, 1>     a_diag,
                         Var<Real, 1>     a_kappa,
                         Real    a_alpha,
                         Real    a_beta,
                         Real    a_dx,
                         int     a_iredBlack)
{
  int sumpt = 0;
  for(int idir = 0; idir < DIM; idir++)
  {
    sumpt += a_pt[idir];
  }
  if(sumpt%2 == a_iredBlack)
  {
    static const Real safety = 1.0;
    Real diagval = a_diag(0);
    Real kappval = a_kappa(0);
    
    Real realdiag = kappval*a_alpha + a_beta*diagval;
//    Real regudiag = a_alpha - 2*DIM*a_beta/(a_dx*a_dx);
    Real lambda = safety/realdiag;
//    Real reglam = safety/regudiag;
    Real phival = a_phi(0);
    Real resval = a_res(0);
//    if(std::abs(lambda) > std::abs(reglam)) lambda= reglam;
    a_phi(0) = phival + lambda*resval;
  }
  return 0;
}
PROTO_KERNEL_END(gsrbResidF, gsrbResid)
/****/
/****/
PROTO_KERNEL_START 
unsigned int  jacobiResidF(int     a_pt[DIM],
                           Var<Real, 1>     a_phi,
                           Var<Real, 1>     a_res,
                           Var<Real, 1>     a_diag,
                           Var<Real, 1>     a_kappa,
                           Real    a_alpha,
                           Real    a_beta,
                           Real    a_dx)
{
    Real phival = a_phi(0);
    Real resval = a_res(0);
    Real kappval = a_kappa(0);
    double lambda = -(a_dx*a_dx)/(16.0);
    a_phi(0) = phival + lambda*resval;
  return 0;
}
PROTO_KERNEL_END(jacobiResidF, jacobiResid)
/****/
#include "Chombo_NamespaceFooter.H"
/****/
#endif
