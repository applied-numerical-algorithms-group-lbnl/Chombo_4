#include "EBParabolicIntegrators.H"
#include "Chombo_NamespaceHeader.H"

//diff comes in holding phinew
//leaves holding (phinew-phiold)/dt

#if 1


///
PROTO_KERNEL_START 
void  getConsDiffF(Proto::Var<Real, 1>     a_diff,
                   Proto::Var<Real, 1>     a_phiold,
                   Real    a_dt)
{
  a_diff(0) = (a_diff(0) - a_phiold(0))/a_dt;
}
PROTO_KERNEL_END(getConsDiffF, getConsDiff)
///this is the conservative way to advance a solution with diffusion.
void 
BaseEBParabolic::
computeDiffusion( EBLevelBoxData<CELL, 1>       &  a_diffusionTerm,
                  const EBLevelBoxData<CELL, 1> &  a_phiold,
                  const Real                    &  a_diffCoef,
                  const Real                    &  a_dt, 
                  const Real                    &  a_tolerance,
                  const unsigned int            &  a_maxIterations)

{
  //this puts phinew into diffterm
  advanceOneStep(a_diffusionTerm, a_phiold, a_diffCoef, 
                 a_dt, a_tolerance, a_maxIterations);
  //this sets diffusion term= (phinew-phiold)/a_dt
  DataIterator dit = m_grids.dataIterator();
  for(unsigned int ibox = 0; ibox < dit.size(); ++ibox)
  {
    auto      & difffab =  a_diffusionTerm[dit[ibox]];
    const auto& phiofab =         a_phiold[dit[ibox]];
    Bx grid = ProtoCh::getProtoBox(m_grids[dit[ibox]]);

    unsigned long long int numflopspt = 3;

    ebforallInPlace(numflopspt, "getConsDiff", getConsDiff, grid, difffab, phiofab, a_dt);
  }
}
///
PROTO_KERNEL_START 
void  setEulerRHSF(Proto::Var<Real, 1>     a_rhs,
                   Proto::Var<Real, 1>     a_phiold,
                   Proto::Var<Real, 1>     a_source,
                   Proto::Var<Real, 1>     a_kappa,
                   Real    a_dt)
{
  //source is presumed to already be kappa weighted
  a_rhs(0) = a_kappa(0)*(a_phiold(0) + a_dt*a_source(0));
}
PROTO_KERNEL_END(setEulerRHSF, setEulerRHS)

///
void 
EBBackwardEuler::
advanceOneStep( EBLevelBoxData<CELL, 1>       &  a_phi,
                const EBLevelBoxData<CELL, 1> &  a_source,
                const Real                    &  a_diffCoef,
                const Real                    &  a_dt,
                const Real                    &  a_tolerance,
                const unsigned int            &  a_maxIterations)
{
  EBLevelBoxData<CELL, 1> & kappa = m_diffusionSolver->getKappa();

  //set rhs = kappa*phi^n + dt*source 
  //the source is assumed to already be kappa weighted
  DataIterator dit = m_grids.dataIterator();
  for(unsigned int ibox = 0; ibox < dit.size(); ++ibox)
  {
    auto       & rhsfab =    m_rhs[dit[ibox]];
    auto       & phifab =    a_phi[dit[ibox]];
    const auto & soufab = a_source[dit[ibox]];
    const auto & kapfab =    kappa[dit[ibox]];

    Bx grid = ProtoCh::getProtoBox(m_grids[dit[ibox]]);

    unsigned long long int numflopspt = 3;

    ebforallInPlace(numflopspt, "setEulerRHS", setEulerRHS, grid, rhsfab, phifab, soufab, kapfab, a_dt);
  }
  
  //solve for new phi
  Real alpha = 1;
  Real beta = -a_dt*a_diffCoef;
  m_diffusionSolver->resetAlphaAndBeta(alpha, beta);
  m_diffusionSolver->solve(a_phi, m_rhs, a_tolerance, a_maxIterations);
}

///
PROTO_KERNEL_START 
void  setCrankNicRHSF(Proto::Var<Real, 1>     a_rhs,
                      Proto::Var<Real, 1>     a_phiold,
                      Proto::Var<Real, 1>     a_source,
                      Proto::Var<Real, 1>     a_kappa,
                      Proto::Var<Real, 1>     a_kappalapphi,
                      Real    a_dt,
                      Real    a_diffCoef)
{
  //source is presumed to already be kappa weighted
  a_rhs(0) = a_kappa(0)*a_phiold(0) + a_dt*a_source(0) + 0.5*a_dt*a_diffCoef*a_kappalapphi(0);
}
PROTO_KERNEL_END(setCrankNicRHSF, setCrankNicRHS)
///
void 
EBCrankNicolson::
advanceOneStep( EBLevelBoxData<CELL, 1>       &  a_phi,
                const EBLevelBoxData<CELL, 1> &  a_source,
                const Real                    &  a_diffCoef,
                const Real                    &  a_dt,
                const Real                    &  a_tolerance,
                const unsigned int            &  a_maxIterations)
{
  EBLevelBoxData<CELL, 1> & kappa = m_diffusionSolver->getKappa();

  //take kappa*laplacian(phiold)
  Real alpha = 0;
  Real beta  = 1;
  m_diffusionSolver->resetAlphaAndBeta(alpha, beta);
  m_diffusionSolver->applyOp(m_kappaLph, a_phi);

  //set rhs = kappa*phi^n + source*dt + (dt/2)*kappa*lapl(phi)
  //the source is assumed to already be kappa weighted
  DataIterator dit = m_grids.dataIterator();
  for(unsigned int ibox = 0; ibox < dit.size(); ++ibox)
  {
    auto       & rhsfab =      m_rhs[dit[ibox]];
    auto       & phifab =      a_phi[dit[ibox]];
    const auto & soufab =   a_source[dit[ibox]];
    const auto & kapfab =      kappa[dit[ibox]];
    const auto & lapfab = m_kappaLph[dit[ibox]];

    Bx grid = ProtoCh::getProtoBox(m_grids[dit[ibox]]);

    unsigned long long int numflopspt = 3;

    ebforallInPlace(numflopspt, "setCrankNicRHS", setCrankNicRHS, grid, rhsfab, phifab, soufab, kapfab, lapfab, a_dt, a_diffCoef);
  }
  
  //solve for new phi
  alpha = 1;
  beta = -0.5*a_dt*a_diffCoef;
  m_diffusionSolver->resetAlphaAndBeta(alpha, beta);
  m_diffusionSolver->solve(a_phi, m_rhs, a_tolerance, a_maxIterations);
}
///
PROTO_KERNEL_START 
void  setTGASrcF(Proto::Var<Real, 1>     a_srct,
                 Proto::Var<Real, 1>     a_input,
                 Real    a_dt)
{
  //source is presumed to already be kappa weighted
  a_srct(0) = a_dt*a_input(0);
}
PROTO_KERNEL_END(setTGASrcF, setTGASrc)

PROTO_KERNEL_START 
void  incrTGARHSF(Proto::Var<Real, 1>     a_phit,
                  Proto::Var<Real, 1>     a_rhst)
{
  a_rhst(0) = a_rhst(0) + a_phit(0);
}
PROTO_KERNEL_END(incrTGARHSF, incrTGARHS)

///
void 
EBTGA::
advanceOneStep( EBLevelBoxData<CELL, 1>       &  a_phi,
                const EBLevelBoxData<CELL, 1> &  a_source,
                const Real                    &  a_diffCoef,
                const Real                    &  a_dt,
                const Real                    &  a_tolerance,
                const unsigned int            &  a_maxIterations)
{

  //this makes m_src = dt*a_src
  DataIterator dit = m_grids.dataIterator();
  for(unsigned int ibox = 0; ibox < dit.size(); ++ibox)
  {
    const auto & inputfab =   a_source[dit[ibox]];
    auto       & srctfab  =     m_srct[dit[ibox]];
    Bx grid = ProtoCh::getProtoBox(m_grids[dit[ibox]]);

    unsigned long long int numflopspt = 1;
    ebforallInPlace(numflopspt, "setTGASrc", setTGASrc, grid, srctfab, inputfab, a_dt);
  }

  Real alpha = 1;
  Real beta  = m_mu4*a_dt*a_diffCoef;
  m_diffusionSolver->resetAlphaAndBeta(alpha, beta); 
  //this makes rhs hold       (k I + dt* mu4 k L) (S)
  //L == kappa*diffcoef*lapl(phi)
  m_diffusionSolver->applyOpNeumann(m_rhst, m_srct); 

  //this makes phit = (k I + mu3*dt kL) phin
  alpha = 1;
  beta  = m_mu3*a_dt*a_diffCoef;
  m_diffusionSolver->resetAlphaAndBeta(alpha, beta);
  m_diffusionSolver->applyOp(m_phit, a_phi);
  //this sets rhst = rhst + phit
  for(unsigned int ibox = 0; ibox < dit.size(); ++ibox)
  {
    const auto & phitfab  =   m_phit[dit[ibox]];
    auto       & rhstfab  =   m_rhst[dit[ibox]];
    Bx grid = ProtoCh::getProtoBox(m_grids[dit[ibox]]);
    unsigned long long int numflopspt = 1;
    ebforallInPlace(numflopspt, "incrTGARHS", incrTGARHS, grid, rhstfab, phitfab);
  }

  alpha = 1;
  beta = -m_mu2*a_dt*a_diffCoef;
  m_diffusionSolver->resetAlphaAndBeta(alpha, beta); 
  //sets phit = (k I - mu2* L)^-1(rhs)
  m_diffusionSolver->solve(m_phit, m_rhst, a_tolerance, a_maxIterations);

  //this uses phit as the RHS (instead of copying it it over to m_rhst)
  alpha = 1;
  beta = -m_mu1*a_dt*a_diffCoef;
  m_diffusionSolver->resetAlphaAndBeta(alpha, beta); 
  m_diffusionSolver->solve(a_phi, m_phit, a_tolerance, a_maxIterations);
}
#endif
///
#include "Chombo_NamespaceFooter.H"
