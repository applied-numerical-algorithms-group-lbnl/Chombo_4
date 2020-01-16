#include "EBParabolicIntegrators.H"
#include "Chombo_NamespaceHeader.H"

//diff comes in holding phinew
//leaves holding (phinew-phiold)/dt

typedef Proto::Var<Real, 1> Sca;

///
PROTO_KERNEL_START 
void  getConsDiffF(Sca     a_diff,
                   Sca     a_phiold,
                   Real    a_dt)
{
  a_diff(0) = (a_diff(0) - a_phiold(0))/a_dt;
}
PROTO_KERNEL_END(getConsDiffF, getConsDiff)
///
PROTO_KERNEL_START 
void  setEulerRHSF(Sca     a_rhs,
                   Sca     a_phiold,
                   Sca     a_source,
                   Sca     a_kappa,
                   Real    a_dt)
{
  //source is presumed to already be kappa weighted
  a_rhs(0) = a_kappa(0)*a_phiold(0) + a_dt*a_source(0);
}
PROTO_KERNEL_END(setEulerRHSF, setEulerRHS)
///

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

#include "Chombo_NamespaceFooter.H"
