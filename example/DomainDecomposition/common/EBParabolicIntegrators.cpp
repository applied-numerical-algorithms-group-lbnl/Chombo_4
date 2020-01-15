#include "EBParabolicIntegrators.H"
#include "Chombo_NamespaceHeader.H"

//diff comes in holding phinew
//leaves holding (phinew-phiold)/dt

typedef Proto::Var<Real, 1> Sca;

PROTO_KERNEL_START 
void  getConsDiffF(Sca     a_diff,
                   Sca     a_phiold,
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
                  const Real                    &  a_dt)
{
  //this puts phinew into diffterm
  advanceOneStep(a_diffusionTerm, a_phiold, a_diffCoef, a_dt);
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

#include "Chombo_NamespaceFooter.H"
