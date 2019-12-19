#include <cmath>
#include <memory>
#include "EBMACProjector.H"
#include "Chombo_NamespaceHeader.H"



  
/// 
EBMACProjector::
EBMACProjector(shared_ptr<EBEncyclopedia<2, Real> >   & a_brit,
               shared_ptr<GeometryService<2> >        & a_geoserv,
               const DisjointBoxLayout                & a_grids,
               const Box                              & a_domain,
               const Real                             & a_dx,
               const IntVect                          & a_nghost)
{
  m_dx         = a_dx;
  m_grids      = a_grids;
  m_domain     = a_domain;
  m_nghost     = a_nghost;
  m_grids      = a_grids;      
  m_domain     = a_domain;     
  m_brit       = a_brit;

  defineInternals();
}


/// advance one time step (via Trebotich et al.) in  an eb context
void 
EBMACProjector::
project(EBLevelFluxData<1>   & a_velo,
        EBLevelFluxData<1>   & a_gphi)
{
}

void 
EBMACProjector::
divergence(EBLevelBoxData<CELL, 1> & a_divu,
           EBLevelFluxData<1>      & a_velo)
{
}

EBLevelBoxData<CELL, 1>& 
EBMACProjector::
getRHSHolder()
{
  return m_rhs;
}
#include "Chombo_NamespaceFooter.H"


