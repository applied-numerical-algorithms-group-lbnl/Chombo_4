#include "EBAdvection.H"
#include "NamespaceHeader.H"
using Proto::Var;
////=================================================
PROTO_KERNEL_START 
void HybridDivergenceF(Var<Real, 1>    a_hybridDiv,
                       Var<Real, 1>    a_kappaConsDiv,
                       Var<Real, 1>    a_nonConsDivF,
                       Var<Real, 1>    a_deltaM,
                       Var<Real, 1>    a_kappa)
{
  a_hybridDiv(0) = a_kappaConsDiv(0) + (1- a_kappa(0))*a_nonConsDivF(0);
  a_deltaM(0)    = (1-a_kappa(0))*(a_kappaConsDiv(0) - a_kappa(0)*a_nonConsDivF(0));
}
PROTO_KERNEL_END(HybridDivergenceF, HybridDivergence)
/// 
EBAdvection::
EBAdvection(shared_ptr<EBEncyclopedia<2, Real> >   & a_dictionaries,
            shared_ptr<GeometryService<2> >        & a_geoserv,
            shared_ptr<EBLevelBoxData<CELL, DIM> > & a_veloCell,
            const DisjointBoxLayout                & a_grids,
            const Box                              & a_domain,
            const Real                             & a_dx,
            const IntVect                          & a_nghostsrc, 
            const IntVect                          & a_nghostdst)
{
  m_grids  = a_grids;
  m_domain = a_domain;
  fillKappa(a_geoserv);
}

void  
EBAdvection::
fillKappa(const shared_ptr<GeometryService<2> >   & a_geoserv)
{
  const shared_ptr<LevelData<EBGraph>  > graphs = a_geoserv->getGraphs(m_domain);
  DataIterator dit = m_grids.dataIterator();
  for(int ibox = 0; ibox < dit.size(); ++ibox)
  {
    Box grid =m_grids[dit[ibox]];
    Bx  grbx = getProtoBox(grid);
    const EBGraph  & graph = (*graphs)[dit[ibox]];
    EBHostData<CELL, Real, 1> hostdat(grbx, graph);
    //fill kappa on the host then copy to the device
    a_geoserv->fillKappa(hostdat, grid, dit[ibox], m_domain);
    // now copy to the device
    EBLevelBoxData<CELL, 1>::copyToDevice(hostdat, m_kappa[dit[ibox]]);
  }
}
///
void
EBAdvection::
kappaConsDiv(EBLevelBoxData<CELL, 1>   & a_scal)
{
}
///
void
EBAdvection::
nonConsDiv(  EBLevelBoxData<CELL, 1>   & a_scal)
{
}
///
void
EBAdvection::
redistribute()
{
}
///
void 
EBAdvection::
advance(EBLevelBoxData<CELL, 1>       & a_phi,
        const  Real                   & a_dt)
{
  //compute kappa div^c F
  kappaConsDiv(a_phi);
  //compute nonconservative divergence = volume weighted ave of div^c
  nonConsDiv(a_phi);

  //advance solution, compute delta M
  DataIterator dit = m_grids.dataIterator();
  for(int ibox = 0; ibox < dit.size(); ibox++)
  {
    
  }

  //redistribute delta M
  redistribute();
    
}
#include "NamespaceFooter.H"

