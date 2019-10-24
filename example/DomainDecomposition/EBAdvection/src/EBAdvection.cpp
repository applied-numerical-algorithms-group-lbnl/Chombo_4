#include "EBAdvection.H"
#include "NamespaceHeader.H"
using Proto::Var;
////
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
////
EBAdvection::
EBAdvection(shared_ptr<EBEncyclopedia<2, Real> >   & a_brit,
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
  m_exchangeCopier.exchangeDefine(m_grids, m_nghostSrc);
  m_grids      = a_grids;      
  m_domain     = a_domain;     
  m_nghostSrc  = a_nghostsrc;
  m_nghostDst  = a_nghostdst;
  m_brit       = a_brit;
  m_veloCell   = a_veloCell;
  defineData(a_geoserv);
  fillKappa(a_geoserv);
  registerStencils();
}
////
void  
EBAdvection::
defineData(shared_ptr<GeometryService<2> >        & a_geoserv)
{
  const shared_ptr<LevelData<EBGraph>  > graphs = a_geoserv->getGraphs(m_domain);
  m_kappa.define(     m_grids, m_nghostSrc, graphs);
  m_kappaDiv.define(  m_grids, m_nghostSrc, graphs);
  m_deltaM.define(    m_grids, m_nghostSrc, graphs);
  m_nonConsDiv.define(m_grids, m_nghostSrc, graphs);
  m_hybridDiv.define( m_grids, m_nghostSrc, graphs);

  m_exchangeCopier.exchangeDefine(m_grids, m_nghostSrc);
}
////
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
  m_kappa.exchange(m_exchangeCopier);
}
////
void  
EBAdvection::
registerStencils()
{
  //volume weighted averaging radius one
  string nobcs("nobcs");
  string averaging("Volume_Weighted_Averaging_rad_1");
  string redistribution("Volume_Weighted_Redistribution_rad_1");
  //false is because I do not need diagonal  weights for any of these stencils
  bool needDiag = false;
  m_brit->m_cellToCell->registerStencil(averaging     , nobcs, nobcs, m_domain, m_domain, needDiag);
  m_brit->m_cellToCell->registerStencil(redistribution, nobcs, nobcs, m_domain, m_domain, needDiag);

  //number of potential cells for redistribution

  for(int idir = 0; idir < DIM; idir++)
  {
    string centInterp("InterpolateToFaceCentroid_");
    string   slopeLow("Slope_Low_");
    string  slopeHigh("Slope_High_");
    centInterp += std::to_string(idir);
    slopeLow   += std::to_string(idir);
    slopeHigh  += std::to_string(idir);
    m_brit->registerFaceStencil(          idir, centInterp, nobcs, nobcs, m_domain, m_domain, needDiag);
    m_brit->m_cellToCell->registerStencil(      slopeLow  , nobcs, nobcs, m_domain, m_domain, needDiag);
    m_brit->m_cellToCell->registerStencil(      slopeHigh , nobcs, nobcs, m_domain, m_domain, needDiag);
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
//    Bx grbx = ProtoCh::getProtoBox(m_grids[dit[ibox]]);
//    ebforallIrreg("HybridDivergence", HybridDivergence, grbx,  
//                  m_hybridDiv[dit[ibox]], 
//                  m_kappaDiv[dit[ibox]], m_nonConsDiv[dit[ibox]], 
//                  m_deltaM[dit[ibox]], m_kappa[dit[ibox]]);
  }

  //redistribute delta M
  redistribute();
    
}
#include "NamespaceFooter.H"

