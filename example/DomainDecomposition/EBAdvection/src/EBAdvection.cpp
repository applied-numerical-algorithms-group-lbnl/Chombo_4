#include "EBAdvection.H"
#include "NamespaceHeader.H"
const string EBAdvection::s_ncdivLabel     = string("Volume_Weighted_Averaging_rad_1"); //this is for the non-conservative div
const string EBAdvection::s_nobcsLabel     = string("no_bcs"); //none of the operators here have eb boundary conditions
const string EBAdvection::s_redistLabel    = string("Volume_Weighted_Redistribution_rad_1"); //for redistribution
const string EBAdvection::s_centInterpLabel= string("InterpolateToFaceCentroid_");
const string EBAdvection::s_slopeLowLabel  = string("Slope_Low_");
const string EBAdvection::s_slopeHighLabel = string("Slope_High_");
using Proto::Var;
////
PROTO_KERNEL_START 
void HybridDivergenceF(Var<Real, 1>    a_hybridDiv,
                       Var<Real, 1>    a_nonConsDivF,
                       Var<Real, 1>    a_deltaM,
                       Var<Real, 1>    a_kappa)
{
  //hybrid divergence comes in holding kappa*consDiv
  Real kappaConsDiv = a_hybridDiv(0);
  a_hybridDiv(0) = kappaConsDiv + (1- a_kappa(0))*a_nonConsDivF(0);
  a_deltaM(0)    = (1-a_kappa(0))*(kappaConsDiv - a_kappa(0)*a_nonConsDivF(0));
}
PROTO_KERNEL_END(HybridDivergenceF, HybridDivergence)
////
PROTO_KERNEL_START 
void AdvanceScalarF(Var<Real, 1>    a_scal,
                    Var<Real, 1>    a_divF,
                    Real            a_dt)
{
  a_scal(0) = a_scal(0) - a_dt*a_divF(0);
}
PROTO_KERNEL_END(AdvanceScalarF, AdvanceScalar)
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
  fillKappa( a_geoserv);
  registerStencils();
}
////
void  
EBAdvection::
defineData(shared_ptr<GeometryService<2> >        & a_geoserv)
{
  m_graphs = a_geoserv->getGraphs(m_domain);
  m_kappa.define(     m_grids, m_nghostSrc, m_graphs);
  m_deltaM.define(    m_grids, m_nghostSrc, m_graphs);
  m_nonConsDiv.define(m_grids, m_nghostSrc, m_graphs);
  m_hybridDiv.define( m_grids, m_nghostSrc, m_graphs);

  m_exchangeCopier.exchangeDefine(m_grids, m_nghostSrc);
}
////
void  
EBAdvection::
fillKappa(shared_ptr<GeometryService<2> >        & a_geoserv)
{
  DataIterator dit = m_grids.dataIterator();
  for(int ibox = 0; ibox < dit.size(); ++ibox)
  {
    Box grid =m_grids[dit[ibox]];
    Bx  grbx = getProtoBox(grid);
    const EBGraph  & graph = (*m_graphs)[dit[ibox]];
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

  //false is because I do not need diagonal  weights for any of these stencils
  bool needDiag = false;
  m_brit->m_cellToCell->registerStencil(s_ncdivLabel , s_nobcsLabel, s_nobcsLabel, m_domain, m_domain, needDiag);
  m_brit->m_cellToCell->registerStencil(s_redistLabel, s_nobcsLabel, s_nobcsLabel, m_domain, m_domain, needDiag);

  //number of potential cells for redistribution

  m_brit->registerFaceStencil(s_centInterpLabel, s_nobcsLabel, s_nobcsLabel, m_domain, m_domain, needDiag);
  for(int idir = 0; idir < DIM; idir++)
  {
    string slopeLow   = s_slopeLowLabel  + std::to_string(idir);
    string slopeHigh  = s_slopeHighLabel + std::to_string(idir);
    m_brit->m_cellToCell->registerStencil(slopeLow  , s_nobcsLabel, s_nobcsLabel, m_domain, m_domain, needDiag);
    m_brit->m_cellToCell->registerStencil(slopeHigh , s_nobcsLabel, s_nobcsLabel, m_domain, m_domain, needDiag);
  }

}
///
void
EBAdvection::
kappaConsDiv(EBLevelBoxData<CELL, 1>   & a_scal)
{
  //coming into this we have the scalar at time = n dt
  // velocity field at cell centers. Leaving, we have filled
  // kappa* div(u scal)
  DataIterator dit = m_grids.dataIterator();
  for(int ibox = 0; ibox < dit.size(); ++ibox)
  {
    Bx   grid   =  ProtoCh::getProtoBox(m_grids[dit[ibox]]);
    Bx  grown   =  grid.grow(1);
    const EBGraph  & graph = (*m_graphs)[dit[ibox]];
    //get face fluxes and interpolate them to centroids
    EBFluxData<Real, 1> centroidFlux(grid , graph);
    EBFluxData<Real, 1>   centerFlux(grown, graph);
    EBFluxStencil<2, Real> stencils =   m_brit->getFluxStencil(s_centInterpLabel, s_nobcsLabel, m_domain, m_domain, ibox);

    // HERE auto& kapdiv =  m_hybridDiv[dit[ibox]];
  }

}
///
void
EBAdvection::
nonConsDiv()
{
  //hybrid div comes in holding kappa*cons_div(F)
  //this makes ncdiv = divF on regular cells
  //and ncdiv = vol_weighted_ave(div) on cut cells
  DataIterator dit = m_grids.dataIterator();
  for(int ibox = 0; ibox < dit.size(); ++ibox)
  {
    auto& ncdiv  = m_nonConsDiv[dit[ibox]];
    auto& kapdiv =  m_hybridDiv[dit[ibox]];
    const auto& stencil = m_brit->m_cellToCell->getEBStencil(s_ncdivLabel,s_nobcsLabel, m_domain, m_domain, ibox);
    bool initToZero = true;
    stencil->apply(ncdiv, kapdiv, initToZero, 1.0);
  }
}
///
void
EBAdvection::
redistribute()
{
  //hybrid div comes in holding kappa*div^c + (1-kappa)div^nc
  //this redistributes delta M into the hybrid divergence.
  DataIterator dit = m_grids.dataIterator();
  for(int ibox = 0; ibox < dit.size(); ++ibox)
  {
    auto& deltaM = m_deltaM   [dit[ibox]];
    auto& hybrid = m_hybridDiv[dit[ibox]];
    const auto& stencil = m_brit->m_cellToCell->getEBStencil(s_redistLabel,s_nobcsLabel, m_domain, m_domain, ibox);
    bool initToZero = false;
    stencil->apply(hybrid, deltaM, initToZero, 1.0);
  }
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
  nonConsDiv();

  //advance solution, compute delta M
  DataIterator dit = m_grids.dataIterator();
  for(int ibox = 0; ibox < dit.size(); ibox++)
  {
    unsigned long long int numflopspt = 7;
    Bx grbx = ProtoCh::getProtoBox(m_grids[dit[ibox]]);
    ebforallInPlace(numflopspt, "HybridDivergence", HybridDivergence, grbx,  
                    m_hybridDiv[dit[ibox]], m_nonConsDiv[dit[ibox]],  
                    m_deltaM[dit[ibox]], m_kappa[dit[ibox]]);
  }
  m_deltaM.exchange(m_exchangeCopier);
  //redistribute delta M
  redistribute();
  for(int ibox = 0; ibox < dit.size(); ibox++)
  {
    Bx grbx = ProtoCh::getProtoBox(m_grids[dit[ibox]]);
    unsigned long long int numflopspt = 2;
    ebforallInPlace(numflopspt, "AdvanceScalar", AdvanceScalar,  grbx,  
                    a_phi[dit[ibox]], m_hybridDiv[dit[ibox]],  a_dt);
  }
    
}
#include "NamespaceFooter.H"

