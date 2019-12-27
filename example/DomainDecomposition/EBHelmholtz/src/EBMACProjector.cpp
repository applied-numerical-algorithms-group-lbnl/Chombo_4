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
  m_dx      = a_dx;

  m_grids   = a_grids;
  m_domain  = a_domain;
  m_nghost  = a_nghost;
  m_grids   = a_grids;      
  m_domain  = a_domain;     
  m_brit    = a_brit;
  m_graphs  = a_geoserv->getGraphs(m_domain);

  defineInternals(a_geoserv);
}
///
void 
EBMACProjector::
defineInternals(shared_ptr<GeometryService<2> >        & a_geoserv)

{
  m_exchangeCopier.exchangeDefine(m_grids, m_nghost);
  m_rhs.define(m_grids, m_nghost, m_graphs);
  m_phi.define(m_grids, m_nghost, m_graphs);
  auto ditch = m_brit->m_cellToCell;
  Real alpha = 0; Real beta = 1; //Poisson's eqn

  m_solver = shared_ptr<EBMultigrid>
    (new EBMultigrid(ditch, a_geoserv, alpha, beta, m_dx, m_grids, 
                     StencilNames::Poisson2, StencilNames::Neumann, StencilNames::Neumann,
                     m_domain, m_nghost, m_nghost));

  registerStencils();
}
////
void  
EBMACProjector::
registerStencils()
{

  //false is because I do not need diagonal  weights for any of these stencils
  bool needDiag = false;
  
  //interpolates from face centers to centroids
  m_brit->registerFaceStencil(StencilNames::InterpToFaceCentroid, StencilNames::NoBC,    StencilNames::NoBC, m_domain, m_domain, needDiag);

  //increment divergence by face difference
  m_brit->registerFaceToCell( StencilNames::DivergeFtoC         , StencilNames::NoBC,    StencilNames::NoBC, m_domain, m_domain, needDiag);

  //face-centered gradient of cell-centered data
  m_brit->registerCellToFace( StencilNames::MACGradient         , StencilNames::Neumann, StencilNames::NoBC, m_domain, m_domain, needDiag);
}
/// 
void 
EBMACProjector::
project(EBLevelFluxData<1>   & a_velo,
        EBLevelFluxData<1>   & a_gphi,
        Real a_tol, unsigned int a_maxiter)
{
  // set rhs = div (vel)
  divergence(m_rhs, a_velo);

  //solve lapl(phi) = rhs
  m_solver->solve(m_phi, m_rhs, a_tol, a_maxiter);
  
  //gphi = grad(phi)
  //v := v - gphi
  DataIterator dit = m_grids.dataIterator();
  int ideb = 0;
  for(int ibox = 0; ibox < dit.size(); ibox++)
  {
    
    Bx   grid   =  ProtoCh::getProtoBox(m_grids[dit[ibox]]);
    //get face fluxes and interpolate them to centroids
    for(unsigned int idir = 0; idir < DIM; idir++)
    {
      bool initToZero = true;
      m_brit->applyCellToFace(StencilNames::MACGradient, StencilNames::NoBC, m_domain,
                              a_gphi[dit[ibox]] ,m_phi[dit[ibox]], idir, ibox, initToZero, 1.0);
    }
    a_velo[dit[ibox]] -= a_gphi[dit[ibox]];
    ideb++;
  }
}
///
void 
EBMACProjector::
divergence(EBLevelBoxData<CELL, 1> & a_divu,
           EBLevelFluxData<1>      & a_velo)
{

  DataIterator dit = m_grids.dataIterator();
  int ideb = 0;
  for(int ibox = 0; ibox < dit.size(); ++ibox)
  {

    Bx   grid   =  ProtoCh::getProtoBox(m_grids[dit[ibox]]);
    const EBGraph  & graph = (*m_graphs)[dit[ibox]];

    //get face fluxes and interpolate them to centroids
    EBFluxData<Real, 1>  centroidFlux(grid, graph);
    EBFluxStencil<2, Real> stencils =
      m_brit->getFluxStencil(StencilNames::InterpToFaceCentroid, StencilNames::NoBC, m_domain, m_domain, ibox);
    EBFluxData<Real,1>& faceCentFlux = a_velo[dit[ibox]];
    stencils.apply(centroidFlux, faceCentFlux, true, 1.0);  //true is to initialize to zero


    auto& kapdiv =  m_rhs[dit[ibox]];
    kapdiv.setVal(0.);
    for(unsigned int idir = 0; idir < DIM; idir++)
    {
      bool initToZero = false;
      m_brit->applyFaceToCell(StencilNames::DivergeFtoC, StencilNames::NoBC, m_domain, kapdiv, centroidFlux,
                              idir, ibox, initToZero, 1.0);
    }
    ideb++;
  }
}
///
EBLevelBoxData<CELL, 1>& 
EBMACProjector::
getRHSHolder()
{
  return m_rhs;
}
///
#include "Chombo_NamespaceFooter.H"


