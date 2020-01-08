#include <cmath>
#include <memory>
#include "EBCCProjector.H"
#include "Chombo_NamespaceHeader.H"
/// 
EBCCProjector::
EBCCProjector(shared_ptr<EBEncyclopedia<2, Real> >   & a_brit,     
              shared_ptr<GeometryService<2> >        & a_geoserv,  
              const DisjointBoxLayout                & a_grids,    
              const Box                              & a_domain,   
              const Real                             & a_dx,       
              const IntVect                          & a_nghost)   
{
  m_macprojector = shared_ptr<EBMACProjector>
    (new EBMACProjector(a_brit,     
                        a_geoserv,  
                        a_grids,    
                        a_domain,   
                        a_dx,       
                        a_nghost));
  registerStencils();
}
////
void  
EBCCProjector::
registerStencils()
{

}
/// 
void 
EBCCProjector::
project(EBLevelBoxData<CELL, DIM>   & a_velo,
        EBLevelBoxData<CELL, DIM>   & a_gphi,
        Real a_tol, unsigned int a_maxiter)
{
  auto & rhs     = m_macprojector->m_rhs;
  auto & brit    = m_macprojector->m_brit;
  auto & phi     = m_macprojector->m_phi;
  auto & solver  = m_macprojector->m_solver;
  // set rhs = div (vel)
  divergence(rhs, a_velo);

  //solve lapl(phi) = rhs
  solver->solve(phi, rhs, a_tol, a_maxiter);
  
  //gphi = grad(phi)
  //v := v - gphi
//  DataIterator dit = m_grids.dataIterator();
//  int ideb = 0;
//  for(int ibox = 0; ibox < dit.size(); ibox++)
//  {
//    
//    Bx   grid   =  ProtoCh::getProtoBox(m_grids[dit[ibox]]);
//    //get face fluxes and interpolate them to centroids
//    for(unsigned int idir = 0; idir < DIM; idir++)
//    {
//      bool initToZero = true;
//      brit->applyCellToFace(StencilNames::MACGradient, StencilNames::NoBC, m_domain,
//                              a_gphi[dit[ibox]] ,m_phi[dit[ibox]], idir, ibox, initToZero, 1.0);
//    }
//    a_velo[dit[ibox]] -= a_gphi[dit[ibox]];
//    ideb++;
//  }
}
///
void 
EBCCProjector::
divergence(EBLevelBoxData<CELL, 1  > & a_divu,
           EBLevelBoxData<CELL, DIM> & a_velo)
{

//  DataIterator dit = m_grids.dataIterator();
//  int ideb = 0;
//  for(int ibox = 0; ibox < dit.size(); ++ibox)
//  {
//
//    Bx   grid   =  ProtoCh::getProtoBox(m_grids[dit[ibox]]);
//    const EBGraph  & graph = (*m_graphs)[dit[ibox]];
//    Bx  grown   =  grid.grow(ProtoCh::getPoint(m_nghost));
//
//    //get face fluxes and interpolate them to centroids
//    EBFluxData<Real, 1>  centroidFlux(grown, graph);
//    EBFluxStencil<2, Real> stencils =
//      m_brit->getFluxStencil(StencilNames::InterpToFaceCentroid, StencilNames::NoBC, m_domain, m_domain, ibox);
//    EBFluxData<Real,1>& faceCentFlux = a_velo[dit[ibox]];
//    stencils.apply(centroidFlux, faceCentFlux, true, 1.0);  //true is to initialize to zero
//
//
//    auto& kapdiv =  m_rhs[dit[ibox]];
//    kapdiv.setVal(0.);
//    for(unsigned int idir = 0; idir < DIM; idir++)
//    {
//      bool initToZero = false;
//      m_brit->applyFaceToCell(StencilNames::DivergeFtoC, StencilNames::NoBC, m_domain, kapdiv, centroidFlux,
//                              idir, ibox, initToZero, 1.0);
//    }
//    ideb++;
//  }
//
//  static bool printed = false;
//  if(!printed)
//  {
//    printed = true;
//    writeEBLevelHDF5(string("divu.hdf5"), m_rhs, m_solver->getKappa(), m_domain, m_graphs);
//  }
}
///
#include "Chombo_NamespaceFooter.H"


