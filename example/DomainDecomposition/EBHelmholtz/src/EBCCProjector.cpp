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
  auto & brit    = m_macprojector->m_brit;
  auto & doma    = m_macprojector->m_domain;
  for(int idir = 0; idir < DIM; idir++)
  {
    string ccGrad   = StencilNames::GradientCtoC  + std::to_string(idir);
    //need to set neumann bcs get gradients of pressure (matches no flow condition)
    brit->m_cellToCell->registerStencil(ccGrad  , StencilNames::Neumann, StencilNames::Neumann, doma, doma, false, Point::Zeroes());

    //dirichlet at domain to get zero normal velocity at domain boundaries
    //grown by one to allow interpolation to face centroids
    brit->registerCellToFace(StencilNames::AveCellToFace, StencilNames::Dirichlet, StencilNames::Neumann, doma, doma, false, Point::Ones());
  }

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
  auto & brit    = m_macprojector->m_brit;
  auto & doma    = m_macprojector->m_domain;
  auto & divu    = m_macprojector->m_rhs;
//  for(int idir = 0; idir < DIM; idir++)
//  {
//    string ccGrad   = StencilNames::GradientCtoC  + std::to_string(idir);
//    //need to set neumann bcs get gradients of pressure (matches no flow condition)
//    brit->m_cellToCell->registerStencil(ccGrad  , StencilNames::Neumann, StencilNames::Neumann, doma, doma, false, Point::Zeroes());
//    brit->m_faceToCell>registerStencil(StencilNames::AveCellToFace, StencilNames::Dirichlet, StencilNames::Neumann, doma, doma, false, Point::Ones());

  auto & grids   = m_macprojector->m_grids;
  auto & graphs  = m_macprojector->m_graphs;
  auto & nghost  = m_macprojector->m_nghost;
  DataIterator dit = grids.dataIterator();
  int ideb = 0;
  for(int ibox = 0; ibox < dit.size(); ++ibox)
  {

    Bx   grid   =  ProtoCh::getProtoBox(grids[dit[ibox]]);
    const EBGraph  & graph = (*graphs)[dit[ibox]];
    Bx  grown   =  grid.grow(ProtoCh::getPoint(nghost));

    //get face fluxes and interpolate them to centroids
    EBFluxData<Real, 1>         centroidv(grown, graph);
    EBFluxData<Real, 1>         facecentv(grown, graph);

    EBBoxData<CELL, Real, DIM>& inputvelo =  a_velo[dit[ibox]];
    //this one was registered by the mac projector
    EBFluxStencil<2, Real> interpsten  = brit->getFluxStencil(StencilNames::InterpToFaceCentroid, StencilNames::NoBC   , doma, doma, ibox);
    bool initZero = 0;
    for(unsigned int idir = 0; idir < DIM; idir++)
    {
      //true is inittozero
      EBBoxData<CELL, Real, 1> velcomp;
      velcomp.define(inputvelo, idir);
      brit->applyCellToFace(StencilNames::AveCellToFace, StencilNames::Neumann, 
                            doma,facecentv, velcomp, idir, ibox,  initZero, 1.0);
    }
    interpsten.apply(centroidv, facecentv, initZero, 1.0);  


    auto& kapdiv = divu[dit[ibox]];
    kapdiv.setVal(0.);
    for(unsigned int idir = 0; idir < DIM; idir++)
    {
      bool initToZero = false;
      //also registered by the mac projector
      brit->applyFaceToCell(StencilNames::DivergeFtoC, StencilNames::NoBC, doma, kapdiv, centroidv,
                            idir, ibox, initToZero, 1.0);
    }
    ideb++;
  }

//  static bool printed = false;
//  if(!printed)
//  {
//    printed = true;
//    writeEBLevelHDF5(string("divu.hdf5"), m_rhs, m_solver->getKappa(), m_domain, m_graphs);
//  }
}
///
#include "Chombo_NamespaceFooter.H"


