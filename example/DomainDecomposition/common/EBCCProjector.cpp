#include <cmath>
#include <memory>
#include "EBCCProjector.H"
#include "Chombo_ParmParse.H"
#include "Chombo_NamespaceHeader.H"
/// 
void
EBCCProjector::
define(shared_ptr<EBEncyclopedia<2, Real> >   & a_brit,     
       shared_ptr<GeometryService<2> >        & a_geoserv,  
       const DisjointBoxLayout                & a_grids,    
       const Box                              & a_domain,   
       const Real                             & a_dx,       
       const IntVect                          & a_nghost,
       const EBIBC                            & a_ebibc)

{
  CH_TIME("EBCCProjector::define");
  m_macprojector = shared_ptr<EBMACProjector>
    (new EBMACProjector(a_brit,     
                        a_geoserv,  
                        a_grids,    
                        a_domain,   
                        a_dx,       
                        a_nghost,
                        a_ebibc));
  registerStencils();
}

////
void  
EBCCProjector::
registerStencils()
{
  CH_TIME("EBCCProjector::registerStencils");
  auto & brit    = m_macprojector->m_brit;
  auto & doma    = m_macprojector->m_domain;
  for(int idir = 0; idir < DIM; idir++)
  {
    //dirichlet at domain to get zero normal velocity at domain boundaries
    //grown by one to allow interpolation to face centroids
    brit->registerCellToFace(StencilNames::AveCellToFace, StencilNames::Neumann, StencilNames::Neumann, doma, doma, false, Point::Ones(2));
    brit->registerFaceToCell(StencilNames::AveFaceToCell, StencilNames::NoBC   , StencilNames::NoBC   , doma, doma);
  }

}
/// 
void 
EBCCProjector::
project(EBLevelBoxData<CELL, DIM>   & a_velo,
        EBLevelBoxData<CELL, DIM>   & a_gphi,
        Real a_tol, unsigned int a_maxiter)
{
  CH_TIME("EBCCProjector::project");
  auto & rhs     = m_macprojector->m_rhs;
  auto & phi     = m_macprojector->m_phi;
  auto & solver  = m_macprojector->m_solver;
  auto & grids   = m_macprojector->m_grids;
  auto & graphs  = m_macprojector->m_graphs;
  
  a_velo.exchange(m_macprojector->m_exchangeCopier);
  
  // set rhs = kappa*div (vel)
  kappaDivU(rhs, a_velo);

  //solve kappa*lapl(phi) = kappa*div(vel)
  solver->solve(phi, rhs, a_tol, a_maxiter);
  
  //v := v - gphi
  DataIterator dit = grids.dataIterator();
  bool useConservativeGradient = false;
  ParmParse pp;
  pp.query("use_conservative_gradient", useConservativeGradient);
  if(useConservativeGradient)
  {
    pout() << "Using a **conservative discretization*** for the CC gradient." << endl;
  }
  else
  {
    pout() << "Using the ***average of MAC gradients*** for the CC gradient." << endl;
  }

  
  for(int ibox = 0; ibox < dit.size(); ++ibox)
  {
    const EBGraph  & graph = (*graphs)[dit[ibox]];
    Bx   grid   =  ProtoCh::getProtoBox(grids[dit[ibox]]);
    auto& phifab =    phi[dit[ibox]];
    auto& gphfab = a_gphi[dit[ibox]];
    auto& velfab = a_velo[dit[ibox]];
    if(useConservativeGradient)
    {
      computeConservativeGradient(gphfab, phifab, graph, grid, ibox);
    }
    else
    {
      computeAverageFaceToCellGradient(gphfab, phifab, graph, grid, ibox);
    }

    velfab -= gphfab;
  }
}
void
EBCCProjector::
computeAverageFaceToCellGradient(EBBoxData<CELL, Real, DIM> & a_gph,
                                 EBBoxData<CELL, Real,   1> & a_phi,
                                 const EBGraph              & a_graph,
                                 const Bx                   & a_grid,
                                 const unsigned int         & a_ibox)
{
  auto & nghost  = m_macprojector->m_nghost;
  auto & doma    = m_macprojector->m_domain;
  auto & brit    = m_macprojector->m_brit;
  Bx  grown   =  a_grid.grow(ProtoCh::getPoint(nghost));

  //get the mac gradient at face centers.
  EBFluxData<Real, 1>         facegrad(grown, a_graph);
  //gphi = grad(phi)
  for(unsigned int idir = 0; idir < DIM; idir++)
  {
    bool initZero = true;
    //registered by the mac projector
    brit->applyCellToFace(StencilNames::MACGradient, StencilNames::NoBC, doma,
                          facegrad, a_phi, idir, a_ibox, initZero, 1.0);
  }
  for(unsigned int idir = 0; idir < DIM; idir++)
  {
    bool initZero = true;
    EBBoxData<CELL, Real, 1> gradComp;
    gradComp.define(a_gph, idir);
    brit->applyFaceToCell(StencilNames::AveFaceToCell, StencilNames::NoBC, doma,
                          gradComp, facegrad,  idir, a_ibox, initZero, 1.0);
  }
  return;
}
void
EBCCProjector::
computeConservativeGradient(EBBoxData<CELL, Real, DIM> & a_gph,
                            EBBoxData<CELL, Real,   1> & a_phi,
                            const EBGraph              & a_graph,
                            const Bx                   & a_grid,
                            const unsigned int         & a_ibox)
{

  auto & nghost  = m_macprojector->m_nghost;
  auto & doma    = m_macprojector->m_domain;
  auto & brit    = m_macprojector->m_brit;
  Bx  grown   =  a_grid.grow(ProtoCh::getPoint(nghost));

  //get phi at face centers.
  EBFluxData<Real, 1>         phiFaceCent(grown, a_graph);
  brit->applyFaceToCell(StencilNames::AveFaceToCell, StencilNames::NoBC, doma,
                        a_phi, phiFaceCent,  a_ibox, true, 1.0);

  //interpolate phi to face centroids
  EBFluxData<Real, 1>         phiCentroid(grown, a_graph);
  //registered by the mac projector
  EBFluxStencil<2, Real> stencils =
      brit->getFluxStencil(StencilNames::InterpToFaceCentroid, StencilNames::NoBC, doma, doma, a_ibox);
  stencils.apply(phiCentroid, phiFaceCent, true, 1.0);  //true is to initialize to zero

  MayDay::Error("not finished");
  return;
}
///
void 
EBCCProjector::
kappaDivU(EBLevelBoxData<CELL, 1  > & a_divu,
          EBLevelBoxData<CELL, DIM> & a_velo)
{
  CH_TIME("EBCCProjector::kappaDivU");
  auto & brit    = m_macprojector->m_brit;
  auto & doma    = m_macprojector->m_domain;
  auto & divu    = m_macprojector->m_rhs;
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
    bool initZero = true;
    brit->applyCellToFace(StencilNames::AveCellToFace, StencilNames::Neumann, 
                          doma,facecentv, inputvelo,  ibox,  initZero, 1.0);

   //interpolate from face centers to centroids
    interpsten.apply(centroidv, facecentv, initZero, 1.0);  

    //brutally enforce flux-based boundary conditions on the hapless data
    m_macprojector->applyVeloBoundaryConditions(centroidv, dit[ibox]);
    auto& kapdiv = divu[dit[ibox]];
    //get kappa*divv
    //also registered by the mac projector
    brit->applyFaceToCell(StencilNames::DivergeFtoC, StencilNames::NoBC, doma, kapdiv, centroidv,
                          ibox, initZero, 1.0);
    ideb++;
  }

}

#include "Chombo_NamespaceFooter.H"


