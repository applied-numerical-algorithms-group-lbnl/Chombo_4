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
       const EBIBC                            & a_ebibc,
       bool a_printStuff)

{
  CH_TIME("EBCCProjector::define");
  if(a_printStuff)
  {
    pout() << "EBCCProjector::define calling EBMACProjector constructor" << endl;
  }
  m_macprojector = shared_ptr<EBMACProjector>
    (new EBMACProjector(a_brit,     
                        a_geoserv,  
                        a_grids,    
                        a_domain,   
                        a_dx,       
                        a_nghost,
                        a_ebibc,
                        a_printStuff));

  //auto & doma    = m_macprojector->m_domain;
  auto & graphs  = m_macprojector->m_graphs;
  auto & grids   = m_macprojector->m_grids;
  m_macGradient.define(grids, a_nghost, graphs);

  if(a_printStuff)
  {
    pout() << "EBCCProjector::define calling register stencils" << endl;
  }
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
  //dirichlet at domain to get zero normal velocity at domain boundaries
  //grown by one to allow interpolation to face centroids
  brit->registerCellToFace(StencilNames::AveCellToFace, StencilNames::Neumann, StencilNames::Neumann, doma, doma, false, Point::Ones(2));
  brit->registerFaceToCell(StencilNames::AveFaceToCell, StencilNames::NoBC   , StencilNames::NoBC   , doma, doma);

  //below here is stuff used in the conservative gradient
  Point ghost = Point::Zeroes();
  bool needDiag = false;
  m_nobcsLabel = StencilNames::NoBC;
  m_ncdivLabel     = StencilNames::NCDivergeRoot + string("1"); //this is for the normalizor 
  brit->m_cellToCell->registerStencil(m_ncdivLabel , m_nobcsLabel, m_nobcsLabel, doma, doma, needDiag);
  brit->m_cellToBoundary->registerStencil(StencilNames::CopyCellValueToCutFace,
                                          m_nobcsLabel,  m_nobcsLabel, doma, doma,
                                          needDiag, ghost);
  for(unsigned int idir = 0; idir < DIM; idir++)
  {
    m_conGradEBFluxName[idir] = StencilNames::CutFaceIncrementToKappaDiv + std::to_string(idir);
    brit->m_boundaryToCell->registerStencil(m_conGradEBFluxName[idir],
                                            m_nobcsLabel, m_nobcsLabel, doma, doma,
                                            needDiag, ghost);
  }
  //extrapolates face values to domain boundaries
  brit->registerFaceStencil(StencilNames::ExtrapToDomainFace,
                            m_nobcsLabel, m_nobcsLabel,  doma, doma, needDiag);
}
///
void
EBCCProjector::
printVectorMax(const EBLevelBoxData<CELL, DIM>& a_velo,
               const string                   & a_prefix) const
{
  Chombo4::pout() << a_prefix << " = (";
  for(int idir = 0; idir < DIM; idir++)
  {
    Real maxvel;
    EBIndex<CELL> maxvof;
    maxvel = a_velo.maxNorm(maxvof, idir);
    Chombo4::pout() << maxvel << "@" << maxvof.m_pt ;
    if(idir < (DIM-1))
    {
      Chombo4:: pout() << ",";
    }
  }
  Chombo4::pout() << a_prefix << ")"  << std::endl;
}
/// 
void 
EBCCProjector::
project(EBLevelBoxData<CELL, DIM>   & a_velo,
        EBLevelBoxData<CELL, DIM>   & a_gphi,
        Real a_tol, unsigned int a_maxiter,
        bool a_printStuff)
{
  CH_TIME("EBCCProjector::project");
  auto & rhs     = m_macprojector->m_rhs;
  auto & phi     = m_macprojector->m_phi;
  auto & solver  = m_macprojector->m_solver;
  auto & grids   = m_macprojector->m_grids;
  
  a_velo.exchange();
  if(a_printStuff)
  {
    printVectorMax(a_velo, string("EBCCProjector::project: |input velo|max"));;
  }
  
  // set rhs = kappa*div (vel)
  kappaDivU(rhs, a_velo, a_printStuff);

  if(a_printStuff)
  {
    EBIndex<CELL> vofmax;
    Real kapDivUMax = rhs.maxNorm(vofmax, 0);
    Chombo4::pout() << "EBCCProjector::project: |kappa DiVU|max  = " << kapDivUMax << "@" << vofmax.m_pt << endl;
  }
    
    
  //solve kappa*lapl(phi) = kappa*div(vel)
  solver->solve(phi, rhs, a_tol, a_maxiter);
  phi.exchange();
  if(a_printStuff)
  {
    EBIndex<CELL> vofmax;
    Real phiMax = phi.maxNorm(vofmax, 0);
    Chombo4::pout() << "EBCCProjector::project: |phi|max  = " << phiMax << "@" << vofmax.m_pt << endl;
  }
    
  //v := v - gphi
  DataIterator dit = grids.dataIterator();
  bool useConservativeGradient = false;
  ParmParse pp;
  pp.query("use_conservative_gradient", useConservativeGradient);
  if(useConservativeGradient)
  {
    pout() << "Using a **conservative discretization*** for the CC gradient." << endl;
    kappaConservativeGradient(a_gphi, phi);
  }
  else
  {
    pout() << "Using the ***average of MAC gradients*** for the CC gradient." << endl;
    averageFaceToCellGradient(a_gphi, phi);
  }
  
  if(a_printStuff)
  {
    printVectorMax(a_gphi, string("EBLevelCCProjector::project:|gph|max"));
  }
  
  //conservative gradient really produces kappa*grad phi so
  //we have to normalize
  if(useConservativeGradient)
  {
    if(a_printStuff)
    {
      printVectorMax(a_gphi, string("EBLevelCCProjector::project:|gph|max before normalization"));
    }
    
    normalizeGradient(a_gphi);

    if(a_printStuff)
    {
      printVectorMax(a_gphi, string("EBLevelCCProjector::project:|gph|max after normalization"));
    }
  }

  for(int ibox = 0; ibox < dit.size(); ++ibox)
  {
    auto& gphfab = a_gphi[dit[ibox]];
    auto& velfab = a_velo[dit[ibox]];
    velfab -= gphfab;
  }
}

///
void EBCCProjector::
normalizeGradient(EBLevelBoxData<CELL, DIM>& a_gphi)
{
  CH_TIME("EBCCProjector::normalizeGradient");
  auto & grids          = m_macprojector->m_grids;
  auto & graphs         = m_macprojector->m_graphs;
  auto & nghost         = m_macprojector->m_nghost;
  auto & doma    = m_macprojector->m_domain;
  auto & brit    = m_macprojector->m_brit;
  
  EBLevelBoxData<CELL, DIM>  kappaGrad(grids, nghost, graphs);
  //right now gphi holds kappa * grad
  a_gphi.copyTo(kappaGrad);
  
  //this makes ncdiv = divF on regular cells
  //and ncdiv = vol_weighted_ave(div) on cut cells
  kappaGrad.exchange();
  
  DataIterator dit = grids.dataIterator();
  for(unsigned int idir = 0; idir < DIM; idir++)
  {
    EBLevelBoxData<CELL, 1> kappComp, gradComp;
    gradComp.define<DIM>(a_gphi   , idir, graphs);
    kappComp.define<DIM>(kappaGrad, idir, graphs);
    for(int ibox = 0; ibox < dit.size(); ++ibox)
    {
      auto& ncdiv  = gradComp[dit[ibox]];
      auto& kapdiv = kappComp[dit[ibox]];

      const auto& stencil = brit->m_cellToCell->getEBStencil(m_ncdivLabel,m_nobcsLabel,
                                                             doma, doma, ibox);
      bool initToZero = true;
      stencil->apply(ncdiv, kapdiv, initToZero, 1.0);
    }
  }
}
///
void
EBCCProjector::
averageFaceToCellGradient(EBLevelBoxData<CELL, DIM>  & a_gph,
                          EBLevelBoxData<CELL,   1>  & a_phi)
{
  auto & doma    = m_macprojector->m_domain;
  auto & grids   = m_macprojector->m_grids;
  auto & brit    = m_macprojector->m_brit;
  
  a_phi.exchange();
  DataIterator dit = grids.dataIterator();

  //first get mac gradient
  for(int ibox = 0; ibox < dit.size(); ++ibox)
  {
    Bx   grid   =  ProtoCh::getProtoBox(grids[dit[ibox]]);
    auto& phifab = a_phi[dit[ibox]];
    //get the mac gradient at face centers.
    auto& facegrad = m_macGradient[dit[ibox]];
    //gphi = grad(phi)
    for(unsigned int idir = 0; idir < DIM; idir++)
    {
      bool initZero = true;
      //registered by the mac projector
      brit->applyCellToFace(StencilNames::MACGradient, StencilNames::NoBC, doma,
                            facegrad, phifab, idir, ibox, initZero, 1.0);
    }
  }
  //exchange mac gradient
  m_macGradient.exchange();
  
  for(int ibox = 0; ibox < dit.size(); ++ibox)
  {
    Bx   grid   =  ProtoCh::getProtoBox(grids[dit[ibox]]);

    auto& gphfab   = a_gph[dit[ibox]];
    auto& facegrad = m_macGradient[dit[ibox]];
    //gphi = grad(phi)
    for(unsigned int idir = 0; idir < DIM; idir++)
    {
      bool initZero = true;
      EBBoxData<CELL, Real, 1> gradComp;
      ///aliasing define to get correct component of the cell centered gradient
      gradComp.define(gphfab, idir);
      brit->applyFaceToCell(StencilNames::AveFaceToCell, StencilNames::NoBC, doma,
                            gradComp, facegrad,  idir, ibox, initZero, 1.0);
    }
  }
}
///
void
EBCCProjector::
averageFaceToCellGradient(EBBoxData<CELL, Real, DIM> & a_gph,
                          EBBoxData<CELL, Real,   1> & a_phi,
                          const EBGraph              & a_graph,
                          const Bx                   & a_grid,
                          const unsigned int         & a_ibox)
{
  bool useStack = true;
  CH_TIME("EBCCProjector::avefacetocell_gradient");
  auto & nghost  = m_macprojector->m_nghost;
  auto & doma    = m_macprojector->m_domain;
  auto & brit    = m_macprojector->m_brit;
  Bx  grown   =  a_grid.grow(ProtoCh::getPoint(nghost));

  //get the mac gradient at face centers.
  EBFluxData<Real, 1>         facegrad(grown, a_graph, useStack);
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
///
void
EBCCProjector::
kappaConservativeGradient(EBLevelBoxData<CELL, DIM>  & a_kappaGrad,
                          EBLevelBoxData<CELL,   1>  & a_phi)
{
  auto & grids   = m_macprojector->m_grids;
  auto & graphs  = m_macprojector->m_graphs;
  
  a_phi.exchange();
  DataIterator dit = grids.dataIterator();
  for(int ibox = 0; ibox < dit.size(); ++ibox)
  {
    const auto & graph = (*graphs)[dit[ibox]];
    Bx   grid   =  ProtoCh::getProtoBox(grids[dit[ibox]]);
    auto& phifab = a_phi[dit[ibox]];
    auto& gphfab = a_kappaGrad[dit[ibox]];
    kappaConservativeGradient(gphfab, phifab, graph, grid, ibox);
  }
}
///
void
EBCCProjector::
kappaConservativeGradient(EBBoxData<CELL, Real, DIM> & a_kappaGrad,
                          EBBoxData<CELL, Real,   1> & a_phi,
                          const EBGraph              & a_graph,
                          const Bx                   & a_grid,
                          const unsigned int         & a_ibox)
{
  CH_TIME("EBCCProjector::kappaConservativeGradient");
  auto & nghost  = m_macprojector->m_nghost;
  auto & doma    = m_macprojector->m_domain;
  auto & brit    = m_macprojector->m_brit;
  auto & grids   = m_macprojector->m_grids;
  Bx  grown   =  a_grid.grow(ProtoCh::getPoint(nghost));

  //get phi at face centers.
  EBFluxData<Real, 1>         phiFaceCent(grown, a_graph);
  brit->applyCellToFace(StencilNames::AveCellToFace, StencilNames::Neumann, 
                        doma, phiFaceCent, a_phi,  a_ibox,  true, 1.0);

  //interpolate phi to face centroids
  EBFluxData<Real, 1>         phiCentroid(grown, a_graph);
  //registered by the mac projector
  EBFluxStencil<2, Real> centroidStencils =
    brit->getFluxStencil(StencilNames::InterpToFaceCentroid, m_nobcsLabel, doma, doma, a_ibox);
  EBFluxStencil<2, Real> domaExtrStencils =
    brit->getFluxStencil(StencilNames::ExtrapToDomainFace  , m_nobcsLabel,  doma, doma, a_ibox);

  //copy that and extrapolate to domain boundary faces
  domaExtrStencils.apply(phiFaceCent, phiFaceCent, true, 1.0);  
  //true is to initialize to zero
  //get phi at centroids
  centroidStencils.apply(phiCentroid, phiFaceCent, true, 1.0);  
  
  EBBoxData<BOUNDARY, Real, 1> ebflux(grown, a_graph);
  //these two are used in the conservative gradient
  auto copystenptr  = brit->m_cellToBoundary->getEBStencil(StencilNames::CopyCellValueToCutFace,
                                                           StencilNames::NoBC, doma, doma, a_ibox);
  copystenptr->apply(ebflux, a_phi, true, 1.0);

  int ideb = 0;
  for(unsigned int idir = 0; idir < DIM; idir++)
  {
    EBBoxData<CELL, Real, 1> gradComp;
    gradComp.define(a_kappaGrad, idir);
    gradComp.setVal(0.); //because this is done via increments
    
    DataIterator dit = grids.dataIterator();
    for(unsigned int ibox = 0; ibox < dit.size(); ibox++)
    {
      //stencil for eb contribution to conservative gradient in this direction
      auto ebdivstenptr =
        brit->m_boundaryToCell->getEBStencil(m_conGradEBFluxName[idir],
                                             StencilNames::NoBC, doma, doma, a_ibox);;

      bool initToZero = false;
      //stencil for face centroid fluxes is registered in the mac projector
      //this not a nested loop because only the normal faces matter here.
      brit->applyFaceToCell(StencilNames::DivergeFtoC, StencilNames::NoBC,
                            doma, gradComp, phiCentroid,
                            idir, ibox, initToZero, 1.0);

      ebdivstenptr->apply(gradComp, ebflux, initToZero, 1.0);  
      ideb++;
    }
  }
  

  return;
}
///
void 
EBCCProjector::
kappaDivU(EBLevelBoxData<CELL, 1  > & a_divu,
          EBLevelBoxData<CELL, DIM> & a_velo,
          bool a_printStuff)
{
  bool useStack = true;
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
    EBFluxData<Real, 1>         centroidv(grown, graph, useStack);
    EBFluxData<Real, 1>         facecentv(grown, graph, useStack);

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


