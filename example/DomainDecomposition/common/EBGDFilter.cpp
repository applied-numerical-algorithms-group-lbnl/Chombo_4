#include <cmath>
#include <memory>
#include "EBGDFilter.H"
#include "Chombo_NamespaceHeader.H"
/// 
void
EBGDFilter::
define(shared_ptr<EBEncyclopedia<2, Real> >   & a_brit,     
       shared_ptr<GeometryService<2> >        & a_geoserv,  
       const DisjointBoxLayout                & a_grids,    
       const Box                              & a_domain,   
       const Real                             & a_dx,       
       const IntVect                          & a_nghost,
       Real a_lambdaScale)

{
  registerStencils();
}

////
void  
EBGDFilter::
registerStencils()
{
  /**
  CH_TIME("EBGDFilter::registerStencils");
  auto & brit    = m_macprojector->m_brit;
  auto & doma    = m_macprojector->m_domain;
  for(int idir = 0; idir < DIM; idir++)
  {
    //dirichlet at domain to get zero normal velocity at domain boundaries
    //grown by one to allow interpolation to face centroids
    brit->registerCellToFace(StencilNames::AveCellToFace, StencilNames::Dirichlet, StencilNames::Neumann, doma, doma, false, Point::Ones());
    brit->registerFaceToCell(StencilNames::AveFaceToCell, StencilNames::NoBC     , StencilNames::NoBC   , doma, doma, false);
  }
  **/

}
/// 
///Phil's filter.   Helps with long term stability in lots of cases.
/**
   This filter helps remove unphysical velocity modes that are numerically 
   divergence-free according  to the cell centered divergence operator.  Most of these modes 
   are very high frequency variations on a checkerboard. 
   We create a face-centered divergence operator and a cell-centered gradient.
   If we were doing a full projection, the operator would look like
   F(v) = (I - G (DG)^-1 D) (v)  
   That is too expensive since we only need to kill very high frequency modes.
   Instead of doing the full solve (the DG-1 bit), we just do one step of Jacoby relaxation.
      Step 1.  Set phi=0, rhs = D(v). 
      Step 2.  jacoby: phi^l = phi^l-1 + lambda(rhs- DG(phi^l)).   
          2a,  Since phi^0 = 0,   phi^1 = lambda*rhs = lambda*D(v)
      Step 3.  Stop the solve (since Jacoby is pretty good at removing high frequency error)
      Step 4.  v := v - G(phi) = v - G(lambda Div(v))

    Lambda is the Jacoby relaxtion coefficient.
 */
void 
EBGDFilter::
filter(EBLevelBoxData<CELL, DIM>   & a_velo)
{
  CH_TIME("EBGDFilter::filter");
}

///
#include "Chombo_NamespaceFooter.H"


