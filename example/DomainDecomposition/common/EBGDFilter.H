#ifndef _EBGDFilter_H__
#define _EBGDFilter_H__
#include <cmath>
#include <memory>
#include "Proto.H"
#include "Chombo_EBDictionary.H"
#include "Chombo_EBChombo.H"
#include "Chombo_EBEncyclopedia.H"
#include "Chombo_EBLevelBoxData.H"
#include "Chombo_EBLevelFluxData.H"
#include "EBMultigrid.H"
#include "EBMACProjector.H"
#include "Chombo_NamespaceHeader.H"

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
class EBGDFilter
{
public:

  EBGDFilter(shared_ptr<EBEncyclopedia<2, Real> >   & a_brit,
             shared_ptr<GeometryService<2> >        & a_geoserv,
             const DisjointBoxLayout                & a_vecgrids,
             const Box                              & a_domain,
             const Real                             & a_dx,
             const IntVect                          & a_nghost,
             Real a_lambdaScale = 0.125)
  {
    define(a_brit,          
           a_geoserv,       
           a_vecgrids,      
           a_domain,        
           a_dx,            
           a_nghost,
           a_lambdaScale);
  }
  
  void define(shared_ptr<EBEncyclopedia<2, Real> >   & a_brit,
              shared_ptr<GeometryService<2> >        & a_geoserv,
              const DisjointBoxLayout                & a_vecgrids,
              const Box                              & a_domain,
              const Real                             & a_dx,
              const IntVect                          & a_nghost,
              Real a_lambdaScale = 0.125);
  
  /// 
  /**
     Input: velo holds input vector field.
     Output: velo holds filtered vector field.
   */
    void     filter(EBLevelBoxData<CELL, DIM>   & a_velo);



protected:

  void registerStencils();


private:

  ///insert usual diatribe about strong vs. weak construction here
  void operator=(const EBGDFilter& a_opin);
  EBGDFilter();

};
#include "Chombo_NamespaceFooter.H"

#endif
