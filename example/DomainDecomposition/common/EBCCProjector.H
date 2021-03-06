#ifndef _EBCCProjector_H__
#define _EBCCProjector_H__
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

///class to approximately project cell centered velocity fields an eb context (as in Trebotich and Graves)
class EBCCProjector
{
public:


  

  EBCCProjector(shared_ptr<EBEncyclopedia<2, Real> >   & a_brit,
                shared_ptr<GeometryService<2> >        & a_geoserv,
                const DisjointBoxLayout                & a_vecgrids,
                const Box                              & a_domain,
                const Real                             & a_dx,
                const IntVect                          & a_nghost,
                const EBIBC                            & a_ebibc)
  {
    define(a_brit,          
           a_geoserv,       
           a_vecgrids,      
           a_domain,        
           a_dx,            
           a_nghost,        
           a_ebibc);
  }
  
  void define(shared_ptr<EBEncyclopedia<2, Real> >   & a_brit,
              shared_ptr<GeometryService<2> >        & a_geoserv,
              const DisjointBoxLayout                & a_vecgrids,
              const Box                              & a_domain,
              const Real                             & a_dx,
              const IntVect                          & a_nghost,
              const EBIBC                            & a_ebibc);
  
  /// Divide  a velocity field into its two Hodge Decomposition bits.
  /**
     Input: velo holds input vector field.
     Output: velo holds divergence-free component
             gphi holds pure gradient component.
   */
  void 
  project(EBLevelBoxData<CELL, DIM>   & a_velo,
          EBLevelBoxData<CELL, DIM>   & a_gphi,
          Real a_tolerance, unsigned int a_maxiterations);


  ///
  /**
     used internally and for  diagnostics.
     Take kappa*divergence of cc velocity field
   */
  void 
  kappaDivU(EBLevelBoxData<CELL, 1  > & a_divu,
            EBLevelBoxData<CELL, DIM> & a_velo);


  ///
  inline EBLevelBoxData<CELL, 1>& getRHSHolder()
  {
    return m_macprojector->m_rhs;
  }
  

  shared_ptr<EBMACProjector> m_macprojector;
protected:
  string m_conGradEBFluxName[DIM];
  string m_ncdivLabel;
  string m_nobcsLabel;
  Copier m_copyToCopier;
  void registerStencils();


  void
  computeAverageFaceToCellGradient(EBBoxData<CELL, Real, DIM> & a_gph,
                                   EBBoxData<CELL, Real,   1> & a_phi,
                                   const EBGraph              & a_graph,
                                   const Bx                   & a_grid,
                                   const unsigned int         & a_ibox);

  void
  kappaConservativeGradient(EBBoxData<CELL, Real, DIM> & a_gph,
                            EBBoxData<CELL, Real,   1> & a_phi,
                            const EBGraph              & a_graph,
                            const Bx                   & a_grid,
                            const unsigned int         & a_ibox);
  
  ///turns kappa grad into grad by volume weighted averaging
  void normalizeGradient(EBLevelBoxData<CELL, DIM>& a_gphi);
private:
  ///insert usual diatribe about strong vs. weak construction here
  void operator=(const EBCCProjector& a_opin);
  EBCCProjector();

};
#include "Chombo_NamespaceFooter.H"

#endif
