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


  
  /// 
  EBCCProjector(shared_ptr<EBEncyclopedia<2, Real> >   & a_brit,
                shared_ptr<GeometryService<2> >        & a_geoserv,
                const DisjointBoxLayout                & a_vecgrids,
                const Box                              & a_domain,
                const Real                             & a_dx,
                const IntVect                          & a_nghost,
                string a_bcname = StencilNames::Neumann)
  {
    string bcnames[2*DIM];
    for(int ivec = 0; ivec < 2*DIM; ivec++)
    {
      bcnames[ivec] = a_bcname;
    }
    define(a_brit,          
           a_geoserv,       
           a_vecgrids,      
           a_domain,        
           a_dx,            
           a_nghost,        
           bcnames);
  }

  EBCCProjector(shared_ptr<EBEncyclopedia<2, Real> >   & a_brit,
                shared_ptr<GeometryService<2> >        & a_geoserv,
                const DisjointBoxLayout                & a_vecgrids,
                const Box                              & a_domain,
                const Real                             & a_dx,
                const IntVect                          & a_nghost,
                string                                   a_bcnames[2*DIM])
  {
    define(a_brit,          
           a_geoserv,       
           a_vecgrids,      
           a_domain,        
           a_dx,            
           a_nghost,        
           a_bcnames);
  }
  
  void define(shared_ptr<EBEncyclopedia<2, Real> >   & a_brit,
              shared_ptr<GeometryService<2> >        & a_geoserv,
              const DisjointBoxLayout                & a_vecgrids,
              const Box                              & a_domain,
              const Real                             & a_dx,
              const IntVect                          & a_nghost,
              string                                   a_bcnames[2*DIM]);
  
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

  void registerStencils();


private:

  ///insert usual diatribe about strong vs. weak construction here
  void operator=(const EBCCProjector& a_opin);
  EBCCProjector();

};
#include "Chombo_NamespaceFooter.H"

#endif
