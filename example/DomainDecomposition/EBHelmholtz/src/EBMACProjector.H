#ifndef _EBMACProjector_H_
#define _EBMACProjector_H_
#include <cmath>
#include <memory>
#include "Proto.H"
#include "Chombo_EBDictionary.H"
#include "Chombo_EBChombo.H"
#include "Chombo_EBEncyclopedia.H"
#include "Chombo_EBLevelBoxData.H"
#include "Chombo_EBLevelFluxData.H"
#include "EBMultigrid.H"
#include "Chombo_NamespaceHeader.H"

///class to advect scalars in an eb context (via Trebotich et al.)
class EBMACProjector
{
public:


  
  /// 
  EBMACProjector(shared_ptr<EBEncyclopedia<2, Real> >   & a_brit,
                 shared_ptr<GeometryService<2> >        & a_geoserv,
                 const DisjointBoxLayout                & a_grids,
                 const Box                              & a_domain,
                 const Real                             & a_dx,
                 const IntVect                          & a_nghost);

  /// Divide  a velocity field into its two Hodge Decomposition bits.
  /**
     Input: velo holds input vector field.
     Output: velo holds divergence-free component
             gphi holds pure gradient component.
   */
  void 
  project(EBLevelFluxData<1>   & a_velo,
          EBLevelFluxData<1>   & a_gphi,
          Real a_tolerance, unsigned int a_maxiterations);


  ///
  /**
     used internally and for  diagnostics.
     Take divergence of mac velocity field
   */
  void 
  divergence(EBLevelBoxData<CELL, 1> & a_divu,
             EBLevelFluxData<1>      & a_velo);

  
  ///
  EBLevelBoxData<CELL, 1>& getRHSHolder();

protected:
  
  void defineInternals(shared_ptr<GeometryService<2> >        & a_geoserv);

  void registerStencils();
  //rule britannica
  shared_ptr<EBEncyclopedia<2, Real> >   m_brit;
  DisjointBoxLayout                      m_grids;
  Box                                    m_domain;
  IntVect                                m_nghost;
  Real m_dx;

  //internals
  shared_ptr<LevelData<EBGraph>  >       m_graphs;
  Copier                                 m_exchangeCopier;
  EBLevelBoxData<CELL, 1>                m_rhs;
  EBLevelBoxData<CELL, 1>                m_phi;

  shared_ptr<EBMultigrid>                m_solver;
private:


  void operator=(const EBMACProjector& a_opin);
  EBMACProjector();

};
#include "Chombo_NamespaceFooter.H"

#endif
