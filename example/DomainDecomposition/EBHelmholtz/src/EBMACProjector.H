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
                 const IntVect                          & a_nghost)
  {
    m_dx     = a_dx;
    m_grids  = a_grids;
    m_domain = a_domain;
    m_nghostSrc  = a_nghostsrc;
    m_nghostDst  = a_nghostdst;
    m_exchangeCopier.exchangeDefine(m_grids, m_nghostSrc);
    m_grids      = a_grids;      
    m_domain     = a_domain;     
    m_brit       = a_brit;

    defineInternals();
  }


  /// advance one time step (via Trebotich et al.) in  an eb context
  void 
  project(EBLevelFluxData<1>   & a_velo,
          EBLevelFluxData<1>   & a_gphi)
  {
  }

  void 
  divergence(EBLevelBoxData<CELL, 1> & a_divu
             EBLevelFluxData<1>      & a_velo)
  {
  }

protected:
  
  void defineInternals()
  {
  }

  //rule britannica
  shared_ptr<EBEncyclopedia<2, Real> >   m_brit;
  DisjointBoxLayout                      m_grids;
  Box                                    m_domain;
  IntVect                                m_nghost;
  Real m_dx;

  //internals
  shared_ptr<LevelData<EBGraph>  >       m_graphs;
  Copier                                 m_exchangeCopier;

private:


  void operator=(const EBMACProjector& a_opin);
  EBMACProjector();

};
#include "Chombo_NamespaceFooter.H"

#endif
