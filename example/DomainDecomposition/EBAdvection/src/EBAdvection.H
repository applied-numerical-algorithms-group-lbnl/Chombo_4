#ifndef _EBAdvection_H_
#define _EBAdvection_H_
#include <cmath>
#include <memory>
#include "Proto.H"
#include "EBDictionary.H"
#include "EBChombo.H"
#include "EBEncyclopedia.H"
#include "EBLevelBoxData.H"
#include "NamespaceHeader.H"

///class to advect scalars in an eb context (via Trebotich et al.)
class EBAdvection
{
public:


  typedef EBStencil<2, Real, CELL , CELL >  celltocell_sten;
  typedef EBStencil<2, Real, XFACE, XFACE>  xfactoxfac_sten;
  typedef EBStencil<2, Real, YFACE, YFACE>  yfactoyfac_sten;
  typedef EBStencil<2, Real, ZFACE, ZFACE>  zfactozfac_sten;

  typedef EBStencil<2, Real, CELL , XFACE>  celltoxfac_sten;
  typedef EBStencil<2, Real, CELL , YFACE>  celltoyfac_sten;
  typedef EBStencil<2, Real, CELL , ZFACE>  celltozfac_sten;

  typedef EBStencil<2, Real, XFACE, CELL >  xfactocell_sten;
  typedef EBStencil<2, Real, YFACE, CELL >  yfactocell_sten;
  typedef EBStencil<2, Real, ZFACE, CELL >  zfactocell_sten;



  /// 
  EBAdvection(shared_ptr<EBEncyclopedia<2, Real> >   & a_brit,
              shared_ptr<GeometryService<2> >        & a_geoserv,
              shared_ptr<EBLevelBoxData<CELL, DIM> > & a_veloCell,
              const DisjointBoxLayout                & a_grids,
              const Box                              & a_domain,
              const Real                             & a_dx,
              const IntVect                          & a_nghostsrc, 
              const IntVect                          & a_nghostdst);

  /// advance one time step (via Trebotich et al.) in  an eb context
  void 
  advance(EBLevelBoxData<CELL, 1>   & a_scal,
          const  Real               & a_dt);

protected:

  //rule britannica
  shared_ptr<EBEncyclopedia<2, Real> >   m_brit;
  DisjointBoxLayout                      m_grids;
  Box                                    m_domain;
  IntVect                                m_nghostSrc;
  IntVect                                m_nghostDst;
  shared_ptr<EBLevelBoxData<CELL, DIM> > m_veloCell;

  EBLevelBoxData<CELL, 1>                m_kappa;
  EBLevelBoxData<CELL, 1>                m_deltaM;
  EBLevelBoxData<CELL, 1>                m_nonConsDiv;
  EBLevelBoxData<CELL, 1>                m_hybridDiv;
  
  Copier                                 m_exchangeCopier;

private:
  void  defineData(shared_ptr<GeometryService<2> >        & a_geoserv);
  void  kappaConsDiv(EBLevelBoxData<CELL, 1>   & a_scal);
  void  nonConsDiv(  EBLevelBoxData<CELL, 1>   & a_scal);
  void redistribute();
  void fillKappa(const shared_ptr<GeometryService<2> >   & a_geoserv);
  void registerStencils();

  void operator=(const EBAdvection& a_opin);
  EBAdvection();

};
#include "NamespaceFooter.H"

#endif
