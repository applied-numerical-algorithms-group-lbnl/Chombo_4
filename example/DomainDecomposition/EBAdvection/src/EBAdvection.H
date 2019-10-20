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

  DisjointBoxLayout m_grids;
  Box               m_domain;
  EBLevelBoxData<CELL    , 1>  m_kappa;
  EBLevelBoxData<CELL    , 1>  m_kappaDiv;
  EBLevelBoxData<BOUNDARY, 1>  m_deltaM;
  EBLevelBoxData<BOUNDARY, 1>  m_nonConsDiv;
  EBLevelBoxData<BOUNDARY, 1>  m_hybridDiv;
  
private:
  void  kappaConsDiv(EBLevelBoxData<CELL, 1>   & a_scal);
  void  nonConsDiv(  EBLevelBoxData<CELL, 1>   & a_scal);
  void redistribute();
  void fillKappa(const shared_ptr<GeometryService<2> >   & a_geoserv);

  void operator=(const EBAdvection& a_opin);
  EBAdvection();

};
#include "NamespaceFooter.H"

#endif
