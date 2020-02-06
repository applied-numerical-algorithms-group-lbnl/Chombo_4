#ifndef _BCGVelAdvect_H_
#define _BCGVelAdvect_H_
#include <cmath>
#include <memory>
#include "Proto.H"
#include "Chombo_EBDictionary.H"
#include "Chombo_EBChombo.H"
#include "Chombo_EBEncyclopedia.H"
#include "EBAdvection.H"
#include "EBMACProjector.H"
#include "Chombo_EBLevelBoxData.H"
#include "Chombo_NamespaceHeader.H"

///class to advect scalars in an eb context (via Trebotich et al.)
class BCGVelAdvect: public EBAdvection
{
public:

  /// 
  BCGVelAdvect(shared_ptr<EBMACProjector>             & a_macproj,
               shared_ptr<EBEncyclopedia<2, Real> >   & a_brit,
               shared_ptr<GeometryService<2> >        & a_geoserv,
               shared_ptr<EBLevelBoxData<CELL, DIM> > & a_veloCell,
               const DisjointBoxLayout                & a_grids,
               const Box                              & a_domain,
               const Real                             & a_dx,
               const IntVect                          & a_nghost)    
  :EBAdvection(a_brit,
               a_geoserv,   
               a_veloCell,  
               a_grids,     
               a_domain,    
               a_dx,        
               a_nghost,    
               a_nghost)
  {
    m_macproj = a_macproj;
    m_nghost  = a_nghost;
    m_macScal.define(m_grids, m_nghost, m_graphs);
    m_macGrad.define(m_grids, m_nghost, m_graphs);
    m_macVect.define(m_grids, m_nghost, m_graphs);
  }


  ///
  virtual void 
  hybridDivergence(EBLevelBoxData<CELL, DIM>& a_divuu,
                   EBLevelBoxData<CELL, DIM>& a_inputVel,
                   const Real               & a_dt);

protected:

  shared_ptr<EBMACProjector>    m_macproj;
  EBLevelFluxData<1>            m_macScal;
  EBLevelFluxData<1>            m_macGrad;
  EBLevelFluxData<DIM>          m_macVect;
  IntVect m_nghost;

  //does hybrid divergence by component--includes resdistribution
  virtual void  assembleDivergence();


  void getAdvectionVelocity(EBLevelBoxData<CELL, DIM>& a_inputVel,
                            const Real               & a_dt);

  void getMACVectorVelocity(EBLevelBoxData<CELL, DIM>& a_inputVel,
                            const Real               & a_dt);

  void operator=(const BCGVelAdvect& a_opin);
  BCGVelAdvect();

};
#include "Chombo_NamespaceFooter.H"

#endif
