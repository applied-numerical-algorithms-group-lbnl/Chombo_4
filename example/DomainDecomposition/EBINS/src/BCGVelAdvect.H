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
               shared_ptr<EBMultigrid>                & a_helmholtz,
               shared_ptr<EBEncyclopedia<2, Real> >   & a_brit,
               shared_ptr<GeometryService<2> >        & a_geoserv,
               shared_ptr<EBLevelBoxData<CELL, DIM> > & a_veloCell,
               const DisjointBoxLayout                & a_grids,
               const Box                              & a_domain,
               const Real                             & a_dx,
               const Real                             & a_viscosity,
               const IntVect                          & a_nghost,
               Real a_tolerance, unsigned int a_maxIter)    
  :EBAdvection(a_brit,
               a_geoserv,   
               a_veloCell,  
               a_grids,     
               a_domain,    
               a_dx,        
               a_nghost,    
               a_nghost)
  {
    m_macproj    = a_macproj;
    m_helmholtz  = a_helmholtz;
    m_nghost     = a_nghost;
    m_viscosity  = a_viscosity;
    m_tol        = a_tolerance;
    m_maxIter    = a_maxIter;
    m_advectionVel.define(m_grids, m_nghost, m_graphs);
    m_macGradient.define(m_grids, m_nghost, m_graphs);
    m_macVelocity.define(m_grids, m_nghost, m_graphs);

    //for tangential velocity correction at faces
    m_brit->registerFaceToCell( StencilNames::TanVelCorrect, StencilNames::NoBC,  StencilNames::NoBC, m_domain, m_domain, false);
  }


  ///
  virtual void 
  hybridDivergence(EBLevelBoxData<CELL, DIM>& a_divuu,
                   EBLevelBoxData<CELL, DIM>& a_inputVel,
                   const Real               & a_dt);

protected:

  shared_ptr<EBMACProjector>    m_macproj;
  shared_ptr<EBMultigrid>       m_helmholtz;
  EBLevelFluxData<1>            m_advectionVel;
  Real m_tol;
  unsigned int m_maxIter;
  EBLevelFluxData<1>            m_macGradient;
  EBLevelFluxData<DIM>          m_macVelocity;
  IntVect m_nghost;
  Real    m_viscosity;
  //does hybrid divergence by component--includes resdistribution
  void  assembleDivergence(EBLevelBoxData<CELL, DIM>& a_divuu, const Real& a_dt);  
  void  copyComp(EBFluxData<Real, 1>&  a_dst,
                 EBFluxData<Real, 1>&  a_src,
                 int a_vecDir);

  void getUpwindState(EBFluxData<Real, 1>&  a_upwindScal,
                      EBFluxData<Real, 1>&  a_faceCentVelo,
                      EBFluxData<Real, 1>&  a_scalLo,
                      EBFluxData<Real, 1>&  a_scalHi);

  void getAdvectionVelocity(EBLevelBoxData<CELL, DIM>& a_inputVel,
                            const Real               & a_dt);

  void getMACVectorVelocity(EBLevelBoxData<CELL, DIM>& a_inputVel,
                            const Real               & a_dt);
  void correctVectorVelocity();

  void assembleFlux(EBFluxData<Real, 1>& a_fcflux,
                    EBFluxData<Real, 1>& a_scalar,
                    EBFluxData<Real, 1>& a_fcvel);

  void operator=(const BCGVelAdvect& a_opin);
  BCGVelAdvect();

};
#include "Chombo_NamespaceFooter.H"

#endif
