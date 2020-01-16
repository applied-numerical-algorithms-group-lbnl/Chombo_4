#ifndef _EBParabolic_H_
#define _EBParabolic_H_
#include <cmath>
#include <memory>
#include "Proto.H"
#include "Chombo_EBDictionary.H"
#include "Chombo_EBChombo.H"
#include "Chombo_EBLevelBoxData.H"
#include "EBMultigrid.H"
#include "Chombo_NamespaceHeader.H"
///
/**
 */
class BaseEBParabolic
{
public:
  ///the boundary conditions for the EBMultigrid have to be correct.  Alpha and Beta do not as they are constantly reset
  BaseEBParabolic(shared_ptr<EBMultigrid>            a_diffusionSolver,
                  shared_ptr<GeometryService<2> >    a_geoserv,
                  const DisjointBoxLayout          & a_grids,
                  const Box                        & a_domain,
                  const IntVect                    & a_nghost)
  {
    m_diffusionSolver = a_diffusionSolver;
    m_geoserv         = a_geoserv;
    m_grids           = a_grids;
    m_domain          = a_domain;
    m_nghost          = a_nghost;
    m_graphs          = a_geoserv->getGraphs(m_domain);

  }

  ///
  virtual ~BaseEBParabolic()
  {
  }

  ///this one will be different for different integrators
  /**
     source term assumed to be already multiplied by kappa (kappa div F for example)
  */
  virtual void advanceOneStep( EBLevelBoxData<CELL, 1>       &  a_phinew,
                               const EBLevelBoxData<CELL, 1> &  a_source,
                               const Real                    &  a_diffCoef,
                               const Real                    &  a_dt,
                               const Real                    &  a_tolerance,
                               const unsigned int            &  a_maxIterations) = 0;

  ///
  /**
     this is the conservative way to advance a solution with diffusion.
     source term assumed to be already multiplied by kappa (kappa div F for example)
  */
  inline
  void computeDiffusion( EBLevelBoxData<CELL, 1>       &  a_diffusionTerm,
                         const EBLevelBoxData<CELL, 1> &  a_phiold,
                         const Real                    &  a_diffCoef,
                         const Real                    &  a_dt,
                         const Real                    &  a_tolerance,
                         const unsigned int            &  a_maxIterations);
                        
protected:
  shared_ptr<EBMultigrid> m_diffusionSolver;  
  DisjointBoxLayout       m_grids;
  Box                     m_domain;
  IntVect                 m_nghost;
  shared_ptr<GeometryService<2> >    m_geoserv;
  shared_ptr<LevelData<EBGraph>  >   m_graphs; 

private:
  
  //outlawed for all the usual reasons
  BaseEBParabolic();
  BaseEBParabolic(const  BaseEBParabolic& a_input);
  void operator= (const  BaseEBParabolic& a_input);
};

///
class EBBackwardEuler: public BaseEBParabolic
{
public:
  ///
  EBBackwardEuler(shared_ptr<EBMultigrid>              a_diffusionSolver,
                  shared_ptr<GeometryService<2> >    a_geoserv,
                  const DisjointBoxLayout          & a_grids,
                  const Box                        & a_domain,
                  const IntVect                    & a_nghost)
    :BaseEBParabolic(a_diffusionSolver, 
                     a_geoserv,         
                     a_grids,           
                     a_domain,          
                     a_nghost)
  {
    m_rhs.define(m_grids, m_nghost, m_graphs);
  }         

  ///
  /** 
      source term assumed to be already multiplied by kappa (kappa div F for example)
  */
  void advanceOneStep( EBLevelBoxData<CELL, 1>       &  a_phi,
                       const EBLevelBoxData<CELL, 1> &  a_source,
                       const Real                    &  a_diffCoef,
                       const Real                    &  a_dt,
                       const Real                    &  a_tolerance,
                       const unsigned int            &  a_maxIterations);

protected:
  
  EBLevelBoxData<CELL, 1>  m_rhs;
private:
  //outlawed for all the usual reasons
  EBBackwardEuler();
  EBBackwardEuler(const  EBBackwardEuler& a_input);
  void operator= (const  EBBackwardEuler& a_input);

};


#include "Chombo_NamespaceFooter.H"
#endif


