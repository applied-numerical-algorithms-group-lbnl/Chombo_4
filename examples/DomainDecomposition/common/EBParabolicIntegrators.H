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
  typedef CH4_Data_Choreography::DistributedData<EBGraph> graph_distrib_t;
 
  ///the boundary conditions for the EBMultigrid have to be correct.  Alpha and Beta do not as they are constantly reset
  BaseEBParabolic(shared_ptr<EBMultigrid>            a_diffusionSolver,
                  shared_ptr<GeometryService<2> >    a_geoserv,
                  const DisjointBoxLayout          & a_grids,
                  const Box                        & a_domain,
                  const IntVect                    & a_nghost,
                  bool a_printStuff = false)
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
  shared_ptr<GeometryService<2> > m_geoserv;
  shared_ptr<graph_distrib_t  >   m_graphs; 

private:
  
//  //outlawed for all the usual reasons
//  BaseEBParabolic();
//  BaseEBParabolic(const  BaseEBParabolic& a_input);
//  void operator= (const  BaseEBParabolic& a_input);
};

///
class EBBackwardEuler: public BaseEBParabolic
{
public:
  ///
  EBBackwardEuler(shared_ptr<EBMultigrid>            a_diffusionSolver,
                  shared_ptr<GeometryService<2> >    a_geoserv,
                  const DisjointBoxLayout          & a_grids,
                  const Box                        & a_domain,
                  const IntVect                    & a_nghost,
                  bool a_printStuff = false)
    :BaseEBParabolic(a_diffusionSolver, 
                     a_geoserv,         
                     a_grids,           
                     a_domain,          
                     a_nghost,
                     a_printStuff)
  {
    if(a_printStuff)
    {
      pout() << "EBBackwardEuler constructor: building rhs" << endl;
    }
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
//private:
//  //outlawed for all the usual reasons
//  EBBackwardEuler();
//  EBBackwardEuler(const  EBBackwardEuler& a_input);
//  void operator= (const  EBBackwardEuler& a_input);
//
};


///
class EBCrankNicolson: public BaseEBParabolic
{
public:
  ///
  EBCrankNicolson(shared_ptr<EBMultigrid>            a_diffusionSolver,
                  shared_ptr<GeometryService<2> >    a_geoserv,
                  const DisjointBoxLayout          & a_grids,
                  const Box                        & a_domain,
                  const IntVect                    & a_nghost,
                  bool a_printStuff = false)
    :BaseEBParabolic(a_diffusionSolver, 
                     a_geoserv,         
                     a_grids,           
                     a_domain,          
                     a_nghost,
                     a_printStuff)
  {
    m_rhs.define(     m_grids, m_nghost, m_graphs);
    m_kappaLph.define(m_grids, m_nghost, m_graphs);
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
  EBLevelBoxData<CELL, 1>  m_kappaLph;
//private:
//  //outlawed for all the usual reasons
//  EBCrankNicolson();
//  EBCrankNicolson(const  EBCrankNicolson& a_input);
//  void operator= (const  EBCrankNicolson& a_input);

};


///
class EBTGA: public BaseEBParabolic
{
public:
  ///
  EBTGA(shared_ptr<EBMultigrid>            a_diffusionSolver,
        shared_ptr<GeometryService<2> >    a_geoserv,
        const DisjointBoxLayout          & a_grids,
        const Box                        & a_domain,
        const IntVect                    & a_nghost,
        bool a_printStuff = false)
    :BaseEBParabolic(a_diffusionSolver, 
                     a_geoserv,         
                     a_grids,           
                     a_domain,          
                     a_nghost,
                     a_printStuff)
  {
    m_rhst.define(     m_grids, m_nghost, m_graphs);
    m_srct.define(     m_grids, m_nghost, m_graphs);
    m_phit.define(     m_grids, m_nghost, m_graphs);

    Real tgaEpsilon = 1.e-12;
#ifdef CH_USE_FLOAT
    tgaEpsilon = sqrt(tgaEpsilon);
#endif
    Real a = 2.0 - sqrt(2.0) - tgaEpsilon;
    m_mu1 = (a - sqrt( a*a - 4.0*a + 2.0))/2.0 ;
    m_mu2 = (a + sqrt( a*a - 4.0*a + 2.0))/2.0 ;
    m_mu3 = (1.0 - a);
    m_mu4 = 0.5 - a;

    Real discr = sqrt(a*a - 4.0*a + 2.0);
    m_r1 = (2.0*a - 1.0)/(a + discr);
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
  Real m_mu1, m_mu2, m_mu3, m_mu4, m_r1;
  
  EBLevelBoxData<CELL, 1>  m_rhst;
  EBLevelBoxData<CELL, 1>  m_srct;
  EBLevelBoxData<CELL, 1>  m_phit;

//private:
//  //outlawed for all the usual reasons
//  EBTGA();
//  EBTGA(const            EBTGA& a_input);
//  void operator= (const  EBTGA& a_input);

};


#include "Chombo_NamespaceFooter.H"
#endif


