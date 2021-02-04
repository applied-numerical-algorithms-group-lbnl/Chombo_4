#ifndef _EBINS_finally__H_
#define _EBINS_finally__H_
#include <cmath>
#include <memory>
#include "Proto.H"
#include "EBAdvection.H"
#include "EBParabolicIntegrators.H"
#include "EBMACProjector.H"
#include "EBCCProjector.H"
#include "BCGVelAdvect.H"
#include "EBIBC.H"
#include "Chombo_NamespaceHeader.H"

///class to run the incompressible algorithm described in Trebotich, et al.
class EBINS
{
public:

  enum ParabolicSolverType {BackwardEuler=0, CrankNicolson, TGA};
  
  ///
  /**Strictly speaking num_species and diffusionCoefficients are 
     not independent arguments since diffuisonCoefficients has
     to be at least a_num_species long.   I could use the vector's size 
     to be num_species but I find that API abstruse.   So the interface is clear,
     albeit slightly overspecified.   I make no apologies.
  */
  EBINS(shared_ptr<EBEncyclopedia<2, Real> >   & a_brit,
        shared_ptr<GeometryService<2> >        & a_geoserv,
        const DisjointBoxLayout                & a_grids,
        const Box                              & a_domain,
        const Real                             & a_dx,
        const Real                             & a_viscosity,
        const IntVect                          & a_nghost,
        ParabolicSolverType a_solver,
        EBIBC               a_ibc,
        int                 a_num_species = 0,  
        vector<Real> a_diffusionCoeffs = vector<Real>()); 

  /// set fixedDt to positive number to force that to be dt
  void 
  run(unsigned int a_max_step,
      Real         a_max_time,
      Real         a_cfl,
      Real         a_fixedDt             = -1.0,
      Real         a_solverTolerance     = 1.0e-10,
      unsigned int a_pressureInitIters   = 1,
      unsigned int a_maxSolverIterations = 27,
      int          a_outputInterval      = 1,
      Real         a_coveredval          = 0.0);



  //rule britannica
  shared_ptr<EBEncyclopedia<2, Real> >   m_brit;
  shared_ptr<GeometryService<2> >        m_geoserv;
  shared_ptr<EBLevelBoxData<CELL, DIM> > m_velo;
  shared_ptr<EBLevelBoxData<CELL, DIM> > m_gphi;
  shared_ptr<EBLevelBoxData<CELL, 1  > > m_scal;
  shared_ptr<EBLevelBoxData<CELL, DIM> > m_sour;
  shared_ptr<EBLevelBoxData<CELL, DIM> > m_divuu;
  
  vector<Real> m_diffusionCoefs;
  vector<shared_ptr<EBLevelBoxData<CELL, 1> > > m_species;

  shared_ptr<EBAdvection>                m_advectOp;
  shared_ptr<EBCCProjector >             m_ccProj;
  shared_ptr<EBMACProjector>             m_macProj;
  shared_ptr<EBMultigrid>                m_helmholtz;
  shared_ptr<BCGVelAdvect>               m_bcgAdvect;
  shared_ptr<BaseEBParabolic>            m_heatSolver;
  shared_ptr<LevelData<EBGraph> >        m_graphs;
  
  DisjointBoxLayout                      m_grids;
  Box                                    m_domain;
  IntVect                                m_nghost;
  Real                                   m_dx;
  Real                                   m_viscosity;
  Copier                                 m_exchangeCopier;
  Copier                                 m_copyCopier;


private:
  
  unsigned int m_step;
  Real         m_time;
  
  bool m_eulerCalc;
  void projectVelocityAndCorrectPressure(Real a_dt,
                                         Real a_tol,                         
                                         unsigned int a_maxIter);
  void getAdvectiveDerivative(Real a_dt, Real a_tol, unsigned int a_maxIter);
  //gets  ustar
  void advanceVelocityEuler(Real a_dt);
  //gets  ustar
  void advanceVelocityNavierStokes(Real a_dt,                          
                                   Real a_tol,                         
                                   unsigned int a_maxSolverIterations);

  Real  computeDt(Real a_cfl) const;
  void initializePressure(Real a_dt,
                          Real a_tol,
                          unsigned int a_numPressureIterations,
                          unsigned int a_maxSolverIterations);
  void advanceVelocityAndPressure(Real a_dt,                          
                                  Real a_tol,                         
                                  unsigned int a_maxSolverIterations);

  void advanceScalar(Real a_dt);

  void outputToFile(unsigned int a_step, Real a_coveredval, Real a_dt, Real a_time) const;
  //forbidden for all the usual reasons
  EBINS();
  EBINS(const EBINS& a_input);
  EBINS& operator=(const EBINS& a_input);

};
#include "Chombo_NamespaceFooter.H"

#endif
