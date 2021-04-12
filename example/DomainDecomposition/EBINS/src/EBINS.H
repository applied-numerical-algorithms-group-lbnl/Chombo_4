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
        ParabolicSolverType                      a_solver,
        EBIBC                                    a_ibc,
        unsigned int                             a_num_species = 0,
        vector<Real> a_diffusionCoeffs = vector<Real>());

  void readDataFromCheckpoint(Real         & a_curTime,
                              unsigned int & a_curStep,
                              const string & a_checkpointName);
  
  /// 
  /**
     Set fixedDt to positive number to force that to be dt.
     Negative intervals turn intervalled stuff off (eg. if checkPointInterval < 0, no checkpointing)
     Negative fixed time step means to use a variable time step.
     This is our weird interface that says sending something non-sensical as an argument
     turns off that bit of functionality.    This quixotic  interface has been standard
     for at least thirty years.  At the very least, it is not over-specified 
     and that can save a lot of code. --dtg 3-12-2021

     Pressure initialization options (only applies if a_startingStep == 0):
     1- Stokes flow initialization -- set gphi = nu*lapl*velo component wise 
        (overrides all if a_stokesFlowInitialization is true)
     2- Standard Treb mode -- initialize pressure with what comes out of 
        the initial projection (v = P(v*), gp = (I-P)(v*).
        For this set stokesFlowInitialization = false, numPressureIterations=0
     3- Fixed point iteration for pressure a-la BCG.  
        For this set stokesFlowInitialization = false, numPressureIterations > 0.
  */
  void 
  run(unsigned int a_max_step,
      Real         a_max_time,
      unsigned int a_startingStep,
      Real         a_startingTime,
      Real         a_cfl,
      Real         a_fixedDt,
      Real         a_solverTolerance, 
      unsigned int a_pressureInitIters,   
      unsigned int a_maxSolverIterations ,
      int          a_plotfileInterval,      // how often to output data
      int          a_checkpointInterval,    // how often to save the state.
      bool         a_stokesFlowInitialization,
      Real         a_coveredval);  //for plot files



  //rule britannica
  shared_ptr<EBEncyclopedia<2, Real> >   m_brit;
  shared_ptr<GeometryService<2> >        m_geoserv;
  shared_ptr<EBLevelBoxData<CELL, DIM> > m_velo;
  shared_ptr<EBLevelBoxData<CELL, DIM> > m_gphi;
  shared_ptr<EBLevelBoxData<CELL, 1  > > m_scal;
  shared_ptr<EBLevelBoxData<CELL, DIM> > m_sour;
  shared_ptr<EBLevelBoxData<CELL, 1  > > m_rhsDiff;
  shared_ptr<EBLevelBoxData<CELL, 1  > > m_sourDiff;
  shared_ptr<EBLevelBoxData<CELL, DIM> > m_divuu;
  
  vector<Real> m_diffusionCoefs;
  vector<shared_ptr<EBLevelBoxData<CELL, 1> > > m_species;

  shared_ptr<EBAdvection>                m_advectOp;
  shared_ptr<EBCCProjector >             m_ccProj;
  shared_ptr<EBMACProjector>             m_macProj;
  shared_ptr<EBMultigrid>                m_helmholtzVelo;
  shared_ptr<EBMultigrid>                m_helmholtzSpec;
  shared_ptr<BCGVelAdvect>               m_bcgAdvect;
  shared_ptr<BaseEBParabolic>            m_heatSolverVelo;
  shared_ptr<BaseEBParabolic>            m_heatSolverSpec;
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
  Real         m_dt;
  
  bool m_eulerCalc;

  void  defineInternals(EBIBC                a_ibc,
                        unsigned int         a_num_species,
                        Real                 a_viscosity,
                        ParabolicSolverType  a_solver);
  
  void projectVelocityAndCorrectPressure(Real a_dt,
                                         Real a_tol,                         
                                         unsigned int a_maxIter);
  void getAdvectiveDerivative(Real a_dt, Real a_tol, unsigned int a_maxIter);
  //gets  ustar
  void advanceVelocityEuler(Real a_dt);
  void getReactionSourceTerm(unsigned int a_ispec);
  
  //gets  ustar
  void advanceVelocityNavierStokes(Real a_dt,                          
                                   Real a_tol,                         
                                   unsigned int a_maxSolverIterations);

  Real  computeDt(Real a_cfl) const;
  void initializePressure(Real         a_dt,
                          Real         a_tol,
                          unsigned int a_maxSolverIterations,
                          bool         a_stokesFlowInitialization,
                          unsigned int a_numPressureIterations);
  
  void advanceVelocityAndPressure(Real         a_dt,                          
                                  Real         a_tol,                         
                                  unsigned int a_maxSolverIterations);

  void advanceScalar(Real a_dt);
  
  void advanceSpecies(Real a_dt,
                      Real         a_tol,    
                      unsigned int a_maxIter);

  void writeCheckpointFile();
  
  void writePlotFile(Real a_coveredval) const;
  //forbidden for all the usual reasons
  EBINS();
  EBINS(const EBINS& a_input);
  EBINS& operator=(const EBINS& a_input);

};
#include "Chombo_NamespaceFooter.H"

#endif
