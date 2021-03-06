#ifndef _EBDarcy_finally__H_
#define _EBDarcy_finally__H_
#include <cmath>
#include <memory>
#include "Proto.H"
#include "EBAdvection.H"
#include "EBParabolicIntegrators.H"
#include "EBMACProjector.H"
#include "EBCCProjector.H"
#include "EBIBC.H"
#include "Chombo_NamespaceHeader.H"

///class to run the incompressible algorithm described in Trebotich, et al.
class EBDarcy
{
public:

  enum ParabolicSolverType {BackwardEuler=0, CrankNicolson, TGA};
  
  ///
  /**
     This application may not make any sense unless you send in exactly the same EBIBC
     that is provided in ../exec/darcy.cpp.   So it really should not be 
     part of the public interface.   It does allow for a lot of code reuse, however,
     so I remain unapologetic.  --dtg
  */
  EBDarcy(shared_ptr<EBEncyclopedia<2, Real> >   & a_brit,
          shared_ptr<GeometryService<2> >        & a_geoserv,
          const DisjointBoxLayout                & a_grids,
          const Box                              & a_domain,
          const Real                             & a_dx,
          const Real                             & a_viscosity,
          const Real                             & a_permeability,
          const Real                             & a_diffusivity,
          const IntVect                          & a_nghost,
          ParabolicSolverType a_solver,
          EBIBC               a_ibc);

  /// set fixedDt to positive number to force that to be dt
  void 
  run(unsigned int a_max_step,
      Real         a_max_time,
      Real         a_cfl,
      Real         a_fixedDt             = -1.0,
      Real         a_solverTolerance     = 1.0e-10,
      unsigned int a_maxSolverIterations = 27,
      int          a_outputInterval      = 1,
      Real         a_coveredval          = 0.0);



  //rule britannica
  shared_ptr<EBEncyclopedia<2, Real> >   m_brit;
  shared_ptr<GeometryService<2> >        m_geoserv;
  shared_ptr<EBLevelBoxData<CELL, DIM> > m_velo;
  shared_ptr<EBLevelBoxData<CELL, DIM> > m_gphi;
  
  shared_ptr<EBLevelBoxData<CELL, 1  > > m_scal;
  shared_ptr<EBLevelBoxData<CELL, 1  > > m_rhs;
  shared_ptr<EBLevelBoxData<CELL, 1  > > m_sour;

  shared_ptr<EBAdvection>                m_advectOp;
  shared_ptr<EBCCProjector >             m_ccProj;
  shared_ptr<EBMultigrid>                m_helmholtz;
  shared_ptr<BaseEBParabolic>            m_heatSolver;
  shared_ptr<LevelData<EBGraph> >        m_graphs;
  
  DisjointBoxLayout                      m_grids;
  Box                                    m_domain;
  IntVect                                m_nghost;
  Real                                   m_dx;

  Real                                   m_viscosity   ; 
  Real                                   m_permeability;
  Real                                   m_diffusivity ;  
  Copier                                 m_exchangeCopier;
  Copier                                 m_copyCopier;
  EBIBC                                  m_ibc;

private:
  
  unsigned int m_step;
  Real         m_time;
  
  Real  computeDt(Real a_cfl) const;

  void advanceScalar(Real a_dt,
                     Real         a_tol,    
                     unsigned int a_maxIter);

  void getReactionSourceTerm();
  void
  initializeVelocity(Real         a_tol,
                     unsigned int a_maxIter);

  void outputToFile(unsigned int a_step, Real a_coveredval, Real a_dt, Real a_time) const;
  //forbidden for all the usual reasons
  EBDarcy();
  EBDarcy(const EBDarcy& a_input);
  EBDarcy& operator=(const EBDarcy& a_input);

};
#include "Chombo_NamespaceFooter.H"

#endif
