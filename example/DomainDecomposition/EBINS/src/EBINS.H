#ifndef _EBINS_finally__H_
#define _EBINS_finally__H_
#include <cmath>
#include <memory>
#include "Proto.H"
#include "EBAdvection.H"
#include "EBParabolicIntegrators.H"
#include "EBMACProjector.H"
#include "EBCCProjector.H"
#include "Chombo_NamespaceHeader.H"

///class to run the incompressible algorithm described in Trebotich, et al.
class EBINS
{
public:

  enum ParabolicSolverType(BackwardEuler=0, CrankNicolson, TGA)
  
  /// 
  EBINS(shared_ptr<EBEncyclopedia<2, Real> >   & a_brit,
        shared_ptr<GeometryService<2> >        & a_geoserv,
        const DisjointBoxLayout                & a_grids,
        const Box                              & a_domain,
        const Real                             & a_dx,
        const Real                             & a_solverTolerance,
        const unsigned int                     & a_maxSolverIterations,
        const Real                             & a_coveredval,
        const IntVect                          & a_nghostdst,
        ParabolicSolverType a_solver);

  /// 
  void 
  run(unsigned int a_max_step,
      Real         a_max_time,
      Real         a_cfl,
      unsigned int a_outputInterval);



  //rule britannica
  shared_ptr<EBEncyclopedia<2, Real> >   m_brit;
  shared_ptr<EBAdvection>                m_advectOp;
  shared_ptr<EBMACProjector>             m_macProj;
  shared_ptr<EBCCProjector >             m_ccProj;
  shared_ptr<BaseEBParabolic>            m_heatSolver;
  shared_ptr<EBLevelBoxData<CELL, DIM> > m_velo;
  shared_ptr<EBLevelBoxData<CELL, DIM> > m_gphi;
  shared_ptr<EBLevelBoxData<CELL, 1  > > m_scal;
#if DIM==2  
  shared_ptr<EBLevelBoxData<CELL, 1  > > m_vort;
#else  
  shared_ptr<EBLevelBoxData<CELL, DIM> > m_vort;
#endif  
  shared_ptr<LevelData<EBGraph> >        m_graphs;
  
  DisjointBoxLayout                      m_grids;
  Box                                    m_domain;
  IntVect                                m_nghost;
  Real                                   m_dx, m_dt, m_time;


private:

  void computeDt(Real a_cfl) const;
  void advanceVelocity();
  void advanceScalar();
  void computeVorticity();
  void outputToFile(unsigned int a_step) const;
  //forbidde for all the usual reasons
  EBINS();
  EBINS(const EBINS& a_input);
  EBINS& operator=(const EBINS& a_input);

};
#include "Chombo_NamespaceFooter.H"

#endif
