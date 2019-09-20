#ifndef _EBMultigrid_H_
#define _EBMultigrid_H_
#include <cmath>
#include <memory>
#include "Proto.H"
#include "EBDictionary.H"
#include "EBChombo.H"
#include "EBLevelBoxData.H"
#include "NamespaceHeader.H"
#if DIM==2
#define MG_NUM_COLORS 4
#else
#define MG_NUM_COLORS 8
#endif

using std::array;
/// Multigrid: solve a FV discretization of Poisson's equation on a rectangle.  periodic boundary conditions.
class EBMultigridLevel
{
public:
  typedef EBStencil<2, Real, CELL, CELL> ebstencil_t;
  typedef shared_ptr<EBDictionary<2, Real, CELL, CELL> > dictionary_t;
  
  ///
  EBMultigridLevel(dictionary_t                      & a_dictionary,
                   const Real                        & a_alpha,
                   const Real                        & a_beta,
                   const Real                        & a_dx,
                   const DisjointBoxLayout           & a_grids,
                   const string                      & a_stenname,
                   const string                      & a_dombcname,
                   const string                      & a_ebbcname,
                   const Box                         & a_domain,
                   const IntVect                     & a_nghostsrc, 
                   const IntVect                     & a_nghostdst);

  ///make coarse from fine 
  EBMultigridLevel(const EBMultigridLevel& a_finerLevel);

  
  ///  
  ~EBMultigridLevel()
  {
  }

  /// 
  void residual(EBLevelBoxData<CELL, 1>       & a_res,
                const EBLevelBoxData<CELL, 1> & a_phi,
                const EBLevelBoxData<CELL, 1> & a_rhs);


  /// 
  void relax(EBLevelBoxData<CELL, 1>      & a_phi,
             const EBLevelBoxData<CELL, 1>& a_rhs);


  /// average down residual to next coarser level.
  void restrictResidual(EBLevelBoxData<CELL, 1>      & a_resc,
                        const EBLevelBoxData<CELL, 1>& a_res);

  /// Piecewise constant interpolation of coarse correction to increment fine solution.
  void prolongIncrement(EBLevelBoxData<CELL, 1>      & a_phiFine,
                        const EBLevelBoxData<CELL, 1>& a_deltaCoarse);


  /// 
  void applyOp(EBLevelBoxData<CELL, 1>       & a_lph,
               const EBLevelBoxData<CELL, 1> & a_phi);

  /// Multigrid v-cycle.
  void vCycle(EBLevelBoxData<CELL, 1>      & a_phi,
              const EBLevelBoxData<CELL, 1>& a_rhs);

  std::shared_ptr<EBMultigridLevel> m_coarser;
private:
  void defineStencils();
  void defineCoarserObjects();


  /// weak construction introduces unnecessary complications
  EBMultigridLevel()
  {
  }

  bool   m_hasCoarser; 

  dictionary_t                      m_dictionary;
  Real                              m_alpha;   
  Real                              m_beta;    
  Real                              m_dx;      
  DisjointBoxLayout                 m_grids;   
  string                            m_stenname;
  string                            m_dombcname;
  string                            m_ebbcname;
  Box                               m_domain;
  EBLevelBoxData<CELL, 1>           m_resid;
  EBLevelBoxData<CELL, 1>           m_relCoef;
  EBLevelBoxData<CELL, 1>           m_lphi;
  EBLevelBoxData<CELL, 1>           m_deltaC;
  EBLevelBoxData<CELL, 1>           m_residC;
  IntVect m_nghostSrc, m_nghostDst;
  Copier                            m_exchangeCopier;
};


///class that outsiders actually call
class EBMultigrid
{
public:
  typedef shared_ptr<EBDictionary<2, Real, CELL, CELL> > dictionary_t;
  static int  s_numSmoothDown ;
  static int  s_numSmoothUp   ;

  /// 
  EBMultigrid(dictionary_t                      & a_dictionary,
              const Real                        & a_alpha,
              const Real                        & a_beta,
              const Real                        & a_dx,
              const DisjointBoxLayout           & a_grids,
              const string                      & a_stenname,
              const string                      & a_dombcname,
              const string                      & a_ebbcname,
              const Box                         & a_domain,
              const IntVect                     & a_nghostsrc, 
              const IntVect                     & a_nghostdst)
  {
    m_finest = std::shared_ptr<EBMultigridLevel>(
      new EBMultigridLevel(a_dictionary, 
                           a_alpha,      
                           a_beta,       
                           a_dx,         
                           a_grids,      
                           a_stenname,   
                           a_dombcname,  
                           a_ebbcname,   
                           a_domain, 
                           a_nghostsrc,
                           a_nghostdst));     
  }
      
  /// 
  void residual(EBLevelBoxData<CELL, 1>       & a_res,
                const EBLevelBoxData<CELL, 1> & a_phi,
                const EBLevelBoxData<CELL, 1> & a_rhs);

  Real norm(const EBLevelBoxData<CELL, 1>       & a_res) const;
  
  /// 
  void applyOp(EBLevelBoxData<CELL, 1>       & a_lph,
               const EBLevelBoxData<CELL, 1> & a_phi);

  /// Multigrid v-cycle.
  void vCycle(EBLevelBoxData<CELL, 1>      & a_phi,
              const EBLevelBoxData<CELL, 1>& a_rhs);

  /// average down residual to next coarser level.
  void restrictResidual(EBLevelBoxData<CELL, 1>      & a_resc,
                        const EBLevelBoxData<CELL, 1>& a_res)
  {
    m_finest->m_coarser->restrictResidual(a_resc, a_res);
  }

  /// Piecewise constant interpolation of coarse correction to increment fine solution.
  void prolongIncrement(EBLevelBoxData<CELL, 1>      & a_phiFine,
                        const EBLevelBoxData<CELL, 1>& a_deltaCoarse)
  {
    m_finest->m_coarser->prolongIncrement(a_phiFine, a_deltaCoarse);
  }

  /// 
  void relax(EBLevelBoxData<CELL, 1>      & a_phi,
             const EBLevelBoxData<CELL, 1>& a_rhs)
  {
    m_finest->relax(a_phi, a_rhs);
  }
private:
  /// weak construction introduces unnecessary complications
  EBMultigrid()
  {
  }

  std::shared_ptr<EBMultigridLevel> m_finest;
};
#include "NamespaceFooter.H"

#endif
