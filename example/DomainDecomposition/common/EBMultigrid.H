#ifndef _EBMultigrid_H_
#define _EBMultigrid_H_
#include <cmath>
#include <memory>
#include "Proto.H"
#include "Chombo_EBDictionary.H"
#include "Chombo_EBChombo.H"
#include "Chombo_EBLevelBoxData.H"
#include "Chombo_NamespaceHeader.H"
#if DIM==2
#define MG_NUM_COLORS 4
#else
#define MG_NUM_COLORS 8
#endif

using std::array;
class EBMultigrid;
class EBRelaxSolver;
/// Multigrid: solve a FV discretization of Poisson's equation on a rectangle.  periodic boundary conditions.
class EBMultigridLevel
{
public:
  friend class EBMultigrid;

  typedef EBStencil<2, Real, CELL, CELL> ebstencil_t;
  typedef shared_ptr<EBDictionary<2, Real, CELL, CELL> > dictionary_t;
  
  ///
  EBMultigridLevel(dictionary_t                      & a_dictionary,
                   shared_ptr<GeometryService<2> >   & a_geoserv,
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
  EBMultigridLevel(const EBMultigridLevel            & a_finerLevel,
                   shared_ptr<GeometryService<2> >   & a_geoserv);

  ///
  void
  preCond(EBLevelBoxData<CELL, 1>       & a_phi,
          const EBLevelBoxData<CELL, 1> & a_rhs);

  
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
             const EBLevelBoxData<CELL, 1>& a_rhs,
             int a_niter);


  /// average down residual to next coarser level.
  void restrictResidual(EBLevelBoxData<CELL, 1>      & a_resc,
                        const EBLevelBoxData<CELL, 1>& a_res);

  /// Piecewise constant interpolation of coarse correction to increment fine solution.
  void prolongIncrement(EBLevelBoxData<CELL, 1>      & a_phiFine,
                        const EBLevelBoxData<CELL, 1>& a_deltaCoarse);


  /// 
  void applyOp(EBLevelBoxData<CELL, 1>       & a_lph,
               const EBLevelBoxData<CELL, 1> & a_phi);

  ///for tga
  void applyOpNeumann(EBLevelBoxData<CELL, 1>       & a_lph,
                      const EBLevelBoxData<CELL, 1> & a_phi);

  /// Multigrid v-cycle.
  void vCycle(EBLevelBoxData<CELL, 1>      & a_phi,
              const EBLevelBoxData<CELL, 1>& a_rhs);

  std::shared_ptr<EBMultigridLevel> m_coarser;

  EBLevelBoxData<CELL, 1> & getKappa()
  {
    return m_kappa;
  }
  
  inline void resetAlphaAndBeta(const Real& a_alpha,
                                const Real& a_beta)
  {
    m_alpha = a_alpha;
    m_beta  = a_beta;
    if(m_hasCoarser)
    {
      m_coarser->resetAlphaAndBeta(a_alpha, a_beta);
    }
  }
private:

  void defineCoarserObjects(shared_ptr<GeometryService<2> >   & a_geoserv);


//need the volume fraction in a data holder so we can evaluate kappa*alpha I 
  void  fillKappa(const shared_ptr<GeometryService<2> >   & a_geoserv,
                  const shared_ptr<LevelData<EBGraph> >   & a_graphs);

  /// weak construction introduces unnecessary complications
  EBMultigridLevel()
  {
  }

  bool   m_hasCoarser; 
  unsigned int m_depth;
  Real                              m_alpha;   
  Real                              m_beta;    
  Real                              m_dx;      
  DisjointBoxLayout                 m_grids;   
  string                            m_stenname;
  string                            m_neumname;
  string                            m_dombcname;
  string                            m_ebbcname;
  dictionary_t                      m_dictionary;
#if DIM==2
  static const unsigned int s_ncolors = 4;  
#else
  static const unsigned int s_ncolors = 8;
#endif
///prolongation has ncolors stencils
  string                            m_prolongationName[s_ncolors];
  string                            m_restrictionName;

  string                            m_nobcname;
  Box                               m_domain;
  EBLevelBoxData<CELL, 1>           m_resid;
  EBLevelBoxData<CELL, 1>           m_residC;
  EBLevelBoxData<CELL, 1>           m_deltaC;
  EBLevelBoxData<CELL, 1>           m_kappa;

  IntVect m_nghostSrc, m_nghostDst;
  Copier                            m_exchangeCopier;
  shared_ptr<EBRelaxSolver>         m_bottomSolver;
};


///class that outsiders actually call
class EBMultigrid
{
public:
  typedef shared_ptr<EBDictionary<2, Real, CELL, CELL> > dictionary_t;
  static int  s_numSmoothDown ;
  static int  s_numSmoothUp   ;
  static bool s_useWCycle     ;

  /// 
  EBMultigrid(dictionary_t                      & a_dictionary,
              shared_ptr<GeometryService<2> >   & a_geoserv,
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
                           a_geoserv,
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

    const shared_ptr<LevelData<EBGraph>  > graphs = a_geoserv->getGraphs(a_domain);
    m_res.define(a_grids, a_nghostsrc  , graphs);
    m_cor.define(a_grids, a_nghostsrc  , graphs);
  }

  ///
  void 
  solve(EBLevelBoxData<CELL, 1>       & a_phi,
        const EBLevelBoxData<CELL, 1> & a_rhs,
        const Real                    & a_tolerance,
        const unsigned int            & a_maxIterations,
        bool a_initToZero = true);

  EBLevelBoxData<CELL, 1> & getKappa()
  {
    return m_finest->getKappa();
  }

  inline void resetAlphaAndBeta(const Real& a_alpha,
                                const Real& a_beta)
  {
    m_finest->resetAlphaAndBeta(a_alpha, a_beta);
  }

  /// 
  void applyOp(EBLevelBoxData<CELL, 1>       & a_lph,
               const EBLevelBoxData<CELL, 1> & a_phi);


  /// for TGA force bcs to be homogeneous neumann
  inline
  void applyOpNeumann(EBLevelBoxData<CELL, 1>       & a_lph,
                      const EBLevelBoxData<CELL, 1> & a_phi)
  {
    m_finest->applyOpNeumann(a_lph, a_phi);
  }

private:

  /// 
  void residual(EBLevelBoxData<CELL, 1>       & a_res,
                const EBLevelBoxData<CELL, 1> & a_phi,
                const EBLevelBoxData<CELL, 1> & a_rhs);

  Real norm(const EBLevelBoxData<CELL, 1>       & a_res) const;
  

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
             const EBLevelBoxData<CELL, 1>& a_rhs,
             int a_niter)
  {
    m_finest->relax(a_phi, a_rhs, a_niter);
  }


  /// weak construction introduces unnecessary complications
  EBMultigrid()
  {
  }

  std::shared_ptr<EBMultigridLevel> m_finest;
  EBLevelBoxData<CELL, 1>           m_res;
  EBLevelBoxData<CELL, 1>           m_cor;

};
#include "Chombo_NamespaceFooter.H"

#endif
