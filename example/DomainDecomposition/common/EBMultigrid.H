#ifndef _EBMultigrid_H_
#define _EBMultigrid_H_
#include <cmath>
#include <memory>
#include "Proto.H"
#include "Chombo_EBDictionary.H"
#include "Chombo_EBChombo.H"
#include "Chombo_EBLevelBoxData.H"
#include "EBPetscSolver.H"
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
/** This not a design that conforms to any pattern I know.. 
    This is not designed with anything but code reuse in mind.   
    Users are strong advised to just use the EBMultigrid API
*/
class EBMultigridLevel
{
public:
  friend class EBMultigrid;

  static int  s_numSmoothDown ;
  static int  s_numSmoothUp   ;
  static bool s_useWCycle     ;
  
  Real                              m_alpha;   
  Real                              m_beta;    
  std::shared_ptr<EBMultigridLevel> m_coarser;
  EBLevelBoxData<CELL, 1>           m_kappa;
  bool                              m_hasCoarser;
  bool                              m_directToBottom;
  EBLevelBoxData<CELL, 1>           m_resid;
  EBLevelBoxData<CELL, 1>           m_deltaC;
  EBLevelBoxData<CELL, 1>           m_residC;
  DisjointBoxLayout                 m_grids;   
  shared_ptr<LevelData<EBGraph> >   m_graphs;
  IntVect                           m_nghost;
  Box                               m_domain;
  unsigned int                      m_depth;
  Real                              m_dx;      
  
  /// Multigrid v-cycle.
  inline void vCycle(EBLevelBoxData<CELL, 1>      & a_phi,
                     const EBLevelBoxData<CELL, 1>& a_rhs)
  {

    PR_TIME("sgmglevel::vcycle");
    if(m_directToBottom)
    {
      bottom_solve(a_phi, a_rhs);
    }
    else
    {
      relax(a_phi,a_rhs, EBMultigridLevel::s_numSmoothDown); 

      if (m_hasCoarser)
      {
        residual(m_resid,a_phi,a_rhs);                      
        //stencils for multilevel objects live with the finer level
        restrictResidual(m_residC,m_resid);
        m_deltaC.setVal(0.);
        m_coarser->vCycle(m_deltaC,m_residC);
        if(EBMultigridLevel::s_useWCycle)
        {
          m_coarser->vCycle(m_deltaC,m_residC);
        }
        prolongIncrement(a_phi,m_deltaC);
      }
      else
      {
        bottom_solve(a_phi, a_rhs);
      }

      relax(a_phi,a_rhs, EBMultigridLevel::s_numSmoothUp);
    }
  }
  ///
  EBMultigridLevel()
  { }

  ///make coarse from fine
  virtual void 
  define(const EBMultigridLevel            & a_finerLevel,
         shared_ptr<GeometryService<2> >   & a_geoserv) =0;

  ///
  virtual void
  preCond(EBLevelBoxData<CELL, 1>       & a_phi,
          const EBLevelBoxData<CELL, 1> & a_rhs) const =0;

  
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
  
  ///  
  virtual ~EBMultigridLevel()
  {
  }

  ///
  virtual void
  residual(EBLevelBoxData<CELL, 1>       & a_res,
           const EBLevelBoxData<CELL, 1> & a_phi,
           const EBLevelBoxData<CELL, 1> & a_rhs,
           bool a_doExchange = true) const=0;


  /// 
  virtual void
  relax(EBLevelBoxData<CELL, 1>      & a_phi,
        const EBLevelBoxData<CELL, 1>& a_rhs,
        int a_niter) const =0;


  /// average down residual to next coarser level.
  virtual void
  restrictResidual(EBLevelBoxData<CELL, 1>      & a_resc,
                   const EBLevelBoxData<CELL, 1>& a_res)=0;

  /// Piecewise constant interpolation of coarse correction to increment fine solution.
  virtual void
  prolongIncrement(EBLevelBoxData<CELL, 1>      & a_phiFine,
                   const EBLevelBoxData<CELL, 1>& a_deltaCoarse)=0;


  /// 
  virtual void
  applyOp(EBLevelBoxData<CELL, 1>       & a_lph,
          const EBLevelBoxData<CELL, 1> & a_phi,
          bool a_doExchange = true) const=0;

  ///for tga
  virtual void
  applyOpNeumann(EBLevelBoxData<CELL, 1>       & a_lph,
                 const EBLevelBoxData<CELL, 1> & a_phi) const=0;


  virtual void
  bottom_solve(EBLevelBoxData<CELL, 1>      &  a_phi,
               const EBLevelBoxData<CELL, 1>&  a_rhs) = 0;
  
  // needed for bicgstab
  inline void create(EBLevelBoxData<CELL, 1>&  a_dataholder) const
  {
    a_dataholder.define(m_grids, m_nghost, m_graphs);
  }
  // needed for bicgstab
  inline void setToZero(EBLevelBoxData<CELL, 1>&  a_dataholder) const
  {
    a_dataholder.setVal(0.);
  }

  // needed for bicgstab
  inline void assignLocal(EBLevelBoxData<CELL, 1>&        a_dst,
                          const EBLevelBoxData<CELL, 1>&  a_src) const
  {
    DataIterator dit = m_grids.dataIterator();
    for(int ibox = 0; ibox < dit.size(); ++ibox)
    {
      a_dst[dit[ibox]].setVal(0.0);
      a_dst[dit[ibox]] += a_src[dit[ibox]] ;
    }
  }
  // needed for bicgstab
  inline Real norm(const EBLevelBoxData<CELL, 1> & a_src) const
  {
    return a_src.maxNorm(0);
  }

  // needed for bicgstab
  inline Real dotProduct(const EBLevelBoxData<CELL, 1> & a_srcOne,
                         const EBLevelBoxData<CELL, 1> & a_srcTwo) const
  {
    EBLevelBoxData<CELL, 1>  dst(m_grids, a_srcOne.ghostVect(), m_graphs);
    DataIterator dit = m_grids.dataIterator();
    for(int ibox = 0; ibox < dit.size(); ++ibox)
    {
      dst[dit[ibox]].setVal(1.0);
      dst[dit[ibox]] *= a_srcOne[dit[ibox]] ;
      dst[dit[ibox]] *= a_srcTwo[dit[ibox]] ;
    }
    return dst.sum(0);
  }
  // needed for bicgstab
  inline void incr(EBLevelBoxData<CELL, 1>       & a_dst,
                   const EBLevelBoxData<CELL, 1> & a_src,
                   Real a_scale) const
  {
    DataIterator dit = m_grids.dataIterator();
    for(int ibox = 0; ibox < dit.size(); ++ibox)
    {
      auto      & dstfab = a_dst[dit[ibox]];
      const auto& srcfab = a_src[dit[ibox]];
      dstfab.incr(srcfab, a_scale);
    }
  }
  // needed for bicgstab
  inline void scale(EBLevelBoxData<CELL, 1>       & a_dst,
                     Real a_scale) const
  {
    DataIterator dit = m_grids.dataIterator();
    for(int ibox = 0; ibox < dit.size(); ++ibox)
    {
      auto      & dstfab = a_dst[dit[ibox]];
      dstfab.scale(a_scale);
    }
  }


};




class EBPoissonOp: public EBMultigridLevel
{
public:
  friend class EBMultigrid;

  typedef EBStencil<2, Real, CELL, CELL> ebstencil_t;
  typedef shared_ptr<EBDictionary<2, Real, CELL, CELL> > dictionary_t;
  
  ///
  EBPoissonOp(dictionary_t                      & a_dictionary,
                   shared_ptr<GeometryService<2> >   & a_geoserv,
                   const Real                        & a_alpha,
                   const Real                        & a_beta,
                   const Real                        & a_dx,
                   const DisjointBoxLayout           & a_grids,
                   const string                      & a_stenname,
                   string                              a_dombcname[2*DIM],
                   const string                      & a_ebbcname,
                   const Box                         & a_domain,
                   const IntVect                     & a_nghost,
                   bool a_directToBottomSolve = false);

  ///make coarse from fine 
  void
  define(const EBMultigridLevel            & a_finerLevel,
         shared_ptr<GeometryService<2> >   & a_geoserv);
  
  ///
  void
  preCond(EBLevelBoxData<CELL, 1>       & a_phi,
          const EBLevelBoxData<CELL, 1> & a_rhs) const;

  
  ///  
  virtual ~EBPoissonOp()
  {
  }

  /// 
  void residual(EBLevelBoxData<CELL, 1>       & a_res,
                const EBLevelBoxData<CELL, 1> & a_phi,
                const EBLevelBoxData<CELL, 1> & a_rhs,
                bool a_doExchange = true) const;


  /// 
  void relax(EBLevelBoxData<CELL, 1>      & a_phi,
             const EBLevelBoxData<CELL, 1>& a_rhs,
             int a_niter) const;


  /// average down residual to next coarser level.
  void restrictResidual(EBLevelBoxData<CELL, 1>      & a_resc,
                        const EBLevelBoxData<CELL, 1>& a_res);

  /// Piecewise constant interpolation of coarse correction to increment fine solution.
  void prolongIncrement(EBLevelBoxData<CELL, 1>      & a_phiFine,
                        const EBLevelBoxData<CELL, 1>& a_deltaCoarse);


  /// 
  void applyOp(EBLevelBoxData<CELL, 1>       & a_lph,
               const EBLevelBoxData<CELL, 1> & a_phi,
               bool a_doExchange = true) const;

  ///for tga
  void applyOpNeumann(EBLevelBoxData<CELL, 1>       & a_lph,
                      const EBLevelBoxData<CELL, 1> & a_phi) const;



  void bottom_solve(EBLevelBoxData<CELL, 1>      &  a_phi,
                    const EBLevelBoxData<CELL, 1>&  a_rhs);
  /// weak construction introduces unnecessary complications
  /** but it sometimes proves useful in object-oriented code 
      nobody should be making one of these outside the EBMultigrid context
      so this unprotected weak construction is probably okay.
      Blame Stroustrup for outlawing virtual constructors. */
  EBPoissonOp()
  {
  }
private:

  void defineCoarserObjects(shared_ptr<GeometryService<2> >   & a_geoserv);


//need the volume fraction in a data holder so we can evaluate kappa*alpha I 
  void  fillKappa(const shared_ptr<GeometryService<2> >   & a_geoserv);


  
  string                            m_stenname;
  string                            m_neumname;
  string                            m_dombcname[2*DIM];
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
  EBLevelBoxData<CELL, 1>           m_diagW;

  Copier                            m_exchangeCopier;
  shared_ptr<EBRelaxSolver>         m_relaxSolver;
#ifdef CH_USE_PETSC  
  shared_ptr<EBPetscSolver<2> >     m_petscSolver;

  void solve_petsc(EBLevelBoxData<CELL, 1>       & a_phi,
                   const EBLevelBoxData<CELL, 1> & a_rhs);
#endif
  void solve_relax(EBLevelBoxData<CELL, 1>       & a_phi,
                   const EBLevelBoxData<CELL, 1> & a_rhs);
  
  void solve_bicgstab(EBLevelBoxData<CELL, 1>       & a_phi,
                      const EBLevelBoxData<CELL, 1> & a_rhs);

  void defineBottomSolvers(shared_ptr<GeometryService<2> >   & a_geoserv);
};


///class that outsiders actually call
class EBMultigrid
{
public:
  typedef shared_ptr<EBDictionary<2, Real, CELL, CELL> > dictionary_t;

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
              const IntVect                     & a_nghost,
              bool a_directToBottom = false)
  {
    string dombc[2*DIM];
    for(int ivec = 0; ivec < 2*DIM; ivec++)
    {
      dombc[ivec] = a_dombcname;
    }
    EBPoissonOp* finestfish = new EBPoissonOp(a_dictionary, 
                                              a_geoserv,
                                              a_alpha,      
                                              a_beta,       
                                              a_dx,         
                                              a_grids,      
                                              a_stenname,   
                                              dombc,
                                              a_ebbcname,   
                                              a_domain, 
                                              a_nghost,
                                              a_directToBottom);
    EBMultigridLevel* baseBelongToUs = static_cast<EBMultigridLevel*>(finestfish);
    m_finest = std::shared_ptr<EBMultigridLevel>(baseBelongToUs);

    const shared_ptr<LevelData<EBGraph>  > graphs = a_geoserv->getGraphs(a_domain);
    m_res.define(a_grids, a_nghost, graphs);
    m_cor.define(a_grids, a_nghost, graphs);
  }


  /// 
  EBMultigrid(dictionary_t                      & a_dictionary,
              shared_ptr<GeometryService<2> >   & a_geoserv,
              const Real                        & a_alpha,
              const Real                        & a_beta,
              const Real                        & a_dx,
              const DisjointBoxLayout           & a_grids,
              const string                      & a_stenname,
              string                              a_dombcname[2*DIM],
              const string                      & a_ebbcname,
              const Box                         & a_domain,
              const IntVect                     & a_nghost,
              bool a_directToBottomSolve = false)
  {
    m_finest = std::shared_ptr<EBPoissonOp>(
      new EBPoissonOp(a_dictionary, 
                           a_geoserv,
                           a_alpha,      
                           a_beta,       
                           a_dx,         
                           a_grids,      
                           a_stenname,   
                           a_dombcname,  
                           a_ebbcname,   
                           a_domain, 
                           a_nghost,
                           a_directToBottomSolve));

    const shared_ptr<LevelData<EBGraph>  > graphs = a_geoserv->getGraphs(a_domain);
    m_res.define(a_grids, a_nghost  , graphs);
    m_cor.define(a_grids, a_nghost  , graphs);
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
               const EBLevelBoxData<CELL, 1> & a_phi,
               bool a_doExchange = true) const;


  /// for TGA force bcs to be homogeneous neumann
  inline
  void applyOpNeumann(EBLevelBoxData<CELL, 1>       & a_lph,
                      const EBLevelBoxData<CELL, 1> & a_phi) const
  {
    m_finest->applyOpNeumann(a_lph, a_phi);
  }

  void residual(EBLevelBoxData<CELL, 1>       & a_res,
                const EBLevelBoxData<CELL, 1> & a_phi,
                const EBLevelBoxData<CELL, 1> & a_rhs,
                bool a_doExchange = true) const;
  
private:
  
  /// 

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
             int a_niter) const
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
