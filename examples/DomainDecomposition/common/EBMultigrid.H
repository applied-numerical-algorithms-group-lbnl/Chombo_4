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
    Users are strongly advised to just use the EBMultigrid API
*/
class EBMultigridLevel
{
public:
  friend class EBMultigrid;

  int  m_numSmooth;
  bool m_useWCycle     ;
  typedef CH4_Data_Choreography::DistributedData<EBGraph> graph_distrib_t;
  
  Real                              m_alpha;   
  Real                              m_beta;    
  std::shared_ptr<EBMultigridLevel> m_coarser;
  EBLevelBoxData<CELL, 1>           m_kappa;
  bool                              m_hasCoarser;
  string                            m_prefix;
  string                            m_bottom_solver;
  bool                              m_direct_to_bottom;
  EBLevelBoxData<CELL, 1>           m_resid;
  EBLevelBoxData<CELL, 1>           m_deltaC;
  EBLevelBoxData<CELL, 1>           m_residC;
  DisjointBoxLayout                 m_grids;   
  shared_ptr<graph_distrib_t>       m_graphs;
  IntVect                           m_nghost;
  Box                               m_domain;
  unsigned int                      m_depth;
  Real                              m_dx;      
  
  shared_ptr<GeometryService<2> >   m_geoserv;

/// Multigrid v-cycle.
  inline void vCycle(EBLevelBoxData<CELL, 1>      & a_phi,
                     const EBLevelBoxData<CELL, 1>& a_rhs,
                     bool a_printStuff = false)
  {

    PR_TIME("sgmglevel::vcycle");
    if(m_direct_to_bottom)
    {
      pout() << "going directly to bottom solve" << endl;
      bottom_solve(a_phi, a_rhs);
    }
    else
    {
      relax(a_phi,a_rhs, m_numSmooth);

      if (m_hasCoarser)
      {
        residual(m_resid,a_phi,a_rhs);                      
        //stencils for multilevel objects live with the coarser
        m_coarser->restrictResidual(m_residC,m_resid);
//begin debug
        {
          Chombo4::pout()  << setprecision(8)
                           << setiosflags(ios::showpoint)
                           << setiosflags(ios::scientific);
          
          EBIndex<CELL> vofmaxP;
          EBIndex<CELL> vofmaxR;
          Real maxP = a_phi.maxNorm(vofmaxP, 0);
          Real maxR = m_residC.maxNorm(vofmaxR, 0);
          Chombo4::pout() << "phi_max  = " << maxP << "@" << vofmaxP.m_pt << endl;
          Chombo4::pout() << "residC_max  = " << maxR << "@" << vofmaxR.m_pt << endl;
        }
//end debug        
        m_deltaC.setVal(0.);
        m_coarser->vCycle(m_deltaC,m_residC);
        if(m_useWCycle)
        {
          m_coarser->vCycle(m_deltaC,m_residC);
        }
        m_coarser->prolongIncrement(a_phi,m_deltaC);
//begin debug
        {
          EBIndex<CELL> vofmaxP;
          EBIndex<CELL> vofmaxD;
          Real maxP = a_phi.maxNorm(vofmaxP, 0);
          Real maxD = m_deltaC.maxNorm(vofmaxD, 0);
          Chombo4::pout() << "phi_max  = " << maxP << "@" << vofmaxP.m_pt << endl;
          Chombo4::pout() << "deC_max  = " << maxD << "@" << vofmaxD.m_pt << endl;
        }
//end debug        
      }
      else
      {
        bottom_solve(a_phi, a_rhs);
      }

      relax(a_phi,a_rhs, m_numSmooth);
    }
  }
  ///
  EBMultigridLevel()
  {
  }


  ///make coarse from fine
  virtual void 
  define(const EBMultigridLevel            & a_finerLevel,
         shared_ptr<GeometryService<2> >   & a_geoserv,
         bool a_printStuff) =0;

  ///
  virtual void
  preCond(EBLevelBoxData<CELL, 1>       & a_phi,
          const EBLevelBoxData<CELL, 1> & a_rhs) const =0;
  
  virtual void
  defineBottomSolvers(shared_ptr<GeometryService<2> >   & a_geoserv, bool a_printStuff) = 0;
  
  
  EBLevelBoxData<CELL, 1> & getKappa()
  {
    return m_kappa;
  }
  
  inline void resetAlphaAndBeta(const Real& a_alpha,
                                const Real& a_beta,
                                bool a_printStuff = false)
  {
    m_alpha = a_alpha;
    m_beta  = a_beta;
      
    if(m_hasCoarser)
    {
      m_coarser->resetAlphaAndBeta(a_alpha, a_beta);
    }
    if(!m_hasCoarser || m_direct_to_bottom)
    {
      this->defineBottomSolvers(m_geoserv, a_printStuff);
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
  typedef CH4_Data_Choreography::DistributedData<EBGraph> graph_distrib_t;
  
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
              string   a_bottom_solver,
              bool     a_direct_to_bottom,
              string a_prefix,
              bool a_useWCycle,
              int  a_numSmooth,
              bool a_printStuff = false);

  ///make coarse from fine 
  void
  define(const EBMultigridLevel            & a_finerLevel,
         shared_ptr<GeometryService<2> >   & a_geoserv,
         bool a_printStuff);
  
  ///
  void
  preCond(EBLevelBoxData<CELL, 1>       & a_phi,
          const EBLevelBoxData<CELL, 1> & a_rhs) const;

  
  ///  
  virtual ~EBPoissonOp()
  {
    m_hasRecoDataAlready  = false;
    m_hasRecoGraphAlready = false;
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

  void restrictResidualOnProc(EBLevelBoxData<CELL, 1>      & a_resc,
                              const EBLevelBoxData<CELL, 1>& a_res);
  void restrictResidualAgglom(EBLevelBoxData<CELL, 1>      & a_resc,
                              const EBLevelBoxData<CELL, 1>& a_res);
  
  void prolongIncrementOnProc(EBLevelBoxData<CELL, 1>      & a_phiFine,
                              const EBLevelBoxData<CELL, 1>& a_deltaCoarse);
  void prolongIncrementAgglom(EBLevelBoxData<CELL, 1>      & a_phiFine,
                              const EBLevelBoxData<CELL, 1>& a_deltaCoarse);

  /**
     For aggregation we sometimes (but not often)  need refined coarse data.
     This data is expensive to create in both memory and time so we have 
     mechanisms in place to only define it when it is needed.   
     That is what this dance is about.   When called the first time, these
     define member data on a refined version of the current grids.   
     It then returns the member data.    
     All subsequent calls just return the member data.
  **/
  shared_ptr<EBLevelBoxData<CELL, 1> >  getDataOnRefinedCoarseLayout()
  {
    PR_TIME("EBPoissonOp::getDataOnRefinedCoarseLayout");

    //we do not need this at every level so only construct reco data when needed
    if(!m_hasRecoDataAlready)
    {
      PR_TIME("data construction bit");
      shared_ptr<graph_distrib_t> graphsReCo = getGraphOnRefinedCoarseLayout();
      DisjointBoxLayout dblReCo = graphsReCo->disjointBoxLayout();

      typedef GraphConstructorFactory< EBBoxData<CELL, Real, 1> > devifactory_t;
      m_dataReCo = shared_ptr<EBLevelBoxData<CELL, 1> >(new EBLevelBoxData<CELL, 1>(dblReCo, m_nghost, graphsReCo));
      m_hasRecoDataAlready = true;
    }
  
    return m_dataReCo;
  }
  
  shared_ptr<graph_distrib_t>& getGraphOnRefinedCoarseLayout()
  {
    PR_TIME("EBPoissonOp::getGraphOnRefinedCoarseLayout");
    if(!m_hasRecoGraphAlready)
    {
      DisjointBoxLayout dblReCo;
      refine(dblReCo, m_grids, 2);
      m_graphReCo = shared_ptr<graph_distrib_t>(new graph_distrib_t(dblReCo, m_nghost, NullConstructorDataFactory<EBGraph >()));
      Box domfine = refine(m_domain, 2);

      DataIterator dit = m_grids.dataIterator();
      //need to fill graph first and exchange because the exchange can change the graph
      // and we want the data to be defined consistent with the final graph
      for(int ibox = 0; ibox < dit.size(); ibox++)
      {
        Box validBox = refine(m_grids[dit[ibox]], 2);
        Box grownBox = grow(validBox, m_nghost);
        Box interBoxCh = grownBox & domfine;
        Bx interBox = getProtoBox(interBoxCh);
        Bx domain   = getProtoBox(domfine);
        Bx validBx  = getProtoBox(validBox);

        (*m_graphReCo)[dit[ibox]].defineAsAllRegular(validBx, interBox,   domfine);
      }
    
      PR_assert(m_geoserv->hasLevel(domfine));
      shared_ptr<graph_distrib_t> graphFine = m_geoserv->getGraphs(domfine);
      graphFine->copyTo(*m_graphReCo);
      m_graphReCo->exchange();
      m_hasRecoGraphAlready = true;
    }
    return m_graphReCo;
  }

  void
  defineCoarserObjects(shared_ptr<GeometryService<2> >   & a_geoserv,
                       bool a_printStuff);
  
  shared_ptr<EBLevelBoxData<CELL, 1> > m_dataReCo;
  shared_ptr<graph_distrib_t>          m_graphReCo;
  bool                                 m_hasRecoDataAlready;
  bool                                 m_hasRecoGraphAlready;


//need the volume fraction in a data holder so we can evaluate kappa*alpha I 
  void  fillKappa(const shared_ptr<GeometryService<2> >   & a_geoserv,
                  bool a_printStuff = false);


  
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
  //for multilevel stuff I am bypassing the dictionary because the layouts often do not conform
  vector<shared_ptr<ebstencil_t> >  m_prolongationStencils[s_ncolors];
  vector<shared_ptr<ebstencil_t> >  m_restrictionStencil;
  
  string                            m_nobcname;
  EBLevelBoxData<CELL, 1>           m_diagW;
  
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

  void defineBottomSolvers(shared_ptr<GeometryService<2> >   & a_geoserv,
                           bool a_printStuff);
private:
  void createAndStoreRestrictionStencil();
  void createAndStoreProlongationStencils();
  void
  getRestrictionStencil(vector<EBIndex<CELL> >                    & a_dstVoFs,                    
                        vector<LocalStencil<CELL, double> >       & a_stencil,                    
                        const Box                                 & a_validCoar,
                        const EBGraph                             & a_graphCoar);

  void
  getProlongationStencil(vector<EBIndex<CELL> >                   & a_dstVoFs,                    
                         vector<LocalStencil<CELL, double> >      & a_stencil,                    
                         const Box                                & a_dstValid,                   
                         const EBGraph                            & a_dstGraph,                      
                         unsigned long                              a_icolor);               
};


///class that outsiders actually call
class EBMultigrid
{
public:
  typedef shared_ptr<EBDictionary<2, Real, CELL, CELL> > dictionary_t;

  /// Constant coefficient, named stencil
  /**
     When this constructor is called, it solves for the named stencil.
     Usually we use this for constant-coefficient Helmholtz.
   */
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
              string   a_bottom_solver,
              bool     a_direct_to_bottom,
              string   a_prefix,
              bool     a_useWCycle,
              int      a_numSmooth,
              bool a_printStuff = false)
  {
    string dombc[2*DIM];
    for(int ivec = 0; ivec < 2*DIM; ivec++)
    {
      dombc[ivec] = a_dombcname;
    }

    if(a_printStuff)
    {
      pout() << "EBMultigrid::EBMultigrid: making finest Poisson  operator" << endl;
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
                                              a_bottom_solver,
                                              a_direct_to_bottom,
                                              a_prefix,
                                              a_useWCycle,
                                              a_numSmooth,
                                              a_printStuff);
    EBMultigridLevel* baseBelongToUs = static_cast<EBMultigridLevel*>(finestfish);
    m_finest = std::shared_ptr<EBMultigridLevel>(baseBelongToUs);

    if(a_printStuff)
    {
      pout() << "EBMultigrid::EBMultigrid: making data holders" << endl;
    }
    auto graphs = a_geoserv->getGraphs(a_domain);
    m_res.define(a_grids, a_nghost, graphs);
    m_cor.define(a_grids, a_nghost, graphs);
    if(a_printStuff)
    {
      pout() << "EBMultigrid::EBMultigrid: leaving" << endl;
    }
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
              string   a_bottom_solver,
              bool     a_direct_to_bottom,
              string   a_prefix,
              bool     a_useWCycle,
              int      a_numSmooth,
              bool a_printStuff = false)
  {
    if(a_printStuff)
    {
      pout() << a_prefix << " EBMultigrid constructor: making finest operator" << endl; 
    }
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
                      a_bottom_solver,
                      a_direct_to_bottom,
                      a_prefix,
                      a_useWCycle,
                      a_numSmooth,
                      a_printStuff));

    if(a_printStuff)
    {
      pout() << a_prefix << " EBMultigrid constructor: making data holders" << endl; 
    }
    
    auto graphs = a_geoserv->getGraphs(a_domain);
    m_res.define(a_grids, a_nghost  , graphs);
    m_cor.define(a_grids, a_nghost  , graphs);
  }
  
  ///
  void 
  solve(EBLevelBoxData<CELL, 1>       & a_phi,
        const EBLevelBoxData<CELL, 1> & a_rhs,
        const Real                    & a_tolerance,
        const unsigned int            & a_maxIterations,
        bool a_initToZero = true,
        bool a_printStuff = false);

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
  inline
  void applyOp(EBLevelBoxData<CELL, 1>       & a_lph,
               const EBLevelBoxData<CELL, 1> & a_phi,
               bool a_doExchange = true) const
  {
    return m_finest->applyOp(a_lph, a_phi, a_doExchange);
  }


  /// for TGA force bcs to be homogeneous neumann
  inline
  void applyOpNeumann(EBLevelBoxData<CELL, 1>       & a_lph,
                      const EBLevelBoxData<CELL, 1> & a_phi) const
  {
    m_finest->applyOpNeumann(a_lph, a_phi);
  }

  inline
  void residual(EBLevelBoxData<CELL, 1>       & a_res,
                const EBLevelBoxData<CELL, 1> & a_phi,
                const EBLevelBoxData<CELL, 1> & a_rhs,
                bool a_doExchange = true) const
  {
    PR_TIME("sgmg::resid");
    return m_finest->residual(a_res, a_phi, a_rhs, a_doExchange);
  }
  
private:
  
  /// 

  Real norm(const EBLevelBoxData<CELL, 1>       & a_res) const;
  

  /// Multigrid v-cycle.
  inline
  void vCycle(EBLevelBoxData<CELL, 1>      & a_phi,
              const EBLevelBoxData<CELL, 1>& a_rhs,
              bool a_printStuff = false)
  {
    PR_TIME("sgmg::vcycle");
    return m_finest->vCycle(a_phi, a_rhs);
  }

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
