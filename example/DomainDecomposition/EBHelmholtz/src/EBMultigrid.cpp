#include "EBMultigrid.H"
#include "Proto.H"
#include "Proto_Timer.H"
#include "NamespaceHeader.H"




int  EBMultigrid::s_numSmoothDown = 2;
int  EBMultigrid::s_numSmoothUp   = 2;

/****/
void
EBMultigridLevel::
applyOp(EBLevelBoxData<CELL, 1>       & a_lph,
        const EBLevelBoxData<CELL, 1> & a_phi)
                    
{
  EBLevelBoxData<CELL, 1>& phi = const_cast<EBLevelBoxData<CELL, 1>&>(a_phi);
  phi.exchange(m_exchangeCopier);
  DataIterator dit = m_grids.dataIterator();
  for(int ibox = 0; ibox < dit.size(); ++ibox)
  {
    string stenname("Second_Order_Poisson");

    shared_ptr<ebstencil_t> stencil = m_dictionary->getEBStencil(stenname, m_ebbcname, ibox);
    //set lphi = beta * div(F)
    stencil->apply(a_lph[dit[ibox]], a_phi[dit[ibox]], true, m_beta);
///HERE
///     ebforall(...);
  }
}
/****/
void
EBMultigrid::
applyOp(EBLevelBoxData<CELL, 1>       & a_lph,
        const EBLevelBoxData<CELL, 1> & a_phi)
{
  return m_finest->applyOp(a_lph, a_phi);
}
/***/
void
EBMultigrid::
residual(EBLevelBoxData<CELL, 1>       & a_res,
         const EBLevelBoxData<CELL, 1> & a_phi,
         const EBLevelBoxData<CELL, 1> & a_rhs)
{
  PR_TIME("sgmg::resid");
  return m_finest->residual(a_res, a_phi, a_rhs);
}
/***/
void
EBMultigrid::
vCycle(EBLevelBoxData<CELL, 1>       & a_phi,
       const EBLevelBoxData<CELL, 1> & a_rhs)
{
  PR_TIME("sgmg::vcycle");
  return m_finest->vCycle(a_phi, a_rhs);
}
/***/
void
EBMultigridLevel::
getMultiColors()
{
// Color offsets are grouped into "red"=even number of nonzeros (first 2^(DIM-1)) 
// and "black= odd number of nonzeros (the rest).
#if DIM==2
  m_colors[0] = Point::Zeros();//(0,0)
  m_colors[1] = Point::Ones();//(1,1)
  m_colors[2] = Point::Zeros() + Point::Basis(1);//(0,1)
  m_colors[3] = Point::Zeros() + Point::Basis(0);//(1,0)
#elif DIM==3
  m_colors[0] = Point::Zeros();//(0,0,0)
  m_colors[1] = Point::Zeros() + Point::Basis(0) + Point::Basis(1);//(1,1,0)
  m_colors[2] = Point::Zeros() + Point::Basis(1) + Point::Basis(2);//(0,1,1)
  m_colors[3] = Point::Zeros() + Point::Basis(0) + Point::Basis(2);//(1,0,1)
  m_colors[4] = Point::Zeros() + Point::Basis(1);//(0,1,0)
  m_colors[5] = Point::Zeros() + Point::Basis(0);//(1,0,0)
  m_colors[6] = Point::Zeros() + Point::Basis(2);//(0,0,1)
  m_colors[7] = Point::Ones();//(1,1,1)
#else
  compiler_error_this_code_is_only_written_for_dim_2_or_3();
#endif
}
/***/
EBMultigridLevel::
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
                 const IntVect                     & a_nghostdst)
{
  m_dictionary = a_dictionary; 
  m_alpha      = a_alpha;      
  m_beta       = a_beta;       
  m_dx         = a_dx;         
  m_grids      = a_grids;      
  m_stenname   = a_stenname;   
  m_dombcname  = a_dombcname;  
  m_ebbcname   = a_ebbcname;   
  m_domain     = a_domain;     
  m_nghostSrc  = a_nghostsrc;
  m_nghostDst  = a_nghostdst;
  defineStencils();

  defineCoarserObjects();
}
/***/
void
EBMultigridLevel::
defineCoarserObjects()
{
  PR_TIME("sgmglevel::defineCoarser");
}
/***/
EBMultigridLevel::
EBMultigridLevel(const EBMultigridLevel& a_finerLevel)
{
  PR_TIME("sgmglevel::constructor");

}
/***/
void
EBMultigridLevel::
defineStencils()
{
  PR_TIME("sgmglevel::definestencils");
  getMultiColors();
  m_exchangeCopier.exchangeDefine(m_grids, m_nghostSrc);
  //register stencil for apply op
  string stenname("Second_Order_Poisson");
  m_dictionary->registerStencil(stenname, m_dombcname, m_ebbcname);
}
/****/
void
EBMultigridLevel::
residual(EBLevelBoxData<CELL, 1>       & a_res,
         const EBLevelBoxData<CELL, 1> & a_phi,
         const EBLevelBoxData<CELL, 1> & a_rhs)
                    
{
}
/****/
void
EBMultigridLevel::
relax(EBLevelBoxData<CELL, 1>       & a_phi,
      const EBLevelBoxData<CELL, 1> & a_rhs)
{
  PR_TIME("sgmglevel::relax");
}
/****/
void
EBMultigridLevel::
restrictResidual(EBLevelBoxData<CELL, 1>       & a_resc,
                 const EBLevelBoxData<CELL, 1> & a_res)
{
  PR_TIME("sgmglevel::restrict");
}
/****/
void
EBMultigridLevel::
prolongIncrement(EBLevelBoxData<CELL, 1>      & a_phi,
                 const EBLevelBoxData<CELL, 1>& a_delta)
{
  PR_TIME("sgmglevel::prolong");
}
/****/
void 
EBMultigridLevel::
vCycle(EBLevelBoxData<CELL, 1>         & a_phi,
       const EBLevelBoxData<CELL, 1>   & a_rhs)
{

  PR_TIME("sgmglevel::vcycle");
  for(int irelax = 0; irelax < EBMultigrid::s_numSmoothDown; irelax++)
  {
    relax(a_phi,a_rhs); //don't do it
  }

  if (m_hasCoarser)
  {
    residual(m_resid,a_phi,a_rhs);                      
    m_coarser->restrictResidual(m_residC,m_resid);
    m_deltaC.setVal(0.);
    m_coarser->vCycle(m_deltaC,m_residC);
    m_coarser->prolongIncrement(a_phi,m_deltaC);
  }

  for(int irelax = 0; irelax < EBMultigrid::s_numSmoothUp; irelax++)
  {
    relax(a_phi,a_rhs);
  }

}
#include "NamespaceFooter.H"
/****/
