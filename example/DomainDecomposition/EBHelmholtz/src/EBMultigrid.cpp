#include "EBMultigrid.H"
#include "Proto.H"
#include "Proto_Timer.H"
#include "NamespaceHeader.H"



int  EBMultigrid::s_numSmoothDown = 2;
int  EBMultigrid::s_numSmoothUp   = 2;

typedef Proto::Var<Real, 1> Sca;

/****/
//lph comes in holding beta*div(F)--leaves holding alpha phi + beta div(F)
PROTO_KERNEL_START 
void  addAlphaPhiF(Sca     a_lph,
                   Sca     a_phi,
                   Real    a_alpha)
{
  Real betadivF = a_lph(0);

  a_lph(0) = a_alpha*a_phi(0) + betadivF;
}
PROTO_KERNEL_END(addAlphaPhiF, addAlphaPhi)

//res comes in holding lphi.   leaves holding res= lphi-rhs
PROTO_KERNEL_START 
void  subtractRHSF(Sca     a_res,
                   Sca     a_phi)
{
  a_res(0) = a_res(0) - a_phi(0);
}
PROTO_KERNEL_END(subtractRHSF, subtractRHS)

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
    shared_ptr<ebstencil_t> stencil = m_dictionary->getEBStencil(m_stenname, m_ebbcname, ibox);
    //set lphi = beta * div(F)
    stencil->apply(a_lph[dit[ibox]], a_phi[dit[ibox]], true, m_beta);
    //this adds alpha*phi (making lphi = alpha*phi + beta*divF)
    unsigned long long int numflopspt = 2;
    Box grid =m_grids[dit[ibox]];
    Bx  grbx = getProtoBox(grid);
    ebforallInPlace(numflopspt, "addAlphaPhi", addAlphaPhi, grbx, a_lph[dit[ibox]], a_phi[dit[ibox]], m_alpha);
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

  m_exchangeCopier.exchangeDefine(m_grids, m_nghostSrc);
  //register stencil for apply op
  m_dictionary->registerStencil(m_stenname, m_dombcname, m_ebbcname);
}
/****/
void
EBMultigridLevel::
residual(EBLevelBoxData<CELL, 1>       & a_res,
         const EBLevelBoxData<CELL, 1> & a_phi,
         const EBLevelBoxData<CELL, 1> & a_rhs)
                    
{
  //this puts lphi into a_res
  applyOp(a_res, a_phi);
  //subtract off rhs so res = lphi - rhs
  DataIterator dit = m_grids.dataIterator();
  for(int ibox = 0; ibox < dit.size(); ++ibox)
  {
    //this adds alpha*phi (making lphi = alpha*phi + beta*divF)
    unsigned long long int numflopspt = 2;
    Box grid = m_grids[dit[ibox]];
    Bx  grbx = getProtoBox(grid);
    ebforallInPlace(numflopspt, "subtractRHS", subtractRHS,  grbx, a_res[dit[ibox]], a_rhs[dit[ibox]]);
  }
}
/****/
void
EBMultigridLevel::
relax(EBLevelBoxData<CELL, 1>       & a_phi,
      const EBLevelBoxData<CELL, 1> & a_rhs)
{
  PR_TIME("sgmglevel::relax");
  residual(m_resid, a_phi, a_rhs);
  //
  DataIterator dit = m_grids.dataIterator();
  for(int ibox = 0; ibox < dit.size(); ++ibox)
  {
    //this adds alpha*phi (making lphi = alpha*phi + beta*divF)
    unsigned long long int numflopspt = 2;
    Box grid = m_grids[dit[ibox]];
    Bx  grbx = getProtoBox(grid);
    ebforallInPlace(numflopspt, "subtractRHS", subtractRHS,  grbx, a_res[dit[ibox]], a_rhs[dit[ibox]]);
  }
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
