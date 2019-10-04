#include "EBMultigrid.H"
#include "Proto.H"
#include "Proto_Timer.H"
#include <sstream>
#include "NamespaceHeader.H"

int  EBMultigrid::s_numSmoothDown = 2;
int  EBMultigrid::s_numSmoothUp   = 2;

typedef Proto::Var<Real, 1> Sca;

/****/
void 
EBMultigrid::
solve(EBLevelBoxData<CELL, 1>       & a_phi,
      const EBLevelBoxData<CELL, 1> & a_rhs,
      const Real                    & a_tol,
      const unsigned int            & a_maxIter)
{
  EBLevelBoxData<CELL, 1>& res =  m_finest->m_resid;

  residual(res, a_phi, a_rhs);
  Real initres = res.maxNorm(0);
  int iter = 0;
  pout() << "EBMultigrid: tol = " << a_tol << ",  max iter = "<< a_maxIter << endl;
  Real resnorm = initres;
  while((iter < a_maxIter) && (resnorm > a_tol*initres))
  {
    pout() << setprecision(3)
           << setiosflags(ios::showpoint)
           << setiosflags(ios::scientific);

    pout() << "EBMultigrid: iter = " << iter << ", |resid| = " << resnorm << endl;

    vCycle(a_phi, a_rhs);

    residual(res, a_phi, a_rhs);
    resnorm = res.maxNorm(0);

    iter++;
  }
}
/****/
//lph comes in holding beta*div(F)--leaves holding alpha phi + beta div(F)
PROTO_KERNEL_START 
void  addAlphaPhiF(Sca     a_lph,
                   Sca     a_phi,
                   Sca     a_kap,
                   Real    a_alpha,
                   Real    a_beta)
{
  //kappa and beta are already in lph
  //kappa because we did not divide by kappa
  //beta was sent in to ebstencil::apply
  Real kappdiv = a_lph(0);
  Real bkdivF  = a_beta*kappdiv;
  Real phival  = a_phi(0);
  Real kapval  = a_kap(0);
  a_lph(0) = a_alpha*phival*kapval + bkdivF;
}
PROTO_KERNEL_END(addAlphaPhiF, addAlphaPhi)


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
    auto& phifab = a_phi[dit[ibox]];
    auto& lphfab = a_lph[dit[ibox]];
    shared_ptr<ebstencil_t> stencil = m_dictionary->getEBStencil(m_stenname, m_ebbcname, ibox);
    //set lphi = kappa* div(F)
    stencil->apply(lphfab, phifab,  true, 1.0);
    //this adds kappa*alpha*phi (making lphi = kappa*alpha*phi + kappa*beta*divF)
    unsigned long long int numflopspt = 3;
    Box grid =m_grids[dit[ibox]];
    Bx  grbx = getProtoBox(grid);
    ebforallInPlace(numflopspt, "addAlphaPhi", addAlphaPhi, grbx, a_lph[dit[ibox]], a_phi[dit[ibox]], m_kappa[dit[ibox]], m_alpha, m_beta);
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
EBMultigridLevel(dictionary_t                            & a_dictionary,
                 shared_ptr<GeometryService<2> >         & a_geoserv,
                 const Real                              & a_alpha,
                 const Real                              & a_beta,
                 const Real                              & a_dx,
                 const DisjointBoxLayout                 & a_grids,
                 const string                            & a_stenname,
                 const string                            & a_dombcname,
                 const string                            & a_ebbcname,
                 const Box                               & a_domain,
                 const IntVect                           & a_nghostsrc, 
                 const IntVect                           & a_nghostdst)
{
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
  m_dictionary = a_dictionary;
  const shared_ptr<LevelData<EBGraph>  > graphs = a_geoserv->getGraphs(m_domain);
  m_resid.define(m_grids, m_nghostSrc  , graphs);
  m_kappa.define(m_grids, IntVect::Zero, graphs);

  defineStencils(a_geoserv, graphs);

  defineCoarserObjects(a_geoserv);
}
/***/
void
EBMultigridLevel::
defineCoarserObjects(shared_ptr<GeometryService<2> >   & a_geoserv)
{
  PR_TIME("sgmglevel::defineCoarser");

  m_hasCoarser = (m_grids.coarsenable(4));
//  if(m_hasCoarser)
//  {
//    m_coarser = shared_ptr<EBMultigridLevel>(new EBMultigridLevel(*this, a_geoserv));
//
//    const shared_ptr<LevelData<EBGraph>  > graphs = a_geoserv->getGraphs(m_coarser->m_domain);
//    m_residC.define(m_coarser->m_grids, m_nghostSrc , graphs);
//    m_deltaC.define(m_coarser->m_grids, m_nghostDst , graphs);
//  }
}
/***/
EBMultigridLevel::
EBMultigridLevel(const EBMultigridLevel            & a_finerLevel,
                 shared_ptr<GeometryService<2> >   & a_geoserv)
{
  PR_TIME("sgmglevel::constructor");

  m_dx         = 2*a_finerLevel.m_dx;         
  m_domain     = coarsen(a_finerLevel.m_domain, 2);      
  coarsen(m_grids, a_finerLevel.m_grids,  2);      

  m_alpha      = a_finerLevel.m_alpha;      
  m_beta       = a_finerLevel.m_beta;       
  m_stenname   = a_finerLevel.m_stenname;   
  m_dombcname  = a_finerLevel.m_dombcname;  
  m_ebbcname   = a_finerLevel.m_ebbcname;   

  m_nghostSrc  = a_finerLevel.m_nghostSrc;
  m_nghostDst  = a_finerLevel.m_nghostDst;
  m_dictionary = a_finerLevel.m_dictionary;

  const shared_ptr<LevelData<EBGraph>  > graphs = a_geoserv->getGraphs(m_domain);
  m_resid.define(m_grids, m_nghostSrc  , graphs);
  m_kappa.define(m_grids, IntVect::Zero, graphs);

  defineStencils(a_geoserv, graphs);

  defineCoarserObjects(a_geoserv);

}
/***/
string convertUInt(unsigned int number)
{
  std::stringstream ss;//create a stringstream
  ss << number;//add number to the stream
  return ss.str();//return a string with the contents of the stream
}
void
EBMultigridLevel::
defineStencils(const shared_ptr<GeometryService<2> >   & a_geoserv,
               const shared_ptr<LevelData<EBGraph> >   & a_graphs)
{
  PR_TIME("sgmglevel::definestencils");

  m_exchangeCopier.exchangeDefine(m_grids, m_nghostSrc);
  //register stencil for apply op
  //true is for need the diagonal wweight
  m_dictionary->registerStencil(m_stenname, m_dombcname, m_ebbcname, true);
  m_nobcname = string("no_bcs");
  m_restrictionName = string("Multigrid_Restriction");
  m_dictionary->registerStencil(m_restrictionName, m_nobcname, m_nobcname, false);
///prolongation has ncolors stencils
  for(unsigned int icolor = 0; icolor < s_ncolors; icolor++)
  {
    string colorstring = "Multigrid_Prolongation_" + convertUInt(icolor);
    m_prolongationName[icolor] = colorstring;
    m_dictionary->registerStencil(colorstring, m_nobcname, m_nobcname, false);
  }
  //need the volume fraction in a data holder so we can evaluate kappa*alpha I 
  fillKappa(a_geoserv, a_graphs);
}
//need the volume fraction in a data holder so we can evaluate kappa*alpha I 
void  
EBMultigridLevel::
fillKappa(const shared_ptr<GeometryService<2> >   & a_geoserv,
          const shared_ptr<LevelData<EBGraph> >   & a_graphs)
{
  Real cellVol = 1;
  for(int idir = 0; idir < DIM; idir++)
  {
    cellVol *= m_dx;
  }
  typedef IndexedMoments<DIM  , 2> IndMomDIM;
  typedef HostIrregData<CELL,      IndMomDIM , 1>  VoluData;
  const shared_ptr<LevelData<VoluData> > volmomld = a_geoserv->getVoluData(m_domain);
  DataIterator dit = m_grids.dataIterator();
  for(int ibox = 0; ibox < dit.size(); ++ibox)
  {
    Box grid =m_grids[dit[ibox]];
    Bx  grbx = getProtoBox(grid);
    const EBGraph  & graph = (*a_graphs)[dit[ibox]];
    EBHostData<CELL, Real, 1> hostdat(grbx, graph);
    //fill kappa on the host then copy to the device
    HostBoxData<        Real, 1>& reghost = hostdat.getRegData();
    HostIrregData<CELL, Real, 1>& irrhost = hostdat.getIrrData();
    const VoluData & volmo = (*volmomld)[dit[ibox]];
    for(auto bit = grbx.begin(); bit != grbx.end();  ++bit)
    {
      if(graph.isRegular(*bit))
      {
        reghost(*bit, 0) = 1.0;
      }
      else if(graph.isCovered(*bit))
      {
        reghost(*bit, 0) = 0.0;
      }
      else
      {
        vector<EBIndex<CELL> > vofs = graph.getVoFs(*bit);
        for(int ivec = 0; ivec < vofs.size(); ivec++)
        {
          const EBIndex<CELL>& vof = vofs[ivec];
          const IndMomDIM&  momspt = volmo(vof, 0);
          double kappavof = momspt[0]/cellVol;

          reghost(*bit, 0) = kappavof;
          irrhost( vof, 0) = kappavof;
        }
      }
    }
    // now copy to the device
    EBLevelBoxData<CELL, 1>::copyToDevice(hostdat, m_kappa[dit[ibox]]);
  }
}
/****/
//res comes in holding lphi.   leaves holding res= rhs-lphi
PROTO_KERNEL_START 
void  subtractRHSF(Sca     a_res,
                   Sca     a_rhs)
{
  a_res(0) = a_rhs(0) - a_res(0);
}
PROTO_KERNEL_END(subtractRHSF, subtractRHS)
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
PROTO_KERNEL_START 
void  gsrbResidF(int     a_pt[DIM],
                 Sca     a_phi,
                 Sca     a_res,
                 Sca     a_diag,
                 Sca     a_kappa,
                 Real    a_alpha,
                 Real    a_beta,
                 Real    a_dx,
                 int     a_iredBlack)
{
  int sumpt = 0;
  for(int idir = 0; idir < DIM; idir++)
  {
    sumpt += a_pt[idir];
  }
  if(sumpt%2 == a_iredBlack)
  {
    static const Real safety = 1.0;
    Real diagval = a_diag(0);
    Real kappval = a_kappa(0);
    
    Real realdiag = kappval*a_alpha + a_beta*diagval;
    Real regudiag = a_alpha + 2*DIM*a_beta*a_dx*a_dx;
    Real lambda = safety/realdiag;
    Real reglam = safety/regudiag;
    Real phival = a_phi(0);
    Real resval = a_res(0);
    if(lambda > reglam) lambda= reglam;
    a_phi(0) = phival + lambda*resval;
  }
}
PROTO_KERNEL_END(gsrbResidF, gsrbResid)
/****/
void
EBMultigridLevel::
relax(EBLevelBoxData<CELL, 1>       & a_phi,
      const EBLevelBoxData<CELL, 1> & a_rhs)
{
  PR_TIME("sgmglevel::relax");
  //
  DataIterator dit = m_grids.dataIterator();
  for(int iredblack = 0; iredblack < 2; iredblack++)
  {
    Point debpt(26, 17);
    EBIndex<CELL> debvof(debpt, 0);

    residual(m_resid, a_phi, a_rhs);
    for(int ibox = 0; ibox < dit.size(); ++ibox)
    {
      shared_ptr<ebstencil_t>                  stencil  = m_dictionary->getEBStencil(m_stenname, m_ebbcname, ibox);
      shared_ptr< EBBoxData<CELL, Real, 1> >   diagptr  = stencil->getDiagonalWeights();

      const       EBBoxData<CELL, Real, 1> &   stendiag = *diagptr;
      Box grid = m_grids[dit[ibox]];
      Bx  grbx = getProtoBox(grid);
      //lambda = safety/diag
      //phi = phi - lambda*(res)
      ///lambda takes floating point to calculate
      //also does an integer check for red/black but I am not sure what to do with that
      auto& phifab =   a_phi[dit[ibox]];
      auto& resfab = m_resid[dit[ibox]];
      unsigned long long int numflopspt = 10;
      ebforallInPlace_i(numflopspt, "gsrbResid", gsrbResid,  grbx, 
                        phifab, resfab, stendiag,
                        m_kappa[dit[ibox]], m_alpha, m_beta, m_dx, iredblack);
      
      numflopspt++;
    }
  }
}
/****/
void
EBMultigridLevel::
restrictResidual(EBLevelBoxData<CELL, 1>       & a_resc,
                 const EBLevelBoxData<CELL, 1> & a_resf)
{
  PR_TIME("sgmglevel::restrict");
  DataIterator dit = m_grids.dataIterator();
  for(int ibox = 0; ibox < dit.size(); ++ibox)
  {
    auto& coarfab = a_resc[dit[ibox]];
    auto& finefab = a_resf[dit[ibox]];
    shared_ptr<ebstencil_t> stencil = m_dictionary->getEBStencil(m_stenname, m_nobcname, ibox);
    //set resc = Ave(resf) (true is initToZero)
    stencil->apply(coarfab, finefab,  true, 1.0);
  }
}
/****/
void
EBMultigridLevel::
prolongIncrement(EBLevelBoxData<CELL, 1>      & a_phi,
                 const EBLevelBoxData<CELL, 1>& a_cor)
{
  PR_TIME("sgmglevel::prolong");
  DataIterator dit = m_grids.dataIterator();
  for(int icolor = 0; icolor < s_ncolors; icolor++)
  {
    for(int ibox = 0; ibox < dit.size(); ++ibox)
    {
      auto& coarfab = a_cor[dit[ibox]];
      auto& finefab = a_phi[dit[ibox]];
      shared_ptr<ebstencil_t> stencil = m_dictionary->getEBStencil(m_prolongationName[icolor], m_nobcname, ibox);
      //phi  = phi + I(correction) (false means do not init to zero)
      stencil->apply(finefab, coarfab,  false, 1.0);
    }
  }
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
    relax(a_phi,a_rhs); 
  }

//  if (m_hasCoarser)
//  {
//    residual(m_resid,a_phi,a_rhs);                      
//    //  stencils live with the destination
//    m_coarser->restrictResidual(m_residC,m_resid);
//    m_deltaC.setVal(0.);
//    m_coarser->vCycle(m_deltaC,m_residC);
//    //   stencils live with the destination
//    prolongIncrement(a_phi,m_deltaC);
//  }

  for(int irelax = 0; irelax < EBMultigrid::s_numSmoothUp; irelax++)
  {
    relax(a_phi,a_rhs);
  }

}
#include "NamespaceFooter.H"
/****/
