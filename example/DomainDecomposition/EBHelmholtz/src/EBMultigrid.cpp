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
                   Sca     a_kap,
                   Real    a_alpha)
{
  //kappa and beta are already in lph
  //kappa because we did not divide by kappa
  //beta was sent in to ebstencil::apply
  Real bkdivF = a_lph(0);
  Real phival  = a_phi(0);
  Real kapval  = a_kap(0);
  a_lph(0) = a_alpha*phival*kapval + bkdivF;
}
PROTO_KERNEL_END(addAlphaPhiF, addAlphaPhi)


/****/
PROTO_KERNEL_START 
void  addAlphaPhiPtF(int a_pt[DIM],
                     Sca     a_lph,
                     Sca     a_phi,
                     Sca     a_kap,
                     Real    a_alpha)
{
  //kappa and beta are already in lph
  //kappa because we did not divide by kappa
  //beta was sent in to ebstencil::apply
  Real bkdivF = a_lph(0);
  Real phival  = a_phi(0);
  Real kapval  = a_kap(0);
  if((a_pt[0] == 26) && (a_pt[1]==16))
  {
    Point debpt(a_pt[0], a_pt[1]);
//    cout << debpt << "addphipt: bkdivf = " << bkdivF << endl;
  }
  a_lph(0) = a_alpha*phival*kapval + bkdivF;
}
PROTO_KERNEL_END(addAlphaPhiPtF, addAlphaPhiPt)
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
    //set lphi = kappa*beta * div(F)
    stencil->apply(a_lph[dit[ibox]], a_phi[dit[ibox]], true, m_beta);
    //this adds kappa*alpha*phi (making lphi = kappa*alpha*phi + kappa*beta*divF)
    unsigned long long int numflopspt = 3;
    Box grid =m_grids[dit[ibox]];
    Bx  grbx = getProtoBox(grid);
    //ebforallInPlace(numflopspt, "addAlphaPhi", addAlphaPhi, grbx, a_lph[dit[ibox]], a_phi[dit[ibox]], m_kappa[dit[ibox]], m_alpha);
    ebforallInPlace_i(numflopspt, "addAlphaPhi", addAlphaPhiPt, grbx, a_lph[dit[ibox]], a_phi[dit[ibox]], m_kappa[dit[ibox]], m_alpha);
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
  const shared_ptr<LevelData<EBGraph>  > graphs = a_geoserv->getGraphs(m_domain);
  m_resid.define(m_grids, m_nghostSrc  , graphs);
  m_kappa.define(m_grids, IntVect::Zero, graphs);

  defineStencils(a_geoserv, graphs);

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
defineStencils(const shared_ptr<GeometryService<2> >   & a_geoserv,
               const shared_ptr<LevelData<EBGraph> >   & a_graphs)
{
  PR_TIME("sgmglevel::definestencils");

  m_exchangeCopier.exchangeDefine(m_grids, m_nghostSrc);
  //register stencil for apply op
  m_dictionary->registerStencil(m_stenname, m_dombcname, m_ebbcname);
  Real cellVol = 1;
  for(int idir = 0; idir < DIM; idir++)
  {
    cellVol *= m_dx;
  }
  //need the volume fraction in a data holder so we can evaluate kappa*alpha I +... HERE
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
//res comes in holding lphi.   leaves holding res= rhs-lphi
PROTO_KERNEL_START 
void  subtractRHSFptF(int     a_pt[DIM],
                      Sca     a_res,
                      Sca     a_rhs)
{
  if((a_pt[0] == 26) && (a_pt[1]==16))
  {
    Point debpt(a_pt[0], a_pt[1]);
//      cout << debpt << "grb: devi diag = " << diagval << endl;
  }
  Real resval = a_res(0);
  Real rhsval = a_rhs(0);
  a_res(0) = a_rhs(0) - a_res(0);
  a_res(0) = rhsval- resval;
}
PROTO_KERNEL_END(subtractRHSFptF, subtractRHSpt)
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
//    ebforallInPlace(numflopspt, "subtractRHS", subtractRHS,  grbx, a_res[dit[ibox]], a_rhs[dit[ibox]]);
    ebforallInPlace_i(numflopspt, "subtractRHSpt", subtractRHSpt,  grbx, a_res[dit[ibox]], a_rhs[dit[ibox]]);
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
    if((a_pt[0] == 26) && (a_pt[1]==16))
    {
      Point debpt(a_pt[0], a_pt[1]);
//      cout << debpt << "grb: devi diag = " << diagval << endl;
    }
    
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
    Point debpt(26, 16);
    EBIndex<CELL> debvof(debpt, 0);
    for(int ibox = 0; ibox < dit.size(); ++ibox)
    {
      auto& phifab = a_phi[dit[ibox]];
      auto& rhsfab = a_rhs[dit[ibox]];
      cout << "gsrb before redblack = " << iredblack << ", phi(26,16) = "<< phifab(debvof, 0) << endl;
    }

    residual(m_resid, a_phi, a_rhs);
    for(int ibox = 0; ibox < dit.size(); ++ibox)
    {
      shared_ptr<ebstencil_t>                  stencil  = m_dictionary->getEBStencil(m_stenname, m_ebbcname, ibox);
      shared_ptr< EBBoxData<CELL, Real, 1> >   diagptr  = stencil->getDiagonalWeights();
//
//      cout << debpt << "rel: devi diag = " << (*diagptr)(debvof, 0) << endl;
      
//      
      const       EBBoxData<CELL, Real, 1> &   stendiag = *diagptr;
      Box grid = m_grids[dit[ibox]];
      Bx  grbx = getProtoBox(grid);
      //lambda = safety/diag
      //phi = phi - lambda*(res)
      ///lambda takes floating point to calculate
      //also does an integer check for red/black but I am not sure what to do with that
      unsigned long long int numflopspt = 10;
      ebforallInPlace_i(numflopspt, "gsrbResid", gsrbResid,  grbx, 
                        a_phi[dit[ibox]], m_resid[dit[ibox]], stendiag,
                        m_kappa[dit[ibox]], m_alpha, m_beta, m_dx, iredblack);
      
      cout << "gsrb after redblack = " << iredblack << ", phi(26,16) = "<< a_phi[dit[ibox]](debvof, 0) << endl;
      numflopspt++;
    }
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

//  if (m_hasCoarser)
//  {
//    residual(m_resid,a_phi,a_rhs);                      
//    m_coarser->restrictResidual(m_residC,m_resid);
//    m_deltaC.setVal(0.);
//    m_coarser->vCycle(m_deltaC,m_residC);
//    m_coarser->prolongIncrement(a_phi,m_deltaC);
//  }

  for(int irelax = 0; irelax < EBMultigrid::s_numSmoothUp; irelax++)
  {
    relax(a_phi,a_rhs);
  }

}
#include "NamespaceFooter.H"
/****/
