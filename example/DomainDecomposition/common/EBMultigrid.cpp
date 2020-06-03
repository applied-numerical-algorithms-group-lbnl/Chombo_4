#include "EBMultigrid.H"
#include "Proto.H"
#include "Proto_Timer.H"
#include "EBRelaxSolver.H"
#include "EBMultigridFunctions.H"
#include <sstream>
#include "BiCGStabSolver.H"
#include "Chombo_ParmParse.H"
#include "Chombo_NamespaceHeader.H"

bool EBMultigrid::s_useWCycle     = true;
int  EBMultigrid::s_numSmoothDown = 4;
int  EBMultigrid::s_numSmoothUp   = 4;

/****/
void 
EBMultigrid::
solve(EBLevelBoxData<CELL, 1>       & a_phi,
      const EBLevelBoxData<CELL, 1> & a_rhs,
      const Real                    & a_tol,
      const unsigned int            & a_maxIter,
      bool a_initToZero)
{
  CH_TIME("EBMultigrid::solve");
  if(a_initToZero)
  {
    a_phi.setVal(0.);
  }

  residual(m_res, a_phi, a_rhs, true);
  Real initres = m_res.maxNorm(0);
  int iter = 0;
  pout() << "EBMultigrid: tol = " << a_tol << ",  max iter = "<< a_maxIter << endl;
  Real resnorm = initres;
  Real resnormold = resnorm;
  while((iter < a_maxIter) && (resnorm > a_tol*initres))
  {
    m_cor.setVal(0.);
    pout() << setprecision(3)
           << setiosflags(ios::showpoint)
           << setiosflags(ios::scientific);

    
    pout() << "EBMultigrid: iter = " << iter << ", |resid| = " << resnorm;
    Real rate = 1;
    if((resnormold > 1.0e-12) && (iter > 0))
    {
      rate = resnormold/resnorm;
      pout() << ", rate = " << rate;
    }
    pout() << endl;
    
    vCycle(m_cor, m_res);

    a_phi += m_cor;
    residual(m_res, a_phi, a_rhs, true);

    resnormold = resnorm;
    resnorm = m_res.maxNorm(0);

    iter++;
  }
  pout() << "EBMultigrid: final     |resid| = " << resnorm << endl;
}

/****/
void
EBMultigridLevel::
applyOp(EBLevelBoxData<CELL, 1>       & a_lph,
        const EBLevelBoxData<CELL, 1> & a_phi,
        bool a_doExchange) const
                    
{
  CH_TIME("EBMultigrid::applyOp");
  CH_assert(  a_lph.ghostVect() == a_phi.ghostVect());
  CH_assert(m_kappa.ghostVect() == a_phi.ghostVect());
  EBLevelBoxData<CELL, 1>& phi = const_cast<EBLevelBoxData<CELL, 1>&>(a_phi);
  if(a_doExchange)
  {
    phi.exchange(m_exchangeCopier);
  }
  
  DataIterator dit = m_grids.dataIterator();
  int ideb  = 0;
  for(int ibox = 0; ibox < dit.size(); ++ibox)
  {
    auto& phifab = a_phi[dit[ibox]];
    auto& lphfab = a_lph[dit[ibox]];
    shared_ptr<ebstencil_t> stencil =
      m_dictionary->getEBStencil(m_stenname, m_ebbcname, m_domain, m_domain, ibox);
    //set lphi = kappa* div(F)
//    Bx lphbox = lphfab.box();
//    Bx phibox = phifab.box();
    stencil->apply(lphfab, phifab,  true, 1.0);
    //this adds kappa*alpha*phi (making lphi = kappa*alpha*phi + kappa*beta*divF)
    unsigned long long int numflopspt = 3;
    Box grid =m_grids[dit[ibox]];
    Bx  grbx = getProtoBox(grid);
    auto& kapfab = m_kappa[dit[ibox]];
#if 1
    ebforallInPlace(numflopspt, "addAlphaPhi", addAlphaPhi, grbx, lphfab, phifab, kapfab, m_alpha, m_beta);
#else
//    ebFastforallInPlace(numflopspt, "addAlphaPhi", addAlphaPhi, grbx, lphfab, phifab, kapfab, m_alpha, m_beta);
    ebFastforallInPlace_i(phifab.inputBox(), numflopspt, "addAlphaPhiPt", addAlphaPhiPt, grbx, lphfab, phifab, kapfab, m_alpha, m_beta);
#endif
    ideb++;
  }
}
/****/
void
EBMultigridLevel::
preCond(EBLevelBoxData<CELL, 1>       & a_phi,
        const EBLevelBoxData<CELL, 1> & a_rhs) const
{
  CH_TIME("EBMultigrid::preCond");
  int maxiter = 27;

  relax(a_phi,a_rhs, maxiter); 
}
/****/
//for tga
void
EBMultigridLevel::
applyOpNeumann(EBLevelBoxData<CELL, 1>       & a_lph,
               const EBLevelBoxData<CELL, 1> & a_phi) const
                    
{
  CH_assert(  a_lph.ghostVect() == a_phi.ghostVect());
  CH_assert(m_kappa.ghostVect() == a_phi.ghostVect());
  CH_TIME("EBMultigrid::applyOpNeumann");
  EBLevelBoxData<CELL, 1>& phi = const_cast<EBLevelBoxData<CELL, 1>&>(a_phi);
  phi.exchange(m_exchangeCopier);
  DataIterator dit = m_grids.dataIterator();
  int ideb  = 0;
  for(int ibox = 0; ibox < dit.size(); ++ibox)
  {
    auto& phifab = a_phi[dit[ibox]];
    auto& lphfab = a_lph[dit[ibox]];
    shared_ptr<ebstencil_t> stencil =
      m_dictionary->getEBStencil(m_neumname, StencilNames::Neumann, m_domain, m_domain, ibox);
    //set lphi = kappa* div(F)
//    Bx lphbox = lphfab.box();
//    Bx phibox = phifab.box();
    stencil->apply(lphfab, phifab,  true, 1.0);
    //this adds kappa*alpha*phi (making lphi = kappa*alpha*phi + kappa*beta*divF)
    unsigned long long int numflopspt = 3;
    Box grid =m_grids[dit[ibox]];
    Bx  grbx = getProtoBox(grid);
    auto& kapfab = m_kappa[dit[ibox]];
#if 1
    ebforallInPlace(numflopspt, "addAlphaPhi", addAlphaPhi, grbx, lphfab, phifab, kapfab, m_alpha, m_beta);
#else
    ebFastforallInPlace(phifab.inputBox(), numflopspt, "addAlphaPhi", addAlphaPhi, grbx, lphfab, phifab, kapfab, m_alpha, m_beta);
#endif
    ideb++;
  }
}
/****/
void
EBMultigrid::
applyOp(EBLevelBoxData<CELL, 1>       & a_lph,
        const EBLevelBoxData<CELL, 1> & a_phi,
        bool a_doExchange) const
{
  return m_finest->applyOp(a_lph, a_phi, a_doExchange);
}
/***/
void
EBMultigrid::
residual(EBLevelBoxData<CELL, 1>       & a_res,
         const EBLevelBoxData<CELL, 1> & a_phi,
         const EBLevelBoxData<CELL, 1> & a_rhs,
         bool a_doExchange) const
{
  PR_TIME("sgmg::resid");
  return m_finest->residual(a_res, a_phi, a_rhs, a_doExchange);
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
                 string                                    a_dombcname[2*DIM],
                 const string                            & a_ebbcname,
                 const Box                               & a_domain,
                 const IntVect                           & a_nghost)
{
  CH_TIME("EBMultigridLevel::define");
  m_depth = 0;

  m_alpha      = a_alpha;      
  m_beta       = a_beta;       
  m_dx         = a_dx;         
  m_grids      = a_grids;      
  m_stenname   = a_stenname;   
  m_neumname   = a_stenname + string("_All_Neumann");   

  for(int ivec = 0; ivec < 2*DIM; ivec++)
  {
    m_dombcname[ivec]  = a_dombcname[ivec];
  }
  m_ebbcname   = a_ebbcname;   
  m_domain     = a_domain;     
  m_nghost     = a_nghost;
  m_dictionary = a_dictionary;

  m_graphs = a_geoserv->getGraphs(m_domain);
  m_resid.define(m_grids, m_nghost, m_graphs);
  m_kappa.define(m_grids, m_nghost, m_graphs);
  m_diagW.define(m_grids, m_nghost, m_graphs);
  
  m_exchangeCopier.exchangeDefine(m_grids, m_nghost);
  //register stencil for apply op
  //true is for need the diagonal wweight
  m_dictionary->registerStencil(m_stenname, m_dombcname, m_ebbcname, m_domain, m_domain, true);

  m_dictionary->registerStencil(m_neumname, StencilNames::Neumann, StencilNames::Neumann, m_domain, m_domain, true);

  //need the volume fraction in a data holder so we can evaluate kappa*alpha I 
  fillKappa(a_geoserv);

  defineCoarserObjects(a_geoserv);
  if(!m_hasCoarser)
  {
    m_bottomSolver = shared_ptr<EBRelaxSolver>(new EBRelaxSolver(this, m_grids, m_graphs, m_nghost));
  }
}
/***/
void
EBMultigridLevel::
defineCoarserObjects(shared_ptr<GeometryService<2> >   & a_geoserv)
{
  PR_TIME("sgmglevel::defineCoarser");

  m_hasCoarser = (m_grids.coarsenable(4));
  if(m_hasCoarser)
  {
    //multilevel operators live with the finer level
    Box coardom = coarsen(m_domain, 2);
    m_nobcname = string("no_bcs");
    m_restrictionName = string("Multigrid_Restriction");
    m_dictionary->registerStencil(m_restrictionName, m_nobcname, m_nobcname, m_domain, coardom, false);
    ///prolongation has ncolors stencils
    for(unsigned int icolor = 0; icolor < s_ncolors; icolor++)
    {
      string colorstring = "PWC_Prolongation_" + std::to_string(icolor);
      m_prolongationName[icolor] = colorstring;
      m_dictionary->registerStencil(colorstring, m_nobcname, m_nobcname, coardom, m_domain, false);
    }
    
    m_coarser = shared_ptr<EBMultigridLevel>(new EBMultigridLevel(*this, a_geoserv));

    const shared_ptr<LevelData<EBGraph>  > graphs = a_geoserv->getGraphs(m_coarser->m_domain);
    m_residC.define(m_coarser->m_grids, m_nghost , graphs);
    m_deltaC.define(m_coarser->m_grids, m_nghost , graphs);
  }
}
/***/
EBMultigridLevel::
EBMultigridLevel(const EBMultigridLevel            & a_finerLevel,
                 shared_ptr<GeometryService<2> >   & a_geoserv)
{
  PR_TIME("sgmglevel::constructor");
  m_depth = a_finerLevel.m_depth + 1;

  m_dx         = 2*a_finerLevel.m_dx;         
  m_domain     = coarsen(a_finerLevel.m_domain, 2);      
  coarsen(m_grids, a_finerLevel.m_grids,  2);      

  m_alpha      = a_finerLevel.m_alpha;      
  m_beta       = a_finerLevel.m_beta;       
  m_stenname   = a_finerLevel.m_stenname;   
  m_neumname   = a_finerLevel.m_neumname;
  for(int ivec = 0; ivec < 2*DIM; ivec++)
  {
    m_dombcname[ivec]  = a_finerLevel.m_dombcname[ivec];
  }
  m_ebbcname   = a_finerLevel.m_ebbcname;   

  m_nghost     = a_finerLevel.m_nghost;
  m_nghost     = a_finerLevel.m_nghost;
  m_dictionary = a_finerLevel.m_dictionary;

  m_graphs = a_geoserv->getGraphs(m_domain);
  m_resid.define(m_grids, m_nghost, m_graphs);
  m_kappa.define(m_grids, m_nghost, m_graphs);
  m_diagW.define(m_grids, m_nghost, m_graphs);

  m_exchangeCopier.exchangeDefine(m_grids, m_nghost);
  //register stencil for apply op
  //true is for need the diagonal wweight
  m_dictionary->registerStencil(m_stenname, m_dombcname, m_ebbcname, m_domain, m_domain, true);

  //should not need the neumann one for coarser levels as TGA only calls it on finest level
  fillKappa(a_geoserv);


  defineCoarserObjects(a_geoserv);
  if(!m_hasCoarser)
  {
    m_bottomSolver = shared_ptr<EBRelaxSolver>(new EBRelaxSolver(this, m_grids, m_graphs, m_nghost));
  }

}
//need the volume fraction in a data holder so we can evaluate kappa*alpha I 
void  
EBMultigridLevel::
fillKappa(const shared_ptr<GeometryService<2> >   & a_geoserv)
{
  CH_TIME("EBMultigridLevel::fillkappa");
  DataIterator dit = m_grids.dataIterator();
  int ideb = 0;
  for(int ibox = 0; ibox < dit.size(); ++ibox)
  {
    auto& kappdat = m_kappa[dit[ibox]];
    auto grid =m_grids[dit[ibox]];
    Bx  grbx = kappdat.inputBox();
    const EBGraph  & graph = (*m_graphs)[dit[ibox]];
    EBHostData<CELL, Real, 1> hostdat(grbx, graph);
    //fill kappa on the host then copy to the device
    a_geoserv->fillKappa(hostdat, grid, dit[ibox], m_domain);
    // now copy to the device
    EBLevelBoxData<CELL, 1>::copyToDevice(hostdat, kappdat);

    shared_ptr<ebstencil_t>                  stencil  = m_dictionary->getEBStencil(m_stenname, m_ebbcname, m_domain, m_domain, ibox);
    shared_ptr< EBBoxData<CELL, Real, 1> >   diagptr  = stencil->getDiagonalWeights();
    auto& diagGhost = m_diagW[dit[ibox]];
    
    Bx gridbx = ProtoCh::getProtoBox(grid);
    size_t numflopspt = 0;
    //this one has to be the old, slow ebforall
    ebforallInPlace(numflopspt, "copyDiag", copyDiag,  gridbx, diagGhost, *diagptr);
    ideb++;
  }
  m_kappa.exchange(m_exchangeCopier);
}
/****/
void
EBMultigridLevel::
residual(EBLevelBoxData<CELL, 1>       & a_res,
         const EBLevelBoxData<CELL, 1> & a_phi,
         const EBLevelBoxData<CELL, 1> & a_rhs,
         bool a_doExchange) const
                    
{
  CH_TIME("EBMultigridLevel::residual");
  //this puts lphi into a_res
  CH_assert(a_res.ghostVect() == a_rhs.ghostVect());
  applyOp(a_res, a_phi, a_doExchange);
  //subtract off rhs so res = lphi - rhs
  DataIterator dit = m_grids.dataIterator();
  int ideb = 0;
  for(int ibox = 0; ibox < dit.size(); ++ibox)
  {
    //this adds alpha*phi (making lphi = alpha*phi + beta*divF)
    unsigned long long int numflopspt = 2;
    auto&       resfab = a_res[dit[ibox]];
    //const auto& phifab = a_phi[dit[ibox]];
    const auto& rhsfab = a_rhs[dit[ibox]];
    Box grid = m_grids[dit[ibox]];
    Bx  grbx = getProtoBox(grid);
#if 1
    ebforallInPlace(numflopspt, "subtractRHS", subtractRHS,  grbx, resfab, rhsfab);
#else
    ebFastforallInPlace(rhsfab.inputBox(), numflopspt, "subtractRHS", subtractRHS,  grbx, resfab, rhsfab);
#endif
    ideb++;
  }
}
/****/
void
EBMultigridLevel::
relax(EBLevelBoxData<CELL, 1>       & a_phi,
      const EBLevelBoxData<CELL, 1> & a_rhs,
      int a_maxiter) const
{
  CH_TIME("EBMultigridLevel::relax");
  CH_assert(a_phi.ghostVect() ==   a_rhs.ghostVect());
  CH_assert(a_phi.ghostVect() == m_kappa.ghostVect());
  CH_assert(a_phi.ghostVect() == m_diagW.ghostVect());
  CH_assert(a_phi.ghostVect() == m_resid.ghostVect());
  //
  ParmParse pp;
  bool do_lazy_relax = false;
  pp.query("do_lazy_relax", do_lazy_relax);
  DataIterator dit = m_grids.dataIterator();
  int ideb = 0;
  for(int iter = 0; iter < a_maxiter; iter++)
  {
    for(int iredblack = 0; iredblack < 2; iredblack++)
    {
      auto & resid = const_cast<EBLevelBoxData<CELL, 1> & >(m_resid);
      bool doExchange = true;
      if(do_lazy_relax)
      {
        doExchange = (iredblack==0);
      }
      residual(resid, a_phi, a_rhs, doExchange);
      {
        CH_TIME("ebforall gsrb bit");
        for(int ibox = 0; ibox < dit.size(); ++ibox)
        {
          shared_ptr<ebstencil_t>                  stencil  = m_dictionary->getEBStencil(m_stenname, m_ebbcname, m_domain, m_domain, ibox);

          Box grid = m_grids[dit[ibox]];
          Bx  grbx = getProtoBox(grid);
          //lambda = safety/diag
          //phi = phi - lambda*(res)
          ///lambda takes floating point to calculate
          //also does an integer check for red/black but I am not sure what to do with that
          auto& phifab   =   a_phi[dit[ibox]];
          auto& resfab   =   resid[dit[ibox]];
          auto& stendiag = m_diagW[dit[ibox]];
//      auto& rhsfab =   a_rhs[dit[ibox]];
          unsigned long long int numflopspt = 10;

#if 1        
          ebforallInPlace_i(numflopspt, "gsrbResid", gsrbResid,  grbx, 
                            phifab, resfab, stendiag,
                            m_kappa[dit[ibox]], m_alpha, m_beta, m_dx, iredblack);
#else           
          ebFastforallInPlace_i(phifab.inputBox(), numflopspt, "gsrbResid", gsrbResid,  grbx, 
                                phifab, resfab, stendiag,
                                m_kappa[dit[ibox]], m_alpha, m_beta, m_dx, iredblack);
#endif
      
          ideb++;
        } //end loop over boxes
        ideb++;
      }
    } //end loop over red and black
    ideb++;
  }// end loop over iteratioons
}
/****/
void
EBMultigridLevel::
restrictResidual(EBLevelBoxData<CELL, 1>       & a_resc,
                 const EBLevelBoxData<CELL, 1> & a_resf)
{
  PR_TIME("sgmglevel::restrict");
  DataIterator dit = m_grids.dataIterator();
  int ideb = 0;
  for(int ibox = 0; ibox < dit.size(); ++ibox)
  {
    auto& coarfab = a_resc[dit[ibox]];
    auto& finefab = a_resf[dit[ibox]];
    //finer level owns the stencil (and the operator)
    Box coardom = coarsen(m_domain, 2);
    shared_ptr<ebstencil_t> stencil =
      m_dictionary->getEBStencil(m_restrictionName, m_nobcname, m_domain, coardom, ibox);
    //set resc = Ave(resf) (true is initToZero)
    stencil->apply(coarfab, finefab,  true, 1.0);
    ideb++;
  }
  ideb++;
}
/****/
void
EBMultigridLevel::
prolongIncrement(EBLevelBoxData<CELL, 1>      & a_phi,
                 const EBLevelBoxData<CELL, 1>& a_cor)
{
  PR_TIME("sgmglevel::prolong");
  //finer level owns the stencil (and the operator)
  Box coardom = coarsen(m_domain, 2);
  DataIterator dit = m_grids.dataIterator();
  for(int icolor = 0; icolor < s_ncolors; icolor++)
  {
    for(int ibox = 0; ibox < dit.size(); ++ibox)
    {
      auto& coarfab = a_cor[dit[ibox]];
      auto& finefab = a_phi[dit[ibox]];
      shared_ptr<ebstencil_t> stencil = m_dictionary->getEBStencil(m_prolongationName[icolor], m_nobcname, coardom, m_domain, ibox);
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
  relax(a_phi,a_rhs, EBMultigrid::s_numSmoothDown); 

  if (m_hasCoarser)
  {
    residual(m_resid,a_phi,a_rhs);                      
    //stencils for multilevel objects live with the finer level
    restrictResidual(m_residC,m_resid);
    m_deltaC.setVal(0.);
    m_coarser->vCycle(m_deltaC,m_residC);
    if(EBMultigrid::s_useWCycle)
    {
      m_coarser->vCycle(m_deltaC,m_residC);
    }
    prolongIncrement(a_phi,m_deltaC);
  }
  else
  {
    typedef BiCGStabSolver<EBLevelBoxData<CELL, 1>, EBMultigridLevel> bicgstab;
    ParmParse pp("bicgstab");
    Real tol   = 1.0e-6;
    Real hang  = 1.0e-8;
    Real small = 1.0e-16;
    int  verb  = 0;
    int  imax  = 0;
    int nrestart = 5;
    pp.query("tol"  , tol);
    pp.query("hang" , hang);
    pp.query("small", small);
    pp.query("imax" , imax);
    pp.query("nrestart", nrestart);
    //the -1.0 is the metric parameter which I do not understand
    //pout() << "calling bicgstab for domain =  " << m_domain << std::endl;
    int status = bicgstab::solve(a_phi, a_rhs, *this, verb, -1.0, tol, hang, small, imax, nrestart);
    if(status != 1)
    {
      pout() << "mild warning: bicgstab returned " << status << std::endl;
    }
  }

  relax(a_phi,a_rhs, EBMultigrid::s_numSmoothUp);

}
#include "Chombo_NamespaceFooter.H"
/****/
