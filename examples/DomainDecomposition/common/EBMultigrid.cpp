#include "EBMultigrid.H"
#include "Proto.H"
#include "base/Proto_Timer.H"
#include "EBRelaxSolver.H"
#include "EBMultigridFunctions.H"
#include <sstream>
#include "BiCGStabSolver.H"
#include "Chombo_ParmParse.H"
#include "Chombo_NamespaceHeader.H"

bool EBMultigridLevel::s_useWCycle     = false;
int  EBMultigridLevel::s_numSmoothDown = 4;
int  EBMultigridLevel::s_numSmoothUp   = 4;

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
EBPoissonOp::
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
    stencil->apply(lphfab, phifab,  true, 1.0);
    //this adds kappa*alpha*phi (making lphi = kappa*alpha*phi + kappa*beta*divF)
    Box grid =m_grids[dit[ibox]];
    Bx  grbx = getProtoBox(grid);
    auto& kapfab = m_kappa[dit[ibox]];

    //pout() << "going into add alphaphi" << endl;
    Bx inputBox = lphfab.inputBox();
    ebforall(inputBox, addAlphaPhi, grbx, lphfab, phifab, kapfab, m_alpha, m_beta);

    ideb++;
  }
}
/****/
void
EBPoissonOp::
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
EBPoissonOp::
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

    stencil->apply(lphfab, phifab,  true, 1.0);
    //this adds kappa*alpha*phi (making lphi = kappa*alpha*phi + kappa*beta*divF)
    Box grid =m_grids[dit[ibox]];
    Bx  grbx = getProtoBox(grid);
    auto& kapfab = m_kappa[dit[ibox]];

    pout() << "going into add alphaphi" << endl;
    Bx inputBox = lphfab.inputBox();
    ebforall(inputBox, addAlphaPhi, grbx, lphfab, phifab, kapfab, m_alpha, m_beta);

    ideb++;
  }
}
/***/
EBPoissonOp::
EBPoissonOp(dictionary_t                            & a_dictionary,
            shared_ptr<GeometryService<2> >         & a_geoserv,
            const Real                              & a_alpha,
            const Real                              & a_beta,
            const Real                              & a_dx,
            const DisjointBoxLayout                 & a_grids,
            const string                            & a_stenname,
            string                                    a_dombcname[2*DIM],
            const string                            & a_ebbcname,
            const Box                               & a_domain,
            const IntVect                           & a_nghost,
            string a_prefix):  EBMultigridLevel()
{
  CH_TIME("EBPoissonOp::define");
  m_prefix = a_prefix;
  getDToB();
  m_depth = 0;
  m_geoserv = a_geoserv;
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
  Point pghost = ProtoCh::getPoint(m_nghost);
  m_dictionary->registerStencil(m_stenname, m_dombcname, m_ebbcname, m_domain, m_domain, true, pghost);

  m_dictionary->registerStencil(m_neumname, StencilNames::Neumann, StencilNames::Neumann, m_domain, m_domain, true, pghost);

  //need the volume fraction in a data holder so we can evaluate kappa*alpha I 
  fillKappa(a_geoserv);

  if(!m_directToBottom)
  {
    defineCoarserObjects(a_geoserv);
  }
  if(!m_hasCoarser || m_directToBottom)
  {
    defineBottomSolvers(a_geoserv);
  }
}
/***/
void
EBPoissonOp::
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
    
    m_coarser = shared_ptr<EBPoissonOp>(new EBPoissonOp());
    m_coarser->define(*this, a_geoserv);

    const shared_ptr<LevelData<EBGraph>  > graphs = a_geoserv->getGraphs(m_coarser->m_domain);
    m_residC.define(m_coarser->m_grids, m_nghost , graphs);
    m_deltaC.define(m_coarser->m_grids, m_nghost , graphs);
  }
}
/***/
void
EBPoissonOp::
define(const EBMultigridLevel            & a_finerLevel,
       shared_ptr<GeometryService<2> >   & a_geoserv)
{
  PR_TIME("sgmglevel::constructor");
  m_depth = a_finerLevel.m_depth + 1;
  m_prefix = a_finerLevel.m_prefix;
  getDToB();
  m_geoserv = a_geoserv;
  m_dx         = 2*a_finerLevel.m_dx;         
  m_domain     = coarsen(a_finerLevel.m_domain, 2);      
  coarsen(m_grids, a_finerLevel.m_grids,  2);      

  m_alpha      = a_finerLevel.m_alpha;      
  m_beta       = a_finerLevel.m_beta;
  EBPoissonOp* finer = (EBPoissonOp*)(& a_finerLevel);
  m_stenname   = finer->m_stenname;   
  m_neumname   = finer->m_neumname;
  for(int ivec = 0; ivec < 2*DIM; ivec++)
  {
    m_dombcname[ivec]  = finer->m_dombcname[ivec];
  }
  m_ebbcname   = finer->m_ebbcname;   
  m_nghost     = finer->m_nghost;
  m_nghost     = finer->m_nghost;
  m_dictionary = finer->m_dictionary;

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

  if(!m_directToBottom)
  {
    defineCoarserObjects(a_geoserv);
  }
  if(!m_hasCoarser || m_directToBottom)
  {
    defineBottomSolvers(a_geoserv);
  }
}

void  
EBPoissonOp::
defineBottomSolvers(shared_ptr<GeometryService<2> >   & a_geoserv)
{
  m_relaxSolver = shared_ptr<EBRelaxSolver>(new EBRelaxSolver(this, m_grids, m_graphs, m_nghost));
#ifdef CH_USE_PETSC

  ParmParse pp(m_prefix.c_str());
  string which_solver("relax");
  pp.query("bottom_solver", which_solver);

  if(which_solver == string("petsc"))
  {
    Point pghost= ProtoCh::getPoint(m_nghost);
    EBPetscSolver<2>* ptrd = 
      (new EBPetscSolver<2>(a_geoserv, m_dictionary, m_graphs, m_grids, m_domain,
                            m_stenname, m_dombcname, m_ebbcname,
                            m_dx, m_alpha, m_beta, pghost));
    m_petscSolver = shared_ptr<EBPetscSolver<2> >(ptrd);
  }
#endif    
}
//need the volume fraction in a data holder so we can evaluate kappa*alpha I 
void  
EBPoissonOp::
fillKappa(const shared_ptr<GeometryService<2> >   & a_geoserv)
{
  CH_TIME("EBPoissonOp::fillkappa");
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
    EBLevelBoxData<CELL, 1>::copyToDevice(kappdat, hostdat);

    shared_ptr<ebstencil_t>                  stencil  = m_dictionary->getEBStencil(m_stenname, m_ebbcname, m_domain, m_domain, ibox);
    shared_ptr< EBBoxData<CELL, Real, 1> >   diagptr  = stencil->getDiagonalWeights();
    auto& diagGhost = m_diagW[dit[ibox]];
    auto& stendiag = *diagptr;
    Bx inputBox = diagGhost.inputBox();

    ebforall(inputBox, copyDiag,  inputBox, diagGhost, stendiag);
    ideb++;
  }
  m_kappa.exchange(m_exchangeCopier);
}
/****/
void
EBPoissonOp::
residual(EBLevelBoxData<CELL, 1>       & a_res,
         const EBLevelBoxData<CELL, 1> & a_phi,
         const EBLevelBoxData<CELL, 1> & a_rhs,
         bool a_doExchange) const
                    
{
  CH_TIME("EBPoissonOp::residual");
  //this puts lphi into a_res
  CH_assert(a_res.ghostVect() == a_rhs.ghostVect());
  applyOp(a_res, a_phi, a_doExchange);
  //subtract off rhs so res = lphi - rhs
  DataIterator dit = m_grids.dataIterator();
  for(int ibox = 0; ibox < dit.size(); ++ibox)
  {
    //this adds alpha*phi (making lphi = alpha*phi + beta*divF)

    auto&       resfab = a_res[dit[ibox]];
    //const auto& phifab = a_phi[dit[ibox]];
    const auto& rhsfab = a_rhs[dit[ibox]];
    Box grid = m_grids[dit[ibox]];
    Bx  grbx = getProtoBox(grid);
    Bx inputBox = resfab.inputBox();
    ebforall(inputBox, subtractRHS,  grbx, resfab, rhsfab);
  }

  a_res.exchange(m_exchangeCopier);
}
/****/
void
EBPoissonOp::
relax(EBLevelBoxData<CELL, 1>       & a_phi,
      const EBLevelBoxData<CELL, 1> & a_rhs,
      int a_maxiter) const
{
  CH_TIME("EBPoissonOp::relax");
  CH_assert(a_phi.ghostVect() ==   a_rhs.ghostVect());
  CH_assert(a_phi.ghostVect() == m_kappa.ghostVect());
  CH_assert(a_phi.ghostVect() == m_diagW.ghostVect());
  CH_assert(a_phi.ghostVect() == m_resid.ghostVect());
  //
  ParmParse pp(m_prefix.c_str());
  bool do_lazy_relax = false;
  bool one_exchange_per_relax = false;
  pp.query("do_lazy_relax", do_lazy_relax);
  pp.query("one_exchange_per_relax", one_exchange_per_relax);
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
        if(doExchange && one_exchange_per_relax)
        {
          doExchange = (iter==0);
        }
      }

      residual(resid, a_phi, a_rhs, doExchange);

      for(int ibox = 0; ibox < dit.size(); ++ibox)
      {

        Box grid = m_grids[dit[ibox]];
        Bx  grbx = getProtoBox(grid);
        //lambda = safety/diag
        //phi = phi - lambda*(res)
        ///lambda takes floating point to calculate
        //also does an integer check for red/black but I am not sure what to do with that
        auto& phifab   =   a_phi[dit[ibox]];
        auto& resfab   =   resid[dit[ibox]];
        auto& stendiag = m_diagW[dit[ibox]];
        auto& kappa    = m_kappa[dit[ibox]];
        Bx  inputBox = phifab.inputBox();
        ebforall_i(inputBox, gsrbResid,  grbx, 
                   phifab, resfab, stendiag,
                   kappa, m_alpha, m_beta, m_dx, iredblack);
      
        ideb++;
      } //end loop over boxes
      //begin debug
      //      a_phi.exchange(m_exchangeCopier);
      //      ideb++;
      //end debug
    } //end loop over red and black
    ideb++;
  }// end loop over iteratioons
}
/****/
void
EBPoissonOp::
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
EBPoissonOp::
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
///
void
EBPoissonOp::
bottom_solve(EBLevelBoxData<CELL, 1>         & a_phi,
             const EBLevelBoxData<CELL, 1>   & a_rhs)
{

  string which_solver("relax");
  ParmParse pp(m_prefix.c_str());

  pp.query("bottom_solver", which_solver);
  if(which_solver == string("bicgstab"))
  {
    solve_bicgstab(a_phi, a_rhs);
  }
  else if(which_solver == string("relax"))
  {
    solve_relax(a_phi, a_rhs);
  }
#ifdef CH_USE_PETSC    
  else if(which_solver == string("petsc"))
  {
    solve_petsc(a_phi, a_rhs);
  }
#endif
  else
  {
    MayDay::Error("bottom solver identifier not found");
  }
}
///
void
EBPoissonOp::
solve_bicgstab(EBLevelBoxData<CELL, 1>         & a_phi,
               const EBLevelBoxData<CELL, 1>   & a_rhs)
{
  ParmParse pp("bicgstab");
  typedef BiCGStabSolver<EBLevelBoxData<CELL, 1>, EBPoissonOp> bicgstab;
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
  pp.query("verbosity", verb);
  //the -1.0 is the metric parameter which I do not understand
  //pout() << "calling bicgstab for domain =  " << m_domain << std::endl;
  int status = bicgstab::solve(a_phi, a_rhs, *this, verb, -1.0, tol, hang, small, imax, nrestart);
  if(status != 1)
  {
    pout() << "mild warning: bicgstab returned " << status << std::endl;
  }
}
///
void
EBPoissonOp::
solve_relax(EBLevelBoxData<CELL, 1>         & a_phi,
            const EBLevelBoxData<CELL, 1>   & a_rhs)
{
  ParmParse pp("relax_solver");
  //typedef BiCGStabSolver<EBLevelBoxData<CELL, 1>, EBPoissonOp> bicgstab;
  Real tol   = 1.0e-6;
  int  imax  = 0;

  pp.query("tol"  , tol);
  pp.query("imax" , imax);

  m_relaxSolver->solve(a_phi, a_rhs, imax, tol);
}
///

#ifdef CH_USE_PETSC
void
EBPoissonOp::
solve_petsc(EBLevelBoxData<CELL, 1>         & a_phi,
            const EBLevelBoxData<CELL, 1>   & a_rhs)
{
  m_petscSolver->solve(a_phi, a_rhs);
}
#endif
#include "Chombo_NamespaceFooter.H"
/****/
