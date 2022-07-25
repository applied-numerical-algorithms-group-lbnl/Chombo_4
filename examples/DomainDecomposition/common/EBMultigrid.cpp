#include "EBMultigrid.H"
#include "Proto.H"
#include "base/Proto_Timer.H"
#include "EBRelaxSolver.H"
#include "EBMultigridFunctions.H"
#include <sstream>
#include "BiCGStabSolver.H"
#include "Chombo_ParmParse.H"
#include "Chombo_NamespaceHeader.H"
/******/
void
EBPoissonOp::
createAndStoreRestrictionStencil()
{
  DataIterator dit = m_grids.dataIterator();
  vector<EBIndex<CELL> >                    dstVoFs;                    
  vector<LocalStencil<CELL, double> >       stenVec;                    
  Point ghostPt = ProtoCh::getPoint(m_nghost);

  auto domCoar = m_domain;
  auto domFine = refine(m_domain,2);
  
  auto dblFine = m_geoserv->getDBL(domFine);
  auto dblCoar = m_geoserv->getDBL(domCoar);
  
  shared_ptr<graph_distrib_t> graphCoar = m_graphs;
  shared_ptr<graph_distrib_t> graphFine;
  if(dblFine.compatible(dblCoar))
  {
    graphFine = m_geoserv->getGraphs(domFine);
  }
  else
  {
    graphFine = getGraphOnRefinedCoarseLayout();
  }
  m_restrictionStencil.resize(dit.size());
  for(int ibox = 0; ibox < dit.size(); ++ibox)
  {

    auto gridCoar = m_grids[dit[ibox]];
    auto gridFine = refine(gridCoar, 2);
    
    Bx  dstValid  = ProtoCh::getProtoBox(gridCoar);
    Bx  srcValid  = ProtoCh::getProtoBox(gridFine);
    Bx  dstDomain = ProtoCh::getProtoBox(domCoar);
    Bx  srcDomain = ProtoCh::getProtoBox(domFine);
    
    const auto& srcGraph = (*graphFine)[dit[ibox]];
    const auto& dstGraph = (*graphCoar)[dit[ibox]];

    getRestrictionStencil(dstVoFs, stenVec,  gridCoar, dstGraph);
    
    m_restrictionStencil[ibox] = shared_ptr<ebstencil_t>
      (new ebstencil_t(dstVoFs,   stenVec,
                       srcGraph,  dstGraph, 
                       srcValid,  dstValid, 
                       srcDomain, dstDomain,
                       ghostPt, ghostPt, false));
  }
}
/******/
void
EBPoissonOp::
createAndStoreProlongationStencils()
{
  auto domCoar = m_domain;
  auto domFine = refine(m_domain,2);
  
  auto dblFine = m_geoserv->getDBL(domFine);
  auto dblCoar = m_geoserv->getDBL(domCoar);
  
  shared_ptr<graph_distrib_t> graphCoar = m_graphs;
  shared_ptr<graph_distrib_t> graphFine;
  if(dblFine.compatible(dblCoar))
  {
    graphFine = m_geoserv->getGraphs(domFine);
  }
  else
  {
    graphFine = getGraphOnRefinedCoarseLayout();
  }
  
  Point ghostPt = ProtoCh::getPoint(m_nghost);
  ///prolongation has ncolors stencils
  for(unsigned int icolor = 0; icolor < s_ncolors; icolor++)
  {
    DataIterator dit = m_grids.dataIterator();
    m_prolongationStencils[icolor].resize(dit.size());
    for(int ibox = 0; ibox < dit.size(); ++ibox)
    {
      auto gridCoar = m_grids[dit[ibox]];
      auto gridFine = refine(gridCoar, 2);
      const EBGraph& srcGraph = (*graphCoar)[dit[ibox]];
      const EBGraph& dstGraph = (*graphFine)[dit[ibox]];
      
      vector<EBIndex<CELL> >                    dstVoFs;                    
      vector<LocalStencil<CELL, double> >       stenVec;                    
      auto srcValid = ProtoCh::getProtoBox(gridCoar);
      auto dstValid = ProtoCh::getProtoBox(gridFine);
      auto srcDomain= ProtoCh::getProtoBox(domCoar);
      auto dstDomain= ProtoCh::getProtoBox(domFine);
      
      getProlongationStencil(dstVoFs, stenVec, dstValid, dstGraph, icolor);

      m_prolongationStencils[icolor][ibox] = shared_ptr<ebstencil_t>
        (new ebstencil_t(dstVoFs, stenVec,
                         srcGraph,  dstGraph, 
                         srcValid,  dstValid, 
                         srcDomain, dstDomain,
                         ghostPt, ghostPt, false));
    }
  }
}
/******/
void
EBPoissonOp::
getRestrictionStencil(vector<EBIndex<CELL> >                    & a_dstVoFs,                    
                      vector<LocalStencil<CELL, double> >       & a_stencil,                    
                      const Box                                 & a_validCoar,
                      const EBGraph                             & a_graphCoar)
{
  a_dstVoFs.resize(0);
  a_stencil.resize(0);
  //this stencil has no span that I can understand so getIrrregLocations not appropriate
  a_dstVoFs = a_graphCoar.getAllVoFs(a_validCoar); 

  double dnumpts = double(s_ncolors);
  double rweight = 1.0/dnumpts;
  a_stencil.resize(a_dstVoFs.size());
  for(int ivof = 0; ivof < a_dstVoFs.size(); ivof++)
  {
    const EBIndex<CELL>&   coarVoF  = a_dstVoFs[ivof];
    vector<EBIndex<CELL> > fineVoFs = a_graphCoar.refine(coarVoF);
    //in multigrid, the rhs is already volume weighted so it is just
    //coarse = (1/ncolors)*(sum(fine))
    for(int ifine = 0; ifine <fineVoFs.size(); ifine++)
    {
      a_stencil[ivof].add(fineVoFs[ifine], rweight);
    }
  }
}

/******/
void
EBPoissonOp::
getProlongationStencil(vector<EBIndex<CELL> >                    & a_dstVoFs,                    
                       vector<LocalStencil<CELL, double> >       & a_stencil,                    
                       const Box                                 & a_dstValid,                   
                       const EBGraph                             & a_dstGraph,                      
                       unsigned long                               a_icolor)               
{
  a_dstVoFs.resize(0);
  a_stencil.resize(0);
  Point colorpt = Proto::EBStencilArchive<CELL, CELL, 2, double>::getColor(a_icolor);

  a_dstVoFs = a_dstGraph.getAllVoFs(a_dstValid); 
  vector<EBIndex<CELL> > coarsePts = a_dstGraph.getAllVoFs(a_dstValid);

  a_stencil.resize(a_dstVoFs.size());
  for(int ivof = 0; ivof < a_dstVoFs.size(); ivof++)
  {
    const EBIndex<CELL>&   fineVoF = a_dstVoFs[ivof];
    EBIndex<CELL>          coarVoF = a_dstGraph.coarsen(fineVoF);
    double weight = 1;

    Point vofpt = fineVoF.m_pt;
    for(int idir = 0; idir < DIM; idir++)
    {
      int mod = (vofpt[idir])%2;
      if(mod != colorpt[idir])
      {
        weight  = 0;
      }
    }
      
    a_stencil[ivof].add(coarVoF, weight);
  }

}


/****/
void 
EBMultigrid::
solve(EBLevelBoxData<CELL, 1>       & a_phi,
      const EBLevelBoxData<CELL, 1> & a_rhs,
      const Real                    & a_tol,
      const unsigned int            & a_maxIter,
      bool a_initToZero,
      bool a_printStuff)
{
  CH_TIME("EBMultigrid::solve");
  if(a_initToZero)
  {
    a_phi.setVal(0.);
  }

  bool printStuff = false;
  bool doExchange = true;
  residual(m_res, a_phi, a_rhs, doExchange, printStuff);
  EBIndex<CELL> vofmax;
  Real initres = m_res.maxNorm(vofmax, 0);
  int iter = 0;
  pout() << "EBMultigrid::solve tol = " << a_tol << ",  max iter = "<< a_maxIter << endl;
  Real resnorm = initres;
  Real resnormold = resnorm;
  while((iter < a_maxIter) && (resnorm > a_tol*initres))
  {
    m_cor.setVal(0.);
    pout() << setprecision(3)
           << setiosflags(ios::showpoint)
           << setiosflags(ios::scientific);

    
    pout() << "EBMultigrid::solve iter = " << iter << ", |resid| = " << resnorm << "@ " << vofmax.m_pt;
    Real rate = 1;
    if((resnormold > 1.0e-12) && (iter > 0))
    {
      rate = resnormold/resnorm;
      pout() << ", rate = " << rate;
    }
    pout() << endl;
    
    {//begin debug
      pout() << "before vcycle m_cor:" << endl;
      DumpArea::dumpLDCell1Area(&m_cor);
      pout() << "before vcycle m_res:" << endl;
      DumpArea::dumpLDCell1Area(&m_res);
    }//end debug
    
    vCycle(m_cor, m_res, a_printStuff);
    
    {//begin debug
      pout() << "relax m_cor:" << endl;
      DumpArea::dumpLDCell1Area(&m_cor);
      pout() << "relax m_res:" << endl;
      DumpArea::dumpLDCell1Area(&m_res);
        
    }//end debug

    a_phi += m_cor;
    residual(m_res, a_phi, a_rhs, doExchange, printStuff);
    
    resnormold = resnorm;
    resnorm = m_res.maxNorm(vofmax, 0);

    iter++;
  }
  pout() << "EBMultigrid:solve: final     |resid| = " << resnorm << "@ " <<vofmax.m_pt << endl;
}

/****/
void
EBPoissonOp::
applyOp(EBLevelBoxData<CELL, 1>       & a_lph,
        const EBLevelBoxData<CELL, 1> & a_phi,
        bool a_doExchange,
        bool a_printStuff) const
                    
{
  CH_TIME("EBMultigrid::applyOp");
  CH_assert(  a_lph.ghostVect() == a_phi.ghostVect());
  CH_assert(m_kappa.ghostVect() == a_phi.ghostVect());
  EBLevelBoxData<CELL, 1>& phi = const_cast<EBLevelBoxData<CELL, 1>&>(a_phi);
  if(a_doExchange)
  {
    phi.exchange();
  }

  DataIterator dit = m_grids.dataIterator();
  int ideb  = 0;
  for(int ibox = 0; ibox < dit.size(); ++ibox)
  {
    auto valid = m_grids[dit[ibox]];
    //begin debug
    IntVect ivdeb(D_DECL(6, 16, 0));
    IntVect ivlod(D_DECL(3, 11, 0));
    IntVect ivhid(D_DECL(9, 21, 0));
    Box debval(ivlod, ivhid);
    Box grownValid = grow(valid, m_nghost);
    debval &= grownValid;
    Point   ptlo = ProtoCh::getPoint(debval.smallEnd());
    Point   pthi = ProtoCh::getPoint(debval.bigEnd());
    bool printStuff = (a_printStuff && (valid.contains(ivdeb)));
    //end debug
    
    auto& phifab = a_phi[dit[ibox]];
    auto& lphfab = a_lph[dit[ibox]];
    shared_ptr<ebstencil_t> stencil =
      m_dictionary->getEBStencil(m_stenname, m_ebbcname, m_domain, m_domain, ibox);
    //set lphi = kappa* div(F)

    if(printStuff) 
    {//begin debug
      pout() << "phifab: " << endl;
      DumpArea::genDumpCell1(&phifab, ptlo, pthi);
      pout() << "lphfab: " << endl;
      DumpArea::genDumpCell1(&phifab, ptlo, pthi);
    }//end debug

    stencil->apply(lphfab, phifab,  true, 1.0);

    if(printStuff) 
    {//begin debug
      pout() << "lphfab: " << endl;
      DumpArea::genDumpCell1(&phifab, ptlo, pthi);
    }//end debug
    
    //this adds kappa*alpha*phi (making lphi = kappa*alpha*phi + kappa*beta*divF)
    Box grid =m_grids[dit[ibox]];
    Bx  grbx = getProtoBox(grid);
    auto& kapfab = m_kappa[dit[ibox]];

    //pout() << "going into add alphaphi" << endl;
    Bx inputBox = lphfab.inputBox();
    ebforall(inputBox, addAlphaPhi, grbx, lphfab, phifab, kapfab, m_alpha, m_beta);

    if(printStuff) 
    {//begin debug
      pout() << "kapfab: " << endl;
      DumpArea::genDumpCell1(&kapfab, ptlo, pthi);
      pout() << "lphfab: " << endl;
      DumpArea::genDumpCell1(&lphfab, ptlo, pthi);
    }//end debug
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
  phi.exchange();
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
            string   a_bottom_solver,
            bool     a_direct_to_bottom,
            string a_prefix,
            bool a_useWCycle,
            int  a_numSmooth,
            bool a_printStuff):EBMultigridLevel()
{
  CH_TIME("EBPoissonOp::define");
  m_prefix = a_prefix;
  m_bottom_solver    = a_bottom_solver    ;
  m_direct_to_bottom = a_direct_to_bottom ;
  m_prefix           = a_prefix           ;
  m_useWCycle        = a_useWCycle        ;
  m_numSmooth        = a_numSmooth        ;    

  m_depth = 0;
  m_hasRecoDataAlready  = false;
  m_hasRecoGraphAlready = false;
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
  
  auto geodbl = m_geoserv->getDBL(m_domain);
  PR_assert(geodbl==m_grids);
  
  if(a_printStuff)
  {
    pout() << "EBPoissonOp constructor: defining internal data" << endl;
  }
  m_graphs = a_geoserv->getGraphs(m_domain);
  m_resid.define(m_grids, m_nghost, m_graphs);
  m_kappa.define(m_grids, m_nghost, m_graphs);
  m_diagW.define(m_grids, m_nghost, m_graphs);
  
  if(a_printStuff)
  {
    pout() << "EBPoissonOp constructor: registering stencil" << endl;
  }
  //register stencil for apply op
  //true is for need the diagonal wweight
  Point pghost = ProtoCh::getPoint(m_nghost);
  m_dictionary->registerStencil(m_stenname, m_dombcname, m_ebbcname, m_domain, m_domain, true, pghost);

  m_dictionary->registerStencil(m_neumname, StencilNames::Neumann, StencilNames::Neumann, m_domain, m_domain, true, pghost);

  if(a_printStuff)
  {
    pout() << "EBPoissonOp constructor: filling volume fraction holder" << endl;
  }
  //need the volume fraction in a data holder so we can evaluate kappa*alpha I 
  fillKappa(a_geoserv, a_printStuff);

  if(!m_direct_to_bottom)
  {
    if(a_printStuff)
    {
      pout() << "EBPoissonOp constructor: making coarser objects" << endl;
    }
    defineCoarserObjects(a_geoserv, a_printStuff);
  }
  if(!m_hasCoarser || m_direct_to_bottom)
  {
    if(a_printStuff)
    {
      pout() << "EBPoissonOp constructor: defining bottom solvers" << endl;
    }
    defineBottomSolvers(a_geoserv, a_printStuff);
  }
  if(a_printStuff)
  {
    pout() << "EBPoissonOp constructor: waiting for other procs to catch up" << endl;
  }
  CH4_SPMD::barrier();
  if(a_printStuff)
  {
    pout() << "EBPoissonOp constructor: leaving" << endl;
  }
}
/***/
void
EBPoissonOp::
defineCoarserObjects(shared_ptr<GeometryService<2> >   & a_geoserv,
                     bool a_printStuff)
{
  PR_TIME("EBPoissonOp::defineCoarser");

  Box coardom = coarsen(m_domain, 2);
  m_hasCoarser = a_geoserv->hasLevel(coardom);
  
  if(m_hasCoarser)
  {
    m_coarser = shared_ptr<EBPoissonOp>(new EBPoissonOp());
    m_coarser->define(*this, a_geoserv, a_printStuff);
    
    auto graphs = a_geoserv->getGraphs(m_coarser->m_domain);
    m_residC.define(m_coarser->m_grids, m_nghost , graphs);
    m_deltaC.define(m_coarser->m_grids, m_nghost , graphs);
  }
}
/***/
void
EBPoissonOp::
define(const EBMultigridLevel            & a_finerLevel,
       shared_ptr<GeometryService<2> >   & a_geoserv,
       bool a_printStuff)
{
  PR_TIME("EBPoissonOp::constructor");
  m_depth            = a_finerLevel.m_depth + 1;

  m_numSmooth        = a_finerLevel.m_numSmooth;
  m_useWCycle        = a_finerLevel.m_useWCycle;
  m_prefix           = a_finerLevel.m_prefix;
  m_bottom_solver    = a_finerLevel.m_bottom_solver;
  m_direct_to_bottom = a_finerLevel.m_direct_to_bottom;
  
  m_hasRecoDataAlready  = false;
  m_hasRecoGraphAlready = false;
  m_geoserv = a_geoserv;
  m_dx         = 2*a_finerLevel.m_dx;         
  m_domain     = coarsen(a_finerLevel.m_domain, 2);      
  m_grids = m_geoserv->getDBL(m_domain);


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

  //register stencil for apply op
  //true is for need the diagonal wweight
  m_dictionary->registerStencil(m_stenname, m_dombcname, m_ebbcname, m_domain, m_domain, true);

  createAndStoreRestrictionStencil();
  createAndStoreProlongationStencils();
  //should not need the neumann one for coarser levels as TGA only calls it on finest level
  fillKappa(a_geoserv, false);

  if(!m_direct_to_bottom)
  {
    defineCoarserObjects(a_geoserv, a_printStuff);
  }
  if(!m_hasCoarser || m_direct_to_bottom)
  {
    defineBottomSolvers(a_geoserv, a_printStuff);
  }
}

void  
EBPoissonOp::
defineBottomSolvers(shared_ptr<GeometryService<2> >   & a_geoserv,
                    bool a_printStuff)
{
  if(a_printStuff)
  {
    pout() << "EBPoissonOp::defineBottomSolver: defining relax solver" << endl;
  }
  m_relaxSolver = shared_ptr<EBRelaxSolver>(new EBRelaxSolver(this, m_grids, m_graphs, m_nghost));
#ifdef CH_USE_PETSC

  ParmParse pp(m_prefix.c_str());
  string which_solver("relax");
  pp.query("bottom_solver", which_solver);

  if(a_printStuff)
  {
    pout() << "EBPoissonOp::defineBottomSolver: defining petsc solver" << endl;
  }
  if(which_solver == string("petsc"))
  {
    Point pghost= ProtoCh::getPoint(m_nghost);
    EBPetscSolver<2>* ptrd = 
      (new EBPetscSolver<2>(a_geoserv, m_dictionary, m_graphs, m_grids, m_domain,
                            m_stenname, m_dombcname, m_ebbcname,
                            m_dx, m_alpha, m_beta, pghost, a_printStuff));
    m_petscSolver = shared_ptr<EBPetscSolver<2> >(ptrd);
  }
#endif    
}
//need the volume fraction in a data holder so we can evaluate kappa*alpha I 
void  
EBPoissonOp::
fillKappa(const shared_ptr<GeometryService<2> >   & a_geoserv,
          bool a_printStuff)
{
  CH_TIME("EBPoissonOp::fillkappa");
  if(a_printStuff)
  {
    pout() << "EBPoissonOp::fillKappa looping through boxes to fill valid data" << endl;
  }
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
  if(a_printStuff)
  {
    pout() << "EBPoissonOp::fillKappa: calling exchange to fill ghost data" << endl;
  }
  m_kappa.exchange();
  if(a_printStuff)
  {
    pout() << "EBPoissonOp::fillKappa: waiting for other procs to catch up" << endl;
  }
  CH4_SPMD::barrier();
}
/****/
void
EBPoissonOp::
residual(EBLevelBoxData<CELL, 1>       & a_res,
         const EBLevelBoxData<CELL, 1> & a_phi,
         const EBLevelBoxData<CELL, 1> & a_rhs,
         bool a_doExchange,
         bool a_printStuff) const
                    
{
  CH_TIME("EBPoissonOp::residual");
  //this puts lphi into a_res
  CH_assert(a_res.ghostVect() == a_rhs.ghostVect());
  //this appears to be necessary.
  a_res.setVal(0.);
  if(a_printStuff)
  {//begin debug
    pout() << "residual a_phi:" << endl;
    DumpArea::dumpAsOneBox(&a_phi, m_geoserv);
    pout() << "residual a_res:" << endl;
    DumpArea::dumpAsOneBox(&a_res, m_geoserv);
    pout() << "residual a_rhs:" << endl;
    DumpArea::dumpAsOneBox(&a_rhs, m_geoserv);
        
  }//end debug
  
  applyOp(a_res, a_phi, a_doExchange, a_printStuff);

  if(a_printStuff)
  {//begin debug
    pout() << "residual a_res:" << endl;
    DumpArea::dumpAsOneBox(&a_res, m_geoserv);
  }//end debug


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

  if(a_printStuff)
  {//begin debug
    pout() << "residual a_res:" << endl;
    DumpArea::dumpAsOneBox(&a_res, m_geoserv);
  }//end debug

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
      resid.setVal(0.);
      bool doExchange = true;
      bool printStuff = false;
      if(do_lazy_relax)
      {
        doExchange = (iredblack==0);
        if(doExchange && one_exchange_per_relax)
        {
          doExchange = (iter==0);
        }
      }

      residual(resid, a_phi, a_rhs, doExchange, printStuff);

//      {//begin debug
//        pout() << "phi: " << endl;
//        DumpArea::dumpAsOneBox(&a_phi,m_geoserv);
//        pout() << "resid: " << endl;
//        DumpArea::dumpAsOneBox(&resid, m_geoserv);
//      }//end debug
      
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
//      {//begin debug
//        pout() << "phi: " << endl;
//        DumpArea::dumpAsOneBox(&a_phi,m_geoserv);
//      }//end debug
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
  PR_TIME("EBPoissonOp::restrict");
  auto dblFine = a_resf.disjointBoxLayout();
  auto dblCoar = a_resc.disjointBoxLayout();
  if(dblFine.compatible(dblCoar))
  {
    restrictResidualOnProc(a_resc, a_resf);
  }
  else
  {
    restrictResidualAgglom(a_resc, a_resf);
  }
}
/****/
void
EBPoissonOp::
restrictResidualAgglom(EBLevelBoxData<CELL, 1>       & a_resc,
                       const EBLevelBoxData<CELL, 1> & a_resf)
{
  PR_TIME("EBPoissonOp::restrictAgglom");
  shared_ptr<EBLevelBoxData<CELL, 1> > resReCo = getDataOnRefinedCoarseLayout();
  a_resf.copyTo(*resReCo);
  restrictResidualOnProc(a_resc, *resReCo);
}
/****/
void
EBPoissonOp::
restrictResidualOnProc(EBLevelBoxData<CELL, 1>       & a_resc,
                       const EBLevelBoxData<CELL, 1> & a_resf)
{
  PR_TIME("EBPoissonOp::restrictOnProc");
  DataIterator dit = m_grids.dataIterator();
  for(int ibox = 0; ibox < dit.size(); ++ibox)
  {
    auto& coarfab = a_resc[dit[ibox]];
    auto& finefab = a_resf[dit[ibox]];
    //finer level owns the stencil (and the operator)
    Box coardom = coarsen(m_domain, 2);
    shared_ptr<ebstencil_t> stencil = m_restrictionStencil[ibox];

    //set resc = Ave(resf) (true is initToZero)
    stencil->apply(coarfab, finefab,  true, 1.0);
  }
}
/****/
void
EBPoissonOp::
prolongIncrement(EBLevelBoxData<CELL, 1>      & a_phi,
                 const EBLevelBoxData<CELL, 1>& a_cor)
{
  PR_TIME("EBPoissonOp::prolong");

  auto dblFine = a_phi.disjointBoxLayout();
  auto dblCoar = a_cor.disjointBoxLayout();
  if(dblFine.compatible(dblCoar))
  {
    prolongIncrementOnProc(a_phi, a_cor);
  }
  else
  {
    prolongIncrementAgglom(a_phi, a_cor);
  }
}
/****/
void
EBPoissonOp::
prolongIncrementOnProc(EBLevelBoxData<CELL, 1>      & a_phi,
                       const EBLevelBoxData<CELL, 1>& a_cor)
{
  PR_TIME("EBPoissonOp::prolong");
  //finer level owns the stencil (and the operator)
  Box coardom = coarsen(m_domain, 2);
  DataIterator dit = m_grids.dataIterator();
  for(int icolor = 0; icolor < s_ncolors; icolor++)
  {
    for(int ibox = 0; ibox < dit.size(); ++ibox)
    {
      auto& coarfab = a_cor[dit[ibox]];
      auto& finefab = a_phi[dit[ibox]];
      shared_ptr<ebstencil_t> stencil = m_prolongationStencils[icolor][ibox];
      //phi  = phi + I(correction) (false means do not init to zero)
      stencil->apply(finefab, coarfab,  false, 1.0);
    }
  }
}
/****/
void
EBPoissonOp::
prolongIncrementAgglom(EBLevelBoxData<CELL, 1>      & a_phi,
                       const EBLevelBoxData<CELL, 1>& a_cor)
{
  PR_TIME("EBPoissonOp::prolongAgglom");
  //finer level owns the stencil (and the operator)
  
  shared_ptr<EBLevelBoxData<CELL, 1> >phiReCo = getDataOnRefinedCoarseLayout();
  a_phi.copyTo(*phiReCo);
  prolongIncrementOnProc(*phiReCo, a_cor);
  phiReCo->copyTo(a_phi);
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
