#include <cmath>
#include <cstdio>
#include <iostream>


#include "EBProto.H"
#include "Chombo_EBLevelBoxData.H"
#include "Chombo_LevelData.H"
#include "Chombo_BaseFab.H"

#include "Chombo_ParmParse.H"
#include "Chombo_LoadBalance.H"
#include "Chombo_ProtoInterface.H"
#include "Chombo_BRMeshRefine.H"
#include "Chombo_GeometryService.H"
#include "Chombo_EBDictionary.H"
#include "Chombo_EBChombo.H"
#include "EBMultigrid.H"
#include "Proto_DebugHooks.H"
#include "DebugFunctions.H"
#include "Hoeb_ExactSolutions.H"
//this one defines HOEB_MAX_ORDER
#include "Hoeb_Utilities.H"
#include <iomanip>


/****/
void
getKappaLphi(EBLevelBoxData<CELL, 1>                                            &  a_klp,
             const EBLevelBoxData<CELL, 1>                                      &  a_phi,
             const shared_ptr<LevelData<EBGraph> >                              &  a_graphs,
             const Chombo4::DisjointBoxLayout                                   &  a_grids,
             const Chombo4::Box                                                 &  a_domain,
             const Real                                                         &  a_dx,
             const shared_ptr<EBDictionary<HOEB_MAX_ORDER, Real, CELL, CELL> >  &  a_dictionary,
             const shared_ptr< GeometryService<HOEB_MAX_ORDER> >                &  a_geoserv)
{

  Chombo4::DataIterator dit = a_grids.dataIterator();
  //register it for every box
  for(unsigned int ibox = 0; ibox < dit.size(); ++ibox)
  {
    string stencilName;
    string ebbcName;
    vector<     EBIndex<CELL>  >          dstVoFs;
    vector<LocalStencil<CELL, Real> >     stencil;
    Proto::Box                            srcValid;
    Proto::Box                            dstValid;
    Proto::Box                            srcDomain;
    Proto::Box                            dstDomain;
    Point                                 srcGhost;
    Point                                 dstGhost;
    bool                                  needDiagonalWeights;
    
    hoeb::
      dharshiLaplStencil(stencilName,        
                         ebbcName,           
                         dstVoFs,            
                         stencil,            
                         srcValid,           
                         dstValid,           
                         srcDomain,          
                         dstDomain,          
                         srcGhost,           
                         dstGhost,
                         needDiagonalWeights,
                         a_geoserv,
                         a_grids,
                         a_domain,
                         a_dx,
                         ibox);
  
    ///registering stencil
    a_dictionary->registerStencil(stencilName,        
                                  ebbcName,           
                                  dstVoFs,            
                                  stencil,            
                                  srcValid,           
                                  dstValid,           
                                  srcDomain,          
                                  dstDomain,          
                                  srcGhost,           
                                  dstGhost,
                                  needDiagonalWeights,
                                  ibox);
    //now apply it
    for(unsigned int ibox = 0; ibox < dit.size(); ++ibox)
    {
      auto      & lphfab = a_klp[dit[ibox]];
      const auto& phifab = a_phi[dit[ibox]];
      auto stencil = a_dictionary->getEBStencil(stencilName, ebbcName, a_domain, a_domain, ibox);
      //set resc = Ave(resf) (true is initToZero)
      stencil->apply(lphfab, phifab,  true, 1.0);
    }
  }
}
/*******/ 
PROTO_KERNEL_START 
void subtractionTractionF(Var<Real, 1>    a_error,
                          Var<Real, 1>    a_klpFtoC,
                          Var<Real, 1>    a_klpCoar)
{
  a_error(0) = a_klpFtoC(0) - a_klpCoar(0);
}
PROTO_KERNEL_END(subtractionTractionF, subtractionTraction)
/****/
void
getKLPhiError(EBLevelBoxData<CELL,   1>                                           &  a_errCoar, 
              const shared_ptr<LevelData<EBGraph> >                               &  a_graphsFine,
              const Chombo4::DisjointBoxLayout                                    &  a_gridsFine,
              const Chombo4::Box                                                  &  a_domFine,
              const Real                                                          &  a_dxFine,
              const shared_ptr<LevelData<EBGraph> >                               &  a_graphsCoar,
              const Chombo4::DisjointBoxLayout                                    &  a_gridsCoar,
              const Chombo4::Box                                                  &  a_domCoar,
              const Real                                                          &  a_dxCoar,
              const shared_ptr<EBDictionary<HOEB_MAX_ORDER, Real, CELL, CELL> >   &  a_dictionary,
              const shared_ptr< GeometryService<HOEB_MAX_ORDER> >                 &  a_geoserv)
{
  IntVect dataGhostIV =   4*IntVect::Unit;
  EBLevelBoxData<CELL,   1>  phiFine(a_gridsFine, dataGhostIV, a_graphsFine);
  EBLevelBoxData<CELL,   1>  phiCoar(a_gridsCoar, dataGhostIV, a_graphsCoar);
  EBLevelBoxData<CELL,   1>  klpFine(a_gridsFine, dataGhostIV, a_graphsFine);
  EBLevelBoxData<CELL,   1>  klpCoar(a_gridsCoar, dataGhostIV, a_graphsCoar);
  EBLevelBoxData<CELL,   1>  klpFtoC(a_gridsCoar, dataGhostIV, a_graphsCoar);
  
  hoeb::fillPhi(phiFine, a_graphsFine, a_gridsFine, a_domFine, a_dxFine, a_geoserv);
  hoeb::fillPhi(phiCoar, a_graphsCoar, a_gridsCoar, a_domCoar, a_dxCoar, a_geoserv);

  getKappaLphi(klpFine, phiFine, a_graphsFine, a_gridsFine, a_domFine, a_dxFine, a_dictionary, a_geoserv);
  getKappaLphi(klpCoar, phiCoar, a_graphsCoar, a_gridsCoar, a_domCoar, a_dxCoar, a_dictionary, a_geoserv);

  hoeb::restrictKappaLphi(klpFtoC, klpFine,
                    a_graphsCoar, a_gridsCoar, a_domCoar, a_dxCoar,
                    a_graphsFine, a_gridsFine, a_domFine, a_dxFine,                    
                    a_dictionary, a_geoserv);

  //error = Ave(klphifine) - klphicoar 
  Chombo4::DataIterator dit = a_gridsCoar.dataIterator();
  for(unsigned int ibox = 0; ibox < dit.size(); ibox++)
  {
    auto& ftocfab =   klpFtoC[dit[ibox]];
    auto& coarfab =   klpCoar[dit[ibox]];
    auto& errfab  = a_errCoar[dit[ibox]];
    auto  inputbx = ftocfab.inputBox();
    auto  validbx = (*a_graphsCoar)[dit[ibox]].validBox();
    ebforall(inputbx, subtractionTraction, validbx, errfab, ftocfab, coarfab);
  }
}

/****/
int
runTest()
{
  Real coveredval = -1;
  int nx      = 32;
  int maxGrid = 32;
    
  ParmParse pp;

  pp.get("nx"        , nx);
  pp.get("max_grid"  , maxGrid);
  pp.get("coveredval", coveredval);         


  pout() << "nx"        << " = " <<  nx         << endl;
  pout() << "max_grid"  << " = " <<  maxGrid    << endl;
  pout() << "coveredval"<< " = " <<  coveredval << endl;         
  

  IntVect domLo = IntVect::Zero;
  IntVect domHi  = (nx - 1)*IntVect::Unit;

  Chombo4::ProblemDomain domain(domLo, domHi);

  Vector<Chombo4::DisjointBoxLayout> vecgrids;
  pout() << "making grids" << endl;
  GeometryService<2>::generateGrids(vecgrids, domain.domainBox(), maxGrid);

  Chombo4::DisjointBoxLayout gridsFine = vecgrids[0];
  Chombo4::DisjointBoxLayout gridsMedi = vecgrids[1];
  Chombo4::DisjointBoxLayout gridsCoar = vecgrids[2];
  IntVect dataGhostIV =   4*IntVect::Unit;
  Point   dataGhostPt = ProtoCh::getPoint(dataGhostIV); 
  int geomGhost = 4;
  RealVect origin = RealVect::Zero();
  
  shared_ptr<BaseIF>    impfunc = hoeb::getImplicitFunction();
  pout() << "defining geometry" << endl;
  Real dx = 1.0/(Real(nx));
  vector<Chombo4::Box>    vecdomain(vecgrids.size(), domain.domainBox());
  vector<Real>   vecdx    (vecgrids.size(), dx);
  for(int ilev = 1; ilev < vecgrids.size(); ilev++)
  {
    vecdomain[ilev] = coarsen(vecdomain[ilev-1], 2);
    vecdx    [ilev] =           2*vecdx[ilev-1];
  }
  
  Real dxFine = vecdx[0];
  Real dxMedi = vecdx[1];
  Real dxCoar = vecdx[2];
  GeometryService<HOEB_MAX_ORDER>* geomptr =
    new GeometryService<HOEB_MAX_ORDER>
    (impfunc, origin, dxFine, domain.domainBox(), vecgrids, geomGhost);
  
  shared_ptr< GeometryService<HOEB_MAX_ORDER> >  geoserv(geomptr);

  pout() << "making dictionary" << endl;

  Chombo4::Box domFine = vecdomain[0];
  Chombo4::Box domMedi = vecdomain[1];
  Chombo4::Box domCoar = vecdomain[2];

  shared_ptr<LevelData<EBGraph> > graphsFine = geoserv->getGraphs(domFine);
  shared_ptr<LevelData<EBGraph> > graphsMedi = geoserv->getGraphs(domMedi);
  shared_ptr<LevelData<EBGraph> > graphsCoar = geoserv->getGraphs(domCoar);

  
  shared_ptr<EBDictionary<HOEB_MAX_ORDER, Real, CELL, CELL> >  dictionary
    (new     EBDictionary<HOEB_MAX_ORDER, Real, CELL, CELL>
     (geoserv, vecgrids, vecdomain, vecdx, dataGhostPt));

  EBLevelBoxData<CELL,   1>  errMedi(gridsMedi, dataGhostIV, graphsMedi);
  EBLevelBoxData<CELL,   1>  errCoar(gridsCoar, dataGhostIV, graphsCoar);
  

  getKLPhiError(errMedi, 
                graphsFine, gridsFine, domFine, dxFine,
                graphsMedi, gridsMedi, domMedi, dxMedi,
                dictionary, geoserv);

  getKLPhiError(errCoar, 
                graphsMedi, gridsMedi, domMedi, dxMedi,
                graphsCoar, gridsCoar, domCoar, dxCoar,
                dictionary, geoserv);


  //Norm!
  Real normMedi = errMedi.maxNorm(0);
  Real normCoar = errCoar.maxNorm(0);

  Real tol = 1.0e-16;
  Real order = 0;
  if((normCoar > tol) && (normMedi > tol))
  {
    order = log(normCoar/normMedi)/log(2.0);
  }
  pout() << "||klphi errMedi||_max = " << normMedi << std::endl;
  pout() << "||klphi errCoar||_max = " << normCoar << std::endl;
  pout() << "Richardson truncation error order for kappa(L(phi))= " << order << std::endl;
  return 0;
}


int main(int a_argc, char* a_argv[])
{
#ifdef CH_USE_PETSC  
  //because of some kind of solipsistic madness, PetscInitialize calls MPI_INIT
   PetscInt ierr = PetscInitialize(&a_argc, &a_argv, "./.petscrc",PETSC_NULL); CHKERRQ(ierr);
#else  
#ifdef CH_MPI
  MPI_Init(&a_argc, &a_argv);
  pout() << "MPI INIT called" << std::endl;
#endif
#endif
  //needs to be called after MPI_Init
  CH_TIMER_SETFILE("trunc.time.table");
  {
    if (a_argc < 2)
    {
      cerr<< " usage " << a_argv[0] << " <input_file_name> " << endl;
      exit(0);
    }
    char* in_file = a_argv[1];
    ParmParse  pp(a_argc-2,a_argv+2,NULL,in_file);
    runTest();
  }

  pout() << "printing time table " << endl;
  CH_TIMER_REPORT();
#ifdef CH_MPI
#ifdef CH_USE_PETSC
  pout() << "about to call petsc Finalize" << std::endl;
  PetscFinalize();
#else  
  pout() << "about to call MPI Finalize" << std::endl;
  MPI_Finalize();
#endif
#endif
  return 0;
}
