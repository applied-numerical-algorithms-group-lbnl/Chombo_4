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
#include <iomanip>

/****/
void
defineStencil(string                                   & a_stencilName,
              string                                   & a_ebbcName,
              vector<EBIndex<CELL> >                   & a_dstVoFs,
              vector<LocalStencil<CELL, Real> >        & a_stencil,
              Stencil<Real>                            & a_regStencilInterior,
              Proto::Box                               & a_regApplyBox,
              Proto::Box                               & a_srcValid,
              Proto::Box                               & a_dstValid,
              Proto::Box                               & a_srcDomain,
              Proto::Box                               & a_dstDomain,
              Proto::Point                             & a_srcGhost,
              Proto::Point                             & a_dstGhost,
              bool                                     & a_irregOnly,
              bool                                     & a_needDiagonalWeights)
{
  ParmParse pp;
    
  int dombc = 1;
  int ebbc  = 1;
  pp.get("domainBC"  , dombc);
  pp.get("EBBC"      , ebbc);
  pout() << "domainBC"  << " = " <<  dombc      << endl;
  pout() << "EBBC"      << " = " <<  ebbc       << endl;
  
  a_needDiagonalWeights = true;
  a_irregOnly           = true;
}
/**
{
  string stencilName;
  string ebbcName;
  vector<     EBIndex<CELL>  >          dstVoFs;
  vector<LocalStencil<CELL, Real> >     stencil;
  Stencil<Real>                         regStencilInterior;
  Proto::Box                            regApplyBox;
  Proto::Box                            srcValid;
  Proto::Box                            dstValid;
  Proto::Box                            srcDomain;
  Proto::Box                            dstDomain;
  Point                                 srcGhost;
  Point                                 dstGhost;
  bool                                  irregOnly;
  bool                                  needDiagonalWeights;

  defineStencil(stencilName,        
                ebbcName,           
                dstVoFs,            
                stencil,            
                regStencilInterior, 
                regApplyBox,        
                srcValid,           
                dstValid,           
                srcDomain,          
                dstDomain,          
                srcGhost,           
                dstGhost,
                irregOnly,
                needDiagonalWeights);
  
  ///registering stencil
  dictionary->registerStencil(stencilName,        
                              ebbcName,           
                              dstVoFs,            
                              stencil,            
                              regStencilInterior, 
                              regApplyBox,        
                              srcValid,           
                              dstValid,           
                              srcDomain,          
                              dstDomain,          
                              srcGhost,           
                              dstGhost,
                              irregOnly,
                              needDiagonalWeights);
  
  
  Chombo4::Box dombox = domain.domainBox();
  pout() << "making data" << endl;
  EBLevelBoxData<CELL,   1>  error(grids, dataGhostIV, graphs);

  fillData(phiexac, lphexac); 


  pout() << "writing to file " << endl;
  
  
  pout() << "exiting " << endl;
  return 0;
}
**/
/****/


//int
//runTest()
//{
//  Real coveredval = -1;
//  int nx      = 32;
//  int maxGrid = 32;
//  Real alpha = 1.0;
//  Real beta = -0.001;
//    
//
//  ParmParse pp;
//
//  pp.get("nx"        , nx);
//  pp.get("max_grid"  , maxGrid);
//  pp.get("alpha"     , alpha);
//  pp.get("beta"      , beta);
//  pp.get("coveredval", coveredval);         
//
//
//  pout() << "nx"        << " = " <<  nx         << endl;
//  pout() << "max_grid"  << " = " <<  maxGrid    << endl;
//  pout() << "alpha"     << " = " <<  alpha      << endl;
//  pout() << "beta"      << " = " <<  beta       << endl;
//  pout() << "coveredval"<< " = " <<  coveredval << endl;         
//  
//
//  IntVect domLo = IntVect::Zero;
//  IntVect domHi  = (nx - 1)*IntVect::Unit;
//
//// EB and periodic do not mix
//  Chombo4::ProblemDomain domain(domLo, domHi);
//
//  Vector<Chombo4::DisjointBoxLayout> vecgrids;
//  pout() << "making grids" << endl;
//  GeometryService<2>::generateGrids(vecgrids, domain.domainBox(), maxGrid);
//
//  Chombo4::DisjointBoxLayout grids = vecgrids[0];
//  grids.printBalance();
//
//  IntVect dataGhostIV =   2*IntVect::Unit;
//  Point   dataGhostPt = ProtoCh::getPoint(dataGhostIV); 
//  int geomGhost = 4;
//  RealVect origin = RealVect::Zero();
//  Real dx = 1.0/nx;
////  Real dx = 1.0;
//  shared_ptr<BaseIF>    impfunc = getImplicitFunction();
////  Bx domainpr = getProtoBox(domain.domainBox());
//
//  pout() << "defining geometry" << endl;
//  GeometryService<2>* geomptr = new GeometryService<2>(impfunc, origin, dx, domain.domainBox(), vecgrids, geomGhost);
////  GeometryService<2>* geomptr = new GeometryService<2>(impfunc, origin, dx, domain.domainBox(), vecgrids[0], geomGhost);
//  shared_ptr< GeometryService<2> >  geoserv(geomptr);
//
//  pout() << "making dictionary" << endl;
//
//  vector<Chombo4::Box>    vecdomain(vecgrids.size(), domain.domainBox());
//  vector<Real>   vecdx    (vecgrids.size(), dx);
//  for(int ilev = 1; ilev < vecgrids.size(); ilev++)
//  {
//    vecdomain[ilev] = coarsen(vecdomain[ilev-1], 2);
//    vecdx    [ilev] =           2*vecdx[ilev-1];
//  }
//  shared_ptr<EBDictionary<2, Real, CELL, CELL> > 
//    dictionary(new EBDictionary<2, Real, CELL, CELL>(geoserv, vecgrids, vecdomain, vecdx, dataGhostPt));
//
//  string stencilName;
//  string ebbcName;
//  vector<     EBIndex<CELL>  >          dstVoFs;
//   vector<LocalStencil<CELL, Real> >     stencil;
//  Stencil<Real>                         regStencilInterior;
//  Proto::Box                            regApplyBox;
//  Proto::Box                            srcValid;
//  Proto::Box                            dstValid;
//  Proto::Box                            srcDomain;
//  Proto::Box                            dstDomain;
//  Point                                 srcGhost;
//  Point                                 dstGhost;
//  bool                                  irregOnly;
//  bool                                  needDiagonalWeights;
//
//  defineStencil(stencilName,        
//                ebbcName,           
//                dstVoFs,            
//                stencil,            
//                regStencilInterior, 
//                regApplyBox,        
//                srcValid,           
//                dstValid,           
//                srcDomain,          
//                dstDomain,          
//                srcGhost,           
//                dstGhost,
//                irregOnly,
//                needDiagonalWeights);
//  
//  ///registering stencil
//  dictionary->registerStencil(stencilName,        
//                              ebbcName,           
//                              dstVoFs,            
//                              stencil,            
//                              regStencilInterior, 
//                              regApplyBox,        
//                              srcValid,           
//                              dstValid,           
//                              srcDomain,          
//                              dstDomain,          
//                              srcGhost,           
//                              dstGhost,
//                              irregOnly,
//                              needDiagonalWeights);
//  
//  
//  Chombo4::Box dombox = domain.domainBox();
//  shared_ptr<LevelData<EBGraph> > graphs = geoserv->getGraphs(dombox);
//
//  pout() << "making data" << endl;
//  EBLevelBoxData<CELL,   1>  phiexac(grids, dataGhostIV, graphs);
//  EBLevelBoxData<CELL,   1>  lphcalc(grids, dataGhostIV, graphs);
//  EBLevelBoxData<CELL,   1>  lphexac(grids, dataGhostIV, graphs);
//  EBLevelBoxData<CELL,   1>  error(grids, dataGhostIV, graphs);
//
//  fillData(phiexac, lphexac); 
//
//
//  pout() << "writing to file " << endl;
//  
//  
//  pout() << "exiting " << endl;
//  return 0;
//}
//



/****/
shared_ptr<BaseIF> 
getImplicitFunction()
{
  RealVect center = 0.5*RealVect::Unit();
  Real radius = 0.1;
  bool inside = false;
  ParmParse pp;
  std::vector<Real> centvec;
  pp.get("radius", radius);
  pp.get("inside", inside);
  pp.getarr("center", centvec, 0, DIM);
  for(int idir = 0; idir < DIM; idir++)
  {
    center[idir] = centvec[idir];
  }
  SimpleSphereIF* sphereptr = new SimpleSphereIF(center, radius, inside);
  shared_ptr<BaseIF> retval(static_cast<BaseIF*>(sphereptr));
  return retval;
}

/****/
void
fillPhi(EBLevelBoxData<CELL, 1>                        &  a_phi,
        const shared_ptr<LevelData<EBGraph> >          &  a_graphs,
        const Chombo4::DisjointBoxLayout               &  a_grids,
        const Chombo4::Box                             &  a_domFine,
        const Real                                     &  a_dx,
        const shared_ptr< GeometryService<2> >         &  a_geoserv)
{
        
}
/****/
void
restrictKappaLphi(EBLevelBoxData<CELL, 1>                                &  a_klpFToC,
                  const EBLevelBoxData<CELL, 1>                          &  a_klpFine,
                  const shared_ptr<LevelData<EBGraph> >                  &  a_graphsFine,
                  const Chombo4::DisjointBoxLayout                       &  a_gridsFine,
                  const Chombo4::Box                                     &  a_domFine,
                  const Real                                             &  a_dxFine,
                  const shared_ptr<LevelData<EBGraph> >                  &  a_graphsCoar,
                  const Chombo4::DisjointBoxLayout                       &  a_gridsCoar,
                  const Chombo4::Box                                     &  a_domCoar,
                  const Real                                             &  a_dxCoar,
                  const shared_ptr<EBDictionary<2, Real, CELL, CELL> >   &  a_dictionary,
                  const shared_ptr< GeometryService<2> >                 &  a_geoserv)
{
        
}
/****/
void
getKappaLphi(EBLevelBoxData<CELL, 1>                                &  a_klp,
             const EBLevelBoxData<CELL, 1>                          &  a_phi,
             const shared_ptr<LevelData<EBGraph> >                  &  a_graphs,
             const Chombo4::DisjointBoxLayout                       &  a_grids,
             const Chombo4::Box                                     &  a_domFine,
             const Real                                             &  a_dx,
             const shared_ptr<EBDictionary<2, Real, CELL, CELL> >   &  a_dictionary,
             const shared_ptr< GeometryService<2> >                 &  a_geoserv)
{
        
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
getKLPhiError(EBLevelBoxData<CELL,   1>                              &  a_errCoar, 
              const shared_ptr<LevelData<EBGraph> >                  &  a_graphsFine,
              const Chombo4::DisjointBoxLayout                       &  a_gridsFine,
              const Chombo4::Box                                     &  a_domFine,
              const Real                                             &  a_dxFine,
              const shared_ptr<LevelData<EBGraph> >                  &  a_graphsCoar,
              const Chombo4::DisjointBoxLayout                       &  a_gridsCoar,
              const Chombo4::Box                                     &  a_domCoar,
              const Real                                             &  a_dxCoar,
              const shared_ptr<EBDictionary<2, Real, CELL, CELL> >   &  a_dictionary,
              const shared_ptr< GeometryService<2> >                 &  a_geoserv)
{
  IntVect dataGhostIV =   4*IntVect::Unit;
  EBLevelBoxData<CELL,   1>  phiFine(a_gridsFine, dataGhostIV, a_graphsFine);
  EBLevelBoxData<CELL,   1>  phiCoar(a_gridsCoar, dataGhostIV, a_graphsCoar);
  EBLevelBoxData<CELL,   1>  klpFine(a_gridsFine, dataGhostIV, a_graphsFine);
  EBLevelBoxData<CELL,   1>  klpCoar(a_gridsCoar, dataGhostIV, a_graphsCoar);
  EBLevelBoxData<CELL,   1>  klpFtoC(a_gridsCoar, dataGhostIV, a_graphsCoar);
  
  fillPhi(phiFine, a_graphsFine, a_gridsFine, a_domFine, a_dxFine, a_geoserv);
  fillPhi(phiCoar, a_graphsCoar, a_gridsCoar, a_domCoar, a_dxCoar, a_geoserv);

  getKappaLphi(klpFine, phiFine, a_graphsFine, a_gridsFine, a_domFine, a_dxFine, a_dictionary, a_geoserv);
  getKappaLphi(klpCoar, phiCoar, a_graphsCoar, a_gridsCoar, a_domCoar, a_dxCoar, a_dictionary, a_geoserv);

  restrictKappaLphi(klpFtoC, klpFine,
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
  Real alpha = 1.0;
  Real beta = -0.001;
    

  ParmParse pp;

  pp.get("nx"        , nx);
  pp.get("max_grid"  , maxGrid);
  pp.get("alpha"     , alpha);
  pp.get("beta"      , beta);
  pp.get("coveredval", coveredval);         


  pout() << "nx"        << " = " <<  nx         << endl;
  pout() << "max_grid"  << " = " <<  maxGrid    << endl;
  pout() << "alpha"     << " = " <<  alpha      << endl;
  pout() << "beta"      << " = " <<  beta       << endl;
  pout() << "coveredval"<< " = " <<  coveredval << endl;         
  

  IntVect domLo = IntVect::Zero;
  IntVect domHi  = (nx - 1)*IntVect::Unit;

// EB and periodic do not mix
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
  
  shared_ptr<BaseIF>    impfunc = getImplicitFunction();
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
  GeometryService<2>* geomptr =
    new GeometryService<2>(impfunc, origin, dxFine, domain.domainBox(), vecgrids, geomGhost);
  shared_ptr< GeometryService<2> >  geoserv(geomptr);

  pout() << "making dictionary" << endl;

  Chombo4::Box domFine = vecdomain[0];
  Chombo4::Box domMedi = vecdomain[1];
  Chombo4::Box domCoar = vecdomain[2];

  shared_ptr<LevelData<EBGraph> > graphsFine = geoserv->getGraphs(domFine);
  shared_ptr<LevelData<EBGraph> > graphsMedi = geoserv->getGraphs(domMedi);
  shared_ptr<LevelData<EBGraph> > graphsCoar = geoserv->getGraphs(domCoar);

  
  shared_ptr<EBDictionary<2, Real, CELL, CELL> >  dictionary
    (new EBDictionary<2, Real, CELL, CELL>(geoserv, vecgrids, vecdomain, vecdx, dataGhostPt));

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
