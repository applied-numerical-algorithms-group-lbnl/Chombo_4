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

#include "DebugFunctions.H"
#include "Hoeb_ExactSolutions.H"
//this one defines HOEB_MAX_ORDER
#include "Hoeb_Utilities.H"
#include "Hoeb_LAPACKMatrix.H"
#include <iomanip>

#ifdef CH_USE_DOUBLE
Real g_tol = 1.0e-12;
#else
Real g_tol = 1.0e-6;
#endif
using std::endl;
using hoeb::LAPACKMatrix;
using Chombo4::pout;
//This is the dirty work of getting functions into the symbol table for debugging.
unsigned int
printNeighborhood(hoeb::Neighborhood<CELL>*  a_blerg)
{
  if(a_blerg != NULL)
  {
    a_blerg->poutAll();
  }
  return 0;
}
unsigned int
printCompositeStencil(hoeb::CompositeStencil*  a_blerg)
{
  if(a_blerg != NULL)
  {
    a_blerg->poutAll();
  }
  return 0;
}

/****/
void
getKappaLphiHomogeneous(EBLevelBoxData<CELL, 1>                                            &  a_klp,
                        const EBLevelBoxData<CELL, 1>                                      &  a_phi,
                        const shared_ptr<hoeb::graph_distrib_t>                              &  a_graphs,
                        const Chombo4::DisjointBoxLayout                                   &  a_grids,
                        const Chombo4::Box                                                 &  a_domain,
                        const Real                                                         &  a_dx,
                        const shared_ptr<EBDictionary<HOEB_MAX_ORDER, Real, CELL, CELL> >  &  a_dictionary,
                        const shared_ptr< GeometryService<HOEB_MAX_ORDER> >                &  a_geoserv,
                        string a_stencilname)
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
    
    if(a_stencilname == string("Devendran"))
    {
      hoeb::
        getHomogeneousDharshiStencil(stencilName,        
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
    }
    else if(a_stencilname == string("Schwartz"))
    {
      hoeb::
        schwartzLaplStencil(stencilName,        
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
    }

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
typedef CH4_Data_Choreography::DistributedData<EBGraph>   graph_distrib_t;
/*******/ 
void
getKappaLphi(EBLevelBoxData<CELL, 1>                                            &  a_klp,
             const EBLevelBoxData<CELL, 1>                                      &  a_phi,
             const shared_ptr<graph_distrib_t>                                  &  a_graphs,
             const Chombo4::DisjointBoxLayout                                   &  a_grids,
             const Chombo4::Box                                                 &  a_domain,
             const Real                                                         &  a_dx,
             const shared_ptr<EBDictionary<HOEB_MAX_ORDER, Real, CELL, CELL> >  &  a_dictionary,
             const shared_ptr< GeometryService<HOEB_MAX_ORDER> >                &  a_geoserv,
             string                                                                a_stencilname)
{

  typedef GraphConstructorFactory<EBHostData<CELL, Real, 1> > hostfactorycell_t;
  typedef CH4_Data_Choreography::DistributedData<EBHostData<CELL, Real, 1> > cell_distrib_t;
  
  IntVect dataGhostIV = a_klp.ghostVect();
  cell_distrib_t  hostklp(a_grids, dataGhostIV, hostfactorycell_t(a_graphs));


  Chombo4::DataIterator dit = a_grids.dataIterator();
  for(unsigned int ibox = 0; ibox < dit.size(); ++ibox)
  {
    auto chomboBox =a_grids[dit[ibox]];
    Proto::Box bx = getProtoBox(chomboBox);
    const auto & graph= (*a_graphs)[dit[ibox]];
    for(auto bit = bx.begin(); bit != bx.end(); ++bit)
    {
      auto pt = *bit;
      auto vofs = graph.getVoFs(pt);
      for(int ivof = 0; ivof < vofs.size(); ivof++)
      {
        const auto& vof = vofs[ivof];
        Real vofval = hoeb::getKappaLphiVoF(vof,
                                            dit[ibox],
                                            a_phi,        
                                            a_graphs,     
                                            a_grids,      
                                            a_domain,     
                                            a_dx,         
                                            a_dictionary, 
                                            a_geoserv,    
                                            a_stencilname);
        auto& klpfab = hostklp[dit[ibox]];
        klpfab(vof, 0) = vofval;
      }
    }
  }
  EBLevelBoxData<CELL,1>::copyToDevice(a_klp, hostklp); 
}
/*******/ 
PROTO_KERNEL_START 
void subtractionTractionF(Var<Real, 1>    a_error,
                          Var<Real, 1>    a_klpFtoC,
                          Var<Real, 1>    a_klpCoar)
{
  Real ftoc  = a_klpFtoC(0);
  Real coar  = a_klpCoar(0);
  Real errr  = ftoc - coar;
//  pout() << ftoc << " - " << coar << " = " << errr << endl;
  a_error(0) = errr;
}
PROTO_KERNEL_END(subtractionTractionF, subtractionTraction)
/****/
void
getKLPhiError(EBLevelBoxData<CELL,   1>                                           &  a_errCoar, 
              const shared_ptr<graph_distrib_t    >                               &  a_graphsFine,
              const Chombo4::DisjointBoxLayout                                    &  a_gridsFine,
              const Chombo4::Box                                                  &  a_domFine,
              const Real                                                          &  a_dxFine,
              const shared_ptr<graph_distrib_t    >                               &  a_graphsCoar,
              const Chombo4::DisjointBoxLayout                                    &  a_gridsCoar,
              const Chombo4::Box                                                  &  a_domCoar,
              const Real                                                          &  a_dxCoar,
              const shared_ptr<EBDictionary<HOEB_MAX_ORDER, Real, CELL, CELL> >   &  a_dictionary,
              const shared_ptr< GeometryService<HOEB_MAX_ORDER> >                 &  a_geoserv,
              string a_stencilname)
{
  ParmParse pp;
  int nghost = 6;
  IntVect dataGhostIV = nghost*IntVect::Unit;
  EBLevelBoxData<CELL,   1>  phiFine(a_gridsFine, dataGhostIV, a_graphsFine);
  EBLevelBoxData<CELL,   1>  phiCoar(a_gridsCoar, dataGhostIV, a_graphsCoar);
  EBLevelBoxData<CELL,   1>  klpFine(a_gridsFine, dataGhostIV, a_graphsFine);
  EBLevelBoxData<CELL,   1>  klpCoar(a_gridsCoar, dataGhostIV, a_graphsCoar);
  EBLevelBoxData<CELL,   1>  klpFtoC(a_gridsCoar, dataGhostIV, a_graphsCoar);
  
  hoeb::fillPhi(phiFine, a_graphsFine, a_gridsFine, a_domFine, a_dxFine, a_geoserv);
  hoeb::fillPhi(phiCoar, a_graphsCoar, a_gridsCoar, a_domCoar, a_dxCoar, a_geoserv);

  getKappaLphi(klpFine, phiFine, a_graphsFine, a_gridsFine, a_domFine, a_dxFine, a_dictionary, a_geoserv, a_stencilname);
  getKappaLphi(klpCoar, phiCoar, a_graphsCoar, a_gridsCoar, a_domCoar, a_dxCoar, a_dictionary, a_geoserv, a_stencilname);

  hoeb::restrictKappaLphi(klpFtoC, klpFine,
                          a_graphsFine, a_gridsFine, a_domFine, a_dxFine,                    
                          a_graphsCoar, a_gridsCoar, a_domCoar, a_dxCoar,
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
void
getICError(EBLevelBoxData<CELL,   1>                                           &  a_errCoar, 
           const shared_ptr<graph_distrib_t>                                   &  a_graphsFine,
           const Chombo4::DisjointBoxLayout                                    &  a_gridsFine,
           const Chombo4::Box                                                  &  a_domFine,
           const Real                                                          &  a_dxFine,
           const shared_ptr<graph_distrib_t>                                   &  a_graphsCoar,
           const Chombo4::DisjointBoxLayout                                    &  a_gridsCoar,
           const Chombo4::Box                                                  &  a_domCoar,
           const Real                                                          &  a_dxCoar,
           const shared_ptr<EBDictionary<HOEB_MAX_ORDER, Real, CELL, CELL> >   &  a_dictionary,
           const shared_ptr< GeometryService<HOEB_MAX_ORDER> >                 &  a_geoserv)
{
  ParmParse pp;
  int nghost = 6;

  IntVect dataGhostIV =   nghost*IntVect::Unit;
  EBLevelBoxData<CELL,   1>  phiFine(a_gridsFine, dataGhostIV, a_graphsFine);
  EBLevelBoxData<CELL,   1>  phiCoar(a_gridsCoar, dataGhostIV, a_graphsCoar);
  EBLevelBoxData<CELL,   1>  phiFtoC(a_gridsCoar, dataGhostIV, a_graphsCoar);
  
  hoeb::fillPhi(phiFine, a_graphsFine, a_gridsFine, a_domFine, a_dxFine, a_geoserv);
  hoeb::fillPhi(phiCoar, a_graphsCoar, a_gridsCoar, a_domCoar, a_dxCoar, a_geoserv);

  hoeb::restrictPhi(phiFtoC, phiFine,
                    a_graphsFine, a_gridsFine, a_domFine, a_dxFine,                    
                    a_graphsCoar, a_gridsCoar, a_domCoar, a_dxCoar,
                    a_dictionary, a_geoserv);

  //error = Ave(phifine) - phicoar 
  Chombo4::DataIterator dit = a_gridsCoar.dataIterator();
  for(unsigned int ibox = 0; ibox < dit.size(); ibox++)
  {
    auto& ftocfab =   phiFtoC[dit[ibox]];
    auto& coarfab =   phiCoar[dit[ibox]];
    auto& errfab  = a_errCoar[dit[ibox]];
    auto  inputbx = ftocfab.inputBox();
    auto  validbx = (*a_graphsCoar)[dit[ibox]].validBox();
    ebforall(inputbx, subtractionTraction, validbx, errfab, ftocfab, coarfab);
  }
}

/****/
unsigned int
runTruncationErrorTest(string a_stencilname)
{
  using Chombo4::pout;
  Real coveredval = -1;
  int nx      = 32;
  int maxGrid = 32;
    
  ParmParse pp;

  pp.get("nx"        , nx);
  pp.get("maxGrid"  , maxGrid);
  pp.get("coveredval", coveredval);         


  pout() << "nx"        << " = " <<  nx         << endl;
  pout() << "max_grid"  << " = " <<  maxGrid    << endl;
  pout() << "coveredval"<< " = " <<  coveredval << endl;         
  

  IntVect domLo = IntVect::Zero;
  IntVect domHi  = (nx - 1)*IntVect::Unit;

  Chombo4::ProblemDomain domain(domLo, domHi);

  vector<Chombo4::DisjointBoxLayout> vecgrids;
  pout() << "making grids" << endl;
  GeometryService<2>::generateGrids(vecgrids, domain.domainBox(), maxGrid);

  Chombo4::DisjointBoxLayout gridsFine = vecgrids[0];
  Chombo4::DisjointBoxLayout gridsMedi = vecgrids[1];
  Chombo4::DisjointBoxLayout gridsCoar = vecgrids[2];
  IntVect dataGhostIV =   6*IntVect::Unit;
  Point   dataGhostPt = ProtoCh::getPoint(dataGhostIV); 
  int geomGhost = 6;
  RealVect origin = RealVect::Zero();
  
  shared_ptr<BaseIF>    impfunc = hoeb::getImplicitFunction();
  pout() << "defining geometry" << endl;
  
  //Real dx = 1.0/(Real(nx));
  Real dx = 32.0/(Real(nx));
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

  shared_ptr<graph_distrib_t> graphsFine = geoserv->getGraphs(domFine);
  shared_ptr<graph_distrib_t> graphsMedi = geoserv->getGraphs(domMedi);
  shared_ptr<graph_distrib_t> graphsCoar = geoserv->getGraphs(domCoar);

  shared_ptr<EBDictionary<HOEB_MAX_ORDER, Real, CELL, CELL> >  dictionary
    (new     EBDictionary<HOEB_MAX_ORDER, Real, CELL, CELL>
     (geoserv, vecgrids, vecdomain, vecdx, dataGhostPt));

  EBLevelBoxData<CELL,   1>  errMedi(gridsMedi, dataGhostIV, graphsMedi);
  EBLevelBoxData<CELL,   1>  errCoar(gridsCoar, dataGhostIV, graphsCoar);
  

  getKLPhiError(errMedi, 
                graphsFine, gridsFine, domFine, dxFine,
                graphsMedi, gridsMedi, domMedi, dxMedi,
                dictionary, geoserv, a_stencilname);

  getKLPhiError(errCoar, 
                graphsMedi, gridsMedi, domMedi, dxMedi,
                graphsCoar, gridsCoar, domCoar, dxCoar,
                dictionary, geoserv, a_stencilname);


  pout() << "writing to file" << endl;
  string fileCoar = a_stencilname + string("errCoar.hdf5");
  string fileMedi = a_stencilname + string("errMedi.hdf5");
  
  errCoar.writeToFileHDF5(fileCoar, 0.0);
  errMedi.writeToFileHDF5(fileMedi, 0.0);
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
  pout() << a_stencilname << " Richardson truncation error order for kappa(L(phi))= " << order << std::endl;
  return 0;
}

/****/
unsigned int
runInitialConditionTest()
{
  using Chombo4::pout;
  Real coveredval = -1;
  int nx      = 32;
  int maxGrid = 32;
    
  ParmParse pp;

  pp.get("nx"        , nx);
  pp.get("maxGrid"  , maxGrid);
  pp.get("coveredval", coveredval);         


  pout() << "nx"        << " = " <<  nx         << endl;
  pout() << "max_grid"  << " = " <<  maxGrid    << endl;
  pout() << "coveredval"<< " = " <<  coveredval << endl;         
  

  IntVect domLo = IntVect::Zero;
  IntVect domHi  = (nx - 1)*IntVect::Unit;

  Chombo4::ProblemDomain domain(domLo, domHi);

  vector<Chombo4::DisjointBoxLayout> vecgrids;
  pout() << "making grids" << endl;
  GeometryService<2>::generateGrids(vecgrids, domain.domainBox(), maxGrid);

  Chombo4::DisjointBoxLayout gridsFine = vecgrids[0];
  Chombo4::DisjointBoxLayout gridsMedi = vecgrids[1];
  Chombo4::DisjointBoxLayout gridsCoar = vecgrids[2];
  IntVect dataGhostIV =   6*IntVect::Unit;
  Point   dataGhostPt = ProtoCh::getPoint(dataGhostIV); 
  int geomGhost = 6;
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

  shared_ptr<graph_distrib_t> graphsFine = geoserv->getGraphs(domFine);
  shared_ptr<graph_distrib_t> graphsMedi = geoserv->getGraphs(domMedi);
  shared_ptr<graph_distrib_t> graphsCoar = geoserv->getGraphs(domCoar);

  
  shared_ptr<EBDictionary<HOEB_MAX_ORDER, Real, CELL, CELL> >  dictionary
    (new     EBDictionary<HOEB_MAX_ORDER, Real, CELL, CELL>
     (geoserv, vecgrids, vecdomain, vecdx, dataGhostPt));

  EBLevelBoxData<CELL,   1>  errMedi(gridsMedi, dataGhostIV, graphsMedi);
  EBLevelBoxData<CELL,   1>  errCoar(gridsCoar, dataGhostIV, graphsCoar);
  

  getICError(errMedi, 
             graphsFine, gridsFine, domFine, dxFine,
             graphsMedi, gridsMedi, domMedi, dxMedi,
             dictionary, geoserv);

  getICError(errCoar, 
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
  pout() << "writing to file" << endl;
  errCoar.writeToFileHDF5(string("IC.errCoar.hdf5"), 0.0);
  errMedi.writeToFileHDF5(string("IC.errMedi.hdf5"), 0.0);
  
  pout() << "||IC errMedi||_max = " << normMedi << std::endl;
  pout() << "||IC errCoar||_max = " << normCoar << std::endl;
  pout() << "Richardson truncation error order for IC= " << order << std::endl;

  return 0;
}

/****/

/****/
template <CENTERING cent>
void
getDevendranFlux(EBLevelBoxData<cent,   1>                                           &  a_deviflux, 
                 const shared_ptr<graph_distrib_t    >                               &  a_graphs,
                 const Chombo4::DisjointBoxLayout                                    &  a_grids,
                 const Chombo4::Box                                                  &  a_dom,
                 const Real                                                          &  a_dx,
                 const shared_ptr< GeometryService<HOEB_MAX_ORDER> >                 &  a_geoserv)
{
  using Chombo4::pout;
  //fill on the host and transfer over to the device rather than go  through the pain and
  //suffering necessary to squeeze this into forall
  IntVect ghost = a_deviflux.ghostVect();
  typedef GraphConstructorFactory<EBHostData<cent, Real, 1> > hostfactoryface_t;
  typedef GraphConstructorFactory<EBHostData<CELL, Real, 1> > hostfactorycell_t;
  typedef CH4_Data_Choreography::DistributedData<EBHostData<CELL, Real, 1> > cell_distrib_t;
  typedef CH4_Data_Choreography::DistributedData<EBHostData<cent, Real, 1> > cent_distrib_t;
  cent_distrib_t  hostflux(a_grids, 1, ghost, hostfactoryface_t(a_graphs));
  int nghost = 6;
  IntVect dataGhostIV = nghost*IntVect::Unit;
  EBLevelBoxData<CELL,   1>  deviphi(a_grids, nghost*IntVect::Unit, a_graphs);
  hoeb::fillPhi(deviphi, a_graphs, a_grids, a_dom, a_dx, a_geoserv);
  
  cell_distrib_t   hostphi(a_grids, 1, ghost, hostfactorycell_t(a_graphs));
  EBLevelBoxData<CELL,   1>::copyToHost(hostphi, deviphi);
  

  int facedir = 0;
  if(cent == XFACE)
  {
    facedir = 0;
  }
  else if(cent == YFACE)
  {
    facedir = 1;
  }
  else if(cent == ZFACE)
  {
    facedir = 2;
  }
  else
  {
    PROTO_ASSERT(false, "bad face type, no donut");
  }
  string dombcname[2*DIM];
  string ebbcname("Dirichlet");
  for(int ivec = 0; ivec < 2*DIM; ivec++)
  {
    dombcname[ivec] = string("Dirichlet");
  }
  
  Chombo4::DataIterator dit = a_grids.dataIterator();
  for(unsigned int ibox = 0; ibox < dit.size(); ibox++)
  {
    hostflux[dit[ibox]].setVal(0.);  //fill all those nasty  covered values with zero
    auto chomboBox =a_grids[dit[ibox]];
    Proto::Box bx = getProtoBox(chomboBox);
    const auto & graph= (*a_graphs)[dit[ibox]];
    for(auto bit = bx.begin(); bit != bx.end(); ++bit)
    {
      auto pt = *bit;
      auto vofs = graph.getVoFs(pt);
      for(int ivof = 0; ivof < vofs.size(); ivof++)
      {
        const auto& vof = vofs[ivof];
        for(SideIterator sit; sit.ok(); ++sit)
        {
          vector<EBIndex<cent> > faces = Proto::getFaces<cent>(vof, sit(), graph);
          for(int iface  = 0; iface < faces.size(); ++iface)
          {
            const auto& face = faces[iface];
            bool dividebyarea = true; //we want the average flux here
            Real fluxpt =
              hoeb::getDevendranFluxFace<cent>(hostphi[dit[ibox]], (*a_graphs)[dit[ibox]],
                                               face, vof,  dombcname, ebbcname, a_geoserv,
                                               a_dom, a_dx, facedir, sit(),
                                               dit[ibox], dividebyarea, false);
            hostflux[dit[ibox]](face, 0) = fluxpt;
          }
        } 
      } 
    }
  }
  EBLevelBoxData<cent, 1>::copyToDevice(a_deviflux, hostflux);
}
/****/
template <CENTERING cent>
void
getDevendranFluxError(EBLevelBoxData<cent,   1>                                           &  a_errCoar, 
                      const shared_ptr<graph_distrib_t>                                   &  a_graphsFine,
                      const Chombo4::DisjointBoxLayout                                    &  a_gridsFine,
                      const Chombo4::Box                                                  &  a_domFine,
                      const Real                                                          &  a_dxFine,
                      const shared_ptr<graph_distrib_t>                               &  a_graphsCoar,
                      const Chombo4::DisjointBoxLayout                                    &  a_gridsCoar,
                      const Chombo4::Box                                                  &  a_domCoar,
                      const Real                                                          &  a_dxCoar,
                      const shared_ptr<EBDictionary<HOEB_MAX_ORDER, Real, cent, cent> >   &  a_dictionary,
                      const shared_ptr< GeometryService<HOEB_MAX_ORDER> >                 &  a_geoserv,
                      bool a_finestLevel)
{
  ParmParse pp;
  int nghost = 6;

  IntVect dataGhostIV =   nghost*IntVect::Unit;
  EBLevelBoxData<cent,   1>  fluxFine(a_gridsFine, dataGhostIV, a_graphsFine);
  EBLevelBoxData<cent,   1>  fluxCoar(a_gridsCoar, dataGhostIV, a_graphsCoar);
  EBLevelBoxData<cent,   1>  fluxFtoC(a_gridsCoar, dataGhostIV, a_graphsCoar);
  
  getDevendranFlux(fluxFine, a_graphsFine, a_gridsFine, a_domFine, a_dxFine, a_geoserv);
  getDevendranFlux(fluxCoar, a_graphsCoar, a_gridsCoar, a_domCoar, a_dxCoar, a_geoserv);

  hoeb::restrictFlux<cent>(fluxFtoC, fluxFine,
                    a_graphsFine, a_gridsFine, a_domFine, a_dxFine,                    
                    a_graphsCoar, a_gridsCoar, a_domCoar, a_dxCoar,
                    a_dictionary, a_geoserv);

  //error = Ave(phifine) - phicoar 
  Chombo4::DataIterator dit = a_gridsCoar.dataIterator();
  for(unsigned int ibox = 0; ibox < dit.size(); ibox++)
  {
    auto& ftocfab =   fluxFtoC[dit[ibox]];
    auto& coarfab =   fluxCoar[dit[ibox]];
    auto& errfab  = a_errCoar[dit[ibox]];
    auto  inputbx = ftocfab.inputBox();
    auto  validbx = (*a_graphsCoar)[dit[ibox]].validBox();
    ebforall(inputbx, subtractionTraction, validbx, errfab, ftocfab, coarfab);
  }

  if(a_finestLevel)
  {
    if(cent == XFACE)
    {
      fluxFine.writeToFileHDF5(string("fluxFine_x.hdf5"), 0.0);
      fluxCoar.writeToFileHDF5(string("fluxCoar_x.hdf5"), 0.0);
      fluxFtoC.writeToFileHDF5(string("fluxFtoC_x.hdf5"), 0.0);
    }
    else if(cent == YFACE)
    {
      fluxFine.writeToFileHDF5(string("fluxFine_y.hdf5"), 0.0);
      fluxCoar.writeToFileHDF5(string("fluxCoar_y.hdf5"), 0.0);
      fluxFtoC.writeToFileHDF5(string("fluxFtoC_y.hdf5"), 0.0);
    }
  }
}
/****/
template<CENTERING cent>
unsigned int
devendranFluxTruncation(string a_prefix)
{
#if 0 
  using Chombo4::pout;
  Real coveredval = -1;
  int nx      = 32;
  int maxGrid = 32;
    
  ParmParse pp;

  pp.get("nx"        , nx);
  pp.get("maxGrid"  , maxGrid);
  pp.get("coveredval", coveredval);         


  pout() << "nx"        << " = " <<  nx         << endl;
  pout() << "max_grid"  << " = " <<  maxGrid    << endl;
  pout() << "coveredval"<< " = " <<  coveredval << endl;         
  

  IntVect domLo = IntVect::Zero;
  IntVect domHi  = (nx - 1)*IntVect::Unit;

  Chombo4::ProblemDomain domain(domLo, domHi);

  vector<Chombo4::DisjointBoxLayout> vecgrids;
  pout() << "making grids" << endl;
  GeometryService<2>::generateGrids(vecgrids, domain.domainBox(), maxGrid);

  Chombo4::DisjointBoxLayout gridsFine = vecgrids[0];
  Chombo4::DisjointBoxLayout gridsMedi = vecgrids[1];
  Chombo4::DisjointBoxLayout gridsCoar = vecgrids[2];
  IntVect dataGhostIV =   6*IntVect::Unit;
  Point   dataGhostPt = ProtoCh::getPoint(dataGhostIV); 
  int geomGhost = 6;
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

  shared_ptr<graph_distrib_t> graphsFine = geoserv->getGraphs(domFine);
  shared_ptr<graph_distrib_t> graphsMedi = geoserv->getGraphs(domMedi);
  shared_ptr<graph_distrib_t> graphsCoar = geoserv->getGraphs(domCoar);

  
  shared_ptr<EBDictionary<HOEB_MAX_ORDER, Real, cent, cent>  >  dictionary
    (new     EBDictionary<HOEB_MAX_ORDER, Real, cent, cent>
     (geoserv, vecgrids, vecdomain, vecdx, dataGhostPt));

  EBLevelBoxData<cent,   1>  errMedi(gridsMedi, dataGhostIV, graphsMedi);
  EBLevelBoxData<cent,   1>  errCoar(gridsCoar, dataGhostIV, graphsCoar);
  

  getDevendranFluxError<cent>(errMedi, 
                              graphsFine, gridsFine, domFine, dxFine,
                              graphsMedi, gridsMedi, domMedi, dxMedi,
                              dictionary, geoserv, true);

  getDevendranFluxError<cent>(errCoar, 
                              graphsMedi, gridsMedi, domMedi, dxMedi,
                              graphsCoar, gridsCoar, domCoar, dxCoar,
                              dictionary, geoserv, false);


  //Norm!
  Real normMedi = errMedi.maxNorm(0);
  Real normCoar = errCoar.maxNorm(0);

  Real tol = 1.0e-16;
  Real order = 0;
  if((normCoar > tol) && (normMedi > tol))
  {
    order = log(normCoar/normMedi)/log(2.0);
  }
  pout() << "writing to file" << endl;
  string fileCoar = a_prefix + string(".errCoar.hdf5");
  string fileMedi = a_prefix + string(".errMedi.hdf5");
  errCoar.writeToFileHDF5(fileCoar, 0.0);
  errMedi.writeToFileHDF5(fileMedi, 0.0);
  
  pout() << "|| " << a_prefix<< "  errMedi ||_max = " << normMedi << std::endl;
  pout() << "|| " << a_prefix<< "  errCoar ||_max = " << normCoar << std::endl;
  pout() << "Richardson truncation error order for Devendran flux = " << order << std::endl;

#endif
  
  return 0;
}
/****/
unsigned int
runTests()
{

  using Chombo4::pout;
  bool runIC, runDCT, runSCT, runDevFluxTE;
  ParmParse pp;

  pp.get("do_initial_condition_test", runIC);
  pp.get("do_devendran_test"        , runDCT);
  pp.get("do_schwartz_test"         , runSCT);
  pp.get("do_devendran_flux_test"   , runDevFluxTE);
  unsigned int errcode = 0;
  if(runDevFluxTE)
  {
    unsigned int errcode1 = devendranFluxTruncation<XFACE>(string("devfluxX")); hoeb::checkError(errcode1, "devendran_flux_truncation_error_x");
    errcode += 10*errcode1;
    unsigned int errcode0 = devendranFluxTruncation<YFACE>(string("devfluxY")); hoeb::checkError(errcode0, "devendran_flux_truncation_error_y");
    errcode += errcode0;
#if DIM==3
    unsigned int errcode2 = devendranFluxTruncation<ZFACE>(string("devfluxZ")); hoeb::checkError(errcode2, "devendran_flux_truncation_error_z");
    errcode += 100*errcode2;
#endif    
  }
  if(runIC)
  {
    unsigned int errcode4 = runInitialConditionTest()     ; hoeb::checkError(errcode4, "initial_condition");
    errcode += 10000*errcode4;
  }
  if(runSCT)
  {
    unsigned int errcode5 = runTruncationErrorTest(string("Schwartz")) ; hoeb::checkError(errcode5, "Schwartz_stencil");
    errcode += 100000*errcode5;
  }
  if(runDCT)
  {
    unsigned int errcode6 = runTruncationErrorTest(string("Devendran")); hoeb::checkError(errcode6, "Devendran_stencil");
    errcode += 1000000*errcode6;
  }
  return errcode;
}
/****/

int main(int a_argc, char* a_argv[])
{
#ifdef CH_USE_PETSC  
  //because of some kind of solipsistic madness, PetscInitialize calls MPI_INIT
  PetscInt ierr = PetscInitialize(&a_argc, &a_argv, "./.petscrc",PETSC_NULL); CHKERRQ(ierr);
#else
  using Chombo4::pout;
#ifdef CH_MPI
  MPI_Init(&a_argc, &a_argv);
  pout() << "MPI INIT called" << std::endl;
#endif
#endif
  using Chombo4::pout;
  //needs to be called after MPI_Init
  CH_TIMER_SETFILE("trunc.time.table");
  {
    printNeighborhood(NULL);
    printCompositeStencil(NULL);
    
    if (a_argc < 2)
    {
      cerr<< " usage " << a_argv[0] << " <input_file_name> " << endl;
      exit(0);
    }
    char* in_file = a_argv[1];
    ParmParse  pp(a_argc-2,a_argv+2,NULL,in_file);
    runTests();
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
