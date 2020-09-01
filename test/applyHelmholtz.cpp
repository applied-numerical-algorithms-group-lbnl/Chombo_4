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


#include "Chombo_EBLevelFluxData.H"

#include <iomanip>
#define PI 3.1419
using Proto::Var;
/****/
/*PROTO_KERNEL_START 
unsigned int  addAlphaPhiF(Var<Real, 1>     a_lph,
                           Var<Real, 1>     a_phi,
                           Var<Real, 1>     a_kap,
                           Real    a_alpha,
                           Real    a_beta)
{
  a_lph(0) = a_alpha*a_phi(0)*a_kap(0)+ a_beta*a_lph(0);
  return 0;
}
PROTO_KERNEL_END(addAlphaPhiF, addAlphaPhi)
*/


PROTO_KERNEL_START 
unsigned int  addAlphaPhiF(Var<Real, 1>     a_lph,
                           Var<Real, 1>     a_phi,
                           Var<Real, 1>     a_kap,
                           Real    a_alpha,
                           Real    a_beta)
{
  a_lph(0) = a_alpha*a_phi(0)*a_kap(0) + a_beta*a_lph(0);
  return 0;
}
PROTO_KERNEL_END(addAlphaPhiF, addAlphaPhi)

/****/
//lph comes in holding beta*div(F)--leaves holding alpha phi + beta div(F)
PROTO_KERNEL_START 
unsigned int  setPhiPtF(int              a_pt[DIM],
                        Var<Real, 1>     a_phi,
                        Var<Real, 1>     a_lph,
                        Real             a_dx)
{
  Real pi = 4.*atan(1.0);
  Real x = Real(a_pt[0])*(a_dx + 0.5);
  Real y = Real(a_pt[1])*(a_dx + 0.5);
  Real val = sin(pi*x*y);
  a_phi(0) = val;
  a_lph(0) = 0;
  return 0;
}
PROTO_KERNEL_END(setPhiPtF, setPhiPt)

int
runTest(int a_argc, char* a_argv[])
{
  typedef EBStencil<2, Real, CELL, CELL> ebstencil_t;
  int nx      = 128;
  int maxGrid = 64;

  Real alpha = 1.0;
  Real beta = -0.001;
    
  Real R = 0.1;
  RealVect ABC = RealVect::Unit();
  RealVect  X0 = RealVect::Unit();
  X0 *= 0.5;
  IntVect domLo = IntVect::Zero;
  IntVect domHi  = (nx - 1)*IntVect::Unit;
  
// EB and periodic do not mix
  ProblemDomain domain(domLo, domHi);

  Vector<DisjointBoxLayout> vecgrids;
  pout() << "making grids" << endl;
  GeometryService<2>::generateGrids(vecgrids, domain.domainBox(), maxGrid);

  DisjointBoxLayout grids = vecgrids[0];
  grids.printBalance();

  IntVect dataGhostIV =   2*IntVect::Unit;
  Point   dataGhostPt = ProtoCh::getPoint(dataGhostIV); 
  int geomGhost = 4;
  RealVect origin = RealVect::Zero();
  Real dx = 1.0/nx;
  shared_ptr<BaseIF>    impfunc(new Proto::SimpleEllipsoidIF(ABC, X0, R, true));

  pout() << "defining geometry" << endl;
  GeometryService<2>* geomptr = new GeometryService<2>(impfunc, origin, dx, domain.domainBox(), vecgrids, geomGhost);
  shared_ptr< GeometryService<2> >  geoserv(geomptr);

  pout() << "making dictionary" << endl;
  vector<Box>    vecdomain(vecgrids.size(), domain.domainBox());
  vector<Real>   vecdx    (vecgrids.size(), dx);
  for(int ilev = 1; ilev < vecgrids.size(); ilev++)
  {
    vecdomain[ilev] = coarsen(vecdomain[ilev-1], 2);
    vecdx    [ilev] =           2*vecdx[ilev-1];
  }
  shared_ptr<EBDictionary<2, Real, CELL, CELL> > 
    dictionary(new EBDictionary<2, Real, CELL, CELL>(geoserv, vecgrids, vecdomain, vecdx, dataGhostPt));
  
  string stenname  = StencilNames::Poisson2;
  string dombcname = StencilNames::Neumann;
  string  ebbcname = StencilNames::Dirichlet;

  Box dombox = domain.domainBox();
  shared_ptr<LevelData<EBGraph> > graphs = geoserv->getGraphs(dombox);

  pout() << "making data" << endl;
  EBLevelBoxData<CELL,   1>  phi(grids, dataGhostIV, graphs);
  EBLevelBoxData<CELL,   1>  lph(grids, dataGhostIV, graphs);

  DataIterator dit = grids.dataIterator();
  pout() << "setting values" << endl;
  for(int ibox = 0; ibox < dit.size(); ibox++)
  {
    EBBoxData<CELL, Real, 1>& phibd = phi[dit[ibox]];
    EBBoxData<CELL, Real, 1>& lphbd = lph[dit[ibox]];
    Bx grbx = phibd.box();
    size_t numflopspt = 0;
    ebforallInPlace_i(numflopspt, "setPhiPt", setPhiPt, grbx, phibd, lphbd, dx);

  }

    Point pghost = ProtoCh::getPoint(dataGhostIV);
    dictionary->registerStencil(stenname, dombcname, ebbcname, dombox, dombox, true, pghost);

  pout() << "applying operator" << endl;
  for(int ibox = 0; ibox < dit.size(); ibox++)
  {
    auto& phifab =   phi[dit[ibox]];
    auto& lphfab =   lph[dit[ibox]];
    Box gridbox  = grids[dit[ibox]];
    Bx grbx = lphfab.box();//ProtoCh::getProtoBox(gridbox);
    //Bx grbx = ProtoCh::getProtoBox(gridbox);
    std::cout << " grbx lo: " << grbx.low() << " high: " << grbx.high() << std::endl;
    shared_ptr<ebstencil_t> stencil = dictionary->getEBStencil(stenname, ebbcname, dombox, dombox, ibox);
//    shared_ptr< EBBoxData<CELL, Real, 1> >   diagptr  = stencil->getDiagonalWeights();
    
   //set lphi = kappa* div(F)
   stencil->apply(lphfab, phifab,  true, 1.0);


  // auto& kapfab = *diagptr;
    //pout() << "going into add alphaphi" << endl;
    Bx inputBox = lphfab.inputBox();
    size_t numflopspt = 0;
//    ebforall(inputBox, addAlphaPhi, grbx, lphfab, phifab, *diagptr, alpha, beta);
    ebforall(inputBox, addAlphaPhi, grbx, lphfab, phifab, *(stencil->getDiagonalWeights()), alpha, beta);
  }
  return 0;
}


int main(int a_argc, char* a_argv[])
{
#ifdef CH_MPI
  MPI_Init(&a_argc, &a_argv);
  pout() << "MPI INIT called" << std::endl;
#endif
  //needs to be called after MPI_Init
  CH_TIMER_SETFILE("ebapply.time.table");

  runTest(a_argc, a_argv);

  pout() << "printing time table " << endl;
  CH_TIMER_REPORT();
  pout() << "applyHelmholtz PASSED " << endl;
#ifdef CH_MPI
  pout() << "about to call MPI Finalize" << std::endl;
  MPI_Finalize();
#endif
  return 0;
}
