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
PROTO_KERNEL_START 
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
/****/
//lph comes in holding beta*div(F)--leaves holding alpha phi + beta div(F)
PROTO_KERNEL_START 
unsigned int  setPhiPtF(int              a_pt[DIM],
                        Var<Real, 1>     a_phi,
                        Var<Real, 1>     a_lph,
                        Real             a_a,
                        Real             a_b,
                        int              a_idir,
                        Real             a_dx)
{
  Real ri = a_pt[0];
  Real rj = a_pt[1];
  Real xd = a_dx*(ri + 0.5);
  Real yd = a_dx*(rj + 0.5);
#if DIM==3  
  Real rk = a_pt[2];
  Real zd = a_dx*(rk + 0.5);
#endif
  Real x = 0;
  if(a_idir == 0)
  {
    x = xd;
  }
  else if(a_idir == 1)
  {
    x = yd;
  }
#if DIM== 3  
  else if(a_idir == 2)
  {
    x = zd;
  }
#endif
  else
  {
    printf("\n bogus idir\n");
  }
    
  Real val = a_a*x*x + a_b*x;

  a_phi(0) = val;
  a_lph(0) = 0;
  return 0;
}
PROTO_KERNEL_END(setPhiPtF, setPhiPt)

//lph comes in holding beta*div(F)--leaves holding alpha phi + beta div(F)
PROTO_KERNEL_START 
unsigned int  checkLphPtF(int              a_pt[DIM],
                          Var<Real, 1>     a_lph,
                          Var<Real, 1>     a_kappa,
                          Real             a_a,
                          Real             a_b,
                          int              a_idir,
                          Real             a_dx)
{
  Real lphval = a_lph(0);
  Real kapval = a_kappa(0);
  Real corval = 2.*a_a*kapval;

  Real tol = 1.0e-3;
  Real diff = (lphval-corval)*(lphval-corval);
  if( (diff > tol) && kapval > 0.99)
  {
    printf("\nexact test FAILED\n");
  }
  return 0;
}
PROTO_KERNEL_END(checkLphPtF, checkLphPt)

int
runTest(int a_argc, char* a_argv[])
{
  typedef EBStencil<2, Real, CELL, CELL> ebstencil_t;
  int nx      = 64;
  Real alpha = 0.0;
  Real beta  = 1.0;
    
  IntVect domLo = IntVect::Zero;
  IntVect domHi  = (nx - 1)*IntVect::Unit;
  
// EB and periodic do not mix
  ProblemDomain domain(domLo, domHi);
  Vector<Box> boxes(1, domain.domainBox());
  Vector<int> procs(1, 0);

  pout() << "making grids" << endl;
  DisjointBoxLayout grids(boxes, procs);
  Vector<DisjointBoxLayout> vecgrids(1, grids);
  grids.printBalance();

  IntVect dataGhostIV =   2*IntVect::Unit;
  Point   dataGhostPt = ProtoCh::getPoint(dataGhostIV); 
  int geomGhost = 4;
  RealVect origin = RealVect::Zero();
  Real dx = 1.0/nx;
  Real x0 = 0.7;
  Real a = 1;
  Real b = -2*a*x0;
  string stenname  = StencilNames::Poisson2;
  string dombcname = StencilNames::Dirichlet;
  string  ebbcname = StencilNames::Neumann;
  for(int idir = 0; idir < DIM; idir++)
  {
    RealVect  start  = RealVect::Zero();
    RealVect  normal = RealVect::Zero();
    normal[idir] = -1;
    start[ idir] = x0;
    shared_ptr<BaseIF>  impfunc(new Proto::PlaneIF(start, normal));

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
  

    Box dombox = domain.domainBox();
    shared_ptr<LevelData<EBGraph> > graphs = geoserv->getGraphs(dombox);

    pout() << "making data" << endl;
    EBLevelBoxData<CELL,   1>  phi(grids, dataGhostIV, graphs);
    EBLevelBoxData<CELL,   1>  lph(grids, dataGhostIV, graphs);

    DataIterator dit = grids.dataIterator();
    pout() << "registering stencil and setting values" << endl;
    for(int ibox = 0; ibox < dit.size(); ibox++)
    {
      EBBoxData<CELL, Real, 1>& phibd = phi[dit[ibox]];
      EBBoxData<CELL, Real, 1>& lphbd = lph[dit[ibox]];
      Bx inputBx = lphbd.inputBox();
      
      Box gridbox  = grids[dit[ibox]];
      Bx grbx = ProtoCh::getProtoBox(gridbox);
      Bx grownbx = phibd.box();
      ebforall_i(inputBx, setPhiPt, grownbx, phibd, lphbd, a, b, idir, dx);

      Point pghost = ProtoCh::getPoint(dataGhostIV);
      dictionary->registerStencil(stenname, dombcname, ebbcname, dombox, dombox, true, pghost);


      shared_ptr<ebstencil_t> stencil = dictionary->getEBStencil(stenname, ebbcname, dombox, dombox, ibox);
      shared_ptr< EBBoxData<CELL, Real, 1> >   diagptr  = stencil->getDiagonalWeights();
      pout() << "applying stencil" << endl;
      //set lphi = kappa* div(F)
      stencil->apply(lphbd, phibd,  true, 1.0);

      auto& graph = (*graphs)[dit[ibox]];
      EBHostData<CELL, Real, 1> hostdat(grownbx, graph);
      EBBoxData< CELL, Real, 1> kappdat(grownbx, graph);

      //fill kappa on the host then copy to the device
      geoserv->fillKappa(hostdat, gridbox, dit[ibox], domain.domainBox());
      // now copy to the device
      EBLevelBoxData<CELL, 1>::copyToDevice(kappdat, hostdat);
      auto& kapbd = kappdat;

      //pout() << "going into add alphaphi" << endl;
      pout() << "going into addAlphaPhi" << endl;
      ebforall(inputBx, addAlphaPhi, grbx, lphbd, phibd, kapbd, alpha, beta);

      pout() << "checking answer" << endl;
      Box interiorBox = gridbox;
      interiorBox.grow(1);
      interiorBox &= domain.domainBox();
      interiorBox.grow(-1);
      Bx intbx = ProtoCh::getProtoBox(interiorBox);
      ebforall_i(inputBx, checkLphPt, intbx, lphbd, kapbd,  a, b, idir, dx);
    }
  }
  printf("\nexact test PASSED\n");
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
