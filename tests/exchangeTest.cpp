#include <cmath>
#include <cstdio>
#include <iostream>
#include <iomanip>
#include <gtest/gtest.h>
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

#define PI 3.1419
using Proto::Var;
/****/

PROTO_KERNEL_START 
unsigned int  setPhiPtF(int              a_pt[DIM],
                        Var<Real, 2>     a_phi,  Real a_dx)
{
  Real x = Real(a_pt[0])*(a_dx + 0.5) + Real(a_pt[1])*(a_dx + 0.8);
  Real y = Real(a_pt[0])*(a_dx + 0.5) - Real(a_pt[1])*(a_dx + 0.8);
  a_phi(0) = x;
  a_phi(1) = y;
  return 0;
}
PROTO_KERNEL_END(setPhiPtF, setPhiPt)

PROTO_KERNEL_START 
unsigned int  verifyPhiPtF(int              a_pt[DIM],
                           Var<Real, 2>     a_phi,  Real a_dx)
{
  Real x = Real(a_pt[0])*(a_dx + 0.5) + Real(a_pt[1])*(a_dx + 0.8);
  Real y = Real(a_pt[0])*(a_dx + 0.5) - Real(a_pt[1])*(a_dx + 0.8);
  Real error_x = a_phi(0) - x;
  Real error_y = a_phi(1) - y;
  if(error_x != 0 || error_y != 0)
    {
      PR_error("Incorrect Exchange Value");
    }
  return 0;
}
PROTO_KERNEL_END(verifyPhiPtF, verifyPhiPt)


Bx convert(const Chombo4::Box& b)
{
  return Bx(b.smallEnd().dataPtr(), b.bigEnd().dataPtr());
}

int runTest()
{
  typedef EBStencil<2, Real, CELL, CELL> ebstencil_t;
  int nx      = 64;
  int maxGrid = 32;

    
  Real R = 0.1;
  RealVect ABC = RealVect::Unit();
  RealVect  X0 = RealVect::Unit();
  X0 *= 0.5;
  IntVect domLo = IntVect::Zero;
  IntVect domHi  = (nx - 1)*IntVect::Unit;
  
// EB and periodic do not mix
  Chombo4::ProblemDomain domain(domLo, domHi);

  Vector<Chombo4::DisjointBoxLayout> vecgrids;
  pout() << "making grids" << endl;
  GeometryService<2>::generateGrids(vecgrids, domain.domainBox(), maxGrid);

  Chombo4::DisjointBoxLayout grids = vecgrids[0];
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


  Chombo4::Box dombox = domain.domainBox();
  shared_ptr<LevelData<EBGraph> > graphs = geoserv->getGraphs(dombox);

  pout() << "making data" << endl;
  EBLevelBoxData<CELL,   2>  phi(grids, dataGhostIV, graphs);

  Chombo4::DataIterator dit = grids.dataIterator();
  pout() << "setting values" << endl;
  for(int ibox = 0; ibox < dit.size(); ibox++)
  {
    EBBoxData<CELL, Real, 2>& phibd = phi[dit[ibox]];
    phibd.setVal(0.0);
    Chombo4::Box validCells = grids[dit[ibox]];
    size_t numflopspt = 0;
    ebforallInPlace_i(numflopspt, "setPhiPt", setPhiPt, convert(validCells), phibd, dx);
  }

  Chombo4::Copier copier(grids, grids, dataGhostIV, true);
  phi.exchange(copier);

  for(int ibox = 0; ibox < dit.size(); ibox++)
    {
      EBBoxData<CELL, Real, 2>& phibd = phi[dit[ibox]];
      Bx d = phibd.box() & convert(dombox);
      size_t flops=0;
      ebforallInPlace_i(flops, "verifyPhiPt", verifyPhiPt, d, phibd, dx);


    }

  
 // string fileq("phi.hdf5");
 // writeEBLevelHDF5<1>(  fileq,  phi, kappa, domain.domainBox(), graphs, coveredval, dx, dt, time);
  
  return 0;
}

TEST(exchangeTest,runTest) {
    int retval = runTest();
    EXPECT_EQ(retval,0);
}

int main(int a_argc, char* a_argv[])
{
  ::testing::InitGoogleTest(&a_argc, a_argv);
#ifdef CH_MPI
  MPI_Init(&a_argc, &a_argv);
  pout() << "MPI INIT called" << std::endl;
#endif
  //needs to be called after MPI_Init
  CH_TIMER_SETFILE("ebapply.time.table");

  int result = RUN_ALL_TESTS();

  pout() << "printing time table " << endl;
  CH_TIMER_REPORT();
  pout() << "exchangeTest PASSED " << endl;
#ifdef CH_MPI
  pout() << "about to call MPI Finalize" << std::endl;
  MPI_Finalize();
#endif
  return result;
}
