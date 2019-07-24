#include <cmath>
#include <cstdio>
#include <iostream>


#include "EBProto.H"
#include "EBLevelBoxData.H"
#include "LevelData.H"
#include "BaseFab.H"

#include "ParmParse.H"
#include "LoadBalance.H"
#include "ProtoInterface.H"
#include "BRMeshRefine.H"
#include "GeometryService.H"
#include "EBDictionary.H"
#include "EBChombo.H"
#include <iomanip>

#define MAX_ORDER 2

using std::cout;
using std::endl;
using std::shared_ptr;

typedef Proto::Box Bx;
using   Proto::Point;
using   Proto::BoxData;
using   Proto::Stencil;
using   ProtoCh::getPoint;
using   ProtoCh::getProtoBox;
using   ProtoCh::getIntVect;
using   ProtoCh::getBox;
using     std::cout;
using     std::endl;
using     std::shared_ptr;
using   Proto::BaseIF;
using   Proto::SimpleEllipsoidIF;
using   Proto::CENTERING;
using   Proto::CELL;
int
runTest(int a_argc, char* a_argv[])
{
  int nx      = 32;
  int maxGrid = 32;
  Real x0 = 0.5;
  Real y0 = 0.5;
  Real z0 = 0.5;
  Real A = 1.0;
  Real B = 1.0;
  Real C = 1.0;
  Real R = 0.25;
  int nIter       = 10;
  int nStream    = 8;
  ParmParse pp;
    
  pp.get("nx"     , nx);
  pp.get("niter"  , nIter);
  pp.get("nstream", nStream);
  pp.get("maxGrid", maxGrid);
  pp.get("x0"     , x0);
  pp.get("y0"     , y0);
  pp.get("z0"     , z0);
  pp.get("A"      , A);
  pp.get("B"      , B);
  pp.get("C"      , C);
  pp.get("R"      , R);         

  pout() << "nx      = " << nx       << endl;
  pout() << "maxGrid = " << maxGrid  << endl;
  pout() << "x0      = " << x0       << endl;
  pout() << "y0      = " << y0       << endl;
  pout() << "z0      = " << z0       << endl;
  pout() << "A       = " << A        << endl;
  pout() << "B       = " << B        << endl;
  pout() << "C       = " << C        << endl;
  pout() << "R       = " << R        << endl;

  pout() << "nIter   = " << nIter    << endl;
  pout() << "nstream = " << nStream  << endl;

#ifdef PROTO_CUDA
  Proto::DisjointBoxLayout::setNumStreams(nStream);
#endif

  RealVect ABC, X0;
  ABC[0] = A;
  ABC[1] = B;
  X0[0] = x0;
  X0[1] = y0;
#if DIM==3
  ABC[2] = C;
  X0[2] = z0;
#endif
  IntVect domLo = IntVect::Zero;
  IntVect domHi  = (nx - 1)*IntVect::Unit;
  constexpr bool is_periodic[] = {true, true, true};
  ProblemDomain domain(domLo, domHi, is_periodic);

  Vector<Box> boxes;
  unsigned int blockfactor = 8;
  domainSplit(domain, boxes, maxGrid, blockfactor);
  Vector<int> procs;
  LoadBalance(procs, boxes);
  DisjointBoxLayout grids(boxes, procs, domain);

  IntVect dataGhostIV =   IntVect::Unit;
  Point   dataGhostPt = getPoint(dataGhostIV); 
  int geomGhost = 2;
  RealVect origin = RealVect::Zero();
  Real dx = 1.0/nx;
  shared_ptr<BaseIF>                       impfunc(new SimpleEllipsoidIF(ABC, X0, R, false));
  Bx domainpr = getProtoBox(domain.domainBox());
  shared_ptr<GeometryService<MAX_ORDER> >  geoserv(new GeometryService<MAX_ORDER>(impfunc, origin, dx, domain.domainBox(), grids, geomGhost));
#if 1
  EBDictionary<2, Real, CELL, CELL> dictionary(geoserv, grids, domain.domainBox(), dataGhostPt, dataGhostPt, dx, true);
  typedef EBStencil<2, Real, CELL, CELL> ebstencil_t;
  string stenname("Second_Order_Poisson");
  string dombcname("Periodic");
  string  ebbcname("Neumann");

  dictionary.registerStencil(stenname, dombcname, ebbcname);

  shared_ptr<LevelData<EBGraph> > graphs = geoserv->getGraphs(domain.domainBox());
  EBLevelBoxData<CELL,   1>  srcData(grids, dataGhostIV, graphs);
  EBLevelBoxData<CELL,   1>  dstData(grids, dataGhostIV, graphs);

  DataIterator dit = grids.dataIterator();
#pragma omp parallel for
  for(int ibox = 0; ibox < dit.size(); ibox++)
  {
    EBBoxData<CELL, Real, 1>& srcebbd = srcData[dit[ibox]];
    EBBoxData<CELL, Real, 1>& dstebbd = dstData[dit[ibox]];
    srcebbd.setVal(0.0);
    dstebbd.setVal(0.0);
  }

  for(int iiter = 0; iiter < nIter; iiter++)
  {    
#pragma omp parallel for
    for(int ibox = 0; ibox < dit.size(); ibox++)
    {
      shared_ptr<ebstencil_t> stencil = dictionary.getEBStencil(stenname, ebbcname, ibox);
      stencil->apply(dstData[dit[ibox]], srcData[dit[ibox]]);
    }
  }
#endif
  return 0;
}


int main(int a_argc, char* a_argv[])
{
#ifdef CH_MPI
  MPI_Init(&a_argc, &a_argv);
#endif
  //needs to be called after MPI_Init
  CH_TIMER_SETFILE("euler.time.table");
  {
    if (a_argc < 2)
    {
      cerr<< " usage " << a_argv[0] << " <input_file_name> " << endl;
      exit(0);
    }
    char* in_file = a_argv[1];
    ParmParse  pp(a_argc-2,a_argv+2,NULL,in_file);
    runTest(a_argc, a_argv);
  }

  CH_TIMER_REPORT();
#ifdef CH_MPI
  MPI_Finalize();
#endif
  return 0;
}
