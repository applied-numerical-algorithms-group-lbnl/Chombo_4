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
#include "Chombo_GeometryService.H"
#include "Chombo_EBDictionary.H"
#include "Chombo_EBChombo.H"
#include "Chombo_EBLevelBoxData.H"
#include "Chombo_EBLevelFluxData.H"
#include "DebugFunctions.H"
#include "HostDebugFunctions.H"
#include "EBMultigridFunctions.H"

#include <iomanip>

#include "Chombo_NamespaceHeader.H"

#define MAX_ORDER 2

using std::cout;
using std::endl;
using std::shared_ptr;


int
runTest(int a_argc, char* a_argv[])
{
  int nx      = 32;
  int maxGrid = 32;
  Real R = 0.25;
  Real C = 0.5;

    
  ParmParse pp;
  pp.get("geom_rad", R);
  pp.get("geom_cen", C);
  pp.get("nx"      , nx);
  pp.get("maxGrid"   , maxGrid);


  pout() << "nx        = " << nx       << endl;
  pout() << "maxGrid   = " << maxGrid  << endl;
  pout() << "R         = " << R        << endl;
  pout() << "C         = " << C        << endl;

  RealVect ABC = RealVect::Unit();
  RealVect X0  = RealVect::Unit();
  X0 *= C;
  IntVect domLo = IntVect::Zero;
  IntVect domHi  = (nx - 1)*IntVect::Unit;

  Box domain(domLo, domHi);
  Vector<Box> boxes(1, domain);
  Vector<int> procs(1, 0);
  DisjointBoxLayout grids(boxes, procs);
  IntVect ivghost =   2*IntVect::Unit;
  Point   ptghost = ProtoCh::getPoint(ivghost);

  int geomGhost = 4;
  RealVect origin = RealVect::Zero();
  Real dx = 1.0/nx;

  bool insideRegular = false;
  pp.get("inside_regular", insideRegular);
                          
  shared_ptr<BaseIF>    impfunc(new Proto::SimpleEllipsoidIF(ABC, X0, R, insideRegular));
  //shared_ptr<BaseIF>    impfunc(new Proto::AllRegularIF());

  pout() << "defining geometry" << endl;
  GeometryService<2>* geomptr = new GeometryService<2>(impfunc, origin, dx, domain, grids, geomGhost);
  shared_ptr< GeometryService<2> >  geoserv(geomptr);

  pout() << "registering stencil" << std::endl;
  string stenname = StencilNames::Poisson2;
  string dombcname = StencilNames::Neumann;
  string  ebbcname = StencilNames::Dirichlet;
  shared_ptr<EBDictionary<2, Real, CELL, CELL> > 
    dictionary(new EBDictionary<2, Real, CELL, CELL>(geoserv, grids, domain, dx, ptghost));
  dictionary->registerStencil(stenname, dombcname, ebbcname, domain, domain, true, ptghost);
  pout() << "done with registering stencil" << endl;
  auto graphs = geoserv->getGraphs(domain);
  pout() << "defining data" << endl;
  EBLevelBoxData<CELL, 1>  kappaDev(grids, ivghost, graphs);
  pout() << "done with data definition"  << endl;
  DataIterator dit = grids.dataIterator();
  for(int ibox = 0; ibox < dit.size(); ++ibox)
  {
    auto& kappdat = kappaDev[dit[ibox]];
    auto grid = grids[dit[ibox]];
    Bx  grbx = kappdat.inputBox();
    const EBGraph  & graph = (*graphs)[dit[ibox]];
    EBHostData<CELL, Real, 1> hostdat(grbx, graph);
    //fill kappa on the host then copy to the device
    pout() << "calling geoserv::fillKappa"  << endl;
    geoserv->fillKappa(hostdat, grid, dit[ibox], domain);
    // now copy to the device
    pout() << "calling geoserv::copyToDevice"  << endl;
    EBLevelBoxData<CELL, 1>::copyToDevice(kappdat, hostdat);

    pout() << "calling getEBStencil"  << endl;
    auto stencil  = dictionary->getEBStencil(stenname, ebbcname, domain, domain, ibox);
    pout() << "calling getDiagonalWeights"  << endl;
    auto  diagptr  = stencil->getDiagonalWeights();
    auto& diagGhost = kappaDev[dit[ibox]];
    auto& stendiag = *diagptr;
    Bx gridbx = ProtoCh::getProtoBox(grid);
    Bx inputBox = diagGhost.inputBox();

    pout() << "calling ebforall"  << endl;
    ebforall(inputBox, copyDiag,  gridbx, diagGhost, stendiag);
    pout() << "out of   ebforall"  << endl;
  }
  pout() << "leaving function"  << endl;
  return 0;
}

#include "Chombo_NamespaceFooter.H"

int main(int a_argc, char* a_argv[])
{
#ifdef CH_MPI
  MPI_Init(&a_argc, &a_argv);
  pout() << "MPI INIT called" << std::endl;
#endif
  //needs to be called after MPI_Init
  CH_TIMER_SETFILE("ebapply.time.table");
  {
    if (a_argc < 2)
    {
      cerr<< " usage " << a_argv[0] << " <input_file_name> " << endl;
      exit(0);
    }
    char* in_file = a_argv[1];
    ParmParse  pp(a_argc-2,a_argv+2,NULL,in_file);
    Chombo4::runTest(a_argc, a_argv);
    pout() << "out of run test"  << endl;
  }

  pout() << "printing time table " << endl;
  CH_TIMER_REPORT();
#ifdef CH_MPI
  pout() << "about to call MPI Finalize" << std::endl;
  MPI_Finalize();
#endif
  return 0;
}
