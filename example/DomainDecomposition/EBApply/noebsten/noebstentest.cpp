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
    
  ParmParse pp;
  pp.get("nx"      , nx);
  pp.get("maxGrid"   , maxGrid);


  pout() << "nx        = " << nx       << endl;
  pout() << "maxGrid   = " << maxGrid  << endl;


  IntVect domLo = IntVect::Zero;
  IntVect domHi  = (nx - 1)*IntVect::Unit;

  Box domain(domLo, domHi);
  Vector<Box> boxes(1, domain);
  Vector<int> procs(1, 0);
  DisjointBoxLayout grids(boxes, procs);
  IntVect ivghost =   2*IntVect::Unit;
  Point   ptghost = ProtoCh::getPoint(ivghost);

  auto sten = Proto::Stencil<double>::Laplacian();
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
  }

  pout() << "printing time table " << endl;
  CH_TIMER_REPORT();
#ifdef CH_MPI
  pout() << "about to call MPI Finalize" << std::endl;
  MPI_Finalize();
#endif
  return 0;
}
