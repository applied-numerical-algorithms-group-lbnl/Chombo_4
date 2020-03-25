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
#include "SetupFunctions.H"
#include <iomanip>

#include "Chombo_NamespaceHeader.H"

#define MAX_ORDER 2

using std::cout;
using std::endl;
using std::shared_ptr;


int
runTest(int a_argc, char* a_argv[])
{
  int nx          = 32;
  ParmParse pp;

  pout() << "defining geometry" << endl;
  shared_ptr<GeometryService<MAX_ORDER> >  geoserv;

  pp.get("nx"        , nx);
  pout() << "nx       = " << nx     << endl;
  Real dx = 1.0/Real(nx);

  Vector<DisjointBoxLayout> vecgrids;
  Vector<Box>               vecdomains;
  Vector<Real> vecdx;
  int whichGeom;

  Real geomCen, geomRad;
  defineGeometry(vecgrids, vecdomains, vecdx, geoserv, geomCen, geomRad, whichGeom, dx, nx);

  auto graphs = m_geoserv->getGraphs(m_domain);

  IntVect dataGhostIV =   4*IntVect::Unit;
  Point   dataGhostPt = ProtoCh::getPoint(dataGhostIV);
  auto grids = vecgrids[0];
  EBLevelBoxData<CELL, 1>  cellDat(grids, dataGhostIV, graphs);
  EBLevelFluxData<1>       fluxDat(grids, dataGhostIV, graphs);
  
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
