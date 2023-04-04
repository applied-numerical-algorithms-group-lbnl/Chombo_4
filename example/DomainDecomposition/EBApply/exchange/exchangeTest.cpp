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
  int nx      = 32;
  int maxGrid = 32;
  Real x0 = 0.5;
  Real y0 = 0.5;
  Real z0 = 0.5;
  Real A = 1.0;
  Real B = 1.0;
  Real C = 1.0;
  Real R = 0.25;
    
  ParmParse pp;
    
  pp.get("nx"        , nx);

  pp.get("maxGrid"   , maxGrid);

  pp.get("x0"        , x0);
  pp.get("y0"        , y0);
  pp.get("z0"        , z0);
  pp.get("A"         , A);
  pp.get("B"         , B);
  pp.get("C"         , C);
  pp.get("R"         , R);         

  pout() << "nx        = " << nx       << endl;
  pout() << "maxGrid   = " << maxGrid  << endl;
  pout() << "x0        = " << x0       << endl;
  pout() << "y0        = " << y0       << endl;
  pout() << "z0        = " << z0       << endl;
  pout() << "A         = " << A        << endl;
  pout() << "B         = " << B        << endl;
  pout() << "C         = " << C        << endl;
  pout() << "R         = " << R        << endl;
  Real blobCen, blobRad;

  pp.get("blob_cen", blobCen);
  pp.get("blob_rad", blobRad);

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

// EB and periodic do not mix
  ProblemDomain domain(domLo, domHi);

  std::vector<DisjointBoxLayout> vecgrids;
  Chombo4::pout() << "making grids" << endl;
  GeometryService<2>::generateGrids(vecgrids, domain.domainBox(), maxGrid);

  DisjointBoxLayout grids = vecgrids[0];
  grids.printBalance();

  IntVect dataGhostIV =   2*IntVect::Unit;

  int geomGhost = 4;
  RealVect origin = RealVect::Zero();
  Real dx = 1.0/nx;
//  Real dx = 1.0;
  bool insideRegular = false;
  pp.get("inside_regular", insideRegular);
                          
  shared_ptr<BaseIF>    impfunc(new Proto::SimpleEllipsoidIF(ABC, X0, R, insideRegular));
//  Bx domainpr = getProtoBox(domain.domainBox());

  Chombo4::pout() << "defining geometry" << endl;
  GeometryService<2>* geomptr = new GeometryService<2>(impfunc, origin, dx, domain.domainBox(), vecgrids, geomGhost);
  shared_ptr< GeometryService<2> >  geoserv(geomptr);

  auto graphs = geoserv->getGraphs(domain.domainBox());

//  Point   dataGhostPt = ProtoCh::getPoint(dataGhostIV);

  DataIterator dit = grids.dataIterator();
  Copier exchangeCopier;
  exchangeCopier.exchangeDefine(grids, dataGhostIV);
  {
    Chombo4::pout() << "Single-variable  data" << endl;
    EBLevelBoxData<CELL, 1>  cellDat(grids, dataGhostIV, graphs);
    EBLevelFluxData<1>       fluxDat(grids, dataGhostIV, graphs);
    
    unsigned int numflops = 0;
    Chombo4::pout() << "initializing valid data" << endl;
    for(int ibox = 0; ibox < dit.size(); ibox++)
    {
      auto& cellfab = cellDat[dit[ibox]];
      auto& fluxfab = fluxDat[dit[ibox]];

      Bx grid = ProtoCh::getProtoBox(grids[dit[ibox]]);
      ebforallInPlace_i(numflops, "IntializeSpot", InitializeSpot,  grid,
                        cellfab, blobCen, blobRad, dx);
      for(int idir = 0; idir < DIM; idir++)
      {
        auto& xfab = *fluxfab.m_xflux;
        auto& yfab = *fluxfab.m_yflux;

        ebforallInPlace_i(numflops, "IntializeSpot", InitializeSpot,  grid,
                          xfab, blobCen, blobRad, dx);
        ebforallInPlace_i(numflops, "IntializeSpot", InitializeSpot,  grid,
                          yfab, blobCen, blobRad, dx);
      }
    }

    Chombo4::pout() << "calling exchange to fill ghost cells" << endl;
    cellDat.exchange();
    fluxDat.exchange();
    Bx domainbx = ProtoCh::getProtoBox(domain.domainBox());
    for(int ibox = 0; ibox < dit.size(); ibox++)
    {
      auto& cellfab = cellDat[dit[ibox]];
      auto& fluxfab = fluxDat[dit[ibox]];
      Bx grid = cellfab.box();
      grid &= domainbx;
      Chombo4::pout() << "checking cell   data" << endl;
      ebforallInPlace_i(numflops, "checkSpot", checkSpot,  grid,
                        cellfab, blobCen, blobRad, dx);
      for(int idir = 0; idir < DIM; idir++)
      {
        auto& xfab = *fluxfab.m_xflux;
        auto& yfab = *fluxfab.m_yflux;

        Bx xgrid = xfab.box();
        Bx ygrid = yfab.box();
        xgrid &= domainbx;
        ygrid &= domainbx;
      
        Chombo4::pout() << "checking x face data" << endl;
        ebforallInPlace_i(numflops, "checkSpot", checkSpot,  xgrid,
                          xfab, blobCen, blobRad, dx);
        Chombo4::pout() << "checking y face data" << endl;
        ebforallInPlace_i(numflops, "checkSpot", checkSpot,  ygrid,
                          yfab, blobCen, blobRad, dx);
      }
    }
  }
  {
    Chombo4::pout() << "DIM variables  data" << endl;
    EBLevelBoxData<CELL, DIM>  cellDat(grids, dataGhostIV, graphs);
    EBLevelFluxData<DIM>       fluxDat(grids, dataGhostIV, graphs);
    
    unsigned int numflops = 0;
    Chombo4::pout() << "initializing valid data" << endl;
    for(int ibox = 0; ibox < dit.size(); ibox++)
    {
      auto& cellfab = cellDat[dit[ibox]];
      auto& fluxfab = fluxDat[dit[ibox]];

      Bx grid = ProtoCh::getProtoBox(grids[dit[ibox]]);
      ebforallInPlace_i(numflops, "IntializeSpotDIM", InitializeSpotDIM,  grid,
                        cellfab, blobCen, blobRad, dx);
      for(int idir = 0; idir < DIM; idir++)
      {
        auto& xfab = *fluxfab.m_xflux;
        auto& yfab = *fluxfab.m_yflux;

        ebforallInPlace_i(numflops, "IntializeSpotDIM", InitializeSpotDIM,  grid,
                          xfab, blobCen, blobRad, dx);
        ebforallInPlace_i(numflops, "IntializeSpotDIM", InitializeSpotDIM,  grid,
                          yfab, blobCen, blobRad, dx);
      }
    }

    Chombo4::pout() << "calling exchange to fill ghost cells" << endl;
    cellDat.exchange();
    fluxDat.exchange();
    Bx domainbx = ProtoCh::getProtoBox(domain.domainBox());
    for(int ibox = 0; ibox < dit.size(); ibox++)
    {
      auto& cellfab = cellDat[dit[ibox]];
      auto& fluxfab = fluxDat[dit[ibox]];
      Bx grid = cellfab.box();
      grid &= domainbx;
      Chombo4::pout() << "checking cell   data" << endl;
      ebforallInPlace_i(numflops, "checkSpotDIM", checkSpotDIM,  grid,
                        cellfab, blobCen, blobRad, dx);
      for(int idir = 0; idir < DIM; idir++)
      {
        auto& xfab = *fluxfab.m_xflux;
        auto& yfab = *fluxfab.m_yflux;

        Bx xgrid = xfab.box();
        Bx ygrid = yfab.box();
        xgrid &= domainbx;
        ygrid &= domainbx;
      
        Chombo4::pout() << "checking x face data" << endl;
        ebforallInPlace_i(numflops, "checkSpotDIM", checkSpotDIM,  xgrid,
                          xfab, blobCen, blobRad, dx);
        Chombo4::pout() << "checking y face data" << endl;
        ebforallInPlace_i(numflops, "checkSpotDIM", checkSpotDIM,  ygrid,
                          yfab, blobCen, blobRad, dx);
      }
    }
  }
  return 0;
}

#include "Chombo_NamespaceFooter.H"

int main(int a_argc, char* a_argv[])
{
#ifdef CH_MPI
  MPI_Init(&a_argc, &a_argv);
  Chombo4::pout() << "MPI INIT called" << std::endl;
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

  Chombo4::pout() << "printing time table " << endl;
  CH_TIMER_REPORT();
#ifdef CH_MPI
  Chombo4::pout() << "about to call MPI Finalize" << std::endl;
  MPI_Finalize();
#endif
  return 0;
}
