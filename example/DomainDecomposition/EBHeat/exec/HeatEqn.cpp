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
#include "Chombo_EBEncyclopedia.H"
#include "Chombo_EBDictionary.H"
#include "Chombo_EBChombo.H"
#include "EBMultigrid.H"
#include "EBParabolicIntegrators.H"
#include "SetupFunctions.H"
#include "Chombo_EBChombo.H"

#include <iomanip>

#define MAX_ORDER 2

using std::cout;
using std::endl;
using std::shared_ptr;
using Proto::Var;
using Proto::SimpleEllipsoidIF;
using Chombo4::DisjointBoxLayout;
using Chombo4::Box;
using Chombo4::DataIterator;
using Chombo4::LevelBoxData;


//=================================================
void 
initializeData(EBLevelBoxData<CELL,   1>   &  a_scalarCell,
               EBLevelBoxData<CELL,   1>   &  a_sourceTerm,
               const Chombo4::DisjointBoxLayout     &  a_grids,
               const Real                  &  a_dx)

{
using Chombo4::DisjointBoxLayout;
using Chombo4::Box;
using Chombo4::DataIterator;
using Chombo4::LevelBoxData;
  DataIterator dit = a_grids.dataIterator();
  Real blobCen, blobRad;
  ParmParse pp;
  pp.get("blob_cen", blobCen);
  pp.get("blob_rad", blobRad);
  Chombo4::pout() << "blob_cen = " << blobCen       << endl;
  Chombo4::pout() << "blob_rad = " << blobRad       << endl;

  Chombo4::pout() << "calling initializespot for scalar rhs " << endl;

  for(int ibox = 0; ibox < dit.size(); ibox++)
  {
    auto& sourcefab = a_sourceTerm[dit[ibox]];
    auto& scalarfab = a_scalarCell[dit[ibox]];
    unsigned long long int numflopsscal = 5*DIM +3;
    Bx scalbox = scalarfab.box();

    ebforallInPlace_i(numflopsscal, "IntializeSpot", InitializeSpot,  scalbox, 
                      sourcefab, scalarfab, blobCen, blobRad, a_dx);
  }
}
///
shared_ptr<EBMultigrid> 
getEBMultigrid(shared_ptr<EBEncyclopedia<2, Real> > & a_brit,
               shared_ptr<GeometryService<2> >      & a_geoserv,
               const Chombo4::DisjointBoxLayout              & a_grids,
               const Real                           & a_dx,
               const IntVect                        & a_dataGhostIV)

{
using Chombo4::DisjointBoxLayout;
using Chombo4::Box;
using Chombo4::DataIterator;
using Chombo4::LevelBoxData;
  ParmParse pp;
  //these get overwritten for heat solves
  Real alpha = 1.; Real beta  = 1.;
  int numSmooth  = 2;
  string stenname = StencilNames::Poisson2;
  string dombcname, ebbcname;
  Box domain = a_grids.physDomain().domainBox();
  pp.get("numSmooth" , numSmooth);         
  Chombo4::pout() << "num smooths = " << numSmooth;

  int dombc, ebbc;
  pp.get("domainBC"  , dombc);
  pp.get("EBBC"      , ebbc);
  if(dombc == 0)
  {
    dombcname = StencilNames::Neumann;
    Chombo4::pout() << "using Neumann BCs at domain" << endl;
  }
  else
  {
    dombcname = StencilNames::Dirichlet;
    Chombo4::pout() << "using Dirichlet BCs at domain" << endl;
  }

  if(ebbc == 0)
  {
    ebbcname = StencilNames::Neumann;
    Chombo4::pout() << "using Neumann BCs at EB" << endl;
  }
  else
  {
    ebbcname = StencilNames::Dirichlet;
    Chombo4::pout() << "using Dirichlet BCs at EB" << endl;
  }
  auto dictionary = a_brit->m_cellToCell;
  shared_ptr<EBMultigrid> 
    solver(new EBMultigrid (dictionary, a_geoserv, alpha, beta, a_dx, a_grids, stenname, dombcname, ebbcname, domain, a_dataGhostIV));

  EBMultigridLevel::s_numSmoothUp   = numSmooth;
  EBMultigridLevel::s_numSmoothDown = numSmooth;

  return solver;

}
///
int
runHeatEqn(int a_argc, char* a_argv[])
{
using Chombo4::DisjointBoxLayout;
using Chombo4::Box;
using Chombo4::DataIterator;
using Chombo4::LevelBoxData;
  Real coveredval = -1;
  int nx          = 32;
  int  max_step   = 10;
  Real max_time   = 1.0;

  int nStream    = 8;
  int outputInterval = -1;
  ParmParse pp;

  pp.get("nStream", nStream);

    
  pp.get("max_step"  , max_step);
  pp.get("max_time"  , max_time);
  pp.get("output_interval", outputInterval);
  pp.get("covered_value", coveredval);

  Chombo4::pout() << "nStream         = " << nStream         << endl;
  Chombo4::pout() << "max_step        = " << max_step        << endl;
  Chombo4::pout() << "max_time        = " << max_time        << endl;
  Chombo4::pout() << "output interval = " << outputInterval  << endl;

#ifdef PROTO_CUDA
  //Proto::DisjointBoxLayout::setNumStreams(nStream);
#endif



  Chombo4::pout() << "defining geometry" << endl;
  shared_ptr<GeometryService<MAX_ORDER> >  geoserv;

  Real diffCoef;
  pp.get("nx"        , nx);

  Chombo4::pout() << "nx       = " << nx     << endl;
  Real dx = 1.0/Real(nx);

  std::vector<Chombo4::DisjointBoxLayout> vecgrids;
  std::vector<Chombo4::Box>               vecdomains;
  std::vector<Real> vecdx;

  defineGeometry(vecgrids, vecdomains, vecdx, geoserv, dx, nx);

  IntVect dataGhostIV =   4*IntVect::Unit;
  Point   dataGhostPt = ProtoCh::getPoint(dataGhostIV); 
  
  
  Chombo4::pout() << "making dictionary" << endl;
  shared_ptr<EBEncyclopedia<2, Real> > 
    brit(new EBEncyclopedia<2, Real>(geoserv, vecgrids, vecdomains, vecdx, dataGhostPt));


  Chombo4::pout() << "inititializing data"   << endl;
  
  Chombo4::Box domain              = vecdomains[0];
  DisjointBoxLayout grids =   vecgrids[0];
  auto graphs = geoserv->getGraphs(domain);
  EBLevelBoxData<CELL,   1>  sourceTerm(grids, dataGhostIV, graphs);
  EBLevelBoxData<CELL,   1>  scalarCell(grids, dataGhostIV, graphs);
  initializeData(scalarCell, sourceTerm,  grids, dx);


  int step = 0; Real time = 0;
  Real dt = dx;
  Chombo4::pout() << "setting the dt = dx" << endl;


  Chombo4::pout() << "defining Helmholtz solver " << endl;
  shared_ptr<EBMultigrid> ebmg =getEBMultigrid(brit, geoserv, grids, dx, dataGhostIV);
  const EBLevelBoxData<CELL, 1> & kappa = ebmg->getKappa();

  if(outputInterval > 0)
  {
    string filep("scal.0.hdf5");
    writeEBLevelHDF5<1>(  filep,  scalarCell, kappa, domain, graphs, coveredval, dx, dt, time);

    string fileq("sourceTerm.hdf5");
    writeEBLevelHDF5<1>(  fileq,  sourceTerm, kappa, domain, graphs, coveredval, dx, dt, time);
  }

  pp.get("diffusion_coefficient", diffCoef);
  Chombo4::pout() << "Diffusion coeficient = " << diffCoef << endl;

  Real tol = 0.00001;
  int  maxIter = 10;

  pp.get("maxIter"   , maxIter);
  pp.get("tolerance" , tol);
  Chombo4::pout() << "tolerance = " << tol    << endl;
  Chombo4::pout() << "maxIter   = " << maxIter    << endl;

  int whichSolver = 0;
  pp.get("which_solver", whichSolver);
  shared_ptr<BaseEBParabolic> heatIntegrator;
  if(whichSolver == 0)
  {
    Chombo4::pout() << "using backward Euler for time integration" << endl;
    heatIntegrator = shared_ptr<BaseEBParabolic>(new EBBackwardEuler(ebmg, geoserv, grids, domain, dataGhostIV));
  }
  else if(whichSolver == 1)
  {
    Chombo4::pout() << "using Crank Nicolson for time integration" << endl;
    heatIntegrator = shared_ptr<BaseEBParabolic>(new EBCrankNicolson(ebmg, geoserv, grids, domain, dataGhostIV));
  }
  else if(whichSolver == 2)
  {
    Chombo4::pout() << "using TGA for time integration" << endl;
    heatIntegrator = shared_ptr<BaseEBParabolic>(new EBTGA(ebmg, geoserv, grids, domain, dataGhostIV));
  }
  else
  {
    Chombo4::MayDay::Error("bogus whichSolver");
  }

  Chombo4::pout() << "running heat equation operator" << endl;
  while((step < max_step) && (time < max_time))
  {
    heatIntegrator->advanceOneStep(scalarCell, sourceTerm, diffCoef, dt, tol, maxIter);

    Chombo4::pout() <<" step = " << step << " time = " << time << " time step = " << dt << endl;
    step++;
    time += dt;

    if((outputInterval > 0) && ( (step%outputInterval == 0) || step == (max_step-1)))
    {
      string filep = string("scal.") + std::to_string(step) + string(".hdf5");
      writeEBLevelHDF5<1>(  filep,  scalarCell, kappa, domain, graphs, coveredval, dx, dt, time);
    }
  }
  Chombo4::pout() << "exiting runAdvection" << endl;
  return 0;
}


int main(int a_argc, char* a_argv[])
{
#ifdef CH_MPI
  MPI_Init(&a_argc, &a_argv);
  Chombo4::pout() << "MPI INIT called" << std::endl;
#endif
  //needs to be called after MPI_Init
  CH_TIMER_SETFILE("heateqn.time.table");
  {
    if (a_argc < 2)
    {
      cerr<< " usage " << a_argv[0] << " <input_file_name> " << endl;
      exit(0);
    }
    char* in_file = a_argv[1];
    ParmParse  pp(a_argc-2,a_argv+2,NULL,in_file);
    runHeatEqn(a_argc, a_argv);
  }

  Chombo4::pout() << "printing time table " << endl;
  CH_TIMER_REPORT();
#ifdef CH_MPI
  Chombo4::pout() << "about to call MPI Finalize" << std::endl;
  MPI_Finalize();
#endif
  return 0;
}


