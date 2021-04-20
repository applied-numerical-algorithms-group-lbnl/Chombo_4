#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include <cstdio>
#include <cstring>
#include <cassert>
#include <cmath>

#include <vector>
#include <memory>

#include <iostream>
#include <fstream>
#include <sstream>

#include "Proto.H"
#include "EulerOp.H"
#include "EulerRK4.H"
#include "Chombo_LevelBoxData.H"
#include "Chombo_LevelData.H"
#include "Chombo_BaseFab.H"
#include "Chombo_BFM.H"

#include "Chombo_ParmParse.H"
#include "Chombo_LoadBalance.H"
#include "Chombo_ProtoInterface.H"
#include "Chombo_BRMeshRefine.H"
#include "Chombo_AMRIO.H"
#include <string>
#include <iostream>
#include <sstream>

#include "Chombo_BFM.H"

#include "Chombo_Box.H"
#include "Chombo_EBChombo.H"

#ifdef _OPENMP
#include<omp.h>
#endif

#define PI 3.141592653589793


typedef ::Proto::Var<Real,DIM> V;
typedef ::Proto::Var<Real,NUMCOMPS> State;
typedef ::Proto::Box Bx;
using   ::Proto::Point;
using   ::Proto::BoxData;
using   ::Proto::Stencil;
using   ::Proto::RK4;
using   ::Proto::Reduction;
using   ::ProtoCh::getPoint;
using   ::ProtoCh::getProtoBox;
using   ::ProtoCh::getIntVect;
using   ::ProtoCh::getBox;
using     std::cout;
using     std::endl;
using     std::shared_ptr;



class RunParams
{
public:
  RunParams()
  {
    gamma    = 1.4;
    nstepmax = 0;
    nx       = 64;
    outinterv= 10;
    tmax     = 2.0;
    domsize  = 1.0;
    cfl      = 0.1;
    numstream= 8;
	dt = 0.0;
    resetDx();
  }

  int numstream;
  int nstepmax;
  int maxgrid;
  int nx;
  int outinterv;
  Real tmax;
  Real domsize;
  Real gamma;
  Real dx;
  Real cfl;
  Real dt;
  void resetDx()
  {
    dx = domsize/nx;
  }

  void coarsen(int a_refrat)
  {
    dx *= a_refrat;
    nx /= a_refrat;;
    
  }
  void print() const
  {

    pout() << "Compressible Euler with periodic bcs." << endl;

    pout() << "parameters: "                           << endl;
    pout() << "nx                  =  "   << nx        << endl;
    pout() << "max_grid            =  "   << maxgrid   << endl;
    pout() << "output interval     =  "   << outinterv << endl;
    pout() << "nstepmax            =  "   << nstepmax  << endl;
    pout() << "tmax                =  "   << tmax      << endl;
    pout() << "domain size         =  "   << domsize   << endl;
    pout() << "CFL number          =  "   << cfl       << endl;
    pout() << "gamma               =  "   << gamma     << endl;
    pout() << "num streams         =  "   << numstream << endl;
    pout() << "Selected dt         =  "   << dt        << endl;
  }
};                  
////////////
void
parseInputs(RunParams& a_params)
{
  ParmParse pp;
  pp.get("nx"                 , a_params.nx);
  pp.get("max_grid"           , a_params.maxgrid);
  pp.get("max_step"           , a_params.nstepmax);
  pp.get("output_interval"    , a_params.outinterv);
  pp.get("max_time"           , a_params.tmax);
  pp.get("domain_size"        , a_params.domsize);
  pp.get("cfl"                , a_params.cfl);
  pp.get("num_streams"        , a_params.numstream);
  pp.get("forced_dt"          , a_params.dt);
  a_params.resetDx();
  a_params.print();
}

PROTO_KERNEL_START 
unsigned int InitializeStateF(State& a_U,
                              const V& a_x)
{
  Real gamma = 1.4;
  Real rho0 = gamma;
  Real p0 = 1.;
  Real umag = .5/sqrt(1.*(DIM));
  Real rho = rho0;
  for (int dir = 0; dir < DIM; dir++)
  {
    rho += .1*rho0*sin(2*2*PI*a_x(0));
  }
  Real p = p0 + (rho - rho0)*gamma*p0/rho0;
  a_U(0) = rho;
  Real ke = 0.;
  for (int dir = 1; dir <= DIM; dir++)
  {
    ke += umag*umag;
    a_U(dir) = rho*umag;
  }
  ke *= .5;
  a_U(NUMCOMPS-1) = p/(gamma-1.0) + rho*ke;
  return 0;
}
PROTO_KERNEL_END(InitializeStateF, InitializeState)

//=================================================================================================
PROTO_KERNEL_START
void iotaFuncF(Point           & a_p,
               V               & a_X,
               Real            a_h)
{
  for (int ii = 0; ii < DIM; ii++)
  {
    a_X(ii) = a_p[ii]*a_h + 0.5*a_h;
  }
}
PROTO_KERNEL_END(iotaFuncF,iotaFunc)

/***/
void 
writeData(int step, const Chombo4::LevelBoxData<NUMCOMPS> & a_U, Real a_time, Real a_dt, Real a_dx, Real a_nx)
{
	#ifdef CH_USE_HDF5
  string filename = std::string("state.step") + std::to_string(step) + "." + std::to_string(DIM) + std::string("d.hdf5");

  //a_U.writeToFileHDF5(filename);
  using Chombo4::DisjointBoxLayout;
  using Chombo4::LevelBoxData;
  using Chombo4::DataIterator;
  Vector<DisjointBoxLayout> vectGrids(1,a_U.getBoxes());
  const Chombo4::Box domain = vectGrids[0].getDomain();
  Vector<LevelData<FArrayBox>* > vectData(1,NULL);
  IntVect ghostVect = NGHOST*IntVect::Unit;
  LevelData<FArrayBox>* level_data = new LevelData<FArrayBox>(vectGrids[0], NUMCOMPS, ghostVect);
  //write to host here
  LevelBoxData<NUMCOMPS>::copyToHost(*level_data, a_U);

  // single level
  Vector<int> refRatio(1,2); 
  int numLevels = 1;

  DataIterator dit = vectGrids[0].dataIterator();
  Real maxRho = 0;
  for(dit.begin(); dit.ok(); ++dit)
    {
      Chombo4::Box b = vectGrids[0][dit];
      BFM1<CH_SPACEDIM>(
		level_data->operator[](dit).dataPtr(),
		b, 
		0, 
		1, 
		b, 
		[&maxRho](Real* v){ maxRho = std::max(*v, maxRho);}
	);
    }
  std::cout<<"max Rho "<<maxRho<<"\n";
  
  Vector<string> vectNames(NUMCOMPS);
  vectNames[0] = "rho";
  vectNames[1] = "momentum_x";
  vectNames[2] = "momentum_y";
  vectNames[3] = "energy";
#if CH_SPACEDIM==3
  vectNames[3] = "momentum_z";
  vectNames[4] = "energy";
#endif
  

  vectData[0] = level_data;
  

  WriteAMRHierarchyHDF5(filename,
                      vectGrids,
                      vectData,
                      vectNames,
                      domain,
                      a_dx,
                      a_dt,
                      a_time,
                      refRatio,
                      numLevels);

   delete level_data;
#endif
}
/***/

void eulerRun(const RunParams& a_params)
{


  CH_TIME("eulerRun");
  Real tstop = a_params.tmax;
  int maxStep  = a_params.nstepmax;
  int nGhost = NGHOST;

#ifdef PROTO_CUDA
//  ::Proto::DisjointBoxLayout::setNumStreams(a_params.numstream);
#endif

  IntVect domLo = IntVect::Zero;
  IntVect domHi  = (a_params.nx - 1)*IntVect::Unit;
  constexpr bool is_periodic[] = {true, true, true};

  Chombo4::ProblemDomain domain(domLo, domHi, is_periodic);

  Vector<Chombo4::Box> boxes;
  unsigned int blockfactor = 8;
  domainSplit(domain, boxes, a_params.maxgrid, blockfactor);
  Vector<int> procs;
  LoadBalance(procs, boxes);
  Chombo4::DisjointBoxLayout grids(boxes, procs, domain);
//  LevelData<FArrayBox> fabdata(grids, NUMCOMPS, 4*IntVect::Unit);
//  fabdata.exchange();

  //this sneaks in parameters not covered by the RK4 template class
  //not the most elegant solution but it works for single level
  EulerOp::s_dx    = a_params.dx;
  EulerOp::s_gamma = a_params.gamma;
  using Chombo4::LevelBoxData;
  shared_ptr<LevelBoxData<NUMCOMPS> > Uptr(new LevelBoxData<NUMCOMPS>(grids, nGhost*IntVect::Unit));
  LevelBoxData<NUMCOMPS> &  U = *Uptr;


  EulerState state(Uptr);
  RK4<EulerState,EulerRK4Op,EulerDX> rk4;
  
  pout() << "before initializestate"<< endl;
  using Chombo4::DataIterator;
  DataIterator dit = grids.dataIterator();

  // This constructor can't be called in the next loop
  Stencil<Real> Lap2nd = Stencil<Real>::Laplacian();

  for(int ibox = 0; ibox < dit.size(); ibox++)
  {
 
    //Box grid = grids[dit[ibox]];
    // Bx valid = getProtoBox(grid);
    Bx grnbx = U[dit[ibox]].box();
    BoxData<Real, DIM> x(grnbx);
    forallInPlace_p(iotaFunc, grnbx, x, EulerOp::s_dx);
    BoxData<Real, NUMCOMPS>& ubd = U[dit[ibox] ];
    
    forallInPlace(InitializeState,grnbx,ubd,x);

    //smooth initial data
    BoxData<Real, NUMCOMPS> lapu = Lap2nd(ubd, 1.0/24.0);
    ubd += lapu;
  }


  Reduction<Real> rxn = state.m_Rxn;
  Real maxwave = EulerOp::maxWave(*state.m_U, rxn);
//  Real dt = .25*a_params.cfl*a_params.dx/maxwave;
  Real dt = a_params.dt;
  pout() << "initial maximum wave speed = " << maxwave << ", dt = "<< dt << endl;

  pout() << "after initializestate"<< endl;
  U.exchange(state.m_exchangeCopier);

  // init Time
  Real time = 0.;

  pout() << "starting time loop"<< endl;
  if(a_params.outinterv > 0)
  {
    writeData(0, U,time,dt,a_params.dx, a_params.nx);
  }
  for (int k = 1;(k <= maxStep) && (time < tstop);k++)
  {
    rxn.reset();
    {
      CH_TIME("rk4_advance");
      rk4.advance(time,dt,state);
    }
    //this was computed during the advance.
    //so the standard trick is to reuse it.
    maxwave = gatherMaxWave(rxn.fetch());

    time += dt;
    
  //  Real dtnew = a_params.cfl*a_params.dx/maxwave; Real dtold = dt;
  //  dt = std::min(1.1*dtold, dtnew);

    pout() <<"nstep = " << k << " time = " << time << ", dt = " << dt << endl;
    if((a_params.outinterv > 0) && (k%a_params.outinterv == 0))
    {
      writeData(k, U,time,dt,a_params.dx, a_params.nx);
    }
  }

}
////////////
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

    RunParams params;
    //std::cout << " parsing inputs " << std::endl;
    parseInputs(params);
    //std::cout << " running euler exmaple " << std::endl;
    eulerRun(params);
  }
  CH_TIMER_REPORT();
#ifdef CH_MPI
  //std::cout << "about to write timer report " << std::endl;  
  MPI_Finalize();
#endif
  return 0;
}
