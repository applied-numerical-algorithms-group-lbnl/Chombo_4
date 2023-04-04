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
#include "MHDsteps.H"
#include "Chombo_LevelBoxData.H"
#include "Chombo_LevelData.H"
#include "Chombo_BaseFab.H"
#include "Chombo_ParmParse.H"
#include "Chombo_LoadBalance.H"
#include "Chombo_ProtoInterface.H"
#include "Chombo_BRMeshRefine.H"
#include "Chombo_AMRIO.H"
#include <string>
#include <iostream>
#include <sstream>

#include "Chombo_Box.H"
#include "Proto_RHScalc.H"
#include "Proto_EulerStep.H"

#define PI 3.141592653589793


typedef ::Proto::Var<Real,DIM> V;
typedef ::Proto::Var<Real,NUMCOMPS> State;
typedef ::Proto::Box Bx;
using   ::Proto::Point;
using   ::Proto::BoxData;
using   ::Proto::Stencil;
using   ::Proto::RK4;
using   ::Proto::RHScalc;
using   ::Proto::EulerStep;
using   ProtoCh::getPoint;
using   ProtoCh::getProtoBox;
using   ProtoCh::getIntVect;
using   ProtoCh::getBox;
using     std::cout;
using     std::endl;
using     std::shared_ptr;

string convertInt(int number)
{
  std::stringstream ss;//create a stringstream
  ss << number;//add number to the stream
  return ss.str();//return a string with the contents of the stream
}

class RunParams
{
public:
  RunParams()
  {
    gamma    = 1.6666666666666666666667;
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

    pout() << "MHD with periodic bcs." << endl;

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
                             const V& a_x,
							 const RunParams& a_params)
{    
    
    double gamma = a_params.gamma;
    double rho = 0.0;
    double p = 0.0;
	double u = 0.0;
	double v = 0.0;
	double w = 0.0;
	double Bx = 0.0;
	double By = 0.0;
	double Bz = 0.0;
	
	
	// //////Modifying parameters for 2D current sheet problem/////
	rho = 1.0;
    p = 0.1;
	u = 0.1 * sin(2*PI*a_x(1));
	if (a_x(0) >= 0.0 && a_x(0) < 0.5){By = 1.0;}
	if (a_x(0) >= 0.5 && a_x(0) < 1.5){By = -1.0;}
	if (a_x(0) >= 1.5 && a_x(0) <= 2.0){By = 1.0;}
	
	// rho = 1.0;
    // p = 0.1;
	// v = 0.1 * sin(2.0*PI*a_x(0));
	// if (a_x(1) >= 0.0 && a_x(1) < 0.5){Bx = 1.0;}
	// if (a_x(1) >= 0.5 && a_x(1) < 1.5){Bx = -1.0;}
	// if (a_x(1) >= 1.5 && a_x(1) <= 2.0){Bx = 1.0;}
	// ////////////////////////////////////////////////////////////
	
	
	//////Modifying parameters for 2D Orszag Tang problem///////
	// Case 1:
	// rho = gamma*((2.0 * gamma)/8.0/PI)*1.0;
    // p = (2.0 * gamma)/8.0/M_PI;
	// u = -sin(2.0*PI*a_x(1));
	// v =  sin(2.0*PI*a_x(0));
	// Bx = -sin(2.0*PI*a_x(1));
	// By =  sin(4.0*PI*a_x(0));
	// Case 2:
	// rho = 1.0;
    // p = 1.0/gamma;
	// u = -sin(2.0*PI*a_x(1));
	// v =  sin(2.0*PI*a_x(0));
	// Bx = -sin(2.0*PI*a_x(1))/gamma;
	// By =  sin(4.0*PI*a_x(0))/gamma;
	////////////////////////////////////////////////////////////
	
	
	//////Modifying parameters for Alfven wave problem///////
	// rho = 1.0;
    // p = 1.0;
	// u =  sin(2.0*PI*a_x(0));
	// Bx = sin(2.0*PI*a_x(0));
	////////////////////////////////////////////////////////////	
	
	double e = p/(gamma-1.0) + rho*(u*u+v*v+w*w)/2.0 + (Bx*Bx+By*By+Bz*Bz)/8.0/PI;
	
#if DIM == 1	
	a_U(0) = rho;      //rho
	a_U(1) = rho*u;    //Momentum-x
	a_U(2) = e;        //Energy
	a_U(3) = Bx;	   //Bx
#endif	
#if DIM == 2
	a_U(0) = rho;      //rho
	a_U(1) = rho*u;    //Momentum-x
	a_U(2) = rho*v;    //Momentum-y 
	a_U(3) = e;        //Energy
	a_U(4) = Bx;	   //Bx
	a_U(5) = By;       //By
#endif
#if DIM == 3
	a_U(0) = rho;      //rho
	a_U(1) = rho*u;    //Momentum-x
	a_U(2) = rho*v;    //Momentum-y 
	a_U(3) = rho*w;    //Momentum-z 
	a_U(4) = e;        //Energy
	a_U(5) = Bx;	   //Bx
	a_U(6) = By;       //By
	a_U(7) = Bz;       //Bz
#endif
	
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
writeData(int step, LevelBoxData<NUMCOMPS> & a_U, Real a_time, Real a_dt, Real a_dx, Real a_nx)
{
#ifdef CH_USE_HDF5
  //string filename = string("state.step") + convertInt(step) + "." + convertInt(a_nx) + "." + convertInt(DIM) + string("d.hdf5");
  string filename = string("state.step") + convertInt(step) + "." + convertInt(DIM) + string("d.hdf5");
  Vector<DisjointBoxLayout> vectGrids(1,a_U.getBoxes());
  const Box domain = vectGrids[0].getDomain();
  Vector<LevelData<FArrayBox>* > vectData(1,NULL);
  IntVect ghostVect = NGHOST*IntVect::Unit;
  LevelData<FArrayBox>* level_data = new LevelData<FArrayBox>(vectGrids[0], NUMCOMPS, ghostVect);
  //write to host here
  LevelBoxData<NUMCOMPS>::copyToHost(*level_data, a_U);

  // single level
  Vector<int> refRatio(1,2); 
  int numLevels = 1;

  Vector<string> vectNames(NUMCOMPS);
  vectNames[0] = "rho";
  vectNames[1] = "momentum_x";
  vectNames[2] = "momentum_y";
  vectNames[3] = "energy";
  vectNames[4] = "B_x";
  vectNames[5] = "B_y";

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
void 
writeData_debug(int step, LevelBoxData<NUMCOMPS> & a_U, Real a_time, Real a_dt, Real a_dx, Real a_nx)
{
#ifdef CH_USE_HDF5
  //string filename = string("debug_state.step") + convertInt(step) + "." + convertInt(a_nx) + "." + convertInt(DIM) + string("d.hdf5");
  string filename = string("debug_state.step") + convertInt(step) + "." + convertInt(DIM) + string("d.hdf5");
  Vector<DisjointBoxLayout> vectGrids(1,a_U.getBoxes());
  const Box domain = vectGrids[0].getDomain();
  Vector<LevelData<FArrayBox>* > vectData(1,NULL);
  IntVect ghostVect = NGHOST*IntVect::Unit;
  LevelData<FArrayBox>* level_data = new LevelData<FArrayBox>(vectGrids[0], NUMCOMPS, ghostVect);
  //write to host here
  LevelBoxData<NUMCOMPS>::copyToHost(*level_data, a_U);

  // single level
  Vector<int> refRatio(1,2); 
  int numLevels = 1;

  Vector<string> vectNames(NUMCOMPS);
  vectNames[0] = "rho";
  vectNames[1] = "momentum_x";
  vectNames[2] = "momentum_y";
  vectNames[3] = "energy";
  vectNames[4] = "B_x";
  vectNames[5] = "B_y";

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
void MHDRun(const RunParams& a_params)
{
  CH_TIME("MHDRun");
  Real tstop = a_params.tmax;
  int maxStep  = a_params.nstepmax;
  int nGhost = NGHOST;

  IntVect domLo = IntVect::Zero;
  IntVect domHi  = (a_params.nx - 1)*IntVect::Unit;
  constexpr bool is_periodic[] = {true, true, true};

  ProblemDomain domain(domLo, domHi, is_periodic);

  Vector<Box> boxes;
  unsigned int blockfactor = 8;
  domainSplit(domain, boxes, a_params.maxgrid, blockfactor);
  Vector<int> procs;
  LoadBalance(procs, boxes);
  DisjointBoxLayout grids(boxes, procs, domain);
//  LevelData<FArrayBox> fabdata(grids, NUMCOMPS, 4*IntVect::Unit);
//  fabdata.exchange();

  //this sneaks in parameters not covered by the RK4 template class
  //not the most elegant solution but it works for single level
  MHDOp::s_dx    = a_params.dx;
  MHDOp::s_gamma = a_params.gamma;

  shared_ptr<LevelBoxData<NUMCOMPS> > Uptr(new LevelBoxData<NUMCOMPS>(grids, nGhost*IntVect::Unit));
  LevelBoxData<NUMCOMPS> &  U = *Uptr;

  bool debug_data = false; // This turns on printing of debug data as well. What data need to be printed is set in a_Rhs in MHDOp::proto_step_test
  bool use_forced_dt = false; // If, true, program will use dt provided in inputs file
  
  shared_ptr<LevelBoxData<NUMCOMPS> > rhsptr(new LevelBoxData<NUMCOMPS>(grids, nGhost*IntVect::Unit));
  LevelBoxData<NUMCOMPS> &  RHS = *rhsptr;
  MHDState state_rhs(rhsptr);
  RHScalc<MHDState, MHDrhsOp, MHDDX> rhscalc;

  
  MHDState state(Uptr);  
  RK4<MHDState, MHDRK4Op, MHDDX> rk4;    
//  EulerStep<MHDState, MHDRK4Op, MHDDX> rk4;    
  EulerStep<MHDState, MHDEulerOp, MHDDX> eulerstep;
  EulerStep<MHDState, MHDViscosityOp, MHDDX> viscositystep;
  
  
  Stencil<Real> Lap2nd = Stencil<Real>::Laplacian();

  pout() << "before initializestate"<< endl;

  DataIterator dit = grids.dataIterator();
#pragma omp parallel for
  for(int ibox = 0; ibox < dit.size(); ibox++)
  {

    Box grid = grids[dit[ibox]];
    Bx valid = getProtoBox(grid);
    Bx grnbx = U[dit[ibox]].box();
    BoxData<Real, DIM> x(grnbx);
    forallInPlace_p(iotaFunc, grnbx, x, MHDOp::s_dx);
    BoxData<Real, NUMCOMPS>& ubd = U[dit[ibox] ];
    
    forallInPlace(InitializeState,grnbx,ubd,x,a_params);

    //smooth initial data
    BoxData<Real, NUMCOMPS> lapu = Lap2nd(ubd, 1.0/24.0);
    ubd += lapu;
  }


  Real maxwave = MHDOp::maxWave(*state.m_U);
  Real dt;  
  if (use_forced_dt == false){
	dt = .25*a_params.cfl*a_params.dx/maxwave;    
  } else {  
	dt = a_params.dt; 
  }
  pout() << "initial maximum wave speed = " << maxwave << ", dt = "<< dt << endl;

  pout() << "after initializestate"<< endl;
  U.exchange(state.m_exchangeCopier);

  Real time = 0.;
  if (debug_data){
    rhscalc.calc(state_rhs,state);
  }
  pout() << "starting time loop"<< endl;
  if(a_params.outinterv > 0)
  {
    writeData(0, U,time,dt,a_params.dx, a_params.nx);
	if (debug_data){
      writeData_debug(0, RHS,time,dt,a_params.dx, a_params.nx);
	}
  }
  for (int k = 1;(k <= maxStep) && (time < tstop);k++)
  {
    {
      CH_TIME("rk4_advance");
      rk4.advance(time,dt,state);
    }
    //this was computed during the advance.
    //so the standard trick is to reuse it.
    maxwave = state.m_velSave;
	// Take step for artificial viscosity
    viscositystep.advance(time,dt,state);
	// Take step for divB term
    eulerstep.advance(time,dt,state);
	
	
    time += dt;
	if (debug_data){
      rhscalc.calc(state_rhs,state);
	}
	if (use_forced_dt == false){
		Real dtnew = a_params.cfl*a_params.dx/maxwave; Real dtold = dt; 
		dt = std::min(1.1*dtold, dtnew);                                
    }
    pout() <<"nstep = " << k << " time = " << time << ", dt = " << dt << endl;
    if((a_params.outinterv > 0) && (k%a_params.outinterv == 0))
    {
      writeData(k, U,time,dt,a_params.dx, a_params.nx);
	  if (debug_data){
        writeData_debug(k, RHS,time,dt,a_params.dx, a_params.nx);
	  }
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
  CH_TIMER_SETFILE("MHD.time.table");

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
    //std::cout << " running MHD exmaple " << std::endl;
    MHDRun(params);
  }
  CH_TIMER_REPORT();
#ifdef CH_MPI
  //std::cout << "about to write timer report " << std::endl;  
  MPI_Finalize();
#endif
  return 0;
}
