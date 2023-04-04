#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

//#define DIM 2


constexpr unsigned int NUMCOMPS=DIM+2;

//#define NUMCOMPS 4

//#define CH_SPACEDIM 2
//#define CH_USE_64 64

//#define CH_USE_DOUBLE TRUE
//#define DCH_USE_HDF5 TRUE

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
unsigned int InitializeStateF(State a_U,
                              V a_x)
{
  for(int icomp = 0; icomp < NUMCOMPS; icomp++)
  {
    Real rightAns = Real(icomp+1)*sin(2.*PI*a_x(0));
    a_U(icomp) = rightAns;
  }
  return 0;
}
PROTO_KERNEL_END(InitializeStateF, InitializeState)


PROTO_KERNEL_START 
unsigned int checkStateF(int a_p[DIM],
                         State a_U,
                         V a_x)
{
  for(int icomp = 0; icomp < NUMCOMPS; icomp++)
  {
    Real rightAns = Real(icomp+1)*sin(2.*PI*a_x(0));
    Real uans = a_U(icomp);
    Real diff = (uans - rightAns)*(uans - rightAns);
    if(diff > 1.0e-5)
    {
      printf("wrong ans at comp %d, point (%d, %d)\n", icomp, a_p[0], a_p[1]);
    }
  }
  return 0;
}
PROTO_KERNEL_END(checkStateF, checkState)

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



void protoRun(const RunParams& a_params)
{

  pout() << "proto version using boxdata" << endl;

  int nGhost = 4;
  using Chombo4::ProblemDomain;
  using Chombo4::DisjointBoxLayout;
  using Chombo4::LevelBoxData;
  using Chombo4::Copier;
  using Chombo4::DataIterator;
  IntVect domLo = IntVect::Zero;
  IntVect domHi  = (a_params.nx - 1)*IntVect::Unit;
  constexpr bool is_periodic[] = {true, true, true};
  cout << "periodic flags: " << is_periodic[0] << ", " << is_periodic[1] << endl;
  ProblemDomain domain(domLo, domHi, is_periodic);

  Vector<Chombo4::Box> boxes;
  unsigned int blockfactor = 8;
  domainSplit(domain, boxes, a_params.maxgrid, blockfactor);
  Vector<int> procs;
  LoadBalance(procs, boxes);
  DisjointBoxLayout grids(boxes, procs, domain);


  shared_ptr<LevelBoxData<NUMCOMPS> > Uptr(new LevelBoxData<NUMCOMPS>(grids, nGhost*IntVect::Unit));
  LevelBoxData<NUMCOMPS> &  U = *Uptr;

  static Copier s_exchangeCopier;
  s_exchangeCopier.exchangeDefine(U.disjointBoxLayout(), U.ghostVect());

  DataIterator dit = grids.dataIterator();


  for(int ibox = 0; ibox < dit.size(); ibox++)
  {
 
    Chombo4::Box grid = grids[dit[ibox]];
    Bx valid = getProtoBox(grid);
    Bx grnbx = U[dit[ibox]].box();
    BoxData<Real, DIM> x(grnbx);
    forallInPlace_p(iotaFunc, grnbx, x, a_params.dx);
    BoxData<Real, NUMCOMPS>& ubd = U[dit[ibox] ];
    
    forallInPlace(InitializeState,valid,ubd,x);
  }
  
  U.exchange(s_exchangeCopier);
  for(int ibox = 0; ibox < dit.size(); ibox++)
  {
 
    Bx grnbx = U[dit[ibox]].box();
    BoxData<Real, DIM> x(grnbx);
    forallInPlace_p(iotaFunc, grnbx, x, a_params.dx);
    BoxData<Real, NUMCOMPS>& ubd = U[dit[ibox] ];
    
    forallInPlace_i(checkState,grnbx,ubd,x);
  }
}
/**/
Real getPerValue(IntVect a_p,
                 Real    a_dx,
                 int     a_icomp)
{
  Real xval = (Real(a_p[0]) + 0.5)*a_dx;
  Real rightAns = Real(a_icomp+1)*sin(2.*PI*xval);
  return rightAns;
}

void fillFAB(FArrayBox& a_fab,
             const Chombo4::Box& a_box,
             Real a_dx)
{
  for(Chombo4::BoxIterator bit(a_box); bit.ok(); ++bit)
  {
    IntVect iv = bit();
    for(int icomp = 0; icomp < a_fab.nComp(); icomp++)
    {
      Real uval = getPerValue(iv, a_dx, icomp);
      a_fab(iv, icomp) = uval;
    }
  }
    
}


void checkFAB(FArrayBox& a_fab,
              const Chombo4::Box& a_box,
              Real a_dx)
{
  for(Chombo4::BoxIterator bit(a_box); bit.ok(); ++bit)
  {
    IntVect iv = bit();
    for(int icomp = 0; icomp < a_fab.nComp(); icomp++)
    {
      Real uval = a_fab(iv, icomp);
      Real corval =  getPerValue(iv, a_dx, icomp);
      Real diff = (uval-corval)*(uval-corval);
      if(diff > 1.0e-5)
      {
        printf("wrong ans at comp %d, point (%d, %d)\n", icomp, iv[0], iv[1]);
      }
    }
  }
}
void chomboRun(const RunParams& a_params)
{

  pout() << "chombo version using FArrayBox" << endl;

  int nGhost = 4;

  using Chombo4::ProblemDomain;
  using Chombo4::DisjointBoxLayout;
  using Chombo4::LevelBoxData;
  using Chombo4::Copier;
  using Chombo4::DataIterator;
  
  IntVect domLo = IntVect::Zero;
  IntVect domHi  = (a_params.nx - 1)*IntVect::Unit;
  constexpr bool is_periodic[] = {true, true, true};
  cout << "periodic flags: " << is_periodic[0] << ", " << is_periodic[1] << endl;
  Chombo4::ProblemDomain domain(domLo, domHi, is_periodic);

  Vector<Chombo4::Box> boxes;
  unsigned int blockfactor = 8;
  domainSplit(domain, boxes, a_params.maxgrid, blockfactor);
  Vector<int> procs;
  LoadBalance(procs, boxes);
  DisjointBoxLayout grids(boxes, procs, domain);
  LevelData<FArrayBox> U(grids, NUMCOMPS, nGhost*IntVect::Unit);
  static Copier s_exchangeCopier;
  s_exchangeCopier.exchangeDefine(U.disjointBoxLayout(), U.ghostVect());

  DataIterator dit = grids.dataIterator();


  for(int ibox = 0; ibox < dit.size(); ibox++)
  {
 
    auto grid = grids[dit[ibox]];
    fillFAB(U[dit[ibox]], grid, a_params.dx);
  }
  
  U.exchange(s_exchangeCopier);
  
  for(int ibox = 0; ibox < dit.size(); ibox++)
  {
 
    auto grnbx = U[dit[ibox]].box();
    checkFAB(U[dit[ibox]], grnbx, a_params.dx);
  }
}

////////////
int main(int a_argc, char* a_argv[])
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
    protoRun(params);
    chomboRun(params);
  return 0;
}
