#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _MHDRK4_H_
#define _MHDRK4_H_

#include "MHDOp.H"
#include "Proto.H"
#include "LevelBoxData.H"


class MHDDX;
using std::shared_ptr;
using Proto::BoxData;

///
class MHDState
{
public:
  MHDState(){};
  ~MHDState(){};
  MHDState(shared_ptr<LevelBoxData<NUMCOMPS> >a_U);
  void increment(const MHDDX& a_DX);
  void copyhere(const MHDDX& a_DX);
  //for max wave speeds
  Real m_velSave;
  shared_ptr<LevelBoxData<NUMCOMPS> > m_U;
  DisjointBoxLayout m_grids;
  Copier m_exchangeCopier;
};

///
class MHDDX
{
public:
  MHDDX(){};
  ~MHDDX(){};
  shared_ptr<LevelBoxData<NUMCOMPS> > m_DU;
  DisjointBoxLayout m_grids;

  void increment(Real& a_weight,const MHDDX& a_incr);
  void init(const MHDState& a_State);
  void operator*=(const Real& a_weight);

};

///

///
class MHDRK4Op
{

public:
  void operator()(MHDDX& a_DX,
                  Real a_time,
                  Real a_dt,
                  MHDState& a_State) const;
  
  static Real maxWave(MHDState& a_State);
};

///

///
class MHDEulerOp
{

public:
  void operator()(MHDDX& a_DX,
                  Real a_time,
                  Real a_dt,
                  MHDState& a_State) const;
};

///
class MHDrhsOp
{

public:
  void operator()(MHDDX& a_DX,
                  MHDState& a_State) const;
};

#endif
