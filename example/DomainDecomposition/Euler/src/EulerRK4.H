#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _EULERRK4_H_
#define _EULERRK4_H_

#include "EulerOp.H"
#include "Proto.H"
#include "Chombo_LevelBoxData.H"
#include "Chombo_NamespaceHeader.H"

class EulerDX;
using std::shared_ptr;
using ::Proto::BoxData;

///
class EulerState
{
public:
  EulerState(){};
  ~EulerState(){};
  EulerState(shared_ptr<LevelBoxData<NUMCOMPS> >a_U);
  void increment(const EulerDX& a_DX);
  //for max wave speeds
  Real m_velSave;
  shared_ptr<LevelBoxData<NUMCOMPS> > m_U;
  DisjointBoxLayout m_grids;
  Copier m_exchangeCopier;
  Reduction<Real,AbsMax> m_Rxn;
};

///
class EulerDX
{
public:
  EulerDX(){};
  ~EulerDX(){};
  shared_ptr<LevelBoxData<NUMCOMPS> > m_DU;
  DisjointBoxLayout m_grids;

  void increment(Real& a_weight,const EulerDX& a_incr);
  void init(const EulerState& a_State);
  void operator*=(const Real& a_weight);

};

///
class EulerRK4Op
{

public:
  void operator()(EulerDX& a_DX,
                  Real a_time,
                  Real a_dt,
                  EulerState& a_State) const;
  
  static Real maxWave(EulerState& a_State);
};

#include "Chombo_NamespaceFooter.H"
#endif
