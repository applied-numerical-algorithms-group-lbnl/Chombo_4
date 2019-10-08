#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "MHDsteps.H"
#include "DataIterator.H"
#include "ProtoInterface.H"

using   Proto::Point;
typedef Proto::Box Bx;

/****/
MHDState::
MHDState(shared_ptr<LevelBoxData<NUMCOMPS> > a_U)
{
  m_U     = a_U;
  m_grids = m_U->disjointBoxLayout();
}
/****/
void 
MHDState::
increment(const MHDDX & a_DX)
{
  CH_TIME("MHDState::increment");
  LevelBoxData<NUMCOMPS> & data  = *m_U;
  LevelBoxData<NUMCOMPS> & delta = *(a_DX.m_DU);
  
  DataIterator dit = m_grids.dataIterator();
#pragma omp parallel
  for(int ibox = 0; ibox < dit.size(); ibox++)
  {
    data[dit[ibox]] += delta[dit[ibox]];
  }
}
/****/
void 
MHDState::
copyhere(const MHDDX & a_DX)
{
  CH_TIME("MHDState::copyhere");
  LevelBoxData<NUMCOMPS> & data  = *m_U;
  LevelBoxData<NUMCOMPS> & delta = *(a_DX.m_DU);
  
  DataIterator dit = m_grids.dataIterator();
#pragma omp parallel
  for(int ibox = 0; ibox < dit.size(); ibox++)
  {
	(*m_U)[dit[ibox]].setVal(0.);
    data[dit[ibox]] += delta[dit[ibox]];
  }
}
/****/


void 
MHDDX::
increment(Real        & a_weight,
          const MHDDX & a_DX)
{
  CH_TIME("MHDDX::increment");
  LevelBoxData<NUMCOMPS>& data  = *m_DU;
  LevelBoxData<NUMCOMPS>& delta = *(a_DX.m_DU);
  
  DataIterator dit = m_grids.dataIterator();
#pragma omp parallel
  for(int ibox = 0; ibox < dit.size(); ibox++)
  {
    BoxData<Real, NUMCOMPS>& incr = delta[dit[ibox]];
    incr *= a_weight;
    data[dit[ibox]] += incr;
  }
}

/****/
void 
MHDDX::
init(const MHDState& a_State)
{
  CH_TIME("MHDDX::init");
  m_grids = a_State.m_grids;
  m_DU = shared_ptr<LevelBoxData<NUMCOMPS> >(new LevelBoxData<NUMCOMPS>(m_grids, a_State.m_U->ghostVect()));

  DataIterator dit = m_grids.dataIterator();
#pragma omp parallel
  for(int ibox = 0; ibox < dit.size(); ibox++)
  {
    (*m_DU)[dit[ibox]].setVal(0.);
  }
}

/****/
void 
MHDDX::
operator*=(const Real& a_weight)
{
  CH_TIME("MHDDX::operator*=");
  DataIterator dit = m_grids.dataIterator();
#pragma omp parallel
  for(int ibox = 0; ibox < dit.size(); ibox++)
  {
    (*m_DU)[dit[ibox]] *= a_weight;
  }
}

/****/
void 
MHDRK4Op::
operator()(MHDDX& a_DX,
           Real a_time,
           Real a_dt,
           MHDState& a_State) const
{
  CH_TIMERS("MHDRKOp::operator()");
  CH_TIMER("defining leveldatas",tdef);
  CH_TIMER("copying to temporary",tcop);
  CH_TIMER("RK arithmetic_no_comm",  trk);

  CH_START(tdef);
  int ncomp  =  a_State.m_U->nComp();
  IntVect gv =  a_State.m_U->ghostVect();
  DisjointBoxLayout grids = a_State.m_grids;
  LevelBoxData<NUMCOMPS>  U_ave(grids, gv);
  LevelBoxData<NUMCOMPS>&  delta = *(a_DX.m_DU);
  CH_STOP(tdef);

  CH_START(tcop);
  Interval interv(0, ncomp-1);
  Copier copier(grids, grids);
  a_State.m_U->copyTo(interv, U_ave, interv, copier);
  CH_STOP(tcop);


  CH_START(trk);
  DataIterator dit = grids.dataIterator();
#pragma omp parallel
  for(int ibox = 0; ibox < dit.size(); ibox++)
  {
    U_ave[dit[ibox]] += delta[dit[ibox]];
  }

  Real velmax = MHDOp::step(*a_DX.m_DU, U_ave);
  a_State.m_velSave = std::max(a_State.m_velSave,velmax);

#pragma omp parallel
  for(int ibox = 0; ibox <  dit.size(); ibox++)
  {
    delta[dit[ibox]] *= a_dt;
  }
  CH_STOP(trk);
}


Real
MHDRK4Op::
maxWave(MHDState& a_State)
{
  CH_TIME("MHDRKOp::maxwave_init");
  MHDDX DX;
  DX.init(a_State);
  Real velmax = MHDOp::step(*(DX.m_DU),*(a_State.m_U));
  return velmax;
}


/****/
void 
MHDEulerOp::
operator()(MHDDX& a_DX,
           Real a_time,
           Real a_dt,
           MHDState& a_State) const
{
  CH_TIMERS("MHDEulerOp::operator()");
  CH_TIMER("defining leveldatas",tdef);
  CH_TIMER("copying to temporary",tcop);
  CH_TIMER("RK arithmetic_no_comm",  trk);

  CH_START(tdef);
  int ncomp  =  a_State.m_U->nComp();
  IntVect gv =  a_State.m_U->ghostVect();
  DisjointBoxLayout grids = a_State.m_grids;
  LevelBoxData<NUMCOMPS>  U_ave(grids, gv);
  LevelBoxData<NUMCOMPS>&  delta = *(a_DX.m_DU);
  CH_STOP(tdef);

  CH_START(tcop);
  Interval interv(0, ncomp-1);
  Copier copier(grids, grids);
  a_State.m_U->copyTo(interv, U_ave, interv, copier);
  CH_STOP(tcop);


  CH_START(trk);
  DataIterator dit = grids.dataIterator();
  //std::cout << "Reached here 4" << std::endl;
#pragma omp parallel
  for(int ibox = 0; ibox < dit.size(); ibox++)
  {
    U_ave[dit[ibox]] += delta[dit[ibox]];
  }

  MHDOp::step2(*a_DX.m_DU, U_ave);

#pragma omp parallel
  for(int ibox = 0; ibox <  dit.size(); ibox++)
  {
    delta[dit[ibox]] *= a_dt;
  }
  CH_STOP(trk);
}



/****/
void 
MHDrhsOp::
operator()(MHDDX& a_DX,
           MHDState& a_State) const
{
  CH_TIMERS("MHDRKOp::operator()");
  CH_TIMER("defining leveldatas",tdef);
  CH_TIMER("copying to temporary",tcop);
  CH_TIMER("RK arithmetic_no_comm",  trk);

  CH_START(tdef);
  int ncomp  =  a_State.m_U->nComp();
  IntVect gv =  a_State.m_U->ghostVect();
  DisjointBoxLayout grids = a_State.m_grids;
  LevelBoxData<NUMCOMPS>  U_ave(grids, gv);
  LevelBoxData<NUMCOMPS>&  delta = *(a_DX.m_DU);
  CH_STOP(tdef);

  CH_START(tcop);
  Interval interv(0, ncomp-1);
  Copier copier(grids, grids);
  a_State.m_U->copyTo(interv, U_ave, interv, copier);
  CH_STOP(tcop);


  CH_START(trk);
  DataIterator dit = grids.dataIterator();
#pragma omp parallel
  for(int ibox = 0; ibox < dit.size(); ibox++)
  {
    U_ave[dit[ibox]] += delta[dit[ibox]];
  }

  Real velmax = MHDOp::step(*a_DX.m_DU, U_ave);

  CH_STOP(trk);
}
