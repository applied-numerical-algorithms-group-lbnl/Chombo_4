#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "EulerRK4.H"
#include "Chombo_DataIterator.H"
#include "Chombo_ProtoInterface.H"
#include "Chombo_NamespaceHeader.H"
//using   ::Proto::Point;


/****/
EulerState::
EulerState(shared_ptr<LevelBoxData<NUMCOMPS> > a_U)
{
  m_U     = a_U;
  m_grids = m_U->disjointBoxLayout();
}
/****/
void 
EulerState::
increment(const EulerDX & a_DX)
{
  CH_TIME("EulerState::increment");
  LevelBoxData<NUMCOMPS> & data  = *m_U;
  LevelBoxData<NUMCOMPS> & delta = *(a_DX.m_DU);
  
  DataIterator dit = m_grids.dataIterator();

  for(int ibox = 0; ibox < dit.size(); ibox++)
  {
    data[dit[ibox]] += delta[dit[ibox]];
  }
}
/****/
void 
EulerDX::
increment(Real        & a_weight,
          const EulerDX & a_DX)
{
  CH_TIME("EulerDX::increment");
  LevelBoxData<NUMCOMPS>& data  = *m_DU;
  LevelBoxData<NUMCOMPS>& delta = *(a_DX.m_DU);
  
  DataIterator dit = m_grids.dataIterator();

 for(int ibox = 0; ibox < dit.size(); ibox++)
  {
    BoxData<Real, NUMCOMPS> incr(delta[dit[ibox]].box());
    delta[dit[ibox]].copyTo(incr); 
    incr *= a_weight;
    data[dit[ibox]] += incr;
  }
}

/****/
void 
EulerDX::
init(const EulerState& a_State)
{
  CH_TIME("EulerDX::init");
  m_grids = a_State.m_grids;

  if(m_DU == NULL)
    m_DU = shared_ptr<LevelBoxData<NUMCOMPS> >(new LevelBoxData<NUMCOMPS>(m_grids, a_State.m_U->ghostVect()));

  if(U_ave == NULL)
    U_ave = shared_ptr<LevelBoxData<NUMCOMPS> >(new LevelBoxData<NUMCOMPS>(m_grids, a_State.m_U->ghostVect()));

  DataIterator dit = m_grids.dataIterator();

  for(int ibox = 0; ibox < dit.size(); ibox++)
  {
    (*m_DU)[dit[ibox]].setVal(0.);
  }
}

/****/
void 
EulerDX::
operator*=(const Real& a_weight)
{
  CH_TIME("EulerDX::operator*=");
  DataIterator dit = m_grids.dataIterator();
//#pragma omp parallel
  for(int ibox = 0; ibox < dit.size(); ibox++)
  {
    (*m_DU)[dit[ibox]] *= a_weight;
  }
}

/****/
void 
EulerRK4Op::
operator()(EulerDX& a_DX,
           Real a_time,
           Real a_dt,
           EulerState& a_State) const
{
  CH_TIMERS("EulerRKOp::operator()");
  CH_TIMER("defining leveldatas",tdef);
  CH_TIMER("copying to temporary",tcop);
  CH_TIMER("RK arithmetic_no_comm",  trk);

  CH_START(tdef);
  int ncomp  =  a_State.m_U->nComp();
  IntVect gv =  a_State.m_U->ghostVect();
  DisjointBoxLayout grids = a_State.m_grids;

  LevelBoxData<NUMCOMPS>&  delta = *(a_DX.m_DU);
  CH_STOP(tdef);

  CH_START(tcop);
  //Interval interv(0, ncomp-1);
  //Copier copier(grids, grids);
  
  LevelBoxData<NUMCOMPS>&  U_ave= *(a_DX.U_ave); 
  U_ave.define(*(a_State.m_U)); // ask why

  CH_STOP(tcop);

  CH_START(trk);
  DataIterator dit = grids.dataIterator();

  for(int ibox = 0; ibox < dit.size(); ibox++)
  {
    U_ave[dit[ibox]] += delta[dit[ibox]];
  }

  EulerOp::step(*a_DX.m_DU, U_ave, a_State.m_Rxn);

  for(int ibox = 0; ibox <  dit.size(); ibox++)
  {
    delta[dit[ibox]] *= a_dt;
  }
  CH_STOP(trk);
}


Real
EulerRK4Op::
maxWave(EulerState& a_State)
{
  CH_TIME("EulerRKOp::maxwave_init");
  EulerDX DX;
  DX.init(a_State);
  Reduction<Real,Op::Abs>& rxn = a_State.m_Rxn;
  rxn.reset(); // initialize device pointer
  EulerOp::step(*(DX.m_DU),*(a_State.m_U), rxn);
  Real velmax = rxn.fetch();
  return velmax;
}
#include "Chombo_NamespaceFooter.H"
