#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _EBLevelFluxData_H____
#define _EBLevelFluxData_H____

#include "Proto.H"
#include "EBProto.H"
#include "Chombo_LevelData.H"
#include "Chombo_AMRIO.H"
#include "Chombo_ProtoInterface.H"
#include "Chombo_EBChombo.H"
#include "Chombo_ProtoFactories.H"
#include "Chombo_LevelBoxData.H"
#include "Chombo_NamespaceHeader.H"


///
/**
*/
template<int  NUMCOMPS>
class EBLevelFluxData: public LevelData<EBFluxData<Real, NUMCOMPS> >
{
public:

  typedef GraphConstructorFactory<EBFluxData<Real, NUMCOMPS> > devifactory_t;

  ///
  EBLevelFluxData(): LevelData<EBFluxData<Real, NUMCOMPS> >()
  {
  }


  ///
  EBLevelFluxData(const DisjointBoxLayout        & a_dp, 
                  const IntVect                  & a_ghost,
                  shared_ptr<LevelData<EBGraph>  > a_graphs)
    : LevelData<EBFluxData<Real, NUMCOMPS> >(a_dp, NUMCOMPS, a_ghost, devifactory_t(a_graphs)), m_graphs(a_graphs)
  {
  }

  ///
  virtual ~EBLevelFluxData()
  {
  }


  ///
  void define(const DisjointBoxLayout& a_dp, 
              const IntVect          & a_ghost,
              shared_ptr<LevelData<EBGraph>  > a_graphs)
  {
    m_graphs = a_graphs;
    LevelData<EBFluxData<Real, NUMCOMPS> >::define(a_dp, NUMCOMPS, a_ghost, devifactory_t(a_graphs));
  }


  void setVal(Real a_value)
  {
    DataIterator dit = this->dataIterator();
    for(int ibox = 0;  ibox < dit.size(); ibox++)
    {
      (*this)[dit[ibox]].setVal(a_value);
    }
  }
  
  ///
  void copyTo(const Interval                      &            a_srcComps,
              LevelData<BoxData<Real, NUMCOMPS> > &            a_dst,
              const Interval                      &            a_dstComps,
              const Copier                        &            a_copier) const
  {
    MayDay::Error("not implemented---not sure if we need this");
  }


  /// 
  virtual void exchange(const Copier& copier)
  {

    MayDay::Error("not implemented---not sure if we need this");
  }

private:
    shared_ptr<LevelData<EBGraph>  > m_graphs;
};


#include "Chombo_NamespaceFooter.H"

#endif
