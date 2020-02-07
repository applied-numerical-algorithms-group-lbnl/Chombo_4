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


class ScalarFluxAliasFactory: public DataFactory< EBFluxData<Real, 1> >
{
public:

  ScalarFluxAliasFactory(LevelData<EBFluxData<Real, DIM> >* a_source, 
                         unsigned int                       a_whichcomp)
  {
    m_source    = a_source;
    m_whichcomp = a_whichcomp;
  }
  virtual ~ScalarFluxAliasFactory()
  {
  }

  /// factory function.  creates a new 'T' object
  /** creates a new 'T' object and returns a pointer to it.
   */
  virtual EBFluxData<Real, 1> * create(const Box& a_box, int ncomps, const DataIndex& a_datInd) const
  {
    EBFluxData<Real, DIM>& vecholder = (*m_source)[a_datInd];
    EBFluxData<Real, 1  >* retval = new EBFluxData<Real, 1  >();
    //aliasing define
    retval->define<DIM>(vecholder, m_whichcomp);
    return retval;
  }

  virtual bool callDelete() const
  {
    return false;
  }

  //return true if copy and construction are thread safe
  virtual bool threadSafe() const
  {
    return false;
  }

private:
  ScalarFluxAliasFactory();
  LevelData<EBFluxData<Real, DIM> >* m_source;
  unsigned int                       m_whichcomp;

};
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

  ///aliasing define
  void define(EBLevelFluxData<DIM>       &     a_input, 
              unsigned int                     a_comp,
              shared_ptr<LevelData<EBGraph>  > a_graphs)
  {
    m_graphs = a_graphs;
    auto dp    = a_input.disjointBoxLayout();
    auto ghost = a_input.ghostVect();
    m_graphs = a_graphs;
    LevelData<EBFluxData<Real, DIM> >* inputPtr = static_cast<LevelData<EBFluxData<Real, DIM> > * >(&a_input);
    LevelData<EBFluxData<Real, NUMCOMPS> >::
      define(dp, NUMCOMPS, ghost, ScalarFluxAliasFactory(inputPtr, a_comp));
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
