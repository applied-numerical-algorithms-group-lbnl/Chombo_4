#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _LEVELBoxDATA_H__
#define _LEVELBoxDATA_H__

#include "LevelData.H"
#include "Proto.H"
#include "NamespaceHeader.H"

using Proto::BoxData;
///
/**
*/
template<int NUMCOMPS>
class LevelBoxData: public LevelData<BoxData<double, NUMCOMPS> >
{
public:

  typedef DefaultDataFactory<BoxData<double, NUMCOMPS> > factory_t;
  ///
  LevelBoxData(): LevelData<BoxData<double, NUMCOMPS> >()
  {
  }

  ///
  LevelBoxData(const DisjointBoxLayout & a_dp, 
               const IntVect           & a_ghost   = IntVect::Zero,
               const DataFactory<T>    & a_factory = factory_t())
    : LevelData<BoxData<double, NUMCOMPS> >(a_dp, NUMCOMPS, a_ghost, a_factory)
  {
  }

  ///
  virtual ~LevelBoxData()
  {
  }

  ///
  virtual void define(const DisjointBoxLayout& a_dp, 
                      const IntVect          & a_ghost = IntVect::Zero,
                      const DataFactory<T>   & a_factory = factory_t())

  {
    LevelData<BoxData<double, NUMCOMPS> >::define(a_dp, NUMCOMPS, a_ghost, a_factory);
  }


  ///
  virtual void copyTo(const Interval&         a_srcComps,
                      LevelBoxData<NUMCOMPS>& a_dst,
                      const Interval&         a_dstComps) const

  {
    LevelData<FArrayBox> hostdataSrc(      m_disjointBoxLayout, NUMCOMPS,        m_ghost, factory_t());
    LevelData<FArrayBox> hostdataDst(a_dst.m_disjointBoxLayout, NUMCOMPS, a_dst.m_ghost, factory_t());
    this->copyToHost( hostdataSrc);
    a_dst->copyToHost(hostdataDst);
    
    hostdataSrc.copyTo(a_srcComps, hostdataDst, a_dstComps);
    
    hostDataDst.copyToDevice(a_dst);
  }


  /// Simplest case -- assumes source and dst have same interval
  virtual void copyTo(LevelBoxData<NUMCOMPS>& a_dst) const
  {
    Interval interv(0, NUMCOMPS-1);
    copyTo(interv, a_dst, interv);
  }

  /*
  ///
  virtual void copyTo(const Interval& srcComps,
                      BoxLayoutData<T>& dst,
                      const Interval& dstComps,
                      const Copier& copier) const;

  /// Simplest case -- assumes source and dst have same interval
  virtual void copyTo(BoxLayoutData<T>& dst,
                      const Copier& copier) const;
  */
  ////////////////////////////// FM MOD START //////////////////////////
   ///
  /** same as copyTo that takes a BoxLayoutData, except that it fills the
    ghost cells of 'dst' with data from 'this' also. USer passes in
    a prebuilt Copier object*/
  virtual void copyTo(const Interval& srcComps,
                      BoxLayoutData<T>& dst,
                      const Interval& dstComps,
                      const Copier& copier,
                      const LDOperator<T>& a_op = LDOperator<T>()) const;

  /// Simplest case -- assumes source and dst have same interval
  virtual void copyTo(BoxLayoutData<T>& dst,
                      const Copier& copier,
                      const LDOperator<T>& a_op = LDOperator<T>()) const;

  ////////////////////////////// FM MOD. END ////////////////////////////
  ///
  /** same as copyTo that takes a BoxLayoutData, except that it fills the
    ghost cells of 'dst' with data from 'this' also. */
  virtual void copyTo(const Interval& srcComps,
                      LevelData<T>& dst,
                      const Interval& dstComps) const;

  /// Simplest case -- assumes source and dst have same interval
  virtual void copyTo(LevelData<T>& dst) const;

   ///
  /** same as copyTo that takes a BoxLayoutData, except that it fills the
    ghost cells of 'dst' with data from 'this' also. USer passes in
    a prebuilt Copier object*/
  virtual void copyTo(const Interval& srcComps,
                      LevelData<T>& dst,
                      const Interval& dstComps,
                      const Copier& copier,
                      const LDOperator<T>& a_op = LDOperator<T>()) const;

  /// Simplest case -- assumes source and dst have same interval
  virtual void copyTo(LevelData<T>& dst,
                      const Copier& copier,
                      const LDOperator<T>& a_op = LDOperator<T>()) const;

  /// Simplest case -- do all components
  virtual void exchange(void);

  /// Simplest case -- do all components. Accepts a pre-built copier.
  virtual void exchange(const Copier& copier);

  /// Accepts an arbitrary component range
  virtual void exchange(const Interval& comps);

  /// The most general case -- can accept an arbitrary component range,
  /// a pre-built Copier object, and an arbitrary accumulation operator.
  virtual void exchange(const Interval& comps,
			const Copier& copier,
			const LDOperator<T>& a_op = LDOperator<T>());

  /// asynchronous exchange start.  load and fire off messages.
  virtual void exchangeBegin(const Copier& copier);
  /// finish asynchronous exchange
  virtual void exchangeEnd();

  virtual void exchangeNoOverlap(const Copier& copier);

  ///
  const IntVect& ghostVect() const
  {
    return m_ghost;
  }

  /**
     \name overidden virtual functions

    These functions will invoke error messages when invoked.
    C++ will ensure that constructors are not called for the
    base class by a user, but a 'define' function has no such protection,
    hence the need to prevent such usage.  A runtime error is not
    a perfect solution...(strong construction gets around this  *sigh*).
    classes that derive from LevelData will have to turn its valid
    defines into runtime errors also and make its own defines.  Thus
    taking over the job of the compiler.
    */
  /*@{*/
  virtual void define(const BoxLayout& dp, int comps,
                      const DataFactory<T>& factory);

  ///
  virtual void define(const BoxLayoutData<T>& da,
                      const DataFactory<T>& factory = DefaultDataFactory<T>());

  ///
  virtual void define(const BoxLayoutData<T>& da, const Interval& comps,
                      const DataFactory<T>& factory = DefaultDataFactory<T>());

  virtual void define(const BoxLayout& deadFunction);

protected:
  DisjointBoxLayout m_disjointBoxLayout;

  IntVect   m_ghost;

  Copier m_exchangeCopier;
};


#include "NamespaceFooter.H"

#endif
