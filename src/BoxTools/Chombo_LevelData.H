#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _Chombo_LEVELDATA_H_
#define _Chombo_LEVELDATA_H_

#include "Chombo_IntVect.H"
#include "Chombo_BoxLayoutData.H"
#include "Chombo_DisjointBoxLayout.H"
#include "Chombo_Copier.H"
#include "Chombo_SPMD.H"
#include "Chombo_NamespaceHeader.H"

///Data over a disjoint union of rectangles

/**
   LevelData is-a BoxLayoutData.  It is built from a DisjointBoxLayout
   instead though, thus making the data in a DisjointBoxLayout uniquely
   defined for all spatial indices defined by the underlying
   DisjointBoxLayout.  It carries within it an inherent concept of
   ghost cells (which the user can choose to be empty).

   Since LevelData is-a BoxLayoutData, all the methods
   required of the template class T for BoxLayoutData are also
   required for LevelData. LevelData does not have any extra requirements.
*/

template <class T>
void aliasLevelData(LevelData<T>& a_alias,
                    LevelData<T>* a_original,
                    const Interval& a_interval);

template<class T> class LevelData : public BoxLayoutData<T>
{
public:
#ifdef PROTO_CUDA  
  //just to shut it up about partially overwritten virtual functions
  using BoxLayoutData<T>::define;
#endif
  ///
  LevelData();

  ///
  LevelData(const DisjointBoxLayout& dp, int comps,
            const IntVect& ghost = IntVect::Zero,
            const DataFactory<T>& a_factory = DefaultDataFactory<T>());

  ///
  virtual ~LevelData();


  ///
  void define(const DisjointBoxLayout& dp, int comps,
              const IntVect& ghost = IntVect::Zero,
              const DataFactory<T>& a_factory = DefaultDataFactory<T>());

  ///
  /**
     copy definer.  this LevelData thrown away and da's LevelData copied
  */
  void define(const LevelData<T>& da,
              const DataFactory<T>& a_factory = DefaultDataFactory<T>());

  ///
  /**
     Copy definer.  'this' LevelData thrown away and da's LevelData copied.
     This LevelData is now defined for the range [0, comps.size()] and filled
     with the data in da from [comps.begin(), comps.end()]
  */
  void define(const LevelData<T>& da, const Interval& comps,
              const DataFactory<T>& a_factory = DefaultDataFactory<T>());

  ///
  void copyTo(const Interval& srcComps,
              BoxLayoutData<T>& dest,
              const Interval& destComps) const;


  ///only works if source and dest have the same disjointboxlayout
  void localCopyTo(const Interval& srcComps,
                   LevelData<T>  & dest,
                   const Interval& destComps) const;

  /// assumes source and dest have same interval
  void localCopyTo(LevelData<T>  & dest) const;


  /// Simplest case -- assumes source and dest have same interval
  void copyTo(BoxLayoutData<T>& dest) const;

  /*
 ///
 void copyTo(const Interval& srcComps,
 BoxLayoutData<T>& dest,
 const Interval& destComps,
 const Copier& copier) const;

 /// Simplest case -- assumes source and dest have same interval
 void copyTo(BoxLayoutData<T>& dest,
 const Copier& copier) const;
  */
  ////////////////////////////// FM MOD START //////////////////////////
  ///
  /** same as copyTo that takes a BoxLayoutData, except that it fills the
      ghost cells of 'dest' with data from 'this' also. USer passes in
      a prebuilt Copier object*/
  void copyTo(const Interval& srcComps,
              BoxLayoutData<T>& dest,
              const Interval& destComps,
              const Copier& copier,
              const LDOperator<T>& a_op = LDOperator<T>()) const;

  /// Simplest case -- assumes source and dest have same interval
  void copyTo(BoxLayoutData<T>& dest,
              const Copier& copier,
              const LDOperator<T>& a_op = LDOperator<T>()) const;

  ////////////////////////////// FM MOD. END ////////////////////////////
  ///
  /** same as copyTo that takes a BoxLayoutData, except that it fills the
      ghost cells of 'dest' with data from 'this' also. */
  void copyTo(const Interval& srcComps,
              LevelData<T>& dest,
              const Interval& destComps) const;

  /// Simplest case -- assumes source and dest have same interval
  void copyTo(LevelData<T>& dest) const;

  ///
  /** same as copyTo that takes a BoxLayoutData, except that it fills the
      ghost cells of 'dest' with data from 'this' also. USer passes in
      a prebuilt Copier object*/
  void copyTo(const Interval& srcComps,
              LevelData<T>& dest,
              const Interval& destComps,
              const Copier& copier,
              const LDOperator<T>& a_op = LDOperator<T>()) const;

  /// Simplest case -- assumes source and dest have same interval
  void copyTo(LevelData<T>& dest,
              const Copier& copier,
              const LDOperator<T>& a_op = LDOperator<T>()) const;

  /// Simplest case -- do all components
  void exchange(void);

  /// Simplest case -- do all components. Accepts a pre-built copier.
  void exchange(const Copier& copier);

  /// Accepts an arbitrary component range
  void exchange(const Interval& comps);

  /// The most general case -- can accept an arbitrary component range,
  /// a pre-built Copier object, and an arbitrary accumulation operator.
  void exchange(const Interval& comps,
                const Copier& copier,
                const LDOperator<T>& a_op = LDOperator<T>());

  /// asynchronous exchange start.  load and fire off messages.
  void exchangeBegin(const Copier& copier);
  /// finish asynchronous exchange
  void exchangeEnd();

  void exchangeNoOverlap(const Copier& copier);

  ///
  const IntVect& ghostVect() const
  {
    return m_ghost;
  }

  /**
     \name overidden  functions

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
  void define(const BoxLayout& dp, int comps,
              const DataFactory<T>& factory);

  ///
  void define(const BoxLayoutData<T>& da,
              const DataFactory<T>& factory = DefaultDataFactory<T>());

  ///
  void define(const BoxLayoutData<T>& da, const Interval& comps,
              const DataFactory<T>& factory = DefaultDataFactory<T>());

  void define(const BoxLayout& deadFunction);
  /*@}*/

  ///
  const DisjointBoxLayout& getBoxes() const
  {
    return m_disjointBoxLayout;
  }

  ///
  const DisjointBoxLayout& disjointBoxLayout() const
  {
    return m_disjointBoxLayout;
  }

  ///
  /** User writes a function with the signature:

      \code
      void myfunction(const Box& box, int n_comps, T& t)
      {
      your code here;
      }
      \endcode

      They can then hand this off to LayoutData::apply, which invokes this
      function for every T.  The argument "box" is the Box (as known to the
      DisjointBoxLayout here) associated with that T and the argument "n_comps"
      is the number of components in this LevelData.

      Your function must not be inline.

      For example:
      \code
      LevelData<FArrayBox> data(layout, 3, IntVect::Unit);
      struct val
      {
      static void set1(const Box& box, int n_comps, const FArrayBox& fab)
      {
      fab.setVal( box.smallEnd(0), box, 0, n_comps );
      }
      };

      data.apply(val::set1);
      \endcode
  */
  void apply( void (*a_Function)(const Box&, int, T&) );

  /** For use with apply( const ApplyFunctor& ) */
  struct ApplyFunctor
  {
    ~ApplyFunctor()
    {
    }

    void operator()( const Box&, int, T& ) const = 0;
  };

/** Like the other apply(), but here the argument is an instance of a class
    derived from LevelData::ApplyFunctor, and which implements ApplyFunctor's
    pure   void operator()(const Box& box, int n_comps, T& t) const.

    Going with an instance of such a class is more convenient if you want
    the thing you pass to apply() to have state.

    For example:
    \code
    class MyFunctor : public LevelData<FArrayBox>::ApplyFunctor
    {
    public:
    MyFunctor( Real x ) : m_x(x)
    {
    ...
    }
    void operator()(const Box& box, int n_comps, FArrayBox& fab) const
    {
    fab.setVal( m_x, box, 0, n_comps );
    }
    private:
    const Real m_x;
    }

    LevelData<FArrayBox> data(layout, 3, IntVect::Unit);
    data.apply( MyFunctor(3.14, 0) );

    \endcode

*/
  void apply( const ApplyFunctor& );

  void degenerate( LevelData<T>& a_to, const SliceSpec& a_ss ) const;

  /// version of degenerate which does strictly local copying
  /** this means that it maintains source ghost cell values in the 
      destination*/
  void degenerateLocalOnly( LevelData<T>& a_to, const SliceSpec& a_ss ) const;

protected:
  DisjointBoxLayout m_disjointBoxLayout;

  IntVect   m_ghost;

  friend void aliasLevelData<T>(LevelData<T>& a_alias,
                                LevelData<T>* a_original,
                                const Interval& a_interval);

  Copier m_exchangeCopier;
};

/// LevelData aliasing function
/**
   @param a_alias aliased LevelData<T> object.  original data in a_alias is destroyed and new aliasing role assumed.
   @param a_original pointer to LevelData<T> that will be aliased.
   @param a_interval range of components of each T in a_original that will be created in the a_alias argument.

   \code

   LevelData<FArrayBox> original(dbl, 4, 2*IntVect::Unit);
   Interval             subcomps(2, 3);
   LevelData<FArrayBox> alias;
   aliasLevelData<FArrayBox>(alias, &original, subcomps);
// component 0 of every FArrayBox in alias references the same data as
// component 2 in original
\endcode

The template class T must have an 'alias' constructor

\code
class A
{
public:
A(const Interval& subcomps, A& original);// alias constructor
...
};

\endcode
*/
template <class T>
void aliasLevelData(LevelData<T>& a_alias, LevelData<T>* a_original,
                    const Interval& a_interval)
{
  AliasDataFactory<T> factory(a_original, a_interval);
  a_alias.define(a_original->disjointBoxLayout(), a_interval.size(), a_original->ghostVect(), factory);
}

// ====== inlined function definitions ================

// As with the BoxLayoutData implementation file.  This file
// just gives the interface and the inline implmentations
// are given in LevelDataI.H

#include "Chombo_NamespaceFooter.H"
#include "Chombo_LevelDataI.H"

#endif
