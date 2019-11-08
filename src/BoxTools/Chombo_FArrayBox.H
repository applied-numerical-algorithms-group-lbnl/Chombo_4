#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _Chombo_FARRAYBOX_H_
#define _Chombo_FARRAYBOX_H_

#ifndef WRAPPER
#include <iostream>

// #include <Pointers.H>
#include "Chombo_Box.H"
#include "Chombo_BaseFab.H"
#endif

#include "Chombo_REAL.H"
#include "Chombo_SPACE.H"
#include "Chombo_NamespaceHeader.H"


///
/**
   Fortran Array Boxes (generally called FABs) are objects constructed
   to interface with arrays in Fortran.  Useful operations can be performed
   upon FABs in C++, and they provide a convenient interface to
   Fortran when it is necessary to retreat into that language for
   doing arithmetic operations when performance matters.

   FArrayBox is derived from BaseFab<Real>.
   FArrayBox adds additional useful capabilities which make sense
   for Real types, such as I/O and L**p norms.

   The C pre-processor macro `CH_SPACEDIM' must be defined to use
   this class.  The internal precision of FArrayBox objects is
   set by defining either `CH_USE_FLOAT' or `CH_USE_DOUBLE'

   This class does NOT provide a copy constructor or assignment operator.
*/
class FArrayBox: public BaseFab<Real>
{
public:

  ///
  /**
     Constructs an invalid FArrayBox with no memory.
  */
  FArrayBox ();

  ///
  /**
     Constructs an initial FArrayBox with the data space allocated but not
     inititialized. a_ncomp is the number of components (variables) at each
     data point in the Box.
  */
  FArrayBox(const Box& a_box,
            int        a_ncomp,
            Real*      a_alias);

  FArrayBox(const Box& a_box,
            int        a_ncomp);
  ///
  /**
     Construct an aliaed FArrayBox.  See BaseFab class for details.
  */
  FArrayBox(const Interval& a_comps,
            FArrayBox&      a_original)
    :
    BaseFab<Real>(a_comps, a_original)
  {}

  ///
  /**
     Defines FArrayBox with the data space allocated but not
     inititialized. a_ncomp is the number of components (variables) at each
     data point in the Box.
  */
  virtual void define(const Box& a_box,
                      int        a_ncomp,
                      Real*      a_alias = NULL)
  {
    BaseFab<Real>::define(a_box, a_ncomp, a_alias);
  }

  ///
  /**
     The (virtual) destructor.
  */
  virtual ~FArrayBox ();

  ///
  /**
     Constructs an 'aliased' BaseFab of the requested interval of the
     argument BaseFab.  This BaseFab does not allocate any memory, but
     sets its data pointer into the memory pointed to by the argument
     BaseFab.  It is the users responsiblity to ensure this aliased
     BaseFab is not used after the original BaseFab has deleted its data ptr
     (resize, define(..) called, or destruction, etc.).

     This aliased BaseFab will also generate side effects (modifying the values
     of data in one will modify the other's data).

     This aliased BaseFab will have a_comps.size() components, starting at zero.
  */
  virtual void define(const Interval& a_comps,
                      FArrayBox&      a_original)
  {
    BaseFab<Real>::define(a_comps, a_original);
  }

  ///
  /**
   * This is here only to make the Intel compiler stop warning about
   * partial override.
   */
  virtual void define(const Interval& a_comps,
              BaseFab<Real>&      a_original)
  {
    BaseFab<Real>::define(a_comps, a_original);
  }

  ///
  /** override resize.  I'm going to try and get rid of this function in the future (bvs)
   */
  virtual void resize(const Box& a_box, 
                      int a_comps,
                      Real* a_alias = NULL)
  {
    BaseFab<Real>::resize(a_box, a_comps, a_alias);
  }



  FArrayBox(FArrayBox&& a_in)=default;
  FArrayBox& operator=(FArrayBox&& a_in)=default;
protected:
  virtual void performCopy(const BaseFab<Real>& a_src,
                           const Box&           a_srcbox,
                           int                  a_srccomp,
                           const Box&           a_destbox,
                           int                  a_destcomp,
                           int                  a_numcomp);


private:
  //
  // These are disallowed.
  //
  //FArrayBox (const FArrayBox&);
  //FArrayBox& operator = (const FArrayBox&);    
};

#include "Chombo_NamespaceFooter.H"
#endif
