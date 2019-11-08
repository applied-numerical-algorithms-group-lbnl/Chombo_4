#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <cfloat>
#include <cmath>
#include <cstring>
#include "Chombo_SPACE.H"
#include "Chombo_FArrayBox.H"
#include "Chombo_MayDay.H"
#include "Chombo_CH_Timer.H"
#include "Chombo_NamespaceHeader.H"
using std::cin;
using std::cout;
using std::cerr;
using std::setw;
using std::setprecision;
using std::ios;
using std::pow;
using std::sqrt;




FArrayBox::FArrayBox()
  :
  BaseFab<Real>()
{
}

FArrayBox::FArrayBox(const Box& a_box,
                     int        a_n,
                     Real*      a_alias)
  :
  BaseFab<Real>(a_box,a_n,a_alias)
{
  // Note: All work is done in BaseFab<Real> constructor
}

FArrayBox::FArrayBox(const Box& a_box,
                     int        a_n)
  :
  BaseFab<Real>(a_box,a_n)
{

}

FArrayBox::~FArrayBox()
{

}

//-----------------------------------------------------------------------

void FArrayBox::performCopy(const BaseFab<Real>& a_src,
                            const Box&           a_srcbox,
                            int                  a_srccomp,
                            const Box&           a_destbox,
                            int                  a_destcomp,
                            int                  a_numcomp)
{
  BaseFab<Real>::performCopy(a_src, a_srcbox, a_srccomp, a_destbox, a_destcomp, a_numcomp);

}
#include "Chombo_NamespaceFooter.H"
