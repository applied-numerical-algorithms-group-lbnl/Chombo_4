#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include <cstring>

#include "BaseFab.H"
#include "BoxIterator.H"
#include "NamespaceHeader.H"

Real BaseFabRealSetVal = BASEFAB_REAL_SETVAL;

template < > int BaseFab<int>::test()
{
  int retbox = testBoxAndComp();
  if (retbox != 0)
    {
      pout() << "testboxandcomp failed" << endl;
      return retbox;
    }

  Box b(IntVect::Zero, IntVect::Unit);
  BaseFab<int> blerg;
  blerg.define(b, 1);
  int val = 4;
  blerg.setVal(val);
  for (BoxIterator bit(b); bit.ok();++bit)
    {
      if (blerg(bit(), 0) != val)
        {
          pout() << "setval or index busted" << endl;
          return -1;
        }
    }

  Box bcop(IntVect::Unit, 2*IntVect::Unit);
  BaseFab<int> bfcopy(bcop, 1);
  bfcopy.setVal(2*val);
  bfcopy.copy(blerg, 0, 0, 1);
  Box binter =bcop;
  binter &= b;
  for (BoxIterator bit(binter); bit.ok();++bit)
    {
      if (bfcopy(bit(), 0) != val)
        {
          pout() << "copy busted" << endl;
          return -2;
        }
    }

  return 0;
}
template < > void BaseFab<Real>::define()
{
  CH_assert(m_nvar > 0);
  CH_assert(m_dptr == 0);
  // CH_assert(m_numpts > 0);
  // petermc, 8 Dec 2010:  Allow empty BaseFab<Real>.
  if (m_numpts == 0) return;
  CH_assert(!m_aliased);

  m_truesize = m_nvar * m_numpts;
  m_dptr     = static_cast<Real*>(malloc(m_truesize * sizeof(Real)));

#ifdef CH_USE_SETVAL
  setVal(BaseFabRealSetVal);
#endif
}

template < > void BaseFab<int>::define()
{
  CH_assert(m_nvar > 0);
  CH_assert(m_dptr == 0);
  // CH_assert(m_numpts > 0);
  // petermc, 10 Apr 2012:  Allow empty BaseFab<int>.
  CH_assert(!m_aliased);


  m_truesize = m_nvar * m_numpts;
  m_dptr     = static_cast<int*>(malloc(m_truesize * sizeof(int)));
}

template < > void BaseFab<Real>::undefine()
{
  if (m_aliased)
  {
    m_dptr = 0;
    return;
  }

  if (m_dptr == 0)
  {
    return;
  }

  free(m_dptr);


  m_dptr = 0;
}

template < > void BaseFab<int>::undefine()
{
  if (m_aliased)
  {
    m_dptr = 0;
    return;
  }

  if (m_dptr == 0)
  {
    return;
  }

  free(m_dptr);

  m_dptr = 0;
}

template < > void BaseFab<Real>::setVal(Real a_val)
{
  if (a_val == 0)
  {
    memset(m_dptr, 0, m_truesize*sizeof(Real));
  }
  else
  {
    Real* end = m_dptr + m_truesize;
    for (Real* v = m_dptr; v<end; v++)
    {
      *v = a_val;
    }
  }
}
#include "NamespaceFooter.H"
