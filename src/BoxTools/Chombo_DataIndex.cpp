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
#include "Chombo_MayDay.H"
#include "Chombo_parstream.H"
#include "Chombo_DataIndex.H"
#include "Chombo_NamespaceHeader.H"

std::ostream&
operator<< (std::ostream&   os,            const LayoutIndex& dit)
{
  os << " (" << dit.m_index << "," << dit.m_datInd << ") " ;
  return os;
}
#include "Chombo_NamespaceFooter.H"
