#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _Chombo_CH_ASSERT_H__
#define _Chombo_CH_ASSERT_H__

#include "Proto.H"

#include "Chombo_BaseNamespaceHeader.H"

#ifndef CH_assert
#define CH_assert(cond) PR_assert(cond)
#endif

#include "Chombo_BaseNamespaceFooter.H"

#endif /* include guard*/
