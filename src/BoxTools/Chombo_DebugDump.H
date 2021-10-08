#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _Chombo_DEBUGDUMP_H_
#define _Chombo_DEBUGDUMP_H_

#include "Chombo_DebugOut.H"
#include "Chombo_NamespaceHeader.H"

void ch4DontYouDareCallMe()
{
  dumpIVS(NULL);
  dumpIVSFAB(NULL);
  dumpBL(NULL);
  dumpLDFLoc(NULL);
  dumpLDFPar(NULL);
  dumpFAB(NULL);
  dumpDBL(NULL);
  dumpBox(NULL);
  dumpVBox(NULL);
  dumpVVBox(NULL);
  dumpLDDBL(NULL);
  dumpBFR(NULL);
  dumpBFI(NULL);
}

#include "Chombo_NamespaceFooter.H"
#endif
