#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _Chombo_FABVIEW_H_
#define _Chombo_FABVIEW_H_

#include "Chombo_AMRIO.H"
#include "Chombo_NodeAMRIO.H"
#include "Chombo_NamespaceHeader.H"

int InitWriteFAB()
{
#ifdef CH_USE_HDF5
  writeFAB(NULL);

  viewFAB(NULL);

  viewBFI(NULL);

  viewBFIV(NULL);

  viewBFRV(NULL);

  viewIVSFAB(NULL);

  viewIVSFABI(NULL);

  viewIVSFABIV(NULL);

  viewIVSFABRV(NULL);

  viewCFAB(NULL);

  browseFAB(NULL);

  browseLevel(NULL);

  writeBFR(NULL);

  viewBFR(NULL);

  writeFABname(NULL, "label");

  writeLevel(NULL);

  viewLevel(NULL);

  viewLevelNoGhost(NULL);

  writeLevelname(NULL, "label");

  writeDBL(NULL);

  viewDBL(NULL);

  writeDBLname(NULL, "label");

  viewLevelBFI(NULL);

  viewIVS(NULL);

  viewVectorBox(NULL);

  writeNFAB(NULL);

  viewNFAB(NULL);

  writeNFABname(NULL, "label");

  writeNodeLevel(NULL);

  viewNodeLevel(NULL);

  writeNodeLevelname(NULL, "label");

#endif // CH_USE_HDF5

  return 1;
}

#include "Chombo_NamespaceFooter.H"
#endif
