#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _Chombo_BRMESHREFINE_H_
#define _Chombo_BRMESHREFINE_H_

#include <climits>
#include <list>

#include <vector>
#include "Chombo_Box.H"

#include "Chombo_REAL.H"

#include "Chombo_Interval.H"
#include "Chombo_NamespaceHeader.H"


///
/**
   Splits domain into a vector of disjoint boxes with
   maximum size maxsize.
   blocking factor is default to one.
   If you make minimum size > 1, then domain must
   be coarsenable and refineable by blockfactor
   (refine(coarsen(domain,blockfactor), minsize) == domain)
   or an error is thrown. This would be defined in
   MeshRefine.H, except that it needs to use a BRMeshRefine object.
   Here a_refineDirs[d] is 1 if refining in dimension d, and 0 if not.
*/
extern void domainSplit(const Box& a_domain, std::vector<Box>& a_vbox,
                        int a_maxSize, int a_blockfactor=1,
                        IntVect a_refineDirs=IntVect::Unit);

///
/**
   Recursive function to enforce max size of boxes in a given direction.
 */
extern void
breakBoxes(std::vector<Box>& a_vboxin,  const int& a_maxSize, const int& a_idir);

#include "Chombo_NamespaceFooter.H"
#endif
