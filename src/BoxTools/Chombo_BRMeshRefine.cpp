#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include <fstream>
#include <iostream>
#include "Chombo_SPMD.H"
#include "Chombo_CH_Timer.H"

// Constants for Berger-Rigoutsos algorithm
#if ! defined(_BR_MIN_INFLECTION_MAG_)
#define _BR_MIN_INFLECTION_MAG_ ( 3 )
#endif

#define MAXBOXES     2000000 // petermc: was 20000
#define P_BUFFERSIZE (MAXBOXES * 8 * CH_SPACEDIM)

// Berger-Rigoutsos Mesh refinement class
// ---------
//  class BRMeshRefine
//
///////////////////////////////////////////////////////////////////////////////

// Include files:

#include "Chombo_BRMeshRefine.H"
#include "Chombo_MayDay.H"
#include "Chombo_parstream.H"
#include "Chombo_NamespaceHeader.H"
//recursive function to enforce max size of boxes in a given direction
void
breakBoxes(Vector<Box>& a_vboxin,  const int& a_maxBoxSize, const int& a_idir)
{
  int nboxes = a_vboxin.size();
  //need to use STL vector for bools.
  using std::vector;
  vector<bool> splitDec(nboxes);
  bool anyBoxesToSplit = false;
  //find out which boxes need to be chopped and in what direction
  for (int ibox = 0; ibox < nboxes; ibox++)
    {
      if ( a_vboxin[ibox].size(a_idir ) > a_maxBoxSize )
        {
          splitDec[ibox] = true;
          anyBoxesToSplit = true;
        }
      else
        {
          splitDec[ibox] = false;
        }
    }
  //if there are no boxes to split, just return
  //otherwise, split all the boxes that need to be
  //split ONCE and then call function recursively
  // and set the return vector to the temporary after
  // the recursion
  if (anyBoxesToSplit)
    {
      Vector<Box> vboxtemp;
      for (int ibox = 0; ibox < nboxes; ibox++)
        {
          Box boxone = a_vboxin[ibox];
          if (splitDec[ibox])
            {
              //int len = (boxone.smallEnd(a_idir) + boxone.bigEnd(a_idir))/2;
              int mid = boxone.smallEnd(a_idir)+a_maxBoxSize ;
              Box boxtwo = boxone.chop(a_idir, mid);
              vboxtemp.push_back(boxone);
              vboxtemp.push_back(boxtwo);
            }
          else
            {
              vboxtemp.push_back(boxone);
            }
        }
      breakBoxes(vboxtemp, a_maxBoxSize, a_idir);
      a_vboxin = vboxtemp;
    }
  return;
}

void
domainSplit(const Box& a_domain, Vector<Box>& a_vbox, int a_maxBoxSize, int a_blockFactor, IntVect a_refineDirs)
{
  CH_TIME("BRMeshRefine::domainSplit");
  a_vbox.resize(0);
  if (a_maxBoxSize == 0)
    {
      a_vbox.push_back(a_domain);
      return;
    }
  int ratio = a_maxBoxSize/a_blockFactor;

  // Set blockFactorDirs[d] = a_blockFactor if a_refineDirs[d] == 1; else 1.
  IntVect blockFactorDirs = (a_blockFactor-1)*a_refineDirs + IntVect::Unit;
  Box d(a_domain);
  d.coarsen(blockFactorDirs);
  if (refine(d, blockFactorDirs) != a_domain)
    {
      MayDay::Error("domainSplit: a_domain not coarsenable by blockingFactor");
    }
  
  /*
    //recursive version here is too cute for very large box counts (bvs)
    // apparently testRThetaZ needs domainSplit to work thsi way.  I suspect the test needs a larger block factor
  a_vbox.push_back(d);
  for (int i=0; i<CH_SPACEDIM; ++i)
    {
      if (a_refineDirs[i] == 1)
        {
          breakBoxes(a_vbox, ratio, i);
        }
    }
   for (int b=0; b<a_vbox.size(); ++b)
    {
      a_vbox[b].refine(blockFactorDirs);
    }
  */

  size_t  count[CH_SPACEDIM+1] = {1};

  for(int i=0; i<CH_SPACEDIM; ++i)
    {
      count[i+1]=count[i];
      if(a_refineDirs[i]==1)
        {
          int c = (d.bigEnd()[i]-d.smallEnd()[i]+ratio)/ratio;
          count[i+1]*=c;
        }
    }
  a_vbox.resize(count[CH_SPACEDIM]);

  //#pragma omp parallel for  
  for (int b=0; b<a_vbox.size(); ++b)
    {
      IntVect Lo(d.smallEnd()), Hi(d.bigEnd());
      for(int i=0; i<CH_SPACEDIM ; i++)
        {
          if(a_refineDirs[i]==1)
            {
              int m = (b/count[i])%(count[i+1]/count[i]);
              Lo[i]+=m*ratio;
              Hi[i]=std::min(d.bigEnd()[i],Lo[i]+ratio-1);
            }
        }
 
      a_vbox[b]=Box(Lo, Hi);
      a_vbox[b].refine(blockFactorDirs);
    }
   
}


#include "Chombo_NamespaceFooter.H"
