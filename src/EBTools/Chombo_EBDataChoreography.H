#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _EB_DATA_CHOREOGRAPHY_H__
#define _EB_DATA_CHOREOGRAPHY_H__

#include "Chombo_DisjointBoxLayout.H"
#include "Chombo_Pool.H"
#include <vector>
#include "Chombo_ProblemDomain.H"
#include "Chombo_CH_Timer.H"
#include "Chombo_NeighborIterator.H"
#include "Chombo_DataIterator.H"

#include <unordered_map>
#include <cstdint>

#include "Chombo_SPMD.H"
#ifdef CH_MPI  
#include "mpi.h"
#endif
using namespace Chombo4;
///
/**
   The Ch4_Data_Choreography space
   (DistributedData/BoxPattern and subclasses) is meant to provide a
   reduced complexity, highly maintainable alternative to the Copier and
   LevelData infrastructure.  If you need high performance or need
   some of the fancier aspects of Copier, you should use the
   standard LevelData/BoxLayoutData/LayoutData/Copier
   infrastructure.

   Periodic boundary conditions are not supported.     

   The applications for which this is intended do not need the optimizations 
   around which LevelData is built and do not need periodic boundary conditions.

   All communication is two phase.
   All communication is done on the host.
   We aggregate all communications between two processors the same way LevelData does.  

   DistributedData holds data over a union of rectangles.
   Each rectangle is a Box in a DisjointBoxLayout.
   There are two communication patterns associated with this data.
   1. copyTo: where one DistributedData writes its data to another 
   over the intersection of their DisjointBoxLayouts.   Ghost data is *not* overwritten.  
   This is  slightly different behavior than standard Chombo3 LevelData.
   2.  exchange: ghost cell information where ghost data from one grid 
   is filled with valid data from a neighboring grid.
     
   DistributedData        --- holds data and manages communcation.
   CommunicationMeta      --- deals with the communication of meta data information
   (also holds the string buffers used to call MPI)
   BoxPattern             --- Manages all the on-processor meta data 
                              (who is talking to whom via what boxes)
   box_interaction_t  the meta data associated with two boxes communciating
   proc_interaction_t the meta data associated with two procs communicating
*/

#include "Chombo_BoxPattern.H"
#include "Chombo_CommunicationMetaData.H"
#include "Chombo_DistributedData.H"
  

#endif

