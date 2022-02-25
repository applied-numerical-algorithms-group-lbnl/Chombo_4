#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _Chombo_Minimal_COPIER_H_
#define _Chombo_Minimal_COPIER_H_

#include "Chombo_DisjointBoxLayout.H"
#include "Chombo_Pool.H"
#include "Chombo_Vector.H"
#include "Chombo_ProblemDomain.H"
#include <unordered_map>
#include <cstdint>

#include "Chombo_SPMD.H"
using namespace Chombo4;
///
/**
   The Ch4_Data_Choreography space
   (HostLevelData/SPMDBuffer/MinimalCopier) is meant to provide a
   reduced complexity, highly maintainable alternative to Copier and
   LevelData infrastructure.  If you need high performance or need
   some of the fancier aspects of Copier, you should use the
   standard LevelData/BoxLayoutData/LayoutData/Copier
   infrastructure.

   Periodic boundary conditions are not supported.     

   All communication is two phase.

   HostLevelData holds data over a union of rectangles.
   Each rectangle is a Box in a DisjointBoxLayout.
   There are two communication patterns associated with this data.
   1. copyTo: where one HostLevelData writes its data to another 
   over the intersection of their DisjointBoxLayouts.   Ghost data is *not* overwritten.  
   This is  slightly different behavior than standard Chombo3 LevelData.
   2.  exchange: ghost cell information where (within a HostData) ghost data from one grid 
   is filled with valid data from a neighboring grid.
     
   HostLevelData --- holds data and manages communcation.
   SPMDBuffer    --- manages buffers associated with two phase communation
   MinimalCopier --- manages the meta data associated with communcation.
*/
namespace CH4_Data_Choreography
{


  ///
  /**
    boxinfo_t
    This is the meta data for one side of a message where the data for one grid 
    is copying over a subset of the data in another grid.
    dblbox --- box of valid data.
    region --- box for copying is a subset of grow(dblbox, nghost).   Can be entirely outside dblbox.
    procid --- process id associated with this side of the message
    datind --- data index associated with this side of the message
   */
  struct boxinfo_t
  {
    Box       m_dblbox;
    Box       m_region;
    int       m_procid;
    DataIndex m_datind;
  };
  
  ///
  /**
     motion_t
     This holds the meta data for both sides of a communication substep
   */
  struct motion_t;
  {
    //region, procid
    boxinfo_t m_src;
    boxinfo_t m_dst;
  };

  ///
  /**
     MinimalCopier class to reduce the complexity of copier and leveldata and make them more maintainable.
     If you need high performance or need some of the fancier aspects of Copier, you should use the standard
     LevelData/BoxLayoutData/LayoutData/Copier infrastructure.  Periodic boundary conditions are not supported.
     All communication is two phase.
  */
  class MinimalCopier
  {
    ///
    /**
       Define for exchange--fills ghost cells around boxes if neighboring boxes are adjacent
    */
    MinimalCopier(const DisjointBoxLayout & a_grids, const IntVect& a_ghost)
    {

    }

    ///Define for copyTo --- does NOT include ghost cells.
    MinimalCopier(const DisjointBoxLayout & a_src,
                  const DisjointBoxLayout & a_dst)
    {
    }
  protected:

    vector<motion_t> m_localMotionPlan;
    vector<motion_t> m_fromMotionPlan;
    vector<motion_t> m_toMotionPlan;

    
  private:
    //Any of this nonsense is asking for trouble with
    //the Cavalier way we are dealing with pointered data.
    //of course, the Roundhead way is worse.
    MinimalCopier();
    MinimalCopier( const MinimalCopier& a_input);
    void operator=(const MinimalCopier& a_input);
  };

}
#endif
