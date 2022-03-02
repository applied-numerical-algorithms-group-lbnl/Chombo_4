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
   (DistributedData/MinimalCopier and subclasses) is meant to provide a
   reduced complexity, highly maintainable alternative to the Copier and
   LevelData infrastructure.  If you need high performance or need
   some of the fancier aspects of Copier, you should use the
   standard LevelData/BoxLayoutData/LayoutData/Copier
   infrastructure.

   Periodic boundary conditions are not supported.     

   All communication is two phase.
   All communication is done on the host.

   DistributedData holds data over a union of rectangles.
   Each rectangle is a Box in a DisjointBoxLayout.
   There are two communication patterns associated with this data.
   1. copyTo: where one HostLevelData writes its data to another 
   over the intersection of their DisjointBoxLayouts.   Ghost data is *not* overwritten.  
   This is  slightly different behavior than standard Chombo3 LevelData.
   2.  exchange: ghost cell information where (within a HostData) ghost data from one grid 
   is filled with valid data from a neighboring grid.
     
   DistributedData        --- holds data and manages communcation.
   CommunicationBuffer    --- manages transient buffers associated with two phase communation
   MinimalCopier          --- manages the meta data associated with communcation.
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
    Box       m_region;
    int       m_procid;
    DataIndex m_datind;
    void 
    define(const Box      & a_region,
           const int      & a_procID,
           const DataIndex& a_datind)
    {
      m_datind = a_datind;
      m_region = a_region;
      m_procID = a_procID;
    }
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
    motion_t(const DataIndex& a_srcInd,
             const DataIndex& a_dstInd,
             const int      & a_srcProcID,
             const int      & a_dstProcID,
             const Box      & a_srcRegion,
             const Box      & a_dstRegion)
    {
      m_src.define(a_srcRegion, a_srcInd, a_srcProcID);
      m_dst.define(a_dstRegion, a_dstInd, a_dstProcID);
    }
             
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
  public:
    ///
    /**
       Define for exchange--fills ghost cells around boxes if neighboring boxes are adjacent
    */
    MinimalCopier(const DisjointBoxLayout & a_grids, const IntVect& a_ghost)
    {
      CH_TIME("MinimalCopier::exchange constructor");

      DataIterator dit = a_grids.dataIterator();
      NeighborIterator nit(a_grids);
      for (dit.begin(); dit.ok(); ++dit)
      {
        const Box& grid = a_grids[dit];
        int gridProcID = procID();
        Box gridGhost(grid);
        gridGhost.grow(a_ghost);

        for (nit.begin(dit()); nit.ok(); ++nit)
        {
          Box neighbor = nit.box();
          int neiProcID = a_grids.procID(nit());
          Box neighborGhost= neighbor;
          neighborGhost.grow(a_ghost);

          Box gridGhostInter(neighbor & gridGhost);
          if (!gridGhostInter.isEmpty())
          {
            //his data copies to my ghost
            motion_t item(DataIndex(nit()), dit(), gridGhostInter, gridGhostInter, neiProcID, myProcID);
            if (srcProcID == dstProcID)
            { // local move
              m_localMotionPlan.push_back(item);
            }
            else
            {
              m_toMotionPlan.push_back(item);
            }
          }
          Bx neighborGhostInter(grid & neighborGhost);
          if (!neighborGhostInter.isEmpty())
          {
            //my data copies to his ghost
            //local case will be taken care of when the data iterator gets to the neighbor box
            If(srcProcID != dstProcID)
            {
              motion_t item(DataIndex(dit(), nit()), neighborGhostInter, neighborGhostInter, myProcID, neiProcID);
              m_fromMotionPlan.push_back(item);
            }
          }
        }
      }
    }

    ///Define for copyTo --- does NOT include ghost cells.
    MinimalCopier(const DisjointBoxLayout & a_src,
                  const DisjointBoxLayout & a_dst)
    {
      CH_TIME("MinimalCopier::exchange constructor");

      LayoutIterator litSrc = a_src.layoutIterator();
      LayoutIterator litDst = a_dst.layoutIterator();
      int myProcID = procID();

      ///first local
      for(litSrc.begin(); litSrc.ok(); ++litSrc)
      {
        int srcProcID = a_src.procID(litSrc());
        Box gridSrc   = a_src[litSrc()];
        for(litDst.begin(); litDst.ok(); ++litDst)
        {
          int dstProcID = a_dst.procID(litDst());
          Box gridDst = a_src[dit()];
          Box intersect(gridSrc & gridDst);
          if(!intersect.isEmpty())
          {
            motion_t item(DataIndex(litSrc()), DataIndex(litDst()), intersect, intersect, srcProcID, dstProcID);
            //if neither processor is == myprocid, nothing to do
            if((srcProcID == myProcID) || (dstProcID == myProcID))
            {
              if((srcProcID == myProcID) && (dstProcID == myProcID))
              {
                m_localMotionPlan.push_back(item);
              }
              else if((srcProcID == myProcID) && (dstProcID != myProcID))
              {
                m_fromMotionPlan.push_back(item);
              }
              else if((srcProcID != myProcID) && (dstProcID == myProcID))
              {
                m_toMotionPlan.push_back(item);
              }
              else
              {
                MayDay::Error("apparently I missed a case");
              }
            } //if(something is on proc)
          } //if there is an intersection
        } //loop over destination boxes
      } //loop over source boxes
    }


    //plan where source == my proc,  dest == my proc
    vector<motion_t> m_localMotionPlan;
    //plan where source != my proc, dest == my proc
    vector<motion_t> m_toMotionPlan;
    //plan where source == my proc, dest != my proc
    vector<motion_t> m_fromMotionPlan;

    
  private:
    //Any of this nonsense is asking for trouble with
    //the Cavalier way we are dealing with pointered data.
    //Of course, the Roundhead way is worse.
    MinimalCopier();
    MinimalCopier( const MinimalCopier& a_input);
    void operator=(const MinimalCopier& a_input);
  };


  struct bufentry_t
  {
    char*    m_buff;
    size_t   m_size;
    motion_t m_plan;
  };
  
  ///
  template<class fabtype_t>
  class DistributedData
  {

    ///datafactory has an ncomps argument that shall be ignored here
    DisributedData(const DisjointBoxLayout & a_grids, const IntVect& a_ghost, const DataFactory<fabtype_t>& a_factory)
    {
      m_grids = a_grids;
      DataIterator dit = a_grids.dataIterator();
      m_data.resize(dit.size(), NULL);
      for(int ibox = 0; ibox < dit.size(); ibox++)
      {
        Box ghosted = a_grids.grow(a_ghost);
        //the zero is for the fake ncomp.   components live in the type system now
        m_data[ibox] = a_factory.create(ghost, 0, dit());
      }
    }

    ~DistributedData()
    {
      for(int ibox = 0; ibox < m_data.size(); ibox++)
      {
        delete m_data[ibox];
        m_data[ibox] = NULL;
      }
    }

    fabtype_t& operator[](const DataIndex& a_dit)
    {
      return m_data(a_dit.datInd());
    }

    const fabtype_t& operator[](const DataIndex& a_dit) const
    {
      return m_data(a_dit.datInd());
    }

    ///
    void exchange(bool a_printStuff = false);
    {
      MinimalCopier copier(m_grids);
      communicate(*this, *this, copier, a_printStuff);
    }
    
    ///
    void copyTo(DistributedData & a_dst,
                bool a_printStuff = false)
    {
      MinimalCopier copier(m_grids, a_dst.m_grids);
      communicate(*this, a_dst, copier, a_printStuff);
    }

  protected:
    //get sizes and stuff
    void  fillBufferEntries(const MinimalCopier& a_pattern,
                            vector<bufentry_t> & a_fromMe,
                            vector<bufentry_t> & a_toMe,
                            bool a_printStuff = false)
    {
      vector<bufentry_t> fromMe(m_copier->m_fromMotionPlan.size());
      vector<bufentry_t> toMe  (m_copier->m_toMotionPlan.size());
      //post sizes
      for(ibuf = 0; ibuf < fromMe.size(); ibuf++)
      {
        const auto& motion = a_pattern.m_fromMotionPlan[ibuf];
        CH_assert(motion.m_src.m_procID == procID());
        //components are dummy arguments
        unsigned long bufsize = (m_data[motion.m_src.m_datInd]).charsize(motion.m_region, 0, 1);
        fromMe[ibuf].m_size = bufsize;
        MPI_SEND(&bufsize, 1, MPI_UNSIGNED_LONG, motion.m_dst.m_procID, motion.m_src.m_procID, comm);
      }

      //get the sizes back
      for(ibuf = 0; ibuf < toMe.size(); ibuf++)
      {
        const auto& motion = a_pattern.m_toMotionPlan[ibuf];
        CH_assert(motion.m_dst.m_procID == procID());
        //components are dummy arguments
        unsigned long bufsize = 0;
        MPI_Recv(&bufsize, 1, MPI_UNSIGNED_LONG, motion.m_dst.m_procID, motion.m_src.m_procID, comm);
        if(a_printStuff)
        {
          pout() << "proc=" << procID() << "will receive " << bufsize << "size buffer from " << motion.m_src.m_procID << endl;
        }
        toMe[ibuf].m_size = bufsize;
      }
      
      if(a_printStuff)
      {
        pout() << "DistributedData::fillBufferEntries(): just before barrier" << endl;
      }
      //do not proceed until everyone has their receives
      MPI_barrier(comm);
      if(a_printStuff)
      {
        pout() << "DistributedData::fillBufferEntries(): just after  barrier" << endl;
      }
    }
    static void communicate(DistributedData       & a_dst,
                            const DistributedData & a_src,
                            const MinimalCopier   & a_pattern,
                            bool a_printStuff = false)
    {
      vector<bufentry_t> fromMe, toMe;
      //includes communication of size information
      fillBufferEntries(a_pattern, fromMe, toMe, a_printStuff);
    }
    //all buffer data is transient.
    
    vector<fabtype_t*> m_data;


    DisjointBoxLayout m_grids;
  private:
    //strong construction brings simplicity.   I like simplicity.
    DistributedData();
    //Copy constuction and assignment are disallowed
    //because this simplifies memory management. I like simplicity.
    DistributedData(const DistributedData& a_input);
    void operator=(const DistributedData& a_input);
      
    
  };
  
}
#endif

