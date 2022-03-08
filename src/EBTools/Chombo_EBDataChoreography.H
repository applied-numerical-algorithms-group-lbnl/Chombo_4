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
#include "Chombo_NeighborIterator.H"
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
   (DistributedData/MinimalCopier and subclasses) is meant to provide a
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

   Messsages are sent one box at at time.    We do not aggregate the communication for a whole processor
   the way LevelData does.   We have two phase communication where we typically only have a few (at most) boxes 
   per box.   This application does not merit the huge increase in complexity that the bit packing
   algorithm in LevelData entails.

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
    Chombo4::Box       m_region;
    int                m_procID;
    Chombo4::DataIndex m_datInd;
    void 
    define(const Chombo4::Box      & a_region,
           const int               & a_procID,
           const Chombo4::DataIndex& a_datind)
    {
      m_datInd = a_datind;
      m_region = a_region;
      m_procID = a_procID;
    }
  };
  
  ///
  /**
     motion_t
     This holds the meta data for both sides of a communication substep
  */
  class motion_t
  {
  public:
    motion_t() {;}
      
    //region, procID
    boxinfo_t m_src;
    boxinfo_t m_dst;
#ifdef CH_MPI    
    MPI_Request m_request;
#endif
    
    motion_t(const Chombo4::DataIndex& a_srcInd,
             const Chombo4::DataIndex& a_dstInd,
             const Chombo4::Box      & a_srcRegion,
             const Chombo4::Box      & a_dstRegion,
             const int               & a_srcProcID,
             const int               & a_dstProcID)
    {
      m_src.define(a_srcRegion, a_srcProcID, a_srcInd);
      m_dst.define(a_dstRegion, a_dstProcID, a_dstInd);
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
    MinimalCopier(const Chombo4::DisjointBoxLayout & a_grids, const IntVect& a_ghost)
    {
      CH_TIME("MinimalCopier::exchange constructor");

      Chombo4::DataIterator dit = a_grids.dataIterator();
      Chombo4::NeighborIterator nit(a_grids);
      for (dit.begin(); dit.ok(); ++dit)
      {
        const Chombo4::Box& grid = a_grids[dit];
        int myProcID = CH4_SPMD::procID();
        Chombo4::Box gridGhost(grid);
        gridGhost.grow(a_ghost);

        for (nit.begin(dit()); nit.ok(); ++nit)
        {
          Chombo4::Box neighbor = nit.box();
          int neiProcID = a_grids.procID(nit());
          Chombo4::Box neighborGhost= neighbor;
          neighborGhost.grow(a_ghost);

          Chombo4::Box gridGhostInter(neighbor & gridGhost);
          if (!gridGhostInter.isEmpty())
          {
            //his data copies to my ghost
            motion_t item(Chombo4::DataIndex(nit()), dit(), gridGhostInter, gridGhostInter, neiProcID, myProcID);
            if (myProcID == neiProcID)
            { // local move
              m_localMotionPlan.push_back(item);
            }
            else
            {
              m_toMotionPlan.push_back(item);
            }
          }
          Chombo4::Box neighborGhostInter(grid & neighborGhost);
          if (!neighborGhostInter.isEmpty())
          {
            //my data copies to his ghost
            //local case will be taken care of when the data iterator gets to the neighbor box
            if(myProcID != neiProcID)
            {
              motion_t item(dit(), Chombo4::DataIndex(nit()), neighborGhostInter, neighborGhostInter, myProcID, neiProcID);
              m_fromMotionPlan.push_back(item);
            }
          }
        }
      }
    }

    ///Define for copyTo --- does NOT include ghost cells.
    MinimalCopier(const Chombo4::DisjointBoxLayout & a_src,
                  const Chombo4::DisjointBoxLayout & a_dst)
    {
      CH_TIME("MinimalCopier::exchange constructor");

      LayoutIterator litSrc = a_src.layoutIterator();
      LayoutIterator litDst = a_dst.layoutIterator();
      int myProcID = CH4_SPMD::procID();

      ///first local
      for(litSrc.begin(); litSrc.ok(); ++litSrc)
      {
        int srcProcID   = a_src.procID(litSrc());
        Chombo4::Box gridSrc   = a_src[litSrc()];
        for(litDst.begin(); litDst.ok(); ++litDst)
        {
          int dstProcID = a_dst.procID(litDst());
          Chombo4::Box gridDst = a_src[litDst()];
          Chombo4::Box intersect(gridSrc & gridDst);
          if(!intersect.isEmpty())
          {
            motion_t item(Chombo4::DataIndex(litSrc()), Chombo4::DataIndex(litDst()), intersect, intersect, srcProcID, dstProcID);
            //if neither processor is == myprocID, nothing to do
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
                Chombo4::MayDay::Error("apparently I missed a case");
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


  class bufentry_t
  {
  public:
    
    bufentry_t()
    {
      m_buff = NULL;
      m_size = 0;
    }
    char*    m_buff;
    size_t   m_size;
    motion_t m_plan;
  };

  ///
  template<class fabtype_t>
  class DistributedData
  {
  public:
    ///datafactory has an ncomps argument that shall be ignored here
    DistributedData(const Chombo4::DisjointBoxLayout     & a_grids,
                    const Chombo4::IntVect               & a_ghost,
                    const Chombo4::DataFactory<fabtype_t>& a_factory)
    {
      CH_TIME("DistributedData constructor");
      m_grids = a_grids;
      m_ghost = a_ghost;
      Chombo4::DataIterator dit = a_grids.dataIterator();
      m_data.resize(dit.size(), NULL);
      for(int ibox = 0; ibox < dit.size(); ibox++)
      {
        Chombo4::Box ghosted = Chombo4::grow(a_grids[dit[ibox]], a_ghost);
        //the zero is for the fake ncomp.   components live in the type system now
        m_data[ibox] = a_factory.create(ghosted, 0, dit());
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

    fabtype_t& operator[](const Chombo4::DataIndex& a_dit)
    {
      return *(m_data[a_dit.datInd()]);
    }

    const fabtype_t& operator[](const Chombo4::DataIndex& a_dit) const
    {
      return *(m_data[a_dit.datInd()]);
    }

    ///
    void exchange(bool a_printStuff = false)
    {
      CH_TIME("DistributedData::exchange");
      MinimalCopier copier(m_grids, m_ghost);
      communicate(*this, copier, *this,  a_printStuff);
    }
    
    ///
    void copyTo(DistributedData & a_dst,
                bool a_printStuff = false)
    {
      CH_TIME("DistributedData::copyTo");
      MinimalCopier copier(m_grids, a_dst.m_grids);
      communicate(*this, copier, a_dst,  a_printStuff);
    }

  private:
    //get sizes and stuff and allocate buffers
    static void  communicateMetaData(MinimalCopier      & a_pattern,
                                     vector<bufentry_t> & a_fromMe,
                                     vector<bufentry_t> & a_toMe,
                                     vector<fabtype_t*> & a_data,
                                     bool a_printStuff = false)
    {
      CH_TIME("DistributedData::communicateMetaData");
#ifdef CH_MPI  
      
      auto comm = CH4_SPMD::Chombo_MPI::comm;

      //post sizes
      for(int ibuf = 0; ibuf < a_fromMe.size(); ibuf++)
      {
        //not constant because I am using the request
        auto& motion  = a_pattern.m_fromMotionPlan[ibuf];
        CH_assert(motion.m_src.m_procID == CH4_SPMD::procID());
        //components are dummy arguments
        const fabtype_t& localfab = *(a_data[motion.m_src.m_datInd.datInd()]);
        unsigned long bufsize = localfab.charsize(motion.m_src.m_region, 0, 1);
        a_fromMe[ibuf].m_size = bufsize;
        a_fromMe[ibuf].m_buff = (char*)(malloc(bufsize));
        auto dstprocID = motion.m_dst.m_procID;
        auto srcprocID = motion.m_src.m_procID;
        auto& request  = motion.m_request;

        if(a_printStuff)
        {
          Chombo4::pout() << "proc = "     << srcprocID
                          << "\t will send   a " << bufsize
                          << "\t sized buffer to   proc " << dstprocID
                          << "\t for region = " << motion.m_src.m_region << endl;
        }
        MPI_Isend(&bufsize, 1, MPI_UNSIGNED_LONG, dstprocID, 0, comm, &request);

      }


      //get the sizes back
      for(int ibuf = 0; ibuf < a_toMe.size(); ibuf++)
      {
        auto& motion = a_pattern.m_toMotionPlan[ibuf];
        CH_assert(motion.m_dst.m_procID == CH4_SPMD::procID());
        //components are dummy arguments
        unsigned long bufsize = 0;
        auto dstprocID = motion.m_dst.m_procID;
        auto srcprocID = motion.m_src.m_procID;
        auto& request  = motion.m_request;

        MPI_Irecv(&bufsize, 1, MPI_UNSIGNED_LONG, srcprocID, MPI_ANY_TAG, comm, &request);
        MPI_Status   status;
        MPI_Wait(&request, &status);
        if(a_printStuff)
        {
          Chombo4::pout() << "proc = "     << dstprocID
                          << "\t will receive a " << bufsize
                          << "\t sized buffer from proc " << srcprocID
                          << "\t for region = " << motion.m_src.m_region << endl;
        }
        a_toMe[ibuf].m_size = bufsize;
        a_toMe[ibuf].m_buff = (char*)(malloc(bufsize));
      }
      
      if(a_printStuff)
      {
        Chombo4::pout() << "DistributedData::communicateMetaData(): just before barrier" << endl;
      }
      //Do not proceed until everyone has their receives.
      //This simplifies matters but should not technically be necessary.
      MPI_Barrier(comm);
      if(a_printStuff)
      {
        Chombo4::pout() << "DistributedData::communicateMetaData(): just after  barrier" << endl;
      }
#endif      
    }

    //actually serialize, send, and receive data
    static void  doRemoteCopies(MinimalCopier      & a_pattern,
                                vector<bufentry_t> & a_fromMe,
                                vector<bufentry_t> & a_toMe,
                                vector<fabtype_t*> & a_data,
                                bool a_printStuff = false)
    {
      CH_TIME("DistributedData::doRemoteCopies");
#ifdef CH_MPI
      if(a_printStuff)
      {
        Chombo4::pout() << "entering DistributedData::doRemoteCopies"  << std::endl;
      }
      auto comm = CH4_SPMD::Chombo_MPI::comm;

      //fill buffers and send them
      for(int ibuf = 0; ibuf < a_fromMe.size(); ibuf++)
      {
        //not constant because I am using the request
        auto& motion  = a_pattern.m_fromMotionPlan[ibuf];
        CH_assert(motion.m_src.m_procID == CH4_SPMD::procID());
        //components are dummy arguments
        const fabtype_t& localfab = *(a_data[motion.m_src.m_datInd.datInd()]);
        auto dstprocID = motion.m_dst.m_procID;
        //auto srcprocID = motion.m_src.m_procID;
        auto& request  = motion.m_request;
        auto  bufsize  = a_fromMe[ibuf].m_size;
        auto* charbuf  = a_fromMe[ibuf].m_buff;
        //encode data into string buffer
        localfab.linearOut(charbuf, motion.m_src.m_region, 0, 1);
        //send string buffer to other proc
        MPI_Isend(charbuf,  bufsize, MPI_BYTE, dstprocID, 0, comm, &request);

      }


      //receive buffers
      for(int ibuf = 0; ibuf < a_toMe.size(); ibuf++)
      {
        auto& motion = a_pattern.m_toMotionPlan[ibuf];
        CH_assert(motion.m_dst.m_procID == CH4_SPMD::procID());
        //components are dummy arguments
        //auto dstprocID = motion.m_dst.m_procID;
        auto srcprocID = motion.m_src.m_procID;
        auto& request  = motion.m_request;
        auto* charbuf = a_toMe[ibuf].m_buff;
        auto  bufsize = a_toMe[ibuf].m_size;

        //ask MPI fill to the buffer 
        MPI_Irecv(charbuf, bufsize, MPI_BYTE, srcprocID, MPI_ANY_TAG, comm, &request);

        //wait until it actually arrives
        MPI_Status   status;
        MPI_Wait(&request, &status);

        //now unwind the buffer into the data
        fabtype_t& localfab = *(a_data[motion.m_dst.m_datInd.datInd()]);
        localfab.linearIn(charbuf, motion.m_dst.m_region, 0, 1);
      }
      
      if(a_printStuff)
      {
        Chombo4::pout() << "DistributedData::communicateMetaData(): just before barrier" << endl;
      }
      //Do not proceed until everyone has their receives.
      //This simplifies matters but should not technically be necessary.
      MPI_Barrier(comm);
      if(a_printStuff)
      {
        Chombo4::pout() << "leaving  DistributedData::doRemoteCopies" << std::endl;
      }
#endif      
    }

    static void freeBuffers(vector<bufentry_t> & a_fromMe,
                            vector<bufentry_t> & a_toMe,
                            bool a_printStuff = false)
    {
      CH_TIME("DistributedData::freeBuffers");
      for(int ibuf = 0; ibuf < a_fromMe.size(); ibuf++)
      {
        free(a_fromMe[ibuf].m_buff);
        a_fromMe[ibuf].m_buff= NULL;
      }
      for(int ibuf = 0; ibuf < a_toMe.size(); ibuf++)
      {
        free(a_toMe[ibuf].m_buff);
        a_toMe[ibuf].m_buff= NULL;
      }



    }
    
    static void doLocalCopies(DistributedData       & a_dst,
                              const DistributedData & a_src,
                              const MinimalCopier   & a_pattern,
                              bool a_printStuff = false)
    {
      CH_TIME("DistributedData::doLocalCopies");
      for(int ibuf = 0; ibuf < a_pattern.m_localMotionPlan.size(); ibuf++)
      {
        //component arguments are artifacts of a more civilized age
        //components are all in the type system now
        const auto& motion = a_pattern.m_localMotionPlan[ibuf];
        const auto& srcdatind = motion.m_src.m_datInd;
        const auto& srcregion = motion.m_src.m_region;
        const auto& dstdatind = motion.m_dst.m_datInd;
        const auto& dstregion = motion.m_dst.m_region;
        CH_assert(srcregion == dstregion);
        Bx regx = ProtoCh::getProtoBox(dstregion);
        unsigned int ico = 0;
        unsigned int ncodst = a_dst[dstdatind].nComp();
        unsigned int ncosrc = a_src[srcdatind].nComp();
        CH_assert(ncodst == ncosrc);
        a_dst[dstdatind].copy(a_src[srcdatind], regx, ico, regx, ico, ncodst);
      }
    }

    
    static void communicate(DistributedData       & a_dst,
                            MinimalCopier         & a_pattern,
                            const DistributedData & a_src,
                            bool a_printStuff = false)
    {
      doLocalCopies(a_dst, a_src, a_pattern, a_printStuff);
#ifdef CH_MPI  
      vector<bufentry_t> fromMe(a_pattern.m_fromMotionPlan.size());
      vector<bufentry_t> toMe  (a_pattern.m_toMotionPlan.size());
      //communication of size information and allocate buffers
      communicateMetaData(a_pattern, fromMe, toMe, a_dst.m_data, a_printStuff);
      //communication of actual data
      doRemoteCopies(a_pattern, fromMe, toMe, a_dst.m_data, a_printStuff);
      //frees up the buffers allocated in communicateMetaData
      freeBuffers(fromMe, toMe, a_printStuff);
#endif      
    }
    //all buffer data is transient.
    
    vector<fabtype_t*> m_data;
    IntVect m_ghost;

    Chombo4::DisjointBoxLayout m_grids;
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

