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
   We aggregate all communications between two processors the same way LevelData does.  

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
     box_interaction_t
     This holds the meta data for both sides of a communication substep
  */
  class box_interaction_t
  {
  public:
    box_interaction_t() {;}
      
    //region, procID
    boxinfo_t m_src;
    boxinfo_t m_dst;
    
    box_interaction_t(const Chombo4::DataIndex& a_srcInd,
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

      m_toMotionPlan.resize(  CH4_SPMD::numProc());
      m_fromMotionPlan.resize(CH4_SPMD::numProc());
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
            box_interaction_t item(Chombo4::DataIndex(nit()), dit(), gridGhostInter, gridGhostInter, neiProcID, myProcID);
            if (myProcID == neiProcID)
            { // local move
              m_localMotionPlan.push_back(item);
            }
            else
            {
              m_toMotionPlan[neiProcID].push_back(item);
            }
          }
          Chombo4::Box neighborGhostInter(grid & neighborGhost);
          if (!neighborGhostInter.isEmpty())
          {
            //my data copies to his ghost
            //local case will be taken care of when the data iterator gets to the neighbor box
            if(myProcID != neiProcID)
            {
              box_interaction_t item(dit(), Chombo4::DataIndex(nit()), neighborGhostInter, neighborGhostInter, myProcID, neiProcID);
              m_fromMotionPlan[neighProcID].push_back(item);
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

      m_toMotionPlan.resize(  CH4_SPMD::numProc());
      m_fromMotionPlan.resize(CH4_SPMD::numProc());
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
            box_interaction_t item(Chombo4::DataIndex(litSrc()), Chombo4::DataIndex(litDst()), intersect, intersect, srcProcID, dstProcID);
            //if neither processor is == myprocID, nothing to do
            if((srcProcID == myProcID) || (dstProcID == myProcID))
            {
              if((srcProcID == myProcID) && (dstProcID == myProcID))
              {
                m_localMotionPlan.push_back(item);
              }
              else if((srcProcID == myProcID) && (dstProcID != myProcID))
              {
                m_fromMotionPlan[dstProcID].push_back(item);
              }
              else if((srcProcID != myProcID) && (dstProcID == myProcID))
              {
                m_toMotionPlan[srcProcID].push_back(item);
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
    vector<box_interaction_t> m_localMotionPlan;

    //plan where source != my proc, dest == my proc
    /// outer loop is per proc 
    vector< vector<box_interaction_t> >m_toMotionPlan;

    //plan where source == my proc, dest != my proc
    /// outer loop is per proc 
    vector< vector<box_interaction_t> >m_fromMotionPlan;

    
  private:
    ///I like strong construction.
    MinimalCopier();
    MinimalCopier( const MinimalCopier& a_input);
    void operator=(const MinimalCopier& a_input);
  };

  /// 
  class proc_interaction_t
  {
  public:
    proc_interaction_t()
    {
      m_srcProcID = -1;
      m_dstProcID = -1;
      m_messageLen   = 0;
      m_metaDataLen  = 0;
      m_buffer       = NULL;
    }

    int                       m_srcProcID;
    int                       m_dstProcID;
    unsigned long             m_messageLen;
    unsigned long             m_metaDataLen;
    char*                     m_buffer;
    vector<unsigned long>     m_boxbufsize;

#ifdef CH_MPI    
    MPI_Request m_dataRequest;
    MPI_Request m_metaRequest;
#endif
             
  };

  ///class to manage meta data for communication between two processors
  class CommunicationMetaData
  {
  public:
    
    CommunicationMetaData(const MinimalCopier& a_pattern,
                          vector<fabtype_t*> & a_data,
                          bool a_printStuff = false)
      
    {
      CH_TIME("CommunicationMetaData constructor");
#ifdef CH_MPI  
      auto comm = CH4_SPMD::Chombo_MPI::comm;

      //post sizes
      for(int iproc = 0; iproc < fromMePattern.size(); iproc++)
      {
        const auto& fromMePattern =  a_pattern.m_fromMe;
        if(fromMePattern[iproc].size() > 0)
        {
          vector<unsigned long> boxbufsize(fromMePattern[iproc].size());

          unsigned long totbufsize = 0;
          for(int ibuf = 0; ibuf < fromMePattern[iproc].size(); ibuf++)
          {
            //not constant because I am using the request
            const auto& motion  = fromMePattern[iproc][ibuf];
            CH_assert(motion.m_src.m_procID == CH4_SPMD::procID());
            const fabtype_t& localfab = *(a_data[motion.m_src.m_datInd.datInd()]);
            unsigned long thisbufsize = localfab.charsize(motion.m_src.m_region, 0, 1);
            boxbufsize[ibuf] = thisbufsize;
            totbufsize += thisbufsize;
          }
          
          shared_ptr<proc_interaction_t> from_pi(new proc_interaction_t());
          from_pi->m_srcProcID  = CH4_SPMD::procID();
          from_pi->m_dstProcID  = iproc;
          from_pi->m_messageLen = totbufsize;
          from_pi->m_buffer     = malloc(totbufsize);
          from_pi->m_metaDataLen= (boxbufsize.size())*sizeof(unsigned long);
          from_pi->m_boxbufsize  = boxbufsize;
          
          m_fromMe.push_back(from_pi);
          if(a_printStuff)
          {
            Chombo4::pout() << "proc = "                   << from_pi.m_srcprocID
                            << "\t will send  " <<  fromMePattern[iproc].size()
                            << "\t boxes with of data amounting to a " <<  totbufsize 
                            << "-sized \t buffer to  proc " << from_pi.m_dstprocID << endl;
          }
          char* sendbuf = (char*)(boxbufsize.data());
          MPI_Isend(sendbuf, from_pi->m_metaDataLen, MPI_BYTE, from_pi.dstprocID, 0, comm, &from_pi.m_metaRequest);
        } //if there is anything to send to this proc
      }   //loop over procs


      for(int iproc = 0; iproc < toMePattern.size(); proc++)
      {
        const auto&   toMePattern =  a_pattern.m_toMe;
        if(toMePattern[iproc].size() > 0)
        {
          vector<unsigned long> boxbufsize(toMePattern[iproc].size());
          shared_ptr<proc_interaction_t> to_pi(new proc_interaction_t());
          to_pi->m_srcProcID   = iproc;
          to_pi->m_dstProcID   = CH4_SPMD::procID();
          to_pi->m_metaDataLen = (boxbufsize.size())*sizeof(unsigned long);
          char* sendbuf = (char*)(boxbufsize.data());
          MPI_Irecv(sendbuf, to_pi->m_metaDataLen, MPI_BYTE, to_pi.srcprocID, MPI_ANY_TAG, comm, &to_pi.m_metaRequest);
          MPI_Status   status;
          MPI_Wait(&to_pi.m_metaRequest, &status);
          
          unsigned long totbufsize = 0;
          for(int ibuf = 0; ibuf < toMePattern[iproc].size(); ibuf++)
          {
            totbufsize += boxbufsize[ibuf];
          }
          to_pi->messageLen = totbufsize;
          to_pi->m_buffer   = malloc(totbufsize);

          if(a_printStuff)
          {
            Chombo4::pout() << "proc = "     << dstprocID
                            << "\t will receive " << toMePattern[iproc].size()
                            << "grids worth of the dat for a total of a " << totbufsize
                            << "-sized buffer from proc " << srcprocID << endl;
          }
      }
      
    }

    //these can be any length up to numProc()
    vector<shared_ptr<proc_interaction_t> > m_toMe;
    vector<shared_ptr<proc_interaction_t> > m_fromMe;
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

    //actually serialize, send, and receive data
    static void  doRemoteCopies(CommunicationMetaData & a_metadataco,
                                vector<fabtype_t*>    & a_data,
                                bool a_printStuff = false)
    {
      CH_TIME("DistributedData::doRemoteCopies");
#if 0
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
      //communication of size information and allocate buffers
      CommunicationMetaData metadataco(a_dst.m_data, pattern, a_printStuff);

      //communication of actual data
      doRemoteCopies(a_pattern, fromMe, toMe, a_dst.m_data, a_printStuff);
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

