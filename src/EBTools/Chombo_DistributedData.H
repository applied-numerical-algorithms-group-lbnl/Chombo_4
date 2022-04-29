#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _Distributed_Data_H__
#define _Distributed_Data_H__

#include "Chombo_DisjointBoxLayout.H"
#include "Chombo_Pool.H"
#include <vector>
#include "Chombo_ProblemDomain.H"
#include "Chombo_CH_Timer.H"
#include "Chombo_NeighborIterator.H"
#include "Chombo_DataIterator.H"
//for datafactory
#include "Chombo_BoxLayoutData.H"
#include "Chombo_BoxPattern.H"
#include "Chombo_CommunicationMetaData.H"
#include <unordered_map>
#include <cstdint>

#include "Chombo_SPMD.H"
#ifdef CH_MPI  
#include "mpi.h"
#endif
using namespace Chombo4;
namespace CH4_Data_Choreography
{


  ///
  /**
     Data over a union of rectangles.   
     CopyTo only copies to valid data (no ghost).
     Exchange does the usual thing (copy from valid data in neighboring cells 
     to overlapping ghost cells).
     Periodic boundary conditions are not supported.
     All communication is two phase so exchange and copyTo can involve as many as 4*N^2 messages
     where N is the number of processors.    
     At each stage, communication between two processors is agreggated into one message.
   */
  template<class fabtype_t>
  class DistributedData
  {
  public:
    ///
    /**
       Strong construction brings simplicity.   I like simplicity.  Sometimes
       I have to use weak construction, anyway to conform to old APIs.
       Let this be a lesson to you.
     */
    DistributedData()
    {
      m_isDefined = false;
    }
    
    ///
    /**
       DataFactory has an ncomps argument that shall be ignored here.  
       It is a historical artifact that future generations will find fascinating.
    */
    DistributedData(const Chombo4::DisjointBoxLayout     & a_grids,
                    const Chombo4::IntVect               & a_ghost,
                    const Chombo4::DataFactory<fabtype_t>& a_factory)
    {
      define(a_grids, a_ghost, a_factory);
    }

    
    ///datafactory has an ncomps argument that shall be ignored here
    void define(const Chombo4::DisjointBoxLayout     & a_grids,
                const Chombo4::IntVect               & a_ghost,
                const Chombo4::DataFactory<fabtype_t>& a_factory)
    {
      CH_TIME("DistributedData constructor");
      m_procID = CH4_SPMD::procID();
      m_isDefined = true;
      m_grids = a_grids;
      m_ghost = a_ghost;
      m_graphConstructed= a_factory.graphConstructor();
      
      Chombo4::DataIterator dit = a_grids.dataIterator();
      m_data.resize(dit.size(), NULL);
      
      for(int ibox = 0; ibox < dit.size(); ibox++)
      {
        Chombo4::Box ghosted = Chombo4::grow(a_grids[dit[ibox]], a_ghost);
        //The fake ncomp (1) is not so fake if you use this on old data holders.
        m_data[ibox] = a_factory.create(ghosted, 1, dit[ibox]);
      }
    }
    
    ///old leveldata functionality that got used a LOT
    Chombo4::DataIterator dataIterator() const
    {
      return m_grids.dataIterator();
    }
    
    virtual ~DistributedData()
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
      BoxPattern pattern(m_grids, m_ghost, a_printStuff);
      communicate(*this, pattern, *this  , a_printStuff, m_graphConstructed);
    }
    
    ///
    void copyTo(DistributedData & a_dst,
                bool a_printStuff = false)
    {
      CH_TIME("DistributedData::copyTo");
      BoxPattern pattern(m_grids, a_dst.m_grids, a_printStuff);
      communicate(*this, pattern, a_dst, a_printStuff, m_graphConstructed);
    }

    ///
    IntVect ghostVect() const
    {
      return m_ghost;
    }

    ///
    Chombo4::DisjointBoxLayout disjointBoxLayout() const
    {
      return m_grids;
    }
    
  private:

    static void communicate(DistributedData       & a_dst,
                            BoxPattern            & a_pattern,
                            const DistributedData & a_src,
                            bool a_printStuff,
                            bool a_graphConstructed)
    {
      CH_TIME("DistributedData::communicate");
#ifdef CH_MPI  
      //communication of size information and allocate buffers
      CommunicationMetaData<fabtype_t> metadataco(a_pattern, a_dst.m_data, a_printStuff);
      postReceives(   a_dst.m_data, metadataco, a_pattern, a_printStuff, a_graphConstructed);
      postSends(      a_dst.m_data, metadataco, a_pattern, a_printStuff, a_graphConstructed);
#endif      
      doLocalCopies(a_dst, a_src, a_pattern, a_printStuff);
#ifdef CH_MPI  
      completePending(a_dst.m_data, metadataco, a_pattern, a_printStuff, a_graphConstructed);
      unpackReceives (a_dst.m_data, metadataco, a_pattern, a_printStuff, a_graphConstructed);
#endif      
    }

#ifdef CH_MPI    
    ///mpi_irecv 
    static void  postReceives(vector<fabtype_t*>               & a_data,
                              CommunicationMetaData<fabtype_t> & a_metacomm,
                              const BoxPattern                 & a_pattern,
                              bool a_printStuff,
                              bool a_graphConstructed)
    {

      CH_TIME("DistributedData::postReceives");
      if(a_printStuff)
      {
        Chombo4::pout() << "DistributedData::postReceives begin" << std::endl;
      }
      for(int iproc = 0; iproc < CH4_SPMD::numProc() ; iproc++)
      {
        auto recvPattern = a_pattern.m_recvPattern[iproc];
        if(recvPattern.size() > 0)
        {
          auto & recv_pi   =  a_metacomm.m_recvPI[iproc];
          const auto&  dstProc     =  recv_pi->m_dstProcID;
          char* charbuf            =  recv_pi->m_buffer;
          const auto&  srcProc     =  recv_pi->m_srcProcID;
          const auto&  totbufsize  =  recv_pi->m_messageLen;

          CH_assert(dstProc == CH4_SPMD::procID());
          
          //ask MPI fill to the buffer 
          auto comm = CH4_SPMD::Chombo_MPI::comm;
          if(a_printStuff)
          {
            Chombo4::pout() << "DistributedData: about to receive " << totbufsize << " from proc " << srcProc << endl;
          }
          MPI_Irecv(charbuf, totbufsize, MPI_BYTE, srcProc, MPI_ANY_TAG, comm,
                    &(recv_pi->m_dataRequest));
        } 
      }
      if(a_printStuff)
      {
        Chombo4::pout() << "DistributedData::postReceives end" << std::endl;
      }
    }
    
    //linearin called here
    static void  unpackReceives(vector<fabtype_t*>               & a_data,
                                CommunicationMetaData<fabtype_t> & a_metacomm,
                                const BoxPattern                 & a_pattern,
                                bool a_printStuff,
                                bool a_graphConstructed)
    {
      CH_TIME("DistributedData::unpackReceives");
      if(a_printStuff)
      {
        Chombo4::pout() << "DistributedData::unpackReceives begin" << std::endl;
      }
      for(int iproc = 0; iproc < CH4_SPMD::numProc() ; iproc++)
      {
        if(a_printStuff)
        {
          Chombo4::pout() << "DistributedData::unpackReceives iproc = "<< iproc << std::endl;
        }
        auto recvPattern = a_pattern.m_recvPattern[iproc];
        
        if(recvPattern.size() > 0)
        {
          auto & recv_pi   =  a_metacomm.m_recvPI[iproc];
          char* charbuf            =  recv_pi->m_buffer;
          const auto&  vecbuf      =  recv_pi->m_boxbufsize;
          char* bufloc = charbuf;
          
          for(int ivec = 0; ivec < vecbuf.size(); ivec++)
          {
            //this is a box_interaction_t
            const auto& motion  = recvPattern[ivec];
            
            CH_assert(motion.m_dst.m_procID == CH4_SPMD::procID());
            fabtype_t& localfab = *(a_data[motion.m_dst.m_datInd.datInd()]);

            int icomp = 0; int ncomp = localfab.nComp(); //to make it fit with old interface

            size_t expectedCharSize = vecbuf[ivec];
            
            if(a_printStuff)
            {
              Chombo4::pout() << "DistributedData::unpackReceives: ivec= "    <<  ivec
                              << ", region = "   << motion.m_dst.m_region
                              << ", srcProc = "   << motion.m_src.m_procID
                              << ", char_size   = " << expectedCharSize << endl;
            }
            localfab.linearIn(bufloc, expectedCharSize, motion.m_dst.m_region, icomp, ncomp, a_printStuff);
            if(a_printStuff)
            {
              Chombo4::pout() << "DistributedData::unpackReceives: made it out of linearin for ivec = " << ivec << endl;
            }
            
            bufloc += vecbuf[ivec];
          } //end loop over sub buffers
          
        }
      }

      if(a_printStuff)
      {
        Chombo4::pout() << "DistributedData::unpackReceives end" << std::endl;
      }
    }

    //mpi_isend called
    static void  postSends(vector<fabtype_t*>               & a_data,
                                CommunicationMetaData<fabtype_t> & a_metacomm,
                                const BoxPattern                 & a_pattern,
                                bool a_printStuff,
                                bool a_graphConstructed)
    {
      CH_TIME("DistributedData::postSends");
      if(a_printStuff)
      {
        Chombo4::pout() << "DistributedData::postSends begin" << std::endl;
      }
      auto comm = CH4_SPMD::Chombo_MPI::comm;
      for(int iproc = 0; iproc < CH4_SPMD::numProc() ; iproc++)
      {
        auto pattern = a_pattern.m_sendPattern[iproc];
        if(pattern.size() > 0)
        {
          auto & send_pi   =  a_metacomm.m_sendPI[iproc];
          const auto&  dstProc     =  send_pi->m_dstProcID;
          char* charbuf            =  send_pi->m_buffer;
          const auto&  srcProc     =  send_pi->m_srcProcID;
          const auto&  totbufsize  =  send_pi->m_messageLen;
          const auto&  vecbuf      =  send_pi->m_boxbufsize;
          CH_assert(srcProc == CH4_SPMD::procID());
          size_t checkbufsize = 0;
          char* bufloc = charbuf;
          for(int ivec = 0; ivec < vecbuf.size(); ivec++)
          {
            //this is a box_interaction_t
            const auto& motion  = pattern[ivec];

            const fabtype_t& localfab = *(a_data[motion.m_src.m_datInd.datInd()]);
            int icomp = 0;  int ncomp = localfab.nComp(); //
            
            size_t charsize = localfab.charsize(motion.m_src.m_region, icomp, ncomp);
            CH_assert(charsize == vecbuf[ivec]);
            checkbufsize += charsize;
          
            if(a_printStuff)
            {
              Chombo4::pout() << "DistributedData::postSends: send:    ivec= "    <<  ivec
                              << ", region = "   << motion.m_src.m_region 
                              << ", dstProc = "   << motion.m_dst.m_procID
                              << ", charsize = "   << charsize << endl;
            }

            size_t expectedSize = charsize;
            localfab.linearOut(bufloc, expectedSize, motion.m_src.m_region, icomp, ncomp, a_printStuff);
            bufloc += charsize;  
            if(a_printStuff)
            {
              Chombo4::pout() << "DistributedData::postSends: send: made it out of linearOut for ivec = " << ivec << endl;
            }
          }
          CH_assert(checkbufsize == totbufsize);

          //now we have the buffer filled so we can send it to the other proc
          if(a_printStuff)
          {
            Chombo4::pout() << "DistributeData: about to send " << totbufsize <<" to proc " << dstProc<< endl;
          }
          MPI_Isend(charbuf,  totbufsize, MPI_BYTE, dstProc, 0, comm, &(send_pi->m_dataRequest));
        } //if pattern.size() > 0
      } //loop over procs
      
      if(a_printStuff)
      {
        Chombo4::pout() << "DistributedData::postSends end  " << std::endl;
      }
    }
    
    //
    static void  completePending(vector<fabtype_t*>              & a_data,
                                CommunicationMetaData<fabtype_t> & a_metacomm,
                                const BoxPattern                 & a_pattern,
                                bool a_printStuff,
                                bool a_graphConstructed)
    {
      CH_TIME("DistributedData::completePending");
      if(a_printStuff)
      {
        Chombo4::pout() << "DistributedData::completePending begin" << std::endl;
      }
      
      const auto& sendPattern =  a_pattern.m_sendPattern;
      for(int iproc = 0; iproc < CH4_SPMD::numProc(); iproc++)
      {
        auto& send_pi = a_metacomm.m_sendPI[iproc];
        if(sendPattern[iproc].size() > 0)
        {
          if(a_printStuff)
          {
            Chombo4::pout() << "Distributed::completePending:begin wait for send iproc = "<< iproc  << endl;
          }
          int resultSend = MPI_Wait(&(send_pi->m_dataRequest),
                                    &(send_pi->m_dataStatus ));
          
          if(a_printStuff)
          {
            Chombo4::pout() << "Distributed::completePending:end   wait for send iproc = "<< iproc  << endl;
          }
          
          if (resultSend != MPI_SUCCESS)
          {
            Chombo4::pout() << "DistributedData::completePending: " 
                            << "WARNING: send MPI returned " << resultSend << endl;
          }
        }
      }

      const auto& recvPattern =  a_pattern.m_recvPattern;
      for(int iproc = 0; iproc < CH4_SPMD::numProc(); iproc++)
      {
        auto& recv_pi = a_metacomm.m_recvPI[iproc];
        if(recvPattern[iproc].size() > 0)
        {

          if(a_printStuff)
          {
            Chombo4::pout() << "Distributed::completePending:begin wait for recv iproc = "<< iproc  << endl;
          }

          int resultRecv = MPI_Wait(&(recv_pi->m_dataRequest),
                                    &(recv_pi->m_dataStatus ));
          if(a_printStuff)
          {
            Chombo4::pout() << "Distributed::completePending:end   wait for recv iproc = "<< iproc  << endl;
          }
          if (resultRecv != MPI_SUCCESS)
          {
            Chombo4::pout() << "DistributeData::completePending: " 
                            << "WARNING: receive MPI returned " << resultRecv << endl;
          }
        }
      }
    }
#endif    //MPI
    ///
    static void doLocalCopiesWithSerialization(DistributedData       & a_dst,
                                               const DistributedData & a_src,
                                               const BoxPattern   & a_pattern,
                                               bool a_printStuff = false)
    {
      CH_TIME("DistributedData::doLocalCopiesWithSerialization");
      for(int ibuf = 0; ibuf < a_pattern.m_localPattern.size(); ibuf++)
      {
        
        const auto& motion = a_pattern.m_localPattern[ibuf];
        const auto& region = motion.m_dst.m_region;
        CH_assert(motion.m_dst.m_region == motion.m_src.m_region);
        
        const auto & srcfab    = a_src[motion.m_src.m_datInd];
        auto       & dstfab    = a_dst[motion.m_dst.m_datInd];
        
        int icomp = 0;
        int ncomp = srcfab.nComp();
        size_t thisbufsize = srcfab.charsize(motion.m_src.m_region, icomp, ncomp);
        char* buffer = (char*)(malloc(thisbufsize));

        srcfab.linearOut(buffer, thisbufsize, region,  icomp, ncomp);
        dstfab.linearIn( buffer, thisbufsize, region,  icomp, ncomp);

        free(buffer);
      }
    }
    ///
    static void doLocalCopies(DistributedData       & a_dst,
                              const DistributedData & a_src,
                              const BoxPattern   & a_pattern,
                              bool a_printStuff = false)
    {
      CH_TIME("DistributedData::doLocalCopies");
      bool goFaster = !a_printStuff;
      if(goFaster)
      {
        //this is probably faster
        doLocalCopiesWithCopy(a_dst, a_src, a_pattern, a_printStuff);
      }
      else
      {
        //this is primarily about debugging serialization routines 
        doLocalCopiesWithSerialization(a_dst, a_src, a_pattern, a_printStuff);
      }
    }

    ///
    static void doLocalCopiesWithCopy(DistributedData       & a_dst,
                                      const DistributedData & a_src,
                                      const BoxPattern   & a_pattern,
                                      bool a_printStuff = false)
    {
      CH_TIME("DistributedData::doLocalCopiesWithCopy");
      for(int ibuf = 0; ibuf < a_pattern.m_localPattern.size(); ibuf++)
      {
        //component arguments are artifacts of a more civilized age
        //components are all in the type system now
        const auto& motion = a_pattern.m_localPattern[ibuf];
        const auto& srcdatind = motion.m_src.m_datInd;
        //const auto& srcregion = motion.m_src.m_region; // unused
        const auto& dstdatind = motion.m_dst.m_datInd;
        const auto& dstregion = motion.m_dst.m_region;
        //CH_assert(srcregion == dstregion);   //just getting the compiler to shut up
        Proto::Box regx = ProtoCh::getProtoBox(dstregion);
        unsigned int ico = 0;
        unsigned int ncodst = a_dst[dstdatind].nComp();
        unsigned int ncosrc = a_src[srcdatind].nComp();
        CH_assert(ncodst == ncosrc);
        a_dst[dstdatind].copy(a_src[srcdatind], regx, ico, regx, ico, ncodst);
      }
    }
    ///
    

    int m_procID;
    ///actual data 
    vector<fabtype_t*> m_data;
    IntVect m_ghost;
    Chombo4::DisjointBoxLayout m_grids;
    bool m_isDefined;
    bool m_graphConstructed;    
  private:
    //Copy constuction and assignment are disallowed
    //because this simplifies memory management. I like simplicity.
    DistributedData(const DistributedData& a_input);
    void operator=( const DistributedData& a_input);
      
    
  };
  
}
#endif

