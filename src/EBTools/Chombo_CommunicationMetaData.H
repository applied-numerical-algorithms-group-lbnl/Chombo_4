#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _EB_Communication_meta_data_H_
#define _EB_Communication_meta_data_H_

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

#include <unordered_map>
#include <cstdint>

#include "Chombo_SPMD.H"
#ifdef CH_MPI  
#include "mpi.h"
#endif
using namespace Chombo4;
namespace CH4_Data_Choreography
{
  /// Information needed for one interaction between two processors.
  class proc_interaction_t
  {
  public:
    proc_interaction_t()
    {
      m_srcProcID    = -1;
      m_dstProcID    = -1;
      m_messageLen   = 0;
      m_metaDataLen  = 0;
      //this bit is important
      m_buffer       = NULL;
    }

    int                m_srcProcID;
    int                m_dstProcID;
    size_t             m_messageLen;
    size_t             m_metaDataLen;
    char*              m_buffer;
    vector<size_t>     m_boxbufsize;
#ifdef CH_MPI    
    MPI_Request m_dataRequest;
    MPI_Request m_metaRequest;

    MPI_Status  m_dataStatus ;
    MPI_Status  m_metaStatus ;
#endif

             
  };

  ///class to manage meta data for communication between processors (reduces to a bunch of proc_interactions)
  template<class fabtype_t>
  class CommunicationMetaData
  {
  public:
    
    CommunicationMetaData(const BoxPattern   & a_pattern,
                          vector<fabtype_t*> & a_data,
                          bool a_printStuff = false)
      
    {
      CH_TIME("CommunicationMetaData constructor");
#ifdef CH_MPI  
      buildMetaData(  a_pattern, a_data, a_printStuff);
      postReceives(   a_pattern, a_data, a_printStuff);
      postSends   (   a_pattern, a_data, a_printStuff);
      completePending(a_pattern, a_data, a_printStuff);
      unpackReceives( a_pattern, a_data, a_printStuff);
#endif
      }
      
#ifdef CH_MPI
    ///build the data you can from the send side
    void
    buildMetaData(const BoxPattern   & a_pattern,
                  vector<fabtype_t*> & a_data,
                  bool a_printStuff = false)
    {
      m_recvPI.resize(CH4_SPMD::numProc());
      m_sendPI.resize(CH4_SPMD::numProc());
      m_numSends = 0;
      m_numRecvs = 0;
      for(int iproc = 0; iproc < CH4_SPMD::numProc(); iproc++)
      {
        auto& recv_pi = m_recvPI[iproc];
        recv_pi = shared_ptr<proc_interaction_t>(new proc_interaction_t());
        recv_pi->m_boxbufsize = vector<size_t>(a_pattern.m_recvPattern[iproc].size());
        if(a_pattern.m_recvPattern[iproc].size() > 0)
        {
          m_numRecvs++;
          recv_pi->m_srcProcID   = iproc;
          recv_pi->m_dstProcID   = CH4_SPMD::procID();
          recv_pi->m_metaDataLen = (a_pattern.m_recvPattern[iproc].size())*sizeof(size_t);
        } // if there is anything coming from this proc

        auto& send_pi = m_sendPI[iproc];
        send_pi = shared_ptr<proc_interaction_t>(new proc_interaction_t());
        send_pi->m_boxbufsize = vector<size_t>(a_pattern.m_sendPattern[iproc].size());
        if(a_pattern.m_sendPattern[iproc].size() > 0)
        {
          m_numSends++;

          size_t totbufsize = 0;
          for(int ibuf = 0; ibuf < a_pattern.m_sendPattern[iproc].size(); ibuf++)
          {
            const auto& motion  = a_pattern.m_sendPattern[iproc][ibuf];
            CH_assert(motion.m_src.m_procID == CH4_SPMD::procID());
            const fabtype_t& localfab = *(a_data[motion.m_src.m_datInd.datInd()]);
            size_t thisbufsize = localfab.charsize(motion.m_src.m_region, 0, 1);
            send_pi->m_boxbufsize[ibuf] = thisbufsize;
            totbufsize += thisbufsize;
          }
          
          send_pi->m_srcProcID  = CH4_SPMD::procID();
          send_pi->m_dstProcID  = iproc;
          send_pi->m_messageLen = totbufsize;
          send_pi->m_buffer     = (char*)(malloc(totbufsize));
          send_pi->m_metaDataLen= (a_pattern.m_sendPattern[iproc].size())*sizeof(size_t);

        } //if there is anything to send to this proc
      }   //loop over procs
   }

    ///let MPI know what messages this processor is expecting
    void
    postReceives(const BoxPattern   & a_pattern,
                 vector<fabtype_t*> & a_data,
                 bool a_printStuff = false)
    {
      CH_TIME("CommunicationMetaData::postReceives");

      for(int iproc = 0; iproc < CH4_SPMD::numProc(); iproc++)
      {
        auto recvPattern = a_pattern.m_recvPattern[iproc];
        if(recvPattern.size() > 0)
        {
          auto& recv_pi = m_recvPI[iproc];
          recv_pi->m_srcProcID   = iproc;
          recv_pi->m_dstProcID   = CH4_SPMD::procID();
          recv_pi->m_boxbufsize.resize(recvPattern.size());
          recv_pi->m_metaDataLen =    (recvPattern.size())*sizeof(size_t);

          char* recvbuf = (char*)(recv_pi->m_boxbufsize.data());
          auto comm = CH4_SPMD::Chombo_MPI::comm;
          if(a_printStuff)
          {
            Chombo4::pout() << "CommunicationMetaData: about to recv " << recv_pi->m_metaDataLen <<" from proc " << recv_pi->m_srcProcID << endl;
          }
          MPI_Irecv(recvbuf, recv_pi->m_metaDataLen, MPI_BYTE, recv_pi->m_srcProcID, MPI_ANY_TAG, comm, &recv_pi->m_metaRequest);
          
        } //loop over processors
      }
    }

    ///get the numbers out of the string that MPI sent us
    void
    unpackReceives(const BoxPattern   & a_pattern,
                   vector<fabtype_t*> & a_data,
                   bool a_printStuff = false)
    {
      CH_TIME("CommunicationMetaData::unpackReceives");

      for(int iproc = 0; iproc < CH4_SPMD::numProc(); iproc++)
      {
        auto recvPattern = a_pattern.m_recvPattern[iproc];
        if(recvPattern.size() > 0)
        {
          auto& recv_pi = m_recvPI[iproc];
          size_t totbufsize = 0;
          for(int ibuf = 0; ibuf < recvPattern.size(); ibuf++)
          {
            totbufsize += recv_pi->m_boxbufsize[ibuf];
            if(a_printStuff)
            {
              Chombo4::pout() << "CommunicationMetaData::unpackReceives: "
                              << "ToMe  : ibuf = "       << ibuf
                              << "  bufSize = "   << recv_pi->m_boxbufsize[ibuf]
                              << "  srcProcID = " << recvPattern[ibuf].m_src.m_procID
                              << "  dstProcID = " << recvPattern[ibuf].m_dst.m_procID 
                              << ", srcRegion = " << recvPattern[ibuf].m_src.m_region
                              << ", dstRegion = " << recvPattern[ibuf].m_dst.m_region 
                              << ", srcValid  = " << recvPattern[ibuf].m_src.m_valid
                              << ", dstValid  = " << recvPattern[ibuf].m_dst.m_valid << std::endl;
            }
          }
          if(a_printStuff)
          {
            Chombo4::pout() << "CommunicationMetaData::unpackReceives: "
                            << "proc = "     << recv_pi->m_dstProcID
                            << "\t will receive " << recvPattern.size()
                            << "\t boxes worth of data for a total of a " << totbufsize
                            << "-sized buffer from proc " << recv_pi->m_srcProcID << endl;
          }
          recv_pi->m_messageLen = totbufsize;
          recv_pi->m_buffer     = (char*)(malloc(totbufsize));
        } //
      }
    }

    /// send messages through MPI
    void
    postSends(const BoxPattern   & a_pattern,
              vector<fabtype_t*> & a_data,
              bool a_printStuff = false)
    {
      CH_TIME("CommunicationMetaData::postSends");
      const auto& sendPattern =  a_pattern.m_sendPattern;
      for(int iproc = 0; iproc < CH4_SPMD::numProc(); iproc++)
      {

        if(sendPattern[iproc].size() > 0)
        {
          auto& send_pi = m_sendPI[iproc];
          size_t totbufsize = 0;
          send_pi->m_boxbufsize.resize(sendPattern[iproc].size());
          for(int ibuf = 0; ibuf < sendPattern[iproc].size(); ibuf++)
          {
            //not constant because I am using the request
            const auto& motion  = sendPattern[iproc][ibuf];
            CH_assert(motion.m_src.m_procID == CH4_SPMD::procID());
            const fabtype_t& localfab = *(a_data[motion.m_src.m_datInd.datInd()]);
            size_t thisbufsize = localfab.charsize(motion.m_src.m_region, 0, 1);

            send_pi->m_boxbufsize[ibuf] = thisbufsize;
            totbufsize += thisbufsize;
            if(a_printStuff)
            {
              Chombo4::pout() << "CommunicationMetaData::postSends:FromMe: ibuf = "       << ibuf
                              << "  bufSize   = "   << thisbufsize
                              << "  srcProcID = " << sendPattern[iproc][ibuf].m_src.m_procID
                              << "  dstProcID = " << sendPattern[iproc][ibuf].m_dst.m_procID 
                              << ", srcRegion = " << sendPattern[iproc][ibuf].m_src.m_region
                              << ", dstRegion = " << sendPattern[iproc][ibuf].m_dst.m_region 
                              << ", srcValid  = " << sendPattern[iproc][ibuf].m_src.m_valid
                              << ", dstValid  = " << sendPattern[iproc][ibuf].m_dst.m_valid << std::endl;
            }
          }

          send_pi->m_srcProcID  = CH4_SPMD::procID();
          send_pi->m_dstProcID  = iproc;
          send_pi->m_messageLen = totbufsize;
          send_pi->m_buffer     = (char*)(malloc(totbufsize));
          send_pi->m_metaDataLen= (sendPattern[iproc].size())*sizeof(size_t);

          
          if(a_printStuff)
          {
            Chombo4::pout() << "CommunicationMetaData::postSends:proc = "                   << send_pi->m_srcProcID
                            << "\t will send    " <<  sendPattern[iproc].size()
                            << "\t boxes worth of data for a total of a " <<  totbufsize 
                            << "-sized  buffer to  proc " << send_pi->m_dstProcID << endl;
          }
          char* sendbuf = (char*)(send_pi->m_boxbufsize.data());
          auto comm = CH4_SPMD::Chombo_MPI::comm;
          if(a_printStuff)
          {
            Chombo4::pout() << "CommunicationMetaData:postSends: about to send " << send_pi->m_metaDataLen <<" to   proc " << send_pi->m_srcProcID << endl;
          }
          MPI_Isend(sendbuf, send_pi->m_metaDataLen, MPI_BYTE, send_pi->m_dstProcID, 0, comm, &(send_pi->m_metaRequest));
        } //if there is anything to send to this proc
      }   //loop over procs
    }


    //call MPI_Waitall for both sends and receives
    void
    completePending(const BoxPattern   & a_pattern,
                    vector<fabtype_t*> & a_data,
                    bool a_printStuff = false)
    {
      CH_TIME("CommunicationMetaData::completePending");
      const auto& sendPattern =  a_pattern.m_sendPattern;
 
      if(a_printStuff)
      {
        Chombo4::pout() << "CommunicationMetaData::completePending: begin "<< endl;
      }
      for(int iproc = 0; iproc < CH4_SPMD::numProc(); iproc++)
      {
        auto& send_pi = m_sendPI[iproc];
        if(sendPattern[iproc].size() > 0)
        {
          if(a_printStuff)
          {
            Chombo4::pout() << "CommunicationMetaData::completePending:send iproc = "<< iproc  << endl;
          }
          int resultSend = MPI_Wait(&(send_pi->m_metaRequest),
                                    &(send_pi->m_metaStatus ));
          if (resultSend != MPI_SUCCESS)
          {
            Chombo4::pout() << "CommunicationMetaData::completePending: " 
                            << "WARNING: send MPI returned " << resultSend << endl;
          }
        }
      }

      const auto& recvPattern =  a_pattern.m_recvPattern;
      for(int iproc = 0; iproc < CH4_SPMD::numProc(); iproc++)
      {
        auto& recv_pi = m_recvPI[iproc];
        if(recvPattern[iproc].size() > 0)
        {
          if(a_printStuff)
          {
            Chombo4::pout() << "CommunicationMetaData::completePending:recv iproc = "<< iproc  << endl;
          }
          int resultRecv = MPI_Wait(&(recv_pi->m_metaRequest),
                                    &(recv_pi->m_metaStatus ));
          if (resultRecv != MPI_SUCCESS)
          {
            Chombo4::pout() << "CommunicationMetaData::completePending: " 
                            << "WARNING: receive MPI returned " << resultRecv << endl;
          }
        }
        if(a_printStuff)
        {
          Chombo4::pout() << "CommunicationMetaData::completePending: end"  << endl;
        }
      }
    }
#endif      

    ///Because we have malloced, so we must free
    ~CommunicationMetaData()
    {
      for(int iinter = 0; iinter < m_recvPI.size(); iinter++)
      {
        auto&  interaction = m_recvPI[iinter];
        if(interaction->m_buffer != NULL)
        {
          free(interaction->m_buffer);
          interaction->m_buffer = NULL;
        }
      }
      m_recvPI.resize(0);

      for(int iinter = 0; iinter < m_sendPI.size(); iinter++)
      {
        auto&  interaction = m_sendPI[iinter];
        if(interaction->m_buffer != NULL)
        {
          free(interaction->m_buffer);
          interaction->m_buffer = NULL;
        }
      }
      m_sendPI.resize(0);
    }
    
    int m_numSends, m_numRecvs;
    //everything length numProc
    vector<shared_ptr<proc_interaction_t> > m_recvPI;
    vector<shared_ptr<proc_interaction_t> > m_sendPI;

  };

  
}
#endif

