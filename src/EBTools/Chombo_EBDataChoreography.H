#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _EB_DATA_CHOREOGRAPHY_blerg__
#define _EB_DATA_CHOREOGRAPHY_blerg__

#include "Chombo_DisjointBoxLayout.H"
#include "Chombo_Pool.H"
#include <vector>
#include "Chombo_ProblemDomain.H"
#include "Chombo_CH_Timer.H"
#include "Chombo_NeighborIterator.H"
#include "Chombo_DataIterator.H"
//for datafactory
#include "Chombo_BoxLayoutData.H"
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
    Chombo4::Box       m_valid;
    int                m_procID;
    Chombo4::DataIndex m_datInd;
    void 
    define(const Chombo4::Box      & a_region,
           const int               & a_procID,
           const Chombo4::DataIndex& a_datind,
           const Chombo4::Box      & a_valid)
    {
      m_datInd = a_datind;
      m_region = a_region;
      m_procID = a_procID;
      m_valid  = a_valid;
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

    //for sorting
    bool operator < (const box_interaction_t& a_input) const
    {
      return m_src.m_region < a_input.m_src.m_region;
    }
    
    box_interaction_t(const Chombo4::DataIndex& a_srcInd,
                      const Chombo4::DataIndex& a_dstInd,
                      const Chombo4::Box      & a_srcRegion,
                      const Chombo4::Box      & a_dstRegion,
                      const int               & a_srcProcID,
                      const int               & a_dstProcID,
                      const Chombo4::Box      & a_srcValid,
                      const Chombo4::Box      & a_dstValid)
    {
      m_src.define(a_srcRegion, a_srcProcID, a_srcInd, a_srcValid);
      m_dst.define(a_dstRegion, a_dstProcID, a_dstInd, a_dstValid);
    }
             
  };

  ///
  /**
     BoxPattern is a class which deals with all the Box intersections associated with copyTo or exchange.
     Periodic boundary conditions are not supported.
  */
  class BoxPattern
  {
  public:
    ///
    /**
       Define for exchange--fills ghost cells around boxes if neighboring boxes are adjacent
    */
    BoxPattern(const Chombo4::DisjointBoxLayout & a_grids,
               const IntVect                    & a_ghost,
               bool a_printStuff = false)
    {
      CH_TIME("BoxPattern::exchange constructor");

      m_recvPattern.resize(CH4_SPMD::numProc());
      m_sendPattern.resize(CH4_SPMD::numProc());
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
            Chombo4::Box srcValid = neighbor;
            Chombo4::Box dstValid = grid;
            box_interaction_t item(Chombo4::DataIndex(nit()), dit(), gridGhostInter, gridGhostInter, neiProcID, myProcID, srcValid, dstValid);
            if (myProcID == neiProcID)
            { // local move
              m_localPattern.push_back(item);
            }
            else
            {
              m_recvPattern[neiProcID].push_back(item);
            }
          }
          Chombo4::Box neighborGhostInter(grid & neighborGhost);
          if (!neighborGhostInter.isEmpty())
          {
            //my data copies to his ghost
            //local case will be taken care of when the data iterator gets to the neighbor box
            if(myProcID != neiProcID)
            {
              Chombo4::Box dstValid = neighbor;
              Chombo4::Box srcValid = grid;
              box_interaction_t item(dit(), Chombo4::DataIndex(nit()), neighborGhostInter, neighborGhostInter, myProcID, neiProcID, srcValid, dstValid);
              m_sendPattern[neiProcID].push_back(item);
            }
          }
        }
      }
      sortAllPatterns();

      
      if(a_printStuff)
      {
        printStuff(string("exchange")); 
      }
    }

    void printStuff(string a_prefix) const
    {
      for(int iproc = 0; iproc < CH4_SPMD::numProc(); iproc++)
      {
        Chombo4::pout() << a_prefix << " BoxPattern constructor:"
                        << " proc = " << iproc 
                        << " m_recvPattern[iproc].size() = " << m_recvPattern[iproc].size()
                        << " m_sendPattern[iproc].size() = " << m_sendPattern[iproc].size() << endl;
                
      }
    }
      
    ///Define for copyTo --- does NOT include ghost cells.
    BoxPattern(const Chombo4::DisjointBoxLayout & a_src,
               const Chombo4::DisjointBoxLayout & a_dst,
               bool a_printStuff = false)
    {
      CH_TIME("BoxPattern::exchange constructor");

      m_recvPattern.resize(CH4_SPMD::numProc());
      m_sendPattern.resize(CH4_SPMD::numProc());
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
            Chombo4::Box dstValid = gridDst;
            Chombo4::Box srcValid = gridSrc;
            box_interaction_t item(Chombo4::DataIndex(litSrc()), Chombo4::DataIndex(litDst()), intersect, intersect, srcProcID, dstProcID, srcValid, dstValid);
            //if neither processor is == myprocID, nothing to do
            if((srcProcID == myProcID) || (dstProcID == myProcID))
            {
              if((srcProcID == myProcID) && (dstProcID == myProcID))
              {
                m_localPattern.push_back(item);
              }
              else if((srcProcID == myProcID) && (dstProcID != myProcID))
              {
                m_sendPattern[dstProcID].push_back(item);
              }
              else if((srcProcID != myProcID) && (dstProcID == myProcID))
              {
                m_recvPattern[srcProcID].push_back(item);
              }
              else
              {
                Chombo4::MayDay::Error("apparently I missed a case");
              }
            } //if(something is on proc)
          } //if there is an intersection
        } //loop over destination boxes
      } //loop over source boxes

      sortAllPatterns();
      if(a_printStuff)
      {
        printStuff(string("copyTo")); 
      }
    }

    void sortAllPatterns()
    {
      //sort everything in sight so that patterns on different procs show up in the same order
      std::sort(m_localPattern.begin(), m_localPattern.end());
      for(int ipat = 0; ipat < m_recvPattern.size(); ipat++)
      {
        std::sort(m_recvPattern[ipat].begin(), m_recvPattern[ipat].end() );
      }
      for(int ipat = 0; ipat < m_sendPattern.size(); ipat++)
      {
        std::sort(m_sendPattern[ipat].begin(), m_sendPattern[ipat].end() );
      }
    }

    //pattern where source == my proc, dest == my proc
    vector<box_interaction_t> m_localPattern;

    //pattern where source != my proc, dest == my proc
    /// outer loop is per proc 
    vector< vector<box_interaction_t> >m_recvPattern;

    //pattern where source == my proc, dest != my proc
    /// outer loop is per proc 
    vector< vector<box_interaction_t> >m_sendPattern;

    
  private:
    ///I like strong construction.
    BoxPattern();
    BoxPattern( const BoxPattern& a_input);
    void operator=(const BoxPattern& a_input);
  };

  /// 
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
          Chombo4::pout() << "CommunicationMetaData: about to recv " << recv_pi->m_metaDataLen <<" from proc " << recv_pi->m_srcProcID << endl;
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
          Chombo4::pout() << "DistributedData: about to receive " << totbufsize << " from proc " << srcProc << endl;
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
        Chombo4::pout() << "DistributedData::unpackReceives iproc = "<< iproc << std::endl;
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
          Chombo4::pout() << "DistributeData: about to send " << totbufsize <<" to proc " << dstProc<< endl;
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
          Chombo4::pout() << "Distributed::completePending:begin wait for send iproc = "<< iproc  << endl;
          int resultSend = MPI_Wait(&(send_pi->m_dataRequest),
                                    &(send_pi->m_dataStatus ));
          Chombo4::pout() << "Distributed::completePending:end   wait for send iproc = "<< iproc  << endl;
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

          Chombo4::pout() << "Distributed::completePending:begin wait for recv iproc = "<< iproc  << endl;
          int resultRecv = MPI_Wait(&(recv_pi->m_dataRequest),
                                    &(recv_pi->m_dataStatus ));
          Chombo4::pout() << "Distributed::completePending:end   wait for recv iproc = "<< iproc  << endl;
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

