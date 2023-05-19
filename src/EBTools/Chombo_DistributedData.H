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
  template<class fabtype_t>
  class MetaDataPair
  {
    public:

    MetaDataPair(const Chombo4::DisjointBoxLayout& a_srcGrids,
                 const Chombo4::DisjointBoxLayout& a_dstGrids,
                 std::vector<fabtype_t *>        & a_data,
                 bool a_printStuff = false)
    {
      m_pattern  = std::shared_ptr<BoxPattern>(
        new BoxPattern(a_srcGrids, a_dstGrids, a_printStuff)
        );

      m_commMeta = std::shared_ptr<CommunicationMetaData<fabtype_t> >(
        new CommunicationMetaData<fabtype_t>(*m_pattern, a_data, a_printStuff)
        );
    }

    MetaDataPair(const Chombo4::DisjointBoxLayout& a_grids,
                 const Chombo4::IntVect          & a_ghost,
                 std::vector<fabtype_t *>        & a_data,
                 bool a_printStuff = false)
    {
      m_pattern  = std::shared_ptr<BoxPattern>(
        new BoxPattern(a_grids, a_ghost, a_printStuff)
        );

      m_commMeta = std::shared_ptr<CommunicationMetaData<fabtype_t> >(
        new CommunicationMetaData<fabtype_t>(*m_pattern, a_data, a_printStuff)
        );
    }
    
    std::shared_ptr<BoxPattern                       >  m_pattern;
    std::shared_ptr<CommunicationMetaData<fabtype_t> >  m_commMeta;
  private:

    MetaDataPair();

      
  };

  
  class exchange_key_t
  {
  public:

    exchange_key_t(Chombo4::DisjointBoxLayout a_grids,
                   Chombo4::IntVect           a_ghost):
      m_grids(a_grids), m_ghost(a_ghost){}
      

    bool operator==(const exchange_key_t& a_input) const
    {
      return ((m_grids == a_input.m_grids)&&
              (m_ghost == a_input.m_ghost));
    }
      
    Chombo4::DisjointBoxLayout m_grids;
    Chombo4::IntVect           m_ghost;

  private:
    exchange_key_t();
      
  };
  

  class copy_to__key_t
  {
  public:

    copy_to__key_t(Chombo4::DisjointBoxLayout a_srcGrids,
                  Chombo4::DisjointBoxLayout a_dstGrids):
      m_srcGrids(a_srcGrids), m_dstGrids(a_dstGrids){}
      

    bool operator==(const copy_to__key_t& a_input) const
    {
      return ((m_srcGrids == a_input.m_srcGrids)&&
              (m_dstGrids == a_input.m_dstGrids));
    }
      
    Chombo4::DisjointBoxLayout m_srcGrids;
    Chombo4::DisjointBoxLayout m_dstGrids;

  private:
    copy_to__key_t();
      
  };

  ///Optimization wherein we archive the meta data for given inputs.   
  template<class fabtype_t>
  class MetaDataCache
  {
  public:
    typedef std::pair<copy_to__key_t  ,std::shared_ptr<MetaDataPair<fabtype_t> > >  copy_to__pair_t;
    typedef std::pair<exchange_key_t ,std::shared_ptr<MetaDataPair<fabtype_t> > > exchange_pair_t;
    
    std::shared_ptr< MetaDataPair<fabtype_t> > getCopyToPair(const Chombo4::DisjointBoxLayout& a_srcGrids,
                                                             const Chombo4::DisjointBoxLayout& a_dstGrids,
                                                             std::vector<fabtype_t *> a_data,
                                                             bool a_printStuff = false)
    {
      typedef std::shared_ptr< MetaDataPair<fabtype_t> > return_type_t;

      return_type_t  retval;
      copy_to__key_t key(a_srcGrids, a_dstGrids);
      //see if already there
      bool found = false;
      for(int ivec = 0; ivec < m_copy_to__cache.size(); ivec++)
      {
        if(m_copy_to__cache[ivec].first == key)
        {
          found = true;
          //already there so we can just return it.
          retval = m_copy_to__cache[ivec].second;
        }
      }
      if(!found)
      {        
        // not already there so make one and save it.
        retval = return_type_t(new MetaDataPair<fabtype_t>(a_srcGrids, a_dstGrids, a_data, a_printStuff) );
        
        copy_to__pair_t  new_entry(key, retval);

        m_copy_to__cache.push_back(new_entry);
      }

      return retval;
    }

    std::shared_ptr< MetaDataPair<fabtype_t> > getExchangePair(const Chombo4::DisjointBoxLayout& a_grids,
                                                               const Chombo4::IntVect&           a_ghost,
                                                               std::vector<fabtype_t *>    a_data,
                                                               bool a_printStuff = false)
    {
      typedef std::shared_ptr< MetaDataPair<fabtype_t> > return_type_t;

      return_type_t  retval;
      exchange_key_t key(a_grids, a_ghost);
      //see if already there
      bool found = false;
      for(int ivec = 0; ivec < m_exchange_cache.size(); ivec++)
      {
        if(m_exchange_cache[ivec].first == key)
        {
          found = true;
          retval = m_exchange_cache[ivec].second;
        }
      }
      if(!found)
      {        
        // not already there
        retval = return_type_t(new MetaDataPair<fabtype_t>(a_grids, a_ghost, a_data, a_printStuff) );
        
        exchange_pair_t  new_entry(key, retval);

        m_exchange_cache.push_back(new_entry);
      }
      return retval;
    }
    
    std::vector<copy_to__pair_t>  m_copy_to__cache;
    std::vector<exchange_pair_t>  m_exchange_cache;
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
                    const Chombo4::DataFactory<fabtype_t>& a_factory,
                    bool                                   a_enableStaticCache = true)
    {
      define(a_grids, a_ghost, a_factory, a_enableStaticCache);
    }

    
    ///datafactory has an ncomps argument that shall be ignored here
    void define(const Chombo4::DisjointBoxLayout     & a_grids,
                const Chombo4::IntVect               & a_ghost,
                const Chombo4::DataFactory<fabtype_t>& a_factory,
                bool                                   a_enableStaticCache = true)
    {
      CH_TIME("DistributedData constructor");
      m_procID = CH4_SPMD::procID();
      m_isDefined = true;
      m_grids = a_grids;
      m_ghost = a_ghost;
      m_graphConstructed = a_factory.graphConstructor();
      m_enableStaticCache = a_enableStaticCache;
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
      if (m_enableStaticCache)
        {
          auto& cache = getMetaDataCache();
          auto metapair_ptr = cache.getExchangePair(m_grids, m_ghost, m_data, a_printStuff);
          if (a_printStuff) {Chombo4::pout() << "DistributedData::exchange: using cache" << endl;}
          BoxPattern& pattern                          = *(metapair_ptr->m_pattern );
          CommunicationMetaData<fabtype_t>& metadataco = *(metapair_ptr->m_commMeta);
          communicate(*this, pattern, metadataco, *this  , a_printStuff, m_graphConstructed);
        }
      else
        {
          if (a_printStuff) {Chombo4::pout() << "DistributedData::exchange: eschewing cache!" << endl;}
          BoxPattern pattern(m_grids, m_ghost, a_printStuff);
          CommunicationMetaData<fabtype_t> metadataco(pattern, m_data, a_printStuff);
          communicate(*this, pattern, metadataco, *this  , a_printStuff, m_graphConstructed);
        }
    }
    
    ///
    void copyTo(DistributedData & a_dst,
                bool a_printStuff = false)
    {
      CH_TIME("DistributedData::copyTo");
      if (m_enableStaticCache)
        {
          auto& cache = getMetaDataCache();
          auto metapair_ptr = cache.getCopyToPair(m_grids, a_dst.m_grids, m_data, a_printStuff);
          if (a_printStuff) {Chombo4::pout() << "DistributedData::copyTo: using cache" << endl;}
          BoxPattern& pattern                          = *(metapair_ptr->m_pattern  );
          CommunicationMetaData<fabtype_t>& metadataco = *(metapair_ptr->m_commMeta);
          communicate(a_dst,  pattern, metadataco, *this, a_printStuff, m_graphConstructed);
        }
      else
        {
          if (a_printStuff) {Chombo4::pout() << "DistributedData::copyTo: eschewing cache!" << endl;}
          BoxPattern pattern(m_grids, a_dst.m_grids, true);
          CommunicationMetaData<fabtype_t> metadataco(pattern,
                                                      m_data, true);
          if (a_printStuff) {Chombo4::pout() << "DistributedData::copyTo: m_grids: " << endl;}
          m_grids.print();
          if (a_printStuff) {Chombo4::pout() << "DistributedData::copyTo: a_dst.m_grids:" << endl;}
          a_dst.m_grids.print();
          communicate(a_dst,  pattern, metadataco, *this, a_printStuff, m_graphConstructed);
        }
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
    
    void printGrids() const
    {
      Chombo4::pout() << "grids = " ;
      m_grids.print();
      Chombo4::pout() << std::endl;
    }
    
  private:

    static void communicate(DistributedData                  & a_dst,
                            BoxPattern                       & a_pattern,
                            CommunicationMetaData<fabtype_t> & a_metadataco,
                            DistributedData                  & a_src,
                            bool a_printStuff,
                            bool a_graphConstructed)
    {
      CH_TIME("DistributedData::communicate");
#ifdef CH_MPI  
      postReceives(   a_dst.m_data, a_metadataco, a_pattern, a_printStuff, a_graphConstructed);
      postSends(      a_src.m_data, a_metadataco, a_pattern, a_printStuff, a_graphConstructed);
#endif      
      doLocalCopies(a_dst, a_src, a_pattern, a_printStuff);
#ifdef CH_MPI  
      completePending(a_dst.m_data, a_metadataco, a_pattern, a_printStuff, a_graphConstructed);
      unpackReceives (a_dst.m_data, a_metadataco, a_pattern, a_printStuff, a_graphConstructed);
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
            if(a_printStuff)
            {
              Chombo4::pout() << "DistributedData::postSends: send:    ivec= "    <<  ivec
                              << ", region = "   << motion.m_src.m_region
                              << ", dstProc = "   << motion.m_dst.m_procID
                              << ", charsize = "   << charsize
                              << ", vecbuf = "   << vecbuf[ivec] << endl;
            }
            CH_assert(charsize == vecbuf[ivec]);
            checkbufsize += charsize;
          
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
      //goFaster = true;
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
        auto      & dstfab = a_dst[dstdatind];
        const auto& srcfab = a_src[srcdatind];
        dstfab.copy(srcfab, regx, ico, regx, ico, ncodst);
      }
    }
    ///

    static MetaDataCache<fabtype_t>& getMetaDataCache()
    {
      static MetaDataCache<fabtype_t>* s_meta_data_ptr = NULL;
      static bool allocated = false;
      if(!allocated)
      {
        s_meta_data_ptr = new MetaDataCache<fabtype_t>();
        allocated = true;
      }
      return *s_meta_data_ptr;
    }
    
    int m_procID;
    ///actual data 
    vector<fabtype_t*> m_data;
    IntVect m_ghost;
    Chombo4::DisjointBoxLayout m_grids;
    bool m_isDefined;
    bool m_enableStaticCache;
    bool m_graphConstructed;

  private:
    //Copy constuction and assignment are disallowed
    //because this simplifies memory management. I like simplicity.
    DistributedData(const DistributedData& a_input);
    void operator=( const DistributedData& a_input);
      
    
  };
  
}
#endif

