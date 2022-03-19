#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _Chombo_BOXLAYOUTDATAI_H_
#define _Chombo_BOXLAYOUTDATAI_H_

#include <cstdlib>
#include <algorithm>
#include <limits.h>
#include <list>
#include "Chombo_CH_OpenMP.H"
#include "Chombo_parstream.H"
#include "Chombo_CH_Timer.H"
#include "Chombo_NamespaceHeader.H"

using std::sort;


template <class T>
T* DefaultDataFactory<T>::create(const Box& box,
                                 int ncomps,
                                 const DataIndex& a_datInd) const
{
  Proto::Box bx = ProtoCh::getProtoBox(box);
  return new T(bx, ncomps);
}

template<class T>
inline bool BoxLayoutData<T>::isDefined() const
{
  return m_isdefined;
}

template <class T>
inline void BoxLayoutData<T>::setVector(const BoxLayoutData<T>& da,
                                        const Interval& srcComps,
                                        const Interval& destComps)
{
  if(&da != this)
  {
    DataIterator it=this->dataIterator();
    int nbox=it.size();
    for(int box=0; box<nbox; box++)
    {
//      this->m_vector[it[box].datInd()]->copy( this->box(it[box]), destComps,
//                                              this->box(it[box]), da[it[box]], srcComps);
      int isrc =  srcComps.begin();
      int idst = destComps.begin();
      int inco = destComps.size();
      ::Proto::Box bx = ProtoCh::getProtoBox(this->box(it[box]));
      this->m_vector[it[box].datInd()]->copy(da[it[box]], bx, isrc, bx, idst, inco);

    }
  }
}

template<class T>
inline void BoxLayoutData<T>::define(const BoxLayoutData<T>& da, const Interval& comps,
                                     const DataFactory<T>& factory)
{
  if (this == &da)
  {
    MayDay::Error("BoxLayoutData<T>::define(const LayoutData<T>& da,.....) called with 'this'");
  }
  CH_assert(comps.size()>0);
  CH_assert(comps.end()<=m_comps);
  //AD: why are the two different 
  //  CH_assert(comps.end()<=da.m_comps);
  CH_assert(comps.begin()>=0);
  this->m_boxLayout = da.boxLayout();

  this->m_comps = comps.size();

  
  Interval dest(0, m_comps-1);
  allocateGhostVector(factory);
  setVector(da, comps, dest);
}

template<class T>
inline void BoxLayoutData<T>::define(const BoxLayout& boxes, int comps,
                                     const DataFactory<T>& factory)
{
  CH_assert(boxes.isClosed());
  this->m_boxLayout = boxes;
  m_comps = comps;

  m_isdefined = true;
  allocateGhostVector(factory);

}

template<class T>
inline void BoxLayoutData<T>::define(const BoxLayout& boxes)
{
  MayDay::Error("BoxLayoutData<T>::define(const BoxLayout& boxes)...needs comps");
}

template <class T>
inline   BoxLayoutData<T>::BoxLayoutData():m_comps(0) ,m_buff(0)
{
  m_isdefined = false;
#ifdef CH_MPI

#endif
}
template<class T>
inline BoxLayoutData<T>::BoxLayoutData(const BoxLayout& boxes, int comps,
                                       const DataFactory<T>& factory)
  :m_comps(comps),m_buff(NULL)
{
  CH_assert(boxes.isClosed());
  this->m_boxLayout = boxes;
  m_isdefined = true;
  allocateGhostVector(factory);
#ifdef CH_MPI

#endif
}

template<class T>
BoxLayoutData<T>::~BoxLayoutData()
{
  CH_TIME("~BoxLayoutData");
  //completePendingSends();
}

template<class T>
inline void BoxLayoutData<T>::define(const BoxLayoutData<T>& da,
                                     const DataFactory<T>& factory)
{
  if (this != &da)
  {
    m_isdefined = da.m_isdefined;
    this->m_boxLayout = da.boxLayout();
    m_comps    = da.nComp();
    Interval srcAnddest(0, m_comps-1);
    allocateGhostVector(factory);
    setVector(da, srcAnddest, srcAnddest);
  }

}
template<class T>
inline void BoxLayoutData<T>::clear()
{
  if (this->m_callDelete == true)
  {
    for (unsigned int i=0; i<this->m_vector.size(); ++i)
    {
      delete this->m_vector[i];
      this->m_vector[i] = NULL;
    }
  }
  m_isdefined = false;
}

template<class T>
inline void BoxLayoutData<T>::allocateGhostVector(const DataFactory<T>& factory, const IntVect& ghost)
{
  //begin debug
  //Chombo4::pout() << "going into BoxLayoutData::allocateGhostVector" << endl;
  //end debug
  
  if (this->m_callDelete == true)
  {
    for (unsigned int i=0; i<this->m_vector.size(); ++i)
    {
      delete this->m_vector[i];
      this->m_vector[i] = NULL;
    }
  }

  this->m_callDelete = factory.callDelete();

  DataIterator it(this->dataIterator()); int nbox=it.size();
  this->m_vector.resize(it.size(), NULL);

  for(int i=0; i<nbox; i++)
  {
    unsigned int index = it[i].datInd();
    Box abox = this->box(it[i]);
    abox.grow(ghost);
    this->m_vector[index] = factory.create(abox, m_comps, it[i]);
    if (this->m_vector[index] == NULL)
    {
      MayDay::Error("OutOfMemory in BoxLayoutData::allocateGhostVector");
    }
  }
  //begin debug
  //Chombo4::pout() << "leaving BoxLayoutData::allocateGhostVector" << endl;
  //end debug
}

template<class T>
inline void BoxLayoutData<T>::apply(void (*a_func)(const Box& box, int comps, T& t))
{
  DataIterator it(this->dataIterator()); int nbox=it.size();

  for(int i=0; i<nbox; i++)
  {
    a_func(this->box(it[i]), m_comps, *(this->m_vector[ it[i].datInd() ]));
  }
}

//======================================================================
template <class T>
AliasDataFactory<T>::AliasDataFactory(BoxLayoutData<T>* a_original, const Interval& interval)
{
  define(a_original, interval);
}

template <class T>
void AliasDataFactory<T>::define(BoxLayoutData<T>* a_original, const Interval& interval)
{
  m_origPointer = a_original;
  m_interval    = interval;
}

template <class T>
T* AliasDataFactory<T>::create(const Box& a_box, int ncomps, const DataIndex& a_dataInd) const
{
  //CH_assert(this->box(a_dataInd) == a_box);
  CH_assert(ncomps = m_interval.size());
  T* rtn = new T(m_interval, m_origPointer->operator[](a_dataInd));
  return rtn;
}

template<class T>
void
BoxLayoutData<T>::
makeItSo(const Interval&   a_srcComps,
         const BoxLayoutData<T>& a_src,
         BoxLayoutData<T>& a_dest,
         const Interval&   a_destComps,
         const Copier&     a_copier,
         const LDOperator<T>& a_op,
         bool a_printStuff) const
{
  if(a_printStuff)
  {
    pout() << "boxlayoutdata::makeitso copier = " << endl;
    a_copier.print();
  }
  if(a_printStuff)
  {
    pout() << "boxlayoutdata::makeitso going into makeitsobegin " << endl;
  }
  {
    CH_TIME("makeItSoBegin");
    makeItSoBegin(a_srcComps, a_src, a_dest, a_destComps, a_copier, a_op, a_printStuff);
  }
  if(a_printStuff)
  {
    pout() << "boxlayoutdata::makeitso going into makeitsolocalcopy " << endl;
  }
  {
    CH_TIME("makeItSoLocalCopy");
    makeItSoLocalCopy(a_srcComps, a_src, a_dest, a_destComps, a_copier, a_op, a_printStuff);
  }
  if(a_printStuff)
  {
    pout() << "boxlayoutdata::makeitso going into makeitsoend " << endl;
  }
  {
    CH_TIME("makeItSoEnd");
    a_dest.makeItSoEnd(a_destComps, a_op, a_printStuff);
  }
}
///
template<class T>
void
BoxLayoutData<T>::
makeItSoBegin(const Interval&   a_srcComps,
              const BoxLayoutData<T>& a_src,
              BoxLayoutData<T>& a_dest,
              const Interval&   a_destComps,
              const Copier&     a_copier,
              const LDOperator<T>& a_op,
              bool a_printStuff) const
{
  // The following five functions are nullOps in uniprocessor mode

#ifdef CH_MPI


  allocateBuffers(a_src,  a_srcComps,
                  a_dest, a_destComps,
                  a_copier,
                  a_op, a_printStuff);  //monkey with buffers, set up 'fromMe' and 'toMe' queues


  writeSendDataFromMeIntoBuffers(a_src, a_srcComps, a_op);

  // If there is nothing to recv/send, don't go into these functions
  // and allocate memory that will not be freed later.  (ndk)
  // The #ifdef CH_MPI is for the m_buff->m_toMe and m_buff->m_fromMe
  {
    CH_TIME("post_messages");
    this->m_buff->numReceives = m_buff->m_toMe.size();

    if (this->m_buff->numReceives > 0)
    {
      postReceivesToMe(); // all non-blocking
    }
  

    this->m_buff->numSends = m_buff->m_fromMe.size();
    if (this->m_buff->numSends > 0)
    {
      postSendsFromMe();  // all non-blocking
    }
  }    
#endif 
}
///
template<class T>
void
BoxLayoutData<T>::
makeItSoLocalCopy(const Interval&   a_srcComps,
                  const BoxLayoutData<T>& a_src,
                  BoxLayoutData<T>& a_dest,
                  const Interval&   a_destComps,
                  const Copier&     a_copier,
                  const LDOperator<T>& a_op,
                  bool a_printStuff) const
{

  CH_TIME("local copying");
  CopyIterator it(a_copier, CopyIterator::LOCAL);  
  int items=it.size();

  for (int n=0; n<items; n++)
  {
    const MotionItem& item = it[n];

    a_op.op(a_dest[item.toIndex], item.fromRegion,
            a_destComps,
            item.toRegion,
            a_src[item.fromIndex],
            a_srcComps);
  }

}

///


template<class T>
void
BoxLayoutData<T>::
makeItSoEnd(const Interval&   a_destComps,
            const LDOperator<T>& a_op,
            bool a_printStuff) 

{
  // Uncomment and Move this out of unpackReceivesToMe()  (ndk)
  completePendingSends(); // wait for sends from possible previous operation

  unpackReceivesToMe(a_destComps, a_op); // nullOp in uniprocessor mode
}



// MPI versions of the above codes.

template<class T>
void BoxLayoutData<T>::completePendingSends() const
{
  CH_TIME("completePendingSends");
#ifdef CH_MPI  
  if (this->m_buff->numSends > 0)
  {
    CH_TIME("MPI_Waitall");
    m_buff->m_sendStatus.resize(this->m_buff->numSends);
    int result = MPI_Waitall(this->m_buff->numSends, &(m_buff->m_sendRequests[0]), &(m_buff->m_sendStatus[0]));
    if (result != MPI_SUCCESS)
    {
      //hell if I know what to do about failed messaging here
    }
  }
  this->m_buff->numSends = 0;
#endif  
}

template<class T>
void BoxLayoutData<T>::allocateBuffers(const BoxLayoutData<T>& a_src,
                                       const Interval& a_srcComps,
                                       const BoxLayoutData<T>& a_dest,
                                       const Interval& a_destComps,
                                       const Copier&   a_copier,
                                       const LDOperator<T>& a_op,
                                       bool a_printStuff) const
{
  
  CH_TIME("MPI_allocateBuffers");
#ifdef CH_MPI  
  m_buff = &(((Copier&)a_copier).m_buffers);
  a_dest.m_buff = m_buff;

  CH_assert(a_srcComps.size() == a_destComps.size());
  if (m_buff->isDefined(a_srcComps.size()) && T::preAllocatable()<2) return;

  if(a_printStuff)
  {
    pout() << " allocate buffers srcComps = " << a_srcComps << ", dest comps = " << a_destComps << endl;
  }
  m_buff->m_ncomps = a_srcComps.size();

  m_buff->m_fromMe.resize(0);
  m_buff->m_toMe.resize(0);
  size_t sendBufferSize = 0;
  size_t recBufferSize  = 0;
  // two versions of code here.  one for preAllocatable T, one not.

  T dummy;
  for (CopyIterator it(a_copier, CopyIterator::FROM); it.ok(); ++it)
  {
    const MotionItem& item = it();
    CopierBuffer::bufEntry b;
    b.item = &item;
    b.size = a_op.size(a_src[item.fromIndex], item.fromRegion, a_srcComps);
    sendBufferSize+=b.size;
    b.procID = item.procID;
    m_buff->m_fromMe.push_back(b);
  }
  sort(m_buff->m_fromMe.begin(), m_buff->m_fromMe.end());
  for (CopyIterator it(a_copier, CopyIterator::TO); it.ok(); ++it)
  {
    const MotionItem& item = it();
    CopierBuffer::bufEntry b;
    b.item = &item;
    if (T::preAllocatable() == 0)
    {
      b.size = a_op.size(dummy, item.fromRegion, a_destComps);
      recBufferSize+=b.size;
    }
    else if (T::preAllocatable() == 1)
    {
      b.size = a_op.size(a_dest[item.toIndex], item.fromRegion, a_destComps);
      recBufferSize+=b.size;
    }
    b.procID = item.procID;
    m_buff->m_toMe.push_back(b);
  }
  sort(m_buff->m_toMe.begin(), m_buff->m_toMe.end());

  if (T::preAllocatable() == 2) // dynamic allocatable, need two pass
  {
    CH_TIME("MPI_ Phase 1 of 2 Phase: preAllocatable==2");
    if (a_printStuff)
    {
      pout()<<"preAllocatable==2\n";
    }

    // in the non-preallocatable case, I need to message the
    // values for the m_buff->m_toMe[*].size
    Vector<unsigned long> fdata;
    Vector<unsigned long> tdata;
    int count = 1;
    int scount = 1;
    if (m_buff->m_toMe.size() > 0)
    {
      tdata.resize(m_buff->m_toMe.size(), ULONG_MAX);
      m_buff->m_receiveRequests.resize(numProc()-1);
      m_buff->m_receiveStatus.resize(numProc()-1);
      MPI_Request* Rptr = &(m_buff->m_receiveRequests[0]);

      unsigned int lastProc = m_buff->m_toMe[0].procID;
      int messageSize = 1;
      unsigned long * dataPtr = &(tdata[0]);
      unsigned int i = 1;

      for (;i<m_buff->m_toMe.size(); ++i)
      {
        CopierBuffer::bufEntry& b = m_buff->m_toMe[i];
        if (b.procID == lastProc)
          messageSize++;
        else
        {

          MPI_Irecv(dataPtr, messageSize, MPI_UNSIGNED_LONG, lastProc,
                    1, CH4_SPMD::Chombo_MPI::comm, Rptr);
          Rptr++;

          lastProc = b.procID;
          messageSize = 1;
          dataPtr = &(tdata[i]);
          count++;
        }
      }

      MPI_Irecv(dataPtr, messageSize, MPI_UNSIGNED_LONG, lastProc,
                1, CH4_SPMD::Chombo_MPI::comm, Rptr );
    }

    if (m_buff->m_fromMe.size() > 0)
    {
      fdata.resize(m_buff->m_fromMe.size());
      fdata[0]=m_buff->m_fromMe[0].size;
      m_buff->m_sendRequests.resize(numProc()-1);
      m_buff->m_sendStatus.resize(numProc()-1);
      MPI_Request* Rptr = &(m_buff->m_sendRequests[0]);

      unsigned int lastProc = m_buff->m_fromMe[0].procID;
      int messageSize = 1;
      unsigned long * dataPtr = &(fdata[0]);
      unsigned int i = 1;
      for (;i<m_buff->m_fromMe.size(); ++i)
      {
        fdata[i]    = m_buff->m_fromMe[i].size;
        CopierBuffer::bufEntry& b = m_buff->m_fromMe[i];
        if (b.procID == lastProc)
          messageSize++;
        else
        {
          MPI_Isend(dataPtr, messageSize, MPI_UNSIGNED_LONG, lastProc,
                    1, CH4_SPMD::Chombo_MPI::comm, Rptr);

          Rptr++;
          lastProc = b.procID;
          messageSize = 1;
          dataPtr = &(fdata[i]);
          scount++;
        }
      }

      MPI_Isend(dataPtr, messageSize, MPI_UNSIGNED_LONG, lastProc,
                1, CH4_SPMD::Chombo_MPI::comm, Rptr);
    }

    if (m_buff->m_toMe.size() > 0)
    {

      int result = MPI_Waitall(count, &(m_buff->m_receiveRequests[0]), &(m_buff->m_receiveStatus[0]));
      if (result != MPI_SUCCESS)
      {
        MayDay::Error("First pass of two-phase communication failed");
      }

      for (unsigned int i=0; i<m_buff->m_toMe.size(); ++i)
      {
        CH_assert(tdata[i] != ULONG_MAX);
        m_buff->m_toMe[i].size = tdata[i];
        recBufferSize+= tdata[i];
      }
    }

    if (m_buff->m_fromMe.size() > 0)
    {

      int result = MPI_Waitall(scount, &(m_buff->m_sendRequests[0]), &(m_buff->m_sendStatus[0]));
      if (result != MPI_SUCCESS)
      {
        MayDay::Error("First pass of two-phase communication failed");
      }

    }
  }

  // allocate send and receveive buffer space.

  if (sendBufferSize > m_buff->m_sendcapacity)
  {
    free((m_buff->m_sendbuffer));

    m_buff->m_sendbuffer = NULL;
    if (a_printStuff)
    {
      pout()<<"malloc send buffer "<<sendBufferSize<<std::endl;
    }
      
    (m_buff->m_sendbuffer) = malloc(sendBufferSize);

    if ((m_buff->m_sendbuffer) == NULL)
    {
      MayDay::Error("Out of memory in BoxLayoutData::allocatebuffers");
    }
    m_buff->m_sendcapacity = sendBufferSize;
  }

  if (recBufferSize > m_buff->m_reccapacity)
  {
    free(m_buff->m_recbuffer);

    m_buff->m_recbuffer = NULL;
    if (a_printStuff)
    {
      pout()<<"malloc receive buffer "<<recBufferSize<<std::endl;
    }

    m_buff->m_recbuffer = malloc(recBufferSize);

    if (m_buff->m_recbuffer == NULL)
    {
      MayDay::Error("Out of memory in BoxLayoutData::allocatebuffers");
    }

    m_buff->m_reccapacity = recBufferSize;
  }
  if(a_printStuff)
  {
    pout()<< "bld::allocateBuffers:" << endl;
    for (int i=0; i<m_buff->m_fromMe.size(); i++)
    pout()<<m_buff->m_fromMe[i].item->toRegion<<"{"<<m_buff->m_fromMe[i].procID<<"}"<<" ";
    pout() <<"::::";
    for (int i=0; i<m_buff->m_toMe.size(); i++)
    pout()<<m_buff->m_toMe[i].item->toRegion<<"{"<<m_buff->m_toMe[i].procID<<"}"<<" ";
    pout() << endl;
  }

  char* nextFree = (char*)(m_buff->m_sendbuffer);
  if (m_buff->m_fromMe.size() > 0)
  {
    for (unsigned int i=0; i<m_buff->m_fromMe.size(); ++i)
    {
      m_buff->m_fromMe[i].bufPtr = nextFree;
      nextFree += m_buff->m_fromMe[i].size;
    }
  }

  nextFree = (char*)m_buff->m_recbuffer;
  if (m_buff->m_toMe.size() > 0)
  {
    for (unsigned int i=0; i<m_buff->m_toMe.size(); ++i)
    {
      m_buff->m_toMe[i].bufPtr = nextFree;
      nextFree += m_buff->m_toMe[i].size;
    }
  }

  // since fromMe and toMe are sorted based on procID, messages can now be grouped
  // together on a per-processor basis.

#endif //mpi
}

template<class T>
void BoxLayoutData<T>::writeSendDataFromMeIntoBuffers(const BoxLayoutData<T>& a_src,
                                                      const Interval&     a_srcComps,
                                                      const LDOperator<T>& a_op) const
{
#ifdef CH_MPI
  CH_TIME("write Data to buffers");


  int isize = m_buff->m_fromMe.size();

  for (unsigned int i=0; i< isize; ++i)
  {
    const CopierBuffer::bufEntry& entry = m_buff->m_fromMe[i];
    a_op.linearOut(a_src[entry.item->fromIndex], entry.bufPtr,
                   entry.item->fromRegion, a_srcComps);
  }
#endif  
}

template<class T>
void BoxLayoutData<T>::postSendsFromMe() const
{
#ifdef CH_MPI  
  CH_TIME("post_Sends");
  // now we get the magic of message coalescence
  // fromMe has already been sorted in the allocateBuffers() step.

  this->m_buff->numSends = m_buff->m_fromMe.size();

  if (this->m_buff->numSends > 1)
  {
    for (unsigned int i=m_buff->m_fromMe.size()-1; i>0; --i)
    {
      if (m_buff->m_fromMe[i].procID == m_buff->m_fromMe[i-1].procID)
      {
        this->m_buff->numSends--;
        m_buff->m_fromMe[i-1].size = m_buff->m_fromMe[i-1].size + m_buff->m_fromMe[i].size;
        m_buff->m_fromMe[i].size = 0;
      }
    }
  }
  m_buff->m_sendRequests.resize(this->m_buff->numSends);
  std::list<MPI_Request> extraRequests;

  protoDeviceSynchronize(T::memTypeAllocation());

  unsigned int next=0;
  long long maxSize = 0;

  for (int i=0; i<this->m_buff->numSends; ++i)
  {
    const CopierBuffer::bufEntry& entry = m_buff->m_fromMe[next];
    char*  buffer = (char*)entry.bufPtr;
    std::size_t bsize = entry.size;
    int idtag=0;
    while (bsize > CH4_SPMD::CH_MAX_MPI_MESSAGE_SIZE)
    {
      extraRequests.push_back(MPI_Request());
      {
        //CH_TIME("MPI_Isend");
        MPI_Isend(buffer, CH4_SPMD::CH_MAX_MPI_MESSAGE_SIZE, MPI_BYTE, entry.procID,
                  idtag, CH4_SPMD::Chombo_MPI::comm, &(extraRequests.back()));
      }
      maxSize = CH4_SPMD::CH_MAX_MPI_MESSAGE_SIZE;
      bsize -= CH4_SPMD::CH_MAX_MPI_MESSAGE_SIZE;
      buffer+=CH4_SPMD::CH_MAX_MPI_MESSAGE_SIZE;
      idtag++;
    }
    {
      //CH_TIME("MPI_Isend");
      MPI_Isend(buffer, bsize, MPI_BYTE, entry.procID,
                idtag, CH4_SPMD::Chombo_MPI::comm, &(m_buff->m_sendRequests[i]));
    }
    maxSize = std::max<long long>(bsize, maxSize);
    ++next;
    while (next < m_buff->m_fromMe.size() && m_buff->m_fromMe[next].size == 0) ++next;
  }
  for (std::list<MPI_Request>::iterator it = extraRequests.begin(); it != extraRequests.end(); ++it)
  {
    m_buff->m_sendRequests.push_back(*it);
  }
  this->m_buff->numSends = m_buff->m_sendRequests.size();

  CH4_SPMD::CH_MaxMPISendSize = std::max<long long>(CH4_SPMD::CH_MaxMPISendSize, maxSize);
#endif  
}

template<class T>
void BoxLayoutData<T>::postReceivesToMe() const
{
  CH_TIME("post_Receives");
#ifndef MPI
  this->m_buff->numReceives = m_buff->m_toMe.size();

  if (this->m_buff->numReceives > 1)
  {
    for (unsigned int i=m_buff->m_toMe.size()-1; i>0; --i)
    {
      if (m_buff->m_toMe[i].procID == m_buff->m_toMe[i-1].procID)
      {
        this->m_buff->numReceives--;
        m_buff->m_toMe[i-1].size += m_buff->m_toMe[i].size;
        m_buff->m_toMe[i].size = 0;
      }

    }
  }
  m_buff->m_receiveRequests.resize(this->m_buff->numReceives);
  std::list<MPI_Request> extraRequests;
  unsigned int next=0;
  long long maxSize = 0;

  for (int i=0; i<this->m_buff->numReceives; ++i)
  {
    const CopierBuffer::bufEntry& entry = m_buff->m_toMe[next];
    char*  buffer = (char*)entry.bufPtr;
    size_t bsize = entry.size;
    int idtag=0;
    while (bsize > CH4_SPMD::CH_MAX_MPI_MESSAGE_SIZE)
    {
      extraRequests.push_back(MPI_Request());
      {
        //CH_TIME("MPI_Irecv");
        MPI_Irecv(buffer, CH4_SPMD::CH_MAX_MPI_MESSAGE_SIZE, MPI_BYTE, entry.procID,
                  idtag, CH4_SPMD::Chombo_MPI::comm, &(extraRequests.back()));
      }
      maxSize = CH4_SPMD::CH_MAX_MPI_MESSAGE_SIZE;
      bsize -= CH4_SPMD::CH_MAX_MPI_MESSAGE_SIZE;
      buffer+=CH4_SPMD::CH_MAX_MPI_MESSAGE_SIZE;
      idtag++;
    }
    {
      //CH_TIME("MPI_Irecv");
      MPI_Irecv(buffer, bsize, MPI_BYTE, entry.procID,
                idtag, CH4_SPMD::Chombo_MPI::comm, &(m_buff->m_receiveRequests[i]));
    }
    ++next;
    maxSize = std::max<long long>(bsize, maxSize);
    while (next < m_buff->m_toMe.size() && m_buff->m_toMe[next].size == 0) ++next;
  }
  for (std::list<MPI_Request>::iterator it = extraRequests.begin(); it != extraRequests.end(); ++it)
  {
    m_buff->m_receiveRequests.push_back(*it);
  }
  this->m_buff->numReceives = m_buff->m_receiveRequests.size();

  CH4_SPMD::CH_MaxMPIRecvSize = std::max<long long>(CH4_SPMD::CH_MaxMPIRecvSize, maxSize);
#endif  
}
///
template<class T>
void
BoxLayoutData<T>::
unpackReceivesToMe(const Interval&      a_destComps,
                   const LDOperator<T>& a_op)
{
#ifdef CH_MPI  
  CH_TIME("unpack_messages");

  if (this->m_buff->numReceives > 0)
  {
    m_buff->m_receiveStatus.resize(this->m_buff->numReceives);
    int result;
    {
      CH_TIME("MPI_Waitall");
      result = MPI_Waitall(this->m_buff->numReceives, &(m_buff->m_receiveRequests[0]),
                           &(m_buff->m_receiveStatus[0]));
    }
    if (result != MPI_SUCCESS)
    {
      //hell if I know what to do about failed messaging here
      //maybe a mayday::warning?
    }

    int isize = m_buff->m_toMe.size();

    for (unsigned int i=0; i< isize; ++i)
    {
      const CopierBuffer::bufEntry& entry = m_buff->m_toMe[i];
      a_op.linearIn(this->operator[](entry.item->toIndex), entry.bufPtr,
                    entry.item->toRegion, a_destComps);
    }

  this->m_buff->numReceives = 0;
  }

#endif  
}

template<class T>
void BoxLayoutData<T>::unpackReceivesToMe_append(LayoutData<Vector<RefCountedPtr<T> > >& a_dest,
                                                 const Interval&   a_destComps,
                                                 int ncomp,
                                                 const DataFactory<T>& factory,

                                                 const LDOperator<T>& a_op) 
{
#ifdef CH_MPI  
  if (this->m_buff->numReceives > 0)
  {
    m_buff->m_receiveStatus.resize(this->m_buff->numReceives);
    int result;
    {
      CH_TIME("MPI_Waitall");
      result = MPI_Waitall(this->m_buff->numReceives, &(m_buff->m_receiveRequests[0]),
                           &(m_buff->m_receiveStatus[0]));
    }
    if (result != MPI_SUCCESS)
    {
      //hell if I know what to do about failed messaging here
    }
    int isize = m_buff->m_toMe.size();

    for (int i=0; i< isize; ++i)
    {
      const CopierBuffer::bufEntry& entry = m_buff->m_toMe[i];
      const MotionItem& item = *(entry.item);
      RefCountedPtr<T> newT( factory.create(item.toRegion, ncomp, item.toIndex) );;

      a_op.linearIn(*newT, entry.bufPtr, item.toRegion, a_destComps);
      a_dest[item.toIndex].push_back(newT);
    }
  }

  this->m_buff->numReceives = 0;
#endif  
}


template <class T>
void BoxLayoutData<T>::generalCopyTo(const BoxLayout& a_destGrids,
                                     LayoutData<Vector<RefCountedPtr<T> > >& a_dest,
                                     const Interval& a_srcComps,
                                     const ProblemDomain& a_domain,
                                     const Copier& a_copier,
                                     const DataFactory<T>& factory) const
{

  CH_assert(T::preAllocatable() == 0);
  a_dest.define(a_destGrids);

  LDOperator<T> a_op;

  int ncomp = a_srcComps.size();
  Interval destComps(0, ncomp-1);
  allocateBuffers(*this,  a_srcComps,
                  *this,  destComps,
                  a_copier, a_op);

  writeSendDataFromMeIntoBuffers(*this, a_srcComps, a_op);

  // If there is nothing to recv/send, don't go into these functions
  // and allocate memory that will not be freed later.  (ndk)
  // The #ifdef CH_MPI is for the m_buff->m_toMe and m_buff->m_fromMe
#ifdef CH_MPI
  this->m_buff->numReceives = m_buff->m_toMe.size();
  if (this->m_buff->numReceives > 0)
  {
    postReceivesToMe(); // all non-blocking
  }

  this->m_buff->numSends = m_buff->m_fromMe.size();
  if (this->m_buff->numSends > 0)
  {
    postSendsFromMe();  // all non-blocking
  }
#endif

  // perform local copy
  CopyIterator it(a_copier, CopyIterator::LOCAL);
  int items=it.size();

//brian says this does not need conditionals because everyone is getting different buffers
  // NOT thread-safe, because of a_dest[item.toIndex].push_back(newT);

  for(int i=0; i<items; ++i)
  {
    const MotionItem& item = it[i];
    RefCountedPtr<T> newT( factory.create(item.toRegion, ncomp, item.toIndex) );

    a_op.op(*newT, item.fromRegion,
            destComps,
            item.toRegion,
            this->operator[](item.fromIndex),
            a_srcComps);
    a_dest[item.toIndex].push_back(newT);
  }
  // }
  // Uncomment and Move this out of unpackReceivesToMe()  (ndk)
  completePendingSends(); // wait for sends from possible previous operation

  unpackReceivesToMe_append(a_dest, destComps, ncomp, factory, a_op); // nullOp in uniprocessor mode
}

template <class T>
void BoxLayoutData<T>::generalCopyTo(const BoxLayout& a_destGrids,
                                     LayoutData<Vector<RefCountedPtr<T> > >& a_dest,
                                     const Interval& a_srcComps,
                                     const ProblemDomain& a_domain,
                                     const DataFactory<T>& factory) const
{
  Copier copier;
  copier.define(this->m_boxLayout, a_destGrids, a_domain, IntVect::Zero);

  generalCopyTo(a_destGrids, a_dest, a_srcComps, a_domain, copier, factory);
}

template <class T>
void BoxLayoutData<T>::addTo(const Interval& a_srcComps,
                             BoxLayoutData<T>& a_dest,
                             const Interval& a_destComps,
                             const ProblemDomain& a_domain) const
{
  Copier copier;
  copier.define(this->m_boxLayout, a_dest.m_boxLayout, a_domain, IntVect::Zero);
  addTo(a_srcComps, a_dest, a_destComps, a_domain, copier);
}

template <class T>
class LDaddOp : public LDOperator<T>
{
public:
  virtual void op(T& dest,
                  const Box& RegionFrom,
                  const Interval& Cdest,
                  const Box& RegionTo,
                  const T& src,
                  const Interval& Csrc) const
  {
    dest.plus(src, RegionFrom, RegionTo, Csrc.begin(), Cdest.begin(), Cdest.size());
  }
  virtual void linearIn(T& arg,  void* buf, const Box& R,
                        const Interval& comps) const
  {
    MayDay::Error("not implemented");
//    Real* buffer = (Real*)buf;
//
//    ForAllXBNNnoindx(Real, arg, R, comps.begin(), comps.size())
//      {
//        argR+=*buffer;
//        buffer++;
//      } EndFor

  }


};

template <class T>
void BoxLayoutData<T>::addTo(const Interval& a_srcComps,
                             BoxLayoutData<T>& a_dest,
                             const Interval& a_destComps,
                             const ProblemDomain& a_domain,
                             const Copier& a_copier) const
{
  CH_TIME("addTo");
  LDaddOp<T> addOp;
  addToBegin(a_srcComps, a_dest, a_destComps, a_copier);
  a_dest.addToEnd(a_destComps);
  //makeItSoBegin(a_srcComps, *this, a_dest, a_destComps, a_copier, addOp);
  //makeItSoLocalCopy(a_srcComps, *this, a_dest, a_destComps, a_copier, addOp);
  //a_dest.makeItSoEnd(a_destComps, addOp);
  //makeItSo(a_srcComps, *this, a_dest, a_destComps, a_copier, addOp);
}


template <class T>
void BoxLayoutData<T>::addToBegin(const Interval& a_srcComps,
                                  BoxLayoutData<T>& a_dest,
                                  const Interval& a_destComps,
                                  const Copier& a_copier) const
{
  CH_TIME("addToBegin");
  LDaddOp<T> addOp;
  makeItSoBegin(a_srcComps, *this, a_dest, a_destComps, a_copier, addOp);
  makeItSoLocalCopy(a_srcComps, *this, a_dest, a_destComps, a_copier, addOp);

}


template <class T>
void BoxLayoutData<T>::addToEnd(const Interval& a_destComps)
{
  CH_TIME("addToEnd");
  LDaddOp<T> addOp;

  makeItSoEnd(a_destComps, addOp);
 
}



#include "Chombo_NamespaceFooter.H"
#endif
