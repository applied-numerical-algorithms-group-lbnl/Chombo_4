#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _Chombo_LEVELBoxDATA_H____
#define _Chombo_LEVELBoxDATA_H____

#include "Chombo_LevelData.H"
#include "Proto.H"
#include "Chombo_AMRIO.H"
#include "Chombo_ProtoInterface.H"
#include "Chombo_NamespaceHeader.H"

using ::Proto::BoxData;
///
/**
*/
template <int NUMCOMPS> class DefaultBoxDataFactory : public DataFactory<BoxData<Real, NUMCOMPS> >
{
public:
  /// factory function.  creates a new 'T' object
  /** creates a new 'T' object and returns a pointer to it.  Responsiblitly
      for calling operator 'delete' on this pointer is passed to the user. */
  virtual BoxData<Real, NUMCOMPS>* create(const Box& box, int ncomps, const DataIndex& a_datInd) const
  {
    ::Proto::Box bx = ProtoCh::getProtoBox(box);
    return new BoxData<Real, NUMCOMPS>(bx, false);
  }
};

///
/**
*/
template<int  NUMCOMPS>
class LevelBoxData: public LevelData<BoxData<Real, NUMCOMPS> >
{

public:

  ///
  LevelBoxData(): LevelData<BoxData<Real, NUMCOMPS> >()
  {
  }


  ///
  LevelBoxData(const DisjointBoxLayout & a_dp, 
               const IntVect           & a_ghost   = IntVect::Zero)
    : LevelData<BoxData<Real, NUMCOMPS> >(a_dp, NUMCOMPS, a_ghost, DefaultBoxDataFactory<NUMCOMPS>())
  {
    if(a_ghost > IntVect::Zero)
	 m_host.define(a_dp, NUMCOMPS, a_ghost);
  }


  ///
  virtual ~LevelBoxData()
  {
  }

  ///
  void define(const DisjointBoxLayout& a_dp, 
              const IntVect          & a_ghost = IntVect::Zero)
  {
    LevelData<BoxData<Real, NUMCOMPS> >::define(a_dp, NUMCOMPS, a_ghost, DefaultBoxDataFactory<NUMCOMPS>());
    if(a_ghost > IntVect::Zero)
	 m_host.define(a_dp, NUMCOMPS, a_ghost);
  }

  void define(const LevelBoxData<NUMCOMPS>& a_lbd)
  {
    LevelData<BoxData<Real, NUMCOMPS> >::define(a_lbd, DefaultBoxDataFactory<NUMCOMPS>());
    if(a_lbd.ghostVect() > IntVect::Zero)
	 m_host.define(a_lbd.disjointBoxLayout(), NUMCOMPS, a_lbd.ghostVect());
  }

  ///
  void copyTo(const Interval                      &            a_srcComps,
              LevelData<BoxData<Real, NUMCOMPS> > &           a_dst,
              const Interval &            a_dstComps,
              const Copier   &            a_copier) 
  {
  //  LevelData<FArrayBox> hostdataSrc(this->disjointBoxLayout(), NUMCOMPS,  this->m_ghost);
  //  LevelData<FArrayBox> hostdataDst(a_dst.disjointBoxLayout(), NUMCOMPS, a_dst.ghostVect());

    // To check
    if(this->m_ghost > IntVect::Zero)
	 m_hostSrc.define(this->disjointBoxLayout(), NUMCOMPS, this->m_ghost);

    // To check
    if(a_dst.ghostVect() > IntVect::Zero)
	 m_hostDst.define(a_dst.disjointBoxLayout(), NUMCOMPS, a_dst.ghostVect());

    copyToHost(m_hostSrc, *this);
    copyToHost(m_hostDst,  a_dst);
    
    m_hostSrc.copyTo(a_srcComps, m_hostDst, a_dstComps, a_copier);
    
    copyToDevice(a_dst, m_hostDst);
  }


  /// 
  virtual void exchange(const Copier& copier)
  {
    CH_TIMERS("LevelBoxData::exchange");
    CH_TIMER("FArrayBox",t0);
    CH_TIMER("copyToHost",t1);
    CH_TIMER("exchange",t2);
    CH_TIMER("copyToDevice",t3);
    CH_START(t0);
// LevelData<FArrayBox> hostdata(this->disjointBoxLayout(), NUMCOMPS, this->m_ghost);
    CH_STOP(t0);
    CH_START(t1);
    copyToHost(m_host, *this);
    CH_STOP(t1);
    CH_START(t2);
    m_host.exchange(copier);
    CH_STOP(t2);
    CH_START(t3);
    copyToDevice(*this, m_host);
    CH_STOP(t3);
  }
    

  void 
  writeToFileHDF5(const string& a_filename) const
  {
#ifdef CH_USE_HDF5
    LevelData<FArrayBox> hostdata(this->disjointBoxLayout(), NUMCOMPS, this->m_ghost);
    copyToHost(hostdata, *this);
    writeLevelname(&hostdata, a_filename.c_str());
#else
    printf("hdf5 not enabled\n");
#endif
  }

  ///
  static 
  void copyToHost(LevelData<FArrayBox>                      & a_hostDat,
                  const LevelData<BoxData<Real, NUMCOMPS> > & a_deviDat)
  {
    DataIterator dit = a_hostDat.dataIterator();
    for(int ibox = 0;  ibox < dit.size(); ibox++)
    {
      FArrayBox & hostfab                     = a_hostDat[dit[ibox]];
      const BoxData<Real, NUMCOMPS> & devifab = a_deviDat[dit[ibox]];
      size_t datSize = sizeof(Real)*devifab.size();
      Real*       hostPtr = hostfab.dataPtr(0);
      const Real* deviPtr = devifab.dataPtr(0);
#ifdef PROTO_CUDA
      //cudaMemcpy(hostPtr, deviPtr, datSize, cudaMemcpyDeviceToHost);
      //cudaStreamSynchronize(Proto::DisjointBoxLayout::getCurrentStream());
      cudaDeviceSynchronize();
      cudaMemcpyAsync(hostPtr, deviPtr, datSize, cudaMemcpyDeviceToHost,Proto::DisjointBoxLayout::getCurrentStream());
      cudaError err = cudaGetLastError();
      if (err != cudaSuccess)
      {
        fprintf(stderr, "cudaGetLastError() failed at %s:%i : %s\n",
                __FILE__, __LINE__, cudaGetErrorString(err));
      }
#else
      memcpy(hostPtr, deviPtr, datSize);
#endif  
    }
  }

  ///
  static 
  void copyToDevice(LevelData<BoxData<Real, NUMCOMPS> > & a_deviDat,
                    const LevelData<FArrayBox>          & a_hostDat) 
  {
    DataIterator dit = a_hostDat.dataIterator();
    for(int ibox = 0;  ibox < dit.size(); ibox++)
    {
      const FArrayBox & hostfab         = a_hostDat[dit[ibox]];
      BoxData<Real, NUMCOMPS> & devifab = a_deviDat[dit[ibox]];
      size_t datSize = sizeof(Real)*devifab.size();
      const Real* hostPtr = hostfab.dataPtr(0);
      Real*       deviPtr = devifab.dataPtr(0);
#ifdef PROTO_CUDA
      cudaDeviceSynchronize();
      cudaMemcpyAsync(deviPtr, hostPtr, datSize, cudaMemcpyHostToDevice, Proto::DisjointBoxLayout::getCurrentStream());

      //cudaMemcpy(deviPtr, hostPtr, datSize, cudaMemcpyHostToDevice);
      cudaError err = cudaGetLastError();
      if (err != cudaSuccess)
      {
        fprintf(stderr, "cudaGetLastError() failed at %s:%i : %s\n",
                __FILE__, __LINE__, cudaGetErrorString(err));
      }
#else
      memcpy(deviPtr, hostPtr, datSize);
#endif  
    }
  }

private:
  LevelData<FArrayBox> m_host;
  LevelData<FArrayBox> m_hostSrc;
  LevelData<FArrayBox> m_hostDst;
};


#include "Chombo_NamespaceFooter.H"

#endif
