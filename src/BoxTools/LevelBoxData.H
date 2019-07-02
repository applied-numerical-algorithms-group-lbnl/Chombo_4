#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _LEVELBoxDATA_H__
#define _LEVELBoxDATA_H__

#include "LevelData.H"
#include "Proto.H"
#include "NamespaceHeader.H"

using Proto::BoxData;
///
/**
*/
template<int NUMCOMPS>
class LevelBoxData: public LevelData<BoxData<REAL, NUMCOMPS> >
{
public:

  typedef DefaultDataFactory<BoxData<REAL, NUMCOMPS> > factory_t;
  typedef BoxLayoutData<BoxData<REAL, NUMCOMPS> > boxlayoutdata_t;
  ///
  LevelBoxData(): LevelData<BoxData<REAL, NUMCOMPS> >()
  {
  }

  ///
  LevelBoxData(const DisjointBoxLayout & a_dp, 
               const IntVect           & a_ghost   = IntVect::Zero,
               const DataFactory<T>    & a_factory = factory_t())
    : LevelData<BoxData<REAL, NUMCOMPS> >(a_dp, NUMCOMPS, a_ghost, a_factory)
  {
  }

  ///
  virtual ~LevelBoxData()
  {
  }

  ///
  virtual void define(const DisjointBoxLayout& a_dp, 
                      const IntVect          & a_ghost = IntVect::Zero,
                      const DataFactory<T>   & a_factory = factory_t())

  {
    LevelData<BoxData<REAL, NUMCOMPS> >::define(a_dp, NUMCOMPS, a_ghost, a_factory);
  }


  ///
  virtual void copyTo(const Interval &            a_srcComps,
                      boxlayoutdata_t &           a_dst,
                      const Interval &            a_dstComps,
                      const Copier   &            a_copier,
                      const LDOperator<T>&        a_op = LDOperator<T>()) const
  {

    LevelData<FArrayBox>     hostdataSrc(      m_disjointBoxLayout, NUMCOMPS,        m_ghost, factory_t());
    BoxLayoutData<FArrayBox> hostdataDst(a_dst.m_boxLayout        , NUMCOMPS,                 factory_t());

    copyToHost( hostdataSrc, *this);
    copyToHost(hostdataDst,  a_dst);
    
    hostdataSrc.copyTo(a_srcComps, hostdataDst, a_dstComps,a_copier, a_op);
    
    copyToDevice(a_dst, hostDataDst);
  }


  /// 
  virtual void exchange(const Copier& copier)
  {
    LevelData<FArrayBox> hostdata(m_disjointBoxLayout, NUMCOMPS,        m_ghost, factory_t());
    copyToHost(hostdata, *this);
    hostdata.exchange(copier);
    copyToDevice(*this, hostdata);
  }
    


private:

  ///
  void copyToHost(BoxLayoutData<FArrayBox>                      & a_hostDat,
                  const BoxLayoutData<BoxData<Real, NUMCOMPS> > & a_deviDat)
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
      cudaMemcpy(hostPtr, deviPtr, datSize, cudaMemcpyDeviceToHost);
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
  void copyToDevice(BoxLayoutData<BoxData<Real, NUMCOMPS> > & a_deviceDat,
                    const BoxLayoutData<FArrayBox>          & a_hostDat)
  {
    DataIterator dit = a_hostDat.dataIterator();
    for(int ibox = 0;  ibox < dit.size(); ibox++)
    {
      FArrayBox & hostfab                     = a_hostDat[dit[ibox]];
      const BoxData<Real, NUMCOMPS> & devifab = a_deviDat[dit[ibox]];
      size_t datSize = sizeof(Real)*devifab.size();
      const Real* hostPtr = hostfab.dataPtr(0);
      Real*       deviPtr = devifab.dataPtr(0);
#ifdef PROTO_CUDA
      cudaMemcpy(deviPtr, hostPtr, datSize, cudaMemcpyHostToDevice);
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
};


#include "NamespaceFooter.H"

#endif
