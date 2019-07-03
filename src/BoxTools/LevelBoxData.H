#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _LEVELBoxDATA_H____
#define _LEVELBoxDATA_H____

#include "LevelData.H"
#include "Proto.H"
#include "NamespaceHeader.H"

using Proto::BoxData;
///
/**
*/

template<int  NUMCOMPS>
class LevelBoxData: public LevelData<BoxData<Real, NUMCOMPS> >
{
public:

  typedef BoxLayoutData<BoxData<Real, NUMCOMPS> > boxlayoutdata_t;

  ///
  LevelBoxData(): LevelData<BoxData<Real, NUMCOMPS> >()
  {
  }

  ///
  LevelBoxData(const DisjointBoxLayout & a_dp, 
               const IntVect           & a_ghost   = IntVect::Zero)
    : LevelData<BoxData<Real, NUMCOMPS> >(a_dp, NUMCOMPS, a_ghost, DefaultBoxDataFactory<NUMCOMPS>())
  {
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
  }


  ///
  void copyTo(const Interval &            a_srcComps,
              boxlayoutdata_t &           a_dst,
              const Interval &            a_dstComps,
              const Copier   &            a_copier) const
  {
    LevelData<FArrayBox>     hostdataSrc(this->disjointBoxLayout(), NUMCOMPS, this->m_ghost);
    BoxLayoutData<FArrayBox> hostdataDst(a_dst.boxLayout()        , NUMCOMPS               );

    copyToHost(hostdataSrc, *this);
    copyToHost(hostdataDst,  a_dst);
    
    hostdataSrc.copyTo(a_srcComps, hostdataDst, a_dstComps,a_copier);
    
    copyToDevice(a_dst, hostdataDst);
  }


  /// 
  virtual void exchange(const Copier& copier)
  {
    LevelData<FArrayBox> hostdata(this->disjointBoxLayout(), NUMCOMPS, this->m_ghost);
    copyToHost(hostdata, *this);
    hostdata.exchange(copier);
    copyToDevice(*this, hostdata);
  }
    



  ///
  static 
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
  static 
  void copyToDevice(BoxLayoutData<BoxData<Real, NUMCOMPS> > & a_deviDat,
                    const BoxLayoutData<FArrayBox>          & a_hostDat) 
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
