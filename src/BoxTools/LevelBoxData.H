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
#include "AMRIO.H"
#include "ProtoInterface.H"
#include "NamespaceHeader.H"

using Proto::BoxData;
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
    Proto::Box bx = ProtoCh::getProtoBox(box);
    return new BoxData<Real, NUMCOMPS>(bx);
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
  void copyTo(const Interval                      &            a_srcComps,
              LevelData<BoxData<Real, NUMCOMPS> > &           a_dst,
              const Interval &            a_dstComps,
              const Copier   &            a_copier) const
  {
    LevelData<FArrayBox> hostdataSrc(this->disjointBoxLayout(), NUMCOMPS,  this->m_ghost);
    LevelData<FArrayBox> hostdataDst(a_dst.disjointBoxLayout(), NUMCOMPS, a_dst.ghostVect());

    copyToHost(hostdataSrc, *this);
    copyToHost(hostdataDst,  a_dst);
    
    hostdataSrc.copyTo(a_srcComps, hostdataDst, a_dstComps,a_copier);
    
    copyToDevice(a_dst, hostdataDst);
  }


  /// 
  virtual void exchange(const Copier& copier)
  {
    LevelData<BoxData<Real, NUMCOMPS> > devdata(this->disjointBoxLayout(), NUMCOMPS, this->m_ghost);
    //copyToHost(hostdata, *this);
    devdata.exchange();
    //copyToDevice(*this, hostdata);
  }
    

  void 
  writeToFileHDF5(const string& a_filename) const
  {
#ifdef CH_USE_HDF5
    LevelData<FArrayBox> hostdata(this->disjointBoxLayout(), NUMCOMPS, this->m_ghost);
    copyToHost(hostdata, *this);
    writeLevelname(&hostdata, a_filename.c_str());
#else
    printf("hdf5 not enabled");
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
