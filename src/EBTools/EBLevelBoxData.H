#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _EBLEVELBoxDATA_H____
#define _EBLEVELBoxDATA_H____

#include "LevelData.H"
#include "Proto.H"
#include "AMRIO.H"
#include "ProtoInterface.H"
#include "EBProto.H"
#include "NamespaceHeader.H"

using Proto::BoxData;
using Proto::CENTERING;
using Proto::EBBoxData;
using Proto::EBHostData;
using Proto::EBGraph;
///
/**
*/
template <CENTERING cent, int NUMCOMPS> class EBBoxDataFactory : public DataFactory<EBBoxData<cent, Real, NUMCOMPS> >
{
public:
  /// factory function.  creates a new 'T' object
  /** creates a new 'T' object and returns a pointer to it.  Responsiblitly
      for calling operator 'delete' on this pointer is passed to the user. */
  virtual EBBoxData<cent, Real, NUMCOMPS>* create(const Box& box, int ncomps, const DataIndex& a_datInd) const
  {
    Proto::Box bx = ProtoCh::getProtoBox(box);
    return new EBBoxData<cent, Real, NUMCOMPS>(bx, (*m_graphs)[a_datInd]);
  }
  
  EBBoxDataFactory(shared_ptr<LevelData<EBGraph>  > a_graphs)
  {
    m_graphs = a_graphs;
  }
    

private:
  EBBoxDataFactory();
  shared_ptr<LevelData<EBGraph>  > m_graphs;
};


template <CENTERING cent, int NUMCOMPS> class EBHostDataFactory : public DataFactory<EBHostData<cent, Real, NUMCOMPS> >
{
public:
  /// factory function.  creates a new 'T' object
  /** creates a new 'T' object and returns a pointer to it.  Responsiblitly
      for calling operator 'delete' on this pointer is passed to the user. */
  virtual EBHostData<cent, Real, NUMCOMPS>* create(const Box& box, int ncomps, const DataIndex& a_datInd) const
  {
    Proto::Box bx = ProtoCh::getProtoBox(box);
    return new EBHostData<cent, Real, NUMCOMPS>(bx, (*m_graphs)[a_datInd]);
  }
  
  EBHostDataFactory(shared_ptr<LevelData<EBGraph>  > a_graphs)
  {
    m_graphs = a_graphs;
  }
    

private:
  EBHostDataFactory();
  shared_ptr<LevelData<EBGraph>  > m_graphs;
};

///
/**
*/
template<CENTERING cent, int  NUMCOMPS>
class EBLevelBoxData: public LevelData<EBBoxData<cent, Real, NUMCOMPS> >
{
public:

  ///
  EBLevelBoxData(): EBLevelData<BoxData<Real, NUMCOMPS> >()
  {
  }


  ///
  LevelBoxData(const DisjointBoxLayout        & a_dp, 
               const IntVect                  & a_ghost,
               shared_ptr<LevelData<EBGraph>  > a_graphs)
    : LevelData<EBBoxData<cent, Real, NUMCOMPS> >(a_dp, NUMCOMPS, a_ghost, EBBoxDataFactory<cent, NUMCOMPS>(a_graphs)), m_graphs(a_graphs)
  {
  }

  ///
  virtual ~LevelBoxData()
  {
  }

  ///
  void define(const DisjointBoxLayout& a_dp, 
              const IntVect          & a_ghost,
               shared_ptr<LevelData<EBGraph>  > a_graphs)
  {
    m_graphs = a_graphs;
    LevelData<EBBoxData<cent, Real, NUMCOMPS> >::define(a_dp, NUMCOMPS, a_ghost, EBBoxDataFactory<cent, NUMCOMPS>(a_graphs));
  }


  ///
  void copyTo(const Interval                      &            a_srcComps,
              LevelData<BoxData<Real, NUMCOMPS> > &            a_dst,
              const Interval                      &            a_dstComps,
              const Copier                        &            a_copier) const
  {
    LevelData<EBHostData<cent, NUMCOMPS> > hostdataSrc(this->disjointBoxLayout(), NUMCOMPS,  this->m_ghost   , EBHostDataFactory<cent, NUMCOMPS>(m_graphs));
    LevelData<EBHostData<cent, NUMCOMPS> > hostdataDst(a_dst.disjointBoxLayout(), NUMCOMPS, a_dst.ghostVect(), EBHostDataFactory<cent, NUMCOMPS>(a_dst.m_graphs));

    copyToHost(hostdataSrc, *this);
    copyToHost(hostdataDst,  a_dst);
    
    hostdataSrc.copyTo(a_srcComps, hostdataDst, a_dstComps,a_copier);
    
    copyToDevice(a_dst, hostdataDst);
  }


  /// 
  virtual void exchange(const Copier& copier)
  {
    LevelData<EBHostData<cent, NUMCOMPS> > hostdata(this->disjointBoxLayout(), NUMCOMPS,  this->m_ghost   , EBHostDataFactory<cent, NUMCOMPS>(m_graphs));
    copyToHost(hostdata, *this);
    hostdata.exchange();
    copyToDevice(*this, hostdata);
  }
    

  void 
  writeToFileHDF5(const string& a_filename) const
  {
//    LevelData<FArrayBox> hostdata(this->disjointBoxLayout(), NUMCOMPS, this->m_ghost);
//    copyToHost(hostdata, *this);
//    writeLevelname(&hostdata, a_filename.c_str());
  }

  ///
  static 
  void copyToHost(LevelData<      EBHostData<cent, NUMCOMPS> > & a_hostDat,
                  const LevelData< EBBoxData<cent, NUMCOMPS> > & a_deviDat)
  {
    LevelData< EBBoxData<cent, NUMCOMPS> > deviDat = const_cast<LevelData< EBBoxData<cent, NUMCOMPS> > &> (a_deviDat);
    DataIterator dit = a_hostDat.dataIterator();
    for(int ibox = 0;  ibox < dit.size(); ibox++)
    {
      EBHostData<cent, NUMCOMPS> & hostfab = a_hostDat[dit[ibox]];
      EBBoxData< cent, NUMCOMPS> & devifab =   deviDat[dit[ibox]];
      for(unsigned int iwhichfab = 0; iwhichfab < 2; iwhichfab++)
      {
        size_t fabsizeHost, fabsizeDevi; 
        Real*  dataptrHost, dataptrDevi;
        hostfab.dataInfo(dataptrHost, fabsizeHost, iwhichfab);
        devifab.dataInfo(dataptrDevi, fabsizeDevi, iwhichfab);

        size_t datSize = sizeof(Real)*fabsizedevi;

#ifdef PROTO_CUDA
        cudaMemcpy(dataptrHost, dataptrDevi, datSize, cudaMemcpyDeviceToHost);
        cudaError err = cudaGetLastError();
        if (err != cudaSuccess)
        {
          fprintf(stderr, "cudaGetLastError() failed at %s:%i : %s\n",
                  __FILE__, __LINE__, cudaGetErrorString(err));
        }
#else
        memcpy(dataptrHost, dataptrDevi, datSize);
#endif  
      }
    }
  }

  ///
  ///
  static 
  void copyToDevice(LevelData<        EBBoxData<cent, NUMCOMPS> > & a_deviDat,
                    const LevelData< EBHostData<cent, NUMCOMPS> > & a_hostDat)
  {
    LevelData< EBBoxData<cent, NUMCOMPS> > deviDat = const_cast<LevelData< EBBoxData<cent, NUMCOMPS> > &> (a_deviDat);
    DataIterator dit = a_hostDat.dataIterator();
    for(int ibox = 0;  ibox < dit.size(); ibox++)
    {
      EBHostData<cent, NUMCOMPS> & hostfab = a_hostDat[dit[ibox]];
      EBBoxData< cent, NUMCOMPS> & devifab =   deviDat[dit[ibox]];
      for(unsigned int iwhichfab = 0; iwhichfab < 2; iwhichfab++)
      {
        size_t fabsizeHost, fabsizeDevi; 
        Real*  dataptrHost, dataptrDevi;
        hostfab.dataInfo(dataptrHost, fabsizeHost, iwhichfab);
        devifab.dataInfo(dataptrDevi, fabsizeDevi, iwhichfab);

        size_t datSize = sizeof(Real)*fabsizedevi;

#ifdef PROTO_CUDA
        cudaMemcpy(dataptrDevi, dataptrHost, datSize, cudaMemcpyHostToDevice);
        cudaError err = cudaGetLastError();
        if (err != cudaSuccess)
        {
          fprintf(stderr, "cudaGetLastError() failed at %s:%i : %s\n",
                  __FILE__, __LINE__, cudaGetErrorString(err));
        }
#else
        memcpy(dataptrDevi, dataptrHost, datSize);
#endif  
      }
    }
  }


  private:
    shared_ptr<LevelData<EBGraph>  > m_graphs;
};


#include "NamespaceFooter.H"

#endif
