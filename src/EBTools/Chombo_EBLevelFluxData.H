#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _EBLevelFluxData_H____
#define _EBLevelFluxData_H____

#include "Proto.H"
#include "EBProto.H"
#include "Chombo_LevelData.H"
#include "Chombo_AMRIO.H"
#include "Chombo_ProtoInterface.H"
#include "Chombo_EBChombo.H"
#include "Chombo_ProtoFactories.H"
#include "Chombo_LevelBoxData.H"
#include "Chombo_NamespaceHeader.H"


template<unsigned int srcComp, unsigned int dstComp>
class ScalarFluxAliasFactory: public DataFactory< EBFluxData<Real, dstComp> >
{
public:

  ScalarFluxAliasFactory(LevelData<EBFluxData<Real, srcComp> >* a_source, 
                         unsigned int                           a_whichcomp)
  {
    m_source    = a_source;
    m_whichcomp = a_whichcomp;
  }
  virtual ~ScalarFluxAliasFactory()
  {
  }

  /// factory function.  creates a new 'T' object
  /** creates a new 'T' object and returns a pointer to it.
   */
  virtual EBFluxData<Real, dstComp> * create(const Box& a_box, int ncomps, const DataIndex& a_datInd) const
  {
    EBFluxData<Real, srcComp>& vecholder = (*m_source)[a_datInd];
    EBFluxData<Real, dstComp>* retval = new EBFluxData<Real, dstComp>();
    //aliasing define
    retval->define(vecholder, m_whichcomp);
    return retval;
  }

  virtual bool callDelete() const
  {
    return false;
  }

  //return true if copy and construction are thread safe
  virtual bool threadSafe() const
  {
    return false;
  }

private:
  ScalarFluxAliasFactory();
  LevelData<EBFluxData<Real, DIM> >* m_source;
  unsigned int                       m_whichcomp;

};
///
/**
*/
template<int  NUMCOMPS>
class EBLevelFluxData: public LevelData<EBFluxData<Real, NUMCOMPS> >
{
public:
  //shutting up the compilers
  using LevelData<EBFluxData<Real, NUMCOMPS> >::define;


  typedef GraphConstructorFactory<EBHostFluxData<Real, NUMCOMPS> > hostfactory_t;
  typedef GraphConstructorFactory<    EBFluxData<Real, NUMCOMPS> > devifactory_t;

  ///
  EBLevelFluxData(): LevelData<EBFluxData<Real, NUMCOMPS> >()
  {
  }


  ///
  EBLevelFluxData(const DisjointBoxLayout        & a_dp, 
                  const IntVect                  & a_ghost,
                  shared_ptr<LevelData<EBGraph>  > a_graphs)
    : LevelData<EBFluxData<Real, NUMCOMPS> >(a_dp, NUMCOMPS, a_ghost, devifactory_t(a_graphs)), m_graphs(a_graphs)
  {
  }

  ///
  virtual ~EBLevelFluxData()
  {
  }


  ///
  virtual
  void define(const DisjointBoxLayout& a_dp, 
              const IntVect          & a_ghost,
              shared_ptr<LevelData<EBGraph>  > a_graphs)
  {
    m_graphs = a_graphs;
    LevelData<EBFluxData<Real, NUMCOMPS> >::define(a_dp, NUMCOMPS, a_ghost, devifactory_t(a_graphs));
  }

  ///aliasing define
  template <unsigned int srcComp>
  void define(EBLevelFluxData<srcComp>   &     a_input, 
              unsigned int                     a_comp,
              shared_ptr<LevelData<EBGraph>  > a_graphs)
  {
    m_graphs = a_graphs;
    auto dp    = a_input.disjointBoxLayout();
    auto ghost = a_input.ghostVect();
    m_graphs = a_graphs;
    LevelData<EBFluxData<Real, srcComp> >* inputPtr = static_cast<LevelData<EBFluxData<Real, DIM> > * >(&a_input);
    ScalarFluxAliasFactory<srcComp, NUMCOMPS> fact(inputPtr, a_comp);
    LevelData<EBFluxData<Real, NUMCOMPS> >::define(dp, NUMCOMPS, ghost, fact);
  }

  void setVal(Real a_value)
  {
    DataIterator dit = this->dataIterator();
    for(int ibox = 0;  ibox < dit.size(); ibox++)
    {
      (*this)[dit[ibox]].setVal(a_value);
    }
  }

  ///
  Real maxNorm(int a_icomp) const
  {
    CH_TIME("EBLevelBoxData::maxNorm");
    LevelData<EBHostFluxData<Real, NUMCOMPS> > hostdata(this->disjointBoxLayout(), NUMCOMPS,  this->m_ghost   , hostfactory_t(m_graphs));

    copyToHost(hostdata, *this);

    DataIterator dit = this->dataIterator();
    DisjointBoxLayout dbl = this->disjointBoxLayout();
    Real maxval = 0; 
    for(int ibox = 0;  ibox < dit.size(); ibox++)
    {
      const Box & grid     =         dbl[dit[ibox]];
      Bx grbx = ProtoCh::getProtoBox(grid);
      Real gridval = hostdata[dit[ibox]].maxNorm(grbx, a_icomp);
      if(gridval > maxval)
      {
        maxval = gridval;
      }
    }

#ifdef CH_MPI    
    Real sendBuf = maxval;
    int result = MPI_Allreduce(&sendBuf, &maxval, 1, MPI_CH_REAL,MPI_MAX, Chombo_MPI::comm);

    if (result != MPI_SUCCESS)
    {
      MayDay::Error("Communication error in EBLevelBoxData::maxNorm");
    }
#endif
    return maxval;
  }
  
  ///
  void copyTo(const Interval                      &            a_srcComps,
              LevelData<BoxData<Real, NUMCOMPS> > &            a_dst,
              const Interval                      &            a_dstComps,
              const Copier                        &            a_copier) const
  {
    LevelData<EBHostFluxData<Real, NUMCOMPS> > hostdataSrc(this->disjointBoxLayout(), NUMCOMPS,  this->m_ghost   , hostfactory_t(m_graphs));
    LevelData<EBHostFluxData<Real, NUMCOMPS> > hostdataDst(a_dst.disjointBoxLayout(), NUMCOMPS, a_dst.ghostVect(), hostfactory_t(a_dst.m_graphs));

    copyToHost(hostdataSrc, *this);
    copyToHost(hostdataDst,  a_dst);
    
    hostdataSrc.copyTo(a_srcComps, hostdataDst, a_dstComps,a_copier);
    
    copyToDevice(a_dst, hostdataDst);
  }


  /// 
  virtual void exchange(const Copier& copier)
  {
    LevelData<EBHostFluxData<Real, NUMCOMPS> > hostdata(this->disjointBoxLayout(), NUMCOMPS,  this->m_ghost   , hostfactory_t(m_graphs));

    copyToHost(hostdata, *this);

    hostdata.exchange(copier);

    copyToDevice(*this, hostdata);
  }
  ////mpi stuff
  static 
  void copyToHost(LevelData<      EBHostFluxData<Real, NUMCOMPS> > & a_hostDat,
                  const LevelData<    EBFluxData<Real, NUMCOMPS> > & a_deviDat)
  {
    LevelData< EBFluxData<Real, NUMCOMPS> >& deviDat = (LevelData< EBFluxData<Real, NUMCOMPS> >&) (a_deviDat);
    DataIterator dit = a_hostDat.dataIterator();
    for(int ibox = 0;  ibox < dit.size(); ibox++)
    {
      auto & hostfab = a_hostDat[dit[ibox]];
      auto & devifab =   deviDat[dit[ibox]];
      copyToHost(hostfab, devifab);
    }
  }


  /// 
  static 
  void copyToHost(EBHostFluxData<Real, NUMCOMPS> & hostfab,
                      EBFluxData<Real, NUMCOMPS> & devifab)
  {
    for(unsigned int idir = 0; idir < DIM; idir++)
    {
      for(unsigned int iwhichfab = 0; iwhichfab < 2; iwhichfab++)
      {
        size_t fabsizeHost, fabsizeDevi; 
        Real*  dataptrHost;
        Real*  dataptrDevi;
        hostfab.dataInfo(dataptrHost, fabsizeHost, iwhichfab, idir);
        devifab.dataInfo(dataptrDevi, fabsizeDevi, iwhichfab, idir);

        size_t datSize = sizeof(Real)*fabsizeDevi;

#ifdef PROTO_CUDA
        protoMemcpy(DEVICE,dataptrHost, dataptrDevi, datSize, protoMemcpyDeviceToHost);
        protoError err = protoGetLastError();
        if (err != protoSuccess)
        {
          fprintf(stderr, "protoGetLastError() failed at %s:%i : %s\n",
                  __FILE__, __LINE__, protoGetErrorString(err));
        }
#else
        memcpy(dataptrHost, dataptrDevi, datSize);
#endif  
      }
    }
  }

  
  ///
  static 
  void copyToDevice(LevelData<          EBFluxData<Real, NUMCOMPS> > & a_deviDat,
                    const LevelData<EBHostFluxData<Real, NUMCOMPS> > & a_hostDat)
  {
    auto & deviDat = const_cast<LevelData<     EBFluxData<Real, NUMCOMPS> > &> (a_deviDat);
    auto & hostDat = const_cast<LevelData< EBHostFluxData<Real, NUMCOMPS> > &> (a_hostDat);
    DataIterator dit = hostDat.dataIterator();
    for(int ibox = 0;  ibox < dit.size(); ibox++)
    {
      auto & hostfab =   hostDat[dit[ibox]];
      auto & devifab =   deviDat[dit[ibox]];
      copyToDevice(hostfab, devifab);
    }
  }

  ///
  static 
  void copyToDevice( EBHostFluxData<Real, NUMCOMPS> & hostfab,
                     EBFluxData<    Real, NUMCOMPS> & devifab)
  {
    for(unsigned int idir = 0; idir < DIM; idir++)
    {
      for(unsigned int iwhichfab = 0; iwhichfab < 2; iwhichfab++)
      {
        size_t fabsizeHost, fabsizeDevi; 
        Real*  dataptrHost, *dataptrDevi;
        hostfab.dataInfo(dataptrHost, fabsizeHost, iwhichfab, idir);
        devifab.dataInfo(dataptrDevi, fabsizeDevi, iwhichfab, idir);

        size_t datSize = sizeof(Real)*fabsizeDevi;

#ifdef PROTO_CUDA
        protoMemcpy(DEVICE,dataptrDevi, dataptrHost, datSize, protoMemcpyHostToDevice);
        protoError err = protoGetLastError();
        if (err != protoSuccess)
        {
          fprintf(stderr, "protoGetLastError() failed at %s:%i : %s\n",
                  __FILE__, __LINE__, protoGetErrorString(err));
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


#include "Chombo_NamespaceFooter.H"

#endif
