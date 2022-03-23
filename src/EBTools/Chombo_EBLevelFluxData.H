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
  typedef CH4_Data_Choreography::DistributedData<EBFluxData<    Real, srcComp> > devi_src_comp_t;
  typedef CH4_Data_Choreography::DistributedData<EBHostFluxData<Real, srcComp> > host_src_comp_t;
  typedef CH4_Data_Choreography::DistributedData<EBFluxData<    Real, dstComp> > devi_dst_comp_t;
  typedef CH4_Data_Choreography::DistributedData<EBHostFluxData<Real, dstComp> > host_dst_comp_t;
  typedef CH4_Data_Choreography::DistributedData<EBGraph> graph_distrib_t;

  ScalarFluxAliasFactory(devi_src_comp_t * a_source, 
                         unsigned int      a_whichcomp)
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

private:
  ScalarFluxAliasFactory();
  devi_src_comp_t* m_source;
  unsigned int     m_whichcomp;

};
///
/**
*/
template<int  NUMCOMPS>
class EBLevelFluxData: public CH4_Data_Choreography::DistributedData<EBFluxData<Real, NUMCOMPS> >
{
public:
  typedef CH4_Data_Choreography::DistributedData<    EBFluxData<Real, NUMCOMPS> > devidistrib_t;
  typedef CH4_Data_Choreography::DistributedData<EBHostFluxData<Real, NUMCOMPS> > hostdistrib_t;
  
  typedef GraphConstructorFactory<    EBFluxData<Real, NUMCOMPS> > devifactory_t;
  typedef GraphConstructorFactory<EBHostFluxData<Real, NUMCOMPS> > hostfactory_t;
  typedef CH4_Data_Choreography::DistributedData<EBGraph> graph_distrib_t;

  ///
  EBLevelFluxData(): devidistrib_t()
  {
  }


  ///
  EBLevelFluxData(const DisjointBoxLayout        & a_dp, 
                  const IntVect                  & a_ghost,
                  shared_ptr<graph_distrib_t>      a_graphs)
    : devidistrib_t(a_dp, a_ghost, devifactory_t(a_graphs)), m_graphs(a_graphs)
  {
  }

  ///
  virtual ~EBLevelFluxData()
  {
  }


  ///
  virtual
  void define(const DisjointBoxLayout    & a_dp, 
              const IntVect              & a_ghost,
              shared_ptr<graph_distrib_t>  a_graphs)
  {
    m_graphs = a_graphs;
    devidistrib_t::define(a_dp, a_ghost, devifactory_t(a_graphs));
  }

  ///aliasing define
  template <unsigned int srcComp>
  void define(EBLevelFluxData<srcComp>   &     a_input, 
              unsigned int                     a_comp,
              shared_ptr<graph_distrib_t>      a_graphs)
  {
    m_graphs = a_graphs;
    auto dp    = a_input.disjointBoxLayout();
    auto ghost = a_input.ghostVect();
    m_graphs = a_graphs;
    typedef CH4_Data_Choreography::DistributedData<EBFluxData< Real, srcComp> > srccomp_devi_t;
    srccomp_devi_t * inputPtr = static_cast<srccomp_devi_t* >(&a_input);

    ScalarFluxAliasFactory<srcComp, NUMCOMPS> fact(inputPtr, a_comp);
    devidistrib_t::define(dp, ghost, fact);
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
  
  
  ///all the cool kids use this copyTo interface
  /**
     copy over box intersections
   */
  void copyTo(EBLevelFluxData<NUMCOMPS> & a_dst, bool a_printStuff = false) const
  {
    hostdistrib_t hostdataSrc(this->disjointBoxLayout(),  this->ghostVect(), hostfactory_t(      m_graphs));
    hostdistrib_t hostdataDst(a_dst.disjointBoxLayout(),  a_dst.ghostVect(), hostfactory_t(a_dst.m_graphs));

    copyToHost(hostdataSrc, *this);
    copyToHost(hostdataDst,  a_dst);
    
    hostdataSrc.copyTo(hostdataDst, a_printStuff);
    
    copyToDevice(a_dst, hostdataDst);
  }


  /// The hip new slick exchange interface.  Use this one.
  virtual void exchange(bool a_printStuff = false)
  {
    if(a_printStuff)
    {
      pout() << "eblevelfluxdata::exchange creating host data" << endl;
    }
    hostdistrib_t hostdata(this->disjointBoxLayout(), this->ghostVect(), hostfactory_t(m_graphs));

    if(a_printStuff)
    {
      pout() << "eblevelfluxdata::exchange copying data to host" << endl;
    }
    copyToHost(hostdata, *this);

    if(a_printStuff)
    {
      pout() << "eblevelfluxdata::exchange calling DistributedData::exchange" << endl;
    }
    hostdata.exchange(a_printStuff);

    if(a_printStuff)
    {
      pout() << "eblevelfluxdata::exchange calling copyToDevice" << endl;
    }
    copyToDevice(*this, hostdata);
    if(a_printStuff)
    {
      pout() << "leaving eblevelfluxdata exchange" << endl;
      
    }
  }
  ////mpi stuff
  static 
  void copyToHost(      hostdistrib_t & a_hostDat,
                  const devidistrib_t & a_deviDat)
  {
    DataIterator dit = a_hostDat.dataIterator();
    for(int ibox = 0;  ibox < dit.size(); ibox++)
    {
      auto & hostfab = a_hostDat[dit[ibox]];
      auto & devifab = a_deviDat[dit[ibox]];
      copyToHost(hostfab, devifab);
    }
  }


  /// 
  static 
  void copyToHost(  EBHostFluxData<Real, NUMCOMPS> & a_hostfab,
                  const EBFluxData<Real, NUMCOMPS> & a_devifab)
  {
    for(unsigned int idir = 0; idir < DIM; idir++)
    {
      for(unsigned int iwhichfab = 0; iwhichfab < 2; iwhichfab++)
      {
        size_t fabsizeHost, fabsizeDevi; 
        Real*  dataptrHost;
        Real*  dataptrDevi;
        a_hostfab.dataInfo(dataptrHost, fabsizeHost, iwhichfab, idir);
        a_devifab.dataInfo(dataptrDevi, fabsizeDevi, iwhichfab, idir);

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
  void copyToDevice(      devidistrib_t & a_deviDat,
                    const hostdistrib_t & a_hostDat)
  {
    DataIterator dit = a_hostDat.dataIterator();
    for(int ibox = 0;  ibox < dit.size(); ibox++)
    {
      const auto & hostfab = a_hostDat[dit[ibox]];
      auto       & devifab = a_deviDat[dit[ibox]];
      copyToDevice(hostfab, devifab);
    }
  }

  ///
  static 
  void copyToDevice(const EBHostFluxData<Real, NUMCOMPS> & a_hostfab,
                          EBFluxData<    Real, NUMCOMPS> & a_devifab)
  {
    for(unsigned int idir = 0; idir < DIM; idir++)
    {
      for(unsigned int iwhichfab = 0; iwhichfab < 2; iwhichfab++)
      {
        size_t fabsizeHost, fabsizeDevi; 
        Real*  dataptrHost, *dataptrDevi;
        a_hostfab.dataInfo(dataptrHost, fabsizeHost, iwhichfab, idir);
        a_devifab.dataInfo(dataptrDevi, fabsizeDevi, iwhichfab, idir);

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
    shared_ptr<graph_distrib_t> m_graphs;
};


#include "Chombo_NamespaceFooter.H"

#endif
