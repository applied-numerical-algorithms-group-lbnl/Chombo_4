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

#include "Proto.H"
#include "EBProto.H"
#include "Chombo_LevelData.H"
#include "Chombo_AMRIO.H"
#include "Chombo_ProtoInterface.H"
#include "Chombo_EBChombo.H"
#include "Chombo_ProtoFactories.H"
#include "Chombo_LevelBoxData.H"
#include "Chombo_SPMD.H"
#include "Proto_EBBoxData.H"
#include "Proto_EBHostData.H"
#include "Chombo_NamespaceHeader.H"

#ifdef CH_MPI
using CH4_SPMD::Chombo_MPI;
#endif
template<unsigned int srcComp, unsigned int dstComp>
class ScalarAliasFactory: public DataFactory< EBBoxData<CELL, Real, dstComp> >
{
public:

  ScalarAliasFactory(LevelData<EBBoxData<CELL, Real, DIM> >* a_source, 
                     unsigned int                            a_whichcomp)
  {
    m_source    = a_source;
    m_whichcomp = a_whichcomp;
  }
  virtual ~ScalarAliasFactory()
  {
  }

  /// factory function.  creates a new 'T' object
  /** creates a new 'T' object and returns a pointer to it.
   */
  virtual EBBoxData<CELL, Real, dstComp> * create(const Box& a_box, int ncomps, const DataIndex& a_datInd) const
  {
    EBBoxData<CELL, Real, srcComp>& vecholder = (*m_source)[a_datInd];
    EBBoxData<CELL, Real, dstComp>* retval = new EBBoxData<CELL, Real, dstComp>();
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
  ScalarAliasFactory();
  LevelData<EBBoxData<CELL, Real, srcComp> >* m_source;
  unsigned int                                m_whichcomp;

};

///
/**
*/
template<CENTERING cent, int  NUMCOMPS>
class EBLevelBoxData: public LevelData<EBBoxData<cent, Real, NUMCOMPS> >
{
public:

  using LevelData<EBBoxData<cent, Real, NUMCOMPS> >::define;
  
#ifdef PROTO_CUDA
  //just shutting up the compiler
//  using LevelData<EBBoxData<cent, Real, NUMCOMPS> >::define;
#endif
  typedef GraphConstructorFactory< EBBoxData<cent, Real, NUMCOMPS> > devifactory_t;
  typedef GraphConstructorFactory<EBHostData<cent, Real, NUMCOMPS> > hostfactory_t;


  ///
  EBLevelBoxData(): LevelData<EBBoxData<cent, Real, NUMCOMPS> >()
  {
  }


  ///
  EBLevelBoxData(const DisjointBoxLayout        & a_dp, 
                 const IntVect                  & a_ghost,
                 shared_ptr<LevelData<EBGraph>  > a_graphs)
    : LevelData<EBBoxData<cent, Real, NUMCOMPS> >(a_dp, NUMCOMPS, a_ghost, devifactory_t(a_graphs)), m_graphs(a_graphs)
  {
    checkIfPeriodic(a_dp);
  }

  ///
  virtual ~EBLevelBoxData()
  {
  }

  ///
  EBLevelBoxData&
  operator+=(const EBLevelBoxData& a_input)
  {
    DisjointBoxLayout dbl = this->disjointBoxLayout();
    DataIterator dit = dbl.dataIterator();
    for(int ibox = 0;  ibox < dit.size(); ibox++)
    {
      (*this)[dit[ibox]] += a_input[dit[ibox]];
    }
    return *this;
  }
  
  ///
  void checkIfPeriodic(const DisjointBoxLayout & a_dp) const
  {
    if(a_dp.physDomain().isPeriodic())
    {
      MayDay::Error("EB applications do not support periodic boundary conditions yet.");
    }
  }

  ///
  virtual
  void define(const DisjointBoxLayout& a_dp, 
              const IntVect          & a_ghost,
              shared_ptr<LevelData<EBGraph>  > a_graphs)
  {
    checkIfPeriodic(a_dp);
    m_graphs = a_graphs;
    LevelData<EBBoxData<cent, Real, NUMCOMPS> >::define(a_dp, NUMCOMPS, a_ghost, devifactory_t(a_graphs));
  }


  template <unsigned int srcComp>
  void define(EBLevelBoxData<cent, srcComp>           & a_input,
              unsigned int                            a_comp,
              shared_ptr<LevelData<EBGraph>  >        a_graphs)
  {
    //for making scalar holders out of vector ones.

    auto dp    = a_input.disjointBoxLayout();
    auto ghost = a_input.ghostVect();
    m_graphs = a_graphs;
    LevelData<EBBoxData<cent, Real, srcComp> >* inputPtr = static_cast<LevelData<EBBoxData<cent, Real, srcComp> > * >(&a_input);
    ScalarAliasFactory<srcComp, NUMCOMPS>  fact(inputPtr, a_comp) ;

    LevelData<EBBoxData<cent, Real, NUMCOMPS> >::
      define(dp, NUMCOMPS, ghost, fact);
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
  void copyTo(const Interval                      &            a_srcComps,
              EBLevelBoxData<cent, NUMCOMPS>      &            a_dst,
              const Interval                      &            a_dstComps,
              const Copier                        &            a_copier) const
  {
    CH_TIME("EBLevelBoxData::copyTo");
    LevelData<EBHostData<cent, Real, NUMCOMPS> > hostdataSrc(this->disjointBoxLayout(), NUMCOMPS,  this->m_ghost   , hostfactory_t(m_graphs));
    LevelData<EBHostData<cent, Real, NUMCOMPS> > hostdataDst(a_dst.disjointBoxLayout(), NUMCOMPS, a_dst.ghostVect(), hostfactory_t(a_dst.m_graphs));

    copyToHost(hostdataSrc, *this);
    copyToHost(hostdataDst,  a_dst);
    
    hostdataSrc.copyTo(a_srcComps, hostdataDst, a_dstComps,a_copier);
    
    copyToDevice(a_dst, hostdataDst);
  }


  /// 
  virtual void exchange(const Copier& copier)
  {
    CH_TIME("EBLevelBoxData::exchange");
    this->LevelData<EBBoxData<cent, Real, NUMCOMPS>>::exchange(copier);
    
/*
    LevelData<EBHostData<cent, Real, NUMCOMPS> > hostdata(this->disjointBoxLayout(), NUMCOMPS,  this->m_ghost   , hostfactory_t(m_graphs));

    copyToHost(hostdata, *this);

    hostdata.exchange(copier);

    copyToDevice(*this, hostdata);
*/

  }

  ///
  Real maxNorm(int a_icomp) const
  {
    CH_TIME("EBLevelBoxData::maxNorm");
    LevelData<EBHostData<cent, Real, NUMCOMPS> > hostdata(this->disjointBoxLayout(), NUMCOMPS,  this->m_ghost   , hostfactory_t(m_graphs));
    copyToHost(     hostdata, *this);
    DataIterator dit = this->dataIterator();
    DisjointBoxLayout dbl = this->disjointBoxLayout();
    Real maxval = 0; 
    for(int ibox = 0;  ibox < dit.size(); ibox++)
    {
      const Box & grid     =         dbl[dit[ibox]];
      Bx grbx = ProtoCh::getProtoBox(grid);
      Real gridval = Proto::maxNormCell(hostdata[dit[ibox]], grbx, a_icomp);
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

  ///for dot products
  Real sum(int a_icomp) const
  {
    CH_TIME("EBLevelBoxData::maxNorm");
    LevelData<EBHostData<cent, Real, NUMCOMPS> > hostdata(this->disjointBoxLayout(), NUMCOMPS,  this->m_ghost   , hostfactory_t(m_graphs));
    copyToHost(     hostdata, *this);
    DataIterator dit = this->dataIterator();
    DisjointBoxLayout dbl = this->disjointBoxLayout();
    Real sumval = 0; 
    for(int ibox = 0;  ibox < dit.size(); ibox++)
    {
      const Box & grid     =         dbl[dit[ibox]];
      const EBGraph& graph = (*m_graphs)[dit[ibox]];
      for(BoxIterator bit(grid); bit.ok(); ++bit)
      {
        Point pt = ProtoCh::getPoint(bit());
        vector<EBIndex<CELL> > vofspt = graph.getVoFs(pt);
        for(int ivec = 0 ; ivec < vofspt.size(); ivec++)
        {
          const EBIndex<CELL>& vof = vofspt[ivec];
          sumval += hostdata[dit[ibox]](vof, a_icomp);
        }
      }
    }

#ifdef CH_MPI    
    Real sendBuf = sumval;
    int result = MPI_Allreduce(&sendBuf, &sumval, 1, MPI_CH_REAL,MPI_SUM, Chombo_MPI::comm);

    if (result != MPI_SUCCESS)
    {
      MayDay::Error("Communication error in EBLevelBoxData::maxNorm");
    }
#endif
    return sumval;
  }

  
  ///
  void 
  writeToFileHDF5(const string& a_filename, Real a_coveredval) const
  {
#ifdef CH_USE_HDF5
    LevelData<FArrayBox>                         hostfabs(this->disjointBoxLayout(), NUMCOMPS, this->m_ghost);
    LevelData<EBHostData<cent, Real, NUMCOMPS> > hostdata(this->disjointBoxLayout(), NUMCOMPS,  this->m_ghost   , hostfactory_t(m_graphs));
    copyToHost(     hostdata, *this);
    copyHostToFAB(  hostfabs, hostdata, a_coveredval);
    writeLevelname(&hostfabs, a_filename.c_str());
#else
    printf("hdf5 not enabled\n");
#endif
  }

  ///
  void 
  writeToCheckPoint(HDF5Handle& a_handle, string a_varname) const
  {
    LevelData<EBHostData<cent, Real, NUMCOMPS> >
      hostdata(this->disjointBoxLayout(), NUMCOMPS,  this->m_ghost   , hostfactory_t(m_graphs));
    copyToHost(     hostdata, *this);
    write(a_handle, hostdata, a_varname.c_str());
  }

  ///
  void 
  readFromCheckPoint(HDF5Handle& a_handle, string a_varname) 
  {
    auto grids = this->disjointBoxLayout();
    Interval interv(0, NUMCOMPS-1);
    LevelData<EBHostData<cent, Real, NUMCOMPS> >
      hostdata(grids, NUMCOMPS, this->m_ghost, hostfactory_t(m_graphs));
    read(a_handle, hostdata, a_varname, grids, interv, false);
    copyToDevice(*this,     hostdata);
  }
  
  void dumpFile(int fileNumber, Real coveredValue) const
  {
    string filename = "tmp"+std::to_string(fileNumber)+".hdf5";
    writeToFileHDF5(filename, coveredValue);
  }

  /// This is for plotfiles only  as multivalued stuff will be lost.
  static 
  void copyHostToFAB(LevelData<FArrayBox>                               & a_fabxDat,
                     const LevelData<EBHostData<cent, Real, NUMCOMPS> > & a_hostDat,
                     const Real& a_coveredval)
  {
    LevelData< EBHostData<cent, Real, NUMCOMPS> >& hostDat = (LevelData< EBHostData<cent, Real, NUMCOMPS> >&) (a_hostDat);
    DisjointBoxLayout dbl = a_fabxDat.disjointBoxLayout();
    DataIterator dit = a_hostDat.dataIterator();
    for(int ibox = 0;  ibox < dit.size(); ibox++)
    {
      EBHostData<cent, Real, NUMCOMPS> & hostdat =    hostDat[dit[ibox]];
      FArrayBox                        & fabxdat =  a_fabxDat[dit[ibox]];
      fabxdat.setVal(a_coveredval);
      Box chbox = dbl[dit[ibox]];
      Proto::Box prbox = ProtoCh::getProtoBox(chbox);
      PointSet allivs(prbox);
      const EBGraph& ebgraph = hostdat.ebgraph();
      
      for(PointSetIterator ivsit(allivs); ivsit.ok(); ++ivsit)
      {
        Point pt = ivsit();
        IntVect iv = ProtoCh::getIntVect(pt);
        if(!ebgraph.isCovered(pt)) //covered data is set with the setVal above
        {
          //might just work for cell centered data
          EBIndex<cent> ebind;
          ebind.m_pt = pt;
          ebind.m_vofIDMe = 0;
          ebind.m_vofIDLo = 0;
          ebind.m_isBoundary = false;
          for(int icomp = 0; icomp < NUMCOMPS; icomp++)
          {
            fabxdat(iv, icomp) = hostdat(ebind, icomp);
          }
        }
      }
    }
  }


  /// this is for mpi stuff
  static 
  void copyToHost(LevelData<      EBHostData<cent, Real, NUMCOMPS> > & a_hostDat,
                  const LevelData< EBBoxData<cent, Real, NUMCOMPS> > & a_deviDat)
  {
    CH_TIME("EBLevelBoxData::copyToHost");
    LevelData< EBBoxData<cent, Real, NUMCOMPS> >& deviDat = (LevelData< EBBoxData<cent, Real, NUMCOMPS> >&) (a_deviDat);
    DataIterator dit = a_hostDat.dataIterator();
    for(int ibox = 0;  ibox < dit.size(); ibox++)
    {
      EBHostData<cent, Real, NUMCOMPS> & hostfab = a_hostDat[dit[ibox]];
      EBBoxData< cent, Real, NUMCOMPS> & devifab =   deviDat[dit[ibox]];
      for(unsigned int iwhichfab = 0; iwhichfab < 2; iwhichfab++)
      {
        size_t fabsizeHost, fabsizeDevi; 
        Real*  dataptrHost;
        Real*  dataptrDevi;
        hostfab.dataInfo(dataptrHost, fabsizeHost, iwhichfab);
        devifab.dataInfo(dataptrDevi, fabsizeDevi, iwhichfab);

        size_t datSize = sizeof(Real)*fabsizeDevi;

#ifdef PROTO_CUDA
        protoMemcpy(dataptrHost, dataptrDevi, datSize, protoMemcpyDeviceToHost);
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
  void copyToHost(EBHostData<cent, Real, NUMCOMPS> & hostfab,
                  EBBoxData< cent, Real, NUMCOMPS> & devifab)
  {

    for(unsigned int iwhichfab = 0; iwhichfab < 2; iwhichfab++)
    {
      size_t fabsizeHost, fabsizeDevi; 
      Real*  dataptrHost;
      Real*  dataptrDevi;
      hostfab.dataInfo(dataptrHost, fabsizeHost, iwhichfab);
      devifab.dataInfo(dataptrDevi, fabsizeDevi, iwhichfab);

      size_t datSize = sizeof(Real)*fabsizeDevi;

#ifdef PROTO_CUDA
      protoMemcpy(dataptrHost, dataptrDevi, datSize, protoMemcpyDeviceToHost);
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
  

  
  ///
  static 
  void copyToDevice(LevelData<        EBBoxData<cent, Real, NUMCOMPS> > & a_deviDat,
                    const LevelData< EBHostData<cent, Real, NUMCOMPS> > & a_hostDat)
  {
    CH_TIME("EBLevelBoxData::copyToDevice");
    LevelData< EBBoxData<cent, Real, NUMCOMPS> >& deviDat = const_cast<LevelData< EBBoxData<cent, Real, NUMCOMPS> >  &> (a_deviDat);
    LevelData<EBHostData<cent, Real, NUMCOMPS> >& hostDat = const_cast<LevelData< EBHostData<cent, Real, NUMCOMPS> > &> (a_hostDat);
    DataIterator dit = a_hostDat.dataIterator();
    for(int ibox = 0;  ibox < dit.size(); ibox++)
    {
      EBHostData<cent, Real, NUMCOMPS> & hostfab =   hostDat[dit[ibox]];
      EBBoxData< cent, Real, NUMCOMPS> & devifab =   deviDat[dit[ibox]];
      copyToDevice(devifab, hostfab);
    }
  }

  ///
  static 
  void copyToDevice( EBBoxData< cent, Real, NUMCOMPS> & devifab,
                     EBHostData<cent, Real, NUMCOMPS> & hostfab)
                     
  {
    for(unsigned int iwhichfab = 0; iwhichfab < 2; iwhichfab++)
    {
      size_t fabsizeHost, fabsizeDevi; 
      Real*  dataptrHost, *dataptrDevi;
      hostfab.dataInfo(dataptrHost, fabsizeHost, iwhichfab);
      devifab.dataInfo(dataptrDevi, fabsizeDevi, iwhichfab);
      if(fabsizeDevi != fabsizeHost)
      {
        MayDay::Error("incompatible host and device sizes in EBLevelBoxData");
      }
      size_t datSize = sizeof(Real)*fabsizeDevi;

#ifdef PROTO_CUDA
      protoMemcpy(dataptrDevi, dataptrHost, datSize, protoMemcpyHostToDevice);
      protoError err = protoGetLastError();
      if (err != protoSuccess)
      {
        fprintf(stderr, "protoGetLastError() failed at %s:%i : %s\n",
                __FILE__, __LINE__, protoGetErrorString(err));
	std::abort();
      }
#else
      memcpy(dataptrDevi, dataptrHost, datSize);
#endif  
    }
  }

  static int 
  checkAnswer(const EBBoxData<cent, Real, NUMCOMPS> & a_calculatedAns,
              const EBBoxData<cent, Real, NUMCOMPS> & a_referenceAns,
              const Bx & a_grid,
              unsigned int a_icomp,
              Real a_tol)
  {
    Bx      inbox = a_calculatedAns.inputBox();
    EBGraph graph = a_calculatedAns.ebgraph();
    CH_assert(inbox == a_referenceAns.inputBox());
  
  
    EBBoxData<cent, Real, NUMCOMPS> difference(inbox, graph);
    difference.setVal(0.);
    difference += a_calculatedAns;
    difference -= a_referenceAns;

    EBHostData<cent, Real, NUMCOMPS> hostdat(inbox, graph);
    copyToHost(hostdat, difference);
    Real maxval = Proto::maxNormCell(hostdat, a_grid, a_icomp);
    int retval = 0;
    if(maxval > a_tol)
    {
      retval = 1;
    }
    return retval;
  }


private:
    shared_ptr<LevelData<EBGraph>  > m_graphs;
};

///write EB file so visit can recognize it as EB
template <int NUMCOMPS>
inline void
writeEBLevelHDF5(const string& a_filename,
                 const EBLevelBoxData<CELL, NUMCOMPS>  & a_data,
                 const EBLevelBoxData<CELL, 1       >  & a_kapp,
                 const Box                             & a_domain,
                 shared_ptr<LevelData<EBGraph>  > a_graphs,
                 Real                             a_coveredval = 0.0,
                 Real a_dx= 1.0, Real a_dt=1.0, Real a_time=1.0)
{

#ifdef CH_USE_HDF5
  typedef GraphConstructorFactory<EBHostData<CELL, Real, NUMCOMPS> > hostfactDat_t;
  typedef GraphConstructorFactory<EBHostData<CELL, Real, 1       > > hostfactKap_t;
  
  DisjointBoxLayout grids = a_data.disjointBoxLayout();
  IntVect dataGh = a_data.ghostVect();
  IntVect kappGh = a_kapp.ghostVect();
  LevelData<EBHostData<CELL, Real, NUMCOMPS> > hostdata(grids, NUMCOMPS, dataGh,  hostfactDat_t(a_graphs));
  LevelData<EBHostData<CELL, Real, 1       > > hostkapp(grids,        1, kappGh,  hostfactKap_t(a_graphs));

  EBLevelBoxData<CELL, NUMCOMPS>::copyToHost(hostdata, a_data);
  EBLevelBoxData<CELL, 1       >::copyToHost(hostkapp, a_kapp);

  LevelData<FArrayBox>    hostfabs(grids, NUMCOMPS+1, IntVect::Zero);
  
  DataIterator dit = grids.dataIterator();
  for(int ibox = 0;  ibox < dit.size(); ibox++)
  {
    EBHostData<CELL, Real, NUMCOMPS> & data =    hostdata[dit[ibox]];
    EBHostData<CELL, Real, 1       > & kapp =    hostkapp[dit[ibox]];
    FArrayBox                        & fabu =    hostfabs[dit[ibox]];

    fabu.setVal(a_coveredval);
    Box chbox = grids[dit[ibox]];
    Proto::Box prbox = ProtoCh::getProtoBox(chbox);
    PointSet allivs(prbox);
    const EBGraph& ebgraph = (*a_graphs)[dit[ibox]];
      
    for(PointSetIterator ivsit(allivs); ivsit.ok(); ++ivsit)
    {
      Point pt = ivsit();
      IntVect iv = ProtoCh::getIntVect(pt);
      if(ebgraph.isCovered(pt))
      {
        fabu(iv, NUMCOMPS) = 0; //kappa = 0
      }
      else
      {
        //might just work for cell centered data
        EBIndex<CELL> ebind;
        ebind.m_pt = pt;
        ebind.m_vofIDMe = 0;
        ebind.m_vofIDLo = 0;
        ebind.m_isBoundary = false;
        for(int icomp = 0; icomp < NUMCOMPS; icomp++)
        {
          Real val = data(ebind, icomp);
          fabu(iv, icomp)  = val;
        }
        fabu(iv, NUMCOMPS) = kapp(ebind, 0);
      }
    }
  }
  Vector<string> names(NUMCOMPS+1);
  for(unsigned int icomp = 0; icomp < NUMCOMPS; icomp++)
  {
    names[icomp] = string("var_") + std::to_string(icomp);
  }
  names[NUMCOMPS] = string("fraction-0");

  Vector<DisjointBoxLayout>     vecgrids(1, grids);
  Vector<LevelData<FArrayBox> *> vecdata(1, &hostfabs);
  int numlev = 1;
  
  WriteAMRHierarchyHDF5(a_filename, vecgrids, vecdata, names, a_domain,
                        a_dx, a_dt, a_time, Vector<int>(1,2), numlev);
#else
    printf("hdf5 not enabled\n");
#endif

}

#include "Chombo_NamespaceFooter.H"

#endif
