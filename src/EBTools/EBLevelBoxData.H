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
#include "EBChombo.H"
#include "ProtoFactories.H"
#include "LevelBoxData.H"
#include "NamespaceHeader.H"


///
/**
*/
template<CENTERING cent, int  NUMCOMPS>
class EBLevelBoxData: public LevelData<EBBoxData<cent, Real, NUMCOMPS> >
{
public:

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
  void checkIfPeriodic(const DisjointBoxLayout & a_dp) const
  {
    if(a_dp.physDomain().isPeriodic())
    {
      MayDay::Error("EB applications do not support periodic boundary conditions yet.");
    }
  }

  ///
  void define(const DisjointBoxLayout& a_dp, 
              const IntVect          & a_ghost,
              shared_ptr<LevelData<EBGraph>  > a_graphs)
  {
    checkIfPeriodic(a_dp);
    m_graphs = a_graphs;
    LevelData<EBBoxData<cent, Real, NUMCOMPS> >::define(a_dp, NUMCOMPS, a_ghost, devifactory_t(a_graphs));
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
              LevelData<BoxData<Real, NUMCOMPS> > &            a_dst,
              const Interval                      &            a_dstComps,
              const Copier                        &            a_copier) const
  {
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

    LevelData<EBHostData<cent, Real, NUMCOMPS> > hostdata(this->disjointBoxLayout(), NUMCOMPS,  this->m_ghost   , hostfactory_t(m_graphs));

    copyToHost(hostdata, *this);

    hostdata.exchange();

    copyToDevice(*this, hostdata);
  }

  ///
  Real maxNorm(int a_icomp) const
  {
    LevelData<EBHostData<cent, Real, NUMCOMPS> > hostdata(this->disjointBoxLayout(), NUMCOMPS,  this->m_ghost   , hostfactory_t(m_graphs));
    copyToHost(     hostdata, *this);
    DataIterator dit = this->dataIterator();
    DisjointBoxLayout dbl = this->disjointBoxLayout();
    Real maxval = 0; 
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
          Real absval = std::abs(hostdata[dit[ibox]](vof, a_icomp));
          if(absval > maxval)
          {
            maxval = absval;
          }
        }
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
  void 
  writeToFileHDF5(const string& a_filename, const Real& a_coveredval) const
  {
    LevelData<FArrayBox>                         hostfabs(this->disjointBoxLayout(), NUMCOMPS, this->m_ghost);
    LevelData<EBHostData<cent, Real, NUMCOMPS> > hostdata(this->disjointBoxLayout(), NUMCOMPS,  this->m_ghost   , hostfactory_t(m_graphs));
    copyToHost(     hostdata, *this);
    copyHostToFAB(  hostfabs, hostdata, a_coveredval);
    writeLevelname(&hostfabs, a_filename.c_str());
  }

  /// this is for hdf5 stuff as multivalued stuff will be lost
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
        if(!ebgraph.isCovered(pt))
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
  void copyToDevice(LevelData<        EBBoxData<cent, Real, NUMCOMPS> > & a_deviDat,
                    const LevelData< EBHostData<cent, Real, NUMCOMPS> > & a_hostDat)
  {
    LevelData< EBBoxData<cent, Real, NUMCOMPS> >& deviDat = const_cast<LevelData< EBBoxData<cent, Real, NUMCOMPS> >  &> (a_deviDat);
    LevelData<EBHostData<cent, Real, NUMCOMPS> >& hostDat = const_cast<LevelData< EBHostData<cent, Real, NUMCOMPS> > &> (a_hostDat);
    DataIterator dit = a_hostDat.dataIterator();
    for(int ibox = 0;  ibox < dit.size(); ibox++)
    {
      EBHostData<cent, Real, NUMCOMPS> & hostfab =   hostDat[dit[ibox]];
      EBBoxData< cent, Real, NUMCOMPS> & devifab =   deviDat[dit[ibox]];
      for(unsigned int iwhichfab = 0; iwhichfab < 2; iwhichfab++)
      {
        size_t fabsizeHost, fabsizeDevi; 
        Real*  dataptrHost, *dataptrDevi;
        hostfab.dataInfo(dataptrHost, fabsizeHost, iwhichfab);
        devifab.dataInfo(dataptrDevi, fabsizeDevi, iwhichfab);

        size_t datSize = sizeof(Real)*fabsizeDevi;

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
