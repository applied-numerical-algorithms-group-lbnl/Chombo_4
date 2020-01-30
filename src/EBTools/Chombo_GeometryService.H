#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif


#ifndef _Chombo_GeometryService_H_
#define _Chombo_GeometryService_H_


#include "Proto.H"
#include "EBProto.H"
#include "Chombo_LevelData.H"
#include "Chombo_AMRIO.H"
#include "Chombo_ProtoInterface.H"
#include "Chombo_EBChombo.H"
#include "Chombo_BRMeshRefine.H"
#include "Chombo_LoadBalance.H"
#include "Chombo_NamespaceHeader.H"

///
/**
   makes geometric description from an implicit function.
*/
template <int order>
class LevelGeometryService
{

public:

  typedef IndexedMoments<DIM  , order> IndMomDIM;
  typedef IndexedMoments<DIM-1, order> IndMomSDMinOne;
  typedef MomentIterator<DIM  , order> MomItDIM;
  typedef MomentIterator<DIM-1, order> MomItSDMinOne;
  typedef IrregNode<order> inode;
  typedef HostIrregData<CELL,      IndMomDIM , 1>  VoluData;
  typedef HostIrregData<BOUNDARY,  IndMomDIM , 1>  EBFaData;
  typedef HostIrregData<XFACE, IndMomSDMinOne, 1>  XFacData;
  typedef HostIrregData<YFACE, IndMomSDMinOne, 1>  YFacData;
  typedef HostIrregData<ZFACE, IndMomSDMinOne, 1>  ZFacData;

  ///
  /**
     basif      :  implicit function
     origin     :  location of the lower left corner of the domain
     dx         :  grid spacing
     domain     :  finest domain in the calculation
     grids      :  layout of the  domain
     ghost      :  maximum number of ghost cells needed in calculation
     maxCoarsen :  maximum level of coarsening 
  */
  LevelGeometryService(const shared_ptr<BaseIF>     &      a_baseIF,
                       const RealVect               &      a_origin,
                       const double                 &      a_dx,
                       const Box                    &      a_domain,
                       const DisjointBoxLayout      &      a_grids,
                       const int                    &      a_ghost);



  LevelGeometryService(const LevelGeometryService   &      a_finerLevel,
                       const DisjointBoxLayout      &      a_grids);

  ///
  /**
     Define the internals of the input ebisRegion.   Public for testing.
  */
  inline void fillGraph(HostBoxData<int>          &   a_regIrregCovered,
                        vector<inode>             &   a_nodes,
                        const Box                 &   a_validRegion,
                        const Box                 &   a_ghostRegion,
                        const Box                 &   a_domain) const;


    
  ///
  inline const shared_ptr<LevelData<EBGraph>  >
  getGraphs() const
  {
    return   m_graphs;
  }

  ///
  inline const shared_ptr<LevelData<VoluData> >
  getVoluData() const
  {
    return m_voludata;
  }

  ///
  inline const shared_ptr<LevelData<EBFaData> >
  getEBFaceData() const
  {
    return m_ebfadata;
  }

  ///
  inline const shared_ptr<LevelData<XFacData> >
  getXFaceData() const
  {
    return m_xfacdata;
  }

  ///
  inline const shared_ptr<LevelData<YFacData> >
  getYFaceData() const
  {
    return m_yfacdata;
  }

  ///
  inline const shared_ptr<LevelData<ZFacData> > 
  getZFaceData() const
  {
    return m_zfacdata;
  }


  /// fill data holder with values of the volume fraction.
  inline void fillKappa(EBHostData<CELL,Real, 1> & a_kappa,
                        const Box                & a_grbx ,
                        const DataIndex          & a_dit) const
  {
    const VoluData & volmo = (*m_voludata)[a_dit];
    const EBGraph  & graph = (*m_graphs)[a_dit];
    HostBoxData<        Real, 1>& reghost = a_kappa.getRegData();
    HostIrregData<CELL, Real, 1>& irrhost = a_kappa.getIrrData();
    for(BoxIterator bit(a_grbx); bit.ok(); ++bit)
    {
      Point pt  = ProtoCh::getPoint(bit());
      if(graph.isRegular(pt))
      {
        reghost(pt, 0) = 1.0;
      }
      else if(graph.isCovered(pt))
      {
        reghost(pt, 0) = 0.0;
      }
      else
      {
        vector<EBIndex<CELL> > vofs = graph.getVoFs(pt);
        for(int ivec = 0; ivec < vofs.size(); ivec++)
        {
          const EBIndex<CELL>& vof = vofs[ivec];
          const IndMomDIM&  momspt = volmo(vof, 0);
          double kappavof = momspt[0]/m_cellVol;

          reghost(pt  , 0) = kappavof;
          irrhost(vof , 0) = kappavof;
        }
      }
    }
  }
  
protected:
  //for coarsening volume moment data
  template <CENTERING cent>
  inline void shiftAndIncrement(IndMomDIM                  &      a_output, 
                                const IndMomDIM            &      a_input,
                                const EBIndex<cent>        &      a_coarVoF,
                                const EBIndex<cent>        &      a_fineVoF,
                                const LevelGeometryService &      a_finerLevel)
  {
    Real coarDx =              m_dx;
    Real fineDx = a_finerLevel.m_dx;
    RealVect coarLoc = indexLoc(a_coarVoF, RealVect::Zero(), coarDx);
    RealVect fineLoc = indexLoc(a_fineVoF, RealVect::Zero(), fineDx);

    IndexTM<Real, SpaceDim> shiftVec;
    for(int idir = 0; idir < SpaceDim; idir++)
    {
      shiftVec[idir] =  fineLoc[idir] - coarLoc[idir];
    }
    IndMomDIM increment = a_input;
    increment.shift(shiftVec);
  
    a_output += increment;
  }

  //for coarsening face moment data
  template <CENTERING cent>
  inline void shiftAndIncrement(IndMomSDMinOne             &      a_output, 
                                const IndMomSDMinOne       &      a_input,
                                const EBIndex<cent>        &      a_coarFace,
                                const EBIndex<cent>        &      a_fineFace,
                                const LevelGeometryService &      a_finerLevel)
  {
    int faceDir = a_coarFace.direction();
    Real coarDx =              m_dx;
    Real fineDx = a_finerLevel.m_dx;
    RealVect coarLoc = indexLoc(a_coarFace, RealVect::Zero(), coarDx);
    RealVect fineLoc = indexLoc(a_fineFace, RealVect::Zero(), fineDx);

    IndexTM<Real, SpaceDim-1> shiftVec;
    int iindex = 0;
    for(int idir = 0; idir < SpaceDim; idir++)
    {
      if(idir != faceDir)
      {
        shiftVec[iindex] =  fineLoc[idir] - coarLoc[idir];
        iindex++;
      }
    }
    IndMomSDMinOne increment = a_input;
    increment.shift(shiftVec);
  
    a_output += increment;
  }

  void 
  coarsenVoFData(const vector< EBIndex<CELL> >  & a_vofs,
                 const LevelGeometryService     & a_finerLevel,
                 const DataIndex                & a_dit);

  void 
  coarsenXFaceData(const vector< EBIndex<XFACE> > & a_faces,
                   const LevelGeometryService     & a_finerLevel,
                   const DataIndex                & a_dit);
  

  void 
  coarsenYFaceData(const vector< EBIndex<YFACE> > & a_faces,
                   const LevelGeometryService     & a_finerLevel,
                   const DataIndex                & a_dit);
                   

#if DIM==3
  void 
  coarsenZFaceData(const vector< EBIndex<ZFACE> > & a_faces,
                   const LevelGeometryService     & a_finerLevel,
                   const DataIndex                & a_dit);
#endif

  void
  defineFinestLevel(const RealVect               & a_origin,
                    const double                 & a_dx );

  void
  defineData();
  void
  coarsenGraph(LevelGeometryService&      a_finerLevel);
  void
  coarsenData( LevelGeometryService&      a_finerLevel);

  
  void
  define(const RealVect               & a_origin,
         const double                 & a_dx );


  bool onBoxBoundary(const Point        & a_iv, 
                     const Box            & a_box,
                     const int            & a_dir,
                     const Side::LoHiSide & a_sd) const;



  //converts a RealVect in physical coordinates to a RealVect in coordinates relative to a cell center
  inline RealVect convert2RelativeCoord(const RealVect& a_rVect)const;

  //converts a IndexTM<double,DIM> in physical coordinates to a
  //RealVect in coordinates relative to a cell center
  inline RealVect convert2RelativeCoord(const IndexTM<double,DIM>& a_rVect)const;

  inline void
  fillNewNode(inode                    &     a_node,
              const PointSet           &     a_ivsIrreg,
              const Point              &     a_iv) const;


  inline shared_ptr<BaseIF> getBaseIF() const
  {
    return m_baseIF;
  }

private:
  inline void
  fixRegularCellsNextToCovered(vector<inode>       & a_nodes, 
                               HostBoxData<int>    & a_regIrregCovered,
                               const Box           & a_validRegion,
                               const Point         & a_iv) const;
  inline void
  getFullNodeWithCoveredFace(inode                   & a_newNode, 
                             const HostBoxData<int>  & a_regIrregCovered,
                             const Point             & a_iv) const;

  ///funny constructor for when grids change and we have to copy data over
  /// to do coarsening
  LevelGeometryService(const Box                               &  a_domain,
                       const DisjointBoxLayout                 &  a_grids,          
                       const shared_ptr<LevelData<EBGraph> >   &  a_graphs,
                       const shared_ptr<LevelData<VoluData> >  &  a_voludata,
                       const shared_ptr<LevelData<EBFaData> >  &  a_ebfadata,
                       const shared_ptr<LevelData<XFacData> >  &  a_xfacdata,       
                       const shared_ptr<LevelData<YFacData> >  &  a_yfacdata,       
                       const shared_ptr<LevelData<ZFacData> >  &  a_zfacdata,       
                       const RealVect                          &  a_origin,         
                       const int                               &  a_order,          
                       const int                               &  a_degreeP,
                       const double                            &  a_dx,             
                       const IntVect                           &  a_geomgrow,       
                       const shared_ptr<BaseIF>                &  a_baseIF);         


  ///calls exchange on everything
  void exchange();

  ///
  Copier m_exchangeCopier;
  
  ///
  Box                                 m_domain;
    
  ///
  DisjointBoxLayout                   m_grids;
                                                
  ///                                         
  shared_ptr<LevelData<EBGraph> >     m_graphs;

  ///
  shared_ptr<LevelData<VoluData> >    m_voludata;

  ///
  shared_ptr<LevelData<EBFaData> >    m_ebfadata;

  ///
  shared_ptr<LevelData<XFacData> >    m_xfacdata;

  ///
  shared_ptr<LevelData<YFacData> >    m_yfacdata;

  ///
  shared_ptr<LevelData<ZFacData> >    m_zfacdata;
    
  //origin                            
  RealVect                            m_origin;
                                      
  int                                 m_order;
  int                                 m_degreeP;
  Real                                m_dx;
  Real                                m_cellVol;
                                      
  IntVect                             m_geomgrow;
                                      
  shared_ptr<BaseIF>                  m_baseIF;         
  

  inline void fillArc(vector<int>                          a_arc[DIM],
                      CutCellMoments<DIM>            &     a_cutCellMoments,
                      const int                      &     a_hilo,
                      const PointSet                 &     a_ivsIrreg,
                      const Point                    &     a_curriv) const;


private:
  //stuff disallowed for all the usual reasons.
  LevelGeometryService();
  LevelGeometryService(const LevelGeometryService& a_workshopin);
  void operator=(const LevelGeometryService& a_workshopin);


};




///
/**
   makes geometric description from an implicit function.
*/
template <int order>
class GeometryService
{

public:

  typedef IndexedMoments<DIM  , order> IndMomDIM;
  typedef IndexedMoments<DIM-1, order> IndMomSDMinOne;
  typedef MomentIterator<DIM  , order> MomItDIM;
  typedef MomentIterator<DIM-1, order> MomItSDMinOne;
  typedef IrregNode<order> inode;
  typedef HostIrregData<CELL,      IndMomDIM, 1>  VoluData;
  typedef HostIrregData<BOUNDARY,  IndMomDIM, 1>  EBFaData;
  typedef HostIrregData<XFACE, IndMomSDMinOne, 1>  XFacData;
  typedef HostIrregData<YFACE, IndMomSDMinOne, 1>  YFacData;
  typedef HostIrregData<ZFACE, IndMomSDMinOne, 1>  ZFacData;



  ///
  /**
     basif      :  implicit function
     origin     :  location of the lower left corner of the domain
     dx         :  grid spacing
     domain     :  finest domain in the calculation
     finestGrids:  layout of the finest domain
     ghost      :  maximum number of ghost cells needed in calculation
  */
  GeometryService(const shared_ptr<BaseIF>          &      a_baseIF,
                  const RealVect                    &      a_origin,
                  const double                      &      a_finestDx,
                  const Box                         &      a_finestDomain,
                  const Vector<DisjointBoxLayout>   &      a_grids,
                  const int                         &      a_ghost)
  {
    int numLevels = a_grids.size();
    m_levels.resize( numLevels);
    m_domains.resize(numLevels);
    m_levels[0]  = shared_ptr<LevelGeometryService<order> >(
      new LevelGeometryService<order>(a_baseIF, a_origin, a_finestDx, a_finestDomain, a_grids[0], a_ghost)
      );
    m_domains[0] = a_finestDomain;
    for(int ilev = 1; ilev < numLevels; ilev++)
    {
      m_levels[ilev]  = shared_ptr<LevelGeometryService<order> >(
        new LevelGeometryService<order>(*m_levels[ilev-1], a_grids[ilev])
        );
      m_domains[ilev] = coarsen(m_domains[ilev-1], 2);
    }
  }

  ///
  /**
     make vector of grids that cover the domain at various coarsenings
     This coarsens grids if it can, rebuilds them if it has to for coarser grids
     maxGridSize: maximum size of boxes on coarser grids
     maxCoarsen : maximum times the service coarsens (-1 means coarsen down to 2x2)
  */
  static void
  generateGrids(Vector<DisjointBoxLayout> & a_grids,
                const Box                 & a_domain,
                int a_maxGridSize =  32,
                int a_maxCoarsen  = -1)
  {
    int numLevels = 1;
    Box domlev = a_domain;
    int numCoarsen = 0;
    int minDomainSize = 2;
    int nref = 2;
    Vector<Box> domains(1, a_domain);
    bool domCheck = domlev.coarsenable(nref*minDomainSize);
    bool maxCheck = (a_maxCoarsen < 0) || (numCoarsen < a_maxCoarsen);
    while( domCheck && maxCheck)
    {
      numLevels++;
      numCoarsen++;
      domlev.coarsen(nref);
      domains.push_back(domlev);

      domCheck = domlev.coarsenable(nref*minDomainSize);
      maxCheck = (a_maxCoarsen < 0) || (numCoarsen < a_maxCoarsen);
    }
    a_grids.resize(numLevels);
    Vector<Box> boxes;
    Vector<int> procs;
    domainSplit(domains[0], boxes, a_maxGridSize);
    LoadBalance(procs, boxes);
    a_grids[0].define(boxes, procs, domains[0]);
    //a_grids[0].print();
    for(int ilev = 1; ilev < numLevels; ilev++)
    {
      if(a_grids[ilev-1].coarsenable(2))
      {
        coarsen(a_grids[ilev], a_grids[ilev-1], 2);
        //a_grids[ilev].print();
      }
      else
      {
        domainSplit(domains[ilev], boxes, a_maxGridSize);
        LoadBalance(procs, boxes);
        a_grids[ilev].define(boxes, procs, domains[ilev]);
      }
    }
  }
                
  ///
  /**
     single level define.
     basif      :  implicit function
     origin     :  location of the lower left corner of the domain
     dx         :  grid spacing
     domain     :  finest domain in the calculation
     finestGrids:  layout of the finest domain
     ghost      :  maximum number of ghost cells needed in calculation
  */
  GeometryService(const shared_ptr<BaseIF>     &      a_baseIF,
                  const RealVect               &      a_origin,
                  const double                 &      a_dx,
                  const Box                    &      a_finestDomain,
                  const DisjointBoxLayout      &      a_grids,
                  const int                    &      a_ghost)
  {
    int numLevels = 1;
    m_levels.resize( numLevels);
    m_domains.resize(numLevels);
    m_domains[0] = a_finestDomain;
    m_levels[0] = shared_ptr<LevelGeometryService<order> >
      (new LevelGeometryService<order>(a_baseIF, a_origin, a_dx, a_finestDomain, a_grids, a_ghost));
  }




  ///
  inline shared_ptr<LevelData<EBGraph>  >
  getGraphs(const Box& a_domain)
  {
    return   (*this)[a_domain]->getGraphs();
  }

  ///
  inline shared_ptr<LevelData<VoluData> >
  getVoluData(const Box& a_domain)
  {
    return   (*this)[a_domain]->getVoluData();
  }

  /// fill a data holder in with volume fractions
  inline void
  fillKappa(EBHostData<CELL,Real, 1>& a_kappa,
            const Box                & a_grid ,
            const DataIndex          & a_dit,
            const Box                & a_domain) const
  {
    
    (*this)[a_domain]->fillKappa(a_kappa, a_grid, a_dit);
  }
      
  ///
  inline shared_ptr<LevelData<EBFaData> >
  getEBFaceData(const Box& a_domain)
  {
    return   (*this)[a_domain]->getEBFaceData();
  }

  ///
  inline shared_ptr<LevelData<XFacData> >
  getXFaceData(const Box& a_domain)
  {
    return   (*this)[a_domain]->getXFaceData();
  }

  ///
  inline shared_ptr<LevelData<YFacData> >
  getYFaceData(const Box& a_domain)
  {
    return   (*this)[a_domain]->getYFaceData();
  }

  ///
  inline shared_ptr<LevelData<ZFacData> >
  getZFaceData(const Box& a_domain)
  {
    return   (*this)[a_domain]->getZFaceData();
  }

  ///
  inline shared_ptr<LevelGeometryService<order> > operator[](const Box& a_domain)
  {
    int ilev = getLevel(a_domain);
    return m_levels[ilev];
  }

    ///
  inline const shared_ptr<LevelGeometryService<order> > operator[](const Box& a_domain) const
  {
    int ilev = getLevel(a_domain);
    return m_levels[ilev];
  }

protected:

  vector<Box>                               m_domains;
  vector<shared_ptr<LevelGeometryService<order> > > m_levels;


private:


  inline const int getLevel(const Box & a_domain) const
  {
    int retval = -1;
    bool found = false;
    for(int ilev = 0; ilev < m_domains.size(); ilev++)
    {
      if(m_domains[ilev] == a_domain)
      {
        retval = ilev;
        found = true;
        break;
      }
    }

    PROTO_ASSERT(found, "did not find asked-for domain");

    return retval;
  }



  GeometryService();

  GeometryService(const GeometryService& a_workshopin);

  void operator=(const GeometryService& a_workshopin);


};
#include "Chombo_GeometryServiceImplem.H"
#include "Chombo_NamespaceFooter.H"
#endif
