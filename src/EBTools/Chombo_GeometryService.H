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
#include "Chombo_EBDataChoreography.H"
#include "Chombo_AMRIO.H"
#include "Chombo_ProtoInterface.H"
#include "Chombo_EBChombo.H"
#include "Chombo_BRMeshRefine.H"
#include "Chombo_LoadBalance.H"
#include <memory>
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
  //this ought to have a boundary type but 
  //in this case boundary and cell have the same semantic
  typedef HostIrregData<CELL,      IndMomDIM , 1>  EBNormalData;
  
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
                        const Box                 &   a_smghRegion,
                        const Box                 &   a_domain) const;

  typedef CH4_Data_Choreography::DistributedData<EBGraph> graph_distrib_t;
  typedef CH4_Data_Choreography::DistributedData<VoluData> volu_distrib_t;
  typedef CH4_Data_Choreography::DistributedData<EBFaData> ebfa_distrib_t;
  typedef CH4_Data_Choreography::DistributedData<XFacData> xfac_distrib_t;
  typedef CH4_Data_Choreography::DistributedData<YFacData> yfac_distrib_t;
  typedef CH4_Data_Choreography::DistributedData<ZFacData> zfac_distrib_t;
  typedef CH4_Data_Choreography::DistributedData<EBNormalData> ebnorm_distrib_t;
    
  ///
  inline DisjointBoxLayout const
  getDBL()
  {
    return m_grids;
  }
  
  ///
  inline const shared_ptr<graph_distrib_t>
  getGraphs() const
  {
    return   m_graphs;
  }

  ///
  inline const shared_ptr<volu_distrib_t>
  getVoluData() const
  {
    return m_voludata;
  }

  
  ///
  inline const shared_ptr<ebfa_distrib_t >
  getEBFaceData() const
  {
    return m_ebfadata;
  }

  ///
  inline const shared_ptr<xfac_distrib_t >
  getXFaceData() const
  {
    return m_xfacdata;
  }

  ///
  inline const shared_ptr<yfac_distrib_t >
  getYFaceData() const
  {
    return m_yfacdata;
  }

  ///
  inline const shared_ptr<zfac_distrib_t > 
  getZFaceData() const
  {
    return m_zfacdata;
  }

  ///Get eb moment in the x direction 
  inline const shared_ptr<ebnorm_distrib_t >
  getEBNormalData_x() const
  {
    return m_ebnormaldata[0];
  }

  ///Get eb moment in the y direction
  inline const shared_ptr<ebnorm_distrib_t>
  getEBNormalData_y() const
  {
    return m_ebnormaldata[1];
  }

  ///Get eb moment in the z direction
  inline const shared_ptr<ebnorm_distrib_t >
  getEBNormalData_z() const
  {
    return m_ebnormaldata[2];
  }

  /// fill data holder with values of the volume fraction.
  inline void fillKappa(EBHostData<CELL,Real, 1> & a_kappa,
                        const Box                & a_grbx ,
                        const DataIndex          & a_dit) const
  {
    const VoluData & volmo = (*m_voludata)[a_dit];
    const EBGraph  & graph = (*m_graphs)[a_dit];
    auto& reghost = a_kappa.getRegData();
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

  /// fill data holder with values of the EBFaceData.
  inline void fillEBFaceData(EBHostData<CELL,Real, 1> & a_ebfacedata,
                             const Box                & a_grbx ,
                             const DataIndex          & a_dit) const
  { 
    const EBFaData & ebfmo = (*m_ebfadata)[a_dit];
    const EBGraph  & graph = (*m_graphs)[a_dit];
    HostBoxData<        Real, 1>& reghost = a_ebfacedata.getRegData();
    HostIrregData<CELL, Real, 1>& irrhost = a_ebfacedata.getIrrData();
    for(BoxIterator bit(a_grbx); bit.ok(); ++bit)
    { 
      Point pt  = ProtoCh::getPoint(bit());
      if(graph.isRegular(pt))
      {
        reghost(pt, 0) = 0.0;
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
          const EBIndex<BOUNDARY> ebf = vof.getCutFace();
          const IndMomDIM&  momspt = ebfmo(ebf, 0);
          double ebfacevof = momspt[0]/m_dx;

          reghost(pt  , 0) = ebfacevof;
          irrhost(vof , 0) = ebfacevof;
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
  coarsenGraph(const LevelGeometryService&      a_finerLevel);
  void
  coarsenData( const LevelGeometryService&      a_finerLevel);

  
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


  ///calls exchange on everything
  void exchange();

  ///
  Box                            m_domain;
    
  ///
  DisjointBoxLayout              m_grids;
                                                
  ///                                         
  shared_ptr<graph_distrib_t >   m_graphs;

  ///
  shared_ptr<volu_distrib_t >    m_voludata;

  ///
  shared_ptr<ebfa_distrib_t >    m_ebfadata;

  ///
  shared_ptr<xfac_distrib_t >    m_xfacdata;

  ///
  shared_ptr<yfac_distrib_t >    m_yfacdata;

  ///
  shared_ptr<zfac_distrib_t >    m_zfacdata;

  ///
  shared_ptr<ebnorm_distrib_t >  m_ebnormaldata[3];
    
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
  typedef HostIrregData<CELL,      IndMomDIM, 1>      VoluData;
  typedef HostIrregData<BOUNDARY,  IndMomDIM, 1>      EBFaData;
  typedef HostIrregData<XFACE, IndMomSDMinOne, 1>     XFacData;
  typedef HostIrregData<YFACE, IndMomSDMinOne, 1>     YFacData;
  typedef HostIrregData<ZFACE, IndMomSDMinOne, 1>     ZFacData;
  typedef HostIrregData<CELL,       IndMomDIM, 1> EBNormalData;
  
  typedef CH4_Data_Choreography::DistributedData<EBGraph>       graph_distrib_t;
  typedef CH4_Data_Choreography::DistributedData<VoluData>       volu_distrib_t;
  typedef CH4_Data_Choreography::DistributedData<EBFaData>       ebfa_distrib_t;
  typedef CH4_Data_Choreography::DistributedData<XFacData>       xfac_distrib_t;
  typedef CH4_Data_Choreography::DistributedData<YFacData>       yfac_distrib_t;
  typedef CH4_Data_Choreography::DistributedData<ZFacData>       zfac_distrib_t;
  typedef CH4_Data_Choreography::DistributedData<EBNormalData> ebnorm_distrib_t;


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
    Vector<Box> boxes;
    Vector<int> procs;
    domainSplit(a_domain, boxes, a_maxGridSize);
    LoadBalance(procs, boxes);
    DisjointBoxLayout currDBL(boxes, procs, a_domain);
    a_grids.push_back(currDBL);
    //a_grids[0].print();
    while(currDBL.coarsenable(4))
    {
      DisjointBoxLayout coarDBL;
      coarsen(coarDBL, currDBL, 2);
      a_grids.push_back(coarDBL);
      currDBL = coarDBL;
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
  inline shared_ptr<graph_distrib_t>
  getGraphs(const Box& a_domain)
  {
    return   (*this)[a_domain]->getGraphs();
  }

  ///
  inline shared_ptr<volu_distrib_t >
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

  /// fill a data holder in with eb face data
  inline void
  fillEBFaceData(EBHostData<CELL,Real, 1>& a_ebfdata,
                 const Box                & a_grid ,
                 const DataIndex          & a_dit,
                 const Box                & a_domain) const
  { 

    (*this)[a_domain]->fillEBFaceData(a_ebfdata, a_grid, a_dit);
  }
      
  ///
  inline shared_ptr<ebfa_distrib_t >
  getEBFaceData(const Box& a_domain)
  {
    return   (*this)[a_domain]->getEBFaceData();
  }

  ///
  inline shared_ptr<xfac_distrib_t >
  getXFaceData(const Box& a_domain)
  {
    return   (*this)[a_domain]->getXFaceData();
  }

  ///
  inline shared_ptr<yfac_distrib_t >
  getYFaceData(const Box& a_domain)
  {
    return   (*this)[a_domain]->getYFaceData();
  }

  ///
  inline shared_ptr<zfac_distrib_t>
  getZFaceData(const Box& a_domain)
  {
    return   (*this)[a_domain]->getZFaceData();
  }

  ///Get eb moment in the x direction
  inline shared_ptr<ebnorm_distrib_t >
  getEBNormalData_x(const Box& a_domain)
  {
    return (*this)[a_domain]->getEBNormalData_x();
  }

  ///Get eb moment in the y direction
  inline shared_ptr<ebnorm_distrib_t >
  getEBNormalData_y(const Box& a_domain)
  {
    return (*this)[a_domain]->getEBNormalData_y();
  }

  inline DisjointBoxLayout getDBL(const Box& a_domain)
  {
    return (*this)[a_domain]->getDBL();
  }
  
  ///Get eb moment in the z direction
  inline shared_ptr<ebnorm_distrib_t >
  getEBNormalData_z(const Box& a_domain)
  {
    return (*this)[a_domain]->getEBNormalData_z();
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


  inline int getLevel(const Box & a_domain) const
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
    if(!found)
    {
      MayDay::Error("did not find asked-for domain");
    }

    return retval;
  }



  GeometryService();

  GeometryService(const GeometryService& a_workshopin);

  void operator=(const GeometryService& a_workshopin);


};
#include "Chombo_GeometryServiceImplem.H"
#include "Chombo_NamespaceFooter.H"
#endif
