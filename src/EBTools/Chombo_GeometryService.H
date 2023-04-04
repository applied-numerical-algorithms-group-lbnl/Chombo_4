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

// #include "Proto.H"
// #include "EBProto.H"
#include "Chombo_EBChombo.H"
#include "Proto_MomentHolder.H"
#include "Chombo_LevelData.H"
#include "Chombo_EBDataChoreography.H"
#include "Chombo_AMRIO.H"
#include "Chombo_ProtoInterface.H"
#include "Chombo_BRMeshRefine.H"
#include "Chombo_LoadBalance.H"
#include "Chombo_ParmParse.H"
#include <memory>
#include "Chombo_NamespaceHeader.H"

// Used in a hack to allow indexing dimensional components again. This save a lot of code duplication higher up
// This shouldn't be needed, but the fix needs to go really low
template<CENTERING>
struct typenameTag {};


///
/**
   makes geometric description from an implicit function.
*/

template <int order>
class LevelGeometryService
{

public:
  typedef Proto::IndexTM<int,2>      Iv2;
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
                       const int                    &      a_ghost,
                       bool a_printStuff = false);




  LevelGeometryService(const LevelGeometryService   &      a_finerLevel,
                       const DisjointBoxLayout      &      a_grids,
                       bool a_printStuff = false);

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
  inline Real dx() const
  {
    return m_dx;
  }
  
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

  const shared_ptr<LevelData<HostIrregData<XFACE, IndMomSDMinOne, 1> > >
  getFaceData_dispatch(typenameTag<XFACE>) const;

  const shared_ptr<LevelData<HostIrregData<YFACE, IndMomSDMinOne, 1> > >
  getFaceData_dispatch(typenameTag<YFACE>) const;

  const shared_ptr<LevelData<HostIrregData<ZFACE, IndMomSDMinOne, 1> > >
  getFaceData_dispatch(typenameTag<ZFACE>) const;

  /// ONE ACCESSOR TO RULE THEM ALL
  template<CENTERING cent>
  const shared_ptr<LevelData<HostIrregData<cent, IndMomSDMinOne, 1> > >
  getFaceData() const
  {
    return getFaceData_dispatch(typenameTag<cent>{});
    //MayDay::Error("invalid face centering");
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

  ///Get eb moment in the
  inline const shared_ptr<LevelData<EBNormalData> >
  getEBNormalData(const int a_dir) const
  {
    CH_assert((a_dir >= 0) && (a_dir < SpaceDim));
    return m_ebnormaldata[a_dir];
  }

  // A struct for holding all the moments,
  // now we can pass just 1 struct instead of 10 arguments and hide the problem dimensionality
  // This only has proto dependencies, so usable in stencils
  const
  Proto::MomentHolder<order>
  getAllMomentData(const DataIndex a_dit) const
    {
      Proto::MomentHolder<order> holder((*(getVoluData()))[a_dit],
                                        (*(getEBFaceData()))[a_dit],
                                        (*(getXFaceData()))[a_dit],
                                        (*(getYFaceData()))[a_dit],
                                        (*(getZFaceData()))[a_dit],
                                        (*(getEBNormalData_x()))[a_dit],
                                        (*(getEBNormalData_y()))[a_dit],
                                        (*(getEBNormalData_z()))[a_dit]);
      return holder;
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


  /// better name---fills with area fraction
  inline void fillEBFaceAlpha(EBHostData<CELL,Real, 1> & a_ebfacedata,
                              const Box                & a_grbx ,
                              const DataIndex          & a_dit) const
  { 
    const EBFaData & ebfmo = (*m_ebfadata)[a_dit];
    const EBGraph  & graph = (*m_graphs)[a_dit];
    HostBoxData<        Real, 1>& reghost = a_ebfacedata.getRegData();
    HostIrregData<CELL, Real, 1>& irrhost = a_ebfacedata.getIrrData();
    Real fullCellArea = D_TERM(1., *m_dx, *m_dx);

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
          const EBIndex<BOUNDARY> ebf = vof.getCutFace();
          const IndMomSDMinOne&  momspt = ebfmo(ebf, 0);
          double ebalpha = momspt[0]/fullCellArea;
          reghost(vof.m_pt, 0) = ebalpha;
          irrhost(vof     , 0) = ebalpha;
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
                    const double                 & a_dx ,
                    bool a_printStuff = false);

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

  //makes levelgeometry service on refined coarse layout.  Used for agglomeration
  LevelGeometryService(const LevelGeometryService   & a_fineLevel,
                       const LevelGeometryService   & a_coarLevel,
                       bool a_printStuff = false);
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
  void exchange(bool a_printStuff = false);

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
template <int order >
inline void
fillFaceAlpha(EBHostData<XFACE, Real, 1>            & a_alpha,
              const Box                             & a_grbx ,
              const DataIndex                       & a_dit,
              const LevelGeometryService<order>     & a_geo)
{
  typedef IndexedMoments<DIM  , order> IndMomDIM;
  typedef IndexedMoments<DIM-1, order> IndMomDMO;
  typedef MomentIterator<DIM  , order> MomItDIM;
  typedef MomentIterator<DIM-1, order> MomItSDMinOne;
  Real dx = a_geo.dx();
  Real fullCellArea = D_TERM(1., *dx, *dx);
  auto dataptr = a_geo.getXFaceData();
  auto graphs  = a_geo.getGraphs();
  const auto & irrdat = (*dataptr)[a_dit];
  const auto & graph  =  (*graphs)[a_dit];
  HostBoxData<         Real, 1>& reghost = a_alpha.getRegData();
  HostIrregData<XFACE, Real, 1>& irrhost = a_alpha.getIrrData();
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
        for(SideIterator sit; sit.ok(); ++sit)
        {
          auto vec_faces = graph.getXFaces(vof, sit());
          for(int iface  = 0; iface < vec_faces.size(); iface++)
          {
            const auto & face = vec_faces[iface];
            const IndMomDMO&  momspt = irrdat(face, 0);
            double facealpha = momspt[0]/fullCellArea;

            reghost(face.m_pt, 0) = facealpha;
            irrhost(face     , 0) = facealpha;
          }
        }
      }
    }
  }
}

template<int order>
const shared_ptr<LevelData<HostIrregData<XFACE, IndexedMoments<DIM-1, order>, 1> > >
LevelGeometryService<order>::getFaceData_dispatch(typenameTag<XFACE>) const
{
  return this->getXFaceData();
}

template<int order>
const shared_ptr<LevelData<HostIrregData<YFACE, IndexedMoments<DIM-1, order>, 1> > >
LevelGeometryService<order>::getFaceData_dispatch(typenameTag<YFACE>) const
{
  return this->getYFaceData();
}

template<int order>
const shared_ptr<LevelData<HostIrregData<ZFACE, IndexedMoments<DIM-1, order>, 1> > >
LevelGeometryService<order>::getFaceData_dispatch(typenameTag<ZFACE>) const
{
  return this->getZFaceData();
}

template <int order >
inline void
fillFaceAlpha(EBHostData<YFACE, Real, 1>        & a_alpha,
              const Box                         & a_grbx ,
              const DataIndex                   & a_dit,
              const LevelGeometryService<order> & a_geo)
{
  typedef IndexedMoments<DIM  , order> IndMomDIM;
  typedef IndexedMoments<DIM-1, order> IndMomDMO;
  typedef MomentIterator<DIM  , order> MomItDIM;
  typedef MomentIterator<DIM-1, order> MomItSDMinOne;
  Real dx = a_geo.dx();
  Real fullCellArea = D_TERM(1., *dx, *dx);
  auto dataptr = a_geo.getYFaceData();
  auto graphs  = a_geo.getGraphs();
  const auto & irrdat = (*dataptr)[a_dit];
  const auto & graph  =  (*graphs)[a_dit];
  HostBoxData<         Real, 1>& reghost = a_alpha.getRegData();
  HostIrregData<YFACE, Real, 1>& irrhost = a_alpha.getIrrData();
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
        for(SideIterator sit; sit.ok(); ++sit)
        {
          auto vec_faces = graph.getYFaces(vof, sit());
          for(int iface  = 0; iface < vec_faces.size(); iface++)
          {
            const auto & face = vec_faces[iface];
            const IndMomDMO&  momspt = irrdat(face, 0);
            double facealpha = momspt[0]/fullCellArea;

            reghost(face.m_pt, 0) = facealpha;
            irrhost(face     , 0) = facealpha;
          }
        }
      }
    }
  }
}

#if DIM==3

template <int order >
inline void
fillFaceAlpha(EBHostData<ZFACE, Real, 1>        & a_alpha,
              const Box                         & a_grbx ,
              const DataIndex                   & a_dit,
              const LevelGeometryService<order> & a_geo)
{
  typedef IndexedMoments<DIM  , order> IndMomDIM;
  typedef IndexedMoments<DIM-1, order> IndMomDMO;
  typedef MomentIterator<DIM  , order> MomItDIM;
  typedef MomentIterator<DIM-1, order> MomItSDMinOne;
  Real dx = a_geo.dx();
  Real fullCellArea = D_TERM(1., *dx, *dx);
  auto dataptr = a_geo.getZFaceData();
  auto graphs  = a_geo.getGraphs();
  const auto & irrdat = (*dataptr)[a_dit];
  const auto & graph  =  (*graphs)[a_dit];
  HostBoxData<         Real, 1>& reghost = a_alpha.getRegData();
  HostIrregData<XFACE, Real, 1>& irrhost = a_alpha.getIrrData();
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
        for(SideIterator sit; sit.ok(); ++sit)
        {
          auto vec_faces = graph.getZFaces(vof, sit());
          for(int iface  = 0; iface < vec_faces.size(); iface++)
          {
            const auto & face = vec_faces[iface];
            const IndMomDMO&  momspt = irrdat(face, 0);
            double facealpha = momspt[0]/fullCellArea;

            reghost(face.m_pt, 0) = facealpha;
            irrhost(face     , 0) = facealpha;
          }
        }
      }
    }
  }
}
#endif

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
                  const std::vector<DisjointBoxLayout>   &      a_grids,
                  const int                         &      a_ghost,
                  bool a_printStuff = false)
  {
    int numLevels = a_grids.size();
    m_levels.resize( numLevels);
    m_domains.resize(numLevels);
    if(a_printStuff)
    {
      pout() << "GeometryService constuctor: making finest level " << endl;
    }
    m_levels[0]  =
      shared_ptr<LevelGeometryService<order> >(
        new LevelGeometryService<order>(a_baseIF, a_origin, a_finestDx,
                                        a_finestDomain, a_grids[0],
                                        a_ghost, a_printStuff));
    m_domains[0] = a_finestDomain;
    for(int ilev = 1; ilev < numLevels; ilev++)
    {
      if(a_printStuff)
      {
        pout() << "GeometryService constuctor: making level << ilev " << endl;
      }
      m_levels[ilev]  = shared_ptr<LevelGeometryService<order> >(
        new LevelGeometryService<order>(*m_levels[ilev-1],
                                        a_grids[ilev],
                                        a_printStuff));
      
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
  generateGrids(std::vector<DisjointBoxLayout> & a_grids,
                const Box                 & a_domain,
                int a_maxGridSize =  32,
                int a_maxCoarsen  = -1)
  {
    std::vector<Box> boxes;
    std::vector<int> procs;
    domainSplit(a_domain, boxes, a_maxGridSize);
    LoadBalance(procs, boxes);
    DisjointBoxLayout dblFine(boxes, procs, a_domain);
    a_grids.push_back(dblFine);

    bool doAgglom = true;
    ParmParse pp("geometry_service");
    pp.query("do_agglomeration", doAgglom);

    //while((currDBL.coarsenable(4)) && ((a_grids.size()-1)!=a_maxCoarsen))
    if(doAgglom)
    {
      int ilev = 0;
      Box domCoar = coarsen(a_domain, 2);
      bool continueCoarsening = (domCoar.coarsenable(2)) && ((a_maxCoarsen < 0) || (ilev < a_maxCoarsen));
      //a_grids[0].print();
      while(continueCoarsening) 
      {
        DisjointBoxLayout dblCoar;
        if(dblFine.coarsenable(4))
        {
          coarsen(dblCoar, dblFine, 2);
        }
        else
        {
          domainSplit(domCoar, boxes, a_maxGridSize);
          LoadBalance(procs, boxes);
          dblCoar = DisjointBoxLayout(boxes, procs, domCoar);
        }
        
        a_grids.push_back(dblCoar);
        
        dblFine = dblCoar;
        ilev ++;
        domCoar = coarsen(domCoar, 2);
        continueCoarsening = (domCoar.coarsenable(2)) && ((a_maxCoarsen < 0) || (ilev < a_maxCoarsen));
      }
    }
    else
    {
      //a_grids[0].print();
      while(dblFine.coarsenable(4))
      {
        DisjointBoxLayout dblCoar;
        coarsen(dblCoar, dblFine, 2);
        a_grids.push_back(dblCoar);
        dblFine = dblCoar;
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

  /// fill a data holder in with eb face area fractions
  inline void
  fillEBFaceAlpha(EBHostData<CELL,Real, 1>& a_ebfdata,
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
  
  ///Get eb moment in the z direction
  inline shared_ptr<ebnorm_distrib_t >
  getEBNormalData_z(const Box& a_domain)
  {
    return (*this)[a_domain]->getEBNormalData_z();
  }

  inline DisjointBoxLayout getDBL(const Box& a_domain)
  {
    return (*this)[a_domain]->getDBL();
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


  bool hasLevel(const Box & a_domain) const
  {
    bool found = false;
    for(int ilev = 0; ilev < m_domains.size(); ilev++)
    {
      if(m_domains[ilev] == a_domain)
      {
        found = true;
        break;
      }
    }

    return found;

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

template <CENTERING cent, int order>
inline shared_ptr<CH4_Data_Choreography::DistributedData<HostIrregData<cent, IndexedMoments<DIM-1, order> , 1> > >
getFaceData(const shared_ptr< GeometryService<order> > & a_geoserv,
            const Box& a_domain, const EBIndex<cent>& a_dummy)
{
  shared_ptr<CH4_Data_Choreography::DistributedData<HostIrregData<cent, IndexedMoments<DIM-1, order> , 1> > >   retval;
  Chombo4::MayDay::Error("did not find template specialization");
  return retval;
}

template <int order  >
inline shared_ptr<CH4_Data_Choreography::DistributedData<HostIrregData<XFACE, IndexedMoments<DIM-1, order> , 1> > >
getFaceData(const shared_ptr< GeometryService<order> >  & a_geoserv,
            const Box& a_domain, const EBIndex<XFACE>& a_dummy)
{
  auto  retval = a_geoserv->getXFaceData(a_domain);;
  return retval;
}
template <int order  >
inline shared_ptr<CH4_Data_Choreography::DistributedData<HostIrregData<YFACE, IndexedMoments<DIM-1, order> , 1> > >
getFaceData(const shared_ptr< GeometryService<order> >  & a_geoserv,
            const Box& a_domain, const EBIndex<YFACE>& a_dummy)
{
  auto  retval = a_geoserv->getYFaceData(a_domain);;
  return retval;
}

#if DIM==3
template <int order  >
inline shared_ptr<CH4_Data_Choreography::DistributedData<HostIrregData<ZFACE, IndexedMoments<DIM-1, order> , 1> > >
getFaceData(const shared_ptr< GeometryService<order> >  & a_geoserv,
            const Box& a_domain, const EBIndex<ZFACE>& a_dummy)
{
  auto  retval = a_geoserv->getZFaceData(a_domain);;
  return retval;
}

#endif
#include "Chombo_GeometryServiceImplem.H"
#include "Chombo_NamespaceFooter.H"
#endif
