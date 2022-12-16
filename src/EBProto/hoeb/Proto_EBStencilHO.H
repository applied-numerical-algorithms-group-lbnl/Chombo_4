#ifndef _Proto_EBStencilHO_H_
#define _Proto_EBStencilHO_H_

#include "EBProto.H"
#include "Proto_EBExactSolutions.H"
//#include "Proto_Normalizor.H"
#include "Chombo_ParmParse.H"

namespace Proto
{
  struct stencilType
  {
    // How many cells around the stencil to grab
    int cellRadius = 2;
    // The weighting power to use
    int powerExponent = -3;
    // Neighbor list type. 0-FaceBased or 1-CellBased
    int neighListType = 1;

    void
    readStencilType(const std::string a_PPprefix)
      {
        ParmParse pp;
        pp.get("radius", cellRadius);
        pp.get("powerE", powerExponent);
        pp.query("neighListType", neighListType);
      }
  };

  template <CENTERING srcCenter, CENTERING dstCenter, int order, typename data_t>
  class EBStencilHO
  {
  public:
    EBStencilHO() {}
    ~EBStencilHO() {}

    typedef IndexedMoments<DIM  , order> IndMomDIM;
    typedef IndexedMoments<DIM-1, order> IndMomSDMinOne;
    typedef MomentIterator<DIM  , order> MomItDIM;
    typedef MomentIterator<DIM-1, order> MomItSDMinOne;
    typedef HostIrregData<CELL,     IndMomDIM  , 1>  VoluData;
    typedef HostIrregData<BOUNDARY, IndMomDIM  , 1>  EBFaData;
    typedef HostIrregData<XFACE, IndMomSDMinOne, 1>  XFacData;
    typedef HostIrregData<YFACE, IndMomSDMinOne, 1>  YFacData;
    typedef HostIrregData<ZFACE, IndMomSDMinOne, 1>  ZFacData;
    typedef HostIrregData<CELL,      IndMomDIM , 1>  EBNormalData;

  };

  ///
  /**
   * EBStencil for high-order methods class. First implementation will get you
   * the corrresponding stencil for Poisson operator.
   */
  template <int order>
  class EBStencilHO<CELL, CELL, order, Real>
  {
  public:
    typedef IndexedMoments<DIM  , order> IndMomDIM;
    typedef IndexedMoments<DIM-1, order> IndMomSDMinOne;
    typedef MomentIterator<DIM  , order> MomItDIM;
    typedef MomentIterator<DIM-1, order> MomItSDMinOne;
    typedef HostIrregData<CELL,      IndMomDIM , 1>  VoluData;
    typedef HostIrregData<BOUNDARY,  IndMomDIM , 1>  EBFaData;
    typedef HostIrregData<XFACE, IndMomSDMinOne, 1>  XFacData;
    typedef HostIrregData<YFACE, IndMomSDMinOne, 1>  YFacData;
    typedef HostIrregData<ZFACE, IndMomSDMinOne, 1>  ZFacData;
    typedef HostIrregData<CELL,      IndMomDIM , 1>  EBNormalData;

    EBStencilHO() {}
    ~EBStencilHO() {}

    inline static vector<EBIndex<CELL> >
    getFlapVoFs(const EBGraph& a_graph,  EBIndex<CELL>& a_start); 

    inline static bool 
    isRegularStencil(const EBIndex<CELL> & a_vof,
                     const int             a_radius,
                     const EBGraph       & a_graph,
                     const VoluData      & a_voludata,
                     const Real          a_dx);

    inline static void 
    printWLSStencil(int                                  radius,
                    const LocalStencil< CELL,  Real> & a_stencil,
                    const Point                        & a_vofPt); 
   
    ///   returns true if there are domain boundary conditions
    inline static bool
    getStencil(vector<EBIndex<CELL> >                    & a_dstVoFs,
               vector<LocalStencil<CELL, Real> >       & a_stencil,
               vector<Real>                            & a_inhoContr,
               Stencil<Real>                           & a_regStencilInterior,
               vector<Stencil<Real> >                    a_regStencilBCS[2*DIM],
               vector<Box>                                 a_BCApplyBoxes[2*DIM],
               Box                                       & a_regApplyBox,
               bool                                      & a_bcOnly,
               string                                      a_stencilName,
               string                                      a_domainBCName[2*DIM],
               string                                      a_ebbcName,
               const Box                                 & a_srcValid,
               const Box                                 & a_dstValid,
               const Box                                 & a_srcDomain,
               const Box                                 & a_dstDomain,
               const Point                               & a_srcGhost,
               const Point                               & a_dstGhost,
               const EBGraph                             & a_srcGraph,
               const EBGraph                             & a_dstGraph,
               const VoluData                            & a_voludata,
               const EBFaData                            & a_ebfadata,
               const XFacData                            & a_xfacdata,
               const YFacData                            & a_yfacdata,
               const ZFacData                            & a_zfacdata,
               const EBNormalData                        & a_ebnormxdata,
               const EBNormalData                        & a_ebnormydata,
               const EBNormalData                        & a_ebnormzdata,
               Real                                      a_dx,
               bool                                        a_periodiconly,
               Point                                       a_dstGrown,
               const vector<shared_ptr<BaseExactSolution<order> > > & a_boundarySolutions = vector<shared_ptr<BaseExactSolution<order> > >(),
               const shared_ptr<BaseExactSolution<order> > & a_boundarySolutionsEB = nullptr);

    ///Gets WLS Poisson stencil
    inline static bool
    getWLSPoissonStencil(vector<EBIndex<CELL> >                    & a_dstVoFs,
                         vector<LocalStencil<CELL, Real> >       & a_stencil,
                         vector<Real>                            & a_inhoContr,
                         Stencil<Real>                           & a_regStencilInterior,
                         vector<Stencil<Real> >                    a_regStencilBC[2*DIM],
                         vector<Box>                                 a_BCApplyBoxes[2*DIM],
                         Box                                       & a_regApplyBox,
                         string                                      a_stencilName,
                         string                                      a_domainBCName[2*DIM],
                         string                                      a_ebbcName,
                         const Box                                 & a_validBox,
                         const Box                                 & a_domain,
                         const Point                               & a_srcGhost,
                         const Point                               & a_dstGhost,
                         const EBGraph                             & a_graph,
                         const VoluData                            & a_voludata,
                         const EBFaData                            & a_ebfadata,
                         const XFacData                            & a_xfacdata,
                         const YFacData                            & a_yfacdata,
                         const ZFacData                            & a_zfacdata,
                         const EBNormalData                        & a_ebnormxdata,
                         const EBNormalData                        & a_ebnormydata,
                         const EBNormalData                        & a_ebnormzdata,
                         Real                                      a_dx,
                         bool                                        a_periodiconly,
                         bool  a_forceAllNeumann,
                         const vector<shared_ptr<BaseExactSolution<order> > > & a_boundarySolutions = vector<shared_ptr<BaseExactSolution<order> > >(),
                         const shared_ptr<BaseExactSolution<order> > & a_boundarySolutionsEB = nullptr);

    ///Gets WLS Advection stencil
    inline static bool
    getWLSAdvectionStencil(vector<EBIndex<CELL> >                    & a_dstVoFs,
                           vector<LocalStencil<CELL, Real> >       & a_stencil,
                           vector<Real>                            & a_inhoContr,
                           Stencil<Real>                           & a_regStencilInterior,
                           vector<Stencil<Real> >                    a_regStencilBC[2*DIM],
                           vector<Box>                                 a_BCApplyBoxes[2*DIM],
                           Box                                       & a_regApplyBox,
                           string                                      a_stencilName,
                           string                                      a_domainBCName[2*DIM],
                           string                                      a_ebbcName,
                           const Box                                 & a_validBox,
                           const Box                                 & a_domain,
                           const Point                               & a_srcGhost,
                           const Point                               & a_dstGhost,
                           const EBGraph                             & a_graph,
                           const VoluData                            & a_voludata,
                           const EBFaData                            & a_ebfadata,
                           const XFacData                            & a_xfacdata,
                           const YFacData                            & a_yfacdata,
                           const ZFacData                            & a_zfacdata,
                           const EBNormalData                        & a_ebnormxdata,
                           const EBNormalData                        & a_ebnormydata,
                           const EBNormalData                        & a_ebnormzdata,
                           Real                                      a_dx,
                           bool                                        a_periodiconly,
                           bool  a_forceAllNeumann);

    ///Gets WLS Divergence stencil
    inline static bool
    getWLSDivergenceCompStencil(vector<EBIndex<CELL> >                    & a_dstVoFs, // could this be tied to IVS or similar?
                                vector<vector<LocalStencil<CELL, Real> > >      & a_stencil, // should be an space dim sized array, tied to IVS
                                vector<Real>                            & a_inhoContr,
                                Stencil<Real>                           & a_regStencilInterior, // actually should use this
                                vector<Stencil<Real> >                    a_regStencilBC[2*DIM], // not used
                                vector<Box>                                 a_BCApplyBoxes[2*DIM], // not used
                                Box                                       & a_regApplyBox,
                                string                                      a_stencilName, // why?
                                vector<string>                                      a_domainBCName[2*DIM], // how is this ordered?
                                string                                      a_ebbcName,
                                const Box                                 & a_validBox,
                                const Box                                 & a_domain,
                                const Point                               & a_srcGhost, // should this matter?
                                const Point                               & a_dstGhost, // should this matter?
                                const EBGraph                             & a_graph,
                                //const LevelGeometryService<order>         & a_levelGeoServ,
                                const VoluData                            & a_voludata,  // Really I am just asking for all the data stored in geometry service on a level
                                const EBFaData                            & a_ebfadata,
                                const XFacData                            & a_xfacdata,
                                const YFacData                            & a_yfacdata,
                                const ZFacData                            & a_zfacdata,
                                const EBNormalData                        & a_ebnormxdata,
                                const EBNormalData                        & a_ebnormydata,
                                const EBNormalData                        & a_ebnormzdata, // last geoserv data
                                Real                                      a_dx, // should be part of apply
                                const vector<shared_ptr<BaseExactSolution<order> > >& a_boundarySolution = vector<shared_ptr<BaseExactSolution<order> > >(), // these are a mess, should be an array or struct to help with organization, similar to BCName. Ideally part of apply
                                const shared_ptr<BaseExactSolution<order> > & a_boundarySolutionsEB = nullptr);

    ///Gets WLS Divergence stencil
    inline static bool
    getWLSGradientStencil(vector<EBIndex<CELL> >                    & a_dstVoFs,
                          vector<vector<LocalStencil<CELL, Real> > >      & a_stencil,
                          vector<vector<Real> >                            & a_inhoContr,
                          Stencil<Real>                           & a_regStencilInterior,
                          vector<Stencil<Real> >                    a_regStencilBC[2*DIM],
                          vector<Box>                                 a_BCApplyBoxes[2*DIM],
                          Box                                       & a_regApplyBox,
                          string                                      a_stencilName,
                          string                                      a_domainBCName[2*DIM],
                          string                                      a_ebbcName,
                          const Box                                 & a_validBox,
                          const Box                                 & a_domain,
                          const Point                               & a_srcGhost,
                          const Point                               & a_dstGhost,
                          const EBGraph                             & a_graph,
                          const VoluData                            & a_voludata,
                          const EBFaData                            & a_ebfadata,
                          const XFacData                            & a_xfacdata,
                          const YFacData                            & a_yfacdata,
                          const ZFacData                            & a_zfacdata,
                          const EBNormalData                        & a_ebnormxdata,
                          const EBNormalData                        & a_ebnormydata,
                          const EBNormalData                        & a_ebnormzdata,
                          Real                                      a_dx,
                          bool                                        a_periodiconly,
                          bool  a_forceAllNeumann,
                          const shared_ptr<BaseExactSolution<order> >& a_boundarySolution = nullptr);

  };

}

#include "Proto_EBStencilHOImplem.H"

#endif
  
