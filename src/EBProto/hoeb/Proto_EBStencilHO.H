
#ifndef _Proto_EBStencilHO_H_
#define _Proto_EBStencilHO_H_


#include "EBProto.H"
//#include "Proto_Normalizor.H"

namespace Proto
{
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
  class EBStencilHO<CELL, CELL, order, double>
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
                     const EBGraph       & a_graph);

    inline static void 
    printWLSStencil(int                                  radius,
                    const LocalStencil< CELL,  double> & a_stencil,
                    const Point                        & a_vofPt); 
   
    ///   returns true if there are domain boundary conditions
    inline static bool
    getStencil(vector<EBIndex<CELL> >                    & a_dstVoFs,
                  vector<LocalStencil<CELL, double> >       & a_stencil,
                  vector<double>                            & a_inhoContr,
                  Stencil<double>                           & a_regStencilInterior,
                  vector<Stencil<double> >                    a_regStencilBCS[2*DIM],
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
                  double                                      a_dx,
                  bool                                        a_periodiconly,
                  Point                                       a_dstGrown);

    ///Gets WLS Poisson stencil
    inline static bool
    getWLSPoissonStencil(vector<EBIndex<CELL> >                    & a_dstVoFs,
                         vector<LocalStencil<CELL, double> >       & a_stencil,
                         vector<double>                            & a_inhoContr,
                         Stencil<double>                           & a_regStencilInterior,
                         vector<Stencil<double> >                    a_regStencilBC[2*DIM],
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
                         double                                      a_dx,
                         bool                                        a_periodiconly,
                         bool  a_forceAllNeumann);

    ///Gets WLS Advection stencil
    inline static bool
    getWLSAdvectionStencil(vector<EBIndex<CELL> >                    & a_dstVoFs,
                           vector<LocalStencil<CELL, double> >       & a_stencil,
                           vector<double>                            & a_inhoContr,
                           Stencil<double>                           & a_regStencilInterior,
                           vector<Stencil<double> >                    a_regStencilBC[2*DIM],
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
                           double                                      a_dx,
                           bool                                        a_periodiconly,
                           bool  a_forceAllNeumann);

  };

}

#include "Proto_EBStencilHOImplem.H"

#endif
  
