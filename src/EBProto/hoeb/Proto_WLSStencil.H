#ifndef _Proto_WLSStencil_H_
#define _Proto_WLSStencil_H_

#include "Proto_WLSCell.H"

namespace Proto
{
  ///some basic infrastructure to support the WLS stencils
  template<int order>  //this is so we can send in the geometry information
  class WLSStencil
  {
    typedef IndexedMoments<DIM  , order> IndMomDIM;
    typedef IndexedMoments<DIM-1, order> IndMomSDMinOne;
    typedef MomentIterator<DIM  , order> MomItDIM;
    typedef MomentIterator<DIM-1, order> MomItSDMinOne;
    typedef HostIrregData<CELL    ,  IndMomDIM , 1>  VoluData;
    typedef HostIrregData<BOUNDARY,  IndMomDIM , 1>  EBFaData;
    typedef HostIrregData<XFACE, IndMomSDMinOne, 1>  XFacData;
    typedef HostIrregData<YFACE, IndMomSDMinOne, 1>  YFacData;
    typedef HostIrregData<ZFACE, IndMomSDMinOne, 1>  ZFacData;
    typedef HostIrregData<CELL    ,  IndMomDIM , 1>  EBNormalData;
  public:

    ///
    //Gets Laplacian WLS flux divergence stencil
    inline void
    getLaplaWLSDivFStencil(LocalStencil< CELL,  Real>        & a_vofStencil,
                           const EBIndex<CELL>                 & a_vof,
                           const EBGraph                       & a_graph,
                           const VoluData                      & a_voludata,
                           const EBFaData                      & a_ebfadata,
                           const XFacData                      & a_xfacdata,
                           const YFacData                      & a_yfacdata,
                           const ZFacData                      & a_zfacdata,
                           const EBNormalData                  & a_ebnormxdata,
                           const EBNormalData                  & a_ebnormydata,
                           const EBNormalData                  & a_ebnormzdata,
                           const BaseWLSFlux<XFACE   , order>  & a_xfaceflux,
                           const BaseWLSFlux<YFACE   , order>  & a_yfaceflux,
                           const BaseWLSFlux<ZFACE   , order>  & a_zfaceflux,
                           const BaseWLSFlux<BOUNDARY, order>  & a_ebflux,
                           string                                a_ebbcName,
                           string                                a_dombcName[2*DIM],
                           Real                                & a_inhoTerm,
                           int                                 & a_neighType,
                           Real                                a_dx);

    //Gets Advection WLS flux divergence stencil
    inline void
    getAdvWLSDivFStencil(LocalStencil< CELL,  Real>        & a_vofStencil,
                         const EBIndex<CELL>                 & a_vof,
                         const EBGraph                       & a_graph,
                         const VoluData                      & a_voludata,
                         const EBFaData                      & a_ebfadata,
                         const XFacData                      & a_xfacdata,
                         const YFacData                      & a_yfacdata,
                         const ZFacData                      & a_zfacdata,
                         const EBNormalData                  & a_ebnormxdata,
                         const EBNormalData                  & a_ebnormydata,
                         const EBNormalData                  & a_ebnormzdata,
                         const ExactVelAdvectionWLSFlux<XFACE   , order>  & a_xfaceflux,
                         const ExactVelAdvectionWLSFlux<YFACE   , order>  & a_yfaceflux,
                         const ExactVelAdvectionWLSFlux<ZFACE   , order>  & a_zfaceflux,
                         const ExactVelAdvectionWLSFlux<BOUNDARY, order>  & a_ebflux,
                         string                                a_ebbcName,
                         string                                a_dombcName[2*DIM],
                         Real                                & a_inhoTerm,
                         int                                 & a_neighType,
                         Real                                a_dx);

    //Gets WLS stencil for each component the gradient of a scalar
    // return a vector of stencils, one per component
    inline void
    getWLSGradStencil(LocalStencil<CELL, Real>    & a_vofStencilx,
                      LocalStencil<CELL, Real>    & a_vofStencily,
                      LocalStencil<CELL, Real>    & a_vofStencilz,
                      const EBIndex<CELL>                 & a_vof,
                      const EBGraph                       & a_graph,
                      const VoluData                      & a_voludata,
                      const EBFaData                      & a_ebfadata,
                      const XFacData                      & a_xfacdata,
                      const YFacData                      & a_yfacdata,
                      const ZFacData                      & a_zfacdata,
                      const EBNormalData                  & a_ebnormxdata,
                      const EBNormalData                  & a_ebnormydata,
                      const EBNormalData                  & a_ebnormzdata,
                      const BaseWLSCell<order>  & a_cellConstruct,
                      string                                a_ebbcName,
                      string                                a_dombcName[2*DIM],
                      vector<Real>                                & a_inhoTerm,
                      int                                 & a_neighType,
                      Real                                a_dx);

        //Gets WLS flux divergence stencil in a direction
    // return a vector of stencils, one per component
    inline void
    getWLSDivFCompStencil(std::array<LocalStencil<CELL, Real>, SpaceDim >& a_vofStencil,
                          const EBIndex<CELL>                 & a_vof,
                          const EBGraph                       & a_graph,
                          const MomentHolder<order>           & a_momentData,
                          const DivWLSFlux<XFACE   , order>  & a_xfaceflux, // Why build these separate?
                          const DivWLSFlux<YFACE   , order>  & a_yfaceflux,
                          const DivWLSFlux<ZFACE   , order>  & a_zfaceflux,
                          const DivWLSFlux<BOUNDARY, order>  & a_ebflux,
                          string                                a_ebbcName,
                          vector<string>                                a_dombcName[2*DIM],
                          Real                                & a_inhoTerm,
                          int                                 a_neighType,
                          Real                                a_dx);

    template<int dir>
    void getWLSDivFCompStencilDir(
      std::array<LocalStencil<CELL, Real>, SpaceDim >& a_vofStencil,
      const EBIndex<CELL>                 & a_vof,
      const EBGraph                       & a_graph,
      const MomentHolder<order>           & a_momentData,
      const DivWLSFlux<XFACE   , order>  & a_xfaceflux, // Why build these separate?
      const DivWLSFlux<YFACE   , order>  & a_yfaceflux,
      const DivWLSFlux<ZFACE   , order>  & a_zfaceflux,
      const DivWLSFlux<BOUNDARY, order>  & a_ebflux,
      string                                a_ebbcName,
      vector<string>                        a_dombcName[2*DIM],
      Real                                & a_inhoTerm,
      int                                 a_neighType,
      Real                                a_dx);

    void
    getDistanceWeightedFluxAvgCoeffs(
      Real&    a_loSideWeight,
      Real&    a_hiSideWeight,
      const RealVect        a_faceCentroid,
      const EBIndex<CELL> & a_voflo,
      const EBIndex<CELL> & a_vofhi,
      const int             a_dir,
      const VoluData       & a_voludata,
      const EBGraph        & a_graph,
      const Real             a_dx);

    // Get a single face flux stencil at a face by a distance weighted average between
    // two cell-centered flux stencils states computed by the neighbor cells sharing the face
    // This was originally designed for the Laplacian stencil, but is generally applicable
    LocalStencil< CELL,  Real>
    getDistanceWeightedFluxAvgStencil(
      LocalStencil< CELL,  Real> & a_fluxMStencil,
      LocalStencil< CELL,  Real> & a_fluxPStencil,
      const RealVect                 a_centroid,
      const EBIndex<CELL>          & a_voflo,
      const EBIndex<CELL>          & a_vofhi,
      const int                      a_dir,
      const VoluData               & a_voludata,
      const EBGraph                & a_graph,
      const Real                     a_dx);

    Real
    getLapinhomogTerm(
      Real                         a_inhomogTermM,
      Real                         a_inhomogTermP,
      const RealVect               a_centroid,
      const EBIndex<CELL>          & a_voflo,
      const EBIndex<CELL>          & a_vofhi,
      const int                      a_dir,
      const VoluData               & a_voludata,
      const EBGraph                & a_graph,
      const Real                     a_dx);

    inline
    LocalStencil< CELL,  Real>
    getAdvFluxStencil(LocalStencil< CELL,  Real> & a_fluxMStencil,
                      LocalStencil< CELL,  Real> & a_fluxPStencil,
                      RealVect                     & a_centroid,
                      const EBIndex<CELL>          & a_voflo,
                      const EBIndex<CELL>          & a_vofhi,
                      const RealVect               & a_v_faceLoc,
                      int                            a_dir,
                      const VoluData               & a_voludata,
                      const EBGraph                & a_graph,
                      Real                         a_dx);

    virtual ~WLSStencil()
    {
    }

    WLSStencil()
    {
    }

    std::map<EBIndex<XFACE>,LocalStencil<CELL, Real> > xfacesFluxStencil;
    std::map<EBIndex<YFACE>,LocalStencil<CELL, Real> > yfacesFluxStencil;
    std::map<EBIndex<ZFACE>,LocalStencil<CELL, Real> > zfacesFluxStencil;
    std::map<EBIndex<XFACE>,Real> xfacesInhTerm;
    std::map<EBIndex<YFACE>,Real> yfacesInhTerm;
    std::map<EBIndex<ZFACE>,Real> zfacesInhTerm;
  };
} // end proto namespace
#include "Proto_WLSStencilImplem.H"
#endif _Proto_WLSStencil_H_
