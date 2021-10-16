
#ifndef _Proto_SecondOrderCell_H_
#define _Proto_SecondOrderCell_H_


#include "EBProto.H"

namespace Proto
{
  template<int order> //order is so we can send in geometry info
  class NeumannEBFlux
  {
  public:
    typedef IndexedMoments<DIM  , order>         IndMomDIM;
    typedef IndexedMoments<DIM-1, order>         IndMomSDMinOne;
    typedef HostIrregData<CELL    ,  IndMomDIM , 1>  VoluData;
    typedef HostIrregData<BOUNDARY,  IndMomDIM , 1>  EBFaData;
    typedef HostIrregData<XFACE, IndMomSDMinOne, 1>  XFacData;
    typedef HostIrregData<YFACE, IndMomSDMinOne, 1>  YFacData;
    typedef HostIrregData<ZFACE, IndMomSDMinOne, 1>  ZFacData;

    NeumannEBFlux()
    {
    }

    virtual ~NeumannEBFlux()
    {
    }

    inline static void
    getEBFluxStencil(LocalStencil<CELL, double>    &  a_facesten,
                     double                        &  a_weight, //for inhomogeneous bcs
                     const EBIndex<BOUNDARY>       &  a_face,
                     const EBGraph                 & a_graph,
                     const VoluData                & a_voludata,
                     const EBFaData                & a_ebfadata,
                     const XFacData                & a_xfacdata,
                     const YFacData                & a_yfacdata,
                     const ZFacData                & a_zfacdata,
                     const double                  & a_dx) 
    {

      //neumann has no stencil
    }
 
  };


  template<int order> //order is so we can send in geometry info
  class DirichletEBFlux
  {
  public:
    typedef IndexedMoments<DIM  , order>         IndMomDIM;
    typedef IndexedMoments<DIM-1, order>         IndMomSDMinOne;
    typedef HostIrregData<CELL    ,  IndMomDIM , 1>  VoluData;
    typedef HostIrregData<BOUNDARY,  IndMomDIM , 1>  EBFaData;
    typedef HostIrregData<XFACE, IndMomSDMinOne, 1>  XFacData;
    typedef HostIrregData<YFACE, IndMomSDMinOne, 1>  YFacData;
    typedef HostIrregData<ZFACE, IndMomSDMinOne, 1>  ZFacData;

    DirichletEBFlux()
    {
    }

    virtual ~DirichletEBFlux()
    {
    }

    inline static void
    getEBFluxStencil(LocalStencil<CELL, double>    &  a_facesten,
                     double                        &  a_weight, //for inhomogeneous bcs
                     const EBIndex<BOUNDARY>       &  a_face,
                     const EBGraph                 & a_graph,
                     const VoluData                & a_voludata,
                     const EBFaData                & a_ebfadata,
                     const XFacData                & a_xfacdata,
                     const YFacData                & a_yfacdata,
                     const ZFacData                & a_zfacdata,
                     const double                  & a_dx); 

    
    inline static void
    calculateWeightingMatrix(RealVect           x0,
                             vector<RealVect>&  xp,
                             vector<RealVect>&  weightMatrix,
                             bool&              detZero);

  };

  ///stuff in common for all second order flux calculations
  template<CENTERING cent, int order>
  class Base2ndOrderFlux
  {
  public:
    typedef IndexedMoments<DIM  , order> IndMomDIM;
    typedef IndexedMoments<DIM-1, order> IndMomSDMinOne;
    typedef MomentIterator<DIM  , order> MomItDIM;
    typedef MomentIterator<DIM-1, order> MomItSDMinOne;
    
    virtual void
    getFaceCenteredFluxStencil(LocalStencil<CELL, double>                   &  a_facesten,
                               const EBIndex<cent>                          &  a_face,
                               const EBGraph                                &  a_graph,
                               const HostIrregData<cent, IndMomSDMinOne, 1> &  a_facdata,
                               string                                          a_dombc,
                               double                                          a_dx) const
    {
      //sometimes need other parts of this class
      PROTO_ASSERT(false, "not implemented");
    }

    void
    getFaceInterpStencil(LocalStencil<cent, double>                   &  a_sten,
                         const EBIndex<cent>                          &  a_face,
                         const EBGraph                                &  a_graph,
                         const HostIrregData<cent, IndMomSDMinOne, 1> &  a_facdata,
                         double                                          a_dx) const;

    void
    getCentroidStencil(LocalStencil<CELL, double>                   &  a_sten,
                       const EBIndex<cent>                          &  a_face,
                       const EBGraph                                &  a_graph,
                       const HostIrregData<cent, IndMomSDMinOne, 1> &  a_facdata,
                       string                                          a_dombc,
                       double                                          a_dx) const;

    
  };
  
  ///second order Poisson flux
  template<CENTERING facecent, int order>
  class Poisson2ndOrder: public Base2ndOrderFlux<facecent, order>
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
  public:
    inline void 
    getFaceCenteredFluxStencil(LocalStencil<CELL, double>                      & a_facesten,
                               const EBIndex<facecent>                         & a_face,
                               const EBGraph                                   & a_graph,
                               const HostIrregData<facecent, IndMomSDMinOne, 1>& a_facedat,
                               string                                            a_dombc,
                               double                                            a_dx) const;
 
    virtual ~Poisson2ndOrder()
    {
    }

    Poisson2ndOrder()
    {
    }
  };
  ///some basic infrastructure to support the second order stuff
  template<int order>  //this is so we can send in the geometry informaion sanely
  class SecondOrderStencil
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
  public:

    ///
    inline static void
    get2ndOrderDivFStencil(LocalStencil< CELL,  double>             & a_vofStencil,
                           const EBIndex<CELL>                      & a_vof,
                           const EBGraph                            & a_graph,
                           const VoluData                           & a_voludata,
                           const EBFaData                           & a_ebfadata,
                           const XFacData                           & a_xfacdata,
                           const YFacData                           & a_yfacdata,
                           const ZFacData                           & a_zfacdata,
                           const Base2ndOrderFlux<XFACE   , order>  & a_xfaceflux,
                           const Base2ndOrderFlux<YFACE   , order>  & a_yfaceflux,
                           const Base2ndOrderFlux<ZFACE   , order>  & a_zfaceflux,
                           string                                  a_dombc[2*DIM],
                           string                                  a_ebbcName,
                           double                                     a_dx);




    inline static void
    getDomainBCs(vector<Stencil<double> >                    a_regStencilBC[2*DIM],       
                 vector<Box>                                 a_BCApplyBoxes[2*DIM],       
                 string                                      a_domainBCName[2*DIM],               
                 const Box                                 & a_validBox,
                 const Box                                 & a_domain,
                 const Point                               & a_srcGhost);
    
    virtual ~SecondOrderStencil()
    {
    }

    SecondOrderStencil()
    {
    }
  };



}
#include "Proto_SecondOrderCellImplem.H"

#endif
  
