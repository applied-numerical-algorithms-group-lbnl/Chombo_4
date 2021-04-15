
#ifndef _Proto_EBStencilArchive_H_
#define _Proto_EBStencilArchive_H_


#include "EBProto.H"
//#include "Proto_Normalizor.H"

namespace Proto
{
  ///
  /**
   */
  template <CENTERING srcCenter, CENTERING dstCenter, int order, typename data_t>
  class EBStencilArchive
  {
  public:
    EBStencilArchive() {}
    ~EBStencilArchive() {}
    
    typedef IndexedMoments<DIM  , order> IndMomDIM;
    typedef IndexedMoments<DIM-1, order> IndMomSDMinOne;
    typedef MomentIterator<DIM  , order> MomItDIM;
    typedef MomentIterator<DIM-1, order> MomItSDMinOne;
    typedef HostIrregData<CELL,     IndMomDIM  , 1>  VoluData;
    typedef HostIrregData<BOUNDARY, IndMomDIM  , 1>  EBFaData;
    typedef HostIrregData<XFACE, IndMomSDMinOne, 1>  XFacData;
    typedef HostIrregData<YFACE, IndMomSDMinOne, 1>  YFacData;
    typedef HostIrregData<ZFACE, IndMomSDMinOne, 1>  ZFacData;

    ///   returns true if there are domain boundary conditions (not all ops have them)
    inline static bool
    getStencil(vector<EBIndex<dstCenter> >               & a_dstVoFs,
               vector<LocalStencil<srcCenter,  data_t> > & a_stencil,
               Stencil<data_t>                           & a_regStencilInterior,
               vector<Stencil<data_t> >                    a_regStencilBCS[2*DIM],
               vector<Box>                                 a_BCApplyBoxes[2*DIM],       
               Box                                       & a_regApplyBox,
               bool                                      & a_bcOnly,
               bool                                      & a_irregOnly,
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
               double                                      a_dx,
               bool                                        a_periodiconly,
               Point                                       a_dstGrown)               
    {
      PROTO_ASSERT(false, "not implemented");
      return false;
    }
 
  };



    ///
  /**
     Stuff to go from faces to faces (same direction)
   */
  template <CENTERING fcent, int order>
  class EBStencilArchive<fcent, fcent, order, double>
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

    EBStencilArchive() {}
    ~EBStencilArchive() {}
    

    ///   returns true if there are domain boundary conditions
    inline static bool
    getStencil(vector<EBIndex<fcent> >                    & a_dstVoFs,
               vector<LocalStencil<fcent, double> >       & a_stencil,
               Stencil<double>                           & a_regStencil,
               vector<Stencil<double> >                    a_regStencilBCS[2*DIM],
               vector<Box>                                 a_BCApplyBoxes[2*DIM],       
               Box                                       & a_regApplyBox,
               bool                                      & a_bcOnly,
               bool                                      & a_irregOnly,
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
               double                                      a_dx,
               bool                                        a_periodiconly, 
               Point                                       a_dstGrown);               

    inline static bool
    getFaceCentroidInterpStencil(vector<EBIndex<fcent > >                  & a_dstFaces,
                                 vector<LocalStencil<fcent, double> >      & a_stencil,
                                 Stencil<double>                           & a_regStencil,
                                 Box                                       & a_regApplyBox,
                                 bool                                      & a_irregOnly,
                                 const Box                                 & a_validBox,
                                 const Box                                 & a_srcDomain,
                                 const Point                               & a_srcGhost,
                                 const Point                               & a_dstGhost,
                                 const XFacData                            & a_xfacdata,
                                 const YFacData                            & a_yfacdata,
                                 const ZFacData                            & a_zfacdata,
                                 const EBGraph                             & a_graph,
                                 const double                              & a_dx);

    inline static bool
    getDomainExtrapolationStencil(vector<EBIndex<fcent > >                  & a_dstFaces,
                                  vector<LocalStencil<fcent, double> >      & a_stencil,
                                  Stencil<double>                           & a_regStencil,
                                  Box                                       & a_regApplyBox,
                                  bool                                      & a_irregOnly,
                                  const Box                                 & a_validBox,
                                  const Box                                 & a_srcDomain,
                                  const Point                               & a_srcGhost,
                                  const Point                               & a_dstGhost,
                                  const XFacData                            & a_xfacdata,
                                  const YFacData                            & a_yfacdata,
                                  const ZFacData                            & a_zfacdata,
                                  const EBGraph                             & a_graph,
                                  const double                              & a_dx);
    

  };


  ///
  /**
     Stuff to go from cells to faces
   */
  template <CENTERING fcent, int order>
  class EBStencilArchive<CELL, fcent, order, double>
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

    EBStencilArchive() {}
    ~EBStencilArchive() {}
    

    ///   returns true if there are domain boundary conditions
    inline static bool
    getStencil(vector<EBIndex<fcent> >                    & a_dstVoFs,
               vector<LocalStencil<CELL, double> >       & a_stencil,
               Stencil<double>                           & a_regStencil,
               vector<Stencil<double> >                    a_regStencilBCS[2*DIM],
               vector<Box>                                 a_BCApplyBoxes[2*DIM],       
               Box                                       & a_regApplyBox,
               bool                                      & a_bcOnly,
               bool                                      & a_irregOnly,
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
               double                                      a_dx,
               bool                                        a_periodiconly, 
               Point                                       a_dstGrown);               


    /// for averaging cell centered stuff to face centers
    inline static bool
    getCellToFaceAveStencil(vector<EBIndex<fcent> >                   & a_dstVoFs,
                            vector<LocalStencil<CELL, double> >       & a_stencil,
                            Stencil<double>                           & a_regStencilInterior,
                            vector<Stencil<double> >                    a_regStencilBCS[2*DIM],
                            vector<Box>                                 a_BCApplyBoxes[2*DIM],       
                            Box                                       & a_regApplyBox,
                            string                                      a_domainBCName[2*DIM],
                            const Box                                 & a_srcValid,
                            const Box                                 & a_srcDomain,
                            const Point                               & a_srcGhost,
                            const Point                               & a_dstGhost,
                            const EBGraph                             & a_dstGraph,
                            double                                      a_dx,
                            Point                                       a_dstGrown);               


    ///
    inline static bool
    getCellToFaceStencil(vector<EBIndex<fcent> >                   & a_dstFaces,
                         vector<LocalStencil<CELL, double> >       & a_stencil,
                         Stencil<double>                           & a_regStencil,
                         Box                                       & a_regApplyBox,
                         const Box                                 & a_validBox,
                         const Box                                 & a_srcDomain,
                         const Point                               & a_srcGhost,
                         const Point                               & a_dstGhost,
                         const EBGraph                             & a_graph,
                         double                                      a_dx,
                         Point                                       a_dstGrown,
                         Side::LoHiSide                              a_side)    ;           

    ///
    inline static bool
    getMACGradientStencil(vector<EBIndex<fcent> >                   & a_dstFaces,
                          vector<LocalStencil<CELL, double> >       & a_stencil,
                          Stencil<double>                           & a_regStencil,
                          vector<Stencil<double> >                    a_regStencilBCS[2*DIM],
                          vector<Box>                                 a_BCApplyBoxes[2*DIM],       
                          Box                                       & a_regApplyBox,
                          const Box                                 & a_validBox,
                          const Box                                 & a_srcDomain,
                          const Point                               & a_srcGhost,
                          const Point                               & a_dstGhost,
                          const EBGraph                             & a_graph,
                          double                                      a_dx,
                          Point                                       a_dstGrown);
  };


  ///
  /**
     Stuff to go from faces to cells
   */
  template <CENTERING fcent, int order>
  class EBStencilArchive<fcent, CELL, order, double>
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

    EBStencilArchive() {}
    ~EBStencilArchive() {}
    

    ///   returns true if there are domain boundary conditions
    inline static bool
    getStencil(vector<EBIndex<CELL> >                     & a_dstVoFs,
               vector<LocalStencil<fcent, double> >       & a_stencil,
               Stencil<double>                           & a_regStencil,
               vector<Stencil<double> >                    a_regStencilBCS[2*DIM],
               vector<Box>                                 a_BCApplyBoxes[2*DIM],       
               Box                                       & a_regApplyBox,
               bool                                      & a_bcOnly,
               bool                                      & a_irregOnly,
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
               double                                      a_dx,
               bool                                        a_periodiconly, 
               Point                                       a_dstGrown);               


    inline static void
    getIrregLocations(vector<EBIndex<CELL> >        & a_dstVoFs,                    
                      const Stencil<double>         & a_regStencil,         
                      const Box                     & a_boundingBox,
                      const EBGraph                 & a_dstGraph,
                      const int                     & a_facedir);
    
    /// for averaging cell centered stuff to face centers
    inline static bool
    getDivergenceStencil(vector<EBIndex<CELL > >                   & a_dstVoFs,
                         vector<LocalStencil<fcent, double> >      & a_stencil,
                         Stencil<double>                           & a_regStencil,
                         Box                                       & a_regApplyBox,
                         const Box                                 & a_validBox,
                         const Box                                 & a_srcDomain,
                         const Point                               & a_srcGhost,
                         const Point                               & a_dstGhost,
                         const VoluData                            & a_voludata,
                         const EBFaData                            & a_ebfadata,
                         const XFacData                            & a_xfacdata,
                         const YFacData                            & a_yfacdata,
                         const ZFacData                            & a_zfacdata,
                         const EBGraph                             & a_graph,
                         double                                      a_dx);


    /// for averaging cell centered stuff to face centers
    inline static bool
    getFaceToCellAveStencil(vector<EBIndex<CELL > >                   & a_dstVoFs,
                            vector<LocalStencil<fcent, double> >      & a_stencil,
                            Stencil<double>                           & a_regStencil,
                            Box                                       & a_regApplyBox,
                            const Box                                 & a_validBox,
                            const Box                                 & a_srcDomain,
                            const Point                               & a_srcGhost,
                            const Point                               & a_dstGhost,
                            const VoluData                            & a_voludata,
                            const EBFaData                            & a_ebfadata,
                            const XFacData                            & a_xfacdata,
                            const YFacData                            & a_yfacdata,
                            const ZFacData                            & a_zfacdata,
                            const EBGraph                             & a_graph,
                            double                                      a_dx);

  };
  ///
  /**
   */
  template <int order>
  class EBStencilArchive<CELL, CELL, order, double>
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

    EBStencilArchive() {}
    ~EBStencilArchive() {}
    
    ///for prolongation and multicolor relaxation
    inline static Point getColor(unsigned long icolor)
    {
      Point retval;
      // Color offsets are grouped into "red"=even number of nonzeros (first 2^(DIM-1)) 
      // and "black= odd number of nonzeros (the rest).
#if DIM==2
      if(icolor == 0)
      {
        retval = Point::Zeros();//(0,0)
      }
      else if(icolor == 1)
      {
        retval = Point::Ones();//(1,1)
      }
      else if(icolor == 2)
      {
        retval = Point::Zeros() + Point::Basis(1);//(0,1)
      }
      else if(icolor == 3)
      {
        retval = Point::Zeros() + Point::Basis(0);//(1,0)
      }
      else
      {
        PROTO_ASSERT(false, "color not found");
      }
#elif DIM==3
      if(icolor == 0)
      {
        retval = Point::Zeros();//(0,0,0)
      }
      else if(icolor == 1)
      {
        retval = Point::Zeros() + Point::Basis(0) + Point::Basis(1);//(1,1,0)
      }
      else if(icolor == 2)
      {
        retval = Point::Zeros() + Point::Basis(1) + Point::Basis(2);//(0,1,1)
      }
      else if(icolor == 3)
      {
        retval = Point::Zeros() + Point::Basis(0) + Point::Basis(2);//(1,0,1)
      }
      else if(icolor == 4)
      {
        retval = Point::Zeros() + Point::Basis(1);//(0,1,0)
      }
      else if(icolor == 5)
      {
        retval = Point::Zeros() + Point::Basis(0);//(1,0,0)
      }
      else if(icolor == 6)
      {
        retval = Point::Zeros() + Point::Basis(2);//(0,0,1)
      }
      else if(icolor == 7)
      {
        retval = Point::Ones();//(1,1,1)
      }
      else
      {
        PROTO_ASSERT(false, "color not found");
      }
#else
      compiler_error_this_code_is_only_written_for_dim_2_or_3();
#endif
      return retval;
    }

    ///   returns true if there are domain boundary conditions
    inline static bool
    getStencil(vector<EBIndex<CELL> >                    & a_dstVoFs,
               vector<LocalStencil<CELL, double> >       & a_stencil,
               Stencil<double>                           & a_regStencilInterior,
               vector<Stencil<double> >                    a_regStencilBCS[2*DIM],
               vector<Box>                                 a_BCApplyBoxes[2*DIM],       
               Box                                       & a_regApplyBox,
               bool                                      & a_bcOnly,
               bool                                      & a_irregOnly,
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
               double                                      a_dx,
               bool                                        a_periodiconly,
               Point                                       a_dstGrown);               

    inline static bool
    getDomainBCStencil(vector<Stencil<double> >                    a_regStencilBC[2*DIM],       
                       vector<Box>                                 a_BCApplyBoxes[2*DIM],
                       bool                                        a_bcOnly,
                       string                                      a_domainBCName[2*DIM],
                       const  Box                                & a_validBox,
                       const Box                                 & a_domain,
                       const Point                               & a_srcGhost);                   
    ///stencil to do noramlize by volume fraction 
    inline static bool
    getNormalizationStencil(vector<EBIndex<CELL> >                    & a_dstVoFs,                    
                            vector<LocalStencil<CELL, double> >       & a_stencil,                    
                            Stencil<double>                           & a_regStencil,         
                            Box                                       & a_regApplyBox,                
                            const Box                                 & a_validBox,                   
                            const Box                                 & a_domain,                     
                            const Point                               & a_srcGhost,                   
                            const Point                               & a_dstGhost,                   
                            const EBGraph                             & a_graph,                      
                            const VoluData                            & a_voludata,                   
                            unsigned long                             & a_radius,
                            double                                      a_dx);

    ///Second order Poisson operator (interpolation to face centroids and all that)
    inline static bool
    getSecondOrderPoissonStencil(vector<EBIndex<CELL> >                    & a_dstVoFs,                    
                                 vector<LocalStencil<CELL, double> >       & a_stencil,                    
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
                                 double                                      a_dx,                         
                                 bool                                        a_periodiconly,
                                 bool  a_forceAllNeumann);               

    ///WLS stencil
    inline static bool
    getWLSStencil(vector<EBIndex<CELL> >                    & a_dstVoFs,
                  vector<LocalStencil<CELL, double> >       & a_stencil,
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
    
    ///Restriction
    /**
       This  produces phi_coarse = Averate(phi_fine).
       If kappa weighting is true, this does volume weighted averaging.
       If kappa weighting is false, this does straight averaging assuming that the inputs
       are already weighted with the volume fraction (as one does with multigrid).
       The source objects are on the finer level, the destination are on the coarser level.
       Refinement ratio of 2 is assumed.    Moment data should be on the finer level.
       Should be used with incrementonly = false;
       The valid box should be on the coarse level
       Returns false (no bcs needed here)
       Dx is on the fine level
     **/
    inline static bool
    getRestrictionStencil(vector<EBIndex<CELL> >                    & a_dstVoFs,                    
                          vector<LocalStencil<CELL, double> >       & a_stencil,                    
                          Stencil<double>                           & a_regStencilInterior,         
                          Box                                       & a_regApplyBox,                
                          const Box                                 & a_srcValid,                   
                          const Box                                 & a_dstValid,                   
                          const Point                               & a_srcGhost,                   
                          const Point                               & a_dstGhost,                   
                          const EBGraph                             & a_srcGraph,                      
                          const EBGraph                             & a_dstGraph,                      
                          const VoluData                            & a_voludata,                   
                          double                                      a_dx,
                          bool                                        a_volumeWeighed)               ;


    ///Restriction
    /**
       This  produces phi_coarse = Averate(phi_fine).
       If kappa weighting is true, this does volume weighted averaging.
       If kappa weighting is false, this does straight averaging assuming that the inputs
       are already weighted with the volume fraction (as one does with multigrid).
       The source objects are on the finer level, the destination are on the coarser level.
       Refinement ratio of 2 is assumed.    Moment data should be on the finer level.
       Should be used with incrementonly = false;
       The valid box should be on the coarse level
       Returns false (no bcs needed here)
       Dx is on the fine level
     **/
    inline static bool
    getRedistributionStencil(vector<EBIndex<CELL> >                    & a_dstVoFs,                    
                             vector<LocalStencil<CELL, double> >       & a_stencil,                    
                             Stencil<double>                           & a_regStencil,
                             Box                                       & a_regApplyBox,                
                             const Box                                 & a_validBox,                   
                             const Box                                 & a_domain,                   
                             const Point                               & a_srcGhost,                   
                             const Point                               & a_dstGhost,                   
                             const EBGraph                             & a_graph,                      
                             const VoluData                            & a_voludata,
                             int  a_redistRad,
                             double a_dx);                   
    
    
    /// slopes
    /**
       This produces low sided (phi_i - phi_i-1) or high sided (phi_i+1 - phi_i)
       depending on whether slopeSign is Lo or Hi
     */
    inline static bool
    getSlopeStencil(vector<EBIndex<CELL> >                    & a_dstVoFs,                    
                    vector<LocalStencil<CELL, double> >       & a_stencil,                    
                    Stencil<double>                           & a_regStencil,
                    vector<Stencil<double> >                    a_regStencilBC[2*DIM],       
                    vector<Box>                                 a_BCApplyBoxes[2*DIM],       
                    Box                                       & a_regApplyBox,                
                    string                                      a_domainBCName[2*DIM],
                    const Box                                 & a_validBox,                   
                    const Box                                 & a_domain,                   
                    const Point                               & a_srcGhost,                   
                    const Point                               & a_dstGhost,                   
                    const EBGraph                             & a_graph,                      
                    const VoluData                            & a_voludata,
                    const Point                               & a_dstGrown,
                    unsigned int  a_slopeDir,
                    Side::LoHiSide a_slopeSign);


    ///Prolongation
    /**
       This  produces phi_fine = piecewise_constant(phi_coar).  Mostly used in multigrid.
       The color says which fine cell color (in the multicolored Gauss Seidel sense) will be changed.  
       Should be used with incrementonly = true;
       The source objects are on the coarser level, the destination are on the finer level.
       Refinement ratio of 2 is assumed.    Moment data is not used.
       The valid box should be on the fine level
       Returns false (no bcs needed here)

     **/
    inline static bool
    getProlongationStencil(vector<EBIndex<CELL> >                    & a_dstVoFs,                    
                           vector<LocalStencil<CELL, double> >       & a_stencil,                    
                           Stencil<double>                           & a_regStencilInterior,         
                           Box                                       & a_regApplyBox,                
                           const Box                                 & a_srcValid,                   
                           const Box                                 & a_dstValid,                   
                           const Point                               & a_srcGhost,                   
                           const Point                               & a_dstGhost,                   
                           const EBGraph                             & a_srcGraph,                      
                           const EBGraph                             & a_dstGraph,                      
                           unsigned long                               a_icolor) ;              

    //need to use a grown set of the irregular cells/faces so that regular and irregular data talk
    inline static void
    getIrregLocations(vector<EBIndex<CELL> >        & a_dstVoFs,                    
                      const Stencil<double>         & a_regStencil,         
                      const Box                     & a_boundingBox,
                      const EBGraph                 & a_dstGraph);
  };

  ///
  /**
     Stuff to go from yface to xfaces
   */
  template <int order>
  class EBStencilArchive<YFACE, XFACE, order, double>
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

    EBStencilArchive() {}
    ~EBStencilArchive() {}
    

    ///   returns true if there are domain boundary conditions
    inline static bool
    getStencil(vector<EBIndex<XFACE> >                    & a_dstVoFs,
               vector<LocalStencil<YFACE, double> >       & a_stencil,
               Stencil<double>                           & a_regStencil,
               vector<Stencil<double> >                    a_regStencilBCS[2*DIM],
               vector<Box>                                 a_BCApplyBoxes[2*DIM],       
               Box                                       & a_regApplyBox,
               bool                                      & a_bcOnly,
               bool                                      & a_irregOnly,
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
               double                                      a_dx,
               bool                                        a_periodiconly, 
               Point                                       a_dstGrown);

  };


  ///
  /**
     Stuff to go from xface to yfaces
   */
  template <int order>
  class EBStencilArchive<XFACE, YFACE, order, double>
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

    EBStencilArchive() {}
    ~EBStencilArchive() {}
    

    ///   returns true if there are domain boundary conditions
    inline static bool
    getStencil(vector<EBIndex<YFACE> >                    & a_dstVoFs,
               vector<LocalStencil<XFACE, double> >       & a_stencil,
               Stencil<double>                           & a_regStencil,
               vector<Stencil<double> >                    a_regStencilBCS[2*DIM],
               vector<Box>                                 a_BCApplyBoxes[2*DIM],       
               Box                                       & a_regApplyBox,
               bool                                      & a_bcOnly,
               bool                                      & a_irregOnly,
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
               double                                      a_dx,
               bool                                        a_periodiconly, 
               Point                                       a_dstGrown);

  };


  ///
  /**
     Stuff to go from xface to zfaces
   */
  template <int order>
  class EBStencilArchive<XFACE, ZFACE, order, double>
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

    EBStencilArchive() {}
    ~EBStencilArchive() {}
    

    ///   returns true if there are domain boundary conditions
    inline static bool
    getStencil(vector<EBIndex<ZFACE> >                    & a_dstVoFs,
               vector<LocalStencil<XFACE, double> >       & a_stencil,
               Stencil<double>                           & a_regStencil,
               vector<Stencil<double> >                    a_regStencilBCS[2*DIM],
               vector<Box>                                 a_BCApplyBoxes[2*DIM],       
               Box                                       & a_regApplyBox,
               bool                                      & a_bcOnly,
               bool                                      & a_irregOnly,
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
               double                                      a_dx,
               bool                                        a_periodiconly, 
               Point                                       a_dstGrown);
  };

    ///
  /**
     Stuff to go from yface to zfaces
   */
  template <int order>
  class EBStencilArchive<YFACE, ZFACE, order, double>
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

    EBStencilArchive() {}
    ~EBStencilArchive() {}
    

    ///   returns true if there are domain boundary conditions
    inline static bool
    getStencil(vector<EBIndex<ZFACE> >                    & a_dstVoFs,
               vector<LocalStencil<YFACE, double> >       & a_stencil,
               Stencil<double>                           & a_regStencil,
               vector<Stencil<double> >                    a_regStencilBCS[2*DIM],
               vector<Box>                                 a_BCApplyBoxes[2*DIM],       
               Box                                       & a_regApplyBox,
               bool                                      & a_bcOnly,
               bool                                      & a_irregOnly,
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
               double                                      a_dx,
               bool                                        a_periodiconly, 
               Point                                       a_dstGrown);
  };

  ///
  /**
     Stuff to go from zface to xfaces
   */
  template <int order>
  class EBStencilArchive<ZFACE, XFACE, order, double>
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

    EBStencilArchive() {}
    ~EBStencilArchive() {}
    

    ///   returns true if there are domain boundary conditions
    inline static bool
    getStencil(vector<EBIndex<XFACE> >                    & a_dstVoFs,
               vector<LocalStencil<ZFACE, double> >       & a_stencil,
               Stencil<double>                           & a_regStencil,
               vector<Stencil<double> >                    a_regStencilBCS[2*DIM],
               vector<Box>                                 a_BCApplyBoxes[2*DIM],       
               Box                                       & a_regApplyBox,
               bool                                      & a_bcOnly,
               bool                                      & a_irregOnly,
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
               double                                      a_dx,
               bool                                        a_periodiconly, 
               Point                                       a_dstGrown);
    
  };

  ///
  /**
     Stuff to go from zface to yfaces
   */
  template <int order>
  class EBStencilArchive<ZFACE, YFACE, order, double>
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

    EBStencilArchive() {}
    ~EBStencilArchive() {}
    

    ///   returns true if there are domain boundary conditions
    inline static bool
    getStencil(vector<EBIndex<YFACE> >                    & a_dstVoFs,
               vector<LocalStencil<ZFACE, double> >       & a_stencil,
               Stencil<double>                           & a_regStencil,
               vector<Stencil<double> >                    a_regStencilBCS[2*DIM],
               vector<Box>                                 a_BCApplyBoxes[2*DIM],       
               Box                                       & a_regApplyBox,
               bool                                      & a_bcOnly,
               bool                                      & a_irregOnly,
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
               double                                      a_dx,
               bool                                        a_periodiconly, 
               Point                                       a_dstGrown);
    
  };

#if 1

  ///
  /**
     Stuff to go from CELL to BOUNDARYs
   */
  template <int order>
  class EBStencilArchive<CELL, BOUNDARY, order, double>
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


    EBStencilArchive() {}
    ~EBStencilArchive() {}
    

    ///   returns true if there are domain boundary conditions
    inline static bool
    getStencil(vector<EBIndex<BOUNDARY> >                & a_dstVoFs,
               vector<LocalStencil<CELL , double> >      & a_stencil,
               Stencil<double>                           & a_regStencil,
               vector<Stencil<double> >                    a_regStencilBCS[2*DIM],
               vector<Box>                                 a_BCApplyBoxes[2*DIM],       
               Box                                       & a_regApplyBox,
               bool                                      & a_bcOnly,
               bool                                      & a_irregOnly,
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
               double                                      a_dx,
               bool                                        a_periodiconly, 
               Point                                       a_dstGrown);
  };

  ///
  /**
     Stuff to go from BOUNDARY to CELLs
   */
  template <int order>
  class EBStencilArchive<BOUNDARY, CELL, order, double>
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

    EBStencilArchive() {}
    ~EBStencilArchive() {}
    

    ///   returns true if there are domain boundary conditions
    inline static bool
    getStencil(vector<EBIndex<CELL> >                    & a_dstVoFs,
               vector<LocalStencil<BOUNDARY, double> >   & a_stencil,
               Stencil<double>                           & a_regStencil,
               vector<Stencil<double> >                    a_regStencilBCS[2*DIM],
               vector<Box>                                 a_BCApplyBoxes[2*DIM],       
               Box                                       & a_regApplyBox,
               bool                                      & a_bcOnly,
               bool                                      & a_irregOnly,
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
               double                                      a_dx,
               bool                                        a_periodiconly, 
               Point                                       a_dstGrown);
  };
#endif
}

#include "Proto_EBStencilArchiveImplem.H"

#endif
  
