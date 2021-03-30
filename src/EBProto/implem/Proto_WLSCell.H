
#ifndef _Proto_WLSCell_H_
#define _Proto_WLSCell_H_

#include "Proto_Normalizor.H"
#include "Chombo_ParmParse.H"
#include "Chombo_ProtoInterface.H"
#include "Chombo_LAPACKMatrix.H"
#include "EBProto.H"

namespace Proto
{
  template<int order>
  class Neighbors
  {
  public:
    typedef IndexedMoments<DIM  , order> IndMomDIM;
    typedef IndexedMoments<DIM-1, order> IndMomSDMinOne;
    typedef MomentIterator<DIM  , order> MomItDIM;
    typedef MomentIterator<DIM-1, order> MomItSDMinOne;

    Neighbors()
    {
      m_volNeighs.resize(0);
      m_volMom.resize(0);
      m_facMom.resize(0);
      m_ebfMom.resize(0);
      m_volLoc.resize(0);
      m_facLoc.resize(0);
      m_ebfLoc.resize(0);
    };

    Neighbors(Vector<EBIndex<CELL> > & volNeighs,
              Vector<IndMomDIM>      & volMom,
              Vector<IndMomSDMinOne> & facMom,
              Vector<IndMomDIM>      & ebfMom,
              Vector<RealVect>       & volLoc,
              Vector<RealVect>       & facLoc,
              Vector<RealVect>       & ebfLoc)
    {
      m_volNeighs = volNeighs;
      m_volMom = volMom;
      m_facMom = facMom;
      m_ebfMom = ebfMom;
      m_volLoc = volLoc;
      m_facLoc = facLoc;
      m_ebfLoc = ebfLoc;
    };


    Vector<EBIndex<CELL> >  m_volNeighs;
    Vector<IndMomDIM>       m_volMom;
    Vector<IndMomSDMinOne>  m_facMom;
    Vector<IndMomDIM>       m_ebfMom;
    Vector<RealVect>        m_volLoc;
    Vector<RealVect>        m_facLoc;
    Vector<RealVect>        m_ebfLoc;
  
  };

  ///stuff in common for all WLS flux calculations
  template<CENTERING cent, int order>
  class BaseWLSFlux
  {
  public:
    typedef IndexedMoments<DIM  , order> IndMomDIM;
    typedef IndexedMoments<DIM-1, order> IndMomSDMinOne;
    typedef MomentIterator<DIM  , order> MomItDIM;
    typedef MomentIterator<DIM-1, order> MomItSDMinOne;
    typedef HostIrregData<CELL    ,  IndMomDIM , 1>  VoluData;
    typedef HostIrregData<BOUNDARY,  IndMomDIM , 1>  EBFaData;
    typedef HostIrregData<XFACE, IndMomSDMinOne, 1>  XFacData;
    typedef HostIrregData<YFACE, IndMomSDMinOne, 1>  YFacData;
    typedef HostIrregData<ZFACE, IndMomSDMinOne, 1>  ZFacData;
    
    vector<EBIndex<CELL> >
    getFlapVoFs(const EBGraph       & a_graph,  
	              EBIndex<CELL> & a_start) const;

    vector<EBIndex<CELL> >
    getVoFsNeighbors(const EBIndex<cent> & a_face,
		     unsigned int          a_radius,
		     Vector<int>         & a_weightDisVol,
		     const EBGraph       & a_graph) const;

    void
    getMomentsInfoFromNeighbors(vector<EBIndex<CELL> > & a_vols,
		                const EBGraph          & a_graph,
                                const VoluData         & a_voludata,
                                const EBFaData         & a_ebfadata,
                                const XFacData         & a_xfacdata,
                                const YFacData         & a_yfacdata,
                                const ZFacData         & a_zfacdata,
		                Vector<IndMomDIM>      & a_volmom, 
				Vector<IndMomDIM>      & a_ebfmom,
                                Vector<IndMomSDMinOne> & a_facmom,
                                Vector<RealVect>       & a_volneighLoc,
                                Vector<RealVect>       & a_ebfneighLoc,
                                Vector<RealVect>       & a_facneighLoc,
				Vector<int>            & a_weightDisVol,
                                Vector<int>            & a_weightDisEB,
                                Vector<int>            & a_weightDisFac,
                                Vector<int>            & a_facneighDir,
				Vector<Side::LoHiSide> & a_facneighSid,
				Vector<RealVect>       & a_ebfneighNorm,
				double                   a_dx) const;

    void
    getShiftedMomentMatrix(Chombo4::LAPACKMatrix  & a_Mvol,
		           RealVect               & a_faceLoc,
			   double                   a_dx,
			   Vector<IndMomDIM>      & a_volmom,
                           Vector<IndMomDIM>      & a_ebfmom,
                           Vector<IndMomSDMinOne> & a_facmom,
                           Vector<RealVect>       & a_volneighLoc,
                           Vector<RealVect>       & a_ebfneighLoc,
                           Vector<RealVect>       & a_facneighLoc,
			   Vector<int>            & a_facneighDir,
			   Vector<Side::LoHiSide> & a_facneighSid,
			   Vector<RealVect>       & a_ebfneighNorm,
			   string                   a_ebbcName,
                           string                   a_dombcName[2*DIM]) const;
		           
    void
    getFaceCenteredStencil(LocalStencil<CELL, double> & a_sten,
                           const EBIndex<cent>        & a_face,
                           const EBGraph              & a_graph,
                           const VoluData             & a_voludata,
                           const EBFaData             & a_ebfadata,
                           const XFacData             & a_xfacdata,
                           const YFacData             & a_yfacdata,
                           const ZFacData             & a_zfacdata,
			   string                       a_ebbcName,
                           string                       a_dombcName[2*DIM],
                           const int                  & a_dir,
			   IndMomSDMinOne             & a_faceMom,
                           double                       a_dx) const;

    void
    getFaceFluxMatrix(Chombo4::LAPACKMatrix    & a_Qshift,
                      const EBIndex<cent>      & a_face,
		      RealVect                 & a_faceLoc,
	              IndMomSDMinOne           & a_faceMom,
                      double                     a_dx) const;

    void
    getEBFluxMatrix(Chombo4::LAPACKMatrix    & a_Qshift,
                    const EBIndex<cent>      & a_face,
		    RealVect                 & a_faceLoc,
		    IndMomDIM                & a_ebmom,
		    RealVect                 & a_normal,
		    vector<IndMomDIM>        & a_ebnorm,
                    double                     a_dx) const;

    void
    getWeightsMatrixFace(Chombo4::LAPACKMatrix  & a_weights,
		         Vector<int>            & a_weightDisVol,
                         Vector<int>            & a_weightDisFac,
                         Vector<int>            & a_weightDisEB) const;

    void
    getStencilWeights(Chombo4::LAPACKMatrix  & a_phiCoef,
		      Chombo4::LAPACKMatrix  & a_Qshift, 
		      Chombo4::LAPACKMatrix  & a_Mvol, 
		      Chombo4::LAPACKMatrix  & a_weights, 
                      double                   a_dx) const;

    void
    getEBFluxStencil(LocalStencil<CELL, double>    &  a_facesten,
                     double                        &  a_weight, //for inhomogeneous bcs
                     const EBIndex<cent>           &  a_face,
                     const EBGraph                 & a_graph,
                     const VoluData                & a_voludata,
                     const EBFaData                & a_ebfadata,
                     const XFacData                & a_xfacdata,
                     const YFacData                & a_yfacdata,
                     const ZFacData                & a_zfacdata,
		     vector<IndMomDIM>             & a_ebnorm,
		     string                          a_ebbcName,
                     string                          a_dombcName[2*DIM],
		     IndMomDIM                     & a_ebmom,
		     RealVect                      & a_normal,
                     const double                  & a_dx) const;

    virtual ~BaseWLSFlux()
    {
    }

    BaseWLSFlux()
    {
    }
    
  };
  
  ///WLS Poisson flux
  template<CENTERING cent, int order>
  class PoissonWLS: public BaseWLSFlux<cent, order>
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

    virtual ~PoissonWLS()
    {
    }

    PoissonWLS()
    {
    }
  };

  ///some basic infrastructure to support the WLS stuff
  template<int order>  //this is so we can send in the geometry informaion sanely
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
    inline static void
    getWLSDivFStencil(LocalStencil< CELL,  double>             & a_vofStencil,
                           const EBIndex<CELL>                 & a_vof,
                           const EBGraph                       & a_graph,
                           const VoluData                      & a_voludata,
                           const EBFaData                      & a_ebfadata,
                           const XFacData                      & a_xfacdata,
                           const YFacData                      & a_yfacdata,
                           const ZFacData                      & a_zfacdata,
 			   const EBNormalData                  & a_ebnormxdata,
                           const EBNormalData                  & a_ebnormydata,
#if DIM==3
                           const EBNormalData                  & a_ebnormzdata,
#endif
                           const BaseWLSFlux<XFACE   , order>  & a_xfaceflux,
                           const BaseWLSFlux<YFACE   , order>  & a_yfaceflux,
                           const BaseWLSFlux<ZFACE   , order>  & a_zfaceflux,
                           const BaseWLSFlux<BOUNDARY, order>  & a_ebflux,
                           string                                a_ebbcName,
                           string                                a_dombcName[2*DIM],
                           double                                a_dx);

    virtual ~WLSStencil()
    {
    }

    WLSStencil()
    {
    }
  };

  template<int order> //order is so we can send in geometry info
  class DirichletBC
  {
  public:
    typedef IndexedMoments<DIM  , order> IndMomDIM;
    typedef IndexedMoments<DIM-1, order> IndMomSDMinOne;
    typedef MomentIterator<DIM  , order> MomItDIM;

    DirichletBC()
    {
    }

    virtual ~DirichletBC()
    {
    }

    inline static void
    getFaceMoments(IndMomDIM      & a_fac2volmom,
                   int              a_faceDir,
                   IndMomSDMinOne & a_facmom);

  };

  template<int order> //order is so we can send in geometry info
  class NeumannBC
  {
  public:
    typedef IndexedMoments<DIM  , order> IndMomDIM;
    typedef IndexedMoments<DIM-1, order> IndMomSDMinOne;
    typedef MomentIterator<DIM  , order> MomItDIM;

    NeumannBC()
    {
    }

    virtual ~NeumannBC()
    {
    }

    inline static void
    getFaceMoments(IndMomDIM      & a_fac2volmom,
                   int              a_faceDir,
                   IndMomSDMinOne & a_facmom);

    inline static void
    getEBMoments(IndMomDIM      & a_ebfmom,
		 RealVect       & a_ebfNorm);

  };

}
#include "Proto_WLSCellImplem.H"

#endif
  
