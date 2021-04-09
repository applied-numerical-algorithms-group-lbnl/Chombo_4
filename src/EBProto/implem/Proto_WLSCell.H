
#ifndef _Proto_WLSCell_H_
#define _Proto_WLSCell_H_

#include "Proto_Normalizor.H"
#include "Chombo_ParmParse.H"
#include "Chombo_ProtoInterface.H"
#include "Chombo_LAPACKMatrix.H"
#include "EBProto.H"

namespace Proto
{
 
  ///class where phi is a monomial in x, y, z
  template<int order> //order is so we can send in geometry info
  class MonomialPhi
  {
  public:
    typedef IndexTM<int,DIM> IndMomDIM;
    ///
    /** monomial is a funciton of x, y, z
        a_power = (p1, p2, p3) 
        a_point = (x0, y0, z0)
        the monomial is (x-x0)^p1*(y-y0)^p2*(z-z0)^p3 */
    MonomialPhi();
  
    ///
    static inline Real getTaylorCoefficient(const RealVect   & a_x0,
                                            const IndMomDIM  & a_p
                                           ) 
    {
      IntVect m_power;
      RealVect m_point;
      Vector<int> powerVec;
      Vector<Real> pointVec;

      ParmParse pp;
      pp.getarr("phi_power",powerVec,0,DIM);
      pp.getarr("phi_point",pointVec,0,DIM);
      for(int idir = 0; idir < DIM; idir++)
      {
        m_point[idir] = pointVec[idir];
        m_power[idir] = powerVec[idir];
      }

      Real retval = 1.0;
      for(int idir = 0; idir < DIM; idir++)
      {
        if(a_p[idir] > m_power[idir])
        {
          return 0.;
        }
        else
        {
          Real prefix = 1.0;
          for(int ipre = 0; ipre < a_p[idir]; ipre++)
          {
            prefix *= (m_power[idir] - ipre);
          }
          retval *= prefix*pow(a_x0[idir]-m_point[idir],m_power[idir]-a_p[idir])/factorial(a_p[idir]);
        }
      }
    
        return retval;
    }

  
    virtual ~MonomialPhi()
    {
    }
  
  };


	 

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
    typedef HostIrregData<CELL    ,  IndMomDIM , 1>  EBNormalData;
    
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
				const EBNormalData     & a_ebnormxdata,
                                const EBNormalData     & a_ebnormydata,
                                const EBNormalData     & a_ebnormzdata,
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
				Vector<Vector<IndMomDIM> > & a_ebfneighNorm,
				Vector<double>         & a_volneighFrac,
                                Vector<double>         & a_ebfneighFrac,
                                Vector<double>         & a_facneighFrac,
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
			   Vector<Vector<IndMomDIM> > & a_ebfneighNorm,
			   Vector<double>         & a_volneighFrac,
                           Vector<double>         & a_ebfneighFrac,
                           Vector<double>         & a_facneighFrac,
			   string                   a_ebbcName,
                           string                   a_dombcName[2*DIM]) const;
		           
    void
    getFaceCenteredStencil(LocalStencil<CELL, double> & a_sten,
                           Real                       &  a_inhomogTerm, //for inhomogeneous bcs
                           const EBIndex<cent>        & a_face,
                           const EBGraph              & a_graph,
                           const VoluData             & a_voludata,
                           const EBFaData             & a_ebfadata,
                           const XFacData             & a_xfacdata,
                           const YFacData             & a_yfacdata,
                           const ZFacData             & a_zfacdata,
                           const EBNormalData         & a_ebnormxdata,
                           const EBNormalData         & a_ebnormydata,
                           const EBNormalData         & a_ebnormzdata,
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
		    vector<IndMomDIM>        & a_ebnorm,
		    double                   & a_bndryArea,
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
                     Real                          &  a_inhomogTerm, //for inhomogeneous bcs
                     const EBIndex<cent>           &  a_face,
                     const EBGraph                 & a_graph,
                     const VoluData                & a_voludata,
                     const EBFaData                & a_ebfadata,
                     const XFacData                & a_xfacdata,
                     const YFacData                & a_yfacdata,
                     const ZFacData                & a_zfacdata,
		     const EBNormalData            & a_ebnormxdata,
                     const EBNormalData            & a_ebnormydata,
                     const EBNormalData            & a_ebnormzdata,
		     vector<IndMomDIM>             & a_ebnorm,
		     string                          a_ebbcName,
                     string                          a_dombcName[2*DIM],
		     IndMomDIM                     & a_ebmom,
		     RealVect                      & a_normal,
		     double                        & a_bndryArea,
                     const double                  & a_dx) const;

    Real
    computeInhomogTerm(const Chombo4::LAPACKMatrix    & a_stencil,
                     const Vector<IndMomSDMinOne>   & a_facmom,
                     const Vector<IndMomDIM>        & a_ebfmom,
                     const int                        a_nvolneigh,
		     const Vector<RealVect>         & a_facneighLoc,
		     const Vector<int>              & a_facneighDir,
                     const Vector<Side::LoHiSide>   & a_facneighSid,
		     const Vector<RealVect>         & a_ebfneighLoc,
	             const Vector<Vector<IndMomDIM> >& a_ebfneighNorm,
		     string                           a_ebbcName,
                     string                           a_dombcName[2*DIM],
                     double                           a_dx) const;

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
                           const EBNormalData                  & a_ebnormzdata,
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
    getEBMoments(IndMomDIM         & a_ebfmom,
		 Vector<IndMomDIM> & a_ebfNorm);

  };

  template<int order> //order is so we can send in geometry info
  class DirichletExactSoln
  {
  public:
    typedef IndexedMoments<DIM  , order> IndMomDIM;
    typedef IndexedMoments<DIM-1, order> IndMomSDMinOne;
    typedef MomentIterator<DIM  , order> MomItDIM;

    DirichletExactSoln()
    {
    }

    virtual ~DirichletExactSoln()
    {
    }

    inline static Real getFaceIntegral(const IndMomSDMinOne  & a_faceMom,
		                       const RealVect        & a_faceLoc,
                                       int                     a_faceDir,
                                       const Real              a_dx);

    inline static Real getEBIntegral(const IndMomDIM  & a_ebfMom,
                                     const RealVect   & a_ebfLoc,
                                     const Real         a_dx); 
  
  };

  template<int order> //order is so we can send in geometry info
  class NeumannExactSoln
  {
  public:
    typedef IndexedMoments<DIM  , order> IndMomDIM;
    typedef IndexedMoments<DIM-1, order> IndMomSDMinOne;
    typedef MomentIterator<DIM  , order> MomItDIM;

    NeumannExactSoln()
    {
    }

    virtual ~NeumannExactSoln()
    {
    }

    inline static Real getFaceIntegral(const IndMomSDMinOne  & a_faceMom,
		                       const RealVect        & a_faceLoc,
                                       int                     a_faceDir,
		                       const Real              a_dx); 

    inline static Real getEBIntegral(const Vector<IndMomDIM> & a_ebfNormalMom,
                                     const RealVect          & a_ebfLoc,
                                     const Real                a_dx);
 
  };

}
#include "Proto_WLSCellImplem.H"

#endif
  
