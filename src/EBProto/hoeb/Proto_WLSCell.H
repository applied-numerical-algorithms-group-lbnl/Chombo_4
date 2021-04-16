
#ifndef _Proto_WLSCell_H_
#define _Proto_WLSCell_H_

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

  //Class that will contain all the information related to neighbors VoFs,
  //boundary faces and eb faces in a certain radius.
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
      m_volMom.resize(0);
      m_facMom.resize(0);
      m_ebfMom.resize(0);
      m_volLoc.resize(0);
      m_facLoc.resize(0);
      m_ebfLoc.resize(0);
      m_weightDisVol.resize(0);
      m_weightDisFac.resize(0);
      m_weightDisEB.resize(0);
      m_volneighFrac.resize(0);
      m_ebfneighFrac.resize(0);
      m_facneighFrac.resize(0);
      m_ebfneighNorm.resize(0);
      m_facneighDir.resize(0);
      m_facneighSid.resize(0);
    };

    Neighbors(Vector<IndMomDIM>      & a_volMom,
              Vector<IndMomSDMinOne> & a_facMom,
              Vector<IndMomDIM>      & a_ebfMom,
              Vector<RealVect>       & a_volLoc,
              Vector<RealVect>       & a_facLoc,
              Vector<RealVect>       & a_ebfLoc,
	      Vector<double>         & a_volneighFrac,
              Vector<double>         & a_ebfneighFrac,
              Vector<double>         & a_facneighFrac,
              Vector<Vector<IndMomDIM> > & a_ebfneighNorm,
              Vector<int>            & a_facneighDir,
              Vector<Side::LoHiSide> & a_facneighSid)
    {
      m_volMom = a_volMom;
      m_facMom = a_facMom;
      m_ebfMom = a_ebfMom;
      m_volLoc = a_volLoc;
      m_facLoc = a_facLoc;
      m_ebfLoc = a_ebfLoc;
      m_volneighFrac = a_volneighFrac;
      m_ebfneighFrac = a_ebfneighFrac;
      m_facneighFrac = a_facneighFrac;
      m_ebfneighNorm = a_ebfneighNorm;
      m_facneighDir = a_facneighDir;
      m_facneighSid = a_facneighSid;
    };

    inline void setWeightDisVol (Vector<int> & a_weightDisVol)
    {
      m_weightDisVol = a_weightDisVol;
    };

    inline void setMomentsInfo (Vector<IndMomDIM>      & a_volMom,
		                Vector<IndMomSDMinOne> & a_facMom,
				Vector<IndMomDIM>      & a_ebfMom,
				Vector<Vector<IndMomDIM> > & a_ebfneighNorm)
    {
      m_volMom = a_volMom;
      m_facMom = a_facMom;
      m_ebfMom = a_ebfMom;
      m_ebfneighNorm.resize(a_ebfneighNorm.size());
      for (int i = 0; i < a_ebfneighNorm.size(); ++i)
      {
        m_ebfneighNorm[i] = a_ebfneighNorm[i];
      }
    };

    inline Vector<Vector<IndMomDIM> > getEBNormMoments ()
    { 
      Vector<Vector<IndMomDIM> > ebfneighNorm;
      ebfneighNorm.resize(m_ebfneighNorm.size());
      for (int i = 0; i < m_ebfneighNorm.size(); ++i)
      {
        ebfneighNorm[i] = m_ebfneighNorm[i];
      }
      return ebfneighNorm;
    };

    inline void setLocInfo (Vector<RealVect> & a_volLoc,
                            Vector<RealVect> & a_facLoc,
                            Vector<RealVect> & a_ebfLoc)
    {
      m_volLoc = a_volLoc;
      m_facLoc = a_facLoc;
      m_ebfLoc = a_ebfLoc;
    };

    inline void setWeightDisFacAndEB (Vector<int> & a_weightDisFac,
		                      Vector<int> & a_weightDisEB)
    {
      m_weightDisFac = a_weightDisFac;
      m_weightDisEB = a_weightDisEB;
    };

    inline void setFracInfo (Vector<double> & a_volneighFrac,
		             Vector<double> & a_ebfneighFrac,
			     Vector<double> & a_facneighFrac)
    {
      m_volneighFrac = a_volneighFrac;
      m_ebfneighFrac = a_ebfneighFrac;
      m_facneighFrac = a_facneighFrac;
    };

    inline void setAdditionalFaceInfo (Vector<int>            & a_facneighDir,
                                       Vector<Side::LoHiSide> & a_facneighSid)
    {
      m_facneighDir = a_facneighDir;
      m_facneighSid = a_facneighSid;
    };

    //Volume moments
    Vector<IndMomDIM>       m_volMom;
    //Face moments
    Vector<IndMomSDMinOne>  m_facMom;
    //EB Faces moments
    Vector<IndMomDIM>       m_ebfMom;
    //EB Normal moments
    Vector<Vector<IndMomDIM> > m_ebfneighNorm;
    //Cell location, will be the cell center
    Vector<RealVect>        m_volLoc;
    //Boundary Face location, will be the face center 
    Vector<RealVect>        m_facLoc;
    //EB Face location, will be the cell center containing the EB Face
    Vector<RealVect>        m_ebfLoc;
    //Contain the manhattan distance from the face 
    //or EB face being evaluated to the cell
    Vector<int>             m_weightDisVol;
    //Contain the manhattan distance from the face 
    //or EB face being evaluated to the boundary face
    Vector<int>             m_weightDisFac;
    //Contain the manhattan distance from the face
    //or EB face being evaluated to the EB Face
    Vector<int>             m_weightDisEB;
    //Volume fraction
    Vector<double>          m_volneighFrac;
    //EB Face area fraction
    Vector<double>          m_ebfneighFrac;
    //Face area fraction
    Vector<double>          m_facneighFrac;
    //Face direction
    Vector<int>             m_facneighDir;
    //Face Side -- lo or hi
    Vector<Side::LoHiSide>  m_facneighSid;
  
  };

  //Main function calls for all WLS flux calculations
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
   
    //Get neighbors VoFs through the faces of a specific cell.
    vector<EBIndex<CELL> >
    getFlapVoFs(const EBGraph       & a_graph,  
	              EBIndex<CELL> & a_start) const;

    //Get neighbors VoFs of a given face within a specific radius.
    //Radius starts at one for the cells sharing the face, and 
    //increases by one for the neighbor cells sharing a face with the original
    //cells of radius one.
    Vector<EBIndex<CELL> >
    getVoFsNeighbors(const EBIndex<cent> & a_face,
		     const int             a_radius,
		     Neighbors<order>    & a_neigh,
		     const EBGraph       & a_graph) const;

    //Gets all the moments and geometric info for the neighbors cells, boundary faces
    //and EB faces. Stores all this info in the Neighbors object.
    void
    getMomentsInfoFromNeighbors(Vector<EBIndex<CELL> > & a_vols,
		                const EBGraph          & a_graph,
                                const VoluData         & a_voludata,
                                const EBFaData         & a_ebfadata,
                                const XFacData         & a_xfacdata,
                                const YFacData         & a_yfacdata,
                                const ZFacData         & a_zfacdata,
				const EBNormalData     & a_ebnormxdata,
                                const EBNormalData     & a_ebnormydata,
                                const EBNormalData     & a_ebnormzdata,
				Neighbors<order>       & a_neigh,
				double                   a_dx) const;

    //Compute the Moments Matrix
    void
    getShiftedMomentMatrix(Chombo4::LAPACKMatrix  & a_Mvol,
		           RealVect               & a_faceLoc,
			   double                   a_dx,
		           Neighbors<order>       & a_neigh,
			   string                   a_ebbcName,
                           string                   a_dombcName[2*DIM]) const;
		           
    //Compute rhs matrix for a face
    void
    getFaceFluxMatrix(Chombo4::LAPACKMatrix    & a_Qshift,
                      const EBIndex<cent>      & a_face,
		      RealVect                 & a_faceLoc,
	              IndMomSDMinOne           & a_faceMom,
	              const double             & a_areaFrac,
                      double                     a_dx) const;

    //Compute rhs matrix for a EB face
    void
    getEBFluxMatrix(Chombo4::LAPACKMatrix    & a_Qshift,
                    const EBIndex<cent>      & a_face,
		    RealVect                 & a_faceLoc,
		    vector<IndMomDIM>        & a_ebnorm,
		    double                   & a_bndryArea,
                    double                     a_dx) const;

    //Compute weight matrix 
    void
    getWeightsMatrixFace(Chombo4::LAPACKMatrix  & a_weights,
		         Neighbors<order>       & a_neigh) const;

    //Compute stencil coefficients matrix 
    void
    getStencilWeights(Chombo4::LAPACKMatrix  & a_phiCoef,
		      Chombo4::LAPACKMatrix  & a_Qshift, 
		      Chombo4::LAPACKMatrix  & a_Mvol, 
		      Chombo4::LAPACKMatrix  & a_weights, 
                      double                   a_dx) const;

    //Compute inhomogeneous term
    Real
    computeInhomogTerm(const Chombo4::LAPACKMatrix     & a_stencil,
                       Neighbors<order>                & a_neigh,
                       const int                         a_nvolneigh,
                       string                            a_ebbcName,
                       string                            a_dombcName[2*DIM],
                       double                            a_dx) const;

    //Compute face centered flux stencil
    void
    getFaceCenteredStencil(LocalStencil<CELL, double> & a_sten,
                           Real                       & a_inhomogTerm, //for inhomogeneous bcs
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
                           const double               & a_areaFrac,
                           double                       a_dx) const;

    //Compute EB face stencil
    void
    getEBFluxStencil(LocalStencil<CELL, double>    & a_facesten,
                     Real                          & a_inhomogTerm, //for inhomogeneous bcs
                     const EBIndex<cent>           & a_face,
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

    virtual ~BaseWLSFlux()
    {
    }

    BaseWLSFlux(int          & a_radius,
		int          & a_powerE,
		bool           a_verbose=false)
    {
      m_radius = a_radius;
      m_powerE = a_powerE;
      m_verbose = a_verbose;
    }

    //Radius needed for neighbors list around a face (Manhattan style)
    int m_radius = 1;
    //Exponent for the distance between the VoFs and the face in the weights matrix
    int m_powerE = -3;
    //Boolean to enable print out statements for debugging
    bool m_verbose = false;
    
  };
  
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
    //Gets WLS flux divergence stencil
    inline static void
    getWLSDivFStencil(LocalStencil< CELL,  double>        & a_vofStencil,
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

  //Class to compute Dirichlet BCs at faces and EB Faces for the Moment Matrix.
  template<int order>
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

  //Class to compute Neumann BCs at faces and EB Faces for the Moment Matrix.
  template<int order> 
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

  //Class to compute Dirichlet BCs at faces and EB Faces. Required for inhomogeneous conditions.
  template<int order> 
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

    //getIntegralPhi
    inline static Real getFaceIntegral(const IndMomSDMinOne  & a_faceMom,
		                       const RealVect        & a_faceLoc,
                                       int                     a_faceDir,
                                       const Real              a_dx);

    //getIntegralPhiAlongEB
    inline static Real getEBIntegral(const IndMomDIM  & a_ebfMom,
                                     const RealVect   & a_ebfLoc,
                                     const Real         a_dx); 
  
  };

  //Class to compute Neumann BCs at faces and EB Faces. Required for inhomogeneous conditions.
  template<int order> 
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

    //getIntegralNormalDerivative
    inline static Real getFaceIntegral(const IndMomSDMinOne  & a_faceMom,
		                       const RealVect        & a_faceLoc,
                                       int                     a_faceDir,
		                       const Real              a_dx); 

    //getIntegralGradPhiDotNormal
    inline static Real getEBIntegral(const Vector<IndMomDIM> & a_ebfNormalMom,
                                     const RealVect          & a_ebfLoc,
                                     const Real                a_dx);
 
  };

}
#include "Proto_WLSCellImplem.H"

#endif
  
