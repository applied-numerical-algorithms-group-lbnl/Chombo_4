
#ifndef _Proto_WLSCell_H_
#define _Proto_WLSCell_H_

#include "Chombo_ProtoInterface.H"
#include "Chombo_LAPACKMatrix.H"
#include "EBProto.H"
#include "Proto_EBExactSolutions.H"

namespace Proto
{
  template<int order> 
  class BaseExactVelocity
  {
  public:
    typedef IndexedMoments<DIM  , order> IndMomDIM;
    typedef IndexedMoments<DIM-1, order> IndMomSDMinOne;
    typedef MomentIterator<DIM  , order> MomItDIM;

    bool locationBased() const
    {
      return m_locationBased;
    }
  
    ///most of our examples are location based
    BaseExactVelocity(bool a_locationBased = true)
    {
      m_locationBased = a_locationBased;
    }
  
    ///
    /**
     */
    virtual IndMomDIM getVelCompTaylorCoeffs(const RealVect   & a_x0,
                                             const int        & a_velcomp
                                            ) const
    {
      Chombo4::MayDay::Error("not implemented--only for location based vels");
      IndMomDIM retval;
      retval.setToZero();
      return retval;
    }
 
    virtual RealVect getVelocityAtPoint(const RealVect& a_x) const
    {
      Chombo4::MayDay::Error("not implemented--only implemented for derived classes");
      RealVect retval;
      return retval;
    }

    virtual void setDx(const Real a_dx)
    {
      m_dx = a_dx;
    }

    virtual Real getDx()
    {
      return m_dx;
    }

    virtual ~BaseExactVelocity()
    {
    }
  
    Real m_dx; 
  protected:
  
    bool m_locationBased;
  
  };
  
  ///
  /**
   */
  template<int order> 
  class ConstantVectorVel: public BaseExactVelocity<order>
  {
  public:
    typedef IndexedMoments<DIM  , order> IndMomDIM;
    typedef IndexedMoments<DIM-1, order> IndMomSDMinOne;
    typedef MomentIterator<DIM  , order> MomItDIM;
 
    ///  
    /**
     */
    ConstantVectorVel(const RealVect & a_value)
      :BaseExactVelocity<order>(true) // this is location based
    {
      m_value = a_value;
    }
  
    ///
    /**
     */
    virtual IndMomDIM getVelCompTaylorCoeffs(const RealVect   & a_x0,
                                             const int        & a_velcomp) const;
  
  
    virtual RealVect getVelocityAtPoint(const RealVect& a_x) const;
  
    virtual ~ConstantVectorVel()
    {
    }
  
  private:
    //weak construction is bad
    ConstantVectorVel();
  
    RealVect m_value;
  };

  ///
  /**
   */
  template<int order> 
  class PoiseuilleVectorVel: public BaseExactVelocity<order>
  {
  public:
    typedef IndexedMoments<DIM  , order> IndMomDIM;
    typedef IndexedMoments<DIM-1, order> IndMomSDMinOne;
    typedef MomentIterator<DIM  , order> MomItDIM;
 
    ///  
    /**
     */
    PoiseuilleVectorVel(const Real & a_centerValue,
                        const Real & a_wallValue,
                        const Real & a_radius,
                        const RealVect & a_center)
      :BaseExactVelocity<order>(true) // this is location based
    {
      m_centerValue = a_centerValue;
      m_wallValue = a_wallValue;
      m_center = a_center;
      m_radius = a_radius;
    }
  
    ///
    /**
     */
    virtual IndMomDIM getVelCompTaylorCoeffs(const RealVect   & a_x0,
                                             const int        & a_velcomp) const;
  
  
    virtual RealVect getVelocityAtPoint(const RealVect& a_x) const;
  
    virtual ~PoiseuilleVectorVel()
    {
    }
  
  private:
    //weak construction is bad
    PoiseuilleVectorVel();
  
    RealVect m_center;
    Real m_radius;
    Real m_centerValue;
    Real m_wallValue;
  };


  ///
  /**
   */
  template<int order> 
  class SolidBodyRotVel: public BaseExactVelocity<order>
  {
  public:
    typedef IndexedMoments<DIM  , order> IndMomDIM;
    typedef IndexedMoments<DIM-1, order> IndMomSDMinOne;
    typedef MomentIterator<DIM  , order> MomItDIM;
  
    ///  
    /**
     */
    SolidBodyRotVel(const Real     & a_omega,
                    const RealVect & a_center)
      :BaseExactVelocity<order>(true) // this is location based
    {
      m_omega  = a_omega;
      m_center = a_center;
    }
  
    ///
    /**
     */
    virtual IndMomDIM getVelCompTaylorCoeffs(const RealVect   & a_x0,
                                             const int        & a_velcomp) const;
  
  
    virtual RealVect getVelocityAtPoint(const RealVect& a_x) const;
  
    virtual ~SolidBodyRotVel()
    {
    }
  
  private:
    //weak construction is bad
    SolidBodyRotVel();
  
    RealVect m_center;
    Real     m_omega;
  
  };
 
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

    static inline Real getTaylorCoefficient(const RealVect   & a_x0,
                                            const IndMomDIM  & a_p,
                                            const IntVect    & a_power,
                                            const RealVect   & a_point
                                           )
    {
      Real retval = 1.0;
      for(int idir = 0; idir < DIM; idir++)
      {
        if(a_p[idir] > a_power[idir])
        {
          return 0.;
        }
        else
        {
          Real prefix = 1.0;
          for(int ipre = 0; ipre < a_p[idir]; ipre++)
          {
            prefix *= (a_power[idir] - ipre);
          }
          retval *= prefix*pow(a_x0[idir]-a_point[idir],a_power[idir]-a_p[idir])/factorial(a_p[idir]);
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
				Vector<Vector<IndMomDIM> > & a_ebfneighNorm,
				Real                     a_dx)
    {
      m_volMom.resize(a_volMom.size());
      Real volScaleFactor  = (D_TERM(a_dx, *a_dx, *a_dx));
      Real areaScaleFactor = (D_TERM(1.0,  *a_dx, *a_dx));
      for(int i = 0; i < a_volMom.size(); ++i)
      {
        IndMomDIM momspt = a_volMom[i];
        for(MomentIterator<DIM, order> iter; iter.ok(); ++iter)
        {
          const IndexTM<int, DIM>& p = iter();
          m_volMom[i][iter()] = momspt[iter()]/(volScaleFactor*POW(a_dx,p));
        }
      }
      m_facMom.resize(a_facMom.size());
      for(int i = 0; i < a_facMom.size(); ++i)
      {
        IndMomSDMinOne momspt = a_facMom[i];
        for(MomentIterator<DIM-1, order> iter; iter.ok(); ++iter)
        {
          const IndexTM<int, DIM-1>& p = iter();
          m_facMom[i][iter()] = momspt[iter()]/(areaScaleFactor*POW(a_dx,p));
        }
      }
      m_ebfMom.resize(a_ebfMom.size());
      for(int i = 0; i < a_ebfMom.size(); ++i)
      {
        IndMomDIM momspt = a_ebfMom[i];
        for(MomentIterator<DIM, order> iter; iter.ok(); ++iter)
        {
          const IndexTM<int, DIM>& p = iter();
          m_ebfMom[i][iter()] = momspt[iter()]/(areaScaleFactor*POW(a_dx,p));
        }
      }
      m_ebfneighNorm.resize(a_ebfneighNorm.size());
      for (int i = 0; i < a_ebfneighNorm.size(); ++i)
      {
        m_ebfneighNorm[i].resize(a_ebfneighNorm[i].size());
	for(int j = 0; j < a_ebfneighNorm[i].size(); ++j)
        {
          IndMomDIM momspt = a_ebfneighNorm[i][j];
          for(MomentIterator<DIM, order> iter; iter.ok(); ++iter)
          {
            const IndexTM<int, DIM>& p = iter();
            m_ebfneighNorm[i][j][iter()] = momspt[iter()]/(areaScaleFactor*POW(a_dx,p));
          }
        }
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
	              const EBIndex<CELL> & a_start) const;

    //Get neighbors VoFs of a given face within a specific radius.
    //Radius starts at one for the cells sharing the face, and 
    //increases by one for the neighbor cells sharing a face with the original
    //cells of radius one.
    Vector<EBIndex<CELL> >
    getVoFsNeighbors(const EBIndex<cent> & a_face,
                     const EBIndex<CELL> & a_vof,
                     const int             a_radius,
                     Neighbors<order>    & a_neigh,
                     const EBGraph       & a_graph,
                     int                   a_neighType) const;

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
                                const Neighbors<order> & a_neigh,
                                int                      a_neighType,
                                string                   a_ebbcName,
                                string                   a_dombcName[2*DIM],
                                double                   a_dx) const;

    //Compute the Moments Matrix
    void
    getShiftedMomentMatrix(Chombo4::LAPACKMatrix  & a_Mvol,
                           RealVect               & a_faceLoc,
                           double                   a_dx,
                           const Neighbors<order> & a_neigh,
                           string                   a_ebbcName,
                           string                   a_dombcName[2*DIM]) const;
		           
    //Compute rhs matrix for a face
    virtual void
    getFaceFluxMatrix(Chombo4::LAPACKMatrix    & a_Qshift,
                      const EBIndex<cent>      & a_face,
                      RealVect                 & a_faceLoc,
                      IndMomSDMinOne           & a_faceMom,
                      double                     a_signNormal,
                      double                     a_dx) const;

    //Compute rhs matrix for a EB face
    virtual void
    getEBFluxMatrix(Chombo4::LAPACKMatrix    & a_Qshift,
                    const EBIndex<cent>      & a_face,
                    RealVect                 & a_faceLoc,
                    Vector<IndMomDIM>        & a_ebnorm,
                    double                     a_dx) const;

    //Compute weight matrix 
    virtual void
    getWeightsMatrixFace(Chombo4::LAPACKMatrix  & a_weights,
                         const Neighbors<order> & a_neigh,
                         RealVect               & a_faceLoc) const;

    //Compute stencil coefficients matrix 
    void
    getStencilWeights(Chombo4::LAPACKMatrix        & a_phiCoef,
                      const Chombo4::LAPACKMatrix  & a_Qshift,
                      const Chombo4::LAPACKMatrix  & a_Mvol,
                      const Chombo4::LAPACKMatrix  & a_weights,
                      double                         a_dx) const;

    //Compute inhomogeneous term
    Real
    computeInhomogTerm(const Chombo4::LAPACKMatrix     & a_stencil,
                       const Neighbors<order>          & a_neigh,
                       const int                         a_nvolneigh,
                       string                            a_ebbcName,
                       string                            a_dombcName[2*DIM],
                       double                            a_dx) const;

    //Compute face centered flux stencil
    void
    getFaceCenteredStencil(LocalStencil<CELL, double> & a_sten,
                           Real                       & a_inhomogTerm, //for inhomogeneous bcs
                           const EBIndex<CELL>        & a_vof,
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
                           const int                    a_dir,
                           IndMomSDMinOne             & a_faceMom,
                           int                        & a_neighType,
                           double                       a_signNormal, 
                           double                       a_dx) const;

    //Compute EB face stencil
    void
    getEBFluxStencil(LocalStencil<CELL, double>    & a_facesten,
                     Real                          & a_inhomogTerm, //for inhomogeneous bcs
                     const EBIndex<CELL>           & a_vof,
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
                     Vector<IndMomDIM>             & a_ebnorm,
                     string                          a_ebbcName,
                     string                          a_dombcName[2*DIM],
                     int                             a_neighType,
                     const double                    a_dx) const;

    virtual ~BaseWLSFlux()
    {
    }

    BaseWLSFlux()
      {}

    BaseWLSFlux(int           a_radius,
                int           a_powerE,
                shared_ptr<BaseExactSolution<order> >  a_exactSol,
                int           a_fluxType=0,
                bool          a_verbose=false)
      {
        m_radius = a_radius;
        m_powerE = a_powerE;
        m_verbose = a_verbose;
        m_exactSol = a_exactSol;
        m_fluxType = a_fluxType;
      }

    BaseWLSFlux(int           a_radius,
                int           a_powerE,
                int           a_exactSolFun,
                bool          a_verbose=false)
      {
        m_radius = a_radius;
        m_powerE = a_powerE;
        m_verbose = a_verbose;
        ParmParse pp;

        if (a_exactSolFun == 1)
          {
            Vector<int> phiPowerX;
            Vector<int> phiPowerY;
            Vector<int> phiPowerZ;
            Vector<Real> phiCoef;
            int num_mono = 0;

            pp.get("phi_num_mono"   , num_mono);
            pp.getarr("phi_power_x",phiPowerX,0,num_mono);
            pp.getarr("phi_power_y",phiPowerY,0,num_mono);
            pp.getarr("phi_power_z",phiPowerZ,0,num_mono);
            pp.getarr("phi_coef",phiCoef,0,num_mono);

            vector<pair<Point, Real> > entries;
            for (int i=0; i<num_mono; i++)
              {
                std::pair<Point, Real> mono;
                mono.first = Point::Zeros();
                mono.first[0] = phiPowerX[i];
                mono.first[1] = phiPowerY[i];
#if DIM==3
                mono.first[2] = phiPowerZ[i];
#endif
                mono.second = phiCoef[i];
                entries.push_back(mono);
              }

            m_exactSol = shared_ptr<BaseExactSolution<order> >(new PolynomialEF<order>(entries));
          }
        else if (a_exactSolFun == 2 || a_exactSolFun == 3)
          {
            Real amplitude;
            pp.get("phi_amp",amplitude);

            IntVect frequency;
            RealVect center;
            Vector<int> freqVec;
            Vector<Real> centerVec;
  
        pp.getarr("phi_freq",freqVec,0,DIM);
        pp.getarr("phi_center",centerVec,0,DIM);
        for(int idir = 0; idir < DIM; idir++)
        {
          frequency[idir] = freqVec[idir];
          center[idir] = centerVec[idir];
        }

        if (a_exactSolFun == 2) m_exactSol = shared_ptr<BaseExactSolution<order> >(new OneDimCosProduct<order>(amplitude,frequency,center));
        if (a_exactSolFun == 3) m_exactSol = shared_ptr<BaseExactSolution<order> >(new OneDimSinProduct<order>(amplitude,frequency,center));
      }
      else if (a_exactSolFun == 4)
      {
        Real rad;
        pp.get("phi_rad", rad);
        Real cen;
        pp.get("phi_cen", cen);

        m_exactSol = shared_ptr<BaseExactSolution<order> >(new SineSphereEF<order>(rad, cen));
      }
      else if (a_exactSolFun == 5)
      {
        Real coeff;
        pp.get("phi_coeff",coeff);

        Real dfac;
        pp.get("phi_dfac",dfac);

        RealVect center;
        Vector<Real> centerVec;

        pp.getarr("phi_center",centerVec,0,DIM);
        for(int idir = 0; idir < DIM; idir++)
        {
          center[idir] = centerVec[idir];
        }

        m_exactSol = shared_ptr<BaseExactSolution<order> >(new GaussianBlob<order>(coeff,dfac,center));
      }
    }

    public:
    //Radius needed for neighbors list around a face (Manhattan style)
    int m_radius = 1;
    //Exponent for the distance between the VoFs and the face in the weights matrix
    int m_powerE = -3;
    //Boolean to enable print out statements for debugging
    bool m_verbose = false;
    //Flux type. 0-Laplacian, 1-Advection
    int m_fluxType = 0;
    //Exact solution class holding different function types
    shared_ptr<BaseExactSolution<order> >  m_exactSol;
    
  };

  //Main function calls for all WLS cell calculations, which follow the same process as the flux calculations
  template<int order>
  class BaseWLSCell : public BaseWLSFlux<CELL, order>
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

    //Compute rhs matrix for a face
    void
    getEvalCellMatrix(Chombo4::LAPACKMatrix    & a_Qshift,
                      const EBIndex<CELL>      & a_cell,
                      RealVect                 & a_cellLoc,
                      IndMomDIM                & a_cellMom,
                      double                     a_dx) const;

    //Compute cell stencil
    void
    getCellStencil(LocalStencil<CELL, double> & a_sten,
                   Real                       & a_inhomogTerm, //for inhomogeneous bcs
                   const EBIndex<CELL>        & a_vof,
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
                   const int                    a_dir,
                   IndMomDIM                  & a_cellMom,
                   double                       a_dx) const;

    virtual ~BaseWLSCell()
    {
    }

    BaseWLSCell()
      {}

    BaseWLSCell(int           a_radius,
                int           a_powerE,
                shared_ptr<BaseExactSolution<order> >  a_exactSol,
                bool          a_verbose=false)
      :BaseWLSFlux<CELL, order>(a_radius, a_powerE, a_exactSol, 0, a_verbose)
      {}

    BaseWLSCell(int           a_radius,
                int           a_powerE,
                int           a_exactSolFun,
                bool          a_verbose=false)
      :BaseWLSFlux<CELL, order>(a_radius, a_powerE, a_exactSolFun, a_verbose)
      {}
  };

  //Main function calls for ExactVelAdvectionWLS flux calculations
  template<CENTERING cent, int order>
  class ExactVelAdvectionWLSFlux: public BaseWLSFlux<cent,order>
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

    using BaseWLSFlux<cent,order>::BaseWLSFlux;

    //Compute rhs matrix for a face
    virtual void
    getFaceFluxMatrix(Chombo4::LAPACKMatrix    & a_Qshift,
                      const EBIndex<cent>      & a_face,
                      RealVect                 & a_faceLoc,
                      IndMomSDMinOne           & a_faceMom,
                      double                     a_signNormal,
                      double                     a_dx) const;

    //Compute rhs matrix for a EB face
    virtual void
    getEBFluxMatrix(Chombo4::LAPACKMatrix    & a_Qshift,
                    const EBIndex<cent>      & a_face,
                    RealVect                 & a_faceLoc,
                    Vector<IndMomDIM>        & a_ebnorm,
                    double                     a_dx) const;

    //Compute weight matrix 
    virtual void
    getWeightsMatrixFace(Chombo4::LAPACKMatrix  & a_weights,
                         const Neighbors<order> & a_neigh,
                         RealVect               & a_faceLoc) const;

    virtual ~ExactVelAdvectionWLSFlux()
    {
    }

    ExactVelAdvectionWLSFlux(int          a_radius,
                             int          a_powerE,
                             int          a_exactSolFun,
                             bool         a_verbose=false)
      :BaseWLSFlux<cent, order>(a_radius,
                                a_powerE,
                                a_exactSolFun,
                                a_verbose)
      {
        this->m_fluxType = 1;
        ParmParse pp;
        int exactVel = 1;
        //Exact solution
        pp.get("exactVel",exactVel);

      if (exactVel == 1)
      {
        RealVect centerVel;
        Real omega;
        pp.get("vel_omega",omega);
  
        Vector<Real> centerVec;
  
        pp.getarr("vel_center",centerVec,0,DIM);
        for(int idir = 0; idir < DIM; idir++)
        {
          centerVel[idir] = centerVec[idir];
        }
  
	m_exactVel = shared_ptr<BaseExactVelocity<order> >(new SolidBodyRotVel<order>(omega, centerVel));
      }
      else if (exactVel == 2)
      {
        RealVect centerVel;
        Vector<Real> velVec;

        pp.getarr("vel_vec",velVec,0,DIM);
        for(int idir = 0; idir < DIM; idir++)
        {
          centerVel[idir] = velVec[idir];
        }
#if DIM==2
        Real mag = pow(pow(centerVel[0],2)+pow(centerVel[1],2),0.5);
#elif DIM==3
        Real mag = pow(pow(centerVel[0],2)+pow(centerVel[1],2)+pow(centerVel[2],2),0.5);
#endif
        centerVel = centerVel/mag;

        m_exactVel = shared_ptr<BaseExactVelocity<order> >(new ConstantVectorVel<order>(centerVel));
      }    
      else if (exactVel == 3)
      {
        RealVect centerVel;
        Vector<Real> centerVec;

        pp.getarr("vel_center",centerVec,0,DIM);
        for(int idir = 0; idir < DIM; idir++)
        {
          centerVel[idir] = centerVec[idir];
        }

        Real centerMagVel;
        Real wallVel;
        Real radiusVel;
        pp.get("centerMagVel",centerMagVel);
        pp.get("wallVel",wallVel);
        pp.get("radiusVel",radiusVel);
 
        m_exactVel = shared_ptr<BaseExactVelocity<order> >(new PoiseuilleVectorVel<order>(centerMagVel, wallVel, radiusVel, centerVel));
      }     
 
    }

    //Exact velocity class holding different function types
    shared_ptr<BaseExactVelocity<order> >  m_exactVel;

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
    //Gets Laplacian WLS flux divergence stencil
    inline void
    getLaplaWLSDivFStencil(LocalStencil< CELL,  double>        & a_vofStencil,
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
                           double                                a_dx);

    //Gets Advection WLS flux divergence stencil
    inline void
    getAdvWLSDivFStencil(LocalStencil< CELL,  double>        & a_vofStencil,
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
                         double                                a_dx);

    //Gets WLS stencil for each component the gradient of a scalar
    // return a vector of stencils, one per component
    inline void
    getWLSGradStencil(LocalStencil<CELL, double>    & a_vofStencilx,
                      LocalStencil<CELL, double>    & a_vofStencily,
                      LocalStencil<CELL, double>    & a_vofStencilz,
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
                      double                                a_dx);

        //Gets WLS flux divergence stencil in a direction
    // return a vector of stencils, one per component
    inline void
    getWLSDivFCompStencil(LocalStencil<CELL, double>    & a_vofStencilx,
                          LocalStencil<CELL, double>    & a_vofStencily,
                          LocalStencil<CELL, double>    & a_vofStencilz,
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
                          double                                a_dx);

    inline
    LocalStencil< CELL,  double>
    getLapFluxStencil(LocalStencil< CELL,  double> & a_fluxMStencil, 
                      LocalStencil< CELL,  double> & a_fluxPStencil, 
                      RealVect                     & a_centroid,
                      const EBIndex<CELL>          & a_voflo, 
                      const EBIndex<CELL>          & a_vofhi, 
                      int                            a_dir,
                      const VoluData               & a_voludata, 
                      const EBGraph                & a_graph, 
                      double                         a_dx);

    inline
    double
    getLapinhomogTerm(double                       & a_inhomogTermM,
                      double                       & a_inhomogTermP, 
                      RealVect                     & a_centroid,
                      const EBIndex<CELL>          & a_voflo, 
                      const EBIndex<CELL>          & a_vofhi,
                      int                            a_dir,
                      const VoluData               & a_voludata,
                      const EBGraph                & a_graph,
                      double                         a_dx);
 
    inline
    LocalStencil< CELL,  double>
    getAdvFluxStencil(LocalStencil< CELL,  double> & a_fluxMStencil, 
                      LocalStencil< CELL,  double> & a_fluxPStencil, 
                      RealVect                     & a_centroid,
                      const EBIndex<CELL>          & a_voflo, 
                      const EBIndex<CELL>          & a_vofhi, 
                      const RealVect               & a_v_faceLoc, 
                      int                            a_dir,
                      const VoluData               & a_voludata, 
                      const EBGraph                & a_graph, 
                      double                         a_dx);

    virtual ~WLSStencil()
    {
    }

    WLSStencil()
    {
    }

    std::map<EBIndex<XFACE>,LocalStencil<CELL, double> > xfacesFluxStencil;
    std::map<EBIndex<YFACE>,LocalStencil<CELL, double> > yfacesFluxStencil;
    std::map<EBIndex<ZFACE>,LocalStencil<CELL, double> > zfacesFluxStencil;
    std::map<EBIndex<XFACE>,double> xfacesInhTerm;
    std::map<EBIndex<YFACE>,double> yfacesInhTerm;
    std::map<EBIndex<ZFACE>,double> zfacesInhTerm;
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
				       const shared_ptr<BaseExactSolution<order> >  & m_exactSol,
                                       const Real              a_dx);

    //getIntegralPhiAlongEB
    inline static Real getEBIntegral(const IndMomDIM  & a_ebfMom,
                                     const RealVect   & a_ebfLoc,
				     const shared_ptr<BaseExactSolution<order> > & m_exactSol,
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
				       const shared_ptr<BaseExactSolution<order> >  & m_exactSol,
		                       const Real              a_dx); 

    //getIntegralGradPhiDotNormal
    inline static Real getEBIntegral(const Vector<IndMomDIM> & a_ebfNormalMom,
                                     const RealVect          & a_ebfLoc,
				     const shared_ptr<BaseExactSolution<order> > & m_exactSol,
                                     const Real                a_dx);
 
  };

}
#include "Proto_WLSCellImplem.H"

#endif
  
