#ifndef _Proto_WLSCell_H_
#define _Proto_WLSCell_H_

#include "Chombo_ProtoInterface.H"
#include "Chombo_LAPACKMatrix.H"
#include "EBProto.H"
#include "Proto_EBExactSolutions.H"
#include "Proto_MomentHolder.H"
#include "Proto_WLSneighbors.H"

namespace Proto
{

  /// A small class to hold parameters about how a stencil should be formed
  struct StencilTraits
  {
    void
    readStencilTraits(const std::string a_PPprefix)
      {
        ParmParse pp;
        pp.get("radius", m_radius);
        pp.get("powerE", m_weightPower);
        pp.query("neighListType", m_neighListType);
      }

    void
    defaultStencilTraits(const int a_order)
      {
        // defaults are an over-estimate but sufficiently large, and the minimal reasonable weighting
        int m_radius = a_order; // how many cells out to grab
        int m_weightPower = -(a_order+1); // WLS weighting power dist^P
        int m_neighListType = 0;
      }

    StencilTraits()
      {
        readStencilTraits("");
      }

    StencilTraits(const int a_order)
      {
        defaultStencilTraits(a_order);
      }

    inline int radius() const
      { return m_radius;}

    inline int weightPower() const
      { return m_weightPower;}

    inline int neighListType() const
      { return m_neighListType;}

    // How many cells around the stencil center to grab
    int m_radius = 0;
    // The weighting power to use
    int m_weightPower = 0;
    // Neighbor list type. 0-FaceBased or 1-CellBased
    int m_neighListType = 0;
    // at what point do small cells stop being relevant
    Real m_smallCellTol = 1e-11; // FIXME
  };

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

    virtual Real getDx() const
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

    //Gets all the moments and geometric info for the neighbors cells, boundary faces
    //and EB faces. Stores all this info in the Neighbors object.
    // DEPRECATED
    void
    getMomentsInfoFromNeighbors(Vector<EBIndex<CELL> > & a_vols,
                                const EBGraph          & a_graph,
                                const VoluData         & a_voludata, // FIXME so we just need the whole geoservg
                                const EBFaData         & a_ebfadata,
                                const XFacData         & a_xfacdata,
                                const YFacData         & a_yfacdata,
                                const ZFacData         & a_zfacdata,
                                const EBNormalData     & a_ebnormxdata,
                                const EBNormalData     & a_ebnormydata,
                                const EBNormalData     & a_ebnormzdata, // last geoserv argument
                                Neighbors<order> & a_neigh,
                                int                      a_neighType,
                                string                   a_ebbcName,
                                string                   a_dombcName[2*DIM], // order is is not obvious
                                Real                   a_dx) const;

    //Compute the Moments Matrix
    void
    getShiftedMomentMatrix(Chombo4::LAPACKMatrix  & a_Mvol,
                           RealVect               & a_faceLoc,
                           Real                     a_dx,
                           const Neighbors<order> & a_neigh,
                           string                   a_ebbcName,
                           string                   a_dombcName[2*DIM]) const;

    //Compute the moments matrix portion for cells in the neighborhood
    void
    getShiftedMomentMatrixVolumes(Chombo4::LAPACKMatrix  & a_Mblock,
                                  RealVect           a_faceLoc,
                                  Real               a_dx,
                                  const Neighbors<order> & a_neigh,
                                  IndexTM<int, 2>   a_start = IndexTM<int, 2>{0, 0}) const;

    //Compute the moments matrix portion for regular faces in the neighborhood
    void
    getShiftedMomentMatrixRegFaces(Chombo4::LAPACKMatrix  & a_Mblock,
                                   RealVect                 a_faceLoc,
                                   Real                     a_dx,
                                   const Neighbors<order> & a_neigh,
                                   string                   a_dombcName[2*DIM],
                                   IndexTM<int, 2>   a_start = IndexTM<int, 2>{0, 0}) const;

    //Compute the moments matrix portion for cut faces in the neighborhood
    void
    getShiftedMomentMatrixEBFaces(Chombo4::LAPACKMatrix  & a_Mblock,
                                  RealVect                 a_faceLoc,
                                  Real                     a_dx,
                                  const Neighbors<order> & a_neigh,
                                  string                   a_ebbcName,
                                  IndexTM<int, 2>   a_start = IndexTM<int, 2>{0, 0}) const;

    //Compute the moments matrix portion for cut faces in the neighborhood
    void
    getShiftedMomentMatrixEBNormalFaces(Chombo4::LAPACKMatrix  & a_Mblock,
                                        RealVect                 a_faceLoc,
                                        Real                     a_dx,
                                        const Neighbors<order> & a_neigh,
                                        int                      a_dir,
                                        string                   a_ebbcName,
                                        IndexTM<int, 2>   a_start = IndexTM<int, 2>{0, 0}) const;

    //Compute the moments matrix portion for cut faces in the neighborhood
    void
    getMomentMatrix(Chombo4::LAPACKMatrix   & a_Mblock,
                    const Vector<IndMomDIM> & a_moments,
                    const Vector<Real>      & a_fraction,
                    IndexTM<int, 2>   a_start = IndexTM<int, 2>{0, 0}) const;

    // Will a block fit
    void
    checkSubMatrixSize(Chombo4::LAPACKMatrix  & a_Mblock,
                       const IndexTM<int, 2>   a_blockSize,
                       const IndexTM<int, 2>   a_start) const;

    void
    shiftMomentsToLoc(IndMomDIM & a_moment,
                      RealVect    a_origLoc,
                      RealVect    a_newLoc,
                      Real        a_dx) const;

    //Compute weight matrix
    virtual void
    getWeightsMatrixFace(Chombo4::LAPACKMatrix  & a_weights,
                         const Neighbors<order> & a_neigh,
                         RealVect               & a_faceLoc) const;

    //Compute weight matrix for a set of distances
    void
    getWeightsMatrixComponent(Chombo4::LAPACKMatrix  & a_weights,
                              const Vector<int>      & a_weightDist,
                              IndexTM<int, 2>         a_start = IndexTM<Real, 2>{0, 0}) const;

    //Compute rhs matrix for a face
    virtual void
    getFaceFluxMatrix(Chombo4::LAPACKMatrix    & a_Qshift,
                      const EBIndex<cent>      & a_face,
                      RealVect                 & a_faceLoc,
                      IndMomSDMinOne           & a_faceMom,
                      Real                     a_signNormal,
                      Real                     a_dx) const;

    //Compute rhs matrix for a EB face
    virtual void
    getEBFluxMatrix(Chombo4::LAPACKMatrix    & a_Qshift,
                    const EBIndex<cent>      & a_face,
                    RealVect                 & a_faceLoc,
                    Vector<IndMomDIM>        & a_ebnorm,
                    Real                     a_dx) const;

    //Compute stencil coefficients matrix 
    void
    computeWLSStencilValues(Chombo4::LAPACKMatrix        & a_phiCoef,
                      const Chombo4::LAPACKMatrix  & a_Qshift,
                      const Chombo4::LAPACKMatrix  & a_Mvol,
                      const Chombo4::LAPACKMatrix  & a_weights) const;

    //Compute inhomogeneous term
    //Soon to be DEPRECATED
    Real
    computeInhomogTerm(const Chombo4::LAPACKMatrix     & a_stencil,
                       const Neighbors<order>          & a_neigh,
                       const int                         a_nvolneigh,
                       string                            a_ebbcName,
                       string                            a_dombcName[2*DIM],
                       Real                            a_dx) const;

    Real
    computeInhomogTermReg(const Chombo4::LAPACKMatrix     & a_stencil,
                          const Neighbors<order>          & a_neigh,
                          const int                         a_nvolneigh,
                          string                            a_dombcName[2*DIM],
                          Real                            a_dx) const;

    Real
    computeInhomogTermEB(const Chombo4::LAPACKMatrix     & a_stencil,
                         const Neighbors<order>          & a_neigh,
                         const int                         a_nvolneigh,
                         string                            a_ebbcName,
                         Real                            a_dx) const;

    //Compute face centered flux stencil
    //// DEPRECATED
    void
    getFaceCenteredStencil(LocalStencil<CELL, Real> & a_sten,
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
                           Real                       a_signNormal,
                           Real                       a_dx) const;

    //Compute face centered flux stencil
    // New version, use this
    void
    getFaceCenteredStencil(LocalStencil<CELL, Real> & a_sten,
                           Real                       & a_inhomogTerm, //for inhomogeneous bcs
                           const EBIndex<CELL>        & a_vof,
                           const EBIndex<cent>        & a_face,
                           const EBGraph              & a_graph,
                           const MomentHolder<order>  & a_momentData,
                           string                       a_ebbcName,
                           string                       a_dombcName[2*DIM],
                           const int                    a_dir,
                           IndMomSDMinOne             & a_faceMom,
                           int                        & a_neighType,
                           Real                       a_signNormal,
                           Real                       a_dx) const;

    //Compute EB face stencil
    // DEPRECATED
    void
    getEBFluxStencil(LocalStencil<CELL, Real>    & a_facesten,
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
                     const Real                    a_dx) const;

    //Compute EB face stencil
    // New version
    void
    getEBFluxStencil(LocalStencil<CELL, Real>    & a_facesten,
                     Real                          & a_inhomogTerm, //for inhomogeneous bcs
                     const EBIndex<CELL>           & a_vof,
                     const EBIndex<cent>           & a_face,
                     const EBGraph                 & a_graph,
                     const MomentHolder<order>  & a_momentData,
                     Vector<IndMomDIM>             & a_ebnorm,
                     string                          a_ebbcName,
                     string                          a_dombcName[2*DIM],
                     int                             a_neighType,
                     const Real                    a_dx) const;

    virtual ~BaseWLSFlux()
    {
    }

    BaseWLSFlux()
      {}

    BaseWLSFlux(int           a_radius,
                int           a_powerE,
                shared_ptr<BaseExactSolution<order> >  a_exactSol,
                shared_ptr<BaseExactSolution<order> >  a_exactSolEB = nullptr,
                int           a_fluxType=0,
                bool          a_verbose=false)
      {
        m_radius = a_radius;
        m_powerE = a_powerE;
        m_verbose = a_verbose;
        m_exactSolLo = a_exactSol;
        m_exactSolHi = a_exactSol;
        if (a_exactSolEB == nullptr)
          m_exactSolEB = a_exactSol;
        else
          m_exactSolEB = a_exactSolEB;
        m_fluxType = a_fluxType;
        CH_assert(false); //DEPRECATED use the updated constructor
      }

    BaseWLSFlux(int           a_radius,
                int           a_powerE,
                shared_ptr<BaseExactSolution<order> >  a_exactSolLo,
                shared_ptr<BaseExactSolution<order> >  a_exactSolHi,
                shared_ptr<BaseExactSolution<order> >  a_exactSolEB,
                int           a_fluxType=0,
                bool          a_verbose=false)
      {
        m_radius = a_radius;
        m_powerE = a_powerE;
        m_verbose = a_verbose;
        m_exactSolLo = a_exactSolLo;
        m_exactSolHi = a_exactSolHi;
        m_exactSolEB = a_exactSolEB;
        m_fluxType = a_fluxType;
        //CH_assert(false); //DEPRECATED use the updated constructor
      }

    // BaseWLSFlux(const StencilTraits          & a_stenTraits,
    //             const EBDomainBC<DIM, order> & a_domainBC
    //             bool          a_verbose=false)
    //   {
    //     m_radius = a_stenTraits.radius();
    //     m_powerE = a_stenTraits.weightPower();
    //     m_verbose = a_verbose;
    //     m_exactSolLo = a_exactSolLo;
    //     m_exactSolHi = a_exactSolHi;
    //     m_exactSolEB = a_exactSolEB;
    //     m_fluxType = 0;
    //   }

    BaseWLSFlux(int           a_radius,
                int           a_powerE,
                int           a_exactSolFun,
                bool          a_verbose=false)
      {
        CH_assert(false); // DEPRECATED THIS FUNCTION IS A HORRIBLE DESIGN, DONT DO THIS!
        m_radius = a_radius;
        m_powerE = a_powerE;
        m_verbose = a_verbose;
        ParmParse pp;
        shared_ptr<BaseExactSolution<order> > exactSol;
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

            exactSol = shared_ptr<BaseExactSolution<order> >(new PolynomialEF<order>(entries));
          }
        else if (a_exactSolFun == 2 || a_exactSolFun == 3)
          {
            Real amplitude;
            pp.get("phi_amp",amplitude);

            RealVect frequency;
            RealVect center;
            Vector<Real> freqVec;
            Vector<Real> centerVec;

            pp.getarr("phi_freq",freqVec,0,DIM);
            pp.getarr("phi_center",centerVec,0,DIM);
            for(int idir = 0; idir < DIM; idir++)
              {
                frequency[idir] = freqVec[idir];
                center[idir] = centerVec[idir];
              }

            if (a_exactSolFun == 2) exactSol = shared_ptr<BaseExactSolution<order> >(new OneDimCosProduct<order>(amplitude,frequency,center));
            if (a_exactSolFun == 3) exactSol = shared_ptr<BaseExactSolution<order> >(new OneDimSinProduct<order>(amplitude,frequency,center));
          }
        else if (a_exactSolFun == 4)
          {
            Real rad;
            pp.get("phi_rad", rad);
            Real cen;
            pp.get("phi_cen", cen);

            exactSol = shared_ptr<BaseExactSolution<order> >(new SineSphereEF<order>(rad, cen));
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

            exactSol = shared_ptr<BaseExactSolution<order> >(new GaussianBlob<order>(coeff,dfac,center));
          }
        m_exactSolLo = exactSol;
        m_exactSolHi = exactSol;
        m_exactSolEB = exactSol;
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
    shared_ptr<BaseExactSolution<order> >  m_exactSolLo;
    shared_ptr<BaseExactSolution<order> >  m_exactSolHi;
    shared_ptr<BaseExactSolution<order> >  m_exactSolEB;

    // Common, flux independent, dials for how to create a stencil
    // Things like the starting neighborhood radius and weighting power
    StencilTraits m_stenTraits;
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
                      Real                     a_dx,
                      const IndexTM<int, DIM> a_deriv = IndexTM<int, DIM>::Zero) const;

    //Compute cell stencil
    void
    getCellStencil(LocalStencil<CELL, Real> & a_sten,
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
                   IndMomDIM                  & a_cellMom,
                   Real                       a_dx,
                   const IndexTM<int, DIM> a_deriv = IndexTM<int, DIM>::Zero) const;

    virtual ~BaseWLSCell()
      {}

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
                      Real                     a_signNormal,
                      Real                     a_dx) const;

    //Compute rhs matrix for a EB face
    virtual void
    getEBFluxMatrix(Chombo4::LAPACKMatrix    & a_Qshift,
                    const EBIndex<cent>      & a_face,
                    RealVect                 & a_faceLoc,
                    Vector<IndMomDIM>        & a_ebnorm,
                    Real                     a_dx) const;

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

  //Main function calls for DivWLS flux calculations
  template<CENTERING cent, int order>
  class DivWLSFlux: public BaseWLSFlux<cent,order>
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
    using BaseWLSFlux<cent,order>::m_powerE;
    using BaseWLSFlux<cent,order>::m_radius;
    using BaseWLSFlux<cent,order>::m_verbose;

    //Compute rhs matrix for a face
    virtual void
    getFaceFluxMatrix(Chombo4::LAPACKMatrix    & a_Qshift,
                      const EBIndex<cent>      & a_face,
                      RealVect                 & a_faceLoc,
                      IndMomSDMinOne           & a_faceMom,
                      Real                     a_signNormal,
                      Real                     a_dx) const override;

    //Compute rhs matrix for a EB face
    virtual void
    getEBFluxMatrix(Chombo4::LAPACKMatrix    & a_Qshift,
                    const EBIndex<cent>      & a_face,
                    RealVect                 & a_faceLoc,
                    Vector<IndMomDIM>        & a_ebnorm,
                    Real                     a_dx) const override;

    //Compute face centered flux stencil, using all components and coupling BCs
    void
    getFaceCenteredStencilFull(std::array<LocalStencil<CELL, Real>, SpaceDim >& a_sten,
                               Real                       & a_inhomogTerm, //for inhomogeneous bcs
                               const EBIndex<CELL>        & a_vof,
                               const EBIndex<cent>        & a_face,
                               const EBGraph              & a_graph,
                               const MomentHolder<order>  & a_momentData,
                               string                       a_ebbcName,
                               string                       a_dombcName[2*DIM],
                               const int                    a_dir,
                               IndMomSDMinOne             & a_faceMom,
                               int                        & a_neighType,
                               Real                       a_signNormal,
                               Real                       a_dx) const;

    //Compute EB face stencil
    // New version
    void
    getEBFluxStencilFull(std::array<LocalStencil<CELL, Real>, SpaceDim >& a_sten,
                         Real                          & a_inhomogTerm, //for inhomogeneous bcs
                         const EBIndex<CELL>           & a_vof,
                         const EBIndex<cent>           & a_face,
                         const EBGraph                 & a_graph,
                         const MomentHolder<order>  & a_momentData,
                         Vector<IndMomDIM>             & a_ebnorm,
                         string                          a_ebbcName,
                         string                          a_dombcName[2*DIM],
                         int                             a_neighType,
                         const Real                    a_dx) const;

    virtual ~DivWLSFlux()
      {}

    DivWLSFlux()
      {}

    DivWLSFlux(int           a_radius,
               int           a_powerE,
               shared_ptr<BaseExactSolution<order> >  a_exactSolLo,
               shared_ptr<BaseExactSolution<order> >  a_exactSolHi,
               shared_ptr<BaseExactSolution<order> >  a_exactSolEB,
               bool          a_verbose=false)
      :BaseWLSFlux<cent, order>(a_radius, a_powerE, a_exactSolLo, a_exactSolHi, a_exactSolEB, 0, a_verbose)
      {}

    DivWLSFlux(const StencilTraits          & a_stenTraits,
               const EBDomainBC<DIM, order> & a_domainBC,
               bool          a_verbose=false)
      {
        m_radius = a_stenTraits.radius();
        m_powerE = a_stenTraits.weightPower();
        m_domainBC = a_domainBC;
      }

    EBDomainBC<DIM, order> m_domainBC;
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

  // template <int order>
  // template <>
  // inline void
  // WLSStencil<order>::
  //   getWLSDivFCompStencilFull<0>(
  //     std::array<LocalStencil<CELL, Real>, SpaceDim >& a_vofStencil,
  //     const EBIndex<CELL>                 & a_vof,
  //     const EBGraph                       & a_graph,
  //     const MomentHolder<order>           & a_momentData,
  //     const DivWLSFlux<XFACE   , order>  & a_xfaceflux, // Why build these separate?
  //     const DivWLSFlux<YFACE   , order>  & a_yfaceflux,
  //     const DivWLSFlux<ZFACE   , order>  & a_zfaceflux,
  //     const DivWLSFlux<BOUNDARY, order>  & a_ebflux,
  //     string                                a_ebbcName,
  //     vector<string>                                a_dombcName[2*DIM],
  //     Real                                & a_inhoTerm,
  //     int                                 a_neighType,
  //     Real                                a_dx)
  //     {
  //       getWLSDivFCompStencilDir<0>(a_vofStencil,
  //                                   a_vof,
  //                                   a_graph,
  //                                   a_momentData,
  //                                   a_xfaceflux,
  //                                   a_yfaceflux,
  //                                   a_zfaceflux,
  //                                   a_ebflux,
  //                                   a_ebbcName,
  //                                   a_dombcName,
  //                                   a_inhoTerm,
  //                                   a_neighType,
  //                                   a_dx);
  //     }

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
  
