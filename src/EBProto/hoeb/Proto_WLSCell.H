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
  /// This exists for the meta data we need about stencil creation, and as an interface
  class StencilAttributes
  {
  public:
    void
    readInputs(const std::string a_PPprefix)
      {
        // FIXME use prefix, and better names
        ParmParse pp;
        pp.get("radius", m_radius);
        pp.get("powerE", m_weightPower);
        pp.query("neighListType", m_neighListType);
      }

    void
    defaultStencilAttributes(const int a_order)
      {
        // defaults are an over-estimate but sufficiently large, and the minimal reasonable weighting
        m_radius = a_order; // how many cells out to grab
        m_weightPower = -(a_order+1); // WLS weighting power dist^P
        m_neighListType = 0;
      }

    StencilAttributes()
      {
        readInputs("");
      }

    StencilAttributes(const int a_order)
      {
        defaultStencilAttributes(a_order);
      }

    inline int radius() const
      { return m_radius;}

    inline int weightPower() const
      { return m_weightPower;}

    inline int neighListType() const
      { return m_neighListType;}

    inline Real smallCellTol() const
      { return m_smallCellTol;}

  private:
    /// Defaults are generally invalid, should be set in the constructor
    // How many cells around the stencil center to grab
    int m_radius = 0;
    // The weighting power to use
    int m_weightPower = 0;
    // Neighbor list type. 0-FaceBased or 1-CellBased
    int m_neighListType = 0;
    // at what point do small cells stop being relevant
    Real m_smallCellTol = 1e-11; // FIXME real aware
  };

  template<int order> 
  class BaseExactVelocity // FIXME a subset of exact solutions, don't reinvent them here
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
  //FIXME this is pulling double duty as both general methods for WLS stencils and the laplacian stencil
  // split into a base class for moment/weight matrix building and inherit a laplacian operator
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

    // FIXME we need to stop using the dombcName strings. Maybe a single comp EBDomRegister is general enough?

    //Gets all the moments and geometric info for the neighbors cells, boundary faces
    //and EB faces. Stores all this info in the Neighbors object.
    // DEPRECATED now lives in neighbors
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
                                Neighbors<order>       & a_neigh,
                                int                      a_neighType,
                                string                   a_ebbcName,
                                string                   a_dombcName[2*DIM], // order is is not obvious
                                Real                     a_dx) const
      {
        CH_assert(false);
      }

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
                                  RealVect                 a_faceLoc,
                                  Real                     a_dx,
                                  const Neighbors<order> & a_neigh,
                                  IndexTM<int, 2>          a_start = IndexTM<int, 2>{0, 0}) const;

    //Compute the moments matrix portion for regular faces in the neighborhood
    void
    getShiftedMomentMatrixRegFaces(Chombo4::LAPACKMatrix  & a_Mblock,
                                   RealVect                 a_faceLoc,
                                   Real                     a_dx,
                                   const Neighbors<order> & a_neigh,
                                   string                   a_dombcName[2*DIM],
                                   IndexTM<int, 2>          a_start = IndexTM<int, 2>{0, 0}) const;

    //Compute the moments matrix portion for cut faces in the neighborhood
    void
    getShiftedMomentMatrixEBFaces(Chombo4::LAPACKMatrix  & a_Mblock,
                                  RealVect                 a_faceLoc,
                                  Real                     a_dx,
                                  const Neighbors<order> & a_neigh,
                                  string                   a_ebbcName,
                                  IndexTM<int, 2>          a_start = IndexTM<int, 2>{0, 0}) const;

    //Compute the moments matrix portion for cut faces in the neighborhood
    void
    getShiftedMomentMatrixEBNormalFaces(Chombo4::LAPACKMatrix  & a_Mblock,
                                        RealVect                 a_faceLoc,
                                        Real                     a_dx,
                                        const Neighbors<order> & a_neigh,
                                        int                      a_dir,
                                        string                   a_ebbcName,
                                        IndexTM<int, 2>          a_start = IndexTM<int, 2>{0, 0}) const;

    // Will a block fit, this is just some asserts
    inline void
    checkSubMatrixSize(Chombo4::LAPACKMatrix  & a_Mblock,
                       const IndexTM<int, 2>   a_blockSize,
                       const IndexTM<int, 2>   a_start) const;

    // Shift moments from a_origLoc ot a_newLoc. Put here to avoid reuse
    // REVIEW point to math
    inline void
    shiftMomentToLoc(IndMomDIM & a_moment,
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
                              IndexTM<int, 2>          a_start = IndexTM<int, 2>{0, 0}) const;

    //Compute rhs matrix for a face
    void
    getFaceOpEvalMatrix(Chombo4::LAPACKMatrix    & a_Op,
                        const EBIndex<cent>      & a_face,
                        const RealVect             a_faceLoc,
                        IndMomSDMinOne           & a_faceMom,
                        const Real                 a_dx,
                        const IndexTM<int, DIM>    a_deriv) const;

    //Compute rhs matrix for a face. Pure virtual here
    virtual void
    getFaceFluxMatrix(Chombo4::LAPACKMatrix    & a_Qshift,
                      const EBIndex<cent>      & a_face,
                      RealVect                 & a_faceLoc,
                      IndMomSDMinOne           & a_faceMom,
                      Real                     a_signNormal,
                      Real                     a_dx) const = 0;

    //Compute rhs matrix for a EB face. Pure virtual here
    virtual void
    getEBFluxMatrix(Chombo4::LAPACKMatrix    & a_Qshift,
                    const EBIndex<cent>      & a_face,
                    RealVect                 & a_faceLoc,
                    Vector<IndMomDIM>        & a_ebnorm,
                    Real                     a_dx) const = 0;

    //Compute the stencil from a weighted least squares system
    void
    computeWLSStencilValues(Chombo4::LAPACKMatrix        & a_stencil,
                            const Chombo4::LAPACKMatrix  & a_Qshift,
                            const Chombo4::LAPACKMatrix  & a_Mvol,
                            const Chombo4::LAPACKMatrix  & a_weights) const;

    virtual ~BaseWLSFlux()
    {
    }

    BaseWLSFlux()
      {}

    BaseWLSFlux(const StencilAttributes& a_stenAttb,
                bool          a_verbose=false)
      {
        m_stenAttb = a_stenAttb;
        m_verbose = a_verbose;
      }

    BaseWLSFlux(int           a_radius,
                int           a_powerE,
                int           a_exactSolFun,
                bool          a_verbose=false)
      {
        CH_assert(false); // DEPRECATED
      }

  public:
    //Boolean to enable print out statements for debugging
    bool m_verbose = false;
    //Flux type. 0-Laplacian, 1-Advection
    int m_fluxType = 0;
    // at what point to small cells become too small to deal with
    Real m_eps = 1e-11; // FIXME should be precision dependent

    // Common, flux independent, dials for how to create a stencil
    // Things like the starting neighborhood radius and weighting power
    // this is the only member data that makes sense in this class
    StencilAttributes m_stenAttb;
  };

  //Main function calls for DivWLS flux calculations
  template<CENTERING cent, int order>
  class LapWLSFlux: public BaseWLSFlux<cent,order>
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
    using BaseWLSFlux<cent,order>::m_stenAttb;
    using BaseWLSFlux<cent,order>::m_verbose;

    //Compute rhs matrix for a face
    virtual void
    getFaceFluxMatrix(Chombo4::LAPACKMatrix    & a_Qshift,
                      const EBIndex<cent>      & a_face,
                      RealVect                 & a_faceLoc,
                      IndMomSDMinOne           & a_faceMom,
                      Real                       a_signNormal,
                      Real                       a_dx) const override;

    //Compute rhs matrix for a EB face
    virtual void
    getEBFluxMatrix(Chombo4::LAPACKMatrix    & a_Qshift,
                    const EBIndex<cent>      & a_face,
                    RealVect                 & a_faceLoc,
                    Vector<IndMomDIM>        & a_ebnorm,
                    Real                       a_dx) const override;

    //Compute inhomogeneous term
    Real
    computeInhomogTerm(const Chombo4::LAPACKMatrix     & a_stencil,
                       const Neighbors<order>          & a_neigh,
                       const int                         a_nvolneigh,
                       const EBDomainBC<1>      & a_domainBC,
                       Real                              a_dx) const;

    Real
    computeInhomogTermReg(const Chombo4::LAPACKMatrix     & a_stencil,
                          const Neighbors<order>          & a_neigh,
                          const int                         a_nvolneigh,
                          const EBDomainBC<1>      & a_domainBC,
                          Real                              a_dx) const;

    Real
    computeInhomogTermEB(const Chombo4::LAPACKMatrix     & a_stencil,
                         const Neighbors<order>          & a_neigh,
                         const int                         a_nvolneigh,
                         const EBDomainBC<1>      & a_domainBC,
                         Real                              a_dx) const;

    //Compute face centered flux stencil
    void
    getFaceCenteredStencil(LocalStencil<CELL, Real> & a_sten,
                           Real                       & a_inhomogTerm, //for inhomogeneous bcs
                           const EBIndex<CELL>        & a_vof,
                           const EBIndex<cent>        & a_face,
                           const EBGraph              & a_graph,
                           const MomentHolder<order>  & a_momentData,
                           const int                    a_dir,
                           IndMomSDMinOne             & a_faceMom,
                           int                        & a_neighType,
                           Real                         a_signNormal,
                           Real                         a_dx) const;

    //Compute EB face stencil
    void
    getEBFluxStencil(LocalStencil<CELL, Real>      & a_facesten,
                     Real                          & a_inhomogTerm, //for inhomogeneous bcs
                     const EBIndex<CELL>           & a_vof,
                     const EBIndex<cent>           & a_face,
                     const EBGraph                 & a_graph,
                     const MomentHolder<order>     & a_momentData,
                     Vector<IndMomDIM>             & a_ebnorm,
                     int                             a_neighType,
                     const Real                      a_dx) const;

    virtual ~LapWLSFlux()
      {}

    LapWLSFlux() = delete;

    LapWLSFlux(const StencilAttributes          & a_stenAttb,
               const EBDomainBC<1>   & a_domainBC,
               bool                           a_verbose=false)
      :BaseWLSFlux<cent, order>(a_stenAttb, a_verbose)
      {
        m_domainBC = a_domainBC;
      }

    EBDomainBC<1> m_domainBC;
  };

  // Methods for evaluating cell averaged gradient using the WLS construction method. same process as the flux calculations
  // TODO eventually the cell averaged reconstruction machinery could be separate from the operator
  template<int order>
  class GradWLSCell : public BaseWLSFlux<CELL, order>
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

    using BaseWLSFlux<CELL, order>::BaseWLSFlux;
    using BaseWLSFlux<CELL, order>::m_stenAttb;
    using BaseWLSFlux<CELL, order>::m_verbose;

    //Compute rhs matrix for a face. Doesn't make sense in the context
    virtual void
    getFaceFluxMatrix(Chombo4::LAPACKMatrix    & a_Qshift,
                      const EBIndex<CELL>      & a_face,
                      RealVect                 & a_faceLoc,
                      IndMomSDMinOne           & a_faceMom,
                      Real                     a_signNormal,
                      Real                     a_dx) const override
      {CH_assert(false)};

    //Compute rhs matrix for a EB face. Doesn't make sense in this context
    virtual void
    getEBFluxMatrix(Chombo4::LAPACKMatrix    & a_Qshift,
                    const EBIndex<CELL>      & a_face,
                    RealVect                 & a_faceLoc,
                    Vector<IndMomDIM>        & a_ebnorm,
                    Real                     a_dx) const override
    {CH_assert(false)};

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
                   const MomentHolder<order>  & a_momentData,
                   IndMomDIM                  & a_cellMom,
                   Real                       a_dx,
                   const IndexTM<int, DIM> a_deriv = IndexTM<int, DIM>::Zero) const;

    //Compute inhomogeneous term
    Real
    computeInhomogTerm(const Chombo4::LAPACKMatrix     & a_stencil,
                       const Neighbors<order>          & a_neigh,
                       const int                         a_nvolneigh,
                       const EBDomainBC<1>      & a_domainBC,
                       Real                              a_dx) const;

    Real
    computeInhomogTermReg(const Chombo4::LAPACKMatrix     & a_stencil,
                          const Neighbors<order>          & a_neigh,
                          const int                         a_nvolneigh,
                          const EBDomainBC<1>      & a_domainBC,
                          Real                              a_dx) const;

    Real
    computeInhomogTermEB(const Chombo4::LAPACKMatrix     & a_stencil,
                         const Neighbors<order>          & a_neigh,
                         const int                         a_nvolneigh,
                         const EBDomainBC<1>      & a_domainBC,
                         Real                              a_dx) const;

    ~GradWLSCell()
      {}

    GradWLSCell() = delete;

    GradWLSCell(const StencilAttributes          & a_stenAttb,
                const EBDomainBC<1>   & a_domainBC,
                bool                           a_verbose=false)
      :BaseWLSFlux<CELL, order>(a_stenAttb, a_verbose)
      {
        m_domainBC = a_domainBC;
      }

    EBDomainBC<1> m_domainBC;
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
                      Real                     a_dx) const override;

    //Compute rhs matrix for a EB face
    virtual void
    getEBFluxMatrix(Chombo4::LAPACKMatrix    & a_Qshift,
                    const EBIndex<cent>      & a_face,
                    RealVect                 & a_faceLoc,
                    Vector<IndMomDIM>        & a_ebnorm,
                    Real                     a_dx) const override;

    //Compute weight matrix 
    virtual void
    getWeightsMatrixFace(Chombo4::LAPACKMatrix  & a_weights,
                         const Neighbors<order> & a_neigh,
                         RealVect               & a_faceLoc) const override;

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
    using BaseWLSFlux<cent,order>::m_stenAttb;
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
                    Real                     a_dx) const override
      {
        CH_assert(false); // Don't need this function for divergence
      }

    //Compute face centered flux stencil, using all components and coupling BCs
    void
    getFaceCenteredStencilFull(std::array<LocalStencil<CELL, Real>, SpaceDim >& a_sten,
                               Real                       & a_inhomogTerm, //for inhomogeneous bcs
                               const EBIndex<CELL>        & a_vof,
                               const EBIndex<cent>        & a_face,
                               const EBGraph              & a_graph,
                               const MomentHolder<order>  & a_momentData,
                               const int                    a_dir,
                               IndMomSDMinOne             & a_faceMom,
                               int                        & a_neighType,
                               Real                       a_signNormal,
                               Real                       a_dx) const;

    //Compute EB face stencil
    void
    getEBFluxStencilFull(std::array<LocalStencil<CELL, Real>, SpaceDim >& a_sten,
                         Real                          & a_inhomogTerm, //for inhomogeneous bcs
                         const EBIndex<CELL>           & a_vof,
                         const EBIndex<cent>           & a_face,
                         const EBGraph                 & a_graph,
                         const MomentHolder<order>  & a_momentData,
                         int                             a_neighType,
                         const Real                    a_dx) const;

    Real
    computeInhomogTermReg(const Chombo4::LAPACKMatrix     & a_stencil,
                          const Neighbors<order>          & a_neigh,
                          const int                         a_dir,
                          const EBDomainBC<DIM>    & a_domainBC,
                          Real                            a_dx) const;

    Real
    computeInhomogTermEB(const Chombo4::LAPACKMatrix     & a_stencil,
                         const Neighbors<order>          & a_neigh,
                         const EBDomainBC<DIM>    & a_domainBC,
                         Real                              a_dx) const;

    // Constructor and destructor
    virtual ~DivWLSFlux()
      {}

    DivWLSFlux() = delete;

    DivWLSFlux(const StencilAttributes          & a_stenAttb,
               const EBDomainBC<DIM> & a_domainBC,
               bool          a_verbose=false)
      :BaseWLSFlux<cent, order>(a_stenAttb, a_verbose)
      {
        m_domainBC = a_domainBC;
      }

    EBDomainBC<DIM> m_domainBC;
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
                                       const shared_ptr<BaseExactSolution >  & m_exactSol,
                                       const Real              a_dx);

    //getIntegralPhiAlongEB
    inline static Real getEBIntegral(const IndMomDIM  & a_ebfMom,
                                     const RealVect   & a_ebfLoc,
                                     const shared_ptr<BaseExactSolution > & m_exactSol,
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
                                       const shared_ptr<BaseExactSolution >  & m_exactSol,
                                       const Real              a_dx,
                                       const int a_evalType = 0);

    //getIntegralGradPhiDotNormal
    inline static Real getEBIntegral(const Vector<IndMomDIM> & a_ebfNormalMom,
                                     const RealVect          & a_ebfLoc,
                                     const shared_ptr<BaseExactSolution > & m_exactSol,
                                     const Real                a_dx);
 
  };

}
#include "Proto_WLSCellImplem.H"

#endif
  
