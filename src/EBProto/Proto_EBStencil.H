
#ifndef _Proto_EBStencil_H_
#define _Proto_EBStencil_H_


#include "Proto_RealVect.H"
#include "Proto_Box.H"
#include "Proto_Point.H"
#include "implem/Proto_LSquares.H"
#include "Proto_BaseIF.H"
#include "implem/Proto_IrregNode.H"
#include "implem/Proto_CutCellMoments.H"
#include "Proto_HostBoxData.H"
#include "implem/Proto_HostIrregData.H"
#include "implem/Proto_PointSet.H"
//#include "implem/Proto_LevelData.H"
#include "Proto_EBGraph.H"
#include "implem/Proto_EBStencilArchive.H"
#include <implem/Proto_StencilImplem.H>
#include <implem/Proto_FusedStencil.H>

using std::shared_ptr;
using std::vector;

namespace Proto
{
  /// class to encapsulate evaluating a stencil for an embedded boundary calculation
  /**
     This is Proto::Stencil generalized to embedded boundaries.   One does not 
     generate these directly, you get it out of your EBDictionary once you have registered it.
     The only part of this API that the user should use are EBStencil::apply and EBStencil::getDiagonalWeights
  */
  template <int order, typename data_t, CENTERING srcCenter, CENTERING dstCenter>
  class EBStencil
  {

  public:
    typedef EBStencilArchive<srcCenter, dstCenter, order, data_t> archive_t;

    typedef IndexedMoments<DIM  , order> IndMomDIM;
    typedef IndexedMoments<DIM-1, order> IndMomSDMinOne;
    typedef MomentIterator<DIM  , order> MomItDIM;
    typedef MomentIterator<DIM-1, order> MomItSDMinOne;
    typedef HostIrregData<CELL    ,  IndMomDIM , 1>  VoluData;
    typedef HostIrregData<BOUNDARY,  IndMomDIM , 1>  EBFaData;
    typedef HostIrregData<XFACE, IndMomSDMinOne, 1>  XFacData;
    typedef HostIrregData<YFACE, IndMomSDMinOne, 1>  YFacData;
    typedef HostIrregData<ZFACE, IndMomSDMinOne, 1>  ZFacData;

    /// Constructor for algorithm developers
    /** 
        If you are designing a stencil and want it executed via EBStencil, use this constructor.
        Users wishing to just used named, existing stencils should use the other constructor.
        In both cases, the user really should get her EBStencil through EBDictionary::registerStencil 
        (which calls this constructor or the other constructor if needed). EBDictionary maintains
        a cache of registered stencils so they never have to be defined twice.

        Only the AggStencil bit will be executed so you have to create a stencil 
        for every point in the box (not just the cut cells)
    */
    EBStencil(const vector<EBIndex<dstCenter> >               & a_dstVoFs,
              const vector<LocalStencil<srcCenter,  data_t> > & a_stencil,
              const EBGraph                                   & a_srcGraph,
              const EBGraph                                   & a_dstGraph,
              const Box                                       & a_srcValid,
              const Box                                       & a_dstValid,
              const Box                                       & a_srcDomain,
              const Box                                       & a_dstDomain,
              const Point                                     & a_srcGhost,
              const Point                                     & a_dstGhost,
              bool                                              a_needDiagonalWeights);
  

    /// Constructor for named, archived stencils.
    /** 
        Users wishing to just used named, existing stencils should use this constructor.
        This just asks the EBStencilArchive for the stencil, so if the stencil is not there, 
        it will be a run time error.
        If you are designing a stencil and want it executed via EBStencil, use the other constructor
        so you can actually just specify what it  is.
        In both cases, the user really should get her EBStencil through EBDictionary::registerStencil 
       (which calls this constructor or the other constructor if needed). EBDictionary maintains
       a cache of registered stencils so they never have to be defined twice.
    */
    EBStencil(string              a_stencilName,
              string              a_domainBCName[2*DIM],
              string              a_ebbcName,
              const Box         & a_srcValid,
              const Box         & a_dstValid,
              const Box         & a_srcDomain,
              const Box         & a_dstDomain,
              const Point       & a_srcGhost,
              const Point       & a_dstGhost,
              const EBGraph     & a_srcGraph,
              const EBGraph     & a_dstGraph,
              const VoluData    & a_voludata,
              const EBFaData    & a_ebfadata,
              const XFacData    & a_xfacdata,
              const YFacData    & a_yfacdata,
              const ZFacData    & a_zfacdata,
              double              a_dx,
              bool                a_periodiconly,
              bool                a_needDiagonalWeights,
              Point               a_dstGrown);
    ~EBStencil()
    {
//      using std::cout;
//      using std::endl;
//      cout << "in ebstencil destructor" << endl;
    }
    ///
    /**
       regular stencil gets applied in the interior of the domain.
       The aggstencil overwrites it near the EB.
       The domainBC stencils fill ghost cells outside the domain.  
       The vector length needs to be the same as the one-sided span of the stencil.
       
       Boundary conditions are applied if m_hasDomainBCs has been set to zero.
       See how  complicated this is?  Use EBDictionary.
    */
    inline
    void apply(EBBoxData<       dstCenter, data_t, 1>  &   a_dst,
               const EBBoxData< srcCenter, data_t, 1>  &   a_src,
               bool   a_initToZero = false,
               double a_scale      = 1)  const;
               

    
    inline
    void defineApplyITZ();
    
    // this forces initToZero = true
    inline
    void applyITZ(EBBoxData<       dstCenter, data_t, 1>  &   a_dst,
                  const EBBoxData< srcCenter, data_t, 1>  &   a_src,
                  double a_scale      = 1,
                  bool   a_initToZero = false) const;

               

    /// get the diagonal weights of the stencil
    /**
       can only be used if it was asked for at construction.
       Useful for relaxation methods.
     */
    shared_ptr<EBBoxData<dstCenter, data_t, 1> > getDiagonalWeights() const
    {
      return m_diagonalWeights;
    }
  protected:
    void
    computeDiagonalWeights(const vector< EBIndex<dstCenter> >              & a_dstVoFs,
                           const vector< LocalStencil<srcCenter, data_t> > & a_stencil,
                           const EBGraph                                   & a_graph,
                           const Box                                       & a_validBox,
                           const Point                                     & a_dstGrown);
    bool m_hasDomainBCs;
    bool m_bcsOnly;
    bool m_irregOnly;
    ///these 
    Stencil<double> m_regStencilInterior;
    //placeholder.  not sure what to do with domain bcs.
    //vector because there can be multiple ghost  cells
    vector<Stencil<double> >  m_regStencilBCS[2*DIM];
    vector<Box>               m_BCApplyBoxes[2*DIM];
    bool                      m_periodic_domain_bcs_only;
    Box                       m_regApplyBox;
    shared_ptr< AggStencil<srcCenter, dstCenter, data_t> >  m_irregStencil;
    shared_ptr<  EBBoxData<dstCenter, data_t, 1> >          m_diagonalWeights;

// disable for debugging
//#ifdef PROTO_CUDA
//    FusedStencil<data_t> sFused;
//#endif

    private:
    EBStencil();
    EBStencil(     const EBStencil<order, data_t, srcCenter, dstCenter>& a_inp);
    void operator=(const EBStencil<order, data_t, srcCenter, dstCenter>& a_inp);

  };
}
#include "implem/Proto_EBStencilImplem.H"
#endif
