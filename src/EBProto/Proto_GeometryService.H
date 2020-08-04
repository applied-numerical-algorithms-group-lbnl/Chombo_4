
#ifndef _Proto_GeometryService_H_
#define _Proto_GeometryService_H_


#include "EBProto.H"
#include "Proto_Box.H"
#include "Proto_Point.H"
#include "implem/Proto_LSquares.H"
#include "Proto_BaseIF.H"
#include "implem/Proto_IrregNode.H"
#include "implem/Proto_IrregData.H"
#include "implem/Proto_CutCellMoments.H"
#include "Proto_HostBoxData.H"
#include "implem/Proto_HostIrregData.H"
#include "implem/Proto_PointSet.H"
using std::shared_ptr;
using std::vector;

namespace Proto
{

  ///Makes geometric description from an implicit function.
  /**
     Given an implicit function and grid information, this cuts 
     the grids and provides moments.
  */
  template <int order>
  class GeometryService
  {

  public:

    typedef IndexedMoments<DIM  , order> IndMomDIM;
    typedef IndexedMoments<DIM-1, order> IndMomSDMinOne;
    typedef MomentIterator<DIM  , order> MomItDIM;
    typedef MomentIterator<DIM-1, order> MomItSDMinOne;
    typedef IrregNode<order> inode;
    typedef HostIrregData<CELL    ,  IndMomDIM , 1>  VoluData;
    typedef HostIrregData<BOUNDARY,  IndMomDIM , 1>  EBFaData;
    typedef HostIrregData<XFACE, IndMomSDMinOne, 1>  XFacData;
    typedef HostIrregData<YFACE, IndMomSDMinOne, 1>  YFacData;
    typedef HostIrregData<ZFACE, IndMomSDMinOne, 1>  ZFacData;

    ///
    /**
       basif      :  implicit function
       origin     :  location of the lower left corner of the domain
       dx         :  grid spacing
       domain     :  finest domain in the calculation
       finestGrids:  layout of the finest domain
       ghost      :  maximum number of ghost cells needed in calculation
       maxCoarsen :  maximum level of coarsening 
    */
    GeometryService(const shared_ptr<BaseIF>     &      a_baseIF,
                    const RealVect               &      a_origin,
                    const double                 &      a_dx,
                    const Box                    &      a_domain,
                    const DisjointBoxLayout      &      a_finestGrids,
                    const Point                  &      a_ghost,
                    int                                 a_maxCoarsen);



    ///
    /**
       Define the internals of the input ebisRegion.   Public for testing.
    */
    inline void fillGraph(HostBoxData<int>          &   a_regIrregCovered,
                          vector<inode>             &   a_nodes,
                          const Box                 &   a_validRegion,
                          const Box                 &   a_ghostRegion,
                          const Box                 &   a_domain) const;


    inline int getLevel(const Box & a_domain)
    {
      int retval = -1;
      bool found = false;
      for(int ilev = 0; ilev < m_domains.size(); ilev++)
      {
        if(m_domains[ilev] == a_domain)
        {
          retval = ilev;
          found = true;
          break;
        }
      }
      PROTO_ASSERT(found, "did not find asked-for domain");
      if(!found)
      {
        std::cout<< "did not find asked-for domain" << std::endl;
        abort();
      }
      return retval;
    }
    
    ///
    inline shared_ptr<LevelData<EBGraph>  >
    getGraphs(const Box& a_domain)
    {
      int lev = getLevel(a_domain);
      return   m_graphs[lev];
    }
    ///
    inline shared_ptr<LevelData<VoluData> >
    getVoluData(const Box& a_domain)
    {
      int lev = getLevel(a_domain);
      return m_voludata[lev];
    }
    ///
    inline shared_ptr<LevelData<EBFaData> >
    getEBFaceData(const Box& a_domain)
    {
      int lev = getLevel(a_domain);
      return m_ebfadata[lev];
    }
    ///
    inline shared_ptr<LevelData<XFacData> >
    getXFaceData(const Box& a_domain)
    {
      int lev = getLevel(a_domain);
      return m_xfacdata[lev];
    }
    ///
    inline shared_ptr<LevelData<YFacData> >
    getYFaceData(const Box& a_domain)
    {
      int lev = getLevel(a_domain);
      return m_yfacdata[lev];
    }
    ///
    inline shared_ptr<LevelData<ZFacData> >
    getZFaceData(const Box& a_domain)
    {
      int lev = getLevel(a_domain);
      return m_zfacdata[lev];
    }

  protected:

    void
    defineFinestLevel(const RealVect               & a_origin,
                      const double                 & a_dx );


    bool onBoxBoundary(const Point        & a_iv, 
                       const Box            & a_box,
                       const int            & a_dir,
                       const Side::LoHiSide & a_sd) const;



    //converts a RealVect in physical coordinates to a RealVect in coordinates relative to a cell center
    inline RealVect convert2RelativeCoord(const RealVect& a_rVect)const;

    //converts a IndexTM<double,DIM> in physical coordinates to a
    //RealVect in coordinates relative to a cell center
    inline RealVect convert2RelativeCoord(const IndexTM<double,DIM>& a_rVect)const;

    inline void
    fillNewNode(inode                    &     a_node,
                const PointSet           &     a_ivsIrreg,
                const Point              &     a_iv) const;


    inline shared_ptr<BaseIF> getBaseIF() const
    {
      return m_baseIF;
    }

  private:
    inline void
    fixRegularCellsNextToCovered(vector<inode>       & a_nodes, 
                                 HostBoxData<int>    & a_regIrregCovered,
                                 const Box           & a_validRegion,
                                 const Point         & a_iv) const;
    inline void
    getFullNodeWithCoveredFace(inode                   & a_newNode, 
                               const HostBoxData<int>  & a_regIrregCovered,
                               const Point             & a_iv) const;

  
  
    ///these vectors are of the length of number of coarsenings

    ///
    vector<Box>                                     m_domains;
    
    ///
    vector<DisjointBoxLayout>                       m_grids;
                                                    
    ///                                             
    vector<shared_ptr<LevelData<EBGraph> >  >       m_graphs;

    ///
    vector<shared_ptr<LevelData<VoluData> > >       m_voludata;

    ///
    vector<shared_ptr<LevelData<EBFaData> > >       m_ebfadata;

    ///
    vector<shared_ptr<LevelData<XFacData> > >       m_xfacdata;

    ///
    vector<shared_ptr<LevelData<YFacData> > >       m_yfacdata;

    ///
    vector<shared_ptr<LevelData<ZFacData> > >       m_zfacdata;
    
    double                m_threshold;

    //origin
    RealVect            m_origin;

    int                 m_order;
    int                 m_degreeP;
    double              m_dx;

    Box m_domain;
    Point m_geomgrow;

    shared_ptr<BaseIF>       m_baseIF;
  

    inline void fillArc(vector<int>                          a_arc[DIM],
                        CutCellMoments<DIM>            &     a_cutCellMoments,
                        const int                      &     a_hilo,
                        const PointSet                 &     a_ivsIrreg,
                        const Point                    &     a_curriv) const;


    //stuff disallowed for all the usual reasons.
  private:
    GeometryService();

    GeometryService(const GeometryService& a_workshopin);

    void operator=(const GeometryService& a_workshopin);


  };
}
#include "implem/Proto_GeometryServiceImplem.H"
#endif
