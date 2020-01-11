
#ifndef _Chombo_EBEBEncyclopedia_H_
#define _Chombo_EBEBEncyclopedia_H_


#include "Proto.H"
#include "EBProto.H"
#include "Chombo_EBDictionary.H"
#include "Chombo_EBChombo.H"
#include "Chombo_GeometryService.H"
#include "Chombo_NamespaceHeader.H"


///
template <int order, typename data_t>
struct EBFluxStencil
{
  shared_ptr<EBStencil< order, data_t, XFACE, XFACE> > m_xsten;
  shared_ptr<EBStencil< order, data_t, YFACE, YFACE> > m_ysten;
  shared_ptr<EBStencil< order, data_t, ZFACE, ZFACE> > m_zsten;

  void apply(EBFluxData<data_t, 1>       & a_dst,
             const EBFluxData<data_t, 1> & a_src,
             bool a_initToZero, data_t a_scale)
  {
    m_xsten->apply(*a_dst.m_xflux, *a_src.m_xflux, a_initToZero, a_scale);
    m_ysten->apply(*a_dst.m_yflux, *a_src.m_yflux, a_initToZero, a_scale);
#if DIM==3    
    m_zsten->apply(*a_dst.m_zflux, *a_src.m_zflux, a_initToZero, a_scale);
#endif    
  }
};
///
/**
   EBEncylcopedia is  collection of EBDictionaries of all useful centerings.
   Nothing fancy.  The argument lists were just getting too long.
   Not all combinations are here but the ones I actually use are here.
   I realize this is not a great use of metaphor.
*/
template <int order, typename data_t>
class EBEncyclopedia
{

public:
  shared_ptr<EBDictionary<order, data_t, CELL    ,CELL    > >    m_cellToCell;
  shared_ptr<EBDictionary<order, data_t, XFACE   ,XFACE   > >    m_xFaceToXFace;
  shared_ptr<EBDictionary<order, data_t, YFACE   ,YFACE   > >    m_yFaceToYFace;
  shared_ptr<EBDictionary<order, data_t, ZFACE   ,ZFACE   > >    m_zFaceToZFace;
                                                            
  shared_ptr<EBDictionary<order, data_t, CELL    ,XFACE   > >    m_cellToXFace;
  shared_ptr<EBDictionary<order, data_t, CELL    ,YFACE   > >    m_cellToYFace;
  shared_ptr<EBDictionary<order, data_t, CELL    ,ZFACE   > >    m_cellToZFace;
  shared_ptr<EBDictionary<order, data_t, CELL    ,BOUNDARY> >    m_cellToBoundary;
                                                            
  shared_ptr<EBDictionary<order, data_t, XFACE   ,CELL    > >    m_xFaceToCell;
  shared_ptr<EBDictionary<order, data_t, YFACE   ,CELL    > >    m_yFaceToCell;
  shared_ptr<EBDictionary<order, data_t, ZFACE   ,CELL    > >    m_zFaceToCell;
  shared_ptr<EBDictionary<order, data_t, BOUNDARY,CELL    > >    m_boundaryToCell;

  //register a stencil cells to face things
  inline void registerCellToFace (string a_stencilName,
                                  string a_domainBCName,
                                  string a_ebbcName,
                                  Box    a_srcDomain,
                                  Box    a_dstDomain,
                                  bool   needDiagonalWeights = false,
                                  Point  a_dstGrown = Point::Zeros())
  {
    m_cellToXFace->registerStencil(a_stencilName,              
                                   a_domainBCName,             
                                   a_ebbcName,                 
                                   a_srcDomain,                
                                   a_dstDomain,                
                                   needDiagonalWeights, 
                                   a_dstGrown);
    m_cellToYFace->registerStencil(a_stencilName,              
                                   a_domainBCName,             
                                   a_ebbcName,                 
                                   a_srcDomain,                
                                   a_dstDomain,                
                                   needDiagonalWeights, 
                                   a_dstGrown);
#if DIM==3      
    m_cellToZFace->registerStencil(a_stencilName,              
                                   a_domainBCName,             
                                   a_ebbcName,                 
                                   a_srcDomain,                
                                   a_dstDomain,                
                                   needDiagonalWeights, 
                                   a_ pdstGrown);
#endif      
  }




  ///  apply all the stencils
  inline void applyCellToFace(string                        a_stenName,
                              string                        a_ebbcName,
                              const Box                 &   a_domain,
                              EBFluxData<Real, 1>       &   a_fluxData,
                              EBBoxData<CELL, Real, 1>  &   a_cellData,
                              int a_ibox, bool a_initToZero, Real a_scale)
  {

    {
      const auto& sten = m_cellToXFace->getEBStencil(a_stenName, a_ebbcName, a_domain, a_domain, a_ibox);
      sten->apply(*a_fluxData.m_xflux, a_cellData, a_initToZero, a_scale);
    }

    {
      const auto& sten = m_cellToYFace->getEBStencil(a_stenName, a_ebbcName, a_domain, a_domain, a_ibox);
      sten->apply(*a_fluxData.m_yflux, a_cellData, a_initToZero, a_scale);

    }
#if DIM==3      

    {
      const auto& sten = m_cellToZFace->getEBStencil(a_stenName, a_ebbcName, a_domain, a_domain, a_ibox);
      sten->apply(*a_fluxData.m_zflux, a_cellData, a_initToZero, a_scale);
    }
#endif
  }


  ///  apply the appropriate stencil
  inline void applyCellToFace(string                        a_stenName,
                              string                        a_ebbcName,
                              const Box                 &   a_domain,
                              EBFluxData<Real, 1>       &   a_fluxData,
                              EBBoxData<CELL, Real, 1>  &   a_cellData,
                              unsigned int  a_idir, int a_ibox, bool a_initToZero, Real a_scale)
  {
    if(a_idir == 0)
    {
      const auto& sten = m_cellToXFace->getEBStencil(a_stenName, a_ebbcName, a_domain, a_domain, a_ibox);
      sten->apply(*a_fluxData.m_xflux, a_cellData, a_initToZero, a_scale);
    }
    else if(a_idir == 1)
    {
      const auto& sten = m_cellToYFace->getEBStencil(a_stenName, a_ebbcName, a_domain, a_domain, a_ibox);
      sten->apply(*a_fluxData.m_yflux, a_cellData, a_initToZero, a_scale);

    }
#if DIM==3      
    else if(a_idir == 2)
    {
      const auto& sten = m_cellToZFace->getEBStencil(a_stenName, a_ebbcName, a_domain, a_domain, a_ibox);
      sten->apply(*a_fluxData.m_zflux, a_cellData, a_initToZero, a_scale);
    }
#endif
    else
    {
      PROTO_ASSERT(false, "bogus idir");
    }
  }


  ///  apply all the stencils
  inline void applyFaceToCell(string                        a_stenName,
                              string                        a_ebbcName,
                              const Box                 &   a_domain,
                              EBBoxData<CELL, Real, 1>  &   a_cellData,
                              EBFluxData<Real, 1>       &   a_fluxData,
                              int a_ibox, bool a_initToZero, Real a_scale)
  {
    bool initToZero = a_initToZero;
    {
      const auto& sten = m_xFaceToCell->getEBStencil(a_stenName, a_ebbcName, a_domain, a_domain, a_ibox);
      sten->apply(a_cellData, *a_fluxData.m_xflux, initToZero, a_scale);
    }
    initToZero = false;
    {
      const auto& sten = m_yFaceToCell->getEBStencil(a_stenName, a_ebbcName, a_domain, a_domain, a_ibox);
      sten->apply(a_cellData, *a_fluxData.m_yflux,  initToZero, a_scale);

    }
#if DIM==3      
    else if(a_idir == 2)
    {
      const auto& sten = m_zFaceToCell->getEBStencil(a_stenName, a_ebbcName, a_domain, a_domain, a_ibox);
      sten->apply(a_cellData, *a_fluxData.m_zflux,  initToZero, a_scale);
    }
#endif

  }


  //  apply the right stencil
  inline void applyFaceToCell(string                        a_stenName,
                              string                        a_ebbcName,
                              const Box                 &   a_domain,
                              EBBoxData<CELL, Real, 1>  &   a_cellData,
                              EBFluxData<Real, 1>       &   a_fluxData,
                              unsigned int a_idir, int a_ibox, bool a_initToZero, Real a_scale)
  {
    if(a_idir == 0)
    {
      const auto& sten = m_xFaceToCell->getEBStencil(a_stenName, a_ebbcName, a_domain, a_domain, a_ibox);
      sten->apply(a_cellData, *a_fluxData.m_xflux, a_initToZero, a_scale);
    }
    else if(a_idir == 1)
    {
      const auto& sten = m_yFaceToCell->getEBStencil(a_stenName, a_ebbcName, a_domain, a_domain, a_ibox);
      sten->apply(a_cellData, *a_fluxData.m_yflux,  a_initToZero, a_scale);

    }
#if DIM==3      
    else if(a_idir == 2)
    {
      const auto& sten = m_zFaceToCell->getEBStencil(a_stenName, a_ebbcName, a_domain, a_domain, a_ibox);
      sten->apply(a_cellData, *a_fluxData.m_zflux,  a_initToZero, a_scale);
    }
#endif
    else
    {
      PROTO_ASSERT(false, "bogus idir");
    }
  }


    //register a stencil for a face-centered thing
  inline void registerFaceToCell (string a_stencilName,
                                  string a_domainBCName,
                                  string a_ebbcName,
                                  Box    a_srcDomain,
                                  Box    a_dstDomain,
                                  bool   needDiagonalWeights = false)
  {
    m_xFaceToCell->registerStencil(a_stencilName,              
                                   a_domainBCName,             
                                   a_ebbcName,                 
                                   a_srcDomain,                
                                   a_dstDomain,                
                                   needDiagonalWeights);
    m_yFaceToCell->registerStencil(a_stencilName,              
                                   a_domainBCName,             
                                   a_ebbcName,                 
                                   a_srcDomain,                
                                   a_dstDomain,                
                                   needDiagonalWeights);
#if DIM==3      
    m_zFaceToCell->registerStencil(a_stencilName,              
                                   a_domainBCName,             
                                   a_ebbcName,                 
                                   a_srcDomain,                
                                   a_dstDomain,                
                                   needDiagonalWeights);
#endif      
  }

  //register a stencil for a face-centered thing
  inline void registerFaceStencil(string a_stencilName,
                                  string a_domainBCName,
                                  string a_ebbcName,
                                  Box    a_srcDomain,
                                  Box    a_dstDomain,
                                  bool   needDiagonalWeights = false)
  {
    m_xFaceToXFace->registerStencil(a_stencilName,              
                                    a_domainBCName,             
                                    a_ebbcName,                 
                                    a_srcDomain,                
                                    a_dstDomain,                
                                    needDiagonalWeights);
    m_yFaceToYFace->registerStencil(a_stencilName,              
                                    a_domainBCName,             
                                    a_ebbcName,                 
                                    a_srcDomain,                
                                    a_dstDomain,                
                                    needDiagonalWeights);
#if DIM==3      
    m_zFaceToZFace->registerStencil(a_stencilName,              
                                    a_domainBCName,             
                                    a_ebbcName,                 
                                    a_srcDomain,                
                                    a_dstDomain,                
                                    needDiagonalWeights);
#endif      
  }


  //since the above does not work, let's do this
  EBFluxStencil<order, data_t>
  getFluxStencil(const string & a_stencilName, 
                 const string & a_ebbcName,
                 const Box    & a_srcDomain,
                 const Box    & a_dstDomain,               
                 const int    & a_boxid)
  {
    EBFluxStencil<order, data_t> retval;
    retval.m_xsten = m_xFaceToXFace->getEBStencil(a_stencilName, a_ebbcName, a_srcDomain, a_dstDomain, a_boxid);
    retval.m_ysten = m_yFaceToYFace->getEBStencil(a_stencilName, a_ebbcName, a_srcDomain, a_dstDomain, a_boxid);
#if DIM==3      
    retval.m_zsten = m_zFaceToZFace->getEBStencil(a_stencilName, a_ebbcName, a_srcDomain, a_dstDomain, a_boxid);
#endif
    return retval;
  }

  
  ///single level constructor
  EBEncyclopedia(const shared_ptr<GeometryService<order> >     & a_geoserv,
                 const DisjointBoxLayout                       & a_grids,
                 const Box                                     & a_domain,
                 data_t a_dx, Point  a_ghost)
  {
    m_cellToCell      = shared_ptr<EBDictionary<order, data_t, CELL,    CELL>     >(new EBDictionary<order, data_t, CELL,    CELL>     (a_geoserv, a_grids, a_domain, a_dx, a_ghost));
    m_xFaceToXFace    = shared_ptr<EBDictionary<order, data_t, XFACE,   XFACE>    >(new EBDictionary<order, data_t, XFACE,   XFACE>    (a_geoserv, a_grids, a_domain, a_dx, a_ghost));
    m_yFaceToYFace    = shared_ptr<EBDictionary<order, data_t, YFACE,   YFACE>    >(new EBDictionary<order, data_t, YFACE,   YFACE>    (a_geoserv, a_grids, a_domain, a_dx, a_ghost));
    m_zFaceToZFace    = shared_ptr<EBDictionary<order, data_t, ZFACE,   ZFACE>    >(new EBDictionary<order, data_t, ZFACE,   ZFACE>    (a_geoserv, a_grids, a_domain, a_dx, a_ghost));
    m_cellToXFace     = shared_ptr<EBDictionary<order, data_t, CELL,    XFACE>    >(new EBDictionary<order, data_t, CELL,    XFACE>    (a_geoserv, a_grids, a_domain, a_dx, a_ghost));
    m_cellToYFace     = shared_ptr<EBDictionary<order, data_t, CELL,    YFACE>    >(new EBDictionary<order, data_t, CELL,    YFACE>    (a_geoserv, a_grids, a_domain, a_dx, a_ghost));
    m_cellToZFace     = shared_ptr<EBDictionary<order, data_t, CELL,    ZFACE>    >(new EBDictionary<order, data_t, CELL,    ZFACE>    (a_geoserv, a_grids, a_domain, a_dx, a_ghost));
    m_cellToBoundary  = shared_ptr<EBDictionary<order, data_t, CELL,    BOUNDARY> >(new EBDictionary<order, data_t, CELL,    BOUNDARY> (a_geoserv, a_grids, a_domain, a_dx, a_ghost));
    m_xFaceToCell     = shared_ptr<EBDictionary<order, data_t, XFACE,   CELL>     >(new EBDictionary<order, data_t, XFACE,   CELL>     (a_geoserv, a_grids, a_domain, a_dx, a_ghost));
    m_yFaceToCell     = shared_ptr<EBDictionary<order, data_t, YFACE,   CELL>     >(new EBDictionary<order, data_t, YFACE,   CELL>     (a_geoserv, a_grids, a_domain, a_dx, a_ghost));
    m_zFaceToCell     = shared_ptr<EBDictionary<order, data_t, ZFACE,   CELL>     >(new EBDictionary<order, data_t, ZFACE,   CELL>     (a_geoserv, a_grids, a_domain, a_dx, a_ghost));
    m_boundaryToCell  = shared_ptr<EBDictionary<order, data_t, BOUNDARY,CELL>     >(new EBDictionary<order, data_t, BOUNDARY,CELL>     (a_geoserv, a_grids, a_domain, a_dx, a_ghost));
  }

  ///multi-level constructor
  EBEncyclopedia(const shared_ptr<GeometryService<order> >   & a_geoserv,
                 const vector<DisjointBoxLayout>             & a_grids,
                 const vector<Box>                           & a_domain,
                 const vector<data_t>                        & a_dx,
                 Point  a_ghost)                      
  {
    m_cellToCell      = shared_ptr<EBDictionary<order, data_t, CELL,    CELL>     >(new EBDictionary<order, data_t, CELL,    CELL>     (a_geoserv, a_grids, a_domain, a_dx, a_ghost));
    m_xFaceToXFace    = shared_ptr<EBDictionary<order, data_t, XFACE,   XFACE>    >(new EBDictionary<order, data_t, XFACE,   XFACE>    (a_geoserv, a_grids, a_domain, a_dx, a_ghost));
    m_yFaceToYFace    = shared_ptr<EBDictionary<order, data_t, YFACE,   YFACE>    >(new EBDictionary<order, data_t, YFACE,   YFACE>    (a_geoserv, a_grids, a_domain, a_dx, a_ghost));
    m_zFaceToZFace    = shared_ptr<EBDictionary<order, data_t, ZFACE,   ZFACE>    >(new EBDictionary<order, data_t, ZFACE,   ZFACE>    (a_geoserv, a_grids, a_domain, a_dx, a_ghost));
    m_cellToXFace     = shared_ptr<EBDictionary<order, data_t, CELL,    XFACE>    >(new EBDictionary<order, data_t, CELL,    XFACE>    (a_geoserv, a_grids, a_domain, a_dx, a_ghost));
    m_cellToYFace     = shared_ptr<EBDictionary<order, data_t, CELL,    YFACE>    >(new EBDictionary<order, data_t, CELL,    YFACE>    (a_geoserv, a_grids, a_domain, a_dx, a_ghost));
    m_cellToZFace     = shared_ptr<EBDictionary<order, data_t, CELL,    ZFACE>    >(new EBDictionary<order, data_t, CELL,    ZFACE>    (a_geoserv, a_grids, a_domain, a_dx, a_ghost));
    m_cellToBoundary  = shared_ptr<EBDictionary<order, data_t, CELL,    BOUNDARY> >(new EBDictionary<order, data_t, CELL,    BOUNDARY> (a_geoserv, a_grids, a_domain, a_dx, a_ghost));
    m_xFaceToCell     = shared_ptr<EBDictionary<order, data_t, XFACE,   CELL>     >(new EBDictionary<order, data_t, XFACE,   CELL>     (a_geoserv, a_grids, a_domain, a_dx, a_ghost));
    m_yFaceToCell     = shared_ptr<EBDictionary<order, data_t, YFACE,   CELL>     >(new EBDictionary<order, data_t, YFACE,   CELL>     (a_geoserv, a_grids, a_domain, a_dx, a_ghost));
    m_zFaceToCell     = shared_ptr<EBDictionary<order, data_t, ZFACE,   CELL>     >(new EBDictionary<order, data_t, ZFACE,   CELL>     (a_geoserv, a_grids, a_domain, a_dx, a_ghost));
    m_boundaryToCell  = shared_ptr<EBDictionary<order, data_t, BOUNDARY,CELL>     >(new EBDictionary<order, data_t, BOUNDARY,CELL>     (a_geoserv, a_grids, a_domain, a_dx, a_ghost));
  }
  
  ///after this is done, you cannot register any more stencils.  This will delete the moment information.
  void  close()
  {
    m_cellToCell    .close();
    m_xFaceToXFace  .close();
    m_yFaceToYFace  .close();
    m_zFaceToZFace  .close();
    m_cellToXFace   .close();
    m_cellToYFace   .close();
    m_cellToZFace   .close();
    m_cellToBoundary.close();
    m_xFaceToCell   .close();
    m_yFaceToCell   .close();
    m_zFaceToCell   .close();
    m_boundaryToCell.close();
  }
private:
  EBEncyclopedia();
  EBEncyclopedia(const EBEncyclopedia<order, data_t>& a_inp);
  void operator=(const EBEncyclopedia<order, data_t>& a_inp);


};

#include "Chombo_NamespaceFooter.H"
#endif
