
#ifndef _Chombo_EBEBEncyclopedia_H_
#define _Chombo_EBEBEncyclopedia_H_


#include "EBDictionary.H"
#include "Proto.H"
#include "EBProto.H"
#include "EBChombo.H"
#include "GeometryService.H"

#include "NamespaceHeader.H"


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

  //register a stencil for a face-centered thing
  inline void registerFaceStencil(int    a_faceDir,
                                  string a_stencilName,
                                  string a_domainBCName,
                                  string a_ebbcName,
                                  Box    a_srcDomain,
                                  Box    a_dstDomain,
                                  bool   needDiagonalWeights = false)
  {
    PR_assert((a_faceDir >= 0) && (a_faceDir < DIM));
    if(a_faceDir == 0)
    {
      m_xFaceToXFace->registerStencil(a_stencilName,              
                                      a_domainBCName,             
                                      a_ebbcName,                 
                                      a_srcDomain,                
                                      a_dstDomain,                
                                      needDiagonalWeights);
    }
    else if(a_faceDir == 1)
    {
      m_yFaceToYFace->registerStencil(a_stencilName,              
                                      a_domainBCName,             
                                      a_ebbcName,                 
                                      a_srcDomain,                
                                      a_dstDomain,                
                                      needDiagonalWeights);
    }
    else if(a_faceDir == 2)
    {
      m_zFaceToZFace->registerStencil(a_stencilName,              
                                      a_domainBCName,             
                                      a_ebbcName,                 
                                      a_srcDomain,                
                                      a_dstDomain,                
                                      needDiagonalWeights);
    }
  }


//  ///get the stencil for a particular box at a particular level of refinement
//  inline template<typename ebstencil_t>
//  shared_ptr< ebstencil_t >
//  getFaceEBStencil(int a_faceDir, 
//                   const string & a_stencilName, 
//                   const string & a_ebbcName,
//                   const Box    & a_srcDomain,
//                   const Box    & a_dstDomain,               
//                   const int    & a_boxid)
//  {
//    PR_assert((a_faceDir >= 0) && (a_faceDir < DIM));
//    if(a_faceDir == 0)
//    {
//      return 
//        m_xFaceToXFace->getEBStencil(a_stencilName,              
//                                     a_domainBCName,             
//                                     a_ebbcName,                 
//                                     a_srcDomain,                
//                                     a_dstDomain,                
//                                     a_boxid);
//    }
//    else if(a_faceDir == 1)
//    {
//      m_yFaceToYFace->getEBStencil(a_stencilName,              
//                                   a_domainBCName,             
//                                   a_ebbcName,                 
//                                   a_srcDomain,                
//                                   a_dstDomain,                
//                                   a_boxid);
//    }
//    else if(a_faceDir == 2)
//    {
//      m_zFaceToZFace->getEBStencil(a_stencilName,              
//                                   a_domainBCName,             
//                                   a_ebbcName,                 
//                                   a_srcDomain,                
//                                   a_dstDomain,                
//                                   a_boxid);
//    }
//  }




  ///single level constructor
  EBEncyclopedia(const shared_ptr<GeometryService<order> >     & a_geoserv,
                 const DisjointBoxLayout                       & a_grids,
                 const Box                                     & a_domain,
                 data_t a_dx, Point  a_srcGhost,  Point  a_dstGhost )
  {
    m_cellToCell      = shared_ptr<EBDictionary<order, data_t, CELL,    CELL>     >(new EBDictionary<order, data_t, CELL,    CELL>     (a_geoserv, a_grids, a_domain, a_dx, a_srcGhost, a_dstGhost));
    m_xFaceToXFace    = shared_ptr<EBDictionary<order, data_t, XFACE,   XFACE>    >(new EBDictionary<order, data_t, XFACE,   XFACE>    (a_geoserv, a_grids, a_domain, a_dx, a_srcGhost, a_dstGhost));
    m_yFaceToYFace    = shared_ptr<EBDictionary<order, data_t, YFACE,   YFACE>    >(new EBDictionary<order, data_t, YFACE,   YFACE>    (a_geoserv, a_grids, a_domain, a_dx, a_srcGhost, a_dstGhost));
    m_zFaceToZFace    = shared_ptr<EBDictionary<order, data_t, ZFACE,   ZFACE>    >(new EBDictionary<order, data_t, ZFACE,   ZFACE>    (a_geoserv, a_grids, a_domain, a_dx, a_srcGhost, a_dstGhost));
    m_cellToXFace     = shared_ptr<EBDictionary<order, data_t, CELL,    XFACE>    >(new EBDictionary<order, data_t, CELL,    XFACE>    (a_geoserv, a_grids, a_domain, a_dx, a_srcGhost, a_dstGhost));
    m_cellToYFace     = shared_ptr<EBDictionary<order, data_t, CELL,    YFACE>    >(new EBDictionary<order, data_t, CELL,    YFACE>    (a_geoserv, a_grids, a_domain, a_dx, a_srcGhost, a_dstGhost));
    m_cellToZFace     = shared_ptr<EBDictionary<order, data_t, CELL,    ZFACE>    >(new EBDictionary<order, data_t, CELL,    ZFACE>    (a_geoserv, a_grids, a_domain, a_dx, a_srcGhost, a_dstGhost));
    m_cellToBoundary  = shared_ptr<EBDictionary<order, data_t, CELL,    BOUNDARY> >(new EBDictionary<order, data_t, CELL,    BOUNDARY> (a_geoserv, a_grids, a_domain, a_dx, a_srcGhost, a_dstGhost));
    m_xFaceToCell     = shared_ptr<EBDictionary<order, data_t, XFACE,   CELL>     >(new EBDictionary<order, data_t, XFACE,   CELL>     (a_geoserv, a_grids, a_domain, a_dx, a_srcGhost, a_dstGhost));
    m_yFaceToCell     = shared_ptr<EBDictionary<order, data_t, YFACE,   CELL>     >(new EBDictionary<order, data_t, YFACE,   CELL>     (a_geoserv, a_grids, a_domain, a_dx, a_srcGhost, a_dstGhost));
    m_zFaceToCell     = shared_ptr<EBDictionary<order, data_t, ZFACE,   CELL>     >(new EBDictionary<order, data_t, ZFACE,   CELL>     (a_geoserv, a_grids, a_domain, a_dx, a_srcGhost, a_dstGhost));
    m_boundaryToCell  = shared_ptr<EBDictionary<order, data_t, BOUNDARY,CELL>     >(new EBDictionary<order, data_t, BOUNDARY,CELL>     (a_geoserv, a_grids, a_domain, a_dx, a_srcGhost, a_dstGhost));
  }

  ///multi-level constructor
  EBEncyclopedia(const shared_ptr<GeometryService<order> >   & a_geoserv,
                 const vector<DisjointBoxLayout>             & a_grids,
                 const vector<Box>                           & a_domain,
                 const vector<data_t>                        & a_dx,
                 Point  a_srcGhost,  Point  a_dstGhost)                      
  {
    m_cellToCell      = shared_ptr<EBDictionary<order, data_t, CELL,    CELL>     >(new EBDictionary<order, data_t, CELL,    CELL>     (a_geoserv, a_grids, a_domain, a_dx, a_srcGhost, a_dstGhost));
    m_xFaceToXFace    = shared_ptr<EBDictionary<order, data_t, XFACE,   XFACE>    >(new EBDictionary<order, data_t, XFACE,   XFACE>    (a_geoserv, a_grids, a_domain, a_dx, a_srcGhost, a_dstGhost));
    m_yFaceToYFace    = shared_ptr<EBDictionary<order, data_t, YFACE,   YFACE>    >(new EBDictionary<order, data_t, YFACE,   YFACE>    (a_geoserv, a_grids, a_domain, a_dx, a_srcGhost, a_dstGhost));
    m_zFaceToZFace    = shared_ptr<EBDictionary<order, data_t, ZFACE,   ZFACE>    >(new EBDictionary<order, data_t, ZFACE,   ZFACE>    (a_geoserv, a_grids, a_domain, a_dx, a_srcGhost, a_dstGhost));
    m_cellToXFace     = shared_ptr<EBDictionary<order, data_t, CELL,    XFACE>    >(new EBDictionary<order, data_t, CELL,    XFACE>    (a_geoserv, a_grids, a_domain, a_dx, a_srcGhost, a_dstGhost));
    m_cellToYFace     = shared_ptr<EBDictionary<order, data_t, CELL,    YFACE>    >(new EBDictionary<order, data_t, CELL,    YFACE>    (a_geoserv, a_grids, a_domain, a_dx, a_srcGhost, a_dstGhost));
    m_cellToZFace     = shared_ptr<EBDictionary<order, data_t, CELL,    ZFACE>    >(new EBDictionary<order, data_t, CELL,    ZFACE>    (a_geoserv, a_grids, a_domain, a_dx, a_srcGhost, a_dstGhost));
    m_cellToBoundary  = shared_ptr<EBDictionary<order, data_t, CELL,    BOUNDARY> >(new EBDictionary<order, data_t, CELL,    BOUNDARY> (a_geoserv, a_grids, a_domain, a_dx, a_srcGhost, a_dstGhost));
    m_xFaceToCell     = shared_ptr<EBDictionary<order, data_t, XFACE,   CELL>     >(new EBDictionary<order, data_t, XFACE,   CELL>     (a_geoserv, a_grids, a_domain, a_dx, a_srcGhost, a_dstGhost));
    m_yFaceToCell     = shared_ptr<EBDictionary<order, data_t, YFACE,   CELL>     >(new EBDictionary<order, data_t, YFACE,   CELL>     (a_geoserv, a_grids, a_domain, a_dx, a_srcGhost, a_dstGhost));
    m_zFaceToCell     = shared_ptr<EBDictionary<order, data_t, ZFACE,   CELL>     >(new EBDictionary<order, data_t, ZFACE,   CELL>     (a_geoserv, a_grids, a_domain, a_dx, a_srcGhost, a_dstGhost));
    m_boundaryToCell  = shared_ptr<EBDictionary<order, data_t, BOUNDARY,CELL>     >(new EBDictionary<order, data_t, BOUNDARY,CELL>     (a_geoserv, a_grids, a_domain, a_dx, a_srcGhost, a_dstGhost));
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

#include "NamespaceFooter.H"
#endif
