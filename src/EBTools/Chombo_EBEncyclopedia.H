
#ifndef _Chombo_EBEBEncyclopedia_H_
#define _Chombo_EBEBEncyclopedia_H_


#include "Proto.H"
#include "EBProto.H"
#include "Chombo_EBDictionary.H"
#include "Chombo_EBChombo.H"
#include "Chombo_GeometryService.H"
#include "Chombo_NamespaceHeader.H"


/// for x to x (same direction stencils)
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


/// for x to y (different direction stencils)
template <int order, typename data_t>
struct EBCrossFaceStencil
{
  shared_ptr<EBStencil< order, data_t, XFACE, YFACE> > m_xtoysten;
  shared_ptr<EBStencil< order, data_t, XFACE, ZFACE> > m_xtozsten;
  
  shared_ptr<EBStencil< order, data_t, YFACE, XFACE> > m_ytoxsten;
  shared_ptr<EBStencil< order, data_t, YFACE, ZFACE> > m_ytozsten;
  
  shared_ptr<EBStencil< order, data_t, ZFACE, XFACE> > m_ztoxsten;
  shared_ptr<EBStencil< order, data_t, ZFACE, YFACE> > m_ztoysten;

  //apply every stencil
  void apply(EBFluxData<data_t, 1>       & a_dst,
             const EBFluxData<data_t, 1> & a_src,
             bool a_initToZero, data_t a_scale)
  {
    m_xtoysten->apply(*a_dst.m_yflux, *a_src.m_xflux, a_initToZero, a_scale);
    m_ytoxsten->apply(*a_dst.m_xflux, *a_src.m_yflux, a_initToZero, a_scale);
    
#if DIM==3    
    m_xtozsten->apply(*a_dst.m_zflux, *a_src.m_xflux, a_initToZero, a_scale);
    m_ytozsten->apply(*a_dst.m_zflux, *a_src.m_yflux, a_initToZero, a_scale);

    m_ztoxsten->apply(*a_dst.m_xflux, *a_src.m_zflux, a_initToZero, a_scale);
    m_ztoysten->apply(*a_dst.m_yflux, *a_src.m_zflux, a_initToZero, a_scale);
#endif    
  }


  //apply one particular stencil
  void apply(EBFluxData<data_t, 1>       & a_dst,
             const EBFluxData<data_t, 1> & a_src,
             unsigned int a_dstDir, unsigned int a_srcDir,
             bool a_initToZero, data_t a_scale)
  {
    if((a_srcDir == 0) && (a_dstDir==1))
    {
      m_xtoysten->apply(*a_dst.m_yflux, *a_src.m_xflux, a_initToZero, a_scale);
    }
    else if((a_srcDir == 1) && (a_dstDir==0))
    {
      m_ytoxsten->apply(*a_dst.m_xflux, *a_src.m_yflux, a_initToZero, a_scale);
    }
#if DIM==3
    else if((a_srcDir == 0) && (a_dstDir==2))
    {
      m_xtozsten->apply(*a_dst.m_zflux, *a_src.m_xflux, a_initToZero, a_scale);
    }
    else if((a_srcDir == 1) && (a_dstDir==2))
    {
      m_ytozsten->apply(*a_dst.m_zflux, *a_src.m_yflux, a_initToZero, a_scale);
    }
    else if((a_srcDir == 2) && (a_dstDir==0))
    {
      m_ztoxsten->apply(*a_dst.m_xflux, *a_src.m_zflux, a_initToZero, a_scale);
    }
    else if((a_srcDir == 2) && (a_dstDir==1))
    {
      m_ztoysten->apply(*a_dst.m_yflux, *a_src.m_zflux, a_initToZero, a_scale);
    }
#endif
    else
    {
      MayDay::Error("ebcrossstencil::apply: bogus input directions");
    }
  }
};
///
/**
   EBEncylcopedia is  collection of EBDictionaries of all useful centerings.
   Nothing fancy.  The argument lists were just getting too long.
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

  shared_ptr<EBDictionary<order, data_t, XFACE   ,YFACE   > >    m_xFaceToYFace;
  shared_ptr<EBDictionary<order, data_t, XFACE   ,ZFACE   > >    m_xFaceToZFace;

  shared_ptr<EBDictionary<order, data_t, YFACE   ,XFACE   > >    m_yFaceToXFace;
  shared_ptr<EBDictionary<order, data_t, YFACE   ,ZFACE   > >    m_yFaceToZFace;

  shared_ptr<EBDictionary<order, data_t, ZFACE   ,XFACE   > >    m_zFaceToXFace;
  shared_ptr<EBDictionary<order, data_t, ZFACE   ,YFACE   > >    m_zFaceToYFace;

  ///register a stencil that goes from cell data  to face data here if it is an archived stencil
  inline void registerCellToFace (string a_stencilName,
                                  string a_domainBCName,
                                  string a_ebbcName,
                                  Box    a_srcDomain,
                                  Box    a_dstDomain,
                                  bool   a_needDiagonalWeights = false,
                                  Point  a_dstGrown = Point::Zeros())
  {
    m_cellToXFace->registerStencil(a_stencilName,              
                                   a_domainBCName,             
                                   a_ebbcName,                 
                                   a_srcDomain,                
                                   a_dstDomain,                
                                   a_needDiagonalWeights, 
                                   a_dstGrown);
    m_cellToYFace->registerStencil(a_stencilName,              
                                   a_domainBCName,             
                                   a_ebbcName,                 
                                   a_srcDomain,                
                                   a_dstDomain,                
                                   a_needDiagonalWeights, 
                                   a_dstGrown);
#if DIM==3      
    m_cellToZFace->registerStencil(a_stencilName,              
                                   a_domainBCName,             
                                   a_ebbcName,                 
                                   a_srcDomain,                
                                   a_dstDomain,                
                                   a_needDiagonalWeights, 
                                   a_dstGrown);
#endif      
  }

  ///
  /**
     This one is for algorithm developers
     If you know what stencil you want and you just want to use EBStencil to execute the stencil, 
     use this option.    
     Only the AggStencil bit will be executed so there has to be an entry
     corresponding to  every cell for which  you want to have an answer.   You have to call this for 
     every box as well.

     In 2D, z arguments do not matter.   I just hate dimension-dependent APIs.
  */
  inline void
  registerCellToFace (string                                            a_stencilName,
                      string                                            a_ebbcName,
                      const std::vector<EBIndex<     XFACE>         > & a_destXFace,
                      const std::vector<LocalStencil<CELL,  data_t> > & a_stenXFace,
                      const std::vector<EBIndex<     YFACE>         > & a_destYFace,
                      const std::vector<LocalStencil<CELL,  data_t> > & a_stenYFace,
                      const std::vector<EBIndex<     ZFACE>         > & a_destZFace,
                      const std::vector<LocalStencil<CELL,  data_t> > & a_stenZFace,
                      const Proto::Box                                & a_srcValid,
                      const Proto::Box                                & a_dstValid,
                      const Proto::Box                                & a_srcDomain,
                      const Proto::Box                                & a_dstDomain,
                      const Point                                     & a_srcGhost,
                      const Point                                     & a_dstGhost,
                      bool                                              a_needDiagonalWeights,
                      unsigned int                                      a_ibox)                
  {
    m_cellToXFace->registerStencil(a_stencilName,              
                                    a_ebbcName,
                                    a_destXFace,
                                    a_stenXFace,
                                    a_srcValid,
                                    a_dstValid,
                                    a_srcDomain,
                                    a_dstDomain,
                                    a_srcGhost,
                                    a_dstGhost,
                                    a_needDiagonalWeights,
                                    a_ibox);
    
    m_cellToYFace->registerStencil(a_stencilName,              
                                    a_ebbcName,
                                    a_destYFace,
                                    a_stenYFace,
                                    a_srcValid,
                                    a_dstValid,
                                    a_srcDomain,
                                    a_dstDomain,
                                    a_srcGhost,
                                    a_dstGhost,
                                    a_needDiagonalWeights,
                                    a_ibox);
#if DIM==3      
    m_cellToZFace->registerStencil(a_stencilName,              
                                    a_ebbcName,
                                    a_destZFace,
                                    a_stenZFace,
                                    a_srcValid,
                                    a_dstValid,
                                    a_srcDomain,
                                    a_dstDomain,
                                    a_srcGhost,
                                    a_dstGhost,
                                    a_needDiagonalWeights,
                                    a_ibox);
#endif      
  }
  



  ///  apply all the stencils
  inline void applyCellToFace(string                                a_stenName,
                              string                                a_ebbcName,
                              const Box                    &        a_domain,
                              EBFluxData<data_t, 1>        &        a_fluxData,
                              EBBoxData<CELL, data_t, DIM> &        a_cellData,
                              int a_ibox, bool a_initToZero, data_t a_scale)
  {

    {
      EBBoxData<CELL, data_t, 1> cellComp;
      unsigned int comp = 0;
      cellComp.define(a_cellData, comp);
      const auto& sten = m_cellToXFace->getEBStencil(a_stenName, a_ebbcName, a_domain, a_domain, a_ibox);
      sten->apply(*a_fluxData.m_xflux, cellComp, a_initToZero, a_scale);
    }

    {
      EBBoxData<CELL, data_t, 1> cellComp;
      unsigned int comp = 1;
      cellComp.define(a_cellData, comp);
      const auto& sten = m_cellToYFace->getEBStencil(a_stenName, a_ebbcName, a_domain, a_domain, a_ibox);
      sten->apply(*a_fluxData.m_yflux, cellComp, a_initToZero, a_scale);

    }
#if DIM==3      
    {
      EBBoxData<CELL, data_t, 1> cellComp;
      unsigned int comp = 2;
      cellComp.define(a_cellData, comp);
      const auto& sten = m_cellToZFace->getEBStencil(a_stenName, a_ebbcName, a_domain, a_domain, a_ibox);
      sten->apply(*a_fluxData.m_zflux, cellComp, a_initToZero, a_scale);
    }
#endif
  }


  ///  apply the appropriate stencil
  inline void applyCellToFace(string                          a_stenName,
                              string                          a_ebbcName,
                              const Box                   &   a_domain,
                              EBFluxData<data_t, 1>       &   a_fluxData,
                              EBBoxData<CELL, data_t, 1>  &   a_cellData,
                              unsigned int  a_idir, int a_ibox, bool a_initToZero, data_t a_scale)
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

  ///  apply the appropriate stencil
  inline void applyCellToFace(string                        a_stenName,
                              string                        a_ebbcName,
                              const Box                 &   a_domain,
                              EBFluxData<data_t, 1>       &   a_fluxData,
                              EBBoxData<CELL, data_t, 1>  &   a_cellData,
                              int a_ibox, bool a_initToZero, data_t a_scale)
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
  

  ///  apply all the stencils
  inline void applyFaceToCell(string                        a_stenName,
                              string                        a_ebbcName,
                              const Box                 &   a_domain,
                              EBBoxData<CELL, data_t, 1>  &   a_cellData,
                              EBFluxData<data_t, 1>       &   a_fluxData,
                              int a_ibox, bool a_initToZero, data_t a_scale)
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
    {
      const auto& sten = m_zFaceToCell->getEBStencil(a_stenName, a_ebbcName, a_domain, a_domain, a_ibox);
      sten->apply(a_cellData, *a_fluxData.m_zflux,  initToZero, a_scale);
    }
#endif

  }


  ///  apply a particular stencil
  inline void applyFaceToCell(string                        a_stenName,
                              string                        a_ebbcName,
                              const Box                 &   a_domain,
                              EBBoxData<CELL, data_t, 1>  &   a_cellData,
                              EBFluxData<data_t, 1>       &   a_fluxData,
                              unsigned int a_idir, int a_ibox, bool a_initToZero, data_t a_scale)
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


  ///apply all the facetocell stencils
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
  //where the faces act only on stuff of the same flavor of face
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


  /// registerFaceStencil version for algorithm developers
  /**
     If you know what stencil you want and you just want to use EBStencil to execute the stencil, 
     use this registerFaceStencil option.    
     Only the AggStencil bit will be executed so there has to be an entry
     corresponding to  every cell for which  you want to have an answer.   You have to call this for 
     every box as well.

     In 2D, z arguments do not matter.   I just hate dimension-dependent APIs.
   */
  inline void
  registerFaceStencil(string                                            a_stencilName,
                      string                                            a_ebbcName,
                      const vector<EBIndex<     XFACE>          >     & a_destXFace,
                      const vector<LocalStencil<XFACE,  data_t> >     & a_stenXFace,
                      const vector<EBIndex<     YFACE>          >     & a_destYFace,
                      const vector<LocalStencil<YFACE,  data_t> >     & a_stenYFace,
                      const vector<EBIndex<     ZFACE>          >     & a_destZFace,
                      const vector<LocalStencil<ZFACE,  data_t> >     & a_stenZFace,
                      const Proto::Box                                & a_srcValid,
                      const Proto::Box                                & a_dstValid,
                      const Proto::Box                                & a_srcDomain,
                      const Proto::Box                                & a_dstDomain,
                      const Point                                     & a_srcGhost,
                      const Point                                     & a_dstGhost,
                      bool                                              a_needDiagonalWeights,
                      unsigned int                                      a_ibox)                
  {
    m_xFaceToXFace->registerStencil(a_stencilName,              
                                    a_ebbcName,
                                    a_destXFace,
                                    a_stenXFace,
                                    a_srcValid,
                                    a_dstValid,
                                    a_srcDomain,
                                    a_dstDomain,
                                    a_srcGhost,
                                    a_dstGhost,
                                    a_needDiagonalWeights,
                                    a_ibox);

    m_yFaceToYFace->registerStencil(a_stencilName,              
                                    a_ebbcName,
                                    a_destYFace,
                                    a_stenYFace,
                                    a_srcValid,
                                    a_dstValid,
                                    a_srcDomain,
                                    a_dstDomain,
                                    a_srcGhost,
                                    a_dstGhost,
                                    a_needDiagonalWeights,
                                    a_ibox);

#if DIM==3
    m_zFaceToZFace->registerStencil(a_stencilName,              
                                    a_ebbcName,
                                    a_destZFace,
                                    a_stenZFace,
                                    a_srcValid,
                                    a_dstValid,
                                    a_srcDomain,
                                    a_dstDomain,
                                    a_srcGhost,
                                    a_dstGhost,
                                    a_needDiagonalWeights,
                                    a_ibox);
#endif    

  }
  
  //register a stencil for a face-centered thing
  //where faces act on different direction faces
  inline void registerCrossFaceStencil(string a_stencilName,
                                       string a_domainBCName,
                                       string a_ebbcName,
                                       Box    a_srcDomain,
                                       Box    a_dstDomain,
                                       bool   needDiagonalWeights = false)
  {
    m_xFaceToYFace->registerStencil(a_stencilName,              
                                    a_domainBCName,             
                                    a_ebbcName,                 
                                    a_srcDomain,                
                                    a_dstDomain,                
                                    needDiagonalWeights);
    m_yFaceToXFace->registerStencil(a_stencilName,              
                                    a_domainBCName,             
                                    a_ebbcName,                 
                                    a_srcDomain,                
                                    a_dstDomain,                
                                    needDiagonalWeights);
#if DIM==3      
    m_zFaceToXFace->registerStencil(a_stencilName,              
                                    a_domainBCName,             
                                    a_ebbcName,                 
                                    a_srcDomain,                
                                    a_dstDomain,                
                                    needDiagonalWeights);
    m_zFaceToYFace->registerStencil(a_stencilName,              
                                    a_domainBCName,             
                                    a_ebbcName,                 
                                    a_srcDomain,                
                                    a_dstDomain,                
                                    needDiagonalWeights);
    m_xFaceToZFace->registerStencil(a_stencilName,              
                                    a_domainBCName,             
                                    a_ebbcName,                 
                                    a_srcDomain,                
                                    a_dstDomain,                
                                    needDiagonalWeights);
    m_yFaceToZFace->registerStencil(a_stencilName,              
                                    a_domainBCName,             
                                    a_ebbcName,                 
                                    a_srcDomain,                
                                    a_dstDomain,                
                                    needDiagonalWeights);
#endif      
  }
  

  //get same face stencils all gathered together
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
  
  //get cross face stencils all gathered together
  EBCrossFaceStencil<order, data_t>
  getCrossFaceStencil(const string & a_stencilName, 
                      const string & a_ebbcName,
                      const Box    & a_srcDomain,
                      const Box    & a_dstDomain,               
                      const int    & a_boxid)
  {
    EBCrossFaceStencil<order, data_t> retval;
    retval.m_xtoysten = m_xFaceToYFace->getEBStencil(a_stencilName, a_ebbcName, a_srcDomain, a_dstDomain, a_boxid);
    retval.m_ytoxsten = m_yFaceToXFace->getEBStencil(a_stencilName, a_ebbcName, a_srcDomain, a_dstDomain, a_boxid);
#if DIM==3      
    retval.m_ztoxsten = m_zFaceToXFace->getEBStencil(a_stencilName, a_ebbcName, a_srcDomain, a_dstDomain, a_boxid);
    retval.m_ztoysten = m_zFaceToYFace->getEBStencil(a_stencilName, a_ebbcName, a_srcDomain, a_dstDomain, a_boxid);
    retval.m_ytozsten = m_yFaceToZFace->getEBStencil(a_stencilName, a_ebbcName, a_srcDomain, a_dstDomain, a_boxid);
    retval.m_xtozsten = m_xFaceToZFace->getEBStencil(a_stencilName, a_ebbcName, a_srcDomain, a_dstDomain, a_boxid);
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

    m_xFaceToYFace    = shared_ptr<EBDictionary<order, data_t, XFACE,   YFACE>    >(new EBDictionary<order, data_t, XFACE,   YFACE>    (a_geoserv, a_grids, a_domain, a_dx, a_ghost));
    m_xFaceToZFace    = shared_ptr<EBDictionary<order, data_t, XFACE,   ZFACE>    >(new EBDictionary<order, data_t, XFACE,   ZFACE>    (a_geoserv, a_grids, a_domain, a_dx, a_ghost));

    m_yFaceToXFace    = shared_ptr<EBDictionary<order, data_t, YFACE,   XFACE>    >(new EBDictionary<order, data_t, YFACE,   XFACE>    (a_geoserv, a_grids, a_domain, a_dx, a_ghost));
    m_yFaceToZFace    = shared_ptr<EBDictionary<order, data_t, YFACE,   ZFACE>    >(new EBDictionary<order, data_t, YFACE,   ZFACE>    (a_geoserv, a_grids, a_domain, a_dx, a_ghost));

    m_zFaceToXFace    = shared_ptr<EBDictionary<order, data_t, ZFACE,   XFACE>    >(new EBDictionary<order, data_t, ZFACE,   XFACE>    (a_geoserv, a_grids, a_domain, a_dx, a_ghost));
    m_zFaceToYFace    = shared_ptr<EBDictionary<order, data_t, ZFACE,   YFACE>    >(new EBDictionary<order, data_t, ZFACE,   YFACE>    (a_geoserv, a_grids, a_domain, a_dx, a_ghost));

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


    m_xFaceToYFace    = shared_ptr<EBDictionary<order, data_t, XFACE,   YFACE>    >(new EBDictionary<order, data_t, XFACE,   YFACE>    (a_geoserv, a_grids, a_domain, a_dx, a_ghost));
    m_xFaceToZFace    = shared_ptr<EBDictionary<order, data_t, XFACE,   ZFACE>    >(new EBDictionary<order, data_t, XFACE,   ZFACE>    (a_geoserv, a_grids, a_domain, a_dx, a_ghost));

    m_yFaceToXFace    = shared_ptr<EBDictionary<order, data_t, YFACE,   XFACE>    >(new EBDictionary<order, data_t, YFACE,   XFACE>    (a_geoserv, a_grids, a_domain, a_dx, a_ghost));
    m_yFaceToZFace    = shared_ptr<EBDictionary<order, data_t, YFACE,   ZFACE>    >(new EBDictionary<order, data_t, YFACE,   ZFACE>    (a_geoserv, a_grids, a_domain, a_dx, a_ghost));

    m_zFaceToXFace    = shared_ptr<EBDictionary<order, data_t, ZFACE,   XFACE>    >(new EBDictionary<order, data_t, ZFACE,   XFACE>    (a_geoserv, a_grids, a_domain, a_dx, a_ghost));
    m_zFaceToYFace    = shared_ptr<EBDictionary<order, data_t, ZFACE,   YFACE>    >(new EBDictionary<order, data_t, ZFACE,   YFACE>    (a_geoserv, a_grids, a_domain, a_dx, a_ghost));


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
