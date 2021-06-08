#ifndef __hoeb_utilities_implem_
#define __hoeb_utilities_implem_
#include "Hoeb_Neighborhood.H"
namespace hoeb
{
  /******/
  inline RealVect
  getVoFLocation(const EBIndex<CELL>& a_vof,
                 Real                 a_dx)
  {
    RealVect retval;
    auto pt = a_vof.m_pt; 
    for(int idir = 0; idir < DIM; idir++)
    {
      retval[idir] = a_dx*(0.5 + Real(pt[idir]));
    }
    return retval;
  }
  /******/
  inline 
  IndexTM<Real, DIM> getIndexDistance(const RealVect& a_vofloc,
                                      const RealVect& a_xbar)
  {
    IndexTM<Real, DIM> distance;
    for(int idir = 0; idir < DIM; idir++)
    {
      distance[idir] = a_vofloc[idir] - a_xbar[idir];
    }
    return distance;
  }
  /******/
  template <unsigned int D, CENTERING cent>
  inline IndexedMoments<D, HOEB_MAX_ORDER>   
  getIndexedMoment(const EBIndex<cent>                                              &   a_index,
                   const HostIrregData<cent, IndexedMoments<D, HOEB_MAX_ORDER>, 1>  &   a_irregfab,
                   const EBGraph                                                    &   a_graph,
                   const Real                                                       &   a_dx)
  {
    IndexedMoments<D, HOEB_MAX_ORDER> retval;
    
    if(a_graph.isRegular(a_index.m_pt))
    {
      retval.setToRegular(a_dx);
    }
    else if(a_graph.isIrregular(a_index.m_pt))
    {
      retval = a_irregfab(a_index, 0);
    }
    else
    {
      PROTO_ASSERT(false, "should not get here by construction");
    }
    return retval;
  }
  /***
  template <unsigned int D, CENTERING cent>
  inline IndexedMoments<D, HOEB_MAX_ORDER>   
  getIndexedMoment(const EBIndex<cent>                                 & a_index,
                   const shared_ptr< GeometryService<HOEB_MAX_ORDER> > & a_geoserv,
                   const Chombo4::Box                                  & a_srcDomain,
                   unsigned int                                          a_ibox,
                   Real                                                  a_dx)
  {
    IndexedMoments<D, HOEB_MAX_ORDER>   retval;
    PROTO_ASSERT(false, "should not get here");
    return retval;
  }

  ***/
  ///template specialization for x
  inline IndexedMoments<DIM-1, HOEB_MAX_ORDER>   
  getIndexedMoment(const EBIndex<XFACE>                                & a_index,
                   const shared_ptr< GeometryService<HOEB_MAX_ORDER> > & a_geoserv,
                   const Chombo4::Box                                  & a_srcDomain,
                   unsigned int                                          a_ibox,
                   Real                                                  a_dx)
  {
    auto dit = a_geoserv->getDBL(a_srcDomain).dataIterator();
    const auto & xfacldptr  =  a_geoserv->getXFaceData(a_srcDomain);
    const auto & xfacdatfab = (*xfacldptr)[dit[a_ibox]];
    const auto & graphsldptr =  a_geoserv->getGraphs(  a_srcDomain);
    const auto & graph     = (*graphsldptr)[dit[a_ibox]];
    return getIndexedMoment<DIM-1,XFACE>(a_index, xfacdatfab, graph, a_dx);
  }
  ///template specialization for y
  inline IndexedMoments<DIM-1, HOEB_MAX_ORDER>   
  getIndexedMoment(const EBIndex<YFACE>                                &   a_index,
                   const shared_ptr< GeometryService<HOEB_MAX_ORDER> > & a_geoserv,
                   const Chombo4::Box                                  & a_srcDomain,
                   unsigned int                                          a_ibox,
                   Real                                                  a_dx)
  {
    auto dit = a_geoserv->getDBL(a_srcDomain).dataIterator();
    const auto & yfacldptr  =  a_geoserv->getYFaceData(a_srcDomain);
    const auto & yfacdatfab = (*yfacldptr)[dit[a_ibox]];
    const auto & graphsldptr =  a_geoserv->getGraphs(  a_srcDomain);
    const auto & graph     = (*graphsldptr)[dit[a_ibox]];
    return getIndexedMoment<DIM-1, YFACE>(a_index, yfacdatfab, graph, a_dx);
  }
  ///template specialization for z
  inline IndexedMoments<DIM-1, HOEB_MAX_ORDER>   
  getIndexedMoment(const EBIndex<ZFACE>                                & a_index,
                   const shared_ptr< GeometryService<HOEB_MAX_ORDER> > & a_geoserv,
                   const Chombo4::Box                                  & a_srcDomain,
                   unsigned int                                          a_ibox,
                   Real                                                  a_dx,
                   int                                                   a_facedir)

  {
    auto dit = a_geoserv->getDBL(a_srcDomain).dataIterator();
    const auto & zfacldptr  =  a_geoserv->getZFaceData(a_srcDomain);
    const auto & zfacdatfab = (*zfacldptr)[dit[a_ibox]];
    const auto & graphsldptr =  a_geoserv->getGraphs(  a_srcDomain);
    const auto & graph     = (*graphsldptr)[dit[a_ibox]];
    return getIndexedMoment<DIM-1, ZFACE>(a_index, zfacdatfab, graph, a_dx);
  }
  
  template <CENTERING cent>
  inline RealVect
  getFaceLocation(const EBIndex<cent>& a_face,
                  Real                 a_dx,
                  int                  a_facedir)
  {
    RealVect retval;
    auto pt = a_face.m_pt; //high size of the face
    //this will work for cell centers or boundary faces as well
    for(int idir = 0; idir < DIM; idir++)
    {
      if(idir == a_facedir)
      {
        retval[idir] = a_dx*Real(pt[idir]);
      }
      else
      {
        retval[idir] = a_dx*(0.5 + Real(pt[idir]));
      }
    }
    return retval;
  }
  ///
  template <unsigned int D>
  IndexTM<int, D-1>
  getFaceIndexFromVol(const  IndexTM<int, D>& a_volind, int a_facedir)
  {
    IndexTM<int, D-1> retval;
    int ifac = 0;
    for(int ivol = 0; ivol < D; ivol++)
    {
      if(ivol != a_facedir)
      {
        retval[ifac] = a_volind[ivol];
        ifac++;
      }
    }
    return retval;
  }
  ///
  template <unsigned int D>
  IndexedMoments<D, HOEB_MAX_ORDER>
  getVolMomFromFacMom(const IndexedMoments<D-1, HOEB_MAX_ORDER>& a_facmom,
                      int a_facedir)
  {
    IndexedMoments<D, HOEB_MAX_ORDER> retval;
    for(MomentIterator< D, HOEB_MAX_ORDER> momit; momit.ok(); ++momit)
    {
      IndexTM<int, D> volpow = momit();
      Real momval;
      if(volpow[a_facedir] > 0) //no derivs in this direction
      {
        momval = 0;
      }
      else
      {
        IndexTM<int, D-1> facpow = getFaceIndexFromVol<DIM>(volpow, a_facedir);
        momval = a_facmom[facpow];
      }
      retval[volpow] = momval;
    }
    return retval;
  }
  /*************/
  inline IndexedMoments<DIM, HOEB_MAX_ORDER>
  getEBIndexedMoment(const EBIndex<BOUNDARY>                             & a_face,
                     const shared_ptr< GeometryService<HOEB_MAX_ORDER> > & a_geoserv,
                     const Chombo4::Box                                  & a_srcDomain,
                     const Chombo4::DataIndex                            & a_datind,
                     const EBGraph                                       & a_graph,
                     Real a_dx)
  {
    IndexedMoments<DIM, HOEB_MAX_ORDER> retval;
    //have to do this because the ebnormal moments have a cell type
    const auto   ebMomLDptr = a_geoserv->getEBFaceData(a_srcDomain);
    const auto & ebmomfab = (*ebMomLDptr)[a_datind];
    retval = getIndexedMoment<DIM, BOUNDARY>(a_face, ebmomfab, a_graph, a_dx);

    return retval;
  }
  /*************/
  inline IndexedMoments<DIM, HOEB_MAX_ORDER>
  getEBNormalMoment(const EBIndex<BOUNDARY>                             & a_face,
                    const shared_ptr< GeometryService<HOEB_MAX_ORDER> > & a_geoserv,
                    const Chombo4::Box                                  & a_srcDomain,
                    const Chombo4::DataIndex                            & a_datind,
                    const EBGraph                                       & a_graph,
                    unsigned int a_normDir, Real a_dx)
  {
    IndexedMoments<DIM, HOEB_MAX_ORDER> retval;
    //have to do this because the ebnormal moments have a cell type
    EBIndex<CELL> vof = a_face.getVoF();
    if(a_normDir == 0)
    {
      const auto   xnormMomLDptr = a_geoserv->getEBNormalData_x(a_srcDomain);
      const auto & xnormmomfab = (*xnormMomLDptr)[a_datind];
      retval = getIndexedMoment<DIM, CELL>(vof, xnormmomfab, a_graph, a_dx);
    }
    else if(a_normDir == 1)
    {
      const auto   ynormMomLDptr = a_geoserv->getEBNormalData_y(a_srcDomain);
      const auto & ynormmomfab = (*ynormMomLDptr)[a_datind];
      retval = getIndexedMoment<DIM, CELL>(vof, ynormmomfab, a_graph, a_dx);
    }
#if DIM==3
    else if(a_normDir == 2)
    {
      const auto   znormMomLDptr = a_geoserv->getEBNormalData_z(a_srcDomain);
      const auto & znormmomfab = (*znormMomLDptr)[a_datind];
      retval = getIndexedMoment<DIM,CELL>(vof, znormmomfab, a_graph, a_dx);
    }
#endif
    else
    {
      PROTO_ASSERT(false, "getshiftedindexedmoment: bogus normdir");
    }
    return retval;
  }
  ////
  inline void
  shiftMomentAndFillRow(LAPACKMatrix                            & a_mat,
                        IndexedMoments<    DIM, HOEB_MAX_ORDER> & a_bcmom,
                        const IndexTM<Real, DIM>                & a_xbardistance,
                        const unsigned int                      & a_currentRow)
  {
                          
    //shift so it has the same xbar
    a_bcmom.shift(a_xbardistance);
    for(MomentIterator<DIM, HOEB_MAX_ORDER> momit; momit.ok(); ++momit)
    {
      auto momind = momit();
      unsigned int currentCol = IndexedMoments<DIM, HOEB_MAX_ORDER>::indexOf(momind);
      a_mat(a_currentRow, currentCol) = a_bcmom[momind];
    }
  }
  /////
  inline void
  fillMatrixRowForEBBC(LAPACKMatrix                                        & a_mat,
                       const shared_ptr< GeometryService<HOEB_MAX_ORDER> > & a_geoserv,
                       const Chombo4::Box                                  & a_srcDomain,
                       const EBIndex<BOUNDARY>                             & a_face,
                       const EBGraph                                       & a_graph,
                       const IndexTM<Real, DIM>                            & a_xbardistance,
                       const Chombo4::DataIndex                            & a_datind,
                       const unsigned int                                  & a_currentRow,
                       int a_facedir, Real a_dx, string a_bcname)
  {
    IndexedMoments<    DIM, HOEB_MAX_ORDER> bcmom;
    if(a_bcname == string("Dirichlet"))
    {
      bcmom =  getEBIndexedMoment(a_face, a_geoserv, a_srcDomain,
                                  a_datind, a_graph, a_dx);
    }
    else if(a_bcname == string("Neumann"))
    {
      bcmom.setVal(0.);
      for(unsigned int normDir = 0; normDir < DIM; normDir++)
      {
        for(MomentIterator<DIM, HOEB_MAX_ORDER> momit; momit.ok(); ++momit)
        {
          auto momind = momit();
          Real bcval = 0;
          
          IndexedMoments<DIM, HOEB_MAX_ORDER> normMom=
            getEBNormalMoment(a_face, a_geoserv, a_srcDomain,
                              a_datind, a_graph, normDir, a_dx);

          if(momind[normDir] > 0)
          {
            int pd = momind[normDir];
            IndexTM<int, DIM> derind =momind - IndexTM<int, DIM>::BASISV(normDir);
            Real derval = normMom[derind];
            bcmom[momind] += pd*derval;
          }
        }
      }
    }
    shiftMomentAndFillRow(a_mat, bcmom, a_xbardistance, a_currentRow);
  }

  ///
  template < CENTERING cent>
  void
  fillMatrixRowForDomBC
  (LAPACKMatrix                                                         & a_mat,
   const HostIrregData<cent, IndexedMoments<DIM-1, HOEB_MAX_ORDER>, 1 > & a_facdatfab,
   const EBIndex<cent>                                                  & a_face,
   const EBGraph                                                        & a_graph,
   const IndexTM<Real, DIM>                                             & a_xbardistance,
   const unsigned int                                                   & a_currentRow,
   int a_facedir, Real  a_dx, string a_bcname)
  {
                       
    IndexedMoments<DIM-1, HOEB_MAX_ORDER>
      facmom = getIndexedMoment<DIM-1, cent>(a_face, a_facdatfab, a_graph, a_dx);
    //need to translate this into higher dim because that is what the matrix sees
    int facedir = 0;
    IndexedMoments<    DIM, HOEB_MAX_ORDER> volmom = getVolMomFromFacMom<DIM>(facmom, facedir);
    IndexedMoments<    DIM, HOEB_MAX_ORDER> bcmom;
    if(a_bcname == string("Dirichlet"))
    {
      bcmom = volmom;
    }
    else if(a_bcname == string("Neumann"))
    {
      for(MomentIterator<DIM, HOEB_MAX_ORDER> momit; momit.ok(); ++momit)
      {
        auto momind = momit();
        Real bcval = 0;
        if(momind[a_facedir] > 0)
        {
          int pd = momind[a_facedir];
          IndexTM<int, DIM> derind =momind - IndexTM<int, DIM>::BASISV(a_facedir);
          Real derval = volmom[derind];
          bcval = pd*derval;
        }
        bcmom[momind] = bcval;
      }
    }
    else
    {
      PROTO_ASSERT(false, "fillmatrixrow: bogus bc name");
    }
    shiftMomentAndFillRow(a_mat, bcmom, a_xbardistance, a_currentRow);
  }
  /******/  
  template <CENTERING cent>
  LAPACKMatrix
  getMomentMatrix(LAPACKMatrix                                        & a_weights,
                  const Neighborhood<cent>                            & a_block,
                  const shared_ptr< GeometryService<HOEB_MAX_ORDER> > & a_geoserv,
                  const Chombo4::Box                                  & a_srcDomain,
                  const RealVect                                      & a_xbar,
                  unsigned int                                          a_ibox,
                  Real                                                  a_dx)
  {
    typedef IndexedMoments<DIM  , HOEB_MAX_ORDER> IndMomDIM;
    typedef IndexedMoments<DIM-1, HOEB_MAX_ORDER> IndMomSDMinOne;
    unsigned int Pcols = IndMomDIM::size();
    unsigned int Nrows = a_block.size();
    LAPACKMatrix retval(Nrows, Pcols);
    a_weights.define(Nrows, Nrows);
    a_weights.setVal(0.);
    auto dit = a_geoserv->getDBL(a_srcDomain).dataIterator();
    
    const auto & graphsldptr =  a_geoserv->getGraphs(  a_srcDomain);
    const auto & volmomldptr =  a_geoserv->getVoluData(a_srcDomain);
    
    const auto & graph     = (*graphsldptr)[dit[a_ibox]];
    const auto & volmomfab = (*volmomldptr)[dit[a_ibox]];

    unsigned int currentRow = 0;
    {
      //volumes
      for(unsigned int ivol = 0; ivol < a_block.m_volumes.size(); ivol++)
      {
        EBIndex<cent> vof    = a_block.m_volumes[ivol];
        RealVect vofloc = hoeb::getVoFLocation(vof, a_dx);
        IndexTM<Real, DIM> distance = getIndexDistance(vofloc, a_xbar);
        IndexedMoments<DIM  , HOEB_MAX_ORDER>
          volmom = getIndexedMoment<DIM, cent>(vof, volmomfab, graph, a_dx);
        volmom.shift(distance);
        for(MomentIterator<DIM, HOEB_MAX_ORDER> momit; momit.ok(); ++momit)
        {
          auto momind = momit();
          unsigned int currentCol = IndMomDIM::indexOf(momind);
          retval(currentRow, currentCol) = volmom[momind];
        }
        a_weights(currentRow, currentRow) = a_block.m_volumesWgt[currentRow];
        currentRow++;
      }
    }
    {
      //now for the  domain x faces
      const auto & xfacldptr  =  a_geoserv->getXFaceData(a_srcDomain);
      const auto & xfacdatfab = (*xfacldptr)[dit[a_ibox]];
      for(unsigned int ivol = 0; ivol < a_block.m_xfacesDiri.size(); ivol++)
      {
        const auto & face = a_block.m_xfacesDiri[ivol];
        RealVect vofloc = getFaceLocation<XFACE>(face, a_dx, 0);
        IndexTM<Real, DIM> distance = getIndexDistance(vofloc, a_xbar);
        fillMatrixRowForDomBC<XFACE>(retval, xfacdatfab, face, graph, distance,
                                     currentRow, 0, a_dx, string("Dirichlet"));

        a_weights(currentRow,currentRow) = a_block.m_xfacesDiriWgt[ivol];
        currentRow++;
      }
      for(unsigned int ivol = 0; ivol < a_block.m_xfacesNeum.size(); ivol++)
      {
        const auto & face = a_block.m_xfacesNeum[ivol];
        RealVect vofloc = getFaceLocation<XFACE>(face, a_dx, 0);
        IndexTM<Real, DIM> distance = getIndexDistance(vofloc, a_xbar);
        fillMatrixRowForDomBC<XFACE>(retval, xfacdatfab, face, graph, distance,
                                     currentRow, 0, a_dx, string("Neumann"));

        a_weights(currentRow,currentRow) = a_block.m_xfacesNeumWgt[ivol];
        currentRow++;
      }
    }
    {
      //now for the  domain y faces
      const auto & yfacldptr  =  a_geoserv->getYFaceData(a_srcDomain);
      const auto & yfacdatfab = (*yfacldptr)[dit[a_ibox]];
      for(unsigned int ivol = 0; ivol < a_block.m_yfacesDiri.size(); ivol++)
      {
        const auto & face = a_block.m_yfacesDiri[ivol];
        RealVect vofloc = getFaceLocation<YFACE>(face, a_dx, 1);
        IndexTM<Real, DIM> distance = getIndexDistance(vofloc, a_xbar);
        fillMatrixRowForDomBC<YFACE>(retval, yfacdatfab, face, graph, distance,
                                     currentRow, 1, a_dx, string("Dirichlet"));
      
        a_weights(currentRow, currentRow) = a_block.m_yfacesDiriWgt[ivol];
        currentRow++;
      }
      for(unsigned int ivol = 0; ivol < a_block.m_yfacesNeum.size(); ivol++)
      {
        const auto & face = a_block.m_yfacesNeum[ivol];
        RealVect vofloc = getFaceLocation<YFACE>(face, a_dx, 1);
        IndexTM<Real, DIM> distance = getIndexDistance(vofloc, a_xbar);
        fillMatrixRowForDomBC<YFACE>(retval, yfacdatfab, face, graph, distance,
                                     currentRow, 1, a_dx, string("Neumann"));
      
        a_weights(currentRow, currentRow) = a_block.m_yfacesNeumWgt[ivol];
        currentRow++;
      }
    }
#if DIM==3
    {
      //now for the  domain z faces
      const auto & zfacldptr  =  a_geoserv->getYFaceData(a_srcDomain);
      const auto & zfacdatfab = (*yfacldptr)[dit[a_ibox]];
      for(unsigned int ivol = 0; ivol < a_block.m_zfacesDiri.size(); ivol++)
      {
        const auto & face = a_block.m_zFacesDiri[ivol];
        RealVect vofloc = getFaceLocation<ZFACE>(face, a_dx, 2);
        IndexTM<Real, DIM> distance = getIndexDistance(vofloc, a_xbar);
        fillMatrixRowForDomBC<ZFACE>(retval, zfacdatfab, face, graph, distance,
                                     currentRow, 2, a_dx, string("Dirichlet"));
      
        a_weights(currentRow,currentRow) = a_block.m_zfacesDiriWgt[ivol];
        currentRow++;
      }
      for(unsigned int ivol = 0; ivol < a_block.m_zfacesNeum.size(); ivol++)
      {
        const auto & face = a_block.m_zFacesNeum[ivol];
        RealVect vofloc = getFaceLocation<ZFACE>(face, a_dx, 2);
        IndexTM<Real, DIM> distance = getIndexDistance(vofloc, a_xbar);
        fillMatrixRowForDomBC<ZFACE>(retval, zfacdatfab, face, graph, distance,
                                     currentRow, 2, a_dx, string("Neumann"));
      
        a_weights(currentRow,currentRow) = a_block.m_zfacesNeumWgt[ivol];
        currentRow++;
      }
    }
#endif
    {
      //Cut faces
      for(unsigned int ivol = 0; ivol < a_block.m_cutfacesDiri.size(); ivol++)
      {
        const auto & face = a_block.m_cutfacesDiri[ivol];
        //-1 gets the cell centered location
        RealVect vofloc = getFaceLocation<BOUNDARY>(face, a_dx, -1);
        IndexTM<Real, DIM> distance = getIndexDistance(vofloc, a_xbar);
        fillMatrixRowForEBBC(retval, a_geoserv, a_srcDomain, face, graph, distance,
                             dit[a_ibox], currentRow, 2, a_dx, 
                             string("Dirichlet"));
      
        a_weights(currentRow,currentRow) = a_block.m_cutfacesDiriWgt[ivol];
        currentRow++;
      }
      for(unsigned int ivol = 0; ivol < a_block.m_cutfacesNeum.size(); ivol++)
      {
        const auto & face = a_block.m_cutfacesNeum[ivol];
        //-1 gets the cell centered location
        RealVect vofloc = getFaceLocation<BOUNDARY>(face, a_dx, -1);
        IndexTM<Real, DIM> distance = getIndexDistance(vofloc, a_xbar);
        fillMatrixRowForEBBC(retval, a_geoserv, a_srcDomain, face, graph, distance,
                             dit[a_ibox], currentRow, 2, a_dx, 
                             string("Neumann"));
      
        a_weights(currentRow,currentRow) = a_block.m_cutfacesNeumWgt[ivol];
        currentRow++;
      }
    }

    return retval;
  }  


  //Moore-Penrose inverse of the moment matrix is calulcated using weighted least squares
  inline LAPACKMatrix
  getMoorePenroseInverse(const LAPACKMatrix & a_M)
  {
    LAPACKMatrix MT = a_M;
    MT.transpose();

    LAPACKMatrix retval;
    multiply(retval, MT, a_M);
    retval.invert();
    return retval;
  }
  ///
  inline      LAPACKMatrix
  getStencilMatrix(const LAPACKMatrix & a_Qmat,
                   const LAPACKMatrix & a_Gmat,
                   const LAPACKMatrix & a_GTGinvMat)
  {
    LAPACKMatrix GTmat = a_Gmat;
    GTmat.transpose();
    //thi9s is (gtg)-1 GT
    LAPACKMatrix gtginvgtMat;
    multiply(gtginvgtMat, a_GTGinvMat, GTmat);
    
    LAPACKMatrix Smat;
    multiply(Smat, a_Qmat, gtginvgtMat);
    return Smat;
  }

   
  ///
  template <CENTERING cent>
  inline LAPACKMatrix
  getDharshiQMatrix(const EBIndex<cent>                                 & a_face,
                    const Neighborhood<CELL>                            & a_block,
                    const shared_ptr< GeometryService<HOEB_MAX_ORDER> > & a_geoserv,
                    const Chombo4::Box                                  & a_srcDomain,
                    const RealVect                                      & a_xbar,
                    unsigned int                                          a_ibox,
                    Real                                                  a_dx,
                    int                                                   a_facedir,
                    const Chombo4::DataIndex                            & a_datind)
  {
    IndexedMoments<DIM-1, HOEB_MAX_ORDER> 
      facmom = getIndexedMoment(a_face, a_geoserv, a_srcDomain, a_ibox, a_dx);
    IndexedMoments<    DIM, HOEB_MAX_ORDER> volmom = getVolMomFromFacMom<DIM>(facmom, a_facedir);

    RealVect facloc = hoeb::getFaceLocation<cent>(a_face, a_dx, a_facedir);
    IndexTM<Real, DIM> distance = getIndexDistance(facloc, a_xbar);
    volmom.shift(distance);

    LAPACKMatrix Qmat;


    for(MomentIterator<DIM, HOEB_MAX_ORDER> momit; momit.ok(); ++momit)
    {
      auto momind = momit();
      Real matval = 0;
      if(momind[a_facedir] > 0)
      {
        int pd = momind[a_facedir];
        IndexTM<int, DIM> derind =momind - IndexTM<int, DIM>::BASISV(a_facedir);
        Real derval = volmom[derind];
        matval = pd*derval;
      }
      unsigned int currentCol = IndexedMoments<DIM, HOEB_MAX_ORDER>::indexOf(momind);
      Qmat(0,currentCol) = matval;
    }
    
    return Qmat;
  }

  ///template specialization to get the eb flux q matrix
  template < >
  inline LAPACKMatrix
  getDharshiQMatrix(const EBIndex<BOUNDARY>                             & a_face,
                    const Neighborhood<CELL>                            & a_block,
                    const shared_ptr< GeometryService<HOEB_MAX_ORDER> > & a_geoserv,
                    const Chombo4::Box                                  & a_srcDomain,
                    const RealVect                                      & a_xbar,
                    unsigned int                                          a_ibox,
                    Real                                                  a_dx,
                    int                                                   a_facedir,
                    const Chombo4::DataIndex                            & a_datind)
  {
    LAPACKMatrix Qmat;
    Qmat.setVal(0.);
    const auto & graphsldptr = a_geoserv->getGraphs(    a_srcDomain);
    const auto & graph = (*graphsldptr)[a_datind];
    for(unsigned int normDir = 0; normDir < DIM; normDir++)
    {
      IndexedMoments<DIM, HOEB_MAX_ORDER> normMom=
        getEBNormalMoment(a_face, a_geoserv, a_srcDomain,
                          a_datind, graph, normDir, a_dx);
        
      for(MomentIterator<DIM, HOEB_MAX_ORDER> momit; momit.ok(); ++momit)
      {
        auto momind = momit();
        if(momind[normDir] > 0)
        {
          int pd = momind[normDir];
          IndexTM<int, DIM> derind =momind - IndexTM<int, DIM>::BASISV(normDir);
          Real derval = normMom[derind];
          unsigned int currentCol = IndexedMoments<DIM, HOEB_MAX_ORDER>::indexOf(momind);
          Qmat(0,currentCol) += pd*derval;
        }
      }
    }
    return Qmat;
  }
  
  ///
  inline 
  LocalStencil<CELL, Real>
  getStencilFromMatrix(const LAPACKMatrix      & a_stenmat,
                       const Neighborhood<CELL>& a_block)
  {
    LocalStencil<CELL, Real> retval;
    for(unsigned int ivol = 0; ivol < a_block.m_volumes.size(); ivol++)
    {
      const auto& vof = a_block.m_volumes[ivol];
      //this depends upon the volumes being the first bits
      Real wgt = a_stenmat(ivol, 0);
      retval.add(vof, wgt);
    }

    return retval;
  }
                       
  ///
  /**
     Return integral F dA  (do NOT divide out the area).
     For EB Fluxes, it is Integral F-n dA with an outward-facing normal.
  */
  template <CENTERING cent>
  LocalStencil<CELL, Real> 
  getDharshiIntFluxDAStencil(const EBIndex<cent>                                 & a_face,
                             const EBIndex<CELL>                                 & a_startingVoF,
                             const std::string                                     a_dombcname[2*DIM],
                             const std::string                                   & a_ebbcname,
                             const shared_ptr< GeometryService<HOEB_MAX_ORDER> > & a_geoserv,
                             const Chombo4::Box                                  & a_srcDomain,
                             unsigned int a_ibox,                              
                             Real a_alpha, Real a_beta, Real a_dx, int a_facedir,
                             Side::LoHiSide a_vofside)
  {
    //facedir == -1 means the EB
    const auto & graphsldptr = a_geoserv->getGraphs(    a_srcDomain);
    auto dit = a_geoserv->getDBL(a_srcDomain).dataIterator();
    const auto & datind = dit[a_ibox];
    const auto & graph = (*graphsldptr)[datind];
    LocalStencil<CELL, Real>  retval;
    //when I know the answer, short circuit all the other stuff and return it
    if(a_facedir==-1)
    {
      if(a_ebbcname  == string("Neumann"))
      {
        return retval;
      }
    }
    else if ((a_facedir >= 0) && a_face.m_isBoundary)
    {
      int iindex = ebp_index(a_facedir, a_vofside);
      if(a_dombcname[iindex] == string("Neumann"))
      {
        return retval;
      }
    }
    else
    {
      //get the neighborhood near the target face
      hoeb::Neighborhood<CELL>  block =
        hoeb::getCellNeighbors<cent>(a_face, graph, a_dombcname, a_ebbcname, a_facedir);
      const EBIndex<cent>& face = a_face;
      Real                 dx = a_dx;
      int                  facedir = a_facedir;
      RealVect xbar =   getFaceLocation(face, dx, facedir);

      //this matrix contains all the moments and boundary condition constraints
      LAPACKMatrix weight;
      LAPACKMatrix Mmat = getMomentMatrix<CELL>(weight, block, a_geoserv,
                                                a_srcDomain, xbar, a_ibox, a_dx);

      LAPACKMatrix Gmat; //G = WM
      multiply(Gmat, weight, Mmat);
      
      //Moore-Penrose inverse of the moment matrix
      //is calulcated using weighted least squares
      //because G is weighted version of M
      LAPACKMatrix GTGinvMat = getMoorePenroseInverse(Gmat);

      //this is the matrix that defines the operator
      LAPACKMatrix Qmat =
        getDharshiQMatrix<cent>
        (a_face, block, a_geoserv, a_srcDomain, xbar, a_ibox, a_dx, a_facedir, datind);

      LAPACKMatrix Smat = getStencilMatrix(Qmat, Gmat, GTGinvMat);
      
      retval = hoeb::getStencilFromMatrix(Smat, block);
    }
    return retval;
  }

}
#endif
