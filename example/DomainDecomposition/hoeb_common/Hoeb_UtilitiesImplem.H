#ifndef __hoeb_utilities_implem_
#define __hoeb_utilities_implem_
#include "Hoeb_Neighborhood.H"
namespace hoeb
{
  /******/
  inline RealVect
  getVoFLocation(const EBIndex<CELL>& a_face,
                 Real                 a_dx)
  {
    RealVect retval;
    auto pt = a_face.m_pt; //high size of the face
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
                   const HostIrregData<cent, IndexedMoments<D, HOEB_MAX_ORDER>, 1> &   a_irregfab,
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
  template <CENTERING cent>
  inline RealVect
  getFaceLocation(const EBIndex<cent>& a_face,
                  Real                 a_dx,
                  int                  a_facedir)
  {
    RealVect retval;
    auto pt = a_face.m_pt; //high size of the face
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

  ///
  template < CENTERING cent>
  void
  fillMatrixRowForDomBC
  (LAPACKMatrix                                                         & a_mat,
   const HostIrregData<cent, IndexedMoments<DIM-1, HOEB_MAX_ORDER>, 1 > & a_facdatfab,
   const EBIndex<cent>                                                  & a_face,
   const EBGraph                                                        & a_graph,
   const IndexTM<Real, DIM>                                             & a_xbardistance,
   const unsigned int                                                     a_currentRow,
   int                                                                    a_facedir,
   Real                                                                   a_dx,
   string                                                                 a_bcname)
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
    //shift so it has the same xbar
    bcmom.shift(a_xbardistance);
    for(MomentIterator<DIM, HOEB_MAX_ORDER> momit; momit.ok(); ++momit)
    {
      auto momind = momit();
      unsigned int currentCol = IndexedMoments<DIM, HOEB_MAX_ORDER>::indexOf(momind);
      a_mat(a_currentRow, currentCol) = bcmom[momind];
    }

  }
  /******/  
  template <CENTERING cent>
  LAPACKMatrix
  getMomentMatrix(const Neighborhood<cent>                            & a_block,
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

        currentRow++;
      }
      for(unsigned int ivol = 0; ivol < a_block.m_xfacesNeum.size(); ivol++)
      {
        const auto & face = a_block.m_xfacesNeum[ivol];
        RealVect vofloc = getFaceLocation<XFACE>(face, a_dx, 0);
        IndexTM<Real, DIM> distance = getIndexDistance(vofloc, a_xbar);
        fillMatrixRowForDomBC<XFACE>(retval, xfacdatfab, face, graph, distance,
                                     currentRow, 0, a_dx, string("Neumann"));

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
      
        currentRow++;
      }
      for(unsigned int ivol = 0; ivol < a_block.m_yfacesNeum.size(); ivol++)
      {
        const auto & face = a_block.m_yfacesNeum[ivol];
        RealVect vofloc = getFaceLocation<YFACE>(face, a_dx, 1);
        IndexTM<Real, DIM> distance = getIndexDistance(vofloc, a_xbar);
        fillMatrixRowForDomBC<YFACE>(retval, yfacdatfab, face, graph, distance,
                                     currentRow, 1, a_dx, string("Neumann"));
      
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
      
        currentRow++;
      }
      for(unsigned int ivol = 0; ivol < a_block.m_zfacesNeum.size(); ivol++)
      {
        const auto & face = a_block.m_zFacesNeum[ivol];
        RealVect vofloc = getFaceLocation<ZFACE>(face, a_dx, 2);
        IndexTM<Real, DIM> distance = getIndexDistance(vofloc, a_xbar);
        fillMatrixRowForDomBC<ZFACE>(retval, zfacdatfab, face, graph, distance,
                                     currentRow, 2, a_dx, string("Neumann"));
      
        currentRow++;
      }
    }
#endif    
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
    const auto & graph = (*graphsldptr)[dit[a_ibox]];
    LocalStencil<CELL, Real>  retval;

    if(a_facedir==-1)
    {
      if(a_ebbcname  == string("Neumann"))
      {
        return retval;
      }
      else
      {
        PROTO_ASSERT(false, "not implemented");
      }
          
    }
    else if ((a_facedir >= 0) && a_face.m_isBoundary)
    {
      int iindex = ebp_index(a_facedir, a_vofside);
      if(a_dombcname[iindex] == string("Neumann"))
      {
        return retval;
      }
      else
      {
        PROTO_ASSERT(false, "not implemented");
      }
    }
    else
    {
      Neighborhood<CELL>  block =
        getCellNeighbors<cent>(a_face, graph, a_dombcname, a_ebbcname, a_facedir);
      const EBIndex<cent>& face = a_face;
      Real                 dx = a_dx;
      int                  facedir = a_facedir;
      RealVect xbar =   getFaceLocation(face, dx, facedir);

      //this matrix contains all the moments and boundary condition constraints
      LAPACKMatrix momMat = getMomentMatrix<CELL>(block, a_geoserv, a_srcDomain, xbar, a_ibox, a_dx);
    }
    return retval;
  }

}
#endif
