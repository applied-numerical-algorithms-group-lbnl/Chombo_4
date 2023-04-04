#ifndef __hoeb_basics__
#define __hoeb_basics__
#include "Proto.H"
#include "Chombo_EBChombo.H"
#include "Chombo_EBLevelBoxData.H"
#include "Chombo_GeometryService.H"
#include "Chombo_EBDictionary.H"

#define HOEB_MAX_ORDER 4

namespace hoeb_basics
{


  ///get location of face center
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
  IndexTM<int, D>
  getVolIndexFromFace(const  IndexTM<int, D-1>& a_facind, int a_facedir)
  {
    IndexTM<int, D> retval;
    int ifac = 0;
    for(int ivol = 0; ivol < D; ivol++)
    {
      if(ivol != a_facedir)
      {
        retval[ivol] = a_facind[ifac];
        ifac++;
      }
      else
      {
        retval[ivol] = 0;
      }
    }
    return retval;
  }
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
  /******/
  template <unsigned int D, CENTERING cent>
  inline IndexedMoments<D, HOEB_MAX_ORDER>   
  getFaceIndexedMoment(const EBIndex<cent>                                              &   a_index,
                       const HostIrregData<cent, IndexedMoments<D, HOEB_MAX_ORDER>, 1>  &   a_irregfab,
                       const EBGraph                                                    &   a_graph,
                       const Real                                                       &   a_dx,
                       int a_facedir)
  {
    IndexedMoments<D, HOEB_MAX_ORDER> retval;
    auto dom = a_graph.getDomain();
    Point regCheckPt = a_index.m_pt;
    if(!dom.contains(regCheckPt))
    {
      regCheckPt -= Point::Basis(a_facedir);
      PR_assert(dom.contains(regCheckPt));
    }
    
    if(a_graph.isRegular(regCheckPt))
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
  /******/

  inline IndexedMoments<DIM, HOEB_MAX_ORDER>   
  getEBIndexedMoment(const EBIndex<BOUNDARY>                                          &   a_index,
                     const HostIrregData<CELL, IndexedMoments<DIM, HOEB_MAX_ORDER>, 1>  &   a_irregfab,
                     const EBGraph                                                    &   a_graph,
                     const Real                                                       &   a_dx)
  {
    IndexedMoments<DIM, HOEB_MAX_ORDER> retval;
    
    if(a_graph.isRegular(a_index.m_pt))
    {
      retval.setToZero();
    }
    else if(a_graph.isIrregular(a_index.m_pt))
    {
      EBIndex<CELL> vof = a_index.getVoF();
      retval = a_irregfab(vof, 0);
    }
    else
    {
      PROTO_ASSERT(false, "should not get here by construction");
    }
    return retval;
  }
  template <CENTERING cent>
  inline IndexedMoments<DIM-1, HOEB_MAX_ORDER>   
  getIndexedMomentFace(const EBIndex<cent>                                 & a_index,
                       const shared_ptr< GeometryService<HOEB_MAX_ORDER> > & a_geoserv,
                       const Chombo4::Box                                  & a_srcDomain,
                       const Chombo4::DataIndex                            & a_datInd,
                       Real                                                  a_dx)
  {
    IndexedMoments<DIM-1, HOEB_MAX_ORDER>   retval;
    Chombo4::MayDay::Error("did not find template specialization");
    return retval;
  }
  ///template specialization for x
  template < >
  inline IndexedMoments<DIM-1, HOEB_MAX_ORDER>   
  getIndexedMomentFace(const EBIndex<XFACE>                                & a_index,
                       const shared_ptr< GeometryService<HOEB_MAX_ORDER> > & a_geoserv,
                       const Chombo4::Box                                  & a_srcDomain,
                       const Chombo4::DataIndex                            & a_datInd,
                       Real                                                  a_dx)
  {
    const auto & xfacldptr  =  a_geoserv->getXFaceData(a_srcDomain);
    const auto & xfacdatfab = (*xfacldptr)[a_datInd];
    const auto & graphsldptr =  a_geoserv->getGraphs(  a_srcDomain);
    const auto & graph     = (*graphsldptr)[a_datInd];
    return getFaceIndexedMoment<DIM-1,XFACE>(a_index, xfacdatfab, graph, a_dx, 0);
  }

  ///template specialization for y
  template < >
  inline IndexedMoments<DIM-1, HOEB_MAX_ORDER>   
  getIndexedMomentFace(const EBIndex<YFACE>                                &   a_index,
                       const shared_ptr< GeometryService<HOEB_MAX_ORDER> > & a_geoserv,
                       const Chombo4::Box                                  & a_srcDomain,
                       const Chombo4::DataIndex                            & a_datInd,
                       Real                                                  a_dx)
  {
    const auto & yfacldptr  =  a_geoserv->getYFaceData(a_srcDomain);
    const auto & yfacdatfab = (*yfacldptr)[a_datInd];
    const auto & graphsldptr =  a_geoserv->getGraphs(  a_srcDomain);
    const auto & graph     = (*graphsldptr)[a_datInd];
    return getFaceIndexedMoment<DIM-1, YFACE>(a_index, yfacdatfab, graph, a_dx, 1);
  }
#if DIM==3  
  ///template specialization for z
  template < >
  inline IndexedMoments<DIM-1, HOEB_MAX_ORDER>   
  getIndexedMomentFace(const EBIndex<ZFACE>                                & a_index,
                       const shared_ptr< GeometryService<HOEB_MAX_ORDER> > & a_geoserv,
                       const Chombo4::Box                                  & a_srcDomain,
                       const Chombo4::DataIndex                            & a_datInd,
                       Real                                                  a_dx)

  {
    const auto & zfacldptr  =  a_geoserv->getZFaceData(a_srcDomain);
    const auto & zfacdatfab = (*zfacldptr)[a_datInd];
    const auto & graphsldptr =  a_geoserv->getGraphs(  a_srcDomain);
    const auto & graph     = (*graphsldptr)[a_datInd];
    return getFaceIndexedMoment<DIM-1, ZFACE>(a_index, zfacdatfab, graph, a_dx, 2);
  }
#endif
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
    if(a_normDir == 0)
    {
      const auto   xnormMomLDptr = a_geoserv->getEBNormalData_x(a_srcDomain);
      const auto & xnormmomfab = (*xnormMomLDptr)[a_datind];
      retval = getEBIndexedMoment(a_face, xnormmomfab, a_graph, a_dx);
    }
    else if(a_normDir == 1)
    {
      const auto   ynormMomLDptr = a_geoserv->getEBNormalData_y(a_srcDomain);
      const auto & ynormmomfab = (*ynormMomLDptr)[a_datind];
      retval = getEBIndexedMoment(a_face, ynormmomfab, a_graph, a_dx);
    }
#if DIM==3
    else if(a_normDir == 2)
    {
      const auto   znormMomLDptr = a_geoserv->getEBNormalData_z(a_srcDomain);
      const auto & znormmomfab = (*znormMomLDptr)[a_datind];
      retval = getEBIndexedMoment(a_face, znormmomfab, a_graph, a_dx);
    }
#endif
    else
    {
      PROTO_ASSERT(false, "getshiftedindexedmoment: bogus normdir");
    }
    return retval;
  }

}

#endif
