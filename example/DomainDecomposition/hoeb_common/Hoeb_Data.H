#ifndef __hoeb_data__
#define __hoeb_data__
#include "Proto.H"
#include "Chombo_EBChombo.H"
#include "Chombo_EBLevelBoxData.H"
#include "Chombo_GeometryService.H"
#include "Chombo_EBDictionary.H"
#include "Hoeb_LAPACKMatrix.H"
#include "Hoeb_Neighborhood.H"


namespace hoeb_data
{

  ///
  /**
     It is a runtime error to send this function a covered index.m_pt.
     If the input index.m_pt is regular, set to regular moments. 
     If the input index.m_pt is irregular, extract the moments from the data holder.
  */
  template <unsigned int D, CENTERING cent, unsigned int GEOMETRY_ORDER>
  inline IndexedMoments<D, GEOMETRY_ORDER>   
  getIndexedMoment(const EBIndex<cent>                                              &   a_index,
                   const HostIrregData<cent, IndexedMoments<D, GEOMETRY_ORDER>, 1>  &   a_irregfab,
                   const EBGraph                                                    &   a_graph,
                   const Real                                                       &   a_dx)
  {
    IndexedMoments<D, GEOMETRY_ORDER> retval;
    
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
  
  
  /*************/
  template <unsigned int GEOMETRY_ORDER>
  inline IndexedMoments<DIM, GEOMETRY_ORDER>
  getEBIndexedMoment(const EBIndex<BOUNDARY>                             & a_face,
                     const shared_ptr< GeometryService<GEOMETRY_ORDER> > & a_geoserv,
                     const Chombo4::Box                                  & a_srcDomain,
                     const Chombo4::DataIndex                            & a_datind,
                     const EBGraph                                       & a_graph,
                     Real a_dx)
  {
    IndexedMoments<DIM, GEOMETRY_ORDER> retval;
    //have to do this because the ebnormal moments have a cell type
    const auto   ebMomLDptr = a_geoserv->getEBFaceData(a_srcDomain);
    const auto & ebmomfab = (*ebMomLDptr)[a_datind];
    retval = getIndexedMoment<DIM, BOUNDARY, GEOMETRY_ORDER>
      (a_face, ebmomfab, a_graph, a_dx);

    return retval;
  }
  /*************/
  template <unsigned int GEOMETRY_ORDER>
  inline IndexedMoments<DIM, GEOMETRY_ORDER>
  getEBNormalMoment(const EBIndex<BOUNDARY>                             & a_face,
                    const shared_ptr< GeometryService<GEOMETRY_ORDER> > & a_geoserv,
                    const Chombo4::Box                                  & a_srcDomain,
                    const Chombo4::DataIndex                            & a_datind,
                    const EBGraph                                       & a_graph,
                    unsigned int a_normDir, Real a_dx)
  {
    IndexedMoments<DIM, GEOMETRY_ORDER> retval;
    //have to do this because the ebnormal moments have a cell type
    EBIndex<CELL> vof = a_face.getVoF();
    if(a_normDir == 0)
    {
      const auto   xnormMomLDptr = a_geoserv->getEBNormalData_x(a_srcDomain);
      const auto & xnormmomfab = (*xnormMomLDptr)[a_datind];
      retval = getIndexedMoment<DIM, CELL, GEOMETRY_ORDER>(vof, xnormmomfab, a_graph, a_dx);
    }
    else if(a_normDir == 1)
    {
      const auto   ynormMomLDptr = a_geoserv->getEBNormalData_y(a_srcDomain);
      const auto & ynormmomfab = (*ynormMomLDptr)[a_datind];
      retval = getIndexedMoment<DIM, CELL, GEOMETRY_ORDER>(vof, ynormmomfab, a_graph, a_dx);
    }
#if DIM==3
    else if(a_normDir == 2)
    {
      const auto   znormMomLDptr = a_geoserv->getEBNormalData_z(a_srcDomain);
      const auto & znormmomfab = (*znormMomLDptr)[a_datind];
      retval = getIndexedMoment<DIM,CELL, GEOMETRY_ORDER>(vof, znormmomfab, a_graph, a_dx);
    }
#endif
    else
    {
      PROTO_ASSERT(false, "getshiftedindexedmoment: bogus normdir");
    }
    return retval;
  }
  

}
///

#endif
