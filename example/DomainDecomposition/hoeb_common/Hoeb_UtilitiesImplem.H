#ifndef __hoeb_utilities_implem_
#define __hoeb_utilities_implem_

namespace hoeb
{
  ///
  /**
     Return integral F dA  (do NOT divide out the area).
     For EB Fluxes, it is Integral F-n dA with an outward-facing normal.
  */
  template <CENTERING cent>
  LocalStencil<CELL, Real> 
  getDharshiIntFluxDAStencil(const EBIndex<cent>                                 & a_face,
                             const EBIndex<CELL>                                 & a_startingVoF,
                             const std::string                                   & a_dombcname,
                             const std::string                                   & a_ebbcname,
                             const shared_ptr< GeometryService<HOEB_MAX_ORDER> > & a_geoserv,
                             Proto::Box                                          & a_srcDomain,
                             unsigned int a_ibox,                              
                             Real a_alpha, Real a_beta, Real a_dx, int a_facedir)
  {
    //facedir == -1 means the EB
    const auto & graphsldptr = a_geoserv->getGraphs(    a_srcDomain);
    auto dit = a_geoserv->getDBL(a_srcDomain).dataIterator();
    const auto & graph = (*graphsldptr)[dit[a_ibox]];
    LocalStencil<CELL, Real>  retval;
    //first let us get the cases where we don't have to do anything.
    if(((a_facedir==-1)         && (a_ebbcname  == string("Neumann")) ||
        ((a_face.m_isBoundary)  && (a_dombcname == string("Neumann")))))
    {
      return retval;
    }
    else
    {
      vector<EBIndex<CELL> > neighborhood = hoeb::getVoFNeighbors(a_startingVoF, graph);
      
      LAPACKMatrix momMat;
    }
    return retval;
  }

}
#endif
