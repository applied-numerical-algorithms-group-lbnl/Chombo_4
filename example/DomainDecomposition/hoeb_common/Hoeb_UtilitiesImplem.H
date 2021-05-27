#ifndef __hoeb_utilities_implem_
#define __hoeb_utilities_implem_
#include "Hoeb_Neighborhood.H"
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
                             const std::string                                    a_dombcname[2*DIM],
                             const std::string                                   & a_ebbcname,
                             const shared_ptr< GeometryService<HOEB_MAX_ORDER> > & a_geoserv,
                             Proto::Box                                          & a_srcDomain,
                             unsigned int a_ibox,                              
                             Real a_alpha, Real a_beta, Real a_dx, int a_facedir,
                             Side::LoHiSide a_vofside)
  {
    //facedir == -1 means the EB
    const auto & graphsldptr = a_geoserv->getGraphs(    a_srcDomain);
    auto dit = a_geoserv->getDBL(a_srcDomain).dataIterator();
    const auto & graph = (*graphsldptr)[dit[a_ibox]];
    LocalStencil<CELL, Real>  retval;

    //first let us get the cases where we don't have to do anything.
    if(((a_facedir==-1)         && (a_ebbcname  == string("Neumann"))))
    {
      return retval;
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
      Neighborhood  block =
        getNeighbors<cent>(a_face, graph, a_dombcname, a_ebbcname, a_facedir);
      //this matrix contains all the moments and boundary condition constraints
      LAPACKMatrix momMat;
    }
    return retval;
  }

}
#endif
