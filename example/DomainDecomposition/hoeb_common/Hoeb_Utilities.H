#ifndef __hoeb_utilities__
#define __hoeb_utilities__
#include "Proto.H"
#include "Chombo_EBChombo.H"
#include "Chombo_EBLevelBoxData.H"
#include "Chombo_GeometryService.H"
#include "Chombo_EBDictionary.H"
#include "Hoeb_LAPACKMatrix.H"
#define HOEB_MAX_ORDER 4

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
                             Real a_alpha, Real a_beta, Real a_dx, int a_facedir)
  {
    //facedir == -1 means the EB
    LocalStencil<CELL, Real>  retval;

    return retval;
  }

  //get the stencil for the laplacian described in Devendran, et al.
  void
  dharshiLaplStencil(string                                              & a_stencilName,
                     string                                              & a_ebbcName,
                     vector<EBIndex<CELL> >                              & a_dstVoFs,
                     vector<LocalStencil<CELL, Real> >                   & a_stencil,
                     Proto::Box                                          & a_srcValid,
                     Proto::Box                                          & a_dstValid,
                     Proto::Box                                          & a_srcDomain,
                     Proto::Box                                          & a_dstDomain,
                     Proto::Point                                        & a_srcGhost,
                     Proto::Point                                        & a_dstGhost,
                     bool                                                & a_needDiagonalWeights,
                     const shared_ptr< GeometryService<HOEB_MAX_ORDER> > & a_geoserv,
                     const Chombo4::DisjointBoxLayout                    & a_grids,
                     const Chombo4::Box                                  & a_domFine,
                     const Real                                          & a_dx,
                     unsigned int                                          a_ibox);

  //return an implicit function based on stuff in the input file
  extern shared_ptr<BaseIF> getImplicitFunction();

  ///
  /**
     Fill in scalar from exactsolutions class
  */
  extern void
  fillPhi(EBLevelBoxData<CELL, 1>                                &  a_phi,
          const shared_ptr<LevelData<EBGraph> >                  &  a_graphs,
          const Chombo4::DisjointBoxLayout                       &  a_grids,
          const Chombo4::Box                                     &  a_domFine,
          const Real                                             &  a_dx,
          const shared_ptr< GeometryService<HOEB_MAX_ORDER> >    &  a_geoserv);

  /****/
  extern void
  restrictKappaLphi(EBLevelBoxData<CELL, 1>                                           &  a_klpFToC,
                    const EBLevelBoxData<CELL, 1>                                     &  a_klpFine,
                    const shared_ptr<LevelData<EBGraph> >                             &  a_graphsFine,
                    const Chombo4::DisjointBoxLayout                                  &  a_gridsFine,
                    const Chombo4::Box                                                &  a_domFine,
                    const Real                                                        &  a_dxFine,
                    const shared_ptr<LevelData<EBGraph> >                             &  a_graphsCoar,
                    const Chombo4::DisjointBoxLayout                                  &  a_gridsCoar,
                    const Chombo4::Box                                                &  a_domCoar,
                    const Real                                                        &  a_dxCoar,
                    const shared_ptr<EBDictionary<HOEB_MAX_ORDER, Real, CELL, CELL> > &  a_dictionary,
                    const shared_ptr< GeometryService<HOEB_MAX_ORDER> >               &  a_geoserv);
}
#endif
