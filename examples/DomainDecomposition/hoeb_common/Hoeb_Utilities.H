#ifndef __hoeb_utilities__
#define __hoeb_utilities__
#include "Proto.H"
#include "Chombo_EBChombo.H"
#include "Chombo_EBLevelBoxData.H"
#include "Chombo_GeometryService.H"
#include "Chombo_EBDictionary.H"
#include "Hoeb_LAPACKMatrix.H"
#include "Hoeb_Neighborhood.H"
#include "Hoeb_CompositeStencil.H"
#define HOEB_MAX_ORDER 4

namespace hoeb
{
  /// get the dreaded M from the neighborhood
  template <CENTERING cent>
  LAPACKMatrix
  getMomentMatrix(const Neighborhood<cent>                            & a_block,
                  const shared_ptr< GeometryService<HOEB_MAX_ORDER> > & a_geoserv,
                  const Chombo4::Box                                  & a_srcDomain,
                  const RealVect                                      & a_xbar,
                  const Chombo4::DataIndex                            & a_dit,
                  Real                                                  a_dx);

  template <CENTERING cent>
  Real
  getDevendranFluxFace(const EBHostData<CELL, Real, 1>                     & a_phi,
                       const EBGraph                                       & a_graph,
                       const EBIndex<cent>                                 & a_face,
                       const EBIndex<CELL>                                 & a_startingVoF,
                       const std::string                                     a_dombcname[2*DIM],
                       const std::string                                   & a_ebbcname,
                       const shared_ptr< GeometryService<HOEB_MAX_ORDER> > & a_geoserv,
                       const Chombo4::Box                                  & a_domain,
                       Real a_dx, int a_facedir,
                       Side::LoHiSide a_vofside,
                       const Chombo4::DataIndex& a_dit,
                       bool a_divideByArea, bool a_printStuff);  //true if you want the average flux.  false if you want int f dA
  
  /// outputs warning message if a_errcode != 0
  void checkError(int a_errcode, string a_prefix);

  Real
  getKappaLphiVoF(const EBIndex<CELL>                                                &  a_vof,
                  const Chombo4::DataIndex                                           &  a_datInd,
                  const EBLevelBoxData<CELL, 1>                                      &  a_phi,
                  const shared_ptr<LevelData<EBGraph> >                              &  a_graphs,
                  const Chombo4::DisjointBoxLayout                                   &  a_grids,
                  const Chombo4::Box                                                 &  a_domain,
                  const Real                                                         &  a_dx,
                  const shared_ptr<EBDictionary<HOEB_MAX_ORDER, Real, CELL, CELL> >  &  a_dictionary,
                  const shared_ptr< GeometryService<HOEB_MAX_ORDER> >                &  a_geoserv,
                  string                                                                a_stencilname);
  ///
  /**
     Return integral F dA  (do NOT divide out the area).
     For EB Fluxes, it is Integral F-n dA with an outward-facing normal.
  */
  
  template <CENTERING cent>
  inline CompositeStencil
  getDevendranIntFluxDACompositeStencil(const EBIndex<cent>                                 & a_face,
                                        const EBIndex<CELL>                                 & a_startingVoF,
                                        const std::string                                     a_dombcname[2*DIM],
                                        const std::string                                   & a_ebbcname,
                                        const shared_ptr< GeometryService<HOEB_MAX_ORDER> > & a_geoserv,
                                        const Chombo4::Box                                  & a_srcDomain,
                                        const Chombo4::DataIndex                            & a_dit,
                                        Real a_dx, int a_facedir,
                                        Side::LoHiSide a_vofside, bool a_printStuff);

  ///
  /**
     Return integral F dA  (do NOT divide out the area).
     For EB Fluxes, it is Integral F-n dA with an outward-facing normal.
     this only returns the bit that acts on cell data.   The other stuff is assumed to be zero
     (homogeneous operator--for use with EBStencil).
  */
  
  template <CENTERING cent>
  inline LocalStencil<CELL, Real> 
  getDharshiIntFluxDAStencil(const EBIndex<cent>                                 & a_face,
                             const EBIndex<CELL>                                 & a_startingVoF,
                             const std::string                                     a_dombcname[2*DIM],
                             const std::string                                   & a_ebbcname,
                             const shared_ptr< GeometryService<HOEB_MAX_ORDER> > & a_geoserv,
                             const Chombo4::Box                                  & a_srcDomain,
                             const Chombo4::DataIndex                            & a_dit,
                             Real a_dx, int a_facedir,
                             Side::LoHiSide a_vofside, bool a_printStuff)
  {
    CompositeStencil composite
      = getDevendranIntFluxDACompositeStencil(a_face,     
                                              a_startingVoF, 
                                              a_dombcname,
                                              a_ebbcname, 
                                              a_geoserv,  
                                              a_srcDomain,
                                              a_dit, a_dx, a_facedir, a_vofside, a_printStuff);
    return composite.m_cellSten;
  }



  ///get the stencil for the laplacian described in Devendran, et al.
  void
  getHomogeneousDharshiStencil(string                                              & a_stencilName,
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
                               const Chombo4::Box                                  & a_domain,
                               const Real                                          & a_dx,
                               unsigned int                                          a_ibox);

  ///get the stencil for the laplacian described in Devendran, et al.
  void
  schwartzLaplStencil(string                                              & a_stencilName,
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
                      const Chombo4::Box                                  & a_domain,
                      const Real                                          & a_dx,
                      unsigned int                                          a_ibox);
  
  shared_ptr<hoeb::BaseExactSolution<HOEB_MAX_ORDER> >   getBaseExactSoltuion();
  
  ///return an implicit function based on stuff in the input file
  shared_ptr<BaseIF> getImplicitFunction();

  ///
  /**
     Fill in scalar from exactsolutions class
  */
  void
  fillPhi(EBLevelBoxData<CELL, 1>                                &  a_phi,
          const shared_ptr<LevelData<EBGraph> >                  &  a_graphs,
          const Chombo4::DisjointBoxLayout                       &  a_grids,
          const Chombo4::Box                                     &  a_domFine,
          const Real                                             &  a_dx,
          const shared_ptr< GeometryService<HOEB_MAX_ORDER> >    &  a_geoserv);

  ///
  /**
   **/
  void
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

  ///
  /**
   **/
  void
  restrictPhi(EBLevelBoxData<CELL, 1>                                           &  a_phiFToC,
              const EBLevelBoxData<CELL, 1>                                     &  a_phiFine,
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

  ///
  /**
   **/
  template<CENTERING cent>
  void
  restrictFlux(EBLevelBoxData<cent, 1>                                          &  a_phiFToC,
              const EBLevelBoxData<cent, 1>                                     &  a_phiFine,
              const shared_ptr<LevelData<EBGraph> >                             &  a_graphsFine,
              const Chombo4::DisjointBoxLayout                                  &  a_gridsFine,
              const Chombo4::Box                                                &  a_domFine,
              const Real                                                        &  a_dxFine,
              const shared_ptr<LevelData<EBGraph> >                             &  a_graphsCoar,
              const Chombo4::DisjointBoxLayout                                  &  a_gridsCoar,
              const Chombo4::Box                                                &  a_domCoar,
              const Real                                                        &  a_dxCoar,
              const shared_ptr<EBDictionary<HOEB_MAX_ORDER, Real, cent, cent> > &  a_dictionary,
               const shared_ptr< GeometryService<HOEB_MAX_ORDER> >               &  a_geoserv);
  
  
}
#include "Hoeb_UtilitiesImplem.H"
#endif
