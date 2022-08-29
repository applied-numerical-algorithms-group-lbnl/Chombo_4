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
  typedef CH4_Data_Choreography::DistributedData<EBGraph>   graph_distrib_t;
  typedef GraphConstructorFactory<EBHostData<CELL, Real, 1> > hostfactorycell_t;
  typedef CH4_Data_Choreography::DistributedData<EBHostData<CELL, Real, 1> > cell_distrib_t;
  
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

  inline Real
  devendranLphiInhomogeneous(const EBIndex<CELL>                                                &  a_vof,
                             const Chombo4::DataIndex                                           &  a_datInd,
                             const cell_distrib_t                                               &  a_hostphi,
                             const shared_ptr<graph_distrib_t >                                 &  a_graphs,
                             const Chombo4::DisjointBoxLayout                                   &  a_grids,
                             const Chombo4::Box                                                 &  a_domain,
                             const Real                                                         &  a_dx,
                             const shared_ptr<EBDictionary<HOEB_MAX_ORDER, Real, CELL, CELL> >  &  a_dictionary,
                             const shared_ptr< GeometryService<HOEB_MAX_ORDER> >                &  a_geoserv)
  {
    std::string  dombcname[2*DIM];
    std::string  ebbcname, dombcone;
    ParmParse pp;
    pp.get("domainBC"  , dombcone);
    pp.get("EBBC"      , ebbcname);
    for(int ivec  = 0; ivec < 2*DIM; ivec++)
    {
      dombcname[ivec] = dombcone;
    }

    bool dividebyarea = false; //we want int flux da here
    
    Real intFdotNDA = 0;
    const auto & graphsldptr =  a_geoserv->getGraphs(  a_domain);
    const auto & volmomldptr =  a_geoserv->getVoluData(a_domain);
    const auto & graph     = (*graphsldptr)[a_datInd];
    const auto & voldat    = (*volmomldptr)[a_datInd];
    const auto & hostphi   =      a_hostphi[a_datInd];
    for (SideIterator sit; sit.ok(); ++sit)
    {
      int isign = sign(sit());
      //begin debug
      int stophere = 0;
//      int hipt = a_domain.bigEnd()[0];
      int mdpt = a_domain.bigEnd()[1]/2;
      bool printStuff = false;
//      if((a_vof.m_pt[0] == hipt) && (a_vof.m_pt[1]== mdpt) && (isign == 1))
      if((a_vof.m_pt[0] == mdpt) && (a_vof.m_pt[1]== mdpt) && (isign == 1))
      {
        stophere = 1;
        printStuff = true;
      }
      //end debug
      {
        vector< EBIndex<XFACE> > xfaces = graph.getXFaces(a_vof, sit());
        for (int iface = 0; iface < xfaces.size(); iface++)
        {
          const auto& face = xfaces[iface];
          int facedir = 0;
          Real fluxpt =
            hoeb::getDevendranFluxFace<XFACE>(hostphi, graph,
                                              face, a_vof,  dombcname, ebbcname, a_geoserv,
                                              a_domain, a_dx, facedir, sit(),
                                              a_datInd, dividebyarea, printStuff);
          intFdotNDA += isign*fluxpt;
        }
      }
      {
        vector< EBIndex<YFACE> > yfaces = graph.getYFaces(a_vof, sit());
        for (int iface = 0; iface < yfaces.size(); iface++)
        {
          const auto& face = yfaces[iface];
          int facedir = 1;
          Real fluxpt =
            hoeb::getDevendranFluxFace<YFACE>(hostphi, graph,
                                              face, a_vof,  dombcname, ebbcname, a_geoserv,
                                              a_domain, a_dx, facedir, sit(),
                                              a_datInd, dividebyarea, false);
          intFdotNDA += isign*fluxpt;
        }
      }
#if DIM==3      
      {
        vector< EBIndex<ZFACE> > zfaces = graph.getZFaces(a_vof, sit());
        for (int iface = 0; iface < zfaces.size(); iface++)
        {
          const auto& face = zfaces[iface];
          int facedir = 2;
          Real fluxpt =
            hoeb::getDevendranFluxFace<ZFACE>(hostphi, graph,
                                              face, a_vof,  dombcname, ebbcname, a_geoserv,
                                              a_domain,  a_dx, facedir, sit(),
                                              a_datInd, dividebyarea, false);
          intFdotNDA += isign*fluxpt;
        }
      }
#endif
    }
    if(graph.isIrregular(a_vof.m_pt))
    {
      EBIndex<BOUNDARY> face = a_vof.getCutFace();
      int facedir = -1;
      Real fluxpt =
        hoeb::getDevendranFluxFace<BOUNDARY>(hostphi, graph,
                                             face, a_vof,  dombcname, ebbcname, a_geoserv,
                                             a_domain,  a_dx, facedir, Side::Lo,
                                             a_datInd, dividebyarea, false);
      intFdotNDA += fluxpt;
    }
      
    //now need to divde by the volume
    IndexedMoments<DIM, HOEB_MAX_ORDER>
      volmom = hoeb_basics::getIndexedMoment<DIM,CELL>(a_vof, voldat, graph, a_dx);
    Real volume = volmom[IndexTM<int, DIM>::Zero];
    Real retval  = 0;
    if(volume > 1.0e-16)
    {
      retval = intFdotNDA/volume;
    }
      
    return retval;

  }
  

  Real
  schwartzLphiInhomogeneous(const EBIndex<CELL>                                                &  a_vof,
                            const Chombo4::DataIndex                                           &  a_datInd,
                            const cell_distrib_t                                               &  a_hostphi,
                            const shared_ptr<graph_distrib_t>                                  &  a_graphs,
                            const Chombo4::DisjointBoxLayout                                   &  a_grids,
                            const Chombo4::Box                                                 &  a_domain,
                            const Real                                                         &  a_dx,
                            const shared_ptr<EBDictionary<HOEB_MAX_ORDER, Real, CELL, CELL> >  &  a_dictionary,
                            const shared_ptr< GeometryService<HOEB_MAX_ORDER> >                &  a_geoserv);
  
  ///
  inline void
  getDevendranKLPInhomogeneous(EBLevelBoxData<CELL, 1>                                            &  a_klp,
                               const EBLevelBoxData<CELL, 1>                                      &  a_phi,
                               const shared_ptr<graph_distrib_t>                                  &  a_graphs,
                               const Chombo4::DisjointBoxLayout                                   &  a_grids,
                               const Chombo4::Box                                                 &  a_domain,
                               const Real                                                         &  a_dx,
                               const shared_ptr<EBDictionary<HOEB_MAX_ORDER, Real, CELL, CELL> >  &  a_dictionary,
                               const shared_ptr< GeometryService<HOEB_MAX_ORDER> >                &  a_geoserv)
  {
    typedef GraphConstructorFactory<EBHostData<CELL, Real, 1> > hostfactorycell_t;
    auto ghost = a_phi.ghostVect();
    cell_distrib_t    hostphi(a_grids,  ghost, hostfactorycell_t(a_graphs));
    cell_distrib_t    hostklp(a_grids,  ghost, hostfactorycell_t(a_graphs));

    EBLevelBoxData<CELL, 1>::copyToHost(hostphi, a_phi);

    Chombo4::DataIterator dit = a_grids.dataIterator();
    for(int ibox = 0; ibox < dit.size(); ibox++)
    {

      Bx grbx = ProtoCh::getProtoBox(a_grids[dit()]);
      auto graph = (*a_graphs)[dit[ibox]];
      auto& hostfab = hostklp[dit[ibox]];
      hostfab.setVal(0.);
      for(Proto::BoxIterator bit = grbx.begin(); bit != grbx.end(); ++bit)
      {
        auto pt   = *bit;
        if(!graph.isCovered(pt))
        {
          auto vofs  = graph.getVoFs(pt);

          for(int ivof = 0; ivof < vofs.size(); ivof++)
          {
            Real klphval =hoeb::
              devendranLphiInhomogeneous(vofs[ivof],
                                         dit[ibox],
                                         hostphi,
                                         a_graphs,
                                         a_grids,
                                         a_domain,
                                         a_dx,
                                         a_dictionary,
                                         a_geoserv);
             hostfab(vofs[ivof], 0) = klphval;
          }
        }
      }
    }
    EBLevelBoxData<CELL, 1>::copyToDevice(a_klp, hostklp);
  }
  ///
  inline void
  getSchwartzKLPInhomogeneous(EBLevelBoxData<CELL, 1>                                            &  a_klp,
                              const EBLevelBoxData<CELL, 1>                                      &  a_phi,
                              const shared_ptr<graph_distrib_t>                                  &  a_graphs,
                              const Chombo4::DisjointBoxLayout                                   &  a_grids,
                              const Chombo4::Box                                                 &  a_domain,
                              const Real                                                         &  a_dx,
                              const shared_ptr<EBDictionary<HOEB_MAX_ORDER, Real, CELL, CELL> >  &  a_dictionary,
                              const shared_ptr< GeometryService<HOEB_MAX_ORDER> >                &  a_geoserv)
  {
    typedef GraphConstructorFactory<EBHostData<CELL, Real, 1> > hostfactorycell_t;
    auto ghost = a_phi.ghostVect();
    cell_distrib_t    hostphi(a_grids,  ghost, hostfactorycell_t(a_graphs));
    cell_distrib_t    hostklp(a_grids,  ghost, hostfactorycell_t(a_graphs));
    EBLevelBoxData<CELL, 1>::copyToHost(hostphi, a_phi);

    Chombo4::DataIterator dit = a_grids.dataIterator();
    for(int ibox = 0; ibox < dit.size(); ibox++)
    {

      Bx grbx = ProtoCh::getProtoBox(a_grids[dit()]);
      hostklp[dit[ibox]].setVal(0.);
      auto graph = (*a_graphs)[dit[ibox]];
      for(Proto::BoxIterator bit = grbx.begin(); bit != grbx.end(); ++bit)
      {
        auto pt   = *bit;
        if(!graph.isCovered(pt))
        {
          auto graph = (*a_graphs)[dit[ibox]];
          auto vofs  = graph.getVoFs(pt);

          for(int ivof = 0; ivof < vofs.size(); ivof++)
          {
            Real klphval = schwartzLphiInhomogeneous(vofs[ivof],
                                                     dit[ibox],
                                                     hostphi,
                                                     a_graphs,
                                                     a_grids,
                                                     a_domain,
                                                     a_dx,
                                                     a_dictionary,
                                                     a_geoserv);

            hostklp[dit[ibox]](vofs[ivof], 0) = klphval;
          }
        }
      }
    }
    EBLevelBoxData<CELL, 1>::copyToDevice(a_klp, hostklp);
  }
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

  ///get the stencil for the laplacian described in Schwartz, et al.
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
          const shared_ptr<graph_distrib_t    >                  &  a_graphs,
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
                    const shared_ptr<graph_distrib_t>                                 &  a_graphsFine,
                    const Chombo4::DisjointBoxLayout                                  &  a_gridsFine,
                    const Chombo4::Box                                                &  a_domFine,
                    const Real                                                        &  a_dxFine,
                    const shared_ptr<graph_distrib_t>                                 &  a_graphsCoar,
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
              const shared_ptr<graph_distrib_t    >                             &  a_graphsFine,
              const Chombo4::DisjointBoxLayout                                  &  a_gridsFine,
              const Chombo4::Box                                                &  a_domFine,
              const Real                                                        &  a_dxFine,
              const shared_ptr<graph_distrib_t    >                             &  a_graphsCoar,
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
              const shared_ptr<graph_distrib_t    >                             &  a_graphsFine,
              const Chombo4::DisjointBoxLayout                                  &  a_gridsFine,
              const Chombo4::Box                                                &  a_domFine,
              const Real                                                        &  a_dxFine,
              const shared_ptr<graph_distrib_t    >                             &  a_graphsCoar,
              const Chombo4::DisjointBoxLayout                                  &  a_gridsCoar,
              const Chombo4::Box                                                &  a_domCoar,
              const Real                                                        &  a_dxCoar,
              const shared_ptr<EBDictionary<HOEB_MAX_ORDER, Real, cent, cent> > &  a_dictionary,
               const shared_ptr< GeometryService<HOEB_MAX_ORDER> >               &  a_geoserv);
  
  
}
#include "Hoeb_UtilitiesImplem.H"
#endif
