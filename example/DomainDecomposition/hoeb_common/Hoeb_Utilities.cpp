#include "Hoeb_Utilities.H"
#include "Hoeb_ExactSolutions.H"
#include "Chombo_GeometryService.H"
#include "Chombo_ParmParse.H"

namespace hoeb
{
  /******/  
  void
  dharshiLaplStencil(string                                 & a_stencilName,
                     string                                 & a_ebbcName,
                     vector<EBIndex<CELL> >                 & a_dstVoFs,
                     vector<LocalStencil<CELL, Real> >      & a_stencil,
                     Stencil<Real>                          & a_regStencilInterior,
                     Proto::Box                             & a_regApplyBox,
                     Proto::Box                             & a_srcValid,
                     Proto::Box                             & a_dstValid,
                     Proto::Box                             & a_srcDomain,
                     Proto::Box                             & a_dstDomain,
                     Proto::Point                           & a_srcGhost,
                     Proto::Point                           & a_dstGhost,
                     bool                                   & a_irregOnly,
                     bool                                   & a_needDiagonalWeights,
                     const shared_ptr<LevelData<EBGraph> >  &  a_graphs,
                     const Chombo4::DisjointBoxLayout       &  a_grids,
                     const Chombo4::Box                     &  a_domFine,
                     const Real                             &  a_dx)

  {
    Chombo4::ParmParse pp;
    
    Real alpha = 1.0;
    Real beta = -0.001;
    int dombc = 1;
    int ebbc  = 1;
    pp.get("domainBC"  , dombc);
    pp.get("EBBC"      , ebbc);
    pp.get("alpha"     , alpha);
    pp.get("beta"      , beta);
    pout() << "domainBC"  << " = " <<  dombc      << endl;
    pout() << "EBBC"      << " = " <<  ebbc       << endl;
  
    a_needDiagonalWeights = true;
    a_irregOnly           = true;
  }

  shared_ptr<BaseIF> getImplicitFunction()
  {
    RealVect center = 0.5*RealVect::Unit();
    Real radius = 0.1;
    bool inside = false;
    ParmParse pp;
    std::vector<Real> centvec;
    pp.get("radius", radius);
    pp.get("inside", inside);
    pp.getarr("center", centvec, 0, DIM);
    for(int idir = 0; idir < DIM; idir++)
    {
      center[idir] = centvec[idir];
    }
    SimpleSphereIF* sphereptr = new SimpleSphereIF(center, radius, inside);
    shared_ptr<BaseIF> retval(static_cast<BaseIF*>(sphereptr));
    return retval;
  }

  /****/
  void
  fillPhi(EBLevelBoxData<CELL, 1>                                &  a_phi,
          const shared_ptr<LevelData<EBGraph> >                  &  a_graphs,
          const Chombo4::DisjointBoxLayout                       &  a_grids,
          const Chombo4::Box                                     &  a_domain,
          const Real                                             &  a_dx,
          const shared_ptr< GeometryService<HOEB_MAX_ORDER> >    &  a_geoserv)
  {
    typedef IndexedMoments<DIM  , HOEB_MAX_ORDER> IndMomDIM;
    typedef HostIrregData<CELL,      IndMomDIM , 1>  VoluData;
    typedef GraphConstructorFactory<EBHostData<CELL, Real, 1> > hostfactory_t;
    
    RealVect center = 0.5*RealVect::Unit();
    Real radius = 0.1;
    bool inside = false;
    ParmParse pp;
    std::vector<Real> centvec;
    pp.get("radius", radius);
    pp.get("inside", inside);
    pp.getarr("center", centvec, 0, DIM);
    for(int idir = 0; idir < DIM; idir++)
    {
      center[idir] = centvec[idir];
    }
    hoeb::SineSphereEF<HOEB_MAX_ORDER> phigen(radius, center);

    LevelData<EBHostData<CELL, Real, 1> > hostdata
      (a_grids, 1,  a_phi.ghostVect(), hostfactory_t(a_graphs));
    
    auto voldatldptr = a_geoserv->getVoluData(a_domain);
    const auto& voldatld = *voldatldptr;
    Chombo4::DataIterator dit = a_grids.dataIterator();
    for(unsigned int ibox = 0; ibox < dit.size(); ibox++)
    {
      const auto& graph  = (*a_graphs)[dit[ibox]];
      const auto& voldat = voldatld[dit[ibox]];
      Bx validbx = graph.validBox();
      for(auto bit = validbx.begin(); bit != validbx.end(); ++bit)
      {
        auto vofs= graph.getVoFs(*bit);
        for(int ivof = 0; ivof < vofs.size(); ++ivof)
        {
          auto& vof = vofs[ivof];
          Real phival = phigen(graph, a_dx, voldat, vof);
        }
      }
    }

    EBLevelBoxData<CELL, 1>::copyToDevice(a_phi, hostdata);
  }
  /****/
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
                    const shared_ptr< GeometryService<HOEB_MAX_ORDER> >               &  a_geoserv)
  {
    std::string  nobcname        = string("no_bcs");
    std::string  restrictionName = string("Multigrid_Restriction");
    string dombc[2*DIM];
    for(unsigned int idom = 0; idom < 2*DIM;  idom++)
    {
      dombc[idom] = nobcname;
    }
    a_dictionary->registerStencil(restrictionName, dombc, nobcname, a_domFine, a_domCoar, false);
    Chombo4::DataIterator dit = a_gridsCoar.dataIterator();
    for(unsigned int ibox = 0; ibox < dit.size(); ++ibox)
    {
      auto      & coarfab = a_klpFToC[dit[ibox]];
      const auto& finefab = a_klpFine[dit[ibox]];
      auto stencil = a_dictionary->getEBStencil(restrictionName, nobcname, a_domFine, a_domCoar, ibox);
      //set resc = Ave(resf) (true is initToZero)
      stencil->apply(coarfab, finefab,  true, 1.0);
    }
  }
} //end namespace hoeb
