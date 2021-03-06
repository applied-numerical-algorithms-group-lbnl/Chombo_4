#include "Hoeb_Utilities.H"
#include "Hoeb_ExactSolutions.H"
#include "Chombo_GeometryService.H"
#include "Chombo_ParmParse.H"

namespace hoeb
{
  
  /******/  
  LocalStencil<CELL, Real> 
  getFullDharshiStencil(const EBIndex<CELL>                                 & a_vof,
                        const std::string                                     a_dombcname[2*DIM],
                        const std::string                                   & a_ebbcname,
                        const shared_ptr< GeometryService<HOEB_MAX_ORDER> > & a_geoserv,
                        Proto::Box                                          & a_srcDomain,
                        unsigned int a_ibox, Real a_alpha, Real a_beta, Real a_dx) 
  {
    LocalStencil<CELL, Real> vofsten;
    //use the age-old trick of building up the stencil by construction
    auto dbl = a_geoserv->getDBL(a_srcDomain);
    auto dit = dbl.dataIterator();
    const auto & graphsldptr = a_geoserv->getGraphs(    a_srcDomain);
    const auto & graph = (*graphsldptr)[dit[a_ibox]];
    for(SideIterator sit; sit.ok(); ++sit)
    {
      int isign = sign(sit());
      {                                                // 
        auto xfaces = graph.getXFaces(a_vof, sit());
        for(unsigned int iface = 0; iface < xfaces.size(); iface++)
        {
          auto face = xfaces[iface];
          LocalStencil<CELL, Real> fluxsten =
            getDharshiIntFluxDAStencil<XFACE>(face,
                                              a_vof,a_dombcname,a_ebbcname,
                                              a_geoserv, a_srcDomain, a_ibox,
                                              a_alpha, a_beta, a_dx, 0, sit());
          fluxsten *= Real(isign);
          vofsten += fluxsten;
        }
      }
      {                                                // 
        auto yfaces = graph.getYFaces(a_vof, sit());
        for(unsigned int iface = 0; iface < yfaces.size(); iface++)
        {
          auto face = yfaces[iface];
          LocalStencil<CELL, Real> fluxsten =
            getDharshiIntFluxDAStencil<YFACE>(face,
                                              a_vof,a_dombcname,a_ebbcname,
                                              a_geoserv, a_srcDomain, a_ibox,
                                              a_alpha, a_beta, a_dx, 1, sit());
          fluxsten *= Real(isign);
          vofsten += fluxsten;
        }
      }
#if DIM==3      
      {                                                // 
        auto zfaces = graph.getZFaces(a_vof, sit());
        for(unsigned int iface = 0; iface < zfaces.size(); iface++)
        {
          auto face = zfaces[iface];
          LocalStencil<CELL, Real> fluxsten =
            getDharshiIntFluxDAStencil<ZFACE>(face,
                                              a_vof, a_dombcname, a_ebbcname,
                                              a_geoserv, a_srcDomain, a_ibox,
                                              a_alpha, a_beta, a_dx, 2, sit());
          fluxsten *= Real(isign);
          vofsten += fluxsten;
        }
      }
#endif
    }
    {
      EBIndex<BOUNDARY> face = a_vof.getCutFace();
      LocalStencil<CELL, Real>
        fluxsten =
        getDharshiIntFluxDAStencil<BOUNDARY>(face,
                                             a_vof,a_dombcname,a_ebbcname,
                                             a_geoserv, a_srcDomain, a_ibox,
                                             a_alpha, a_beta, a_dx, -1, Side::Invalid);
      vofsten += fluxsten;
    }
    //need to divide by dx^d to get kappa*lapl(phi)
    Real cellvolinv = 1;
    for(int idir = 0; idir < DIM; idir++)
    {
      cellvolinv /= a_dx;
    }
    vofsten *= cellvolinv;
    
    return vofsten;
  }

  /******/
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
                     const Chombo4::Box                                  & a_domain,
                     const Real                                          & a_dx,
                     unsigned int                                          a_ibox)

  {
    /* just in case a fit of madness causes me to use geometric multigrid */
    a_needDiagonalWeights = true;
    a_srcDomain = ProtoCh::getProtoBox(a_domain);
    a_dstDomain = ProtoCh::getProtoBox(a_domain);
    Chombo4::ParmParse pp;
    
    Real alpha = 1.0;
    Real beta = -0.001;
    string dombcname;
    int nghost;
    pp.get("num_ghost_cells", nghost);
    pp.get("domainBC"  , dombcname);
    pp.get("EBBC"      , a_ebbcName);
    pp.get("alpha"     , alpha);
    pp.get("beta"      , beta);
    pout() << "domainBC"  << " = " <<  dombcname      << endl;
    pout() << "EBBC"      << " = " <<  a_ebbcName     << endl;
    
    a_stencilName = string("Dharshi_Laplacian");

    const auto & graphsldptr = a_geoserv->getGraphs(    a_srcDomain);
    const auto & voldatldptr = a_geoserv->getVoluData(  a_srcDomain);
    const auto & ebfdatldptr = a_geoserv->getEBFaceData(a_srcDomain);
    const auto & xfadatldptr = a_geoserv->getXFaceData( a_srcDomain);
    const auto & yfadatldptr = a_geoserv->getYFaceData( a_srcDomain);
    const auto & zfadatldptr = a_geoserv->getZFaceData( a_srcDomain);
    const auto & dbl         = a_geoserv->getDBL(a_srcDomain);
    auto dit = dbl.dataIterator();

    const auto & graph = (*graphsldptr)[dit[a_ibox]];
    a_srcValid  = graph.validBox();
    a_dstValid  = graph.validBox();
    a_srcGhost  = Point::Ones(nghost);
    a_dstGhost  = Point::Ones(nghost);
    a_needDiagonalWeights = true;
    
    string dombcarray[2*DIM];
    for(int ivec = 0; ivec < 2*DIM; ivec++)
    {
      dombcarray[ivec] = dombcname;
    }
    for(auto bit = a_dstValid.begin(); bit != a_dstValid.end(); ++bit)
    {
      auto vofs = graph.getVoFs(*bit);
      for(unsigned int ivof = 0; ivof < vofs.size(); ivof++)
      {
        LocalStencil<CELL, Real> vofsten = getFullDharshiStencil(vofs[ivof],
                                                                 dombcarray, a_ebbcName,
                                                                 a_geoserv, a_srcDomain, a_ibox,
                                                                 alpha, beta, a_dx);
        a_dstVoFs.push_back(vofs[ivof]);
        a_stencil.push_back(vofsten);
      }
    }
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
