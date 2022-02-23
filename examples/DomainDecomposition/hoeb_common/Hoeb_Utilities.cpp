#include "Hoeb_Utilities.H"
#include "Hoeb_ExactSolutions.H"
#include "implem/Proto_SecondOrderCell.H"
#include "Chombo_GeometryService.H"
#include "Chombo_ParmParse.H"

namespace hoeb
{
  
  Real
  devendranLphiInhomogeneous(const EBIndex<CELL>                                                &  a_vof,
                             const Chombo4::DataIndex                                           &  a_datInd,
                             const LevelData<EBHostData<CELL, Real, 1> >                        &  a_hostphi,
                             const shared_ptr<LevelData<EBGraph> >                              &  a_graphs,
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
        for (int iface = 0; iface < yfaces.size(); iface++)
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
                            const LevelData<EBHostData<CELL, Real, 1> >                        &  a_hostphi,
                            const shared_ptr<LevelData<EBGraph> >                              &  a_graphs,
                            const Chombo4::DisjointBoxLayout                                   &  a_grids,
                            const Chombo4::Box                                                 &  a_domain,
                            const Real                                                         &  a_dx,
                            const shared_ptr<EBDictionary<HOEB_MAX_ORDER, Real, CELL, CELL> >  &  a_dictionary,
                            const shared_ptr< GeometryService<HOEB_MAX_ORDER> >                &  a_geoserv)
  {
        
    const auto & graphsldptr = a_geoserv->getGraphs(    a_domain);
    const auto & voldatldptr = a_geoserv->getVoluData(  a_domain);
    const auto & ebfdatldptr = a_geoserv->getEBFaceData(a_domain);
    const auto & xfadatldptr = a_geoserv->getXFaceData( a_domain);
    const auto & yfadatldptr = a_geoserv->getYFaceData( a_domain);
    const auto & zfadatldptr = a_geoserv->getZFaceData( a_domain);
    const auto & dbl         = a_geoserv->getDBL(       a_domain);


    typedef IndexedMoments<DIM  , HOEB_MAX_ORDER> IndMomDIM;
    typedef IndexedMoments<DIM-1, HOEB_MAX_ORDER> IndMomSDMinOne;
    typedef MomentIterator<DIM  , HOEB_MAX_ORDER> MomItDIM;
    typedef MomentIterator<DIM-1, HOEB_MAX_ORDER> MomItSDMinOne;
    
    const auto & graph       = (*graphsldptr)[a_datInd];
    const auto & voludata    = (*voldatldptr)[a_datInd];
    const auto & ebfadata    = (*ebfdatldptr)[a_datInd];
    const auto & xfacdata    = (*xfadatldptr)[a_datInd];
    const auto & yfacdata    = (*yfadatldptr)[a_datInd];
    const auto & zfacdata    = (*zfadatldptr)[a_datInd];

    Proto::Poisson2ndOrder<XFACE, HOEB_MAX_ORDER> xflux;
    Proto::Poisson2ndOrder<YFACE, HOEB_MAX_ORDER> yflux;
    Proto::Poisson2ndOrder<ZFACE, HOEB_MAX_ORDER> zflux;
    ParmParse pp;
    string   ebbcname;
    Real alpha, beta;
    pp.get("EBBC"      , ebbcname);
    pp.get("alpha"     , alpha);
    pp.get("beta"      , beta);
    string dombcname("Neumann");
    string dombcarray[2*DIM];
    for(int ivec = 0; ivec < 2*DIM; ivec++)
    {
      dombcarray[ivec] = dombcname;
    }
        
    LocalStencil<CELL, Real> vofsten;
    Proto::SecondOrderStencil<HOEB_MAX_ORDER>::
      get2ndOrderDivFStencil(vofsten,
                             a_vof,
                             graph,
                             voludata,
                             ebfadata,
                             xfacdata, yfacdata, zfacdata,
                             xflux, yflux, zflux, dombcarray, ebbcname, 
                             a_dx);
    
    Real divfval = 0;
    for(unsigned int isten = 0; isten < vofsten.size(); isten++)
    {
      const auto& entry = vofsten.m_entries[isten];
      const auto& weight = entry.m_weight;
      const auto& vof    = entry.m_vof;
      Real phival = a_hostphi[a_datInd](vof, 0);
      divfval += weight*phival;
    }
    if(ebbcname == string("Dirichlet"))
    {
      shared_ptr<hoeb::BaseExactSolution<HOEB_MAX_ORDER> > exact = getBaseExactSoltuion();
      EBIndex<BOUNDARY> cutface = a_vof.getCutFace();
      //-1 is for the EB
      RealVect ebloc = hoeb_basics::getFaceLocation<BOUNDARY>(cutface, a_dx, -1);
      //true is to divide by the area
      Real phival = exact->evaluateDIM<BOUNDARY>(graph, a_dx, ebfadata, cutface, true);
      Real inhomogcontrib = vofsten.m_ebbcWeight*phival;
      divfval += inhomogcontrib;
    }
    else
    {
      Chombo4::MayDay::Warning("schwartz neumann bit not implemented yet");
    }
    
    Real phistart = a_hostphi[a_datInd](a_vof, 0);
    Real retval   = alpha*phistart + beta*divfval;
    return retval;
  }
  /*******/
  ///inhomogeneous 
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
                  string                                                                a_stencilName)
  {
    IntVect ghost = a_phi.ghostVect();
    typedef GraphConstructorFactory<EBHostData<CELL, Real, 1> > hostfactorycell_t;
    LevelData<EBHostData<CELL, Real, 1> >    hostphi(a_grids, 1, ghost, hostfactorycell_t(a_graphs));
    EBLevelBoxData<CELL,   1>::copyToHost(hostphi, a_phi);

    Real retval = 0;
    if(a_stencilName == string("Schwartz"))
    {
      retval =   schwartzLphiInhomogeneous(a_vof,
                                           a_datInd,
                                           hostphi,
                                           a_graphs,
                                           a_grids,
                                           a_domain,
                                           a_dx,
                                           a_dictionary,
                                           a_geoserv);
    }
    else if(a_stencilName == string("Devendran"))
    {
      retval =   devendranLphiInhomogeneous(a_vof,
                                           a_datInd,
                                           hostphi,
                                           a_graphs,
                                           a_grids,
                                           a_domain,
                                           a_dx,
                                           a_dictionary,
                                            a_geoserv);
    }
    else
    {
      Chombo4::MayDay::Error("bogus stencil name");
    }
    return retval;
  }
  /******/
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
                      unsigned int                                          a_ibox)

  {
    /* geometric multigrid is not a great idea here*/
    using Chombo4::pout;
    a_needDiagonalWeights = false;
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
    
    a_stencilName = string("Schwartz_Laplacian");

    const auto & graphsldptr = a_geoserv->getGraphs(    a_srcDomain);
    const auto & voldatldptr = a_geoserv->getVoluData(  a_srcDomain);
    const auto & ebfdatldptr = a_geoserv->getEBFaceData(a_srcDomain);
    const auto & xfadatldptr = a_geoserv->getXFaceData( a_srcDomain);
    const auto & yfadatldptr = a_geoserv->getYFaceData( a_srcDomain);
    const auto & zfadatldptr = a_geoserv->getZFaceData( a_srcDomain);
    const auto & dbl         = a_geoserv->getDBL(a_srcDomain);
    auto dit = dbl.dataIterator();


    typedef IndexedMoments<DIM  , HOEB_MAX_ORDER> IndMomDIM;
    typedef IndexedMoments<DIM-1, HOEB_MAX_ORDER> IndMomSDMinOne;
    typedef MomentIterator<DIM  , HOEB_MAX_ORDER> MomItDIM;
    typedef MomentIterator<DIM-1, HOEB_MAX_ORDER> MomItSDMinOne;
    
    const auto & graph       = (*graphsldptr)[dit[a_ibox]];
    auto & voludata    = (*voldatldptr)[dit[a_ibox]];
    const auto & ebfadata    = (*ebfdatldptr)[dit[a_ibox]];
    const auto & xfacdata    = (*xfadatldptr)[dit[a_ibox]];
    const auto & yfacdata    = (*yfadatldptr)[dit[a_ibox]];
    const auto & zfacdata    = (*zfadatldptr)[dit[a_ibox]];
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
    string dombc = dombcname;
    string ebbc = a_ebbcName;
    for(auto bit = a_dstValid.begin(); bit != a_dstValid.end(); ++bit)
    {
      auto vofs = graph.getVoFs(*bit);
      for(unsigned int ivof = 0; ivof < vofs.size(); ivof++)
      {
        LocalStencil<CELL, Real> vofsten;

        Proto::Poisson2ndOrder<XFACE, HOEB_MAX_ORDER> xflux;
        Proto::Poisson2ndOrder<YFACE, HOEB_MAX_ORDER> yflux;
        Proto::Poisson2ndOrder<ZFACE, HOEB_MAX_ORDER> zflux;
        
        Proto::SecondOrderStencil<HOEB_MAX_ORDER>::
          get2ndOrderDivFStencil(vofsten,
                                 vofs[ivof],
                                 graph,
                                 voludata,
                                 ebfadata,
                                 xfacdata, yfacdata, zfacdata,
                                 xflux, yflux, zflux, dombcarray, ebbc, 
                                 a_dx);
        
        a_dstVoFs.push_back(vofs[ivof]);
        a_stencil.push_back(vofsten);
      }
    }
  }
  
  /******/  
  void checkError(int a_errcode, string a_prefix)
  {
    if(a_errcode != 0)
    {
      Chombo4::pout() << "warning: " << a_prefix << " output error code = "  << a_errcode << endl;
    }
  }
  /******/  
  LocalStencil<CELL, Real> 
  getHomogeneousDharshiStencil(const EBIndex<CELL>                                 & a_vof,
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
                                              a_geoserv, a_srcDomain, dit[a_ibox],
                                              a_dx, 0, sit(), false);
          //begin debug
          Point debpt = 15*Point::Ones();
          if(a_vof.m_pt == debpt)
          {
            Chombo4::pout() << "face = " << face.m_pt << endl;
            fluxsten.print();
          }
          //end debug
          
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
                                              a_geoserv, a_srcDomain, dit[a_ibox],
                                              a_dx, 1, sit(), false);
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
                                              a_geoserv, a_srcDomain, dit[a_ibox],
                                              a_dx, 2, sit(), false);
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
                                             a_geoserv, a_srcDomain, dit[a_ibox],
                                             a_dx, -1, Side::Invalid, false);
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
                               unsigned int                                          a_ibox)

  {
    /* geometric multigrid is not a great idea here*/
    a_needDiagonalWeights = false;
    a_srcDomain = ProtoCh::getProtoBox(a_domain);
    a_dstDomain = ProtoCh::getProtoBox(a_domain);
    Chombo4::ParmParse pp;
    using Chombo4::pout;
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
        LocalStencil<CELL, Real> vofsten =
          getHomogeneousDharshiStencil(vofs[ivof],
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
    using Proto::BaseIF;
    shared_ptr<BaseIF>  retval;
    ParmParse pp;
    string which_geom;
    pp.get("which_geom", which_geom);
    using Chombo4::pout;
    if(which_geom == string("sphere"))
    {
      RealVect center = 0.5*RealVect::Unit();
      Real radius = 0.1;
      bool inside = false;
      std::vector<Real> centvec;
      pp.get("radius", radius);
      pp.get("inside", inside);
      pp.getarr("center", centvec, 0, DIM);
      for(int idir = 0; idir < DIM; idir++)
      {
        center[idir] = centvec[idir];
      }
      Proto::SimpleSphereIF* sphereptr = new Proto::SimpleSphereIF(center, radius, inside);
      pout() << "sphere implicit function with radius = " << radius << ", center = " << center << ", and inside = " << inside << endl;
      retval = shared_ptr<BaseIF>(static_cast<BaseIF*>(sphereptr));
    }
    else if(which_geom == string("all_regular"))
    {
      pout() << "all regular geometry" << endl;
      retval = shared_ptr<BaseIF>(new Proto::AllRegularIF());
    }
    else if(which_geom == string("plane"))
    {
      using Proto::PlaneIF;
      pout() << "plane implicit function" << endl;
      RealVect normal, startPt;
      vector<double> v_norm, v_start;
      pp.getarr("geom_normal", v_norm, 0, DIM);
      pp.getarr("geom_start_pt", v_start, 0, DIM);
      for(int idir = 0; idir < DIM; idir++)
      {
        normal[ idir] = v_norm[ idir];
        startPt[idir] = v_start[idir];
        pout() << "normal ["<< idir << "] = " << normal [idir]  << endl;
        pout() << "startPt["<< idir << "] = " << startPt[idir]  << endl;
      }
      retval = shared_ptr<BaseIF>(new PlaneIF(startPt, normal));
    }
    else
    {
      Chombo4::MayDay::Error("bogus geometry");
    }
    return retval;
  }

  /****/
  shared_ptr<hoeb::BaseExactSolution<HOEB_MAX_ORDER> > 
  getBaseExactSoltuion()
  {
    string whichphi;
    //this gets called all the time so only print out the diagnostics once
    static bool printedonce = false;
    using Chombo4::pout;
    ParmParse pp;
    pp.get("which_phi", whichphi);
    shared_ptr<hoeb::BaseExactSolution<HOEB_MAX_ORDER> >  retval;
    if(whichphi == string("SineSphere"))
    {
      RealVect center = 0.5*RealVect::Unit();
      Real radius = 0.1;
      bool inside = false;
      std::vector<Real> centvec;
      pp.get("radius", radius);
      pp.get("inside", inside);
      pp.getarr("center", centvec, 0, DIM);
      for(int idir = 0; idir < DIM; idir++)
      {
        center[idir] = centvec[idir];
      }
      if(!printedonce)
        pout() << "exact solution = sinesphere with radius = " << radius << ", center = " << center << endl;
      retval =  shared_ptr<hoeb::BaseExactSolution<HOEB_MAX_ORDER> >( new hoeb::SineSphereEF<HOEB_MAX_ORDER>(radius, center));;
    }
    else if(whichphi == string("Polynomial"))
    {
      using std::vector;
      using std::pair;
      vector<pair<Point, Real> > entries;
      int num_terms;
      pp.get("polynomial_num_terms", num_terms);
      entries.resize(num_terms);
      if(!printedonce)
        pout() << "exact solution = polynomial" << endl;
      for(int iterm = 0; iterm < num_terms; iterm++)
      {
        string coef_str   = string("polynomial_coefficient_") + to_string(iterm);
        string powers_str = string("polynomial_powers_")      + to_string(iterm);
        vector<int> entry_powers;
        Real        entry_coef;
        pp.getarr(powers_str.c_str(), entry_powers, 0, DIM);
        pp.get(     coef_str.c_str(), entry_coef);
        pair<Point, Real> entry;
        entry.second = entry_coef;
        for(int idir = 0; idir < DIM; idir++)
        {
          entry.first[idir] = entry_powers[idir];
        }
        entries.push_back(entry);
        if(!printedonce)
          pout() << "iterm = " << iterm << ", power = " << entry.first << ", weight = " << entry.second << endl;
      }
      hoeb::PolynomialEF<HOEB_MAX_ORDER>* pointerpoly = (new hoeb::PolynomialEF<HOEB_MAX_ORDER>(entries));
      BaseExactSolution<HOEB_MAX_ORDER> * pointerbase = static_cast<BaseExactSolution<HOEB_MAX_ORDER>* >(pointerpoly);
      
      retval = shared_ptr<hoeb::BaseExactSolution<HOEB_MAX_ORDER> > (pointerbase);
    }
    else
    {
      Chombo4::MayDay::Error("string for whichphi unrecognized");
    }

    printedonce = true;
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
    
    shared_ptr<hoeb::BaseExactSolution<HOEB_MAX_ORDER> > phigen =
      getBaseExactSoltuion();
    
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
          //true is to divide by the area
          Real phival = (*phigen)(graph, a_dx, voldat, vof, true);
          hostdata[dit[ibox]](vof, 0) = phival;
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

  /****/
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
                   const shared_ptr< GeometryService<HOEB_MAX_ORDER> >               &  a_geoserv)
  {
    std::string  nobcname        = string("no_bcs");
    std::string  restrictionName = string("Rho_Restriction");
    string dombc[2*DIM];
    for(unsigned int idom = 0; idom < 2*DIM;  idom++)
    {
      dombc[idom] = nobcname;
    }
    a_dictionary->registerStencil(restrictionName, dombc, nobcname, a_domFine, a_domCoar, false);
    Chombo4::DataIterator dit = a_gridsCoar.dataIterator();
    for(unsigned int ibox = 0; ibox < dit.size(); ++ibox)
    {
      auto      & coarfab = a_phiFToC[dit[ibox]];
      const auto& finefab = a_phiFine[dit[ibox]];
      auto stencil = a_dictionary->getEBStencil(restrictionName, nobcname, a_domFine, a_domCoar, ibox);
      //set resc = Ave(resf) (true is initToZero)
      stencil->apply(coarfab, finefab,  true, 1.0);
    }
  }
} //end namespace hoeb
