#ifndef _Hoeb_MAC_PROJECTOR_H_
#define _Hoeb_MAC_PROJECTOR_H_
#include <cmath>
#include <memory>
#include "Proto.H"
#include "Chombo_EBDictionary.H"
#include "Chombo_EBChombo.H"
#include "Chombo_EBEncyclopedia.H"
#include "Chombo_EBLevelBoxData.H"
#include "Chombo_EBLevelFluxData.H"
#include "Hoeb_PetscSolver.H"
#include "Hoeb_Neighborhood.H"
#include "Hoeb_Utilities.H"
#include "EBIBC.H"

namespace hoeb
{

  ///class to project face-centered velocity fields.
  /**
     for face-valued fields, subtracts off the pure gradient part of a vector field
     Homogeneous boundary conditions for phi and vel are required for now.
     as there is real work to do to make certain those boundary conditions
     are implmented with appropriate accuracy.
  **/

  class Hoeb_MAC_Projector
  {
  public:
    typedef CH4_Data_Choreography::DistributedData<EBGraph>                  graph_distrib_t;
    
    Hoeb_MAC_Projector(shared_ptr<EBEncyclopedia<2, Real> >                      & a_brit,
                       shared_ptr<GeometryService<2> >                           & a_geoserv,
                       const DisjointBoxLayout                                   & a_grids,
                       const Box                                                 & a_domain,
                       const Real                                                & a_dx,
                       const EBIBC                                               & a_ebibc,
                       const IntVect                                             & a_nghost,
                       bool                                                        a_printStuff = false)
    {
      m_brit       = a_brit     ;
      m_geoserv    = a_geoserv  ;
      m_grids      = a_grids    ;
      m_domain     = a_domain   ;
      m_dx         = a_dx       ;
      m_ebibc      = a_ebibc    ;
      m_nghost     = a_nghost   ;      
      m_printStuff = a_printStuff;
      m_graphs     = a_geoserv.getGraphs(a_domain);
      a_ebibc.projectionStencilStrings(m_domBCNamePoisson);
      m_EBBCNamePoisson = StencilNames::Neumann;
      m_EBBCNameGrad = string("no_bcs");
      m_stenNameGrad = string("hoeb_mac_gradient");
      defineInternals();
    }

    virtual ~Hoeb_MAC_Projector()
    {
    }

    /// project advectionvel onto its divergence-free subspace.
    /**
       Input: 
       Vector field on faces U* = U^V + U^P
       div(U^V) = 0
       Nabla x U^p = 0;
       U^p = grad(phi)
       Output:
       advectionVel holds rotational component of vector field = U^V
       m_gphi holds the pure gradient = grad(phi) = U^p
       Alters input in place.   
    */
    void 
    project(shared_ptr< EBLevelFluxData<1> > & a_advectionVel,
            bool a_printStuff = false) const
    {
      computeKappaDivU(a_advectionVel, a_printStuff);

      m_poissonSolver->solve(*m_phi, *m_kappaDivU);

      computeGradient(*m_phi);

      (*a_advectionVel) -= (*m_gradient);
    }


    ///
    /**
       used internally and for  diagnostics.
       compute and return kappa*divergence of mac velocity field
    */
    shared_ptr<EBLevelBoxData<CELL, 1> >
    kappaDivergence(shared_ptr<EBLevelFluxData<1> > & a_velo,
                    bool a_printStuff = false)
    {
      computeKappaDiv(a_velo, a_printStuff);
      return m_kappaDivU;
    }
    
    ///
    /**
       used internally and for  diagnostics.
       compute and return mac gradient of cell-averaged phi
    */
    shared_ptr<EBLevelFluxData<1> > 
    gradient(shared_ptr<EBLevelBoxData<CELL, 1> > & a_phi,
             bool a_printStuff = false)
    {
      computeGradient(a_phi, a_printStuff);
      return m_gradient;
    }
  
    ///input at construction
    shared_ptr<EBEncyclopedia<2, Real> >   m_brit;
    shared_ptr<GeometryService<2> >        m_geoserv;
    DisjointBoxLayout                      m_grids;
    Box                                    m_domain;
    Real                                   m_dx;
    EBIBC                                  m_ebibc;
    IntVect                                m_nghost;
    bool                                   m_printStuff;

    ///junk drawer--some of which can be reused (hence the public)
    shared_ptr< EBLevelBoxData<CELL, 1> >  m_kappaDivU;
    shared_ptr< EBLevelBoxData<CELL, 1> >  m_phi;
    shared_ptr< EBLevelFluxData<     1> >  m_gradient;
    shared_ptr<graph_distrib_t>            m_graphs;
    shared_ptr<Hoeb_PetscSolver<4> >       m_solver;
    string m_domBCNamePoisson[2*DIM];
    string m_EBBCNamePoisson;
    string m_stenNameGrad;
    string m_EBBCNameGrad;
    private:
    
    ///this stuff is too context-dependent to allow into the wild.
    void defineInternals()
    {

      m_kappaDivU       = shared_ptr< EBLevelBoxData<CELL, 1> >  (new EBLevelBoxData<CELL, 1> (m_grids, m_nghost, m_graphs));
      m_phi             = shared_ptr< EBLevelBoxData<CELL, 1> >  (new EBLevelBoxData<CELL, 1> (m_grids, m_nghost, m_graphs));
      m_gradient        = shared_ptr< EBLevelFluxData<     1> >  (new EBLevelFluxData<     1> (m_grids, m_nghost, m_graphs));

      auto cellDict= m_brit->m_cellToCell;
      Real alpha = 0; Real beta = 1; //poisson
      m_solver = shared_ptr< Hoeb_PetscSolver<4> >
        (new Hoeb_PetscSolver<4>(m_geoserv, cellDict,    m_graphs, m_grids, m_domain,
                                 m_domBCNamePoisson, m_EBBCNamePoisson,
                                 m_dx, alpha, beta,  m_nghost, a_printStuff));

      registerStencils();
    }

    void registerStencils()
    {
      //increment divergence by face difference
      m_brit->registerFaceToCell(StencilNames::DivergeFtoC         ,
                                 StencilNames::NoBC,
                                 StencilNames::NoBC,
                                 m_domain, m_domain, false); // false is for needDiagonal

      //face-centered gradient
      DataIterator dit = m_grids.dataIterator();
      for(int ibox = 0; ibox  < dit.size(); ibox++)
      {
        bool needDiag = false;  
        vector<EBIndex<     XFACE>                destXFace;
        vector<LocalStencil<CELL,  data_t>  >     stenXFace;
        vector<EBIndex<     YFACE>                destYFace;
        vector<LocalStencil<CELL,  data_t>  >     stenYFace;
        vector<EBIndex<     ZFACE>                destZFace;
        vector<LocalStencil<CELL,  data_t>  >     stenZFace;
        Proto::Box                                srcValid;
        Proto::Box                                dstValid;
        Proto::Box                                srcDomain;
        Proto::Box                                dstDomain;
        Point                                     srcGhost;
        Point                                     dstGhost;
        
        getGradientStencilComponents(destXFace,
                                     stenXFace,
                                     destYFace,
                                     stenYFace,
                                     destZFace,
                                     stenZFace,
                                     srcValid,
                                     dstValid,
                                     srcDomain,
                                     dstDomain,
                                     srcGhost,
                                     dstGhost,
                                     dit[ibox],
                                     ibox);
        
        m_brit->registerCellToFace(m_stenNameGradient, m_ebbcNameGradient,
                                   destXFace,stenXFace,
                                   destYFace, stenYFace,
                                   destZFace, stenZFace,
                                   srcValid,  dstValid,
                                   srcDomain, dstDomain,
                                   srcGhost,  dstGhost,
                                   needDiag,  ibox);                
      }
      
    }

    //build stencils for higher order gradient
    void
    getGradientStencilComponents(vector<EBIndex<     XFACE>               & a_destXFace,
                                 vector<LocalStencil<CELL,  data_t> >     & a_stenXFace,
                                 vector<EBIndex<     YFACE>               & a_destYFace,
                                 vector<LocalStencil<CELL,  data_t> >     & a_stenYFace,
                                 vector<EBIndex<     ZFACE>               & a_destZFace,
                                 vector<LocalStencil<CELL,  data_t> >     & a_stenZFace,
                                 Proto::Box                               & a_srcValid,
                                 Proto::Box                               & a_dstValid,
                                 Proto::Box                               & a_srcDomain,
                                 Proto::Box                               & a_dstDomain,
                                 Point                                    & a_srcGhost,
                                 Point                                    & a_dstGhost,
                                 const DataIndex                          & a_dit,
                                 unsigned int                               a_ibox)
    {
      const auto& graph = (*m_graphs)[a_dit];
      a_dstValid  = ProtoCh::getProtoBox(m_grids[a_dit]);
      a_srcValid  = a_dstValid;
      a_srcDomain = ProtoCh::getProtoBox(m_domain);
      a_dstDomain = ProtoCh::getProtoBox(m_domain);
      a_srcGhost =  ProtoCh::getPoint(m_nghost);
      a_dstGhost =  ProtoCh::getPoint(m_nghost);
      a_destXFace = graph.getAllFaces<XFACE>(a_dstValid);
      a_destYFace = graph.getAllFaces<YFACE>(a_dstValid);
#if DIM==3        
      a_destZFace = graph.getAllFaces<ZFACE>(a_dstValid);
#else
      a_dstZFace.resize(0);
#endif

      getGradientStencil<XFACE>(a_stenXFace,  a_destXFace, a_dit, 0);
      getGradientStencil<YFACE>(a_stenYFace,  a_destYFace, a_dit, 1);
#if DIM==3                                   
      getGradientStencil<ZFACE>(a_stenZFace,  a_destZFace, a_dit, 2);
#else
      a_stenZFace.resize(0);
#endif      
    }

    template<CENTERING cent>
    getGradientStencil(vector<LocalStencil<CELL, data_t> > & a_stencils,
                       const vector<EBIndex<cent> >        & a_faces,
                       const DataIndex                     & a_dit,
                       int a_facedir)
    {

      const auto& graph = (*m_graphs)[a_dit];
      a_stencils.resize(a_faces.size());
      string dobcName[2*DIM];
      string ebbcName = StencilNames::Neumann;
      m_ebbc.projectionStencilStrings(dobcName);
      for(int iface = 0; iface < faces.size(); iface++)
      {
        const auto& face = a_faces[iface];
        //vofside is used to short circuit all this stuff for neumann boundary conditions
        Side::LoHiSide vofside = Side::Lo;
        if(face.m_isBoundary && !m_domain.contains(face.m_pt))
        {
          vofsize = Side::Hi;
        }
        //this is the whole stencil, including inhomogeneous bits.   
        // rest of it into a data holder.
        CompositeStencil intFluxDA_sten = getDevendranIntFluxDACompositeStencil(face, dobcName, ebbcName, m_geoserv, m_domain, a_dit, m_dx, a_facedir, vofside, a_printStuff);
        a_stencils[iface] = intFluxDA_sten.m_cellSten;
      }
    }
    
    void
    computeKappaDiv(shared_ptr<EBLevelFluxData<1> > & a_velo,
                    bool a_printStuff = false)
    {

      CH_TIME("Hoeb_EBMACProjector::kappaDivU");
      if(a_printStuff)
      {
        pout() << "ebmacproj::kappaDivU: exchanging velocity input" << endl;
      }
      a_velo.exchange(a_printStuff);

      if(a_printStuff)
      {
        pout() << "applying velocity boundary conditions and applying divergence stencils" << endl;
      }
      DataIterator dit = m_grids.dataIterator();

      for(int ibox = 0; ibox < dit.size(); ++ibox)
      {

        Bx   grid   =  ProtoCh::getProtoBox(m_grids[dit[ibox]]);
        const EBGraph  & graph = (*m_graphs)[dit[ibox]];
        Bx  grown   =  grid.grow(ProtoCh::getPoint(m_nghost));

        const EBFluxData<Real, 1>&  centroidFlux(grown, graph, useStack);
        EBFluxStencil<2, Real> stencils = 
          m_brit->getFluxStencil(StencilNames::InterpToFaceCentroid, StencilNames::NoBC, m_domain, m_domain, ibox);
        EBFluxData<Real,1>& flux = (*a_velo)[dit[ibox]]; 

        //I think this is harder in hoeb land.  for now, let us stay in homogeneous land for this operator
        //applyVeloBoundaryConditions(flux, dit[ibox]);

        auto& kapdiv =  (*m_kappaDivU)[dit[ibox]];
        kapdiv.setVal(0.);
        for(unsigned int idir = 0; idir < DIM; idir++)
        {
          bool initToZero = false;
          m_brit->applyFaceToCell(StencilNames::DivergeFtoC, StencilNames::NoBC, m_domain, kapdiv, flux,
                                  idir, ibox, initToZero, 1.0);
        }
      }
      
    }
    
    void
    computeGradient(shared_ptr<EBLevelBoxData<CELL, 1> > & a_phi,
                    bool a_printStuff = false)
    {


      //face-centered gradient
      DataIterator dit = m_grids.dataIterator();
      for(int ibox = 0; ibox  < dit.size(); ibox++)
      {
        bool needDiag = false;  
        vector<EBIndex<     XFACE>                destXFace;
        vector<LocalStencil<CELL,  data_t>  >     stenXFace;
        vector<EBIndex<     YFACE>                destYFace;
        vector<LocalStencil<CELL,  data_t>  >     stenYFace;
        vector<EBIndex<     ZFACE>                destZFace;
        vector<LocalStencil<CELL,  data_t>  >     stenZFace;
        Proto::Box                                srcValid;
        Proto::Box                                dstValid;
        Proto::Box                                srcDomain;
        Proto::Box                                dstDomain;
        Point                                     srcGhost;
        Point                                     dstGhost;
        
        getGradientStencilComponents(destXFace,
                                     stenXFace,
                                     destYFace,
                                     stenYFace,
                                     destZFace,
                                     stenZFace,
                                     srcValid,
                                     dstValid,
                                     srcDomain,
                                     dstDomain,
                                     srcGhost,
                                     dstGhost,
                                     dit[ibox],
                                     ibox);
        
        m_brit->registerCellToFace(m_stenNameGradient, m_ebbcNameGradient,
                                   destXFace,stenXFace,
                                   destYFace, stenYFace,
                                   destZFace, stenZFace,
                                   srcValid,  dstValid,
                                   srcDomain, dstDomain,
                                   srcGhost,  dstGhost,
                                   needDiag,  ibox);                
      }
    }
    

    
    //disallowed because I like strong construction
    Hoeb_MAC_Projector();

  };

}

#endif
