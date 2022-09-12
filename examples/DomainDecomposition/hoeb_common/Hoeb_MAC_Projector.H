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
     Homogeneous boundary conditions for phi  are required for now
     as there is real work to do to make certain those boundary conditions
     are implmented with appropriate accuracy.
  **/
  template <int order>
  class Hoeb_MAC_Projector
  {
  public:
    typedef CH4_Data_Choreography::DistributedData<EBGraph>                  graph_distrib_t;
    
    Hoeb_MAC_Projector(shared_ptr<EBEncyclopedia<order, Real> >                  & a_brit,
                       shared_ptr<GeometryService<order> >                       & a_geoserv,
                       const Chombo4::DisjointBoxLayout                          & a_grids,
                       const Chombo4::Box                                        & a_domain,
                       const Real                                                & a_dx,
                       const EBIBC                                               & a_ebibc,
                       const IntVect                                             & a_nghost,
                       bool                                                        a_printStuff = false)
    {
      CH_TIME("Hoeb_MACProjector::constructor");
      m_brit       = a_brit     ;
      m_geoserv    = a_geoserv  ;
      m_grids      = a_grids    ;
      m_domain     = a_domain   ;
      m_dx         = a_dx       ;
      m_ebibc      = a_ebibc    ;
      m_nghost     = a_nghost   ;      
      m_printStuff = a_printStuff;
      m_graphs     = a_geoserv->getGraphs(a_domain);
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
       VELOCITY BOUNDARY CONDITIONS ARE NOT ENFORCED at the domain boundaries.
       Have the velocity at domain boundary faces set correctly coming into the 
       projection.
    */
    void 
    project(shared_ptr< EBLevelFluxData<1> > & a_advectionVel,
            bool a_printStuff = false) 
    {
      CH_TIME("Hoeb_MACProjector::project");
      computeKappaDiv(a_advectionVel, a_printStuff);
      if(a_printStuff)
      {
        Real error = m_kappaDivU->maxNorm(0);
        Chombo4::pout() << "Hoeb_MACProjector:: project: divergence befo projection = " << error;
      }

      m_solver->solve(*m_phi, *m_kappaDivU);

      computeGradient(m_phi, a_printStuff);

      (*a_advectionVel) -= (*m_gradient);
      if(a_printStuff)
      {
        computeKappaDiv(a_advectionVel, a_printStuff);
        Real error = m_kappaDivU->maxNorm(0);
        Chombo4::pout() << "Hoeb_MACProjector:: project: divergence after projection = " << error;
      }
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
    shared_ptr<EBEncyclopedia<order, Real> >   m_brit;
    shared_ptr<GeometryService<order> >        m_geoserv;
    Chombo4::DisjointBoxLayout             m_grids;
    Chombo4::Box                           m_domain;
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

      CH_TIME("Hoeb_MACProjector::defineInternals");
      m_kappaDivU       = shared_ptr< EBLevelBoxData<CELL, 1> >  (new EBLevelBoxData<CELL, 1> (m_grids, m_nghost, m_graphs));
      m_phi             = shared_ptr< EBLevelBoxData<CELL, 1> >  (new EBLevelBoxData<CELL, 1> (m_grids, m_nghost, m_graphs));
      m_gradient        = shared_ptr< EBLevelFluxData<     1> >  (new EBLevelFluxData<     1> (m_grids, m_nghost, m_graphs));

      auto cellDict= m_brit->m_cellToCell;
      Real alpha = 0; Real beta = 1; //poisson
      bool printStuff  = false;
      Point ptghost = ProtoCh::getPoint(m_nghost);
      m_solver = shared_ptr< Hoeb_PetscSolver<4> >
        (new Hoeb_PetscSolver<4>(m_geoserv, cellDict,    m_graphs, m_grids, m_domain,
                                 m_domBCNamePoisson, m_EBBCNamePoisson,
                                 m_dx, alpha, beta,  ptghost, printStuff));

      registerStencils();
    }

    void registerStencils()
    {
      CH_TIME("Hoeb_MACProjector::registerStencil");
      //increment divergence by face difference
      m_brit->registerFaceToCell(StencilNames::DivergeFtoC         ,
                                 StencilNames::NoBC,
                                 StencilNames::NoBC,
                                 m_domain, m_domain, false); // false is for needDiagonal

      //face-centered gradient
      Chombo4::DataIterator dit = m_grids.dataIterator();
      for(int ibox = 0; ibox  < dit.size(); ibox++)
      {
        bool needDiag = false;  
        vector<EBIndex<     XFACE>          >     destXFace;
        vector<LocalStencil<CELL,  Real  >  >     stenXFace;
        vector<EBIndex<     YFACE>          >     destYFace;
        vector<LocalStencil<CELL,  Real  >  >     stenYFace;
        vector<EBIndex<     ZFACE>          >     destZFace;
        vector<LocalStencil<CELL,  Real  >  >     stenZFace;
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
        
        m_brit->registerCellToFace(m_stenNameGrad,
                                   m_EBBCNameGrad,
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
    inline void
    getGradientStencilComponents(vector<EBIndex<     XFACE>         >     & a_destXFace,
                                 vector<LocalStencil<CELL,  Real  > >     & a_stenXFace,
                                 vector<EBIndex<     YFACE>         >     & a_destYFace,
                                 vector<LocalStencil<CELL,  Real  > >     & a_stenYFace,
                                 vector<EBIndex<     ZFACE>         >     & a_destZFace,
                                 vector<LocalStencil<CELL,  Real  > >     & a_stenZFace,
                                 Proto::Box                               & a_srcValid,
                                 Proto::Box                               & a_dstValid,
                                 Proto::Box                               & a_srcDomain,
                                 Proto::Box                               & a_dstDomain,
                                 Point                                    & a_srcGhost,
                                 Point                                    & a_dstGhost,
                                 const Chombo4::DataIndex                 & a_dit,
                                 unsigned int                               a_ibox)
    {
      const auto& graph = (*m_graphs)[a_dit];
      a_dstValid  = ProtoCh::getProtoBox(m_grids[a_dit]);
      a_srcValid  = a_dstValid;
      a_srcDomain = ProtoCh::getProtoBox(m_domain);
      a_dstDomain = ProtoCh::getProtoBox(m_domain);
      a_srcGhost =  ProtoCh::getPoint(m_nghost);
      a_dstGhost =  ProtoCh::getPoint(m_nghost);
      a_destXFace = getAllFaces<XFACE>(a_dstValid, graph);
      a_destYFace = getAllFaces<YFACE>(a_dstValid, graph);
#if DIM==3                                             
      a_destZFace = getAllFaces<ZFACE>(a_dstValid, graph);
#else
      a_destZFace.resize(0);
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
    inline void
    getGradientStencil(vector<LocalStencil<CELL, Real  > > & a_stencils,
                       const vector<EBIndex<cent> >        & a_faces,
                       const Chombo4::DataIndex            & a_dit,
                       int a_facedir)
    {

      const auto& graph = (*m_graphs)[a_dit];
      auto grid   = m_grids[a_dit];
      auto gridbx = ProtoCh::getProtoBox(m_grids[a_dit]);
      gridbx = gridbx.grow(1);
      gridbx &= ProtoCh::getProtoBox(m_domain);

      a_stencils.resize(a_faces.size());
      string dobcName[2*DIM];
      string ebbcName = StencilNames::Neumann;
      m_ebibc.projectionStencilStrings(dobcName);
      //need area fractions to go between integrals and average values
      EBHostData<cent, Real  , 1> faceAlpha(gridbx, graph);
      auto geolev  = (*m_geoserv)[m_domain];
      Chombo4::fillFaceAlpha<cent, order>(faceAlpha, grid, a_dit, *geolev);
      
      for(int iface = 0; iface < a_faces.size(); iface++)
      {
        const auto& face = a_faces[iface];
        //vofside is used to short circuit all this stuff for neumann boundary conditions
        Side::LoHiSide vofside = Side::Lo;
        if(face.m_isBoundary && !m_domain.contains(face.m_pt))
        {
          vofside = Side::Hi;
        }
        //This is the whole stencil, including inhomogeneous bits.   
        // For inhomogeneous (phi) boundary conditions, the inhomogeneous contribution
        //would have to be archived here and dragged out again for operator evaluation.
        //Since the application on which I am working does not require inhomogeneous pressure
        //boundary conditions, this bit is left as exercise to the reader.
        bool printStuff = false;
        CompositeStencil intFluxDA_composite_sten =
          getDevendranIntFluxDACompositeStencil(face, dobcName, ebbcName, m_geoserv, m_domain, a_dit, m_dx, a_facedir, vofside, printStuff);

        //since the above sten gives int(grad phi) dA, we need to divide out the area to get <grad phi>
        //if there is no real area, set the stencil to zero.
        Real   faceArea = faceAlpha(face, 0);
        auto gradSten = intFluxDA_composite_sten.m_cellSten;
        Real   invAreaFactor = 0;
        if(faceArea > 1.0e-10)
        {
          invAreaFactor = 1.0/faceArea;
        }
        gradSten *= invAreaFactor;
        a_stencils[iface] = gradSten;
      }
    }
    
    void
    computeKappaDiv(shared_ptr<EBLevelFluxData<1> > & a_velo,
                    bool a_printStuff = false)
    {

      CH_TIME("Hoeb_EBMACProjector::kappaDivU");
      if(a_printStuff)
      {
        Chombo4::pout() << "ebmacproj::kappaDivU: exchanging velocity input" << endl;
      }
      a_velo->exchange(a_printStuff);

      if(a_printStuff)
      {
        Chombo4::pout() << "applying velocity boundary conditions and applying divergence stencils" << endl;
      }
      Chombo4::DataIterator dit = m_grids.dataIterator();

      for(int ibox = 0; ibox < dit.size(); ++ibox)
      {

        EBFluxData<Real,1>& flux = (*a_velo)[dit[ibox]]; 
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
      Chombo4::DataIterator dit = m_grids.dataIterator();
      for(int ibox = 0; ibox  < dit.size(); ibox++)
      {
        bool needDiag = false;  
        vector<EBIndex<     XFACE>          >     destXFace;
        vector<LocalStencil<CELL,  Real  >  >     stenXFace;
        vector<EBIndex<     YFACE>          >     destYFace;
        vector<LocalStencil<CELL,  Real  >  >     stenYFace;
        vector<EBIndex<     ZFACE>          >     destZFace;
        vector<LocalStencil<CELL,  Real  >  >     stenZFace;
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
        
        m_brit->registerCellToFace(m_stenNameGrad,
                                   m_EBBCNameGrad,
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
