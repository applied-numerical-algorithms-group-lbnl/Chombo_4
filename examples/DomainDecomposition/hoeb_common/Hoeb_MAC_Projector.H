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
#include "EBIBC.H"

namespace hoeb
{

///class to compute kappa*div_c(u s);
  class Hoeb_MAC_Projector
  {
  public:
    typedef CH4_Data_Choreography::DistributedData<EBGraph> graph_distrib_t;
    
    /// 
    Hoeb_MAC_Projector(shared_ptr<EBEncyclopedia<2, Real> >   & a_brit,
                       shared_ptr<GeometryService<2> >        & a_geoserv,
                       const DisjointBoxLayout                & a_grids,
                       const Box                              & a_domain,
                       const Real                             & a_dx,
                       const EBIBC                            & a_ebibc,
                       const IntVect                          & a_nghost,
                       bool                                     a_printStuff = false)
    :BaseMACProjector(shared_ptr<EBEncyclopedia<2, Real> >   &  a_brit,
                       shared_ptr<GeometryService<2> >        & a_geoserv,
                       const DisjointBoxLayout                & a_grids,
                       const Box                              & a_domain,
                       const Real                             & a_dx,
                       const EBIBC                            & a_ebibc,
                       const IntVect                          & a_nghost,
                       bool                                     a_printStuff)
    {
      m_brit      = a_brit     ;
      m_geoserv   = a_geoserv  ;
      m_grids     = a_grids    ;
      m_domain    = a_domain   ;
      m_dx        = a_dx       ;
      m_ebibc     = a_ebibc    ;
      m_nghost    = a_nghost   ;      
      m_printStuff = a_printStuff;
      m_graphs = a_geoserv.getGraphs(a_domain);
      a_ebibc.projectionStencilStrings(m_domBCName);
      m_EBBCName = StencilNames::Neumann;
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

      m_poissonSolver->modifyRHSForInhomogeneity(*m_kappaDivU, a_printStuff);
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
    string m_domBCName[2*DIM];
    string m_EBBCName;
    private:
    
    ///this stuff is too context-dependent to allow into the wild.
    void defineInternals()
    {

      auto cellDict= m_brit->m_cellToCell;
      m_kappaDivU       = shared_ptr< EBLevelBoxData<CELL, 1> >  (new EBLevelBoxData<CELL, 1> (m_grids, m_nghost, m_graphs));
      m_phi             = shared_ptr< EBLevelBoxData<CELL, 1> >  (new EBLevelBoxData<CELL, 1> (m_grids, m_nghost, m_graphs));
      m_gradient        = shared_ptr< EBLevelFluxData<     1> >  (new EBLevelFluxData<     1> (m_grids, m_nghost, m_graphs));
      Real alpha = 0; Real beta = 1; //poisson
      m_solver = shared_ptr< Hoeb_PetscSolver<4> >
        (new Hoeb_PetscSolver<4>(m_geoserv, cellDict,    m_graphs, m_grids, m_domain,
                                 m_domBCName, m_EBBCName, m_dx, alpha, beta,  m_nghost,
                                 a_printStuff));
      
      //increment divergence by face difference
      m_brit->registerFaceToCell( StencilNames::DivergeFtoC         , StencilNames::NoBC,    StencilNames::NoBC, m_domain, m_domain, needDiag);
    }

    void
    computeKappaDiv(shared_ptr<EBLevelFluxData<1> > & a_velo,
                    bool a_printStuff = false)
    {
      
    }


    void
    computeGradient(shared_ptr<EBLevelBoxData<CELL, 1> > & a_phi,
                    bool a_printStuff = false)
    {
      MayDay::Error("not implemented");
    }

    void
    defineSolver()
    {
      MayDay::Error("not implemented");
    }
    
    
    //disallowed because I like strong construction
    Hoeb_MAC_Projector();

  };

}

#endif
