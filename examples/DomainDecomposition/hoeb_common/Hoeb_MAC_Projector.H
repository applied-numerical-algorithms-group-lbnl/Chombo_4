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
#include "EBIBC.H"
namespace Chombo4
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
                       bool a_printStuff = false)
    {
      Chombo4::MayDay::Error("not implemented");
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
         gphi holds the pure gradient = grad(phi) = U^p
         Alters inputs in place.   
     */
    void 
    project(shared_ptr< EBLevelFluxData<1 > > & a_advectionVel,
            shared_ptr< EBLevelFluxData<1 > > & a_gphi,
            bool a_printStuff = false) const
    {
      Chombo4::MayDay::Error("not implemented");
    }


    ///
    /**
       used internally and for  diagnostics.
       Take kappa*divergence of mac velocity field
    */
    shared_ptr<EBLevelBoxData<CELL, 1> > & a_divu,
    kappaDivU(shared_ptr<EBLevelFluxData<1> > & a_velo,
              bool a_printStuff = false)
    {
      computeKappaDiv(a_velo, a_printStuff);
      return m_kappaDivU;
    }
    
  private:  
\
  
    ///input at construction
    shared_ptr<EBEncyclopedia<2, Real> >   m_brit;
    shared_ptr<GeometryService<2> >        m_geoserv;
    DisjointBoxLayout                      m_grids;
    Box                                    m_domain;
    Real                                   m_dx;
    EBIBC                                  m_ebibc;
    IntVect                                m_nghost;
    bool                                   m_printStuff;

    ///junk drawer 
    shared_ptr< EBLevelBoxData<CELL, 1> >  m_kappaDivU;
    shared_ptr<graph_distrib_t>            m_graphs;

    void
    computeKappaDiv(shared_ptr<EBLevelFluxData<1> > & a_velo,
                    bool a_printStuff = false)
    {
      MayDay::Error("not implemented");
    }
    //disallowed because I like strong construction
    Hoeb_MAC_Projector();

  };

}

#endif
