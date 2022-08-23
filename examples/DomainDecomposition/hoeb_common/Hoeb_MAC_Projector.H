#ifndef _Hoeb_MAC_PROJECTOR_H_
#define _Hoeb_MAC_PROJECTOR_H_
#include <cmath>
#include <memory>
#include "Proto.H"
#include "Chombo_EBDictionary.H"
#include "Chombo_EBChombo.H"
#include "Chombo_EBEncyclopedia.H"
#include "Chombo_EBLevelBoxData.H"
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
                       bool a_printStuff = false);

    virtual ~Hoeb_MAC_Projector()
    {
    }

    /// project advectionvel onto its divergence-free subspace.
    /**
       alters input in place.   
     */
    void 
    project(shared_ptr< EBLevelFluxData<1 >  & a_advectionVel) const;


    ///need outside access to this so I can avoid remaking it.
    shared_ptr< EBLevelBoxData<CELL, 1> > m_kappa;

  private:  

  
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
    shared_ptr< EBLevelBoxData<CELL, 1> >  m_kappaConsDiv;
    shared_ptr<graph_distrib_t>            m_graphs;

  
  
    //disallowed because I like strong construction
    Hoeb_MAC_Projector();

  };

}

#endif
