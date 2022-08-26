#ifndef _Hoeb_Advection_H_
#define _Hoeb_Advection_H_
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
  class Hoeb_Advection
  {
  public:

    /// 
    Hoeb_Advection(shared_ptr<EBEncyclopedia<2, Real> >   & a_brit,
                   shared_ptr<GeometryService<2> >        & a_geoserv,
                   const DisjointBoxLayout                & a_grids,
                   const Box                              & a_domain,
                   const Real                             & a_dx,
                   const EBIBC                            & a_ebibc,
                   const IntVect                          & a_nghost,
                   bool a_printStuff = false);

    virtual ~Hoeb_Advection()
    {
    }


    /// use 3rd order in time tvd to advance solution
    void 
    advance(EBLevelBoxData<CELL, 1>          & a_scal,
            shared_ptr< EBLevelFluxData<1> > & a_advectionVel,
            const  Real                      & a_dt) const;


    ///need outside access to this so I can avoid remaking
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

  
  
    ///disallowed because I like strong construction
    Hoeb_Advection();

  };
}

#endif
