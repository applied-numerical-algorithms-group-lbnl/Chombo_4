#ifndef _Hoeb_Advection_H_
#define _Hoeb_Advection_H_
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
  class Hoeb_Advection
  {
  public:
    typedef CH4_Data_Choreography::DistributedData<EBGraph>                   graph_distrib_t;

    /// 
    Hoeb_Advection(shared_ptr<EBEncyclopedia<2, Real> >   & a_brit,
                   shared_ptr<GeometryService<2> >        & a_geoserv,
                   const DisjointBoxLayout                & a_grids,
                   const Box                              & a_domain,
                   const Real                             & a_dx,
                   const EBIBC                            & a_ebibc,
                   const IntVect                          & a_nghost,
                   bool a_printStuff = false)
    {
      MayDay::Error("not implemented");
    }
    

    virtual ~Hoeb_Advection()
    {
    }


    ///
    /**
       Return kappa*D^C(U phi) + (1-kappa)*D^NC(u phi)  + R(u phi)
       where R(u phi) is the divergence redistribution operator
    **/
    
    shared_ptr<EBLevelBoxData<CELL, 1> >
    hybridDivergence(shared_ptr< EBLevelBoxData<CELL, 1> > & a_scal,
                     shared_ptr< EBLevelFluxData<1    >  > & a_advectionVel)
    {

      zeroRedistributionMass();
      computeHybridDivergence();
      redistributeDivergence();
      return m_hybridDiv;
    }


    ///need outside access to this so I can avoid remaking
    shared_ptr< EBLevelBoxData<CELL, 1> > m_kappa;

  private:  

  
    void  zeroRedistributionMass()  {Chombo4::MayDay::Error("not implemented");}
    void  computeHybridDivergence() {Chombo4::MayDay::Error("not implemented");}
    void  redistributeDivergence()  {Chombo4::MayDay::Error("not implemented");}
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
    shared_ptr< EBLevelBoxData<CELL, 1> >  m_nonConsDiv;
    shared_ptr< EBLevelBoxData<CELL, 1> >  m_hybridDiv;
    shared_ptr< EBLevelBoxData<CELL, 1> >  m_hybridDivergence;
    
    shared_ptr<graph_distrib_t>            m_graphs;

  
  
    ///disallowed because I like strong construction
    Hoeb_Advection();

  };
}

#endif
