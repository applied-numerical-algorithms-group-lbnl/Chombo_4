#include "Hoeb_Advection.H"
#include "EBAdvectionFunctions.H"
#include "EBMACProjector.H"
#include "Chombo_ParmParse.H"
namespace Chombo4
{
/// 
  Hoeb_Conservative_Advection_Op::
  Hoeb_Advection(shared_ptr<EBEncyclopedia<2, Real> >   & a_brit,
                 shared_ptr<GeometryService<2> >        & a_geoserv,
                 const DisjointBoxLayout                & a_grids,
                 const Box                              & a_domain,
                 const Real                             & a_dx,
                 const EBIBC                            & a_ebibc,
                 const IntVect                          & a_nghost,
                 bool a_printStuff = false)
  {
    m_brit       = a_brit      ;
    m_geoserv    = a_geoserv   ;
    m_grids      = a_grids     ;
    m_domain     = a_domain    ;
    m_dx         = a_dx        ;
    m_ebibc      = a_ebibc     ;
    m_nghost     = a_nghost    ;
    m_printStuff = a_printStuff;

    
    m_graphs = a_geoserv->getGraphs(m_domain);
    m_kappaConsDiv  = shared_ptr<EBLevelBoxData<CELL, 1>( new EBLevelBoxData<CELL, 1>(m_grids, m_nghost, m_graphs));

  }


  void
  Hoeb_Conservative_Advection_Op::
  advance(EBLevelBoxData<CELL, 1>          & a_scal,
          shared_ptr< EBLevelFluxData<1> > & a_advectionVel,
          const  Real                      & a_dt) const
  {
    MayDay::Error("not implemented");
  }
} //end namespace chombo4

