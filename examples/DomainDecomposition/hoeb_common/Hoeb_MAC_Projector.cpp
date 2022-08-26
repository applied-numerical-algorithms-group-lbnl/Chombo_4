#include <cmath>
#include <memory>
#include "Proto.H"
#include "Chombo_EBDictionary.H"
#include "Chombo_EBChombo.H"
#include "Chombo_EBEncyclopedia.H"
#include "Chombo_EBLevelBoxData.H"
#include "EBIBC.H"
#include "Hoeb_MAC_Projector.H"
namespace Chombo4
{

  Hoeb_MAC_Projector::
  Hoeb_MAC_Projector(shared_ptr<EBEncyclopedia<2, Real> >   & a_brit,
                     shared_ptr<GeometryService<2> >        & a_geoserv,
                     const DisjointBoxLayout                & a_grids,
                     const Box                              & a_domain,
                     const Real                             & a_dx,
                     const EBIBC                            & a_ebibc,
                     const IntVect                          & a_nghost,
                     bool a_printStuff = false)
  {
  }


  void 
  Hoeb_MAC_Projector::
  project(shared_ptr< EBLevelFluxData<1 >  & a_advectionVel) const
  {
  }

}

