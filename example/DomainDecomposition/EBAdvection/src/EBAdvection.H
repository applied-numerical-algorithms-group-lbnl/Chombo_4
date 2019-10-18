#ifndef _EBAdvection_H_
#define _EBAdvection_H_
#include <cmath>
#include <memory>
#include "Proto.H"
#include "EBDictionary.H"
#include "EBChombo.H"
#include "EBEncyclopedia.H"
#include "EBLevelBoxData.H"
#include "NamespaceHeader.H"

///class to advect scalars in a eb context
class EBAdvection
{
public:

  /// 
  EBAdvection(shared_ptr<EBEncyclopedia<2, Real> >   & a_dictionaries,
              shared_ptr<GeometryService<2> >        & a_geoserv,
              shared_ptr<EBLevelBoxData<CELL, DIM> > & a_veloCell,
              const DisjointBoxLayout                & a_grids,
              const Box                              & a_domain,
              const Real                             & a_dx,
              const IntVect                          & a_nghostsrc, 
              const IntVect                          & a_nghostdst);

  /// advance (via Trebotich et al.) in  an eb context
  void 
  advance(EBLevelBoxData<CELL, 1>   & a_scal,
          const  Real               & a_dt);

private:

  /// weak construction introduces unnecessary complications
  EBAdvection()
  {
  }


};
#include "NamespaceFooter.H"

#endif
