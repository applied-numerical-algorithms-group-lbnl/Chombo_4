#ifndef _crunch_interface_H_
#define _crunch_interface_H_
#include <cmath>
#include <memory>
#include "Proto.H"
#include "Chombo_EBDictionary.H"
#include "Chombo_EBChombo.H"
#include "Chombo_EBEncyclopedia.H"
#include "EBAdvection.H"
#include "EBMACProjector.H"
#include "Chombo_EBLevelBoxData.H"
#include "Chombo_NamespaceHeader.H"

///class to advect scalars in an eb context (via Trebotich et al.)
class CrunchInterface
{
public:

  /// 
  CrunchInterface(shared_ptr<EBEncyclopedia<2, Real> >   & a_brit,
                  shared_ptr<GeometryService<2> >        & a_geoserv,
                  const DisjointBoxLayout                & a_grids,
                  const Box                              & a_domain);


  ///
  void
  getReactionRates(vector<shared_ptr<EBLevelBoxData<CELL, 1> > > a_species,
                   vector<shared_ptr<EBLevelBoxData<CELL, 1> > > a_reactionRates);

  
private:
};
#include "Chombo_NamespaceFooter.H"

#endif
