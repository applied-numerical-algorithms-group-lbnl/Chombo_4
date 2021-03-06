#include "CrunchInterface.H"
#include "Chombo_NamespaceHeader.H"
///class to advect scalars in an eb context (via Trebotich et al.)
CrunchInterface::
CrunchInterface(shared_ptr<EBEncyclopedia<2, Real> >   & a_brit,
                shared_ptr<GeometryService<2> >        & a_geoserv,
                const DisjointBoxLayout                & a_grids,
                const Box                              & a_domain)
{
}


void
CrunchInterface::
getReactionRates(vector<shared_ptr<EBLevelBoxData<CELL, 1> > > a_species,
                 vector<shared_ptr<EBLevelBoxData<CELL, 1> > > a_reactionRates)
{
}

  
#include "Chombo_NamespaceFooter.H"

