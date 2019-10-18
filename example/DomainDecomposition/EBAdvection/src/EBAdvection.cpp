#include "EBAdvection.H"
#include "NamespaceHeader.H"
/// 
EBAdvection::
EBAdvection(shared_ptr<EBEncyclopedia<2, Real> >   & a_dictionaries,
            shared_ptr<GeometryService<2> >        & a_geoserv,
            shared_ptr<EBLevelBoxData<CELL, DIM> > & a_veloCell,
            const DisjointBoxLayout                & a_grids,
            const Box                              & a_domain,
            const Real                             & a_dx,
            const IntVect                          & a_nghostsrc, 
            const IntVect                          & a_nghostdst)
{
}

///
void 
EBAdvection::
advance(EBLevelBoxData<CELL, 1>       & a_phi,
        const  Real                   & a_dt)
{
}

#include "NamespaceFooter.H"

