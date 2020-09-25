#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include <cmath>
using std::pow;

#include "Chombo_BoxLayoutData.H"
#include "Chombo_FluxBox.H"
//Thus did the sins of petsc through the unholy union of mpi.h with petsc.h
//reflect sevenfold unto the generations of developers
#ifdef CH_MPI
#ifdef CH_USE_PETSC
#include "petsc.h"
#else
#include "mpi.h"
#endif
#endif
#include "Chombo_NamespaceHeader.H"

int LinearizationTest = 0;


template <>
BaseFab<int>* DefaultDataFactory<BaseFab<int> >::create(const Box& box,
                                                        int ncomps,
                                                        const DataIndex& a_datInd) const
{
  return new BaseFab<int>(box, ncomps);
}

template <>
FArrayBox* DefaultDataFactory<FArrayBox>::create(const Box& box, int ncomps,
                                          const DataIndex& a_datInd) const
{
  return new FArrayBox(box, ncomps);
}

FABAliasDataFactory::FABAliasDataFactory(const LayoutData<Real*>& aliases)
{
  define(aliases);
}

void FABAliasDataFactory::define(const LayoutData<Real*>& aliases)
{
  aliasPtrs.define(aliases.boxLayout());
  DataIterator dit = aliases.dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
    {
      aliasPtrs[dit()] = aliases[dit()];
    }
}

FArrayBox* FABAliasDataFactory::create(const Box& box, int ncomps, const DataIndex& a_datInd) const
{
/* Brian thinks this check shouldn't be done.
  const Box& b = aliasPtrs.box(a_datInd);
  if (b != box)
    {
      MayDay::Error("Aliased data holder dimensions do not match across LevelData const.");
    }
*/
  FArrayBox* rtn = new FArrayBox(box, ncomps, aliasPtrs[a_datInd]);
  return rtn;
}

//--FABAliasFlBxDataFactory

FABAliasFlBxDataFactory::FABAliasFlBxDataFactory(
  BoxLayoutData<FluxBox>* a_original,
  const Interval&         a_interval,
  const int               a_dir)
  :
  m_origPointer(a_original),
  m_interval(a_interval),
  m_dir(a_dir)
{
}

void
FABAliasFlBxDataFactory::define(
  BoxLayoutData<FluxBox>* a_original,
  const Interval&         a_interval,
  const int               a_dir)
{
  m_origPointer = a_original;
  m_interval    = a_interval;
  m_dir         = a_dir;
}

FArrayBox*
FABAliasFlBxDataFactory::create(const Box&       a_box,
                                int              a_ncomps,
                                const DataIndex& a_dataInd) const
{
  CH_assert(a_ncomps = m_interval.size());
  FluxBox& origFlBx = m_origPointer->operator[](a_dataInd);
  return new FArrayBox(m_interval, origFlBx[m_dir]);
}

//--FaceFabDataFactory

FaceFabDataFactory::FaceFabDataFactory(const int a_dir)
  :
  m_dir(a_dir)
{
}

void FaceFabDataFactory::define(const int a_dir)
{
  m_dir = a_dir;
}

FArrayBox* FaceFabDataFactory::create(const Box      & a_box   ,
                                      int              a_nComps,
                                      const DataIndex& a_dataInd) const
{
  Box edgeBoxDir(surroundingNodes(a_box, m_dir));
  FArrayBox* faceBox = new FArrayBox(edgeBoxDir, a_nComps);
  
  return faceBox;
}

#include "Chombo_NamespaceFooter.H"
