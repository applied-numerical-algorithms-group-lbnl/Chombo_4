#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _EBPoissonPetscSolver_H_
#define _EBPoissonPetscSolver_H_

#ifdef CH_USE_PETSC
#include "petsc.h"
#include "petscmat.h"
#include "petscksp.h"
#include "petscviewer.h"
#endif

#include "Chombo_EBLevelBoxData.H"

#ifdef CH_USE_PETSC

/// Framework to solve an elliptic equation using PETsC (PETsC?  petsC?)
/**
   This is an an adaptation to the new reality of Mark Adams' PetscSolver framework.
   If it were any cooler, he would have to add another d to his name.  
   dtg 
   9-18-2020
*/
template <int order>
class EBPetscSolver
{
public:

  EBPetscSolver(const shared_ptr<GeometryService<order> >  & a_geoserv,
                const DisjointBoxLayout                    & a_grids,
                const Box                                  & a_domain,
                string a_stencilName,
                string a_domainBCName[2*DIM],
                string a_ebbcName,
                Real a_dxPoint,   Point a_ghost)
  {
    m_geoserv =    a_geoserv;
    m_grids   =    a_grids;
    m_domain  =    a_domain;
    defineGIDS();
  }

  virtual void solve(EBLeveBoxData<CELL, 1>       & a_phi,
                     const EBLevelBoxData<CELL, 1>& a_rhs )
  {
  }
  
private:
  shared_ptr<GeometryService<order> >  m_geoserv;
  DisjointBoxLayout                    m_grids;
  Box                                  m_domain;
  
  //index space--allows us to to matrix to solution space
  LevelData<EBHostData<CELL, int, NUMCOMPS> > m_gids;
  void defineGIDS()
  {
  }
  
  Mat m_mat;
  void *m_ctx; // pointer for nonlnear solver call backs

  Vec m_xx, m_rr, m_bb;
  SNES m_snes;
  KSP m_ksp;
  PetscInt m_defined;
  PetscErrorCode (*m_function)(SNES,Vec,Vec,void*);
};

#endif   //if petsc is defined

#endif 
