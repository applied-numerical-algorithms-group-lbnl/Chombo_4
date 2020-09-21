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
  typedef GraphConstructorFactory< EBBoxData<CELL, int, 1> > devifactint_t;
  typedef GraphConstructorFactory<EBHostData<CELL, int, 1> > hostfactint_t;
  typedef EBDictionary<order, Real, CELL, CELL>                dictionary_t;

  EBPetscSolver(const shared_ptr<GeometryService<order> >  & a_geoserv,
                const shared_ptr<dictionary_t           >  & a_ebdictionary,
                const shared_ptr<LevelData<EBGraph> >      & a_graphs,
                const DisjointBoxLayout                    & a_grids,
                const Box                                  & a_domain,
                string a_stencilName,
                string a_domainBCName[2*DIM],
                string a_ebbcName,
                Real a_dxPoint,   Point a_ghost)
  {
    m_geoserv      =    a_geoserv;
    m_ebdictionary =    a_ebdictionary;
    m_graphs       =    a_graphs;
    m_grids        =    a_grids;
    m_domain       =    a_domain;
    m_stencilName  =    a_stencilName;
    for(int iface = 0; iface < 2*DIM; iface++)
    {
      m_domainBCName[iface] = a_domainBCName[iface];
    }
    m_ebbcName = a_ebbcName;
    
    m_ivghost      =    ProtoCh::getIntVect(a_ghost);
    m_ptghost      =    a_ghost;
    defineGIDS();
  }

  virtual void solve(EBLevelBoxData<CELL, 1>       & a_phi,
                     const EBLevelBoxData<CELL, 1> & a_rhs )
  {
  }
  
private:
  
  //index space--allows us to to matrix to solution space
  LevelData<EBHostData<CELL, int, 1> > m_gids;
  void defineGIDS()
  {
    m_gids.define(m_grids, 1, m_ivghost, hostfactint_t(m_graphs));
    DataIterator dit = m_grids.dataIterator();
    //get the number of points on each proc.
    int numpointsThisProc = 0;
    for(int ibox = 0; ibox < dit.size(); ibox++)
    {
      auto graph = (*m_graphs)[dit[ibox]];
      auto grid  =     m_grids[dit[ibox]];
      Bx  grbx = ProtoCh::getProtoBox(grid);
      for(auto bit = grbx.begin(); bit != grbx.end(); ++bit)
      {
        auto vofs = graph.getVoFs(*bit);
        numpointsThisProc += vofs.size();
      }
    }
    //decide which location maps to the first one of this proc
#ifdef CH_MPI    
    Vector<int> numPtsAllProc;
    gather(numPtsAllProc, numPtsThisProc, 0);
    broadcast(numptsAllProc, 0);
    int startgid = 0;
    for(int iproc = 0; iproc < procID(); iproc++)
    {
      startgid += numptsAllProc[iproc];
    }
#else
    int startgid = 0;
#endif
    int curgid = startgid;
    for(int ibox = 0; ibox < dit.size(); ibox++)
    {
      auto graph = (*m_graphs)[dit[ibox]];
      auto grid  =     m_grids[dit[ibox]];
      Bx  grbx = ProtoCh::getProtoBox(grid);
      for(auto bit = grbx.begin(); bit != grbx.end(); ++bit)
      {
        auto vofs = graph.getVoFs(*bit);
        for(int ivof = 0; ivof < vofs.size(); ivof++)
        {
          m_gids[dit[ibox]](vofs[ivof], 0) = curgid;
          curgid++;
        }
      }
    }
  }


  shared_ptr<GeometryService<order> >  m_geoserv;
  shared_ptr< dictionary_t>            m_ebdictionary;
  shared_ptr<LevelData<EBGraph> >      m_graphs;
  DisjointBoxLayout                    m_grids;
  Box                                  m_domain;
  string m_stencilName;
  string m_domainBCName[2*DIM];
  string m_ebbcName;

  Point   m_ptghost;
  IntVect m_ivghost;
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
