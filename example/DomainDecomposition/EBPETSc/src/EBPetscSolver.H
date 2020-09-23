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
  typedef GraphConstructorFactory< EBBoxData<CELL, int, 1> >  devifactint_t;
  typedef GraphConstructorFactory<EBHostData<CELL, int, 1> >  hostfactint_t;
  typedef GraphConstructorFactory<EBHostData<CELL,Real, 1> > hostfactreal_t;
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
    //create map of locations in space to matrix row.
    defineGIDS();
    //create space for a matrix and the necessary vectors.
    createMatrixAndVectors();
    //put actual values into the matrix
    formMatrix();
    setupSolver();
  }

  virtual void solve(EBLevelBoxData<CELL, 1>       & a_phi,
                     const EBLevelBoxData<CELL, 1> & a_rhs )
  {

    
    auto dbl = m_grids;
    int nc = 1;
    PetscErrorCode ierr;
    const PetscInt nc = a_rhs.nComp();
#ifdef CH_MPI
    MPI_Comm wcomm = Chombo_MPI::comm;
#else
    MPI_Comm wcomm = PETSC_COMM_SELF;
#endif
    
    // form RHS -- m_bb
    Real idx2 = 1.e0/(this->m_dx*this->m_dx);
    Real addV;


    // this is an interface in case some operations are needed on the rhs
    // the default does nothing
    ierr = VecSetOption(m_bb,VEC_IGNORE_OFF_PROC_ENTRIES,PETSC_TRUE);CHKERRQ(ierr);
    ierr = VecSetOption(m_xx,VEC_IGNORE_OFF_PROC_ENTRIES,PETSC_TRUE);CHKERRQ(ierr);
    // add X and B from Chombo to PETSc and add stuff for EB to B
    ierr = VecSet( m_xx, 0.);CHKERRQ(ierr);
    ierr = VecSet( m_bb, 0.);CHKERRQ(ierr);

    chomboToPetsc(m_bb, a_rhs);
    chomboToPetsc(m_xx, a_phi);
    // solve
#ifdef CH_MPI
    //with all the host to device and device to host stuff going on, this is probably
    //unnessary.   The peTsC team, however, tells me that the communication between
    //host and device will unnecessary soon so I left the barrier in the code.
    MPI_Barrier(Chombo_MPI::comm);
#endif
    ierr = KSPSolve( m_ksp, m_bb, m_xx );CHKERRQ(ierr);

    // put solution into output
    ierr = petscToChombo( a_phi, m_xx );CHKERRQ(ierr);
    a_phi.exchange();

    return 0;
  }

  
private:
  int chomboToPetsc(Vec                           & a_dst,
                    const EBLevelBoxData<CELL, 1> & a_deviceSrc)
  {
    PetscErrorCode ierr;

    //copy data to host.  The pEtsC team has promised to make this
    //step unnecessary by providing the pointer to the vector that lives
    //on the device
    hostfactreal_t factory(m_graphs);
    LevelData<EBHostData<CELL, , 1> > hostSrc(m_grids, 1, a_deviceSrc.ghost(), factory);
    EBLevelBoxData<CELL, 1>::copyToHost(hostSrc, a_deviceSrc);
    
    DataIterator dit = m_grids.dataIterator();
    for (int ibox=0;ibox< dit.size(); ibox++)
    {
      auto& hostfab =   hostSrc[dit[ibox]];
      auto& graph = (*m_graphs)[dit[ibox]];
      auto& grid  =     m_grids[dit[ibox]];
      Bx grbx = ProtoCh::getProtoBox(grid);
      for(auto bit = grbx.begin(); bit != grbx.end(); ++bit)
      {
        auto pt = *bit;
        auto vofs = m_graph.getVoFs(pt);
        for(int ivof = 0; ivof < vofs.size(); ivof++)
        {
          auto vof = vofs[ivof];
          PetscInt ki = m_gids[dit[ibox]](vof, 0);
          Real v = hostfab(vof, 0);
          ierr = VecSetValues(a_dst,1,&ki,&v,INSERT_VALUES);CHKERRQ(ierr);
        }
      }
    }//dit
    ierr = VecAssemblyBegin(a_dst );  CHKERRQ(ierr);
    ierr = VecAssemblyEnd(  a_dst );  CHKERRQ(ierr);
    return 0;
  }

  int getStencilComponents(vector< EBIndex<CELL> >             & a_dstVoFs, 
                           vector< LocalStencil<CELL, REAL> >  & a_wstencil,
                           Stencil<Real>                       & regStencilInterior)
  {
    return 0;
  }

  int petscToChombo(EBLevelBoxData<CELL, 1> & a_deviceDst,
                    const Vec               & a_src)
                    
  {
  PetscErrorCode ierr;
  const PetscScalar *arr;
  ierr = VecGetArrayRead(a_src,&arr);  CHKERRQ(ierr);

  //copy data to host.  The PeTsC team has promised to make this
  //step unnecessary by providing the pointer to the vector that lives
  //on the device
  hostfactreal_t factory(m_graphs);
  LevelData<EBHostData<CELL, , 1> > hostDst(m_grids, 1, a_deviceSrc.ghost(), factory);
    
  DataIterator dit = m_grids.dataIterator();
  for (int ibox=0;ibox< dit.size(); ibox++)
  {
    auto& hostfab =   hostSrc[dit[ibox]];
    auto& graph = (*m_graphs)[dit[ibox]];
    auto& grid  =     m_grids[dit[ibox]];
    Bx grbx = ProtoCh::getProtoBox(grid);
    for(auto bit = grbx.begin(); bit != grbx.end(); ++bit)
    {
      auto pt = *bit;
      auto vofs = m_graph.getVoFs(pt);
      for(int ivof = 0; ivof < vofs.size(); ivof++)
      {
        auto vof = vofs[ivof];
        PetscInt ki = m_gids[dit[ibox]](vof, 0);
        hostfab(vof, 0) = arr[ki];
      }
    }//dit

    ///copy data back to device.   Soon, the pEtSc team will render this unnecessary.
    EBLevelBoxData<CELL, 1>::copyToDevice(a_deviceDst, hostDst);
    return 0;
  }

  
  virtual void formMatrix()
  {
    
    {
      char str[256];
      strcpy (str,"-");
      strcat (str,m_prestring);
#if PETSC_VERSION_GE(3,6,0)
      strcat (str,"pc_gamg_square_graph 20");
#else
      strcat (str,"pc_gamg_square_graph true");
#endif
#if PETSC_VERSION_GE(3,7,0)
      ierr = PetscOptionsInsertString(PETSC_NULL,str);CHKERRQ(ierr);
#else
      ierr = PetscOptionsInsertString(str);CHKERRQ(ierr);
#endif
      PetscInt  ierr;
      Real nil = 1;
      DataIterator dit = m_grids.dataIterator();
      for(int ibox = 0; ibox < dit.size(): ibox++)
      {
        bool periodic = false;
        vector< EBIndex<CELL> >            dstVoFs;
        vector< LocalStencil<CELL, REAL> > wstencil;
        Stencil<Real>    regStencilInterior;
        getStencilComponents(dstVoFs, wstencil, regStencilInterior, dit[ibox]);
        
        auto& graph = (*m_graphs)[dit[ibox]];
        auto& grid  =     m_grids[dit[ibox]];
        Bx grbx = ProtoCh::getProtoBox(grid);
        for(auto bit = grbx.begin(); bit != grbx.end(); ++bit)
        {
          auto pt = *bit;
          auto vofs = graph.getVoFs(pt);
      for(int ivof = 0; ivof < vofs.size() ivof++)
      {
        auto vof = vofs[ivof];
        int irow  = m_gids[dit[ibox]](vof, 0);
        //the chombo version of this has domain bcs meshed in here.
        //the original mark adams code was hard coded for Neumann boundary
        //conditions,
                  VoFStencil stenc;
                  //get the stencil for this point

                  m_op->getVoFStencil(stenc, vof, dit());
                  for(int ivof = 0; ivof < stenc.size(); ivof++)
                    {
                      const VolIndex& stenvof = stenc.vof(ivof);
                      const Real    & stenwei = stenc.weight(ivof);
                      int icol = m_gids[dit()](stenvof.gridIndex(), 0);  //this will also break with multicells
                      PetscInt irowpet = irow;
                      PetscInt icolpet = icol;
                      ierr = MatSetValues(a_mat,1,&irowpet,1,&icolpet,&stenwei,INSERT_VALUES);
                      CHKERRQ(ierr);
                    }//end loop over stencil
                } //end loop over vofs in cell
            }// end if(notcovered)
        } //box iterator over box because petscsolver is not too smart
    }//data iterator

  ierr = MatAssemblyBegin(a_mat,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(a_mat,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

  return 0;

      
    }
    //HERE use strings to get stencil for this particular operator.
  }
  virtual void setupSolver()
  {
    // create solvers
    PetscBool ism = PETSC_FALSE;
#if PETSC_VERSION_GE(3,7,0)
    PetscOptionsGetBool(PETSC_NULL,m_prestring,"-ksp_monitor",&ism,PETSC_NULL);
#else 
    PetscOptionsGetBool(m_prestring,"-ksp_monitor",&ism,PETSC_NULL);
#endif
    // create the KSP so that we can set KSP parameters
    KSPCreate( wcomm, &m_ksp );
    if ( strlen(m_prestring) > 0 )
    {
      ierr = KSPSetOptionsPrefix( m_ksp, m_prestring );    CHKERRQ(ierr);
    }
    ierr = KSPSetFromOptions(m_ksp);CHKERRQ(ierr);
    if (ism)
    {
      ierr = KSPMonitorSet(m_ksp,ksp_monitor_pout,PETSC_NULL,PETSC_NULL); CHKERRQ(ierr);
    }
#if PETSC_VERSION_GE(3,5,0)
    ierr = KSPSetOperators(m_ksp,m_mat,m_mat);CHKERRQ(ierr);
#else
    ierr = KSPSetOperators(m_ksp,m_mat,m_mat,SAME_NONZERO_PATTERN);CHKERRQ(ierr);
#endif
    ierr = KSPSetInitialGuessNonzero(m_ksp, m_nz_init_guess ? PETSC_TRUE : PETSC_FALSE );CHKERRQ(ierr);

    //there is some stuff here about blocksize but it appears to be optional.
  }
  
  virtual int getNNZPerRow() const
  {
    return 1000;
  }
  
  PetscInt createMatrixAndVectors()
  {
    // create matrix
    PetscInt nnzrow = getNNZPerRow();
    PetscInt *d_nnz=PETSC_NULL, *o_nnz=PETSC_NULL;
#ifdef CH_MPI
    MPI_Comm wcomm = Chombo_MPI::comm;
#else
    MPI_Comm wcomm = PETSC_COMM_SELF;
#endif
    int nc = 1;
    PetscInt ierr;
    ierr = MatCreate(wcomm,&m_mat);CHKERRQ(ierr);
    ierr = MatSetSizes(m_mat,m_NN,m_NN,PETSC_DECIDE,PETSC_DECIDE);CHKERRQ(ierr);
    ierr = MatSetBlockSize(m_mat,nc);CHKERRQ(ierr);
    ierr = MatSetType(m_mat,MATAIJ);CHKERRQ(ierr);
    ierr = MatSeqAIJSetPreallocation(m_mat,nnzrow, d_nnz);CHKERRQ(ierr);
    ierr = MatMPIAIJSetPreallocation(m_mat,nnzrow, d_nnz, nnzrow/2, o_nnz);CHKERRQ(ierr);
    ierr = MatSetOption(m_mat,MAT_NEW_NONZERO_ALLOCATION_ERR,PETSC_FALSE) ;CHKERRQ(ierr);
    ierr = MatSetFromOptions( m_mat ); CHKERRQ(ierr);
      
    if ( d_nnz )
    {
      ierr = PetscFree( d_nnz );  CHKERRQ(ierr);
      ierr = PetscFree( o_nnz );  CHKERRQ(ierr);
    }

    // create vectors
    ierr = VecCreate( wcomm, &m_bb ); CHKERRQ(ierr);
    ierr = VecSetFromOptions( m_bb ); CHKERRQ(ierr);
    ierr = VecSetSizes( m_bb, m_NN, PETSC_DECIDE ); CHKERRQ(ierr);
    ierr = VecDuplicate( m_bb, &m_rr ); CHKERRQ(ierr);
    ierr = VecDuplicate( m_bb, &m_xx ); CHKERRQ(ierr);

    return 0;
  }
  
  void defineGIDS()
  {
    m_gids.define(m_grids, 1, IntVect::Zero, hostfactint_t(m_graphs));
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
    gather(m_numPtsAllProc, numPtsThisProc, 0);
    broadcast(m_numptsAllProc, 0);
    int startgid = 0;
    for(int iproc = 0; iproc < procID(); iproc++)
    {
      startgid += numptsAllProc[iproc];
    }
#else
    int startgid = 0;
#endif
    m_gid0     = startgid;
    m_NN       = numpointsThisProc;
    
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
  //index space--allows us to to matrix to solution space
  LevelData<EBHostData<CELL, int, 1> > m_gids;
  PetscInt m_gid0;
  PetscInt m_NN;


  shared_ptr<GeometryService<order> >  m_geoserv;
  shared_ptr< dictionary_t>            m_ebdictionary;
  shared_ptr<LevelData<EBGraph> >      m_graphs;
  DisjointBoxLayout                    m_grids;
  Box                                  m_domain;
  string m_stencilName;
  string m_domainBCName[2*DIM];
  string m_ebbcName;

  Vector<int> m_numPtsAllProc;
  Point   m_ptghost;
  IntVect m_ivghost;
  Mat m_mat;
  void *m_ctx; // pointer for nonlnear solver call backs

  Vec m_xx, m_rr, m_bb;
  SNES m_snes;
  KSP m_ksp;
  PetscInt m_defined;
};

#endif   //if petsc is defined

#endif 
