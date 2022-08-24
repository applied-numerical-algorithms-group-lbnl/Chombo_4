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
#include "petscvec.h"
#include "petscksp.h"
#include "petscviewer.h"
#include "Chombo_EBChombo.H"
#include "Chombo_EBLevelBoxData.H"

namespace hoeb
{
/// Framework to solve an elliptic equation using PETSc (PeTSC? PETsC? PETsC?)
/**
   This is an an adaptation to the new reality of Mark Adams' PetscSolver framework.
   The current reality is we need dharshi's laplacian  and I c
   dtg 
   8-24-2022
*/
  template <int order>
  class Hoeb_PetscSolver
  {
  public:
    typedef GraphConstructorFactory< EBBoxData<CELL, int, 1> >  devifactint_t;
    typedef GraphConstructorFactory<EBHostData<CELL, int, 1> >  host_fact_int_t;
    typedef GraphConstructorFactory<EBHostData<CELL,Real, 1> >  host_fact_real_t;
    typedef EBDictionary<order, Real, CELL, CELL>                dictionary_t;
    typedef CH4_Data_Choreography::DistributedData<EBBoxData< CELL, Real, 1> > devi_distrib_real_t;
    typedef CH4_Data_Choreography::DistributedData<EBHostData<CELL, Real, 1> > host_distrib_real_t;
    typedef CH4_Data_Choreography::DistributedData<EBBoxData< CELL, int , 1> > devi_distrib_int_t;
    typedef CH4_Data_Choreography::DistributedData<EBHostData<CELL, int , 1> > host_distrib_int_t;
    typedef CH4_Data_Choreography::DistributedData<EBGraph>                   graph_distrib_t;

    Hoeb_PetscSolver(const shared_ptr<GeometryService<order> >  & a_geoserv,
                     const shared_ptr<dictionary_t           >  & a_ebdictionary,
                     const shared_ptr<graph_distrib_t>          & a_graphs,
                     const Chombo4::DisjointBoxLayout           & a_grids,
                     const Chombo4::Box                         & a_domain,
                     string a_stencilName,
                     string a_domainBCName[2*DIM],
                     string a_ebbcName,
                     Real a_dx,   Real a_alpha, Real a_beta, Point a_ghost,
                     bool a_printStuff = false)
    {
      CH_TIME("Hoeb_PetscSolver::Hoeb_PetscSolver");
      m_geoserv      =    a_geoserv;
      m_ebdictionary =    a_ebdictionary;
      m_graphs       =    a_graphs;
      m_grids        =    a_grids;
      m_domain       =    a_domain;
      m_stencilName  =    a_stencilName;
      m_dx           =    a_dx;
      m_alpha        =    a_alpha;
      m_beta         =    a_beta;
      m_prestring[0] = '\0';
      for(int iface = 0; iface < 2*DIM; iface++)
      {
        m_domainBCName[iface] = a_domainBCName[iface];
      }
      m_ebbcName = a_ebbcName;
    
      m_ivghost      =    ProtoCh::getIntVect(a_ghost);
      m_ptghost      =    a_ghost;
    
      if(a_printStuff)
      {
        Chombo4::pout() << "Hoeb_PetscSolver::Hoeb_PetscSolver: create map of locations in space to matrix row." << endl;
      }
    
      defineGIDS(a_printStuff);

      if(a_printStuff)
      {
        Chombo4::pout() << "Hoeb_PetscSolver::Hoeb_PetscSolver: create space for a matrix and the necessary vectors.." << endl;
      }

      createMatrixAndVectors(a_printStuff);

      if(a_printStuff)
      {
        Chombo4::pout() << "Hoeb_PetscSolver::Hoeb_PetscSolver:put actual values into the matrix" << endl;
      }
      //
      formMatrix(a_printStuff);

      if(a_printStuff)
      {
        Chombo4::pout() << "Hoeb_PetscSolver::Hoeb_PetscSolver:calling setupSolver" << endl;
      }
      setupSolver(a_printStuff);

      if(a_printStuff)
      {
        Chombo4::pout() << "Hoeb_PetscSolver::Hoeb_PetscSolver:leaving" << endl;
      }
    }

    virtual
    PetscErrorCode
    solve(EBLevelBoxData<CELL, 1>       & a_phi,
          const EBLevelBoxData<CELL, 1> & a_rhs,
          bool a_printStuff = false)
    {

    
      CH_TIME("Hoeb_PetscSolver::solve");
      auto dbl = m_grids;

      PetscErrorCode ierr;
    
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
      MPI_Comm wcomm = Chombo_MPI::comm;
      //with all the host to device and device to host stuff going on, this is probably
      //unnessary.   The peTsC team, however, tells me that the communication between
      //host and device will unnecessary soon so I left the barrier in the code.
      MPI_Barrier(wcomm);
#endif
      //begin debug
      //std::cout << "debug::just before KSPSolve for proc = " << CH4_SPMD::procID() << std::endl;
      //Chombo4::pout() << "debug::just before KSPSolve for proc = " << CH4_SPMD::procID() << std::endl;
      //MatView(m_mat, PETSC_VIEWER_STDOUT_WORLD);

      //std::cout << "b vector" << std::endl;
      //VecView(m_bb , PETSC_VIEWER_STDOUT_WORLD);
      //std::cout << "x vector" << std::endl;
      //VecView(m_xx  , PETSC_VIEWER_STDOUT_WORLD);
      //end debug
      {
        CH_TIME("actual ksp solve");
        ierr = KSPSolve( m_ksp, m_bb, m_xx );CHKERRQ(ierr);
      }
      //begin debug
      //std::cout << "debug::just after  KSPSolve for proc = " << CH4_SPMD::procID() << std::endl;
      //std::cout << "x vector" << std::endl;
      //VecView(m_xx  , PETSC_VIEWER_STDOUT_WORLD);
      //end debug

      // put solution into output
      ierr = petscToChombo( a_phi, m_xx );CHKERRQ(ierr);
      a_phi.exchange();
      return 0;
    }

    void setInitialGuessNonzero( bool b = true )
    {
      m_nz_init_guess = b;
    }

  
  private:
    int chomboToPetsc(Vec                           & a_dst,
                      const EBLevelBoxData<CELL, 1> & a_deviceSrc)
    {
      CH_TIME("Hoeb_PetscSolver::chomboToPetsc");
      PetscErrorCode ierr;
      //copy data to host.
      host_fact_real_t factory(m_graphs);
      host_distrib_real_t  hostSrc(m_grids,  a_deviceSrc.ghostVect(), factory);
      EBLevelBoxData<CELL, 1>::copyToHost(hostSrc, a_deviceSrc);
    
      Chombo4::DataIterator dit = m_grids.dataIterator();
      for (int ibox=0;ibox< dit.size(); ibox++)
      {
        auto& hostfab =   hostSrc[dit[ibox]];
        auto& graph = (*m_graphs)[dit[ibox]];
        auto  grid  =     m_grids[dit[ibox]];
        Bx grbx = ProtoCh::getProtoBox(grid);
        for(auto bit = grbx.begin(); bit != grbx.end(); ++bit)
        {
          auto pt = *bit;
          auto vofs = graph.getVoFs(pt);
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

    int petscToChombo(EBLevelBoxData<CELL, 1> & a_deviceDst,
                      const Vec               & a_src)
    {
      CH_TIME("Hoeb_PetscSolver::petscToChombo");
      using Chombo4::DataIterator;
      PetscErrorCode ierr;
      const PetscScalar *arr;
      ierr = VecGetArrayRead(a_src,&arr);  CHKERRQ(ierr);


      PetscInt localStart,  localEnd, localSize;
      VecGetOwnershipRange( a_src ,&localStart,&localEnd);
      VecGetLocalSize(a_src ,&localSize);

      host_fact_real_t factory(m_graphs);
      host_distrib_real_t hostDst(m_grids, a_deviceDst.ghostVect(), factory);

      DataIterator dit = m_grids.dataIterator();
      for (int ibox=0;ibox< dit.size(); ibox++)
      {
        auto& hostfab =   hostDst[dit[ibox]];
        auto& graph = (*m_graphs)[dit[ibox]];
        auto  grid  =     m_grids[dit[ibox]];
        Bx grbx = ProtoCh::getProtoBox(grid);
        for(auto bit = grbx.begin(); bit != grbx.end(); ++bit)
        {
          auto pt = *bit;
          auto vofs = graph.getVoFs(pt);
          for(int ivof = 0; ivof < vofs.size(); ivof++)
          {
            auto vof = vofs[ivof];
            PetscInt ki = m_gids[dit[ibox]](vof, 0);
            PetscInt localIndex = ki - localStart;
            hostfab(vof, 0) = arr[localIndex];
	  
          }
        }
      }

      ///copy data back to device. 
      EBLevelBoxData<CELL, 1>::copyToDevice(a_deviceDst, hostDst);


      //end debug
      return 0;
    }

  
    virtual PetscInt formMatrix(bool a_printStuff = false)
    {
    
      CH_TIME("Hoeb_PetscSolver::formMatrix");
      {
        using Chombo4::DataIterator;
        char str[256];
        strcpy (str,"-");
        strcat (str,m_prestring);
#if PETSC_VERSION_GE(3,6,0)
        strcat (str,"pc_gamg_square_graph 20");
#else
        strcat (str,"pc_gamg_square_graph true");
#endif
        PetscInt  ierr;
#if PETSC_VERSION_GE(3,7,0)
        ierr = PetscOptionsInsertString(PETSC_NULL,str);CHKERRQ(ierr);
#else
        ierr = PetscOptionsInsertString(str);CHKERRQ(ierr);
#endif
        DataIterator dit = m_grids.dataIterator();
        for(int ibox = 0; ibox < dit.size(); ibox++)
        {
          Bx srcValid = ProtoCh::getProtoBox(m_grids[dit[ibox]]);
          Bx dstValid = ProtoCh::getProtoBox(m_grids[dit[ibox]]);
          Bx srcDom = ProtoCh::getProtoBox(m_domain);
          Bx dstDom = ProtoCh::getProtoBox(m_domain);

          vector< EBIndex<CELL> >            dstVoFs;
          vector< Proto::LocalStencil<CELL, Real> > wstencil;
          Stencil<Real>    regStencilInterior;
          bool needDiag = false;
          hoeb::getHomogeneousDharshiStencil(m_stencilName, m_ebbcName, dstVoFs, wstencil, srcValid, dstValid,
                                             srcDom, dstDom, m_ptghost, m_ptghost, needDiag, m_geoserv, m_grids, m_domain, m_dx, ibox);
                                             
          //now create the matrix
          for(int ivof = 0; ivof < dstVoFs.size(); ivof++)
          {
            auto vof    =   dstVoFs[ivof];
            auto& stenc =  wstencil[ivof];
            int irow = m_gids[dit[ibox]](vof, 0); 
            for(int ivof = 0; ivof < stenc.size(); ivof++)
            {
              auto& stenvof = stenc.m_entries[ivof].m_vof;
              auto& stenwgt = stenc.m_entries[ivof].m_weight;
              Real weight = m_beta*stenwgt;
              if(stenvof == vof)
              {
                weight += m_alpha;
              }
              int icol = m_gids[dit[ibox]](stenvof, 0);  
              PetscInt irowpet = irow;
              PetscInt icolpet = icol;
              ierr = MatSetValues(m_mat,1,&irowpet,1,&icolpet,&weight,INSERT_VALUES);
            } 
          }
        }

        ierr = MatAssemblyBegin(m_mat,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
        ierr = MatAssemblyEnd(m_mat,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

      }
      return 0;
    }
    static PetscErrorCode ksp_monitor_pout(KSP ksp, PetscInt it, PetscReal rnorm  ,void *ctx)
    {
      using Chombo4::pout;
      Chombo4::pout() << "      KSP:: iteration = " << it << " residual norm = " <<  rnorm << std::endl;
      return 0;
    }
    virtual PetscInt setupSolver(bool a_printStuff = false)
    {
      CH_TIME("Hoeb_PetscSolver::setupSolver");
      // create solvers
      PetscBool ism = PETSC_FALSE;
#if PETSC_VERSION_GE(3,7,0)
      PetscOptionsGetBool(PETSC_NULL,m_prestring,"-ksp_monitor",&ism,PETSC_NULL);
#else 
      PetscOptionsGetBool(m_prestring,"-ksp_monitor",&ism,PETSC_NULL);
#endif
#ifdef CH_MPI
      MPI_Comm wcomm = Chombo_MPI::comm;
#else
      MPI_Comm wcomm = PETSC_COMM_SELF;
#endif
      // create the KSP so that we can set KSP parameters
      KSPCreate( wcomm, &m_ksp );
      PetscInt  ierr;
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
      //ierr = KSPSetInitialGuessNonzero(m_ksp, m_nz_init_guess ? PETSC_TRUE : PETSC_FALSE );CHKERRQ(ierr);

      { // coordinates
        PC pc; 
        PetscInt sz,ind,bs,n,m;
#if PETSC_VERSION_LT(3,4,0) & PETSC_VERSION_RELEASE
        const PCType type;
#else
        PCType type;
#endif  
        ierr = KSPGetPC( m_ksp, &pc );     CHKERRQ(ierr);
        ierr = PCGetType( pc, &type );    CHKERRQ(ierr);
        ierr = MatGetBlockSize( m_mat, &bs );               CHKERRQ( ierr );
        if ( strcmp(type,PCGAMG) == 0 && bs > 1 )
        {
          PetscReal    *coords;
          ierr = MatGetLocalSize( m_mat, &m, &n );  CHKERRQ(ierr);
          sz = CH_SPACEDIM*(m/bs);
          ierr = PetscMalloc( (sz+1)*sizeof(PetscReal), &coords ); CHKERRQ(ierr);
          ind  = 0;//just for parity check with sz
          for (Chombo4::DataIterator  dit = m_grids.dataIterator(); dit.ok() ; ++dit )
          {
            Chombo4::Box grid = m_grids[dit()];
            for (Chombo4::BoxIterator bit(grid); bit.ok(); bit.next())
            {
              IntVect iv = bit(); // coordinate in any scaled, shifted, rotated frame.
              for (PetscInt k=0; k< DIM; k++)
              {
                ind = ind + 1;
                coords[ind] = (PetscReal)iv[k];
              }
            }
          }
          CH_assert(ind==sz);
          ierr = PCSetCoordinates( pc, CH_SPACEDIM, sz/CH_SPACEDIM, coords ); CHKERRQ(ierr);
          ierr = PetscFree( coords );  CHKERRQ(ierr);
        }
      }

      // print_memory_line("After AMG set up");
      //there is some stuff here about blocksize but it appears to be optional.
      return 0;
    }
  
    virtual int getNNZPerRow() const
    {
      return 1000;
    }
  
    PetscInt createMatrixAndVectors(bool a_printStuff = false)
    {
      CH_TIME("Hoeb_PetscSolver::createMatrixAndVectors");
      // create matrix
      //apparently we are in supportnnzexact land
      PetscInt nnzrow = 0;
#ifdef CH_MPI
      MPI_Comm wcomm = Chombo_MPI::comm;
#else
      MPI_Comm wcomm = PETSC_COMM_SELF;
#endif
      int nc = 1;
      PetscInt ierr;
      ierr = MatCreate(wcomm,&m_mat);CHKERRQ(ierr);
      ierr = MatSetOptionsPrefix(m_mat,"");CHKERRQ(ierr);

      //Chombo4::Chombo4::pout() << "debug:: matsetsize m_NN = " << m_NN << std::endl;
      ierr = MatSetSizes(m_mat,m_NN,m_NN,PETSC_DECIDE,PETSC_DECIDE);CHKERRQ(ierr);
      ierr = MatSetBlockSize(m_mat,nc);CHKERRQ(ierr);
      ierr = MatSetType(m_mat,MATAIJ);CHKERRQ(ierr);

      ierr = MatSetFromOptions( m_mat ); CHKERRQ(ierr);
      PetscInt *d_nnz=PETSC_NULL, *o_nnz=PETSC_NULL;

      PetscInt nglobal;
      PetscInt NN = m_NN;
#ifdef CH_MPI
      PetscInt nn = NN;
      MPI_Datatype mtype;
      PetscDataTypeToMPIDataType(PETSC_INT,&mtype);
      MPI_Allreduce(&nn,&nglobal,1,mtype,MPI_SUM,wcomm);
#else
      nglobal = NN;
#endif    

      ierr = PetscMalloc( (m_NN+1)*sizeof(PetscInt), &d_nnz ); CHKERRQ(ierr);
      ierr = PetscMalloc( (m_NN+1)*sizeof(PetscInt), &o_nnz ); CHKERRQ(ierr);
      for (PetscInt kk=0;kk<NN;kk++) d_nnz[kk] = o_nnz[kk] = 0;

      // fix bounds
      for (PetscInt kk=0;kk<NN;kk++) 
      {
        if (d_nnz[kk] > NN) d_nnz[kk] = NN;
        if (o_nnz[kk] > nglobal-NN) o_nnz[kk] = nglobal-NN;
      }
      nnzrow = 0;

      //Chombo4::pout() << "debug:: matmpiaijsetpreallocation nnzrow = "  << nnzrow << ", d_nnz = " << d_nnz << ", o_nnz = " << o_nnz << endl;
      ierr = MatSeqAIJSetPreallocation(m_mat,nnzrow, d_nnz);CHKERRQ(ierr);
      ierr = MatMPIAIJSetPreallocation(m_mat,nnzrow, d_nnz, nnzrow/2, o_nnz);CHKERRQ(ierr);
      ierr = MatSetOption(m_mat,MAT_NEW_NONZERO_ALLOCATION_ERR,PETSC_FALSE) ;CHKERRQ(ierr);

      
      if ( d_nnz )
      {
        ierr = PetscFree( d_nnz );  CHKERRQ(ierr);
        ierr = PetscFree( o_nnz );  CHKERRQ(ierr);
      }

      // create vectors
      ierr = MatCreateVecs(m_mat,&m_bb,&m_xx);CHKERRQ(ierr);
      //ierr = VecCreate( wcomm, &m_bb ); CHKERRQ(ierr);
      //ierr = VecSetFromOptions( m_bb ); CHKERRQ(ierr);
      //ierr = VecSetSizes( m_bb, m_NN, PETSC_DECIDE ); CHKERRQ(ierr);
      ierr = VecDuplicate( m_bb, &m_rr ); CHKERRQ(ierr);
      //ierr = VecDuplicate( m_bb, &m_xx ); CHKERRQ(ierr);

      return 0;
    }
  
    void defineGIDS(bool a_printStuff = false)
    {
      CH_TIME("Hoeb_PetscSolver::defineGIDS");
      using Chombo4::DataIterator;
      IntVect gidghost = 4*IntVect::Unit;
      m_gids.define(m_grids,  gidghost, host_fact_int_t(m_graphs));
      DataIterator dit = m_grids.dataIterator();
      //get the number of points on each proc.
      int numPtsThisProc = 0;
      for(int ibox = 0; ibox < dit.size(); ibox++)
      {
        auto graph = (*m_graphs)[dit[ibox]];
        auto grid  =     m_grids[dit[ibox]];
        Bx  grbx = ProtoCh::getProtoBox(grid);
        for(auto bit = grbx.begin(); bit != grbx.end(); ++bit)
        {
          auto vofs = graph.getVoFs(*bit);
          numPtsThisProc += vofs.size();
        }
      }

      //decide which location maps to the first one of this proc

#ifdef CH_MPI
      std::vector<int> numPtsAllProcs(CH4_SPMD::numProc());

      if(a_printStuff)
      {
        Chombo4::pout() << "Hoeb_PetscSolver::defineGids: doing gather numptsThisProc = " << numPtsThisProc << endl;
      }
      MPI_Gather(&numPtsThisProc, 1, MPI_INT, &numPtsAllProcs[0], 1, MPI_INT, 0, Chombo_MPI::comm);
    
      if(a_printStuff)
      {
        Chombo4::pout() << "Hoeb_PetscSolver::defineGids: after gather " << endl;
        if(CH4_SPMD::procID() == 0)
        {
          Chombo4::pout() << "Hoeb_PetscSolver::defineGids: numPtsAllProcs = " << endl;
          for(int iproc = 0; iproc < numPtsAllProcs.size(); iproc++)
          {
            Chombo4::pout() << "Hoeb_PetscSolver::defineGids: iproc = " << iproc << ", numpts[iproc] = " << numPtsAllProcs[iproc] << endl;
          }
        }
      }

      if(a_printStuff)
      {
        Chombo4::pout() << "Hoeb_PetscSolver::defineGids:before broadcast " << endl;
      }
      MPI_Bcast(numPtsAllProcs.data(), numPtsAllProcs.size(), MPI_INT, 0, Chombo_MPI::comm);

      if(a_printStuff)
      {
        Chombo4::pout() << "Hoeb_PetscSolver::defineGids:after broadcast " << endl;
        for(int iproc = 0; iproc < numPtsAllProcs.size(); iproc++)
        {
          Chombo4::pout() << "Hoeb_PetscSolver::defineGids: iproc = " << iproc << ", numpts[iproc] = " << numPtsAllProcs[iproc] << endl;
        }
      }
      PetscInt startgid = 0;
      PetscInt totalNumPts = 0; 
      for(int iproc = 0; iproc < Chombo4::procID(); iproc++)
      {
        startgid +=   numPtsAllProcs[iproc];
      }
      for(int iproc = 0; iproc < Chombo4::numProc(); iproc++)
      {
        totalNumPts += numPtsAllProcs[iproc];
      }
      using Chombo4::pout;
      using Chombo4::procID;
      Chombo4::pout()  << "procID = " << procID() << ", numPtsThisProc = " << numPtsThisProc << ", startgid = " << startgid <<", totalNumPts = " << totalNumPts <<  endl;
#else
      int startgid = 0;
#endif
      m_gid0     = startgid;
      m_NN       = numPtsThisProc;
    
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

      if(a_printStuff)
      {
        Chombo4::pout() << "Hoeb_PetscSolver::defineGids: calling exchange on gids" << endl;
      }
      m_gids.exchange(a_printStuff);
      if(a_printStuff)
      {
        Chombo4::pout() << "Hoeb_PetscSolver::defineGids: leaving" << endl;
      }
    }
    //index space--allows us to to matrix to solution space
    host_distrib_int_t  m_gids;
    PetscInt m_gid0;
    PetscInt m_NN;
    bool m_nz_init_guess;


    shared_ptr<GeometryService<order> >  m_geoserv;
    shared_ptr< dictionary_t>            m_ebdictionary;
    shared_ptr<graph_distrib_t>          m_graphs;
    Chombo4::DisjointBoxLayout           m_grids;
    Chombo4::Box                         m_domain;
    Real                                 m_dx;
  
    string m_stencilName;
    string m_domainBCName[2*DIM];
    string m_ebbcName;

    std::vector<int> m_numPtsAllProc;
    Point   m_ptghost;
    IntVect m_ivghost;
    Mat m_mat;
    void *m_ctx; // pointer for nonlnear solver call backs

    Vec m_xx, m_rr, m_bb;
    SNES m_snes;
    KSP m_ksp;
    PetscInt m_defined;
    char m_prestring[32];
    Real m_alpha, m_beta;
  };
}
#endif

#endif 
