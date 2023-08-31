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
#include "petscmat.h"
#include "Chombo_EBChombo.H"
#include "Chombo_EBLevelBoxData.H"

/// Stuff petsc frameworks have in common
/**
*/
template <int order>
class EBPetscSolver
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
  typedef CH4_Data_Choreography::DistributedData<EBGraph>                    graph_distrib_t;


  EBPetscSolver(const shared_ptr<GeometryService<order> >  & a_geoserv,
                const shared_ptr<dictionary_t           >  & a_ebdictionary,
                const shared_ptr<graph_distrib_t>          & a_graphs,
                const Chombo4::DisjointBoxLayout           & a_grids,
                const Chombo4::Box                         & a_domain,
                string                                       a_stencilName,
                string                                       a_domainBCName[2*DIM],
                string                                       a_ebbcName,
                Real                                         a_dx,
                Real                                         a_alpha,
                Real                                         a_beta,
                Point                                        a_ghost,
                bool                                         a_printStuff = false)

  {
    m_stencilName = a_stencilName;
    for(int iface = 0; iface < 2*DIM; iface++)
    {
      m_domainBCName[iface] = a_domainBCName[iface];
    }
    m_ebbcName = a_ebbcName;
  
    baseDefine(a_geoserv,           
               a_ebdictionary,      
               a_graphs,            
               a_grids,             
               a_domain,            
               a_dx,                
               a_alpha,             
               a_beta,              
               a_ghost,             
               a_printStuff);
  

  }
  virtual ~EBPetscSolver()
  {
  }

  void 
  baseDefine(const shared_ptr<GeometryService<order> >  & a_geoserv,
                    const shared_ptr<dictionary_t           >  & a_ebdictionary,
                    const shared_ptr<graph_distrib_t>          & a_graphs,
                    const Chombo4::DisjointBoxLayout           & a_grids,
                    const Chombo4::Box                         & a_domain,
                    Real a_dx,   Real a_alpha, Real a_beta, Point a_ghost,
                    bool a_printStuff = false)
  {
    CH_TIME("EBPetscSolver::baseDefine");
    m_geoserv      =    a_geoserv;
    m_ebdictionary =    a_ebdictionary;
    m_graphs       =    a_graphs;
    m_grids        =    a_grids;
    m_domain       =    a_domain;
    m_dx           =    a_dx;
    m_alpha        =    a_alpha;
    m_beta         =    a_beta;
    m_prestring[0] = '\0';
    
    m_ivghost      =    ProtoCh::getIntVect(a_ghost);
    m_ptghost      =    a_ghost;
    
    if(a_printStuff)
    {
      Chombo4::pout() << "BaseEBPetscSolver: create map of locations in space to matrix row." << endl;
    }
    
    defineGIDS(a_printStuff);

    if(a_printStuff)
    {
      Chombo4::pout() << "BaseEBPetscSolver:: create space for a matrix and the necessary vectors.." << endl;
    }

    createMatrixAndVectors(a_printStuff);

    if(a_printStuff)
    {
      Chombo4::pout() << "BaseEBPetscSolver:: put actual values into the matrix" << endl;
    }
    //

    ///this is where the problem-dependent stuff lives
    if(a_printStuff)
    {
      Chombo4::pout() << "BaseEBPetscsolver:calling setupSolver" << endl;
    }
    setupSolver(a_printStuff);

  }


  ///apply the operator to get Y = A X
  virtual
  PetscErrorCode
  applyOp(EBLevelBoxData<CELL, 1>       & a_lph,
          const EBLevelBoxData<CELL, 1> & a_phi,
          bool a_printStuff = false)
  {
    CH_TIME("BaseEBPetscSolver::applyOp");

    VecZeroEntries(m_xx);
    VecZeroEntries(m_bb);

    chomboToPetsc(m_xx, a_phi);

//    if(a_printStuff)
//    {
//      Chombo4::pout() << "EBPetscSolver::applyOp before mult: m_xx = " << endl;
//      VecView(m_xx , PETSC_VIEWER_STDOUT_WORLD);
//      Chombo4::pout() << "EBPetscSolver::applyOp before mult : m_mat = " << endl;
//      MatView(m_mat, PETSC_VIEWER_STDOUT_WORLD);
//    }

    PetscErrorCode err = MatMult(m_mat, m_xx, m_bb);

//      Chombo4::pout() << "EBPetscSolver::applyOp after mult: m_xx = " << endl;
//    VecView(m_xx , PETSC_VIEWER_STDOUT_WORLD);
//      Chombo4::pout() << "EBPetscSolver::applyOp after mult: m_bb = " << endl;
//    VecView(m_bb, PETSC_VIEWER_STDOUT_WORLD);

    petscToChombo(a_lph, m_bb);
    return err;
  }
  
  ///solve for X in AX = rhs
  virtual
  PetscErrorCode
  solve(EBLevelBoxData<CELL, 1>       & a_phi,
        const EBLevelBoxData<CELL, 1> & a_rhs,
        bool a_printStuff = false)
  {

    
    CH_TIME("BaseEBPetscSolver::solve");
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
    //unnessary.   The PETSc team, however, tells me that the communication between
    //host and device will unnecessary soon so I left the barrier in the code.
    MPI_Barrier(wcomm);
#endif

    {
      CH_TIME("actual ksp solve");
      if(a_printStuff)
      {
        PetscInt nRows = 4586;
        PetscInt nCols = 4586;
        MatGetSize(m_mat, &nRows, &nCols);
        Chombo4::pout() << "EBPetscSolver::solve: matrix size = ( " << nRows << " , " << nCols << " )" << endl;
        PetscInt bbSize = 4586;
        PetscInt xxSize = 4586;
        VecGetSize(m_bb, &bbSize);
        VecGetSize(m_xx, &xxSize);
        Chombo4::pout() << "EBPetscSolver::solve: bb size = " << bbSize << endl;
        Chombo4::pout() << "EBPetscSolver::solve: xx size = " << xxSize << endl;
        
      }
      ierr = KSPSolve( m_ksp, m_bb, m_xx );CHKERRQ(ierr);
    }

    // put solution into output
    ierr = petscToChombo( a_phi, m_xx );CHKERRQ(ierr);
    a_phi.exchange();
    return 0;
  }

  void setInitialGuessNonzero( bool b = true )
  {
    m_nz_init_guess = b;
  }

  
protected:
  int chomboToPetsc(Vec                           & a_dst,
                    const EBLevelBoxData<CELL, 1> & a_deviceSrc)
  {
    CH_TIME("BaseEBPetscSolver::chomboToPetsc");
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
    CH_TIME("BaseEBPetscSolver::petscToChombo");
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

  
  static PetscErrorCode ksp_monitor_pout(KSP ksp, PetscInt it, PetscReal rnorm  ,void *ctx)
  {
    using Chombo4::pout;
    Chombo4::pout() << "      KSP:: iteration = " << it << " residual norm = " <<  rnorm << std::endl;
    return 0;
  }
  
  PetscInt setupSolver(bool a_printStuff = false)
  {
    CH_TIME("BaseEBPetscSolver::setupSolver");
    // create solvers
    PetscInt ierr = 4586;
    ierr = this->formMatrix(a_printStuff); CHKERRQ(ierr);
    
#ifdef CH_MPI
    MPI_Comm wcomm = Chombo_MPI::comm;
#else
    MPI_Comm wcomm = PETSC_COMM_SELF;
#endif
    PetscInt irr = 4586;
    // create the KSP so that we can set KSP parameters
    ierr = KSPCreate( wcomm, &m_ksp );;CHKERRQ(ierr);

    ierr = KSPSetOperators(m_ksp,m_mat,m_mat);CHKERRQ(ierr);

    return 0;
  }
  
  PetscInt createMatrixAndVectors(bool a_printStuff = false)
  {
    CH_TIME("BaseEBPetscSolver::createMatrixAndVectors");
    // create matrix
    //apparently we are in supportnnzexact land
    PetscInt nnzrow = 0;
#ifdef CH_MPI
    MPI_Comm wcomm = Chombo_MPI::comm;
#else
    MPI_Comm wcomm = PETSC_COMM_SELF;
#endif
//    int nc = 1;
    PetscInt ierr;
    int nnzperrow = 1000;
    ierr = MatCreate(wcomm,&m_mat);CHKERRQ(ierr);
    ierr = MatSetSizes(m_mat,m_numPtsThisProc,m_numPtsThisProc,m_numPtsAllProcs,m_numPtsAllProcs);CHKERRQ(ierr);
    MatSetFromOptions(m_mat);
#ifdef CH_MPI    
    MatMPIAIJSetPreallocation(m_mat, nnzperrow, NULL, nnzperrow,     NULL);
#else
    MatSeqAIJSetPreallocation(m_mat, nnzperrow, NULL);
#endif
                              
    VecCreate(wcomm, &m_bb);
    VecSetSizes(m_bb, PETSC_DECIDE, m_numPtsAllProcs);
    VecSetFromOptions(m_bb);
    ierr = VecDuplicate( m_bb, &m_xx);CHKERRQ(ierr);
    ierr = VecDuplicate( m_bb, &m_rr ); CHKERRQ(ierr);


    return 0;
  }
  
  void defineGIDS(bool a_printStuff = false)
  {
    CH_TIME("BaseEBPetscSolver::defineGIDS");
    using Chombo4::DataIterator;
    IntVect gidghost = 6*IntVect::Unit;
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

    Chombo4::pout() << "BaseEBPetscSolver::defineGids: numptsThisProc = " << numPtsThisProc << endl;
    //decide which location maps to the first one of this proc

#ifdef CH_MPI
    std::vector<int> numPtsAllProcs(CH4_SPMD::numProc());

    if(a_printStuff)
    {
      Chombo4::pout() << "BaseEBPetscSolver::defineGids: doing gather numptsThisProc = " << numPtsThisProc << endl;
    }
    MPI_Gather(&numPtsThisProc, 1, MPI_INT, &numPtsAllProcs[0], 1, MPI_INT, 0, Chombo_MPI::comm);
    
    if(a_printStuff)
    {
      Chombo4::pout() << "BaseEBPetscSolver::defineGids: after gather " << endl;
      if(CH4_SPMD::procID() == 0)
      {
        Chombo4::pout() << "BaseEBPetscSolver::defineGids: numPtsAllProcs = " << endl;
        for(int iproc = 0; iproc < numPtsAllProcs.size(); iproc++)
        {
          Chombo4::pout() << "EBPetscSolver::defineGids: iproc = " << iproc << ", numpts[iproc] = " << numPtsAllProcs[iproc] << endl;
        }
      }
    }

    if(a_printStuff)
    {
      Chombo4::pout() << "BaseEBPetscSolver::defineGids:before broadcast " << endl;
    }
    MPI_Bcast(numPtsAllProcs.data(), numPtsAllProcs.size(), MPI_INT, 0, Chombo_MPI::comm);

    if(a_printStuff)
    {
      Chombo4::pout() << "BaseEBPetscSolver::defineGids:after broadcast " << endl;
      for(int iproc = 0; iproc < numPtsAllProcs.size(); iproc++)
      {
        Chombo4::pout() << "BaseEBPetscSolver::defineGids: iproc = " << iproc << ", numpts[iproc] = " << numPtsAllProcs[iproc] << endl;
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
    Chombo4::pout()  << "BaseEBPetscSolver: procID = " << procID() << ", numPtsThisProc = " << numPtsThisProc << ", startgid = " << startgid <<", totalNumPts = " << totalNumPts <<  endl;
#else
    int startgid = 0;
    PetscInt totalNumPts = numPtsThisProc;
#endif
    m_gid0     = startgid;
    m_numPtsThisProc       = numPtsThisProc;
    m_numPtsAllProcs       = totalNumPts;
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
      Chombo4::pout() << "EBPetscSolver::defineGids: calling exchange on gids" << endl;
    }
    m_gids.exchange();
    if(a_printStuff)
    {
      Chombo4::pout() << "EBPetscSolver::defineGids: leaving" << endl;
    }
  }
  //index space--allows us to to matrix to solution space
  host_distrib_int_t  m_gids;
  PetscInt m_gid0;
  PetscInt m_numPtsThisProc;
  PetscInt m_numPtsAllProcs;
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


  ///
  void getLocalRegularStencil(Proto::LocalStencil<CELL, Real> & a_localSten,
                              const EBIndex<CELL>             & a_vof,
                              const Proto::Stencil<Real >     & a_regSten)
  {
    a_localSten.m_entries.resize(0);
    auto startpt = a_vof.m_pt;
    auto coefs   = a_regSten.coefs();
    auto offsets = a_regSten.offsets();
    Bx domainbx = ProtoCh::getProtoBox(this->m_domain);
    for(int isten =  0; isten < coefs.size(); isten++)
    {
      Point stenpt = startpt + offsets[isten];
      if(domainbx.contains(stenpt))
      {
        EBIndex<CELL> stenvof;
        stenvof.defineAsRegular(stenpt);
        a_localSten.add(stenvof, coefs[isten]);
      }
      else
      {
        //compensate for neumann or dirichlet bcs
        //here is where we diverge from canon.
        //The original just ignored stuff outside the domain.
        for(int idir = 0; idir < DIM; idir++)
        {
          for(SideIterator sit; sit.ok(); ++sit)
          {
            int isign = sign(sit());
            Point justOutside = startpt + isign*Point::Basis(idir);
            if(stenpt == justOutside)
            {
              int index = ebp_index(idir, sit());
              string bcstr = this->m_domainBCName[index];
              if(bcstr == StencilNames::Dirichlet)
              {
                a_localSten.add(a_vof, -coefs[isten]);
              }
              else if(bcstr == StencilNames::Neumann)
              {
                a_localSten.add(a_vof,  coefs[isten]);
              }
            }
          }
        }
      }
    }
  }
  int getStencilComponents(vector< EBIndex<CELL> >             & a_dstVoFs, 
                           vector< Proto::LocalStencil<CELL, Real> >  & a_wstencil,
                           Proto::Stencil<Real>                       & a_regStencilInterior,
                           PointSet                            & a_allIrregPoints,
                           const Chombo4::DataIndex                     & a_dit)
  {
    CH_TIME("EBPetscSolver::getStencilComponents");
    typedef Proto::EBStencilArchive<CELL, CELL, order, Real> archive_t;
    using Chombo4::Box;
    typedef IndexedMoments<DIM  , order> IndMomDIM;
    typedef IndexedMoments<DIM-1, order> IndMomSDMinOne;
    typedef HostIrregData<CELL    ,  IndMomDIM , 1>  VoluData;
    typedef HostIrregData<BOUNDARY,  IndMomDIM , 1>  EBFaData;
    typedef HostIrregData<XFACE, IndMomSDMinOne, 1>  XFacData;
    typedef HostIrregData<YFACE, IndMomSDMinOne, 1>  YFacData;
    typedef HostIrregData<ZFACE, IndMomSDMinOne, 1>  ZFacData;
    auto voldatpld = this->m_geoserv->getVoluData(  this->m_domain);
    auto ebfdatpld = this->m_geoserv->getEBFaceData(this->m_domain);
    auto xfadatpld = this->m_geoserv->getXFaceData( this->m_domain);
    auto yfadatpld = this->m_geoserv->getYFaceData( this->m_domain);
    auto zfadatpld = this->m_geoserv->getZFaceData( this->m_domain);

    auto& voludata =(*voldatpld)[a_dit];
    auto& ebfadata =(*ebfdatpld)[a_dit];
    auto& xfacdata =(*xfadatpld)[a_dit];
    auto& yfacdata =(*yfadatpld)[a_dit];
    auto& zfacdata =(*zfadatpld)[a_dit];
    auto& graph    =(*(this->m_graphs))[a_dit];
    Box valid = this->m_grids[a_dit];
    Bx  valbx = ProtoCh::getProtoBox(valid);
    vector<Proto::Stencil<Real> >    regStencilBCS[2*DIM];
    vector<Bx >               BCApplyBoxes[2*DIM];
    bool                      bcsOnly = false;
    Bx                        regApplyBox;
    bool                      irregOnly = false;
    Bx domainbx = ProtoCh::getProtoBox(this->m_domain);
    archive_t::getStencil(a_dstVoFs, a_wstencil, a_regStencilInterior, 
                          regStencilBCS, BCApplyBoxes, regApplyBox, bcsOnly, irregOnly,
                          this->m_stencilName,  this->m_domainBCName, this->m_ebbcName,
                          valbx, valbx,  domainbx, domainbx,
                          this->m_ptghost,  this->m_ptghost,    
                          graph,  graph,
                          voludata,  ebfadata,    
                          xfacdata,  yfacdata, zfacdata,    
                          this->m_dx, false, Point::Zeros());

    a_allIrregPoints.makeEmpty();
    for(int ipt = 0; ipt < a_dstVoFs.size(); ipt++)
    {
      a_allIrregPoints |= a_dstVoFs[ipt].m_pt;
    }

    return 0;
  }


  virtual PetscInt formMatrix(bool a_printStuff = false)
  {
    
    PetscInt ierr = 4586;
    CH_TIME("EBPetscSolver::formMatrix");
    {
      Chombo4::DataIterator dit = this->m_grids.dataIterator();
      for(int ibox = 0; ibox < dit.size(); ibox++)
      {
        vector< EBIndex<CELL> >            dstVoFs;
        vector< Proto::LocalStencil<CELL, Real> > wstencil;
        Stencil<Real>    regStencilInterior;
        PointSet irregSet;
        //use strings to get stencil for this particular operator.
        getStencilComponents(dstVoFs, wstencil, regStencilInterior, irregSet, dit[ibox]);
        
        auto& graph = (*this->m_graphs)[dit[ibox]];
        auto  grid  =     this->m_grids[dit[ibox]];
        Bx grbx = ProtoCh::getProtoBox(grid);
        //add stencils for cells that are in the regular stencil regime to the irregular
        //set.   This includes adjusting for boundary condtions
        for(auto bit = grbx.begin(); bit != grbx.end(); ++bit)
        {
          auto pt = *bit;
          if(graph.isCovered(pt))
          {
            //do nothing
          }
          else if(!irregSet.contains(pt))
          {
            EBIndex<CELL> vof;
            vof.defineAsRegular(pt);
            LocalStencil<CELL, Real> regSten;
            getLocalRegularStencil(regSten,  vof, regStencilInterior);
            dstVoFs.push_back(vof);
            wstencil.push_back(regSten);
          }
        }
        if(a_printStuff)
        {
          int nrow = 4586;
          int ncol = 4586;
          MatGetSize(this->m_mat, &nrow, &ncol);
          Chombo4::pout() << "EBPetcSolver::formMatrix: matrix size = ( " << nrow << " . " << ncol << ")" << endl;
        }
        //now create the matrix
        for(int idst = 0; idst < dstVoFs.size(); idst++)
        {
          auto vof    =   dstVoFs[idst];
          auto& stenc =  wstencil[idst];
          int irow = this->m_gids[dit[ibox]](vof, 0); 
          for(int jsrc = 0; jsrc < stenc.size(); jsrc++)
          {
            auto& stenvof = stenc.m_entries[jsrc].m_vof;
            auto& stenwgt = stenc.m_entries[jsrc].m_weight;
            Real weight = (this->m_beta)*stenwgt;
            if(stenvof == vof)
            {
              weight += this->m_alpha;
            }
            int icol = (this->m_gids)[dit[ibox]](stenvof, 0);  
            PetscInt irowpet = irow;
            PetscInt icolpet = icol;
            ierr = MatSetValues(this->m_mat,1,&irowpet,1,&icolpet,&weight,INSERT_VALUES);
            CHKERRQ(ierr);
          } 
        }
      }
      ierr = MatAssemblyBegin(this->m_mat,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
      ierr = MatAssemblyEnd(  this->m_mat,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
    }
    return 0;
  }
};

#endif

#endif 
