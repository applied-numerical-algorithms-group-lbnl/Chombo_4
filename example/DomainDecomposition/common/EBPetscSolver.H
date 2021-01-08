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
#include "Chombo_EBChombo.H"
#include "Chombo_EBLevelBoxData.H"

/// Framework to solve an elliptic equation using PETSc (PeTSC? PETsC? PETsC?)
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
                const Chombo4::DisjointBoxLayout                    & a_grids,
                const Chombo4::Box                                  & a_domain,
                string a_stencilName,
                string a_domainBCName[2*DIM],
                string a_ebbcName,
                Real a_dx,   Real a_alpha, Real a_beta, Point a_ghost)
  {
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
    //create map of locations in space to matrix row.
    defineGIDS();
    //create space for a matrix and the necessary vectors.
    createMatrixAndVectors();
    //put actual values into the matrix
    formMatrix();
    setupSolver();
  }

  virtual
  PetscErrorCode
  solve(EBLevelBoxData<CELL, 1>       & a_phi,
        const EBLevelBoxData<CELL, 1> & a_rhs )
  {

    
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
    ierr = KSPSolve( m_ksp, m_bb, m_xx );CHKERRQ(ierr);

    // put solution into output
    ierr = petscToChombo( a_phi, m_xx );CHKERRQ(ierr);
    Chombo4::Copier styx;
    styx.exchangeDefine(m_grids, m_ivghost);
    a_phi.exchange(styx);
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
    PetscErrorCode ierr;

    //copy data to host.  The pEtsC team has promised to make this
    //step unnecessary by providing the pointer to the vector that lives
    //on the device
    hostfactreal_t factory(m_graphs);
    LevelData<EBHostData<CELL, Real, 1> > hostSrc(m_grids, 1, a_deviceSrc.ghostVect(), factory);
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

  void getLocalRegularStencil(Proto::LocalStencil<CELL, Real> & a_localSten,
                              const EBIndex<CELL>      & a_vof,
                              const Proto::Stencil<Real >     & a_regSten)
  {
    a_localSten.m_entries.resize(0);
    auto startpt = a_vof.m_pt;
    auto coefs   = a_regSten.coefs();
    auto offsets = a_regSten.offsets();
    Bx domainbx = ProtoCh::getProtoBox(m_domain);
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
              string bcstr = m_domainBCName[index];
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
    typedef Proto::EBStencilArchive<CELL, CELL, order, Real> archive_t;
    using Chombo4::Box;
    typedef IndexedMoments<DIM  , order> IndMomDIM;
    typedef IndexedMoments<DIM-1, order> IndMomSDMinOne;
    typedef HostIrregData<CELL    ,  IndMomDIM , 1>  VoluData;
    typedef HostIrregData<BOUNDARY,  IndMomDIM , 1>  EBFaData;
    typedef HostIrregData<XFACE, IndMomSDMinOne, 1>  XFacData;
    typedef HostIrregData<YFACE, IndMomSDMinOne, 1>  YFacData;
    typedef HostIrregData<ZFACE, IndMomSDMinOne, 1>  ZFacData;
    shared_ptr< LevelData< VoluData >  > voldatpld = m_geoserv->getVoluData(  m_domain);
    shared_ptr< LevelData< EBFaData >  > ebfdatpld = m_geoserv->getEBFaceData(m_domain);
    shared_ptr< LevelData< XFacData >  > xfadatpld = m_geoserv->getXFaceData( m_domain);
    shared_ptr< LevelData< YFacData >  > yfadatpld = m_geoserv->getYFaceData( m_domain);
    shared_ptr< LevelData< ZFacData >  > zfadatpld = m_geoserv->getZFaceData( m_domain);

    auto& voludata =(*voldatpld)[a_dit];
    auto& ebfadata =(*ebfdatpld)[a_dit];
    auto& xfacdata =(*xfadatpld)[a_dit];
    auto& yfacdata =(*yfadatpld)[a_dit];
    auto& zfacdata =(*zfadatpld)[a_dit];
    auto& graph    =(*m_graphs )[a_dit];
    Box valid = m_grids[a_dit];
    Bx  valbx = ProtoCh::getProtoBox(valid);
    vector<Proto::Stencil<Real> >    regStencilBCS[2*DIM];
    vector<Bx >               BCApplyBoxes[2*DIM];
    bool                      bcsOnly = false;
    Bx                        regApplyBox;
    Bx domainbx = ProtoCh::getProtoBox(m_domain);
    archive_t::getStencil(a_dstVoFs, a_wstencil, a_regStencilInterior, 
                          regStencilBCS, BCApplyBoxes, regApplyBox, bcsOnly,
                          m_stencilName,  m_domainBCName, m_ebbcName,
                          valbx, valbx,  domainbx, domainbx,
                          m_ptghost,  m_ptghost,    
                          graph,  graph,
                          voludata,  ebfadata,    
                          xfacdata,  yfacdata, zfacdata,    
                          m_dx, false, Point::Zeros());

    a_allIrregPoints.makeEmpty();
    for(int ipt = 0; ipt < a_dstVoFs.size(); ipt++)
    {
      a_allIrregPoints |= a_dstVoFs[ipt].m_pt;
    }

    return 0;
  }

  int petscToChombo(EBLevelBoxData<CELL, 1> & a_deviceDst,
                    const Vec               & a_src)
                    
  {
    using Chombo4::DataIterator;
    PetscErrorCode ierr;
    const PetscScalar *arr;
    ierr = VecGetArrayRead(a_src,&arr);  CHKERRQ(ierr);

    //copy data to host.  The PeTsC team has promised to make this
    //step unnecessary by providing the pointer to the vector that lives
    //on the device
    hostfactreal_t factory(m_graphs);
    LevelData<EBHostData<CELL, Real, 1> > hostDst(m_grids, 1, a_deviceDst.ghostVect(), factory);
    
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
          hostfab(vof, 0) = arr[ki];
        }
      }
    }

    ///copy data back to device.   Soon, the pEtSc team will render this unnecessary.
    EBLevelBoxData<CELL, 1>::copyToDevice(a_deviceDst, hostDst);
    return 0;
  }

  
  virtual PetscInt formMatrix()
  {
    
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
        vector< EBIndex<CELL> >            dstVoFs;
        vector< Proto::LocalStencil<CELL, Real> > wstencil;
        Stencil<Real>    regStencilInterior;
        PointSet irregSet;
        //use strings to get stencil for this particular operator.
        getStencilComponents(dstVoFs, wstencil, regStencilInterior, irregSet, dit[ibox]);
        
        auto& graph = (*m_graphs)[dit[ibox]];
        auto  grid  =     m_grids[dit[ibox]];
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
            ierr = MatSetValues(m_mat,1,&irowpet,1,&icolpet,&stenwgt,INSERT_VALUES);
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
    pout() << "      KSP:: iteration = " << it << " residual norm = " <<  rnorm << std::endl;
    return 0;
  }
  virtual PetscInt setupSolver()
  {
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
    ierr = KSPSetInitialGuessNonzero(m_ksp, m_nz_init_guess ? PETSC_TRUE : PETSC_FALSE );CHKERRQ(ierr);

    //there is some stuff here about blocksize but it appears to be optional.
    return 0;
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
    using Chombo4::DataIterator;
    IntVect gidghost = 4*IntVect::Unit;
    m_gids.define(m_grids, 1, gidghost, hostfactint_t(m_graphs));
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
    CH4_SPMD::gather<int>(m_numPtsAllProc, numPtsThisProc, 0);
    CH4_SPMD::broadcast<Vector<int> >(  m_numPtsAllProc, 0);
    int startgid = 0;
    for(int iproc = 0; iproc < procID(); iproc++)
    {
      startgid += m_numPtsAllProc[iproc];
    }
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
    Chombo4::Copier exchangeCopier;
    exchangeCopier.exchangeDefine(m_grids, gidghost);
    m_gids.exchange(exchangeCopier);
  }
  //index space--allows us to to matrix to solution space
  LevelData<EBHostData<CELL, int, 1> > m_gids;
  PetscInt m_gid0;
  PetscInt m_NN;
  bool m_nz_init_guess;


  shared_ptr<GeometryService<order> >  m_geoserv;
  shared_ptr< dictionary_t>            m_ebdictionary;
  shared_ptr<LevelData<EBGraph> >      m_graphs;
  Chombo4::DisjointBoxLayout                    m_grids;
  Chombo4::Box                                  m_domain;
  Real                                 m_dx;
  
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
  char m_prestring[32];
  Real m_alpha, m_beta;
};

#endif

#endif 
