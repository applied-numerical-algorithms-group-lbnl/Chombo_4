#ifndef __Chombo_EBCM_PETSC_Framework_
#define __Chombo_EBCM_PETSC_Framework_

#ifdef CH_USE_PETSC
#include "petsc.h"
#include "petscmat.h"
#include "petscvec.h"
#include "petscksp.h"
#include "petscviewer.h"
#include "petscmat.h"
/**
   This is where MetaDataLevel lives.  It made sense at the time.
   I should use this opportunity to apologize to the world for the many 
   awful choices I have made in this life in relation to file names and organization.  
   I regret all the consternation I have caused.
**/
#include "Chombo_EBCM_Graph.H"
#include "Chombo_EBLevelBoxData.H"

/// An EBCM interface to the magical world of PETSc,  where the evil bits of linear algebra are done for us.
/**
   I plan to use this for everything.  Geometric multigrid requires
   more infrastructure than I am willing to write so I will use PETSc
   for solves and truncation error tests and I will be hooking in
   SLEPc to do the eigenvalue analysis. --dtg
**/
namespace EBCM
{

  template <int ebcm_order>
  class PETSc_Framework
  {
  public:
    typedef Chombo4::Box                            ch_box;
    typedef Chombo4::DisjointBoxLayout              ch_dbl;
    typedef Chombo4::BoxIterator                    ch_bit;
    typedef Chombo4::IntVect                        ch_iv;
    typedef Chombo4::DataIterator                   ch_dit;
    
    typedef ::Mat                                   petsc_mat;
    typedef Proto::Point                            pr_pt;
    typedef EBCM::MetaDataLevel<ebcm_order>         ebcm_meta;
    typedef EBCM::HostLevelData<int, 1, ebcm_order> ebcm_data;
    typedef EBCM::SubVolumeVector<      ebcm_order> ebcm_subvol_vec;

    ///
    /**
       This is a wrapper around PETSc's  ::Mat object.  The specfics of the stencil 
       are computed in the derived class.   All data is public because
       this design is cleaner than hiding behind access functions.
    **/
    class Base_PETSc_Op
    {
    public:
      shared_ptr<petsc_mat> m_mat_ptr;
      shared_ptr<ebcm_data> m_row_map_ptr; 
      shared_ptr<ebcm_meta> m_meta_ptr;
      int                   m_ghost;
      //each box in data iterator we get total points
      std::vector<int>      m_numPtsPerBox;
      int                   m_numPtsThisProc;
      int                   m_numPtsAllProcs;
      int                   m_startPtThisProc;
      
      virtual ~Base_PETSc_Op()
      {;}
  
      Base_PETSc_Op(const shared_ptr< ebcm_meta  >               & a_meta_ptr,
                    const int                                    & a_ghost,
                    bool                                           a_print = false)
      {
        CH_TIME("EBPetscSolver::EBPetscSolver");
        m_meta_ptr     =    a_meta_ptr;
        m_ghost        =    a_ghost;
    
        if(a_print)
        {
          Chombo4::pout() << "Base_Petsc_Op: creating map of locations in space to matrix row." << endl;
        }
        ///this would be convenient to keep around but it must be temporary as it contains all kinds of moment data that
        ///will no longer be necessary  once the matrix is defined.
        vector<shared_ptr<ebcm_subvol_vec> > a_all_valid_vols;
        createMatrixRowMap(a_all_valid_vols, a_print);

        if(a_print)
        {
          Chombo4::pout() << "Base_Petsc_Op:: creating the matrix" << endl;
        }

        createOperatorPetscMatrix(a_all_valid_vols, a_print);
      }

    private:
      /// banning weak construction but assignment and copy constructors should work fine.
      Base_PETSc_Op();
      
    protected:
      ///makes m_row_map_ptr (operator uses this to map between volumes and petsc matrix rows)
      void createMatrixRowMap(vector<shared_ptr<ebcm_subvol_vec> > a_all_valid_vols, 
                              bool a_print = false)
      {
        m_row_map_ptr = shared_ptr<ebcm_data>(new ebcm_data(m_meta_ptr, m_ghost));
        const auto& graphs = *(m_meta_ptr->m_graphs);
        const auto& grids  =  (m_meta_ptr->m_grids);
        auto dit = grids.dataIterator();

        a_all_valid_vols.resize(dit.size());
        
        //get all the info on points on this proc.
        int numPtsThisProc = 0;
        vector<int> numPtsPerBox(dit.size(), 0)
        for(int ibox =0; ibox < dit.size(); ibox++)
        {
          const auto& valid =  grids[dit[ibox]];
          const auto& graph = graphs[dit[ibox]];
          //valid is correct here.  ghost is handled via exchange.
          a_all_valid_vols[ibox]= shared_ptr<ebcm_subvol_vec>(new ebcm_subvol_vec(graph, valid, a_print));
          numPtsThisProc     += a_all_valid_vols[ibox]->size();
          numPtsPerBox[ibox]  = a_all_valid_vols[ibox]->size();
        }

        Chombo4::pout() << "Base_Petsc_Op::createMatrixRowMap: numptsThisProc = " << numPtsThisProc << endl;

#ifdef CH_MPI
        std::vector<int> numPtsAllProcs(CH4_SPMD::numProc());
        MPI_Gather(&numPtsThisProc, 1, MPI_INT, &numPtsAllProcs[0], 1, MPI_INT, 0, Chombo_MPI::comm);
    
        MPI_Bcast(numPtsAllProcs.data(), numPtsAllProcs.size(), MPI_INT, 0, Chombo_MPI::comm);

        int totalNumPts = 0;
        for(int iproc = 0; iproc < numPtsAllProcs.size(); iproc++)
        {
          totalNumPts += numPtsAllProcs[iproc];
        }
        //decide which location maps to the first one of this proc
        int startgid = 0;
        for(int iproc = 0; iproc <  procID(); iproc++)
        {
          startgid += numPtsAllProcs[iproc];
        }
        Chombo4::pout()  << "Base_Petsc_Op::getMap: procID = " << procID() << ", numPtsThisProc = " << numPtsThisProc << ", startgid = " << startgid <<", totalNumPts = " << totalNumPts <<  endl;
#else
        int startgid = 0;
        PetscInt totalNumPts = numPtsThisProc;
#endif
        m_startPtThisProc      = startgid;
        m_numPtsThisProc       = numPtsThisProc;
        m_numPtsAllProcs       = totalNumPts;
        m_numPtsPerBox         = numPtsPerBox;
        auto& mapdata = *(m_row_map_ptr->m_data);
        int curgid = startgid;
        for(int ibox = 0; ibox < dit.size(); ibox++)
        {
          const auto& volvec = *(a_all_valid_vols[ibox]);
          auto      & mapfab =      mapdata[dit[ibox]];
          for(int ivec = 0; ivec < volvec.size(); ivec++)
          {
            const auto& volume = volvec[ivec];
            mapfab(volume.m_pt, 0) = curgid;
            for(int icell = 0; icell < volume.m_cells.size(); icell++)
            {
              const auto& cell = volume.m_cells[icell];
              mapfab(cell, 0) = curgid;
            }
            curgid++;
          }
        }

        Chombo4::pout() << "Base_Petsc_Op::getMap calling exchange" << endl;

        m_row_map_ptr->exchange(a_print);
      }//end function createMatrixRowMap

      ///makes m_mat_ptr (operator uses this to map between volumes and petsc matrix rows)
      void createOperatorPetscMatrix(vector<shared_ptr<ebcm_subvol_vec> > &  a_all_valid_vols,
                                     bool a_print = false)
      {
        CH_TIME("Base_Petsc_Op::createOperatorPetscMatrix");
        // create matrix
        PetscInt nnzrow = 0;
#ifdef CH_MPI
        MPI_Comm wcomm = Chombo_MPI::comm;
#else
        MPI_Comm wcomm = PETSC_COMM_SELF;
#endif
//    int nc = 1;
        PetscInt ierr;
        ierr = MatCreate(wcomm,&m_mat);CHKERRQ(ierr);
        ierr = MatSetOptionsPrefix(m_mat,"");CHKERRQ(ierr);

        //Chombo4::Chombo4::pout() << "debug:: matsetsize m_numPtsThisProc = " << m_numPtsThisProc << std::endl;
        ierr = MatSetSizes(m_mat,m_numPtsThisProc,m_numPtsThisProc,m_numPtsAllProcs,m_numPtsAllProcs);CHKERRQ(ierr);
//    ierr = MatSetBlockSize(m_mat,nc);CHKERRQ(ierr);
        ierr = MatSetType(m_mat,MATDENSE);CHKERRQ(ierr);
#ifdef CH_MPI
        ierr = MatMPIDenseSetPreallocation(m_mat, NULL);CHKERRQ(ierr);
#else
        PetscInt matsize = m_numPtsThisProc*m_numPtsThisProc;
        PetscScalar* matspace = new PetscScalar[matsize];
        ierr = MatSeqDenseSetPreallocation(m_mat, matspace);CHKERRQ(ierr);
#endif    

        ierr = MatSetFromOptions( m_mat ); CHKERRQ(ierr);


        formMatrix(a_all_valid_vols, a_print);
        return 0;
      }

      struct local_stencil_t
      {
        pr_pt                         m_start
        vector< pair<pr_pt, double> > m_stencil;
      }
      virtual void formMatrix(const vector<shared_ptr<ebcm_subvol_vec> > &  a_all_valid_vols,
                              bool a_print = false)
      {
    
        CH_TIME("Base_Petsc_Op::formMatrix");
        
        int currentRow = m_startPtThisProc;
        ch_dit dit = this->m_grids.dataIterator();
        for(int ibox = 0; ibox < dit.size(); ibox++)
        {
          const auto& volvec = *(a_all_valid_vols[ibox]);
          for(int ivol = 0; ivol < volvec.size(); ivol++)
          {
            shared_ptr<eigen_mat> stenmat =  getStencilMatrix(a_all_valid_vols, currentRow, ibox, ivol, a_print);
            currentRow++;
          }
        }
                int numPtsThisBox = m_numPtsPerBox[ibox];
          for(int ipt = 0; ipt < numPtsThisBox; ipt++)
          {
        vector< EBIndex<ELL> >            dstVoFs;
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
          } 
        }
      }

      ierr = MatAssemblyBegin(this->m_mat,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
      ierr = MatAssemblyEnd(  this->m_mat,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

      
    }
    return 0;
  }
};


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

    }; //end class Base_Petsc_Op
    
    ///test function for example/EBCM/hoeb_petsc_test
    static  void run_hoeb_petsc_tests(shared_ptr< ebcm_meta > & a_meta)
    {
      int nghost = 4;
      typedef Base_PETSc_Op   ebcm_petsc_op;
      ebcm_petsc_op   testOperator(a_meta, nghost, false);

    }
  }; //end petsc_framwork class
}
#endif

#endif 
