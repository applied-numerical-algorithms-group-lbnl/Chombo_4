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
   Chombo_EBCM_Graph.H is where MetaDataLevel lives.  It made sense at the time.
   I should use this opportunity to apologize to the world for the many 
   awful choices I have made in this life in relation to file names and organization.  
   I regret all the consternation I have caused.
**/
#include "Chombo_EBCM_Graph.H"
#include "Chombo_EBLevelBoxData.H"

/// An EBCM interface to the magical world of PETSc,  where the evil bits of distributed linear algebra are done for us.
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
    typedef EBCM::EBCM_Graph<ebcm_order>            ebcm_graph;
    typedef EBCM::HostLevelData<int, 1, ebcm_order> ebcm_data;
    typedef EBCM::SubVolumeVector<      ebcm_order> ebcm_subvol_vec;
    typedef EBCM::Algorithm_Framework<  ebcm_order> ebcm_framework;

    ///
    struct local_stencil_t: public std::vector< pair< pr_pt, double > >
    {
    public:
 
      ///
      /**
         adds each member of input
       */
      void
      operator+= (const local_stencil_t& a_input)
      {
        for(int ivec = 0; ivec < a_input.size(); ivec++)
        {
          add(a_input[ivec]);
        }
      }
      
      ///
      void clear()
      {
        m_stencil.resize(0);
      }

      /// 
      /**
         Add to the vector if the point is not already in the vector.
         Otherwise, just increment the weight.
       */
      void add(const pair<pr_pt, double>& a_entry)
      {
        bool alreadyhere = false;
        for (int ivof = 0; ivof < this->size(); ivof++)
        {
          const pr_pt& curpt (*this)[ivec].first;
          if(curpt == a_entry.first)
          {
            alreadyhere = true;
            curpt.second += a_entry.second;
          }
        }
        if (!alreadyhere)
        {
          this->push_back(a_entry);
        }
      }

      ///scale the stencil by a constant
      void operator*= (const double& a_scale)
      {
        for (int ivof = 0; ivof < this->size(); ivof++)
        {
          (*this)[ivec].second *= a_scale;
        }
      }
        
    }; //end struct local_stencil_t
    
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
        createMatrixRowMap(a_print);

        if(a_print)
        {
          Chombo4::pout() << "Base_Petsc_Op:: creating the matrix" << endl;
        }

        createOperatorPetscMatrix(a_print);
      }

      ///this gets done by the derived class
      virtual shared_ptr<local_stencil_t> getLocalStencil(const pr_pt& a_pt, const ebcm_graph& a_graph, bool a_print = false) = 0;


    private:
      /// banning weak construction but assignment and copy constructors should work fine.
      Base_PETSc_Op();
      
    protected:
      ///makes m_row_map_ptr (operator uses this to map between volumes and petsc matrix rows)
      PetscInt createMatrixRowMap(bool a_print = false)
      {
        m_row_map_ptr = shared_ptr<ebcm_data>(new ebcm_data(m_meta_ptr, m_ghost));
        const auto& graphs = *(m_meta_ptr->m_graphs);
        const auto& grids  =  (m_meta_ptr->m_grids);
        auto dit = grids.dataIterator();

        //get all the info on points on this proc.
        int numPtsThisProc = 0;
        vector<int> numPtsPerBox(dit.size(), 0);
        for(int ibox =0; ibox < dit.size(); ibox++)
        {
          const auto& valid =  grids[dit[ibox]];
          const auto& graph = graphs[dit[ibox]];
          //valid is correct here.  ghost is handled via exchange.

          ebcm_subvol_vec subvollocal(graph, valid, a_print);
          numPtsThisProc     += subvollocal.size();
          numPtsPerBox[ibox]  = subvollocal.size();
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

        //get all the info on points on this proc into the map
        for(int ibox =0; ibox < dit.size(); ibox++)
        {
          const auto& valid =  grids[dit[ibox]];
          const auto& graph = graphs[dit[ibox]];
          ebcm_subvol_vec volvec(graph, valid, a_print);
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
        
        return 0;
      }//end function createMatrixRowMap

      ///makes m_mat_ptr (operator uses this to map between volumes and petsc matrix rows)
      PetscInt createOperatorPetscMatrix(bool a_print = false)
      {
        CH_TIME("Base_Petsc_Op::createOperatorPetscMatrix");
        // create matrix
        PetscInt nnzrow = 0;
#ifdef CH_MPI
        MPI_Comm wcomm = Chombo_MPI::comm;
#else
        MPI_Comm wcomm = PETSC_COMM_SELF;
#endif
        PetscInt ierr;
        ierr = MatCreate(wcomm, & (*m_mat_ptr));CHKERRQ(ierr);
        ierr = MatSetOptionsPrefix(*m_mat_ptr,"");CHKERRQ(ierr);

        ierr = MatSetSizes(*m_mat_ptr,m_numPtsThisProc,m_numPtsThisProc,m_numPtsAllProcs,m_numPtsAllProcs);CHKERRQ(ierr);
        ierr = MatSetType( *m_mat_ptr,MATDENSE);CHKERRQ(ierr);
#ifdef CH_MPI
        ierr = MatMPIDenseSetPreallocation(*m_mat_ptr, NULL);CHKERRQ(ierr);
#else
        PetscInt matsize = m_numPtsThisProc*m_numPtsThisProc;
        PetscScalar* matspace = new PetscScalar[matsize];
        ierr = MatSeqDenseSetPreallocation(*m_mat_ptr, matspace);CHKERRQ(ierr);
#endif    

        ierr = MatSetFromOptions( *m_mat_ptr ); CHKERRQ(ierr);

        //sets the actual values in *m_mat_ptr
        ierr = formMatrix(a_print); CHKERRQ(ierr);
        return 0;
      }

      ///sets the actual values in *m_mat_ptr
      PetscInt formMatrix(bool a_print = false)
      {
    
        CH_TIME("Base_Petsc_Op::formMatrix");
        
        const auto& graphs = *(m_meta_ptr->m_graphs);
        const auto& grids  =  (m_meta_ptr->m_grids);
        auto dit = grids.dataIterator();

        //get all the info on points on this proc.
        int numPtsThisProc = 0;
        vector<int> numPtsPerBox(dit.size(), 0);
        for(int ibox =0; ibox < dit.size(); ibox++)
        {
          const auto& valid =  grids[dit[ibox]];
          const auto& graph = graphs[dit[ibox]];
          ebcm_subvol_vec volvec(graph, valid, a_print);
          for(int ivol = 0; ivol < volvec.size(); ivol++)
          {
            const auto& volu = volvec[ivol];
            shared_ptr<local_stencil_t>  local = getLocalStencil(volu.m_pt, graph);
            PetscInt irow = (*m_row_map_ptr)(local->m_pt, 0);
            for(int isten = 0; isten < sten.size(); isten++)
            {
              const auto& stenpair = local->m_stencil[isten];
              PetscInt icol = (*m_row_map_ptr)(stenpair.first, 0);
              double weight = stenpair.second;
              /// This is PETSc for
              /// (*m_mat)(irow, icol) = weight;
              MatSetValues(this->m_mat,1,&irow,1,&icol,&weight,INSERT_VALUES);
            }
          }
        }
        return 0;
      }

    }; //end class Base_Petsc_Op

    ///
    /**
      Class to force compilation of all this lovely code.
      It is useful but unusable. 
      If a class could be a zen koan, this might be it.
    **/
    class DummyDerivedOp: public Base_PETSc_Op
    {
    public:
      DummyDerivedOp(const shared_ptr< ebcm_meta  >               & a_meta_ptr,
                     const int                                    & a_ghost,
                     bool                                           a_print = false)
        :Base_PETSc_Op(a_meta_ptr, a_ghost, a_print)
      {
      }

      virtual ~DummyDerivedOp()
      {
      }
      
      virtual shared_ptr<local_stencil_t> getLocalStencil(const pr_pt& a_pt, const ebcm_graph& a_graph, bool a_print = false)
      {
        shared_ptr<local_stencil_t> retval(new local_stencil_t());
        return retval;
      }
    };
    

    ///
    template<CENTERING cent>
    local_stencil_t
    getHelmholtzFlux()
    {
    }
    
    ///
    /**
      Class which forms a helmholtz operator.
      This creates an operator of the form
      L(phi) = alpha phi + beta * div(F(phi))
    **/
    class HelmholtzOp: public Base_PETSc_Op
    {
    public:

      typedef Proto::Side         pr_sd;
      typedef Proto::Box          pr_box;
      
      double                                 m_alpha;
      double                                 m_beta;
      int                                    m_weight_power;
      
      HelmholtzOp(const double                                 & a_alpha,
                  const double                                 & a_beta,
                  const int                                    & a_weight_power,
                  const shared_ptr< ebcm_meta  >               & a_meta_ptr,
                  const int                                    & a_ghost,
                  bool                                           a_print = false)
        :Base_PETSc_Op(a_meta_ptr, a_ghost, a_print)
      {
        m_alpha        = a_alpha;       
        m_beta         = a_beta;        
        m_weight_power = a_weight_power;        
      }

      virtual ~HelmholtzOp()
      {
      }
      
      virtual shared_ptr<local_stencil_t> getLocalStencil(const pr_pt& a_pt, const ebcm_graph& a_graph, bool a_print = false)
      {
        shared_ptr<local_stencil_t> retval(new local_stencil_t());
        pr_bx pointbx(a_pt, a_pt);
        pr_bx grbx = pointbox.grow(pr_pt::Ones(m_nghost));
        
        ebcm_volu volume = a_graph.getVolumeCoveringPoint(a_pt);
        for(
        
        return retval;
      }

    };
    
    ///test function for example/EBCM/hoeb_truncation
    static  void run_hoeb_truncation_tests()
    {
      
      int nghost = 4;
      shared_ptr< ebcm_meta  > ebcm;
      ch_dbl                   grids;
      double                   dx;
      ch_probdom               domain;
      ebcm_framework::makeMergedGeometry(ebcm, grids, dx, domain, true, true);
      
      ///testing compilation of the base framework
      typedef DummyDerivedOp   ebcm_zen_koan;
      ebcm_zen_koan  number4586(ebcm, nghost, false);

      ///Helmholtz operator declaration
      double alpha     = 4586;
      double beta      = 4586;
      int weight_power = 4586;
      ParmParse pp("trunc_test");
      pp.get("alpha", alpha);
      pp.get("beta" , beta);
      pp.get("weight_power", weight_power);
      
      typedef HelmholtzOp   ebcm_helm_op;
      ebcm_helm_op   helmholtzOperator(alpha, beta, weight_power, ebcm, nghost, false);
      
    }
  }; //end petsc_framwork class
}
#endif

#endif 
