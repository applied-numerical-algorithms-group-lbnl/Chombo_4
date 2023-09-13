#ifndef __Chombo_EBCM_PETSC_Framework_
#define __Chombo_EBCM_PETSC_Framework_

#ifdef CH_USE_PETSC
#include "petsc.h"
#include "petscmat.h"
#include "petscvec.h"
#include "petscksp.h"
#include "petscviewer.h"
#include "petscmat.h"
/*
  Chombo_EBCM_Graph.H is where MetaDataLevel lives.  It made sense at the time.
  I should use this opportunity to apologize to the world for the many 
  awful choices I have made in this life in relation to file names and organization.  
  I regret all the consternation I have caused.
*/
#include "Chombo_EBCM_Graph.H"
#include "Chombo_EBLevelBoxData.H"

/// 
/**
   An EBCM interface to the magical world of PETSc, 
   wherein the evil bits of distributed linear algebra are done for us.
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

    typedef Chombo4::Box                                   ch_box;
    typedef Chombo4::DisjointBoxLayout                     ch_dbl;
    typedef Chombo4::BoxIterator                           ch_bit;
    typedef Chombo4::IntVect                               ch_iv;
    typedef Chombo4::DataIterator                          ch_dit;
                                                           
    typedef ::Mat                                          PETSc_mat;
    typedef ::Vec                                          PETSc_vec;
    
    typedef  EBCM::MetaDataLevel<             ebcm_order>  ebcm_meta;
    typedef  EBCM::EBCM_Graph<                ebcm_order>  ebcm_graph;
    typedef  EBCM::HostLevelData<int,    1,   ebcm_order>  ebcm_int_sca_data;
    typedef  EBCM::HostLevelData<double, 1,   ebcm_order>  ebcm_dou_sca_data;
    typedef  EBCM::HostLevelData<int,    DIM, ebcm_order>  ebcm_int_vec_data;
    typedef  EBCM::HostLevelData<double, DIM, ebcm_order>  ebcm_dou_vec_data;
    typedef  EBCM::SubVolumeVector<           ebcm_order>  ebcm_subvol_vec;
    typedef  EBCM::EBCM_Volu<                 ebcm_order>  ebcm_volu;
    typedef  EBCM::Algorithm_Framework<       ebcm_order>  algo_framework;
    
    typedef Proto::Point                                   pr_pt;
    typedef Proto::IndexedMoments<DIM  ,      ebcm_order > pr_mom_dim;
    typedef Proto::IndexedMoments<DIM-1,      ebcm_order > pr_mom_dmo;
    typedef Proto::IndexTM<double, DIM>                    pr_itm_r_dim;
    typedef Proto::IndexTM<int , DIM>                      pr_itm_i_dim;
    typedef Proto::IndexTM<double, DIM-1>                  pr_itm_r_dmo;
    typedef Proto::IndexTM<int , DIM-1>                    pr_itm_i_dmo;
    
    typedef EBCM::neighborhood       <  ebcm_order>       ebcm_neighborhood;
    /// Class to replace old vofstencil.  
    /**
       This will have to be generalized for multiple vofs per cell.
     **/
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
        this->resize(0);
      }

      /// 
      /**
         Append the entry to the vector if the point is not already in the vector.
         Otherwise, just increment the weight.
      */
      void add(const pair<pr_pt, double>& a_entry)
      {
        bool alreadyhere = false;
        for (int ivof = 0; ivof < this->size(); ivof++)
        {
          const pr_pt& curpt = (*this)[ivof].first;
          if(curpt == a_entry.first)
          {
            alreadyhere = true;
            (*this)[ivof].second += a_entry.second;
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
          (*this)[ivof].second *= a_scale;
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
      shared_ptr<PETSc_mat> m_mat_ptr;
      shared_ptr<PETSc_vec> m_xx_ptr, m_rr_ptr, m_bb_ptr;
      shared_ptr<ebcm_int_sca_data> m_row_map_ptr; 
      shared_ptr<ebcm_meta> m_meta_ptr;
      int                   m_ghost;
      //each box in data iterator we get total points
      std::vector<int>      m_numPtsPerBox;
      int                   m_numPtsThisProc;
      int                   m_numPtsAllProcs;
      int                   m_startPtThisProc;
      bool                  m_defined;
      virtual ~Base_PETSc_Op()
      {;}
  
      Base_PETSc_Op(const shared_ptr< ebcm_meta  >               & a_meta_ptr,
                    const int                                    & a_ghost,
                    bool                                           a_leaveUndefined = false,  //for testing
                    bool                                           a_print          = false)  //for testing
      {
        CH_TIME("EBPETScSolver::EBPETScSolver");
        m_meta_ptr     =    a_meta_ptr;
        m_ghost        =    a_ghost;
        m_defined = !a_leaveUndefined;
        if(!a_leaveUndefined)
        {
    
          if(a_print)
          {
            Chombo4::pout() << "Base_PETSc_Op: creating map of locations in space to matrix row." << endl;
          }
          createMatrixRowMap(a_print);

          if(a_print)
          {
            Chombo4::pout() << "Base_PETSc_Op:: creating the matrix" << endl;
          }

          defineMatrixAndVectors(a_print);
        }
      }

      ///copy chombo data to PETSc vector
      void chomboToPETSc(shared_ptr<PETSc_vec>         a_phi_PETSc,
                         shared_ptr<ebcm_dou_sca_data> a_phi_chombo,
                         bool a_print= false) const
      {
        const auto& grids   =   m_meta_ptr->m_grids;
        const auto& graphs  = *(m_meta_ptr->m_graphs);
        const auto& mapdata = *(m_row_map_ptr->m_data);
        auto     dit = grids.dataIterator();
        auto& PETSc_vec = *a_phi_PETSc;

        for(int ibox = 0; ibox < dit.size(); ibox++)
        {
          const auto & chombo_fab = (*a_phi_chombo)[dit[ibox]];
          const auto & valid  =   grids[dit[ibox]];
          const auto & graph  =  graphs[dit[ibox]];
          const auto & mapfab = mapdata[dit[ibox]];
          ebcm_subvol_vec volu_vec(graph, valid, a_print);
          for(int ivec = 0; ivec < volu_vec.size(); ivec++)
          {
            PetscInt global_vec_row  =    mapfab( volu_vec[ivec].m_pt, 0);
            double        chombo_val = chombo_fab(volu_vec[ivec].m_pt, 0);
            VecSetValues(PETSc_vec, 1, &global_vec_row, &chombo_val, INSERT_VALUES); 
          }
        }
        /**
           PETSc apparently requires this when the values change.
           The documentation says it assembles the vectors.   
           (something to do with values on correct processors).
           Seems like a good idea.
        **/
        VecAssemblyBegin(PETSc_vec); 
        VecAssemblyEnd(  PETSc_vec); 
      }
      
      ///copy PETSc vector to chombo data
      void PETScToChombo(shared_ptr<ebcm_dou_sca_data> a_phi_chombo,
                         shared_ptr<PETSc_vec>         a_phi_PETSc,
                         bool a_print= false) const
      {

        const auto& grids   =   m_meta_ptr->m_grids;
        const auto& graphs  = *(m_meta_ptr->m_graphs);
        const auto& mapdata = *(m_row_map_ptr->m_data);
        auto     dit = grids.dataIterator();

        
        auto& PETSc_vec = *a_phi_PETSc;
        const PetscScalar *raw_array;
        VecGetArrayRead(PETSc_vec, &raw_array);
        for(int ibox = 0; ibox < dit.size(); ibox++)
        {
          auto       & chombo_fab = (*a_phi_chombo)[dit[ibox]];
          const auto & valid  =   grids[dit[ibox]];
          const auto & graph  =  graphs[dit[ibox]];
          const auto & mapfab = mapdata[dit[ibox]];
          ebcm_subvol_vec volu_vec(graph, valid, a_print);
          for(int ivec = 0; ivec < volu_vec.size(); ivec++)
          {
            PetscInt global_vec_row = mapfab(volu_vec[ivec].m_pt, 0);
            PetscInt localIndex = global_vec_row - m_startPtThisProc;
            double PETSc_val = raw_array[localIndex];
            for(int icell = 0; icell < volu_vec[ivec].m_cells.size(); icell++)
            {
              const auto& cell = volu_vec[ivec].m_cells[icell];
              chombo_fab( cell, 0) = PETSc_val;
            }
          }
        }
      }
      
      /// compute a_lph = L(phi)
      void applyOp(shared_ptr<ebcm_dou_sca_data> a_lph,
                   shared_ptr<ebcm_dou_sca_data> a_phi,
                   bool a_print= false) const
      {

        CH_TIME("BaseEBPETScSolver::applyOp");

        VecZeroEntries(*m_xx_ptr);
        VecZeroEntries(*m_bb_ptr);

        chomboToPETSc(m_xx_ptr, a_phi);

        if(a_print)
        {
          Chombo4::pout() << "Base_PETSc_Op::applyOp before mult: m_xx = " << endl;
          VecView(*m_xx_ptr , PETSC_VIEWER_STDOUT_WORLD);
          /**
             This one produces a very large amount of output.
             Turn it on if you need to but you should probably
             route stdout to  a log file. 
          **/
          /**
          Chombo4::pout() << "Base_PETSc_op::applyOp before mult : m_mat = " << endl;
          MatView(*m_mat_ptr, PETSC_VIEWER_STDOUT_WORLD);
          **/
        }

        MatMult(*m_mat_ptr, *m_xx_ptr, *m_bb_ptr); 

        if(a_print)
        {
          Chombo4::pout() << "Base_PETSc_Op::applyOp after mult: m_bb = " << endl;
          VecView(*m_bb_ptr, PETSC_VIEWER_STDOUT_WORLD);
        }
    
        PETScToChombo(a_lph, m_bb_ptr);


      }
      ///this gets done by the derived class
      virtual shared_ptr<local_stencil_t> getLocalStencil(const pr_pt& a_pt, const ebcm_graph& a_graph, bool a_print = false) = 0;


    private:
      /// banning weak construction but assignment and copy constructors should work fine.
      Base_PETSc_Op();
      
    protected:
      ///makes m_row_map_ptr (operator uses this to map between volumes and PETSc matrix rows)
      PetscInt createMatrixRowMap(bool a_print = false)
      {
        m_row_map_ptr = shared_ptr<ebcm_int_sca_data>(new ebcm_int_sca_data(m_meta_ptr, m_ghost));
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

        Chombo4::pout() << "Base_PETSc_Op::createMatrixRowMap: numptsThisProc = " << numPtsThisProc << endl;

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
        Chombo4::pout()  << "Base_PETSc_Op::getMap: procID = " << procID() << ", numPtsThisProc = " << numPtsThisProc << ", startgid = " << startgid <<", totalNumPts = " << totalNumPts <<  endl;
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

        Chombo4::pout() << "Base_PETSc_Op::getMap calling exchange" << endl;

        m_row_map_ptr->exchange(a_print);
        
        return 0;
      }//end function createMatrixRowMap

      ///makes m_mat_ptr and some compatible vectors.   The matrix gets its valued filled.
      PetscInt defineMatrixAndVectors(bool a_print = false)
      {
        CH_TIME("Base_PETSc_Op::defineMatrixAndVectors");
        // create matrix and vectors
        m_mat_ptr = shared_ptr<PETSc_mat>(new PETSc_mat());
        m_bb_ptr  = shared_ptr<PETSc_vec>(new PETSc_vec());
        m_xx_ptr  = shared_ptr<PETSc_vec>(new PETSc_vec());
        m_rr_ptr  = shared_ptr<PETSc_vec>(new PETSc_vec());
        PetscInt nnz_per_row = 1000;
#ifdef CH_MPI
        MPI_Comm wcomm = Chombo_MPI::comm;
#else
        MPI_Comm wcomm = PETSC_COMM_SELF;
#endif
        PetscInt ierr;
        ierr = MatCreate(wcomm, & (*m_mat_ptr));CHKERRQ(ierr);
        ierr = MatSetSizes(*m_mat_ptr,m_numPtsThisProc,m_numPtsThisProc,m_numPtsAllProcs,m_numPtsAllProcs);CHKERRQ(ierr);
        ierr = MatSetFromOptions(*m_mat_ptr);

#ifdef CH_MPI
        ierr = MatMPIAIJSetPreallocation(*m_mat_ptr, nnz_per_row, NULL, nnz_per_row, NULL);CHKERRQ(ierr);
#else
        ierr = MatSeqAIJSetPreallocation(*m_mat_ptr, nnz_per_row, NULL); CHKERRQ(ierr);
#endif    

        VecCreate(wcomm, &(*m_bb_ptr));
        VecSetSizes(*m_bb_ptr, PETSC_DECIDE, m_numPtsAllProcs);
        VecSetFromOptions(  *m_bb_ptr);
        ierr = VecDuplicate(*m_bb_ptr, &(*m_xx_ptr)); CHKERRQ(ierr);
        ierr = VecDuplicate(*m_bb_ptr, &(*m_rr_ptr)); CHKERRQ(ierr);

        //sets the actual values in *m_mat_ptr
        ierr = formMatrix(a_print); CHKERRQ(ierr);
        return 0;
      }

      ///sets the actual values in *m_mat_ptr
      PetscInt formMatrix(bool a_print = false)
      {
    
        CH_TIME("Base_PETSc_Op::formMatrix");
        
        const auto& graphs = *(m_meta_ptr->m_graphs);
        const auto& grids  =  (m_meta_ptr->m_grids);
        auto dit = grids.dataIterator();

        //get all the info on points on this proc.
        int numPtsThisProc = 0;
        vector<int> numPtsPerBox(dit.size(), 0);
        auto & stencil_matrix = *(this->m_mat_ptr);
        for(int ibox =0; ibox < dit.size(); ibox++)
        {
          const auto& valid =  grids[dit[ibox]];
          const auto& graph = graphs[dit[ibox]];
          ebcm_subvol_vec volvec(graph, valid, a_print);
          const auto& mapfab = (*m_row_map_ptr)[dit[ibox]];
          for(int ivol = 0; ivol < volvec.size(); ivol++)
          {
            const auto& volu = volvec[ivol];
            shared_ptr<local_stencil_t>  local = getLocalStencil(volu.m_pt, graph);
            PetscInt                     irow  =          mapfab(volu.m_pt, 0);

            for(int isten = 0; isten < local->size(); isten++)
            {
              const auto& stenpair = (*local)[isten];
              PetscInt icol = mapfab(stenpair.first, 0);
              double weight =        stenpair.second;

              /// This is PETSc for
              /// (*m_mat)(irow, icol) = weight;
              MatSetValues(stencil_matrix,1,&irow,1,&icol,&weight,INSERT_VALUES);
            }
          }
        }
        
        /**
           PETSc apparently requires this step.
           The documentation just says it assembles the matrix, making sure 
           it has the right values on the right processors.
         **/
        MatAssemblyBegin(stencil_matrix,MAT_FINAL_ASSEMBLY);
        MatAssemblyEnd(  stencil_matrix,MAT_FINAL_ASSEMBLY);
        
        return 0;
      }

      }; //end class Base_PETSc_Op

    /// 
    /**
       Class to force compilation of all this lovely code. It is useful but unusable.  
       If a class can be a zen koan, this may qualify.
    **/
    class DummyDerivedOp: public Base_PETSc_Op
    {
    public:
      DummyDerivedOp(const shared_ptr< ebcm_meta  >               & a_meta_ptr,
                     const int                                    & a_ghost,
                     bool                                           a_print = false)
        :Base_PETSc_Op(a_meta_ptr, a_ghost, true, a_print) //the true is to leave the basePETScop undefinde
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
    

    ///Returns a local stencil for \int (\grad \phi) \cdot (\ehat^facedir)
    /**
       Get stencil for integrated gradient in a  coordinate direction. 
       This is the same as the average gradient * face area.
       \int (\grad \phi) \cdot (\ehat^d) dA = A_f < \frac{\partial \phi}{\partial x_d} >
    */
    static local_stencil_t
    getStencilFromMat(const eigen_mat                                & a_matrix,
                      shared_ptr< EBCM::neighborhood< ebcm_order > > & a_locality)
    {
      int nrows_volumes   = a_locality->m_weig_volumes.size();
      local_stencil_t retval;
      const auto& volvec = *a_locality->m_volumes;
      for( int irow = 0; irow < nrows_volumes; irow++)
      {
        pair<pr_pt, double> entry;
        const ebcm_volu& volume = volvec[irow];
        entry.first = volume.m_pt;
        entry.second = a_matrix(irow, 0);
        retval.add(entry);
      }
      return retval;
    }


    /**
       quick test to get this stuff to compile
    **/ 
    static  void dummy_petsc_test()
    {

#ifdef CH_USE_PETSC
      shared_ptr< ebcm_meta  >  meta_ptr;
      int  nghost = 4586;    
      bool print  = false;
      DummyDerivedOp(meta_ptr, nghost, print);
#endif      
    }

  }; //end PETSc_framwork class
}
#endif

#endif 
