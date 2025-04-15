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
#include "Chombo_EBCM_ExactSolutions.H"
#include "Chombo_EBCM_Debug.H"
#include "Chombo_EBCM_HostLevelData.H"
#include "Chombo_EBCM_FluxLevelData.H"

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

                                                            
    typedef ::Mat                                           petsc_mat;
    typedef ::Vec                                           petsc_vec;
    typedef ::PetscInt                                      petsc_int;
    typedef ::PetscScalar                                   petsc_sca;
    typedef Proto::Point                                    pr_pt;
    typedef Proto::HostBoxData<double, 1>                   pr_fab_t;
    typedef Proto::Box                                      pr_bx;
    typedef Chombo4::Box                                    ch_box;
    typedef Chombo4::DisjointBoxLayout                      ch_dbl;
    typedef Chombo4::BoxIterator                            ch_bit;
    typedef Chombo4::IntVect                                ch_iv;
    typedef Chombo4::DataIterator                           ch_dit;

    typedef Chombo4::GeometryService<         ebcm_order >  ch_geoserv;
    typedef  EBCM::Algorithm_Framework<       ebcm_order >  algo_framework;
    typedef  EBCM::neighborhood<              ebcm_order >  ebcm_neighborhood;
    typedef  EBCM::neig_dom_face<             ebcm_order >  ebcm_dom_face;
    typedef  EBCM::neig_cut_face<             ebcm_order >  ebcm_cut_face;
    typedef  EBCM::MetaDataLevel<             ebcm_order >  ebcm_meta;
    typedef  EBCM::Algorithm_Framework<       ebcm_order >  ebcm_framework;
    typedef  EBCM::EBCM_Graph<                ebcm_order >  ebcm_graph;
    typedef  EBCM::HostLevelData<int,    1,   ebcm_order >  ebcm_int_sca_data;
    typedef  EBCM::HostLevelData<double, 1,   ebcm_order >  ebcm_dou_sca_data;
    typedef  EBCM::HostLevelData<int,    DIM, ebcm_order >  ebcm_int_vec_data;
    typedef  EBCM::HostLevelData<double, DIM, ebcm_order >  ebcm_dou_vec_data;
    typedef  EBCM::FluxLevelData<int,    1,   ebcm_order >  ebcm_int_sca_flux;
    typedef  EBCM::FluxLevelData<double, 1,   ebcm_order >  ebcm_dou_sca_flux;
    typedef  EBCM::FluxLevelData<int,    DIM, ebcm_order >  ebcm_int_vec_flux;
    typedef  EBCM::FluxLevelData<double, DIM, ebcm_order >  ebcm_dou_vec_flux;

    typedef  EBCM::SubVolumeVector<           ebcm_order >  ebcm_subvol_vec;
    typedef  EBCM::EBCM_Volu<                 ebcm_order >  ebcm_volu;

    
    typedef Proto::IndexedMoments<DIM  ,      ebcm_order >  pr_mom_dim;
    typedef Proto::IndexedMoments<DIM-1,      ebcm_order >  pr_mom_dmo;
    typedef Proto::IndexTM<double,DIM                    >  pr_itm_r_dim;
    typedef Proto::IndexTM<int   ,DIM                    >  pr_itm_i_dim;
    typedef Proto::IndexTM<double,DIM-1                  >  pr_itm_r_dmo;
    typedef Proto::IndexTM<int   ,DIM-1                  >  pr_itm_i_dmo;
    
    typedef EBCM_Exact::BaseExactSolution<    ebcm_order >  ebcm_exac;



    static int dummy_func()
    {
      return 0;
    }

    /// think vofstencil
    class local_stencil_t
    {
    public:
      local_stencil_t()
      {
        clear();
      }
      std::vector< std::pair<pr_pt, double> > m_entries;
      double m_inhomogeneous_bit;
      ///
      const std::pair<pr_pt, double> & operator[](int a_ivec) const
      {
        return m_entries[a_ivec];
      }

      ///
      std::pair<pr_pt, double> &      operator[](int a_ivec) 
      {
        return m_entries[a_ivec];
      }
      
      double apply(const pr_fab_t&   a_phi_exac) const
                  
      {
        double retval = 0;
        for(int ivec = 0; ivec < m_entries.size(); ivec++)
        {
          const pr_pt &  loc = m_entries[ivec].first;
          const double&  wgt = m_entries[ivec].second;
          const double&  phi = a_phi_exac(loc, 0);
          retval += wgt*phi;
        }
        retval += m_inhomogeneous_bit;
        return retval;
      }
      
      ///
      void
      operator+= (const local_stencil_t& a_input)
      {
        this->m_inhomogeneous_bit += a_input.m_inhomogeneous_bit;
        for(int ivec = 0; ivec < a_input.m_entries.size(); ivec++)
        {
          add(a_input.m_entries[ivec]);
        }
      }
      
      ///
      void clear()
      {
        m_entries.resize(0);
        m_inhomogeneous_bit = 0;//needs to be zero
      }

      /// 
      /**
         Append the entry to the vector if the point is not already in the vector.
         Otherwise, just increment the weight.
      */
      void add(const std::pair<pr_pt, double>& a_entry)
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
          this->m_entries.push_back(a_entry);
        }
      }

      ///
      int size() const
      {
        return m_entries.size();
      }
      
      ///scale the stencil by a constant
      void operator*= (const double& a_scale)
      {
        this->m_inhomogeneous_bit *= a_scale;
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
    template <class linear_operator_t>
    class PETSc_Interface
    {
    public:
      ///
      /**
         This is m_inhomog_bit, the amount the operator will be incremented with 
         due to inhomogeneous boundary condtions.
         What goes here is up to the operator.   
         Go ahead and setVal(0.) because you think you will never need inhomogeneous BCs.   
         If you do that, to test the software, you need lots of functions that exactly 
         conform to your homogeneous boundary conditions.    Maxima hell awaits.   
         Rather than revisit my condo in that circle of hell,  
         I am going to just allow inhomogenous BCs for my operators this time around the 
         development spiral.
         dtg 
         10-9-2024
      **/
      shared_ptr<ebcm_dou_sca_data> m_inhomog_bit;
      shared_ptr<linear_operator_t> m_oper_ptr;
      shared_ptr<petsc_mat> m_mat_ptr;
      shared_ptr<ebcm_meta> m_meta_ptr;
      shared_ptr<petsc_vec> m_lph_ptr;
      shared_ptr<petsc_vec> m_phi_ptr;
      shared_ptr<ebcm_int_sca_data> m_row_map_ptr; 
      int                   m_ghost;
      //each box in data iterator we get total points
      std::vector<int>      m_numPtsPerBox;
      int                   m_numPtsThisProc;
      int                   m_numPtsAllProcs;
      int                   m_startPtThisProc;
      
      PETSc_Interface(const shared_ptr< linear_operator_t  >       & a_oper_ptr,
                      const shared_ptr< ebcm_meta  >               & a_meta_ptr,
                      const int                                    & a_ghost)
      {
        CH_TIME("EBPetscSolver::EBPetscSolver");
        m_oper_ptr     =    a_oper_ptr;
        m_meta_ptr     =    a_meta_ptr;
        m_ghost        =    a_ghost;
    
        //Chombo4::pout() << "Base_Petsc_Op: creating map of locations in space to matrix row." << endl;
        createMatrixRowMap();
        
        //Chombo4::pout() << "Base_Petsc_Op:: creating the matrix" << endl;

        createOperatorPetscMatrix();
      }


      ///
      void  addInhomogeneousBit(shared_ptr<ebcm_dou_sca_data> a_lph_ptr)
      {
        CH_TIME("PETSc_Framework::addInhomogeneousBit");
        (*a_lph_ptr) += (*m_inhomog_bit);
      }

      
      ///YADF
      void writeDetails(shared_ptr<ebcm_dou_sca_data> a_phi_ptr,
                        string   a_identifier) const
      {
        bool killing_me_smalls = true; //this is a debugging routine
        ch_iv lo_coar(D_DECL( 1, 14, 4586));
        ch_iv hi_coar(D_DECL( 1, 15, 4586));
        ch_iv lo_fine(D_DECL( 1, 30, 4586));
        ch_iv hi_fine(D_DECL( 1, 31, 4586));
        ch_iv lo_disp = lo_fine;
        ch_iv hi_disp = hi_fine;
        if(a_identifier == string("coar"))
        {
          lo_disp = lo_coar;
          hi_disp = hi_coar;
        }
        ch_box disp_box(lo_disp, hi_disp);
        int flux_dir = 0;
        Side::LoHiSide flux_side = Side::Hi;
        const auto& graphs = *(m_meta_ptr->m_graphs);
        const auto& grids  =  (m_meta_ptr->m_grids);
        auto dit = grids.dataIterator();
        
        for(int ibox =0; ibox < dit.size(); ibox++)
        {
          const auto& valid           =            grids[dit[ibox]];
          ch_box intersection = valid & disp_box;
          if(!intersection.isEmpty())
          {
            const auto& graph           =           graphs[dit[ibox]];
            const auto& phi_fab         =     (*a_phi_ptr)[dit[ibox]];

            ebcm_subvol_vec volvec(graph, intersection,  killing_me_smalls);
            for(int ivol = 0; ivol < volvec.size(); ivol++)
            {
              const auto& volu = volvec[ivol];

              //gets stencil from first face in that direction yadf
              shared_ptr<local_stencil_t>  local =
                m_oper_ptr->getVoluSideFluxStencil(volu, graph, flux_dir, flux_side, &phi_fab, killing_me_smalls);
              double inhomog = local->m_inhomogeneous_bit;
              double flux_val = inhomog;
              for(int ivec = 0; ivec < local->m_entries.size(); ivec++)
              {
                pr_pt      sten_iv  = (local->m_entries)[ivec].first;
                double     wght     = (local->m_entries)[ivec].second;
                double     phi_val  = phi_fab(sten_iv, 0);
                double     flux_incr = phi_val*wght;
                flux_val += flux_incr;
              } //end loop over stencil
              Chombo4::pout() << "writeDetails for " << volu.m_pt 
                              << ": inhomog = " << inhomog
                              << ", flux    = " << flux_val
                              << std::endl;
            } //end loop over volumes
          } //end if(intersection exists)
        } //end loop over boxes
      }//end function writeDetails
      ///
      void applyOpPrimitives(shared_ptr<ebcm_dou_sca_data> a_lph_ptr,
                             shared_ptr<ebcm_dou_sca_data> a_phi_ptr,
                             shared_ptr<ebcm_dou_sca_flux> a_flx_ptr,
                             string                        a_identifier) const
      {
        CH_TIME("helmholtz::applyOpPrimitives");

        const auto& graphs = *(m_meta_ptr->m_graphs);
        const auto& grids  =  (m_meta_ptr->m_grids);
        auto dit = grids.dataIterator();
        
        for(int ibox =0; ibox < dit.size(); ibox++)
        {
          const auto& valid           =            grids[dit[ibox]];
          const auto& graph           =           graphs[dit[ibox]];
          const auto& phi_fab         =     (*a_phi_ptr)[dit[ibox]];
          auto      & lph_fab         =     (*a_lph_ptr)[dit[ibox]];
          //auto      & inhomog_bit_fab = (*m_inhomog_bit)[dit[ibox]];
          bool killing_me_smalls = false;
          //begin debug
          //killing_me_smalls = true;
          //end debug
          ebcm_subvol_vec volvec(graph, valid, killing_me_smalls);
          for(int ivol = 0; ivol < volvec.size(); ivol++)
          {
            const auto& volu = volvec[ivol];
            //full  operator
            {
              shared_ptr<local_stencil_t>  local =
                m_oper_ptr->getLocalStencil(volu, graph);
              double lph_val = local->m_inhomogeneous_bit;
              for(int ivec = 0; ivec < local->m_entries.size(); ivec++)
              {
                pr_pt      sten_iv  = (local->m_entries)[ivec].first;
                double     wght     = (local->m_entries)[ivec].second;
                double     phi_val  = phi_fab(sten_iv, 0);
                double     lph_incr = phi_val*wght;
                lph_val += lph_incr;
              } //end loop over stencil
              lph_fab(volu.m_pt, 0) = lph_val;
            }
            //fluxes 
            for(int idir = 0; idir < DIM; idir++)
            {
              auto      & flx_mpi        =     (*a_flx_ptr).m_data[idir];
              auto      & flx_fab        =       (*flx_mpi)[dit[ibox]];
              //gets stencil from first face in that direction
              shared_ptr<local_stencil_t>  local =
                m_oper_ptr->getVoluSideFluxStencil(volu, graph, idir, Side::Hi, NULL, false);
              double flux_val = local->m_inhomogeneous_bit;
              for(int ivec = 0; ivec < local->m_entries.size(); ivec++)
              {
                pr_pt      sten_iv  = (local->m_entries)[ivec].first;
                double     wght     = (local->m_entries)[ivec].second;
                double     phi_val  = phi_fab(sten_iv, 0);
                double     flux_incr = phi_val*wght;
                flux_val += flux_incr;
              } //end loop over stencil
              flx_fab(volu.m_pt, 0) = flux_val;
            } //end loop over 
          } //end loop voer volvec
        } //end loop over boxes
      } //end function applyOpPrimitives
      
      ///
      void applyOp(shared_ptr<ebcm_dou_sca_data> a_lph_ptr,
                   shared_ptr<ebcm_dou_sca_data> a_phi_ptr,
                   bool                          a_homogeneousOnly) 
      {
        CH_TIME("PETSc_Framework::applyOp");
        //this stuff makes lph = laplacian
        VecZeroEntries(*m_phi_ptr);
        VecZeroEntries(*m_lph_ptr);

        chomboToPetsc(*m_phi_ptr, *a_phi_ptr            );
        MatMult(      *m_mat_ptr, *m_phi_ptr, *m_lph_ptr);
        petscToChombo(*a_lph_ptr, *m_lph_ptr            );
        if(!a_homogeneousOnly)
        {
          addInhomogeneousBit(a_lph_ptr);
        }
      }

      ///
      /** 
          These two functions translate between the PETSc data holder (Vec)
          and the Chombo data holder, ebcm_dou_sca_data for short.
      **/


      ///
      void chomboToPetsc(petsc_vec                          & a_petsc_data,
                         const ebcm_dou_sca_data            & a_chomb_data)
      {
        CH_TIME("BasePetsc::chomboToPetsc");
        
        const auto& graphs = *(m_meta_ptr->m_graphs);
        const auto& grids  =  (m_meta_ptr->m_grids);
        auto dit = grids.dataIterator();
        for(int ibox =0; ibox < dit.size(); ibox++)
        {
          const auto& valid     =           grids[dit[ibox]];
          const auto& graph     =          graphs[dit[ibox]];
          const auto& chomb_fab =    a_chomb_data[dit[ibox]];
          const auto& row_map   =(*m_row_map_ptr)[dit[ibox]];

          /// this gets all volumes in valid
          ebcm_subvol_vec subvol_vec(graph, valid, false);
          int ivar = 0; //scalar
          for(int ivec = 0; ivec < subvol_vec.size(); ivec++)
          {
            const auto& volu    = subvol_vec[ivec];
            double    data_val  = chomb_fab(volu.m_pt, 0);
            petsc_int petsc_loc =   row_map(volu.m_pt, 0);

            VecSetValues(a_petsc_data, 1, &petsc_loc, &data_val, INSERT_VALUES);
          } // end loop over volumes in valid
        }//    end loop over boxes
        ///Not sure what they do but the HOEB version had them.   Petsc can be inscrutable.
        VecAssemblyBegin(a_petsc_data);
        VecAssemblyEnd(  a_petsc_data);
      }

      ///
      void petscToChombo(ebcm_dou_sca_data                  & a_chomb_data,
                         const petsc_vec                    & a_petsc_data)
                         
      {
        CH_TIME("BasePetsc::petscToChombo");
        
        const auto& graphs = *(m_meta_ptr->m_graphs);
        const auto& grids  =  (m_meta_ptr->m_grids);
        const petsc_sca * raw_array;
        VecGetArrayRead(a_petsc_data, &raw_array);

        petsc_int local_start,  local_end, local_size;
        VecGetOwnershipRange( a_petsc_data ,&local_start,&local_end);
        VecGetLocalSize(      a_petsc_data ,&local_size);

        auto dit = grids.dataIterator();
        for(int ibox =0; ibox < dit.size(); ibox++)
        {
          const auto& valid     =           grids[dit[ibox]];
          const auto& graph     =          graphs[dit[ibox]];
          const auto& row_map   =(*m_row_map_ptr)[dit[ibox]];
          auto      & chomb_fab =    a_chomb_data[dit[ibox]];

          /// this gets all volumes in valid
          ebcm_subvol_vec subvol_vec(graph, valid, false);
          int ivar = 0; //scalar
          for(int ivec = 0; ivec < subvol_vec.size(); ivec++)
          {
            const auto& volu = subvol_vec[ivec];
            petsc_int global_petsc_loc =   row_map(volu.m_pt, 0);
            petsc_int      local_index = global_petsc_loc - local_start;
            double data_val            = raw_array[local_index];
            for(int ivec = 0; ivec < volu.m_cells.size(); ivec++)
            {
              chomb_fab(volu.m_cells[ivec], ivar)      =  data_val ;
            }
          } //end loop over volumes in box.
        }   //loop over boxes in grid
      } //end function petscToChombo

    private:
      /// banning weak construction but assignment and copy constructors should work fine.
      PETSc_Interface();
      
    protected:
      /**
         This makes m_row_map_ptr 
         (operator uses this to map between volumes and petsc matrix rows)
      **/
      PetscInt createMatrixRowMap()
      {
        CH_TIME("BasePetsc::createMatrixRowMap");
        m_row_map_ptr = shared_ptr<ebcm_int_sca_data>(new ebcm_int_sca_data(m_meta_ptr, m_ghost));
        m_inhomog_bit = shared_ptr<ebcm_dou_sca_data>(new ebcm_dou_sca_data(m_meta_ptr, m_ghost));
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

          ebcm_subvol_vec subvollocal(graph, valid, false);
          numPtsThisProc     += subvollocal.size();
          numPtsPerBox[ibox]  = subvollocal.size();
        }

        //Chombo4::pout() << "Base_Petsc_Op::createMatrixRowMap: numptsThisProc = " << numPtsThisProc << endl;

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
        for(int iproc = 0; iproc <  CH4_SPMD::procID(); iproc++)
        {
          startgid += numPtsAllProcs[iproc];
        }
        //Chombo4::pout()  << "Base_Petsc_Op::getMap: procID = " << CH4_SPMD::procID() << ", numPtsThisProc = " << numPtsThisProc << ", startgid = " << startgid <<", totalNumPts = " << totalNumPts <<  endl;
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
          ebcm_subvol_vec volvec(graph, valid, false);
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

        //Chombo4::pout() << "Base_Petsc_Op::getMap calling exchange" << endl;

        m_row_map_ptr->exchange(false);
        
        return 0;
      }//end function createMatrixRowMap

      ///
      /**
         Makes m_mat_ptr (operator matrix) 
         and 
         m_phi_ptr, m_lph_ptr, vectors to hold petsc copies of solution, rhs.
      **/
      PetscInt createOperatorPetscMatrix()
      {
        CH_TIME("Base_Petsc_Op::createOperatorPetscMatrix");
        // create matrix
//        PetscInt nnzrow = 0;
#ifdef CH_MPI
        MPI_Comm wcomm = Chombo_MPI::comm;
#else
        MPI_Comm wcomm = PETSC_COMM_SELF;
#endif
        PetscInt ierr;
        m_mat_ptr = shared_ptr<petsc_mat>(new petsc_mat());
        m_phi_ptr = shared_ptr<petsc_vec>(new petsc_vec());
        m_lph_ptr = shared_ptr<petsc_vec>(new petsc_vec());
        
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
        ierr = formMatrix();CHKERRQ(ierr);
        return 0;
      }

      ///sets the actual values in *m_mat_ptr
      PetscInt formMatrix()
      {
    
        CH_TIME("Base_Petsc_Op::formMatrix");
        
        const auto& graphs = *(m_meta_ptr->m_graphs);
        const auto& grids  =  (m_meta_ptr->m_grids);
        auto dit = grids.dataIterator();

        //get all the info on points on this proc.
        //int numPtsThisProc = 0;
        vector<int> numPtsPerBox(dit.size(), 0);
        for(int ibox =0; ibox < dit.size(); ibox++)
        {
          const auto& valid           =            grids[dit[ibox]];
          const auto& graph           =           graphs[dit[ibox]];
          auto      & inhomog_bit_fab = (*m_inhomog_bit)[dit[ibox]];
          bool killing_me_smalls = false;
          //begin debug
          //killing_me_smalls = true;
          //end debug
          ebcm_subvol_vec volvec(graph, valid, killing_me_smalls);
          for(int ivol = 0; ivol < volvec.size(); ivol++)
          {
            const auto& volu = volvec[ivol];
            shared_ptr<local_stencil_t>  local =
              m_oper_ptr->getLocalStencil(volu, graph);
            double inhombit = local->m_inhomogeneous_bit;
            inhomog_bit_fab(volu.m_pt, 0)   = inhombit;
            //int local_size = local->size();
            const auto    & maplev = *(m_row_map_ptr->m_data);
            const auto    & mapfab = maplev[dit[ibox]];
            PetscInt irow = mapfab(volu.m_pt, 0);
            for(int isten = 0; isten < local->size(); isten++)
            {
              const auto& stenpair = (*local)[isten];
              PetscInt icol = mapfab(stenpair.first, 0);
              double weight = stenpair.second;
              /// This is PETSc for
              /// (*m_mat)(irow, icol) = weight;
              MatSetValues(*(this->m_mat_ptr),1,&irow,1,&icol,&weight,INSERT_VALUES);
            }
          }
        }
        PetscInt ierr;
        ierr = MatAssemblyBegin(*(this->m_mat_ptr),MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
        ierr = MatAssemblyEnd(  *(this->m_mat_ptr),MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
        ///now for the vectors
#ifdef CH_MPI
        MPI_Comm wcomm = Chombo_MPI::comm;
#else
        MPI_Comm wcomm = PETSC_COMM_SELF;
#endif
        Vec& phiVec = *m_phi_ptr;
        Vec& lphVec = *m_lph_ptr;
        VecCreate(   wcomm, &phiVec                  );
        VecSetSizes(         phiVec, m_numPtsThisProc, m_numPtsAllProcs);
        VecSetFromOptions(   phiVec                  );

        ierr = VecAssemblyBegin(phiVec );  CHKERRQ(ierr);
        ierr = VecAssemblyEnd(  phiVec );  CHKERRQ(ierr);
        
        VecDuplicate(        phiVec, &lphVec         );
        ierr = VecAssemblyBegin(lphVec );  CHKERRQ(ierr);
        ierr = VecAssemblyEnd(  lphVec );  CHKERRQ(ierr);
        
        return 0;
      } //end function formMatrix

    }; //end class petscinterface


    /**
       sets d = facedir.
       (1) Takes face moments (in d) and translates to volume moments.
       (2) It shifts these moments to xbar from locality
       (3) resizes output matrix to (1, pr_mom_dim::size())
       (4) matrix_value(0, currentCol) = pd*derval/momzero

       momzero == the face area. 
       pd      ==  
       matval = 0; 
       if( momzero > eps)
       {
       pr_itm_i_dim derind = momind - pr_itm_i_dim::BASISV(a_facedir);
       Real derval = volmom[derind];
       matval = pd*derval/momzero;
       }
       retval(0, curCol) = matval
    **/
    template <CENTERING cent>
    static shared_ptr<eigen_mat>
    getHelmholtzQMatrix( const EBCM_Face<cent, ebcm_order>  & a_face,
                         shared_ptr<ebcm_neighborhood>        a_locality,
                         const int                          & a_facedir,
                         const double                       & a_dx)
    {
      pr_mom_dmo facmom =  a_face.m_facmom;
      pr_mom_dim volmom = ebcm_basics::getVolMomFromFacMom(facmom, a_facedir);
      //double momzero    = facmom[pr_itm_i_dmo::Zero];

      pr_itm_r_dim xbar = a_locality->getXBar(a_face, a_dx);
      volmom.shift(xbar);

      unsigned int Pcols = pr_mom_dim::size();
      shared_ptr<eigen_mat> retval(new eigen_mat(1, Pcols));

      //Real area = a_face.m_facmom[pr_itm_i_dmo::Zero];
      for(MomentIterator<DIM, ebcm_order> momit; momit.ok(); ++momit)
      {
        auto p_mom = momit();
        double matval = 0;
        int pd = p_mom[a_facedir];
        if(pd == 1) //pd > 1 means x-xbar^p = 0
        {
          pr_itm_i_dim p_minu_ed = p_mom - pr_itm_i_dim::BASISV(a_facedir);
          double  m_p_min_ed = volmom[p_minu_ed];
          matval = pd*m_p_min_ed;
          /**
          if(area > 1.0e-11)
          {
            matval /= area;
          }
          **/
        }
        unsigned int currentCol = pr_mom_dim::indexOf(p_mom);
        (*retval)(0,currentCol) = matval;
      }
    
      return retval;
    }///end function getHelmholtzQMatrix


    ///
    static double
    getDomainBCValue(const ebcm_dom_face     & a_domain_face,
                     const string            & a_dom_bc,
                     const ebcm_exac         & a_bc_exact,
                     const ebcm_neighborhood & a_neigh,
                     bool a_killing_me_smalls)
    {
      double retval = 4586;
      double area = a_domain_face.m_facmom[pr_itm_i_dmo::Zero];
      if(     a_dom_bc == string("dirichlet"))
      {
        retval = a_bc_exact.getIntPhiDADomainFace(a_domain_face, a_neigh, a_killing_me_smalls);
      }
      else if(a_dom_bc == string("neumann"))
      {
        retval = a_bc_exact.getIntGradPhiDotNDADomainFace(a_domain_face, a_neigh);
      }
      else
      {
        pout() << "getDomainBCValue knows not of this dom_bc = " << a_dom_bc << endl;
        ch_mayday::Error();
      }
      //we need values not integrals
      if(area > 1.0e-12)
      {
        retval /= area;
      }
      if(a_killing_me_smalls)
      {
        Chombo4::pout() << "getDomainBCValue: area = " << area << ", dom_bc = " << a_dom_bc << ", bcval = " << retval << endl;
      }
      return retval;
    }//end function getDomBCValue

    ///
    static double
    getEBBCValue(const ebcm_cut_face     & a_cut_face,
                 const string            & a_eb_bc,
                 const ebcm_exac         & a_bc_exact,
                 const ebcm_neighborhood & a_neigh,
                 const double            & a_dx)
    {
      double retval = 4586;

      if(     a_eb_bc == string("dirichlet"))
      {
        retval = a_bc_exact.getIntPhiDA(a_cut_face, a_neigh);
      }
      else if(a_eb_bc == string("neumann"))
      {
        retval = a_bc_exact.getIntGradPhiDotNDA(a_cut_face, a_neigh, a_dx);
      }
      else
      {
        pout() << "getEBBCValue knows not of this eb_bc = " << a_eb_bc << endl;
        ch_mayday::Error();
      }
      
      return retval;
    } //end function getEBBCValue

    ///Returns a local stencil for \int (\grad \phi) \cdot (\ehat^facedir)
    /**
       Get stencil for integrated gradient in a  coordinate direction. 
       This is the same as the average gradient * face area.
       \int (\grad \phi) \cdot (\ehat^d) dA = A_f < \frac{\partial \phi}{\partial x_d} >
    */
    static local_stencil_t
    getStencilFromMat(const eigen_mat                                & a_matrix,
                      shared_ptr< EBCM::neighborhood< ebcm_order > > & a_locality,
                      const string                                   & a_dom_bc,
                      const string                                   & a_eb_bc,
                      const ebcm_exac                                & a_bc_exact,
                      const double                                   & a_dx)
                      
    {
      int nrows_volumes   = a_locality->m_weig_volumes.size();
      local_stencil_t retval;
      const auto& volvec = *a_locality->m_volumes;
      int mat_row = 0;
      for( int irow = 0; irow < nrows_volumes; irow++)
      {
        pair<pr_pt, double> entry;
        const ebcm_volu& volume = volvec[irow];
        entry.first = volume.m_pt;
        entry.second = a_matrix(0, mat_row);
        retval.add(entry);
        mat_row++;
      }
      //Inhomogeneous bit
      ///we have m_bc_exact to get values or derivatives. 
      retval.m_inhomogeneous_bit = 0;
      const auto &  dom_face_vec =*(a_locality->m_dom_faces);
      for(int idom = 0; idom < dom_face_vec.size(); idom++)
      {
        double      dom_weig     = a_matrix(0, mat_row);
        const auto& dom_face     = dom_face_vec[idom];
        double      dom_face_val = getDomainBCValue(dom_face, a_dom_bc, a_bc_exact, *a_locality, false);
        retval.m_inhomogeneous_bit += dom_weig*dom_face_val;
        mat_row++;
      }
      const auto &  cut_face_vec = *(a_locality->m_cut_faces);
      for(int icut = 0; icut < cut_face_vec.size(); icut++)
      {
        double      cut_weig     = a_matrix(0, mat_row);
        const auto& cut_face     = cut_face_vec[icut];
        double      cut_face_val = getEBBCValue(cut_face, a_eb_bc,  a_bc_exact, *a_locality, a_dx);
        retval.m_inhomogeneous_bit += cut_weig*cut_face_val;
        mat_row++;
      }

      return retval;
    } //end function getStencilFromMat

    ///
    /**
       Fill a row vector with exact values of the solution and bc values.
       This is the Phi that we usually factor out.
       
       Used to evaluate intermediate values for debugging purposes.
       Phi = {phi_v, phi_db, phi_eb}^T
       phi_v = Taylor series over volumes: phi_v = M C
       phi_db = Taylor series over domain faces: phi_db = M C 
       (if bcs are Neumann, phi_db is the gradient compoenent in the face normal direction).
       phi_eb = Taylor series over cut faces: phi_eb = M C 
       (if bcs are Neumann, phi_db = grad phi dot n_eb).

     **/
    static shared_ptr<eigen_mat>
    getSolutionRowVector(const shared_ptr<ebcm_neighborhood> & a_locality,
                         const string                        & a_dom_bc,
                         const string                        & a_eb_bc,
                         const ebcm_exac                     & a_bc_exact,
                         const double                        & a_dx,
                         bool a_killing_me_smalls)
    {
      CH_TIME("PETSc_Framework::getSolutionRowVector");
      int n_rows_volumes   = a_locality->m_weig_volumes.size();
      int n_rows_dom_faces = a_locality->m_weig_dom_faces.size();
      int n_rows_cut_faces = a_locality->m_weig_cut_faces.size();
      int n_rows           = n_rows_volumes + n_rows_dom_faces + n_rows_cut_faces;

      if(a_killing_me_smalls)
      {
        Chombo4::pout() << "getSolnRowVec: "
                        << " nvol = " << n_rows_volumes       
                        << ",ndom = " << n_rows_dom_faces
                        << ",ncut = " << n_rows_cut_faces
                        << ",ntot = " << n_rows << endl;
      }
      shared_ptr<eigen_mat> soln_vec_p(new eigen_mat(n_rows, 1));
      soln_vec_p->setVal(4586.);
      int irow_global = 0;
      for(int irow = 0; irow < n_rows_volumes; irow++)
      {
        const ebcm_volu& volu  =   (*(a_locality->m_volumes ))[irow];
        double phi_val = a_bc_exact.getCellAveSolution(volu);
        (*soln_vec_p)(irow_global,0) = phi_val;
        if(a_killing_me_smalls)
        {
          Chombo4::pout() << "getSolnRowVec: (vol) phi[" << irow_global << "]=" << phi_val << endl;
        }
        irow_global++;
      }

      for(int irow = 0; irow < n_rows_dom_faces; irow++)
      {
        const ebcm_dom_face& face = (*(a_locality->m_dom_faces ))[irow];
        double phi_val  = getDomainBCValue(face, a_dom_bc, a_bc_exact, *a_locality, a_killing_me_smalls);
        (*soln_vec_p)(irow_global,0) = phi_val;
        if(a_killing_me_smalls)
        {
          Chombo4::pout() << "getSolnRowVec: (dom) phi[" << irow_global << "]=" << phi_val << endl;
        }
        irow_global++;
      }

      for(int irow = 0; irow < n_rows_cut_faces; irow++)
      {
        const ebcm_cut_face& face = (*(a_locality->m_cut_faces ))[irow];
        double phi_val  = 4586.;
        double area = face.m_ebfmom[pr_itm_i_dim::Zero];
        if(a_eb_bc == string("neumann"))
        {
          phi_val = a_bc_exact.getIntGradPhiDotNDACutFace(face, *a_locality, a_dx);
        }
        else if(a_eb_bc == string("dirchlet"))
        {
          phi_val = a_bc_exact.getIntPhiDACutFace(face, *a_locality);
        }
        else
        {
          ch_mayday::Error("PETScFramework::getSolutionVector::bogus a_dom_bc");
        }
        phi_val /= area;
        (*soln_vec_p)(irow_global,0) = phi_val;
        if(a_killing_me_smalls)
        {
          Chombo4::pout() << "getSolnRowVec: (cut) phi[" << irow_global << "]=" << phi_val << endl;
        }
        irow_global++;
      }
      return soln_vec_p;
   }
    

    ///
    template<CENTERING cent>
    static local_stencil_t
    getHelmholtzFDotEhatStencil(const EBCM_Face<cent, ebcm_order>  & a_face,
                                const int                          & a_stenrad,
                                const int                          & a_weightpower,
                                const EBCM::EBCM_Graph<ebcm_order> & a_graph,
                                const string                       & a_dom_bc,
                                const string                       & a_eb_bc,
                                const int                          & a_facedir,
                                const double                       & a_dx,
                                const ebcm_exac                    & a_bc_exact,
                                const pr_fab_t                     * a_phi_exac_ptr, //can be null
                                bool a_killing_me_smalls)
    {
      CH_TIME("PETSc_Framework::getHelmholtzFDotEhatStencil");
      /** 
          At domain faces, Neumann bcs do not get a stencil but they can have an inhomogeneous value.
      **/
      local_stencil_t retval;

      shared_ptr< EBCM::neighborhood< ebcm_order > > locality;
      shared_ptr< eigen_mat > Amat_p, Wmat_p, Mmat_p, WMTmat_p;
      EBCM::Algorithm_Framework<ebcm_order>::
        getSystemInfo(Amat_p, Wmat_p, Mmat_p, WMTmat_p, locality,
                      a_face, a_graph, a_stenrad, a_weightpower,
                      a_dom_bc, a_eb_bc, false);
      
      if((a_dom_bc == string("neumann")) && (a_face.domainBoundary()))
      {
        if(a_killing_me_smalls)
        {
          Chombo4::pout() << "getHelmFDEHStencil: neumann bc" << endl;
        }
        retval = local_stencil_t();
        retval.clear();
        retval.m_inhomogeneous_bit = 0;

        double flux_val = a_bc_exact.getIntGradPhiDotNDA(a_face, a_dx);
        retval.m_inhomogeneous_bit += flux_val;
      }
      else
      {    
        /**
           From the document (Lawrence Berkeley National Laboratory Technical Report LBNL-2001537):  
           Q : Q matrix (Unfortunately, I did not give it a fancy mathy name in the document.
           It is the one that converts  Taylor coefficients to fluxes).   
           S : stencil matrix.  
           P : phi vector.  
           W : weighting matrix.
           M : moment matrix.
           A : system matrix =  ((WM)^T(WM))
           (S P) = <\nabla \psi> = Q C= (Q ((WM)^T(WM))^{-1} (WM)^T  P) Q A^{-1} WM^T P.
           S = Q C = (Q ((WM)^T(WM))^{-1} (WM)^T  ) = Q A^{-1} (WM)^T.
        **/
        //if(a_killing_me_smalls)
        //{
        //  Chombo4::pout() << "getHelmFDEHStencil: doing the Taylor dance" << endl;
        //}


        eigen_mat Ainvmat = *Amat_p;
        Ainvmat.invert();
        shared_ptr<eigen_mat> Qmat_p = getHelmholtzQMatrix(a_face, locality, a_facedir, a_dx);
      
        eigen_mat WMTW;
        multiply( WMTW, *WMTmat_p, *Wmat_p);
        eigen_mat AinvWMTW;
        multiply( AinvWMTW, Ainvmat, WMTW);
        eigen_mat SmatFFS;
        multiply( SmatFFS, *Qmat_p, AinvWMTW);

        retval = getStencilFromMat(SmatFFS, locality, a_dom_bc, a_eb_bc, a_bc_exact, a_dx);

//        ///begin debug
        static int s_num_prints = 0;
        static int s_max_prints = 4;
        if(a_killing_me_smalls && (s_num_prints < s_max_prints))
        {
          s_num_prints++; 
          int num_volumes = locality->m_volumes->size();
          int num_powers  = pr_mom_dim::size();
          pout() << "getHemmFDEHStencil: num_volumes = " << num_volumes << ", num_powers = " << num_powers << endl;

          /**
             Create a matrix with the solution in it so I can 
             understand what which line of my 
             obviously flawless code 
             is betraying me.
          **/
          shared_ptr<eigen_mat> Phi_p = getSolutionRowVector(locality, a_dom_bc, a_eb_bc, a_bc_exact, a_dx, a_killing_me_smalls);
//          pout() << "getHemmFDEHStencil: phi_p FFS = " << endl;
//          Phi_p->poutAll();
//          eigen_mat Smat = SmatFFS;


//          eigen_mat FmatFFS;
//          multiply( FmatFFS, SmatFFS, *Phi_p);
//          pout() << "getHemmFDEHStencil: FFS SmatFFS = " << endl;
//          SmatFFS.poutAll();
//          pout() << "getHemmFDEHStencil: FFS variables: FmatFFS = " << endl;
//          FmatFFS.poutAll();

//          eigen_mat Fmat;
//          multiply( Fmat, Smat, *Phi_p);
//          eigen_mat WPhi;
//          multiply( WPhi, *Wmat_p, *Phi_p);
//          eigen_mat WMTWPhi;
//          multiply( WMTWPhi, *WMTmat_p, WPhi);
//          eigen_mat taylor_coef;  
//          multiply( taylor_coef, Ainvmat, WMTWPhi);
//
//          ///
//          eigen_mat MC; //should give back phi
//          multiply( MC , * Mmat_p, taylor_coef);
//
//          eigen_mat QC; //should give the flux
//          multiply(QC, *Qmat_p, taylor_coef);
//
//          eigen_mat sten_phi;
//          multiply(sten_phi, Smat, *Phi_p);
//          Chombo4::pout() << "Smat PhiExact should be a scalar = int F dA:  "<< endl;
//          sten_phi.poutAll();
//          
//          Chombo4::pout() << "Flux matrix (should be a scalar = int F dA):  "<< endl;
//          Fmat.poutAll();
//          Chombo4::pout() << "Taylor coefficients: "<< endl;
//          taylor_coef.poutAll();
//
//          Chombo4::pout() << "MC (should give back Phi): "<< endl;
//          MC.poutAll();
//
//          Chombo4::pout() << "QC (should give back  int F dA):" << endl;
//          QC.poutAll();

        }
        /// end  debug
        
      }
      return retval;
    } //end function getHelmholtzFDotEhatStencil
    ///
    /**
       loops through d
       (2) It shifts these moments to xbar from locality
       (3) resizes output matrix to (1, pr_mom_dim::size())
       (4) matrix_value(0, currentCol) = pd*derval/momzero

       momzero == the face area. 
       pd      ==  
       matval = 0; 
       if( momzero > eps)
       {
       pr_itm_i_dim derind = momind - pr_itm_i_dim::BASISV(a_facedir);
       Real derval = volmom[derind];
       matval = pd*derval/momzero;
       }
       retval(0, curCol) = matval
    **/
    static shared_ptr<eigen_mat>
    getHelmholtzEBQMatrix( const IrregFace< ebcm_order >      & a_face,
                           shared_ptr<ebcm_neighborhood>        a_locality,
                           const double                       & a_dx)
    {
      CH_TIME("PETSc_Framework::getHelmholtzEBQMatrix");
      
      double momzero    = a_face.m_ebfmom[pr_itm_i_dim::Zero];
      unsigned int Pcols = pr_mom_dim::size();
      
      shared_ptr<eigen_mat> retval(new eigen_mat(1, Pcols));
      const pr_itm_r_dim xbar_face = a_locality->getXBarEB(a_face, a_dx);
      retval->setVal(0.); // needed
      for(int idir = 0; idir < DIM; idir++)
      {
        pr_mom_dim normal_moment;
        if(idir == 0)
        {
          normal_moment = a_face.m_norm_mom_x;
        }
        else if(idir == 1)
        {
          normal_moment = a_face.m_norm_mom_y;
        }
        else if(idir == 2)
        {
          normal_moment = a_face.m_norm_mom_z;
        }
        else
        { //can only be wetware error.
          ch_mayday::Error("logic error 4586");
        }
        normal_moment.shift(xbar_face);

        for(MomentIterator<DIM, ebcm_order> momit; momit.ok(); ++momit)
        {
          auto momind = momit();
          double matval = 0;
          if(momind[idir] > 0)
          {
            int pd = momind[idir];
            pr_itm_i_dim derind = momind - pr_itm_i_dim::BASISV(idir);
            Real derval = normal_moment[derind];
            if(std::abs(momzero) > 1.0e-10)
            {
              matval = pd*derval/momzero;
            }
          }
          unsigned int currentCol = pr_mom_dim::indexOf(momind);
          (*retval)(0,currentCol) += matval;//note: increments
        } //end loop over moments
      }   //end loop over directions
      return retval;

    } // end function getHelmholtzEBQMatrix
    
    ///Returns a local stencil for \int (\grad \phi) \cdot (\ehat^facedir)
    /**
       Get stencil for integrated gradient in a  coordinate direction. 
       This is the same as the average gradient * face area.
       \int (\grad \phi) \cdot (\ehat^d) dA = A_f < \frac{\partial \phi}{\partial x_d} >
    */
    static local_stencil_t
    getIrregHelmholtzFDotNStencil(const IrregFace<ebcm_order>        & a_face,
                                  const int                          & a_stenrad,
                                  const int                          & a_weightpower,
                                  const EBCM::EBCM_Graph<ebcm_order> & a_graph,
                                  const string                       & a_dom_bc,
                                  const string                       & a_eb_bc,
                                  const double                       & a_dx,
                                  const ebcm_exac                    & a_bc_exact,
                                  bool                                 a_killing_me_smalls)
    {
      CH_TIME("PETSc_Framework::getIrregHelmholtzFDotNStencil");
      /**
         There is no stencil for Neumann EBBC but it can have an inhomogeneous value.
      **/
      local_stencil_t retval;
      shared_ptr< EBCM::neighborhood< ebcm_order > > locality;
      shared_ptr<eigen_mat> Amat_p, Wmat_p, Mmat_p, WMTmat_p;
      algo_framework::getEBSystemInfo(Amat_p, Wmat_p, Mmat_p, WMTmat_p, locality,
                                      a_face, a_graph,
                                      a_stenrad, a_weightpower,
                                      a_dom_bc, a_eb_bc, false);
      ///Irreg faces are always boundary faces (so do not have to check)
      if(a_eb_bc == string("neumann")) 
      {
        retval = local_stencil_t();
        retval.clear();
        retval.m_inhomogeneous_bit = 0;
        double flux_val = a_bc_exact.getIntGradPhiDotNDA(a_face, a_dx);
        retval.m_inhomogeneous_bit += flux_val;
      }
      else
      {
        eigen_mat WMTmat  = *WMTmat_p;
        eigen_mat Ainvmat = *Amat_p;
        Ainvmat.invert();
      
        shared_ptr<eigen_mat> Qmat_p = getHelmholtzEBQMatrix(a_face, locality,  a_dx);
        eigen_mat AinvWMTmat;
        multiply( AinvWMTmat, Ainvmat, WMTmat);
        eigen_mat Smat;
        multiply( Smat, *Qmat_p, AinvWMTmat);
        retval = getStencilFromMat(Smat, locality, a_dom_bc, a_eb_bc, a_bc_exact, a_dx);
      }
      return retval;
    }// end function getIrregHelmholtzFDotNStencil
    
    ///
    /**
       Class which forms a helmholtz operator.
       This creates an operator of the form
       L(phi) = alpha phi + beta * div(F(phi))
    **/
    class HelmholtzOp
    {
    public:

      double                m_alpha;
      double                m_beta;
      int                   m_weight_power;
      int                   m_sten_rad;
      int                   m_ghost;
      string                m_dom_bc;
      string                m_eb_bc;
      string                m_identifier; // debugging tool.
      ///for inhomog
      shared_ptr<ebcm_exac> m_bc_exact;
      shared_ptr<ebcm_meta> m_meta_ptr;
        
      /// 
      /**
         alpha        = identity coeff
         beta         = coeff of laplcian
         ghost        = neighbor ghost size (grow the face box by this much)
         weight_power = goes to neighborhood to weight vofs.
         meta         = ebcm's very data intensive meta data
      **/
      HelmholtzOp(const double                                 & a_alpha,
                  const double                                 & a_beta,
                  const int                                    & a_sten_rad,
                  const int                                    & a_weight_power,
                  const shared_ptr< ebcm_meta  >               & a_meta_ptr,
                  shared_ptr<ebcm_exac> a_bc_exact,
                  string a_dom_bc, string a_eb_bc, int a_ghost, string a_identifier)
      {
        
        m_alpha        = a_alpha;       
        m_beta         = a_beta;
        m_sten_rad     = a_sten_rad;
        m_weight_power = a_weight_power;
        m_meta_ptr     = a_meta_ptr;
        m_bc_exact     = a_bc_exact;
        m_dom_bc       = a_dom_bc;
        m_eb_bc        = a_eb_bc;
        m_ghost        = a_ghost;
        m_identifier   = a_identifier;
      }

      ///this is the function that the driver code calls
      shared_ptr<local_stencil_t>
      getLocalStencil(const ebcm_volu & a_volu,
                      const ebcm_graph& a_graph)
      {
        shared_ptr<local_stencil_t> retval =
          privateVolumeStencil(a_volu,          
                               m_sten_rad,     
                               m_weight_power, 
                               a_graph,       
                               m_dom_bc,      
                               m_eb_bc);
        return retval;
      }

      ///in theory, just for debugging.  only grabs the last stencil on that side in the case of multiple faces.
      shared_ptr<local_stencil_t>
      getVoluSideFluxStencil(const ebcm_volu  & a_volu,
                             const ebcm_graph & a_graph,
                             const int        & a_idir,
                             Side::LoHiSide     a_side,
                             const pr_fab_t   * a_phi_exac_ptr, //can be null
                             bool               a_killing_me_smalls)
      {
        double  dx = a_volu.m_dx;
        int isign  = sign(a_side);
        
        shared_ptr<local_stencil_t> retval(new local_stencil_t());
        if(a_idir == 0)
        {
          if(a_killing_me_smalls)
          {
            Chombo4::pout() << "getVoluFluxStencil: getting x flux stencil " << endl;
          }
          for(int iface = 0; iface < a_volu.m_xfaces.size(); iface++)
          {
            const auto& face = a_volu.m_xfaces[iface];
            bool on_vol_exterior  = false;
            double face_sign = EBCM::faceSign<XFACE, ebcm_order>(on_vol_exterior, a_volu, face, 0, false);
            if(! on_vol_exterior)
            {
              ch_mayday::Error("volume not consistent in x");
            }
            if(face_sign*isign > 0)
            {
              *retval = getFluxStencil<XFACE>(face, a_graph, dx, a_phi_exac_ptr, a_killing_me_smalls);
            }
          }
        }
        else if(a_idir == 1)
        {
          if(a_killing_me_smalls)
          {
            Chombo4::pout() << "getVoluFluxStencil: getting y flux stencil " << endl;
          }
          for(int iface = 0 ;iface < a_volu.m_yfaces.size(); iface++)
          {
            const auto& face = a_volu.m_yfaces[iface];
            bool on_vol_exterior  = false;
            double face_sign = EBCM::faceSign(on_vol_exterior, a_volu, face, 1, false);
            if(! on_vol_exterior)
            {
              ch_mayday::Error("volume not consistent in y");
            }
            if(face_sign*isign > 0)
            {
              *retval = getFluxStencil<YFACE>(face, a_graph, dx, a_phi_exac_ptr, a_killing_me_smalls);
            }
          }
        }
        else if (a_idir == 2)
        {
          if(a_killing_me_smalls)
          {
            Chombo4::pout() << "getVoluFluxStencil: getting z flux stencil " << endl;
          }
          for(int iface = 0 ; iface < a_volu.m_zfaces.size(); iface++)
          {
            const auto& face = a_volu.m_zfaces[iface];
            bool on_vol_exterior  = false;
            double face_sign = EBCM::faceSign(on_vol_exterior, a_volu, face, 2, false);
            if(! on_vol_exterior)
            {
              ch_mayday::Error("volume not consistent in z");
            }
            if(face_sign*isign > 0)
            {
              *retval = getFluxStencil<ZFACE>(face, a_graph, dx, a_phi_exac_ptr, a_killing_me_smalls);
            }
          }
          MayDay::Warning("This diagnostic only written for x and y");
        }
        return retval;
      }

      template <CENTERING cent>
      local_stencil_t
      getFluxStencil(const EBCM_Face<cent, ebcm_order> & a_face,
                     const ebcm_graph                  & a_graph,
                     const double                      & a_dx, 
                     const pr_fab_t                    *  a_phi_exac_ptr, //can be null
                     bool a_killing_me_smalls) const
      {
        int face_dir = a_face.getFaceDir();
        local_stencil_t flux = getHelmholtzFDotEhatStencil(a_face,
                                                           m_sten_rad,     
                                                           m_weight_power, 
                                                           a_graph,       
                                                           m_dom_bc,      
                                                           m_eb_bc,       
                                                           face_dir, a_dx,
                                                           *m_bc_exact,
                                                           a_phi_exac_ptr, 
                                                           a_killing_me_smalls);
        return flux;
      
      }

    private:
      ///this is where the canal locks are dug in the basement
      shared_ptr<local_stencil_t>
      privateVolumeStencil(const ebcm_volu                    & a_volu, 
                           const int                          & a_stenrad,
                           const int                          & a_weightpower,
                           const EBCM::EBCM_Graph<ebcm_order> & a_graph,
                           const string                       & a_dom_bc,
                           const string                       & a_eb_bc)
      {
        
        shared_ptr<local_stencil_t> retval(new local_stencil_t());
        bool killing_me_smalls = false;

        {
          for(int iface = 0; iface < a_volu.m_xfaces.size(); iface++)
          {
            const auto& face =         a_volu.m_xfaces[iface];
            bool on_volu_exterior = false;
            double fluxsign = faceSign(on_volu_exterior, a_volu, face, 0, killing_me_smalls);
            if(!on_volu_exterior)
            {
              ch_mayday::Error("privateVolumeStencil: bogus face in x");
            }
            local_stencil_t flux = getHelmholtzFDotEhatStencil(face,
                                                               a_stenrad,     
                                                               a_weightpower, 
                                                               a_graph,       
                                                               a_dom_bc,      
                                                               a_eb_bc,       
                                                               0, a_volu.m_dx,
                                                               *m_bc_exact,
                                                               NULL,
                                                               killing_me_smalls);
            flux *= fluxsign;
            (*retval) += flux;
          }
        }
        {
          for(int iface = 0; iface < a_volu.m_yfaces.size(); iface++)
          {
            const auto& face =          a_volu.m_yfaces[iface];
            bool on_volu_exterior = false;
            double fluxsign  = faceSign(on_volu_exterior, a_volu, face, 1, killing_me_smalls);
            if(!on_volu_exterior)
            {
              ch_mayday::Error("privateVolumeStencil: bogus face in y");
            }
            local_stencil_t flux = getHelmholtzFDotEhatStencil(face,
                                                               a_stenrad,     
                                                               a_weightpower, 
                                                               a_graph,       
                                                               a_dom_bc,      
                                                               a_eb_bc,       
                                                               1, a_volu.m_dx,
                                                               *m_bc_exact,
                                                               NULL,
                                                               killing_me_smalls);

                                                               
            flux *= fluxsign;
            (*retval) += flux;
          }
        }
        {
          for(int iface = 0; iface < a_volu.m_zfaces.size(); iface++)
          {
            const auto& face =         a_volu.m_zfaces[iface];
            bool on_volu_exterior = false;
            double fluxsign = faceSign(on_volu_exterior, a_volu, face, 2, killing_me_smalls);
            if(!on_volu_exterior)
            {
              ch_mayday::Error("privateVolumeStencil: bogus face in z");
            }
            local_stencil_t flux = getHelmholtzFDotEhatStencil(face,
                                                               a_stenrad,     
                                                               a_weightpower, 
                                                               a_graph,       
                                                               a_dom_bc,      
                                                               a_eb_bc,
                                                               2, a_volu.m_dx,
                                                               *m_bc_exact,
                                                               NULL,
                                                               killing_me_smalls);
                                                               
            flux *= fluxsign;
            (*retval) += flux;
          }
        }
        {
          for(int iface = 0; iface < a_volu.m_ifaces.size(); iface++)
          {
            const auto& face     = a_volu.m_ifaces[iface];
            local_stencil_t flux = getIrregHelmholtzFDotNStencil(face,
                                                                 a_stenrad,     
                                                                 a_weightpower, 
                                                                 a_graph,       
                                                                 a_dom_bc,      
                                                                 a_eb_bc,
                                                                 a_volu.m_dx,
                                                                 *m_bc_exact,
                                                                 killing_me_smalls);
                                                                 
            (*retval) += flux;
          }
        }

        //now that we have int( f dot n) dA so we have to divide by volume magnitude
        double volmag = a_volu.m_volmom[pr_itm_i_dim::Zero];

        (*retval) *= (1./volmag);
        //now retval = laplacian  so we multiply by beta and add alpha I
        (*retval) *= m_beta;
        local_stencil_t ident_term;
        ident_term.clear();
        ident_term.m_inhomogeneous_bit = 0;
        ident_term.m_entries.resize(1);
        ident_term.m_entries[0].first  = a_volu.m_pt;
        ident_term.m_entries[0].second = m_alpha;

        (*retval) += ident_term;
        return retval;
      }

    }; //end class HelmholtzOp
    /**
       just parse stuff.  This decouples geometry generation with parmparse, allowing me to make varying refinements.
    **/ 
    static void parseFineGridParams(int            & a_nx,
                                    double         & a_dx,
                                    ch_probdom     & a_domain,
                                    int            & a_maxGrid,
                                    bool           & a_mergeSmallCells)
    {
      ParmParse pp("parseFineGridParams");

      pp.get("nx"             , a_nx);
      pp.get("maxGrid"        , a_maxGrid);
      pp.get("mergeSmallCells", a_mergeSmallCells);

      ch_iv domLo = ch_iv::Zero;
      ch_iv domHi  = (a_nx - 1)*ch_iv::Unit;
      a_domain = ch_probdom(domLo, domHi);
      a_dx = 1./(double(a_nx));
      
    }
    /**
       Make geometry at specified refinment.
       This decouples geometry generation with parmparse, allowing me to make varying refinements.
    **/ 
    static
    shared_ptr< ebcm_meta  >
    makeGeometryAtSpecificNX( const int            & a_nx,
                              const double         & a_dx,
                              const ch_probdom     & a_domain,
                              const int            & a_maxGrid,
                              const bool           & a_mergeSmallCells)
    {
      vector<ch_dbl> vecgrids;
      //Chombo4::pout() << "making grids" << endl;
      GeometryService<2>::generateGrids(vecgrids, a_domain.domainBox(), a_maxGrid);

      int geomGhost = 6;
      shared_ptr<BaseIF>    impfunc = ebcm_framework::getImplicitFunction();
      //Chombo4::pout() << "defining geometry in EB land" << endl;

      pr_rv origin = pr_rv::Zero();
      shared_ptr< ch_geoserv > geoserv
        (new ch_geoserv(impfunc, origin, a_dx, a_domain.domainBox(), vecgrids, geomGhost));

      int ghost = 4;
      shared_ptr< ebcm_meta  >
        metaDataPtr(new ebcm_meta(geoserv, a_domain.domainBox(), a_dx, ghost,
                                  a_mergeSmallCells, false));
      return metaDataPtr;
    }

    ///
    static void
    setExact( shared_ptr<ebcm_dou_sca_data> a_phi_exac_ptr,
              shared_ptr<ebcm_dou_sca_data> a_lph_exac_ptr,
              shared_ptr< ebcm_meta  >      a_meta_ptr,
              shared_ptr< ebcm_exac  >      a_exac_ptr,
              double a_alpha, double a_beta)
    {
      const auto& graphs = *(a_meta_ptr->m_graphs);
      const auto& grids  =  (a_meta_ptr->m_grids);
      auto dit = grids.dataIterator();
      for(int ibox =0; ibox < dit.size(); ibox++)
      {
        const auto& valid  =             grids[dit[ibox]];
        const auto& graph  =            graphs[dit[ibox]];
        auto& phi_fab      = (*a_phi_exac_ptr)[dit[ibox]];
        auto& lph_fab      = (*a_lph_exac_ptr)[dit[ibox]];

        /// this gets all volumes in valid
        ebcm_subvol_vec subvol_vec(graph, valid, false);
        //int ivar = 0; //scalar
        for(int ivec = 0; ivec < subvol_vec.size(); ivec++)
        {
          const auto& volu = subvol_vec[ivec];
//begin debug
//      double deb_loc_x = 0.5, deb_eps = 0.05;
//      int deb_stop = 0;
//      int    i_loc = volu.m_pt[0];
//      double x_loc = volu.m_dx*(double(i_loc) + 0.5);
//      if((x_loc > (deb_loc_x - deb_eps)) && (x_loc < (deb_loc_x + deb_eps)))
//      {
//        deb_stop = 1;
//      }
//end debug

          
          double phi_val   = a_exac_ptr->getCellAveSolution( volu);
          double lap_val   = a_exac_ptr->getCellAveLaplacian(volu);
          double lph_val   = a_alpha*phi_val + a_beta*lap_val;
          for(int ivec = 0; ivec < volu.m_cells.size(); ivec++)
          {
            phi_fab(volu.m_cells[ivec], 0) = phi_val;
            lph_fab(volu.m_cells[ivec], 0) = lph_val;
          }
        }
      }
    } // end function setExact
          
    static shared_ptr<ebcm_exac>
    getExactSolution(int a_exact_flag)
    {
      shared_ptr<ebcm_exac> retval;
      static bool printed_once = false;
      if(a_exact_flag == 0)
      {
        if(!printed_once) pout() << "exact solution = Zero" << endl;
        retval = shared_ptr<ebcm_exac>(new EBCM_Exact::ZeroSolution<ebcm_order>() );
      }
      else if(a_exact_flag==1)
      {
        using std::vector;
        using std::pair;
        vector<pair<Point, Real> > entries;
        int num_terms = -4586;

        ParmParse pp("polynomial_exact");
        pp.get("num_terms", num_terms);
        entries.resize(num_terms);
        if(!printed_once)  pout() << "exact solution = polynomial" << endl;
        for(int iterm = 0; iterm < num_terms; iterm++)
        {
          string coef_str   = string("coeff_") + to_string(iterm);
          string powers_str = string("power_") + to_string(iterm);
          vector<int> entry_powers;
          Real        entry_coef;
          pp.getarr(powers_str.c_str(), entry_powers, 0, DIM);
          pp.get(     coef_str.c_str(), entry_coef);
          pair<Point, Real> entry;
          entry.second = entry_coef;
          for(int idir = 0; idir < DIM; idir++)
          {
            entry.first[idir] = entry_powers[idir];
          }
          entries[iterm] = entry;
          if(!printed_once) pout() << "iterm = " << iterm << ", power = " << entry.first << ", weight = " << entry.second << endl;
        }//end loop over terms
        EBCM_Exact::PolynomialEF<ebcm_order>* ebcm_poly_ptr = new EBCM_Exact::PolynomialEF<ebcm_order>(entries);
        
        retval = shared_ptr<ebcm_exac>( ebcm_poly_ptr);
      }
      else if(a_exact_flag==2)
      {
        using std::vector;
        using std::pair;
        vector<pair<Point, Real> > entries;

        if(!printed_once)
        {
          pout() << "Exact solution is a polynomial of all the powers with unity coefficients." << endl;
          pout() << " (cp = 1/p! for all p <= (max_polynomial_order+1))." << endl;
        }
        ///Land-ho but for moments.
        typedef  Proto::MomentIterator<DIM ,ebcm_order+1> pr_momit_ho;
        typedef  Proto::IndexedMoments<DIM ,ebcm_order+1> pr_indmom_ho;
        int num_terms = pr_indmom_ho::size();
        entries.resize(num_terms);
        for(pr_momit_ho momit_ho; momit_ho.ok(); ++momit_ho)
        {
          pr_itm_i_dim power_itm = momit_ho();
          pr_pt power_pt;
          for(int idir = 0; idir < DIM; idir++)
          {
            power_pt[idir]  = power_itm[idir];
          }
          double  d_fact = double(power_itm.factorial());
          double  coeff = 1./d_fact;  
          int vec_ind = pr_indmom_ho::indexOf(power_itm);
          entries[vec_ind].first  = power_pt;
          entries[vec_ind].second = coeff;
          if(!printed_once)
          {
            pout() << "   power[" << vec_ind << "] = " << power_pt
                   << ", weight[" << vec_ind << "] = " << coeff  << endl;
          }
        }

        EBCM_Exact::PolynomialEF<ebcm_order>* ebcm_poly_ptr = new EBCM_Exact::PolynomialEF<ebcm_order>(entries);
        
        retval = shared_ptr<ebcm_exac>( ebcm_poly_ptr);
      }
      else if(a_exact_flag==3)
      {
        ParmParse pp("sine_exact");
        if(!printed_once)
        {
          pout() << "Exact solution = One dimensional sine." << endl;
        }
        int    frequency;
        int    direction;
        double amplitude;
        pp.get("frequency",frequency);
        pp.get("direction",direction);
        pp.get("amplitude",amplitude);

        EBCM_Exact::OneDimSineEF<ebcm_order>* ebcm_poly_ptr =
          new EBCM_Exact::OneDimSineEF<ebcm_order>(direction,
                                                   amplitude,
                                                   frequency);
        
        retval = shared_ptr<ebcm_exac>( ebcm_poly_ptr);
      }
      else
      {
        ch_mayday::Error("getExactSolution: bogus exact flag");
      }

      printed_once = true;
      return retval;
    }

    ///
    template<CENTERING cent>
    static double
    getFluxError(const EBCM_Face<cent, ebcm_order>      &   a_face,
                 const ebcm_graph                       &   a_graph,
                 const ebcm_volu                        &   a_volu,
                 const shared_ptr<HelmholtzOp>          &   a_helm_ptr,
                 const pr_fab_t                         &   a_phi_exac,
                 bool                                       a_killing_me_smalls)
    {
      local_stencil_t flux_stencil =
        a_helm_ptr->getFluxStencil(a_face, a_graph, a_volu.m_dx, &a_phi_exac, a_killing_me_smalls);

      double calc_val = flux_stencil.apply(a_phi_exac);
      double exac_val = a_helm_ptr->m_bc_exact->getIntGradPhiDotNDA(a_face, a_volu.m_dx);
      double face_area = a_face.m_facmom[0];
      if(face_area > 1.0e-10)
      {
        calc_val /= face_area;
        exac_val /= face_area;
      }
      return std::abs(calc_val - exac_val);
    }
    ///
    static void      
    getFluxErrorVolume(double                                     a_flux_error[DIM],
                       const ebcm_graph                       &   a_graph,
                       const ebcm_volu                        &   a_volu,
                       const shared_ptr<HelmholtzOp>          &   a_helm_ptr,
                       const pr_fab_t                         &   a_phi_exac,
                       bool a_killing_me_smalls)
    {
      for(int idir = 0; idir < DIM; idir++)
      {
        a_flux_error[idir] = -4586;
      }
      for(int iface = 0; iface < a_volu.m_xfaces.size(); iface++)
      {
        double face_error =
          getFluxError<XFACE>(a_volu.m_xfaces[iface],
                              a_graph, a_volu, a_helm_ptr, a_phi_exac,
                              a_killing_me_smalls);
        a_flux_error[0] = std::max(a_flux_error[0], std::abs(face_error));
      }
      for(int iface = 0; iface < a_volu.m_yfaces.size(); iface++)
      {

        double face_error =
          getFluxError<YFACE>(a_volu.m_yfaces[iface],
                              a_graph, a_volu, a_helm_ptr, a_phi_exac,
                              a_killing_me_smalls);
        a_flux_error[1] = std::max(a_flux_error[1], std::abs(face_error));
      }

#if DIM==3
      for(int iface = 0; iface < a_volu.m_zfaces.size(); iface++)
      {
        double face_error = getFluxError<ZFACE>(a_volu.m_zfaces[iface],
                                                a_graph, a_volu, a_helm_ptr, a_phi_exac,
                                                a_killing_me_smalls);
        a_flux_error[2] = std::max(a_flux_error[2], std::abs(face_error));
      }
#endif      
      
    }
    ///
    static void
    getMaxAbsFluxError( double                     a_max_flux_error[DIM],
                        shared_ptr< ebcm_meta  >   a_meta,
                        string                     a_file_prefix)
    {
      ParmParse pp("main");
      double alpha     = 4586.;
      double  beta     = 4586.;
      int sten_rad     = 4586 ;
      int weight_power = 4586 ;
      int exact_flag   = 4586 ;
      int output_hdf5  = 4586;
      string     dom_bc("4586");
      string      eb_bc("4586");

      pp.get("alpha"       , alpha   );
      pp.get("beta"        , beta    );
      pp.get("sten_rad"    , sten_rad);
      pp.get("weight_power", weight_power);
      pp.get("dom_bc"      , dom_bc);
      pp.get("eb_bc"       , eb_bc);
      pp.get("exact_flag"  , exact_flag);
      pp.get("output_hdf5" , output_hdf5);
      int num_ghost = 4;
      
      shared_ptr<ebcm_dou_sca_data> phi_exac_ptr(new ebcm_dou_sca_data(a_meta, num_ghost));
      shared_ptr<ebcm_dou_sca_data> lph_exac_ptr(new ebcm_dou_sca_data(a_meta, num_ghost));
      shared_ptr<ebcm_exac> exact_sol = getExactSolution(exact_flag);
      setExact(phi_exac_ptr, lph_exac_ptr, a_meta, exact_sol, alpha, beta);

      shared_ptr<HelmholtzOp> helm_ptr(new HelmholtzOp(alpha, beta, sten_rad,
                                                       weight_power, 
                                                       a_meta, 
                                                       exact_sol,
                                                       dom_bc, eb_bc, num_ghost, a_file_prefix));
      const auto& graphs = *(a_meta->m_graphs);
      const auto& grids  =  (a_meta->m_grids);
      auto dit = grids.dataIterator();

      for(int ibox =0; ibox < dit.size(); ibox++)
      {
        const auto& valid  =             grids[dit[ibox]];
        const auto& graph  =            graphs[dit[ibox]];
        auto& phi_fab      =   (*phi_exac_ptr)[dit[ibox]];

        /// this gets all volumes in valid
        ebcm_subvol_vec subvol_vec(graph, valid, false);

        for(int idir = 0; idir < DIM; idir++)
        {
          a_max_flux_error[idir] = -4586.;
        }
        
        for(int ivec = 0; ivec < subvol_vec.size(); ivec++)
        {
          const auto& volu = subvol_vec[ivec];
          double             flux_error[DIM];
          bool killing_me_smalls = false;
          getFluxErrorVolume(flux_error, graph, volu, helm_ptr, phi_fab, killing_me_smalls);
          for(int idir = 0; idir < DIM; idir++)
          {
            a_max_flux_error[idir] = std::max(a_max_flux_error[idir], std::abs(flux_error[idir]));
          }
        } //end loop over volumes
      } //end loop over boxes
    }
    ///
    static
    void
    getWindow( pr_pt& a_lo_coar,
               pr_pt& a_hi_coar,
               pr_pt& a_lo_fine,
               pr_pt& a_hi_fine)
    {
      static bool printed = false;
// lower left
//      if(!printed)
//      {
//        pout() << "displaying LOWER LEFT" << endl;
//        printed = true;
//      }      
//      pr_pt  lo_coar(D_DECL(0 ,  0, 0));
//      pr_pt  hi_coar(D_DECL(3,   3, 0));
//      pr_pt  lo_fine(D_DECL(0 ,  0, 0));
//      pr_pt  hi_fine(D_DECL(3,   3, 0));

//  upper middle      
//      if(!printed)
//      {
//        pout() << "displaying TOP MIDDLE" << endl;
//        printed = true;
//      }
//      pr_pt  lo_coar(D_DECL( 6, 12, 0));
//      pr_pt  hi_coar(D_DECL( 9, 15, 0));
//      pr_pt  lo_fine(D_DECL(14, 28, 0));
//      pr_pt  hi_fine(D_DECL(17, 31, 0));

//  upper left
      if(!printed)
      {
        pout() << "displaying TOP LEFT" << endl;
        printed = true;
      }
      pr_pt  lo_coar(D_DECL( 0, 12, 0));
      pr_pt  hi_coar(D_DECL( 3, 15, 0));
      pr_pt  lo_fine(D_DECL( 0, 28, 0));
      pr_pt  hi_fine(D_DECL( 3, 31, 0));
      
      
      a_lo_coar = lo_coar;
      a_hi_coar = hi_coar;
      a_lo_fine = lo_fine;
      a_hi_fine = hi_fine;
    }
    ///
    static
    void
    writeWindow (string                         a_file_prefix,
                 string                         a_data_label,
                 shared_ptr<ebcm_dou_sca_data>  a_data)
    {

      pr_pt lo_coar, hi_coar, lo_fine, hi_fine;
      getWindow(lo_coar, hi_coar, lo_fine, hi_fine);
      pr_pt lo_disp = lo_fine;
      pr_pt hi_disp = hi_fine;
      if(a_file_prefix == string("coar"))
      {
        lo_disp = lo_coar;
        hi_disp = hi_coar;
      }
      
      string prefix = a_file_prefix + string("_") + a_data_label;
      EBCM_Debug::DataWindow<double, 1, ebcm_order>::
        writeWindowLevel(&(*a_data), lo_disp, hi_disp, prefix);

    }

    ///
    static
    void
    writeWindow (string                         a_file_prefix,
                 string                         a_data_label,
                 shared_ptr<ebcm_dou_sca_flux>  a_data)
    {

      pr_pt lo_coar, hi_coar, lo_fine, hi_fine;
      getWindow(lo_coar, hi_coar, lo_fine, hi_fine);
      pr_pt lo_disp = lo_fine;
      pr_pt hi_disp = hi_fine;
      if(a_file_prefix == string("coar"))
      {
        lo_disp = lo_coar;
        hi_disp = hi_coar;
      }

      for(int idir = 0; idir < DIM; idir++)
      {
        string prefix = a_file_prefix + string("_") + a_data_label;
        EBCM_Debug::DataWindow<double, 1, ebcm_order>::
          writeFluxWindowLevel(&(*a_data), lo_disp, hi_disp, idir, prefix);
      }

    }
    static
    void
    writeStuffToFile(string                        a_file_prefix,
                     shared_ptr<ebcm_dou_sca_data> a_phi_exac_ptr,
                     shared_ptr<ebcm_dou_sca_data> a_lph_exac_ptr,
                     shared_ptr<ebcm_dou_sca_data> a_lph_calc_ptr,
                     shared_ptr<ebcm_dou_sca_data> a_err_calc_ptr)
    {
      ParmParse pp("main");
      int     output_hdf5;
      pp.get("output_hdf5" , output_hdf5);
      if(output_hdf5 != 0)
      {
        string phi_exac_file_name = a_file_prefix + string("_phi_exac.hdf5");
        string lph_exac_file_name = a_file_prefix + string("_lph_exac.hdf5");
        string lph_calc_file_name = a_file_prefix + string("_lph_calc.hdf5");
        string err_calc_file_name = a_file_prefix + string("_err_calc.hdf5");


        a_phi_exac_ptr->writeToHDF5(phi_exac_file_name);
        a_lph_exac_ptr->writeToHDF5(lph_exac_file_name);
        a_lph_calc_ptr->writeToHDF5(lph_calc_file_name);
        a_err_calc_ptr->writeToHDF5(err_calc_file_name);
      }
    }
    ///
    static
    void
    getMaxAbsLphiError( double                   & a_max_error,
                        shared_ptr< ebcm_meta  >   a_meta,
                        string                     a_file_prefix)
    {
      ParmParse pp("main");
      double alpha     = 4586.;
      double  beta     = 4586.;
      int sten_rad     = 4586 ;
      int weight_power = 4586 ;
      int exact_flag   = 4586 ;
      string     dom_bc("4586");
      string      eb_bc("4586");

      pp.get("alpha"       , alpha   );
      pp.get("beta"        , beta    );
      pp.get("sten_rad"    , sten_rad);
      pp.get("weight_power", weight_power);
      pp.get("dom_bc"      , dom_bc);
      pp.get("eb_bc"       , eb_bc);
      pp.get("exact_flag"  , exact_flag);
      int num_ghost = 4;
      
      shared_ptr<ebcm_dou_sca_data> phi_exac_ptr(new ebcm_dou_sca_data(a_meta, num_ghost));
      shared_ptr<ebcm_dou_sca_data> lph_exac_ptr(new ebcm_dou_sca_data(a_meta, num_ghost));
      shared_ptr<ebcm_dou_sca_data> lph_calc_ptr(new ebcm_dou_sca_data(a_meta, num_ghost));
      shared_ptr<ebcm_dou_sca_data> err_calc_ptr(new ebcm_dou_sca_data(a_meta, num_ghost));

      shared_ptr<ebcm_exac> exact_sol = getExactSolution(exact_flag);
      setExact(phi_exac_ptr, lph_exac_ptr, a_meta, exact_sol, alpha, beta);
      shared_ptr<HelmholtzOp> helm_ptr(new HelmholtzOp(alpha, beta, sten_rad,
                                                       weight_power, 
                                                       a_meta, 
                                                       exact_sol,
                                                       dom_bc, eb_bc, num_ghost, a_file_prefix));
                                                          

      ///using the proper petsc stuff
      PETSc_Interface<HelmholtzOp> discreteOp(helm_ptr, a_meta, num_ghost);

#if 0
      ///using the proper petsc stuff
      bool homogeneous_only = false;
      discreteOp.applyOp(lph_calc_ptr, phi_exac_ptr, homogeneous_only);
#else
      ///using stuff I think I actually understand
      discreteOp.writeDetails(phi_exac_ptr, a_file_prefix);
      return;
      shared_ptr<ebcm_dou_sca_flux> flx_calc_ptr(new ebcm_dou_sca_flux(a_meta, num_ghost));      
      discreteOp.applyOpPrimitives(lph_calc_ptr, phi_exac_ptr, flx_calc_ptr, a_file_prefix);
#endif
      //settting  error to lph_calc - lph_exac
      double acoef =  1;
      double bcoef = -1;
      err_calc_ptr->setToAxBy(acoef, *lph_calc_ptr, bcoef, *lph_exac_ptr);

      pr_pt pt_loc;
      double max_val;
      int i_co;
      err_calc_ptr->getAbsMaxLocationAndValue(pt_loc, max_val, i_co);
      a_max_error = max_val;

      writeWindow(a_file_prefix, string("phi_exac"), phi_exac_ptr);
      writeWindow(a_file_prefix, string("flx_calc"), flx_calc_ptr);
      writeWindow(a_file_prefix, string("lph_exac"), lph_exac_ptr);
      writeWindow(a_file_prefix, string("lph_calc"), lph_calc_ptr);
      writeWindow(a_file_prefix, string("err_calc"), err_calc_ptr);

      writeStuffToFile(a_file_prefix, phi_exac_ptr, lph_exac_ptr, lph_calc_ptr, err_calc_ptr);
      
      return;
    }
    ///
    static void
    getMaxAbsTruncation(double                   & a_lphi_err_max,
                        double                     a_flux_err_max[DIM],
                        shared_ptr< ebcm_meta  >   a_meta,
                        string                     a_file_prefix)
    {
      getMaxAbsFluxError(a_flux_err_max,   a_meta, a_file_prefix);
      getMaxAbsLphiError(a_lphi_err_max,   a_meta, a_file_prefix);
    }

    ///
    /**
       test function for example/EBCM/hoeb_truncation
    **/ 
    static  void run_hoeb_truncation_tests()
    {

      int           nx_fine = 4586;
      double        dx_fine = 4586.;
      ch_probdom domain_fine;
      ///these two have to be the same for fine and coar      
      int maxGrid = 4586;      bool mergeSmallCells = true;   
      parseFineGridParams(nx_fine, dx_fine, domain_fine, maxGrid, mergeSmallCells);
      int    nx_coar = nx_fine/2;
      double dx_coar = dx_fine*2.;
      ch_probdom domain_coar = domain_fine;
      domain_coar.coarsen(2);
      shared_ptr< ebcm_meta  > meta_fine = makeGeometryAtSpecificNX(nx_fine, dx_fine, domain_fine, maxGrid, mergeSmallCells);
      shared_ptr< ebcm_meta  > meta_coar = makeGeometryAtSpecificNX(nx_coar, dx_coar, domain_coar, maxGrid, mergeSmallCells);
      double max_lphi_err_fine = -4586;
      double max_lphi_err_coar = -4586;
      double max_flux_err_fine[DIM];
      double max_flux_err_coar[DIM];
      
      getMaxAbsTruncation(max_lphi_err_fine, max_flux_err_fine, meta_fine, string("fine")); 
      getMaxAbsTruncation(max_lphi_err_coar, max_flux_err_coar, meta_coar, string("coar"));


//      for(int itest = 0; itest < (DIM+1); itest++)
//      {
//        string which_test("4586");
//        double max_err_coar = -4586;
//        double max_err_fine = -4586;
//        if(itest == DIM)
//        {
//          which_test = string("L_phi");
//          max_err_coar = max_lphi_err_coar;
//          max_err_fine = max_lphi_err_fine;
//        }
//        else
//        {
//          which_test = string("Grad_phi[")  + to_string(itest) + string("]");
//          max_err_coar = max_flux_err_coar[itest];
//          max_err_fine = max_flux_err_fine[itest];
//        }
//
//        Chombo4::pout() << which_test<< ":err_coar_max=" << max_err_coar << endl;
//        Chombo4::pout() << which_test<< ":err_fine_max=" << max_err_fine << endl;
//        if(max_err_fine > 1.0e-12) //err >= 0  by construction
//        {
//          double order = (log(max_err_coar/max_err_fine))/log(2.0);
//          Chombo4::pout() << which_test << ":trunc_order=" << order << endl;
//        }
//        else
//        {
//          Chombo4::pout() <<  which_test << ":Fine error is too small for safe order computation." << endl;
//        }
//
//      } //end loop over tests
    } //  end function run_truncation_error_tests
  }; //   end petsc_framwork class
} //      end namespace EBCM

#endif //petsc
#endif  //include guard
