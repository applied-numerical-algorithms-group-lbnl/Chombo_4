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
    typedef Proto::Box                                      pr_bx;
    typedef Chombo4::Box                                    ch_box;
    typedef Chombo4::DisjointBoxLayout                      ch_dbl;
    typedef Chombo4::BoxIterator                            ch_bit;
    typedef Chombo4::IntVect                                ch_iv;
    typedef Chombo4::DataIterator                           ch_dit;

    typedef Chombo4::GeometryService<         ebcm_order >  ch_geoserv;
    typedef  EBCM::Algorithm_Framework<       ebcm_order >  algo_framework;
    typedef  EBCM::neighborhood       <       ebcm_order >  ebcm_neighborhood;

    typedef  EBCM::neig_dom_face                            ebcm_dom_face;
    typedef  EBCM::neig_cut_face                            ebcm_cut_face;
    
    typedef  EBCM::MetaDataLevel<             ebcm_order >  ebcm_meta;
    typedef  EBCM::Algorithm_Framework<       ebcm_order >  ebcm_framework;
    typedef  EBCM::EBCM_Graph<                ebcm_order >  ebcm_graph;
    typedef  EBCM::HostLevelData<int,    1,   ebcm_order >  ebcm_int_sca_data;
    typedef  EBCM::HostLevelData<double, 1,   ebcm_order >  ebcm_dou_sca_data;
    typedef  EBCM::HostLevelData<int,    DIM, ebcm_order >  ebcm_int_vec_data;
    typedef  EBCM::HostLevelData<double, DIM, ebcm_order >  ebcm_dou_vec_data;
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
         Having visited that circle already,  I am going to just allow inhomogenous BCs for my operators this time.
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
      
      virtual ~PETSc_Interface()
      {;}
  
      PETSc_Interface(const shared_ptr< linear_operator_t  >       & a_oper_ptr,
                      const shared_ptr< ebcm_meta  >               & a_meta_ptr,
                      const int                                    & a_ghost)
      {
        CH_TIME("EBPetscSolver::EBPetscSolver");
        m_oper_ptr     =    a_oper_ptr;
        m_meta_ptr     =    a_meta_ptr;
        m_ghost        =    a_ghost;
    
        Chombo4::pout() << "Base_Petsc_Op: creating map of locations in space to matrix row." << endl;
        createMatrixRowMap();
        
        Chombo4::pout() << "Base_Petsc_Op:: creating the matrix" << endl;

        createOperatorPetscMatrix();
      }


      ///
      void  addInhomogeneousBit(shared_ptr<ebcm_dou_sca_data> a_lph_ptr)
      {
        (*a_lph_ptr) += (*m_inhomog_bit);
      }
      
      ///
      void applyOp(shared_ptr<ebcm_dou_sca_data> a_lph_ptr,
                   shared_ptr<ebcm_dou_sca_data> a_phi_ptr,
                   bool                          a_homogeneousOnly) 
      {
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

        Chombo4::pout() << "Base_Petsc_Op::getMap calling exchange" << endl;

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
        PetscInt nnzrow = 0;
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
        int numPtsThisProc = 0;
        vector<int> numPtsPerBox(dit.size(), 0);
        for(int ibox =0; ibox < dit.size(); ibox++)
        {
          const auto& valid           =            grids[dit[ibox]];
          const auto& graph           =           graphs[dit[ibox]];
          auto      & inhomog_bit_fab = (*m_inhomog_bit)[dit[ibox]];
          ebcm_subvol_vec volvec(graph, valid, false);
          for(int ivol = 0; ivol < volvec.size(); ivol++)
          {
            const auto& volu = volvec[ivol];
            shared_ptr<local_stencil_t>  local =
              m_oper_ptr->getLocalStencil(volu.m_pt, graph);
            double inhombit = local->m_inhomogeneous_bit;
            inhomog_bit_fab(volu.m_pt, 0)   = inhombit;
            int local_size = local->size();
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
        VecDuplicate(        phiVec, &lphVec         );
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
      double momzero    = facmom[pr_itm_i_dmo::Zero];

      pr_itm_r_dim xbar = a_locality->getXBar(a_face, a_dx);
      volmom.shift(xbar);

      unsigned int Pcols = pr_mom_dim::size();
      shared_ptr<eigen_mat> retval(new eigen_mat(1, Pcols));


      for(MomentIterator<DIM, ebcm_order> momit; momit.ok(); ++momit)
      {
        auto momind = momit();
        double matval = 0;
        if(momind[a_facedir] > 0)
        {
          int pd = momind[a_facedir];
          pr_itm_i_dim derind = momind - pr_itm_i_dim::BASISV(a_facedir);
          Real derval = volmom[derind];
          if(std::abs(momzero) > 1.0e-10)
          {
            matval = pd*derval/momzero;
          }
        }
        unsigned int currentCol = pr_mom_dim::indexOf(momind);
        (*retval)(0,currentCol) = matval;
      }
    
      return retval;
    }///end function getHelmholtzQMatrix


    ///
    static double
    getDomainBCValue(const ebcm_dom_face     & a_domain_face,
                     const string            & a_dom_bc,
                     const ebcm_exac         & a_bc_exact,
                     const ebcm_neighborhood & a_neigh)
    {
      double retval = 4586;
      
      if(     a_dom_bc == string("Dirichlet"))
      {
        retval = a_bc_exact.getIntPhiDA(a_domain_face, a_neigh);
      }
      else if(a_dom_bc == string("Neumann"))
      {
        retval = a_bc_exact.getIntGradPhiDotNDA(a_domain_face, a_neigh);
      }
      else
      {
        pout() << "getDomainBCValue knows not of this dom_bc = " << a_dom_bc << endl;
        ch_mayday::Error();
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

      if(     a_eb_bc == string("Dirichlet"))
      {
        retval = a_bc_exact.getIntPhiDA(a_cut_face, a_neigh);
      }
      else if(a_eb_bc == string("Neumann"))
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
      for( int irow = 0; irow < nrows_volumes; irow++)
      {
        pair<pr_pt, double> entry;
        const ebcm_volu& volume = volvec[irow];
        entry.first = volume.m_pt;
        entry.second = a_matrix(irow, 0);
        retval.add(entry);
      }
      //Inhomogeneous bit
      ///we have m_bc_exact to get values or derivatives. 
      retval.m_inhomogeneous_bit = 0;
      const auto &  dom_face_vec =*(a_locality->m_dom_faces);
      const auto &  dom_weig_vec = (a_locality->m_weig_dom_faces);
      for(int idom = 0; idom < dom_face_vec.size(); idom++)
      {
        double      dom_weig     = dom_weig_vec[idom];
        const auto& dom_face     = dom_face_vec[idom];
        double      dom_face_val = getDomainBCValue(dom_face, a_dom_bc, a_bc_exact, *a_locality);
        retval.m_inhomogeneous_bit += dom_weig*dom_face_val;
      }
      const auto &  cut_face_vec = *(a_locality->m_cut_faces);
      const auto &  cut_weig_vec =   a_locality->m_weig_cut_faces;
      for(int icut = 0; icut < cut_face_vec.size(); icut++)
      {
        double      cut_weig     = cut_weig_vec[icut];
        const auto& cut_face     = cut_face_vec[icut];
        double      cut_face_val = getEBBCValue(cut_face, a_eb_bc,  a_bc_exact, *a_locality, a_dx);
        retval.m_inhomogeneous_bit += cut_weig*cut_face_val;
      }

      return retval;
    } //end function getStencilFromMat

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
                                const ebcm_exac                    & a_bc_exact)
    {
      /** 
          At domain faces, Neumann bcs do not get a stencil but they can have an inhomogeneous value.
      **/
      local_stencil_t retval;
      shared_ptr< EBCM::neighborhood< ebcm_order > > locality;
      shared_ptr< eigen_mat > Amat_p =
        EBCM::Algorithm_Framework<ebcm_order>::
        getSystemMatrix(locality, a_face, a_graph,
                        a_stenrad, a_weightpower,
                        a_dom_bc, a_eb_bc, false);
      if((a_dom_bc == string("Neumann")) && (a_face.domainBoundary()))
      {
        retval = local_stencil_t();
        retval.clear();
        retval.m_inhomogeneous_bit = 0;

        double flux_val = a_bc_exact.getIntGradPhiDotNDA(a_face, *locality, a_facedir, a_dx);
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
           A: system matrix =  ((WM)^T(WM))
           (S P) = <\nabla \psi> = Q C= (Q ((WM)^T(WM))^{-1} (WM)^T  P) Q A^{-1} WM^T P.
           S = Q C = (Q ((WM)^T(WM))^{-1} (WM)^T  ) = Q A^{-1} (WM)^T.
        **/

        shared_ptr<eigen_mat> Wmat_p = algo_framework::getWeightMatrix(locality,                   false);
        shared_ptr<eigen_mat> Mmat_p = algo_framework::getMomentMatrix(locality, a_dom_bc, a_eb_bc,false);

        eigen_mat WMmat;
        multiply(WMmat, *Wmat_p, *Mmat_p);
        eigen_mat WMTmat = WMmat;
        WMTmat.transpose();

        eigen_mat Ainvmat = *Amat_p;
        Ainvmat.invert();
      
        shared_ptr<eigen_mat> Qmat_p = getHelmholtzQMatrix(a_face, locality, a_facedir, a_dx);
        eigen_mat AinvWMTmat;
        multiply( AinvWMTmat, Ainvmat, WMTmat);
        eigen_mat Smat;
        multiply( Smat, *Qmat_p, AinvWMTmat);
        retval = getStencilFromMat(Smat, locality, a_dom_bc, a_eb_bc, a_bc_exact, a_dx);
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
                                  const ebcm_exac                    & a_bc_exact)
    {
        /**
           There is no stencil for Neumann EBBC but it can have an inhomogeneous value.
        **/
      local_stencil_t retval;
      shared_ptr< EBCM::neighborhood< ebcm_order > > locality;
      shared_ptr< eigen_mat > Amat_p =
        algo_framework::getEBSystemMatrix(locality, a_face, a_graph,
                                          a_stenrad, a_weightpower,
                                          a_dom_bc, a_eb_bc, false);
      ///Irreg faces are always boundary faces (so do not have to check)
      if(a_eb_bc == string("Neumann")) 
      {
        retval = local_stencil_t();
        retval.clear();
        retval.m_inhomogeneous_bit = 0;
        double flux_val = a_bc_exact.getIntGradPhiDotNDA(a_face, *locality, a_dx);
        retval.m_inhomogeneous_bit += flux_val;
      }
      else
      {
        shared_ptr< eigen_mat > Wmat_p = algo_framework::getWeightMatrix(  locality, false);
        shared_ptr< eigen_mat > Mmat_p = algo_framework::getMomentMatrix(  locality, a_dom_bc, a_eb_bc, false);

        eigen_mat WMmat;
        multiply(WMmat, *Wmat_p, *Mmat_p);
        eigen_mat WMTmat = WMmat;
        WMTmat.transpose();

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
                  string a_dom_bc, string a_eb_bc, int a_ghost)
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
      }

      virtual ~HelmholtzOp()
      {
      }

      ///this is the function that the driver code calls
      shared_ptr<local_stencil_t> getLocalStencil(const pr_pt& a_pt, const ebcm_graph& a_graph)
      {
        shared_ptr<local_stencil_t> retval =
          privateStencil(a_pt,          
                         m_sten_rad,     
                         m_weight_power, 
                         a_graph,       
                         m_dom_bc,      
                         m_eb_bc);
        return retval;
      }

    private:
      ///this is where the canal locks are dug in the basement
      virtual shared_ptr<local_stencil_t> privateStencil(const pr_pt                         & a_pt, 
                                                         const int                          & a_stenrad,
                                                         const int                          & a_weightpower,
                                                         const EBCM::EBCM_Graph<ebcm_order> & a_graph,
                                                         const string                       & a_dom_bc,
                                                         const string                       & a_eb_bc)
      {
        
        shared_ptr<local_stencil_t> retval(new local_stencil_t());
        pr_bx pointbx(a_pt, a_pt);
        pr_bx grbx = pointbx.grow(pr_pt::Ones(m_ghost));
        
        auto volume = a_graph.getVolumeCoveringPoint(a_pt);
        auto dx = volume.m_dx;
        {
          for(int iface = 0; iface < volume.m_xfaces.size(); iface++)
          {
            const auto& face = volume.m_xfaces[iface];
            double fluxsign = faceSign(volume, face);

            local_stencil_t flux = getHelmholtzFDotEhatStencil(face,
                                                               a_stenrad,     
                                                               a_weightpower, 
                                                               a_graph,       
                                                               a_dom_bc,      
                                                               a_eb_bc,       
                                                               0, dx,
                                                               *m_bc_exact);
            flux *= fluxsign;
            (*retval) += flux;
          }
        }
        {
          for(int iface = 0; iface < volume.m_yfaces.size(); iface++)
          {
            const auto& face = volume.m_yfaces[iface];
            double fluxsign = faceSign(volume, face);
            local_stencil_t flux = getHelmholtzFDotEhatStencil(face,
                                                               a_stenrad,     
                                                               a_weightpower, 
                                                               a_graph,       
                                                               a_dom_bc,      
                                                               a_eb_bc,       
                                                               1, dx,
                                                               *m_bc_exact);

                                                               
            flux *= fluxsign;
            (*retval) += flux;
          }
        }
        {
          for(int iface = 0; iface < volume.m_zfaces.size(); iface++)
          {
            const auto& face = volume.m_zfaces[iface];
            double fluxsign = faceSign(volume, face);
            local_stencil_t flux = getHelmholtzFDotEhatStencil(face,
                                                               a_stenrad,     
                                                               a_weightpower, 
                                                               a_graph,       
                                                               a_dom_bc,      
                                                               a_eb_bc,
                                                               2, dx,
                                                               *m_bc_exact);
                                                               
            flux *= fluxsign;
            (*retval) += flux;
          }
        }
        {
          for(int iface = 0; iface < volume.m_ifaces.size(); iface++)
          {
            const auto& face = volume.m_ifaces[iface];
            local_stencil_t flux = getIrregHelmholtzFDotNStencil(face,
                                                                 a_stenrad,     
                                                                 a_weightpower, 
                                                                 a_graph,       
                                                                 a_dom_bc,      
                                                                 a_eb_bc,
                                                                 dx,
                                                                 *m_bc_exact);
                                                                 
            (*retval) += flux;
          }
        }

        //now that we have int( f dot n) dA so we have to divide by volume magnitude
        double volmag = volume.m_volmom[pr_itm_i_dim::Zero];

        (*retval) *= (1./volmag);
        //now retval = laplacian  so we multiply by beta and add alpha I
        (*retval) *= m_beta;
        local_stencil_t ident_term;
        ident_term.clear();
        ident_term.m_inhomogeneous_bit = 0;
        ident_term.m_entries.resize(1);
        ident_term.m_entries[0].first  = a_pt;
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
      Chombo4::pout() << "making grids" << endl;
      GeometryService<2>::generateGrids(vecgrids, a_domain.domainBox(), a_maxGrid);

      int geomGhost = 6;
      shared_ptr<BaseIF>    impfunc = ebcm_framework::getImplicitFunction();
      Chombo4::pout() << "defining geometry in EB land" << endl;

      pr_rv origin = pr_rv::Zero();
      shared_ptr< ch_geoserv > geoserv
        (new ch_geoserv(impfunc, origin, a_dx, a_domain.domainBox(), vecgrids, geomGhost));

      int ghost = 0;
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
        int ivar = 0; //scalar
        for(int ivec = 0; ivec < subvol_vec.size(); ivec++)
        {
          const auto& volu = subvol_vec[ivec];
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
      ParmParse pp("polynomial_exact");
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
          entries.push_back(entry);
          if(!printed_once) pout() << "iterm = " << iterm << ", power = " << entry.first << ", weight = " << entry.second << endl;
        }//end loop over terms
        retval = shared_ptr<ebcm_exac>( new EBCM_Exact::PolynomialEF<ebcm_order>(entries));
      }
      else
      {
        ch_mayday::Error("getExactSolution: bogus exact flag");
      }

      printed_once = true;
      return retval;
    }
    ///
    static double
    getMaxAbsTruncation(shared_ptr< ebcm_meta  >   a_meta)
    {
      ParmParse pp("main");
      double alpha     = 4586.;
      double  beta     = 4586.;
      int sten_rad     = 4586 ;
      int weight_power = 4586;
      int exact_flag   = 4586;
      string dom_bc("4586");
      string  eb_bc("4586");

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
                                                       dom_bc, eb_bc, num_ghost));
                                                          
      
      PETSc_Interface<HelmholtzOp> discreteOp(helm_ptr, a_meta, num_ghost);

      bool homogeneous_only = false;
      discreteOp.applyOp(lph_calc_ptr, phi_exac_ptr, homogeneous_only);

      //settting  error to lph_calc - lph_exac
      double acoef =  1;
      double bcoef = -1;
      err_calc_ptr->setToAxBy(acoef, *lph_calc_ptr, bcoef, *lph_exac_ptr);
      double retval = err_calc_ptr->maxNorm();
      return retval;
    }

    ///
    /**
       test function for example/EBCM/hoeb_truncation
    **/ 
    static  void run_hoeb_truncation_tests()
    {

#ifdef CH_USE_PETSC
      static const int num_ghost = 4;  
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

      double max_err_fine = getMaxAbsTruncation(meta_fine); 
      double max_err_coar = getMaxAbsTruncation(meta_coar);

      Chombo4::pout() << "hoeb_truncation: Lmax(coarse error) =  " << max_err_coar << endl;
      Chombo4::pout() << "hoeb_truncation: Lmax(fine   error) =  " << max_err_fine << endl;
      if(max_err_fine > 1.0e-12) //err >= 0  by construction
      {
        double order = (log(max_err_coar/max_err_fine))/log(2.0);
        Chombo4::pout() << "hoeb_truncation: order of accuracy = " <<  order  << endl;
      }
      else
      {
        Chombo4::pout() << "hoeb_truncation: Fine error is too small for safe order computation." << endl;
      }
#endif
    }

  }; //end petsc_framwork class
}
#endif

#endif 
