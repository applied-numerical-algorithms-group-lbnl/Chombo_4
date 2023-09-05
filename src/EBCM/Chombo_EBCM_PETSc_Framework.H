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
    typedef Chombo4::Box                             ch_box;
    typedef Chombo4::DisjointBoxLayout               ch_dbl;
    typedef Chombo4::BoxIterator                     ch_bit;
    typedef Chombo4::IntVect                         ch_iv;
    typedef Chombo4::DataIterator                    ch_dit;
                                                   
    typedef ::Mat                                    petsc_mat;
    typedef Proto::Point                             pr_pt;
    typedef EBCM::MetaDataLevel<ebcm_order>          ebcm_meta;
    typedef EBCM::EBCM_Graph<ebcm_order>             ebcm_graph;
    typedef EBCM::HostLevelData<int, 1, ebcm_order > ebcm_data;
    typedef EBCM::SubVolumeVector<      ebcm_order > ebcm_subvol_vec;
    typedef EBCM::Algorithm_Framework<  ebcm_order > ebcm_framework;
    typedef Chombo4::GeometryService<   ebcm_order > ch_geoserv;

    static int dummy_func()
    {
      return 0;
    }
    struct local_stencil_t
    {
      pr_pt                           m_start;
      vector< pair< pr_pt, double > > m_stencil;
    };
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
            for(int isten = 0; isten < local->size(); isten++)
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
                              const bool           & a_mergeSmallCells,
                              const bool           & a_print)
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
                                  a_mergeSmallCells, a_print));
      return metaDataPtr;
    }

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
      bool print = false;
      shared_ptr< ebcm_meta  > ebcm_fine = makeGeometryAtSpecificNX(nx_fine, dx_fine, domain_fine, maxGrid, mergeSmallCells, print);
      shared_ptr< ebcm_meta  > ebcm_coar = makeGeometryAtSpecificNX(nx_coar, dx_coar, domain_coar, maxGrid, mergeSmallCells, print);

      /**
      typedef DummyDerivedOp   ebcm_petsc_op;
      ebcm_petsc_op   testOperator(ebcm, num_ghost, false);
      shared_ptr<ebcm_data> phi_exac_ptr(new ebcm_data(m_meta_ptr, num_ghost));
      shared_ptr<ebcm_data> lph_exac_ptr(new ebcm_data(m_meta_ptr, num_ghost));
      shared_ptr<ebcm_data> lph_calc_ptr(new ebcm_data(m_meta_ptr, num_ghost));
      **/
#endif
    }
  }; //end petsc_framwork class
}
#endif

#endif 
