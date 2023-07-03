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
    typedef EBCM::MetaDataLevel<ebcm_order>         ebcm_meta;
    typedef ::Mat                                   petsc_mat;
    typedef Proto::Point                            pr_pt;
    typedef Chombo4::IntVect                        ch_iv;
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
    
        createMatrixRowMap(a_print);

        if(a_print)
        {
          Chombo4::pout() << "Base_Petsc_Op:: creating the matrix" << endl;
        }

        createOperatorPetscMatrix(a_print);
      }

    private:
      /// banning weak construction but assignment and copy constructors should work fine.
      Base_PETSc_Op();
      
      ///makes m_row_map_ptr (operator uses this to map between volumes and petsc matrix rows)
      void createMatrixRowMap(bool a_print = false)
      {
        m_row_map_ptr = shared_ptr<ebcm_data>(new ebcm_data(m_meta_ptr, m_ghost));
        const auto& graphs = *(m_meta_ptr->m_graphs);
        const auto& grids  =  (m_meta_ptr->m_grids);
        auto dit = grids.dataIterator();

        vector<shared_ptr<ebcm_subvol_vec> > all_valid_vols(dit.size());
        
        //get the number of points on each proc.
        int numPtsThisProc = 0;
        for(int ibox =0; ibox < dit.size(); ibox++)
        {
          const auto& valid =  grids[dit[ibox]];
          const auto& graph = graphs[dit[ibox]];
          //valid is correct here.  ghost is handled via exchange.
          all_valid_vols[ibox]= shared_ptr<ebcm_subvol_vec>(new ebcm_subvol_vec(graph, valid, a_print));
          numPtsThisProc += all_valid_vols[ibox]->size();
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
        auto& mapdata = *(m_row_map_ptr->m_data);
        int curgid = startgid;
        for(int ibox = 0; ibox < dit.size(); ibox++)
        {
          const auto& volvec = *(all_valid_vols[ibox]);
          auto      & mapfab =      mapdata[dit[ibox]];
          for(int ivec = 0; ivec < volvec.size(); ivec++)
          {
            const auto& volume = volvec[ivec];
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
      void createOperatorPetscMatrix(bool a_print = false)
      {
      } //end function createOperatorPetscMatrix
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
