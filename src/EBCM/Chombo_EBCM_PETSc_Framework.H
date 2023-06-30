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
    typedef EBCM::MetaDataLevel<ebcm_order>         ebcm_meta;
    typedef ::Mat                                   petsc_mat;
    typedef Proto::Point                            pr_pt;
    typedef Chombo4::IntVect                         ch_iv;
    typedef EBCM::HostLevelData<int, 1, ebcm_order> ebcm_data;

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
      pr_pt                 m_ptghost;
      ch_iv                 m_ivghost;
    
      virtual ~Base_PETSc_Op()
      {;}
  
      Base_PETSc_Op(const shared_ptr< ebcm_meta  >               & a_meta_ptr,
                    const pr_pt                                  & a_ghost,
                    bool                                           a_print = false)
      {
        CH_TIME("EBPetscSolver::EBPetscSolver");
        m_meta_ptr     =    a_meta_ptr;
        m_ivghost      =    ProtoCh::getIntVect(a_ghost);
        m_ptghost      =    a_ghost;
    
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

      void createOperatorPetscMatrix(bool a_print = false)
      {
      }
      void createMatrixRowMap(bool a_print = false)
      {
      }
    }; //end class Base_Petsc_Op
    
    ///test function for example/EBCM/hoeb_petsc_test
    static  void run_hoeb_petsc_tests(shared_ptr< ebcm_meta > & a_meta)
    {
      pr_pt ghost = Point::Ones(4);
      typedef Base_PETSc_Op   ebcm_petsc_op;
      ebcm_petsc_op   testOperator(a_meta, ghost);

    }
  }; //end petsc_framwork class
}
#endif

#endif 
