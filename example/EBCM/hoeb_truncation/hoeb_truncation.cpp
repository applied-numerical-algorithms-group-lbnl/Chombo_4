#include <cmath>
#include <cstdio>
#include <iostream>

#include "EBProto.H"
#include "Chombo_EBLevelBoxData.H"
#include "Chombo_LevelData.H"
#include "Chombo_BaseFab.H"

#include "Chombo_ParmParse.H"
#include "Chombo_LoadBalance.H"
#include "Chombo_ProtoInterface.H"
#include "Chombo_BRMeshRefine.H"
#include "Chombo_GeometryService.H"
#include "Chombo_EBDictionary.H"
#include "Chombo_EBChombo.H"
#include "Chombo_EBCM_Graph.H"
#include "Chombo_EBCM_HostLevelData.H"
#include "Chombo_EBCM_Algorithm_Framework.H"
#include "Chombo_EBCM_PETSc_Framework.H"
#include "Chombo_EBCM_Helmholtz.H"
#include "Chombo_EigenMatrix.H"
#include "Chombo_EBCM_GDB_Utilities.H"
#include <iomanip>


/****/


int main(int a_argc, char* a_argv[])
{
#ifdef CH_USE_PETSC  
  //because of some kind of solipsistic madness, PetscInitialize calls MPI_INIT
  PetscInt ierr = PetscInitialize(&a_argc, &a_argv, "./.petscrc",PETSC_NULL); CHKERRQ(ierr);
  Chombo4::pout() << "PetscInitialized called" << std::endl;
#else  
#ifdef CH_MPI
  MPI_Init(&a_argc, &a_argv);
  Chombo4::pout() << "MPI_Init called" << std::endl;
#endif
#endif

  CH_TIMER_SETFILE("ebcm_framework_time_table");
  {
    
    if (a_argc < 2)
    {
      cerr<< " usage " << a_argv[0] << " <input_file_name> " << endl;
      return -2;
    }
    char* in_file = a_argv[1];
    //We need to declare the first  parmparse this way
    //(the reason for this is lost in time).
    ParmParse  ppdecl(a_argc-2,a_argv+2,NULL,in_file);

    int order;
    ParmParse pp("main");
    pp.get("polynomial_order", order);
    Chombo4::pout() <<  "main: running tests for polynomial order = " << order << endl;

    if(order == 6)
    {
      EBCM::PETSc_Framework<   6>::dummy_petsc_test();   
      EBCM::Elliptic_Framework<6>::run_hoeb_truncation_tests();
    }
    else if(order == 5)
    {
      EBCM::PETSc_Framework<   5>::dummy_petsc_test();   
      EBCM::Elliptic_Framework<5>::run_hoeb_truncation_tests();
    }
    else if(order == 4)
    {
      EBCM::PETSc_Framework<   4>::dummy_petsc_test();   
      EBCM::Elliptic_Framework<4>::run_hoeb_truncation_tests();
    }
    else if(order == 3)
    {
      EBCM::PETSc_Framework<   3>::dummy_petsc_test();   
      EBCM::Elliptic_Framework<3>::run_hoeb_truncation_tests();
    }
    else if(order == 2)
    {
      EBCM::PETSc_Framework<   2>::dummy_petsc_test();   
      EBCM::Elliptic_Framework<2>::run_hoeb_truncation_tests();
    }
    else if(order == 1)
    {
      EBCM::PETSc_Framework<   1>::dummy_petsc_test();   
      EBCM::Elliptic_Framework<1>::run_hoeb_truncation_tests();
    }
    else
    {
      Chombo4::pout() << "main: Doh! unknown order = " << order << endl;
      return -1;
    }
  }

  Chombo4::pout() << "main: printing time table " << endl;
  CH_TIMER_REPORT();

#ifdef CH_MPI
#ifdef CH_USE_PETSC
  Chombo4::pout() << "main: Woo hoo! (about to call PetscFinalize)" << std::endl;
  PetscFinalize();
#else  
  Chombo4::pout() << "main: Woo hoo! (about to call MPI_Finalize" << std::endl;
  MPI_Finalize();
#endif
#endif
  
  return 0;
}
