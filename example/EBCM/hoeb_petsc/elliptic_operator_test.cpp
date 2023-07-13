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
#include "Chombo_EigenMatrix.H"

#include <iomanip>


/****/


int main(int a_argc, char* a_argv[])
{
#ifdef CH_MPI
  MPI_Init(&a_argc, &a_argv);
  Chombo4::pout() << "MPI INIT called" << std::endl;
#endif

  using Chombo4::pout;
  CH_TIMER_SETFILE("ebcm_framework.time_table");
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
    Chombo4::pout() <<  "Running hoeb_petsc_test for polynomial order = " << order << endl;

    if(order == 6)
    {
      shared_ptr< EBCM::MetaDataLevel< 6 > >             ebcm_meta;      
      EBCM::Algorithm_Framework<6>::run_hoeb_petsc_tests(ebcm_meta);
    }
    else if(order == 5)
    {
      shared_ptr< EBCM::MetaDataLevel< 5 > >             ebcm_meta;      
      EBCM::Algorithm_Framework<5>::run_hoeb_petsc_tests(ebcm_meta);
    }
    if(order == 4)
    {
      shared_ptr< EBCM::MetaDataLevel< 4 > >             ebcm_meta;      
      EBCM::Algorithm_Framework<4>::run_hoeb_petsc_tests(ebcm_meta);
    }
    else if(order == 3)
    {
      shared_ptr< EBCM::MetaDataLevel< 3 > >             ebcm_meta;      
      EBCM::Algorithm_Framework<3>::run_hoeb_petsc_tests(ebcm_meta);
    }
    else if(order == 2)
    {
      shared_ptr< EBCM::MetaDataLevel< 2 > >             ebcm_meta;      
      EBCM::Algorithm_Framework<2>::run_hoeb_petsc_tests(ebcm_meta);
    }
    else if(order == 1)
    {
      shared_ptr< EBCM::MetaDataLevel< 1 > >             ebcm_meta;      
      EBCM::Algorithm_Framework<1>::run_hoeb_petsc_tests(ebcm_meta);
    }
    else
    {
      Chombo4::pout() << "main: Doh! unknown order = " << order << endl;
      return -1;
    }
  }

  Chombo4::pout() << "Printing time table " << endl;
  CH_TIMER_REPORT();
#ifdef CH_MPI
  Chombo4::pout() << "Woo hoo! \n About to call MPI Finalize" << std::endl;
  MPI_Finalize();
#endif
  return 0;
}
