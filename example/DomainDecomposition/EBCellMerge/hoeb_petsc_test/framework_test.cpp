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
#include "EBMultigrid.H"

#include "DebugFunctions.H"
#include "Hoeb_ExactSolutions.H"

#include "Chombo_LAPACKMatrix.H"
#include <iomanip>


typedef Chombo4::Box                                   ch_box;
typedef Chombo4::DisjointBoxLayout                     ch_dbl;
typedef Chombo4::ProblemDomain                         ch_probdom;
typedef Chombo4::DataIterator                          ch_dit;
typedef Chombo4::IntVect                               ch_iv;
typedef Chombo4::MayDay                                ch_mayday;
//typedef Chombo4::LAPACKMatrix                          ch_mat;
typedef Proto::RealVect                                pr_rv;
typedef Proto::IndexTM<Real, DIM>                      pr_itm_r_dim;
typedef Proto::IndexTM<int , DIM>                      pr_itm_i_dim;
typedef Proto::IndexTM<Real, DIM-1>                    pr_itm_r_dmo;
typedef Proto::IndexTM<int , DIM-1>                    pr_itm_i_dmo;
typedef Proto::BaseIF                                  pr_baseif;
typedef ch_eigen::Matrix                               eigen_mat;

/****/

int main(int a_argc, char* a_argv[])
{
#ifdef CH_MPI
  MPI_Init(&a_argc, &a_argv);
  Chombo4::pout() << "MPI INIT called" << std::endl;
#endif

  using Chombo4::pout;
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
    Chombo4::pout() <<  "Running test for polynomial order = " << order << endl;
    if(order == 6)
    {
      EBCM::Algorithm_Framework<6>::runTests();
    }
    else if(order == 5)
    {
      EBCM::Algorithm_Framework<5>::runTests();
    }
    if(order == 4)
    {
      EBCM::Algorithm_Framework<4>::runTests();
    }
    else if(order == 3)
    {
      EBCM::Algorithm_Framework<3>::runTests();
    }
    else if(order == 2)
    {
      EBCM::Algorithm_Framework<2>::runTests();
    }
    else if(order == 1)
    {
      EBCM::Algorithm_Framework<1>::runTests();
    }
    else
    {
      Chombo4::pout() << "main: unknown order = " << order << endl;
      return -1;
    }
  }

  Chombo4::pout() << "printing time table " << endl;
  CH_TIMER_REPORT();
#ifdef CH_MPI
  Chombo4::pout() << "about to call MPI Finalize" << std::endl;
  MPI_Finalize();
#endif
  return 0;
}
