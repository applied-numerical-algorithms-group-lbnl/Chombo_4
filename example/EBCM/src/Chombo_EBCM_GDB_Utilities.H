#ifndef __EBCM_GDB_DBGFunctionss__
#define __EBCM_GDB_DBGFunctionss__

///  Chombo_EBCM_Graph.H is where MetaDataLevel lives.  It made sense at the time.
#include "Chombo_EBCM_Graph.H"
#include "Chombo_EBCM_PETSc_Framework.H"
#include "Chombo_EBCM_ExactSolutions.H"
#include "Chombo_EBLevelBoxData.H"

///
/**
   This is a bunch of stuff that can be called to print data
   (with varying degrees of neatness) to screen.  

   Everything here goes through raw pointers of POD because that is all
   that many versions of gdb will take.   
**/
namespace EBCM
{
  typedef ch_eigen::Matrix                               eigen_mat;

  ///print matrix to screen
  inline 
  void dumpEigen(const eigen_mat*           a_data)
  {
    if(a_data != NULL)
    {
      a_data->poutAll();
    }
  }

  template<int ebcm_order>
  class GDB_Debug_Framework
  {
  public:
    typedef EBCM::HostLevelData<double,   1, ebcm_order >                  ebcm_dou_sca_data;
    typedef EBCM::HostLevelData<double, DIM, ebcm_order >                  ebcm_dou_dim_data;
    ///write scalar data to a file
    static void writeDoubleScalar(const ebcm_dou_sca_data* a_data)
    {
      if(a_data != NULL)
      {
        string filename("ebcm_double_scalar_debug_data.hdf5");
        a_data->writeToHDF5(filename);
      }
    }
    ///write vector data to a file
    static void writeDoubleVector(const ebcm_dou_sca_data* a_data)
    {
      if(a_data != NULL)
      {
        string filename("ebcm_double_vector_debug_data.hdf5");
        a_data->writeToHDF5(filename);
      }
    }

    /// this exists to get these functions into the symbol table 
    static void allTheFuncsNULL()
    {
      writeDoubleScalar(NULL);
      writeDoubleVector(NULL);
    }
  }; //end gdb debug framework
}//end class EBCM
#endif
