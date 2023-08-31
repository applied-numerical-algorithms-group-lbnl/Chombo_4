#ifndef _PROTO_WRITERECTMDARRAY_H_
#define _PROTO_WRITERECTMDARRAY_H_
#include "Proto_BoxData.H"
#include "Proto_VisitWriter.H"

#include <cstring>
#include <cstdio>

namespace Proto 
{

/** @defgroup boxdata_writer BoxData Writer */
/** @addtogroup boxdata_writer */
/*@{*/

/// Write Scalar BoxData
/**
    Write's the data in a BoxData<double> to a default named .vtk file.

    \param a_array    The BoxData to be printed.
    \param dx         (Optional) Physical spacing between points in a_array. (default: 1)
*/
inline const char* WriteBoxData(BoxData<double>& a_array,
                    double dx = 1.0);

/// Write Scalar BoxData with Name
/**
    Write's the data in a BoxData<double> to a given file name. The given filename should include the relevent extension. (e.g. name.vtk)

    \param a_filename Name of output file
    \param a_array    The BoxData to be printed.
    \param dx         (Optional) Physical spacing between points in a_array. (default: 1)
*/
inline void WriteBoxData(const char*                a_filename,
                   BoxData<double>& a_array,
                   double dx = 1.0);


/// Write Non-Scalar BoxData with Name
/**
    Write's the data in a BoxData<T,C,D,E> to a given file name. The given filename should include the relevent extension. (e.g. name.vtk)

    \param a_filename Name of output file
    \param a_array    The BoxData to be printed.
    \param a_varnames Array of strings corresponding to the variable names of a_array.
    \param a_corner   DIM long array of doubles defining the desired origin
    \param a_h        Physical spacing between points
*/
  template<class T, unsigned int C, MemType MEMTYPE, unsigned char D, unsigned char E>
inline void WriteBoxData(const char*           a_filename,
             BoxData<T, C,MEMTYPE, D, E>&       a_array,
             const char* const          a_varnames[C*D*E],
             const double*              a_corner,
             const double&              a_h);

/// Write Non-Scalar BoxData to Existing File
/**
    Write's the data in a BoxData<T,C,D,E> to an existing file.

    \param a_filename Name of output file
    \param a_array    The BoxData to be printed.
    \param a_varnames Array of strings corresponding to the variable names of a_array.
    \param a_corner   DIM long array of doubles defining the desired origin
    \param a_h        Physical spacing between points
*/
  template<class T, unsigned int C,MemType MEMTYPE, unsigned char D, unsigned char E>
inline void WriteBoxData(FILE*                 a_fp,
                         BoxData<T,C,MEMTYPE,D,E>&          a_array,
             const char* const          a_varnames[C*D*E],
             const double*              a_corner,
             const double&              a_h);

// End of Doxygen documentation module
/*@}*/
using std::string;
//convenience function, turn off iter in filename by setting a_iter < 0
template <unsigned int n>
inline void WriteData(BoxData<double,n> &    a_state,
               const   int       &    a_iter, 
               const   double    &    a_dx,
               const   string    &    a_varnamebase,
               const   string    &    a_filenamebase)
{
  double origin[DIM];
  for (int ii = 0; ii < DIM; ii++)
  {
    origin[ii] = 0.0;
  }

  char basename[1024];
  if(a_iter >= 0) 
  {
    sprintf(basename,"%s.%06d",a_filenamebase.c_str(),a_iter);
  }
  else
  {
    sprintf(basename,"%s",a_filenamebase.c_str());
  }
    
  char* velnames[n];
  for(int idir = 0; idir < n; idir++)
  {
    velnames[idir] = new char[1024];
    sprintf(velnames[idir],"%s%06d",a_varnamebase.c_str(),idir);
  }

  WriteBoxData(basename,a_state,velnames,origin,a_dx);

  for(int idir = 0; idir < n; idir++)
  {
    delete[] velnames[idir];
  }

}  

template <int n> inline void viewData(BoxData<double, n>* a_statePtr)
{
  WriteData<n>(*a_statePtr, -1, 1.0, string("var"), string("tempdata"));
  system("visit -o tempdata.vtk");
}

#include "implem/Proto_WriteBoxDataImplem.H"

} // end Proto namespace
#endif 
