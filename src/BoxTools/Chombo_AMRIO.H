#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif


#ifndef _Chombo_AMRIO_H_
#define _Chombo_AMRIO_H_

#include "Chombo_REAL.H"

#ifndef WRAPPER
#include <string>
#include "Chombo_LevelData.H"
#include "Chombo_REAL.H"
#include <vector>
#include "Chombo_FArrayBox.H"

#include "Chombo_RealVect.H"
#include "Chombo_NamespaceHeader.H"
using std::string;
using std::istream;
#endif

#ifdef CH_USE_HDF5
///
/**
   Writes hierarchy of levels in HDF5 format.  Only available if the
   preprocessor macro HDF5 is defined at compilation.

   {\bf Arguments:}\\
 filename  :  file to output to.\\
 a_vectGrids : grids at each level.\\
 a_vectData :  data at each level.\\
 a_vectNames:  names of variables.\\
 a_domain :  domain at coarsest level.\\
 a_dx     :  grid spacing at coarsest level.\\
 a_dt     :  time step at coarsest level.\\
 a_time     :  time.\\
 a_vectRatio :  refinement ratio at all levels
 (ith entry is refinement ratio between levels i and i + 1).\\
 a_numLevels :  number of levels to output.\\
This is blocking.

*/
void
WriteAMRHierarchyHDF5(const string& filename,
                      const std::vector<DisjointBoxLayout>& a_vectGrids,
                      const std::vector<LevelData<FArrayBox>* > & a_vectData,
                      const std::vector<string>& a_vectNames,
                      const Box& a_domain,
                      const Real& a_dx,
                      const Real& a_dt,
                      const Real& a_time,
                      const std::vector<int>& a_vectRatio,
                      const int& a_numLevels);

///
/**
   Writes hierarchy of levels in HDF5 format.  Only available if the
   preprocessor macro HDF5 is defined at compilation.

   {\bf Arguments:}\\
 filename  :  file to output to.\\
 a_vectGrids : grids at each level.\\
 a_vectData :  data at each level.\\
 a_vectNames:  names of variables.\\
 a_domain :  domain at coarsest level.\\
 a_dx     :  grid spacing in each direction at coarsest level.\\
 a_dt     :  time step at coarsest level.\\
 a_time     :  time.\\
 a_vectRatio :  refinement ratio in each direction at all levels
 (ith entry is refinement ratio in each direction between levels i and i + 1).\\
 a_numLevels :  number of levels to output.\\
This is blocking.

*/
void
WriteAnisotropicAMRHierarchyHDF5(
    const string& filename,
    const std::vector<DisjointBoxLayout>& a_vectGrids,
    const std::vector<LevelData<FArrayBox>* > & a_vectData,
    const std::vector<string>& a_vectNames,
    const Box& a_domain,
    const RealVect& a_dx,
    const Real& a_dt,
    const Real& a_time,
    const std::vector<IntVect>& a_vectRatios,
    const int& a_numLevels);

///
/**
   Writes hierarchy of levels in HDF5 format.  Only available if the
   preprocessor macro HDF5 is defined at compilation.

   {\bf Arguments:}\\
handle  :  handle to output to.\\
 a_vectGrids : grids at each level.\\
 a_vectData :  data at each level.\\
 a_vectNames:  names of variables.\\
 a_domain :  domain at coarsest level.\\
 a_dx     :  grid spacing at coarsest level.\\
 a_dt     :  time step at coarsest level.\\
 a_time     :  time.\\
 a_vectRatio :  refinement ratio at all levels
 (ith entry is refinement ratio between levels i and i + 1).\\
 a_numLevels :  number of levels to output.\\
This is not blocking.

*/
void
WriteAMRHierarchyHDF5(HDF5Handle& handle,
                      const std::vector<DisjointBoxLayout>& a_vectGrids,
                      const std::vector<LevelData<FArrayBox>* > & a_vectData,
                      const std::vector<string>& a_vectNames,
                      const Box& a_domain,
                      const Real& a_dx,
                      const Real& a_dt,
                      const Real& a_time,
                      const std::vector<int>& a_vectRatio,
                      const int& a_numLevels);

//
/**
   Writes hierarchy of levels in HDF5 format.  Only available if the
   preprocessor macro HDF5 is defined at compilation.

   {\bf Arguments:}\\
handle  :  handle to output to.\\
 a_vectGrids : grids at each level.\\
 a_vectData :  data at each level.\\
 a_vectNames:  names of variables.\\
 a_domain :  domain at coarsest level.\\
 a_dx     :  grid spacing in each direction at coarsest level.\\
 a_dt     :  time step at coarsest level.\\
 a_time     :  time.\\
 a_vectRatio :  refinement ratio in each direction at all levels
 (ith entry is refinement ratio between levels i and i + 1).\\
 a_numLevels :  number of levels to output.\\
This is not blocking.

*/
void
WriteAnisotropicAMRHierarchyHDF5(
    HDF5Handle& handle,
    const std::vector<DisjointBoxLayout>& a_vectGrids,
    const std::vector<LevelData<FArrayBox>* > & a_vectData,
    const std::vector<string>& a_vectNames,
    const Box& a_domain,
    const RealVect& a_dx,
    const Real& a_dt,
    const Real& a_time,
    const std::vector<IntVect>& a_vectRatios,
    const int& a_numLevels);

///
/**
   Writes hierarchy of levels in HDF5 format.  Only available if the
   preprocessor macro HDF5 is defined at compilation.

   {\bf Arguments:}\\
 filename  :  file to output to.\\
 a_vectGrids : grids at each level.\\
 a_vectData :  data at each level.\\
 a_domain :  domain at coarsest level.\\
 a_vectRatio :  refinement ratio at all levels
 (ith entry is refinement ratio between levels i and i + 1).\\
 a_numLevels :  number of levels to output.\\
This is blocking.
*/
void
WriteAMRHierarchyHDF5(const string& filename,
                      const std::vector<DisjointBoxLayout>& a_vectGrids,
                      const std::vector<LevelData<FArrayBox>* > & a_vectData,
                      const Box& a_domain,
                      const std::vector<int>& a_vectRatio,
                      const int& a_numLevels);

///
/**
   Writes hierarchy of levels in HDF5 format.  Only available if the
   preprocessor macro HDF5 is defined at compilation.

   {\bf Arguments:}\\
 handle  :  handle to output to.\\
 a_vectGrids : grids at each level.\\
 a_vectData :  data at each level.\\
 a_domain :  domain at coarsest level.\\
 a_vectRatio :  refinement ratio at all levels
 (ith entry is refinement ratio between levels i and i + 1).\\
 a_numLevels :  number of levels to output.\\
This is not blocking.
*/
void
WriteAMRHierarchyHDF5(HDF5Handle& handle,
                      const std::vector<DisjointBoxLayout>& a_vectGrids,
                      const std::vector<LevelData<FArrayBox>* > & a_vectData,
                      const Box& a_domain,
                      const std::vector<int>& a_vectRatio,
                      const int& a_numLevels);

#ifndef WRAPPER
///
/**
   Reads hierarchy of levels in HDF5 format.  Only available if the
   preprocessor macro HDF5 is defined at compilation.  Returns 0 on success.

   {\bf Arguments:}\\
filename  :  file to input from.\\
a_vectGrids : grids at each level.\\
a_vectData :  data at each level.\\
a_vectNames:  names of variables.\\
a_domain :  domain at coarsest level.\\
a_dx     :  grid spacing at coarsest level.\\
a_dt     :  time step at coarsest level.\\
a_time     :  time.\\
a_vectRatio :  refinement ratio at all levels
(ith entry is refinement ratio between levels i and i + 1).\\
a_numLevels :  number of levels to read.\\

{\bf Returns:} \\
status code with values:\\
0: success \\
-1: bogus number of levels \\
-2: bogus number of components \\
-3: error in readlevel \\
-4: file open failed \\
This is blocking
*/
int
ReadAMRHierarchyHDF5(const string& filename,
                     std::vector<DisjointBoxLayout>& a_vectGrids,
                     std::vector<LevelData<FArrayBox>* > & a_vectData,
                     std::vector<string>& a_vectNames,
                     Box& a_domain,
                     Real& a_dx,
                     Real& a_dt,
                     Real& a_time,
                     std::vector<int>& a_refRatio,
                     int& a_numLevels);

///
/**
   Reads hierarchy of levels in HDF5 format.  Only available if the
   preprocessor macro HDF5 is defined at compilation.  Returns 0 on success.

   {\bf Arguments:}\\
handle  :  handle to input from.\\
a_vectGrids : grids at each level.\\
a_vectData :  data at each level.\\
a_vectNames:  names of variables.\\
a_domain :  domain at coarsest level.\\
a_dx     :  grid spacing at coarsest level.\\
a_dt     :  time step at coarsest level.\\
a_time     :  time.\\
a_vectRatio :  refinement ratio at all levels
(ith entry is refinement ratio between levels i and i + 1).\\
a_numLevels :  number of levels to read.\\

{\bf Returns:} \\
status code with values:\\
0: success \\
-1: bogus number of levels \\
-2: bogus number of components \\
-3: error in readlevel \\
-4: file open failed \\
This is not blocking
*/
int
ReadAMRHierarchyHDF5(HDF5Handle& handle,
                     std::vector<DisjointBoxLayout>& a_vectGrids,
                     std::vector<LevelData<FArrayBox>* > & a_vectData,
                     std::vector<string>& a_vectNames,
                     Box& a_domain,
                     Real& a_dx,
                     Real& a_dt,
                     Real& a_time,
                     std::vector<int>& a_refRatio,
                     int& a_numLevels);

///
/**
   Reads hierarchy of levels in HDF5 format.  Only available if the
   preprocessor macro HDF5 is defined at compilation.  Returns 0 on success.

   {\bf Arguments:}\\
filename  :  file to input from.\\
a_vectGrids : grids at each level.\\
a_vectData :  data at each level.\\
a_domain :  domain at coarsest level.\\
a_vectRatio :  refinement ratio at all levels
(ith entry is refinement ratio between levels i and i + 1).\\
a_numLevels :  number of levels to read.\\

{\bf Returns:} \\
status code with values:\\
0: success \\
-1: bogus number of levels \\
-2: bogus number of components \\
-3: error in readlevel \\
-4: file open failed \\
This is blocking.
*/
int
ReadAMRHierarchyHDF5(const string& filename,
                     std::vector<DisjointBoxLayout>& a_vectGrids,
                     std::vector<LevelData<FArrayBox>* > & a_vectData,
                     Box& a_domain,
                     std::vector<int>& a_refRatio,
                     int& a_numLevels);

///
/**
   Reads hierarchy of levels in HDF5 format.  Only available if the
   preprocessor macro HDF5 is defined at compilation.  Returns 0 on success.

   {\bf Arguments:}\\
handle  :  handle to input from.\\
a_vectGrids : grids at each level.\\
a_vectData :  data at each level.\\
a_domain :  domain at coarsest level.\\
a_vectRatio :  refinement ratio at all levels
(ith entry is refinement ratio between levels i and i + 1).\\
a_numLevels :  number of levels to read.\\

{\bf Returns:} \\
status code with values:\\
0: success \\
-1: bogus number of levels \\
-2: bogus number of components \\
-3: error in readlevel \\
-4: file open failed \\
This is not blocking.
*/
int
ReadAMRHierarchyHDF5(HDF5Handle& handle,
                     std::vector<DisjointBoxLayout>& a_vectGrids,
                     std::vector<LevelData<FArrayBox>* > & a_vectData,
                     Box& a_domain,
                     std::vector<int>& a_refRatio,
                     int& a_numLevels);

///
/**
   Reads hierarchy of levels in ANISOTROPIC HDF5 format.  Only available if the
   preprocessor macro HDF5 is defined at compilation.  Returns 0 on success.

   {\bf Arguments:}\\
filename  :  file to input from.\\
a_vectGrids : grids at each level.\\
a_vectData :  data at each level.\\
a_vectNames:  names of variables.\\
a_domain :  domain at coarsest level.\\
a_dx     :  grid spacing at coarsest level.\\
a_dt     :  time step at coarsest level.\\
a_time     :  time.\\
a_vectRatio :  refinement ratio at all levels
(ith entry is refinement ratio between levels i and i + 1).\\
a_numLevels :  number of levels to read.\\

{\bf Returns:} \\
status code with values:\\
0: success \\
-1: bogus number of levels \\
-2: bogus number of components \\
-3: error in readlevel \\
-4: file open failed \\
This is blocking
*/
int
ReadAnisotropicAMRHierarchyHDF5(const string& filename,
                     std::vector<DisjointBoxLayout>& a_vectGrids,
                     std::vector<LevelData<FArrayBox>* > & a_vectData,
                     std::vector<string>& a_vectNames,
                     Box& a_domain,
                     RealVect& a_dx,
                     Real& a_dt,
                     Real& a_time,
                     std::vector<IntVect>& a_refRatio,
                     int& a_numLevels);

///
/**
   Reads hierarchy of levels in ANISOTROPIC HDF5 format.  Only available if the
   preprocessor macro HDF5 is defined at compilation.  Returns 0 on success.

   {\bf Arguments:}\\
handle  :  handle to input from.\\
a_vectGrids : grids at each level.\\
a_vectData :  data at each level.\\
a_vectNames:  names of variables.\\
a_domain :  domain at coarsest level.\\
a_dx     :  grid spacing at coarsest level.\\
a_dt     :  time step at coarsest level.\\
a_time     :  time.\\
a_vectRatio :  refinement ratio at all levels
(ith entry is refinement ratio between levels i and i + 1).\\
a_numLevels :  number of levels to read.\\

{\bf Returns:} \\
status code with values:\\
0: success \\
-1: bogus number of levels \\
-2: bogus number of components \\
-3: error in readlevel \\
-4: file open failed \\
This is not blocking
*/
int
ReadAnisotropicAMRHierarchyHDF5(HDF5Handle& handle,
                     std::vector<DisjointBoxLayout>& a_vectGrids,
                     std::vector<LevelData<FArrayBox>* > & a_vectData,
                     std::vector<string>& a_vectNames,
                     Box& a_domain,
                     RealVect& a_dx,
                     Real& a_dt,
                     Real& a_time,
                     std::vector<IntVect>& a_refRatio,
                     int& a_numLevels);

///
/**
   Reads hierarchy of levels in ANISOTROPIC HDF5 format.  Only available if the
   preprocessor macro HDF5 is defined at compilation.  Returns 0 on success.

   {\bf Arguments:}\\
filename  :  file to input from.\\
a_vectGrids : grids at each level.\\
a_vectData :  data at each level.\\
a_domain :  domain at coarsest level.\\
a_vectRatio :  refinement ratio at all levels
(ith entry is refinement ratio between levels i and i + 1).\\
a_numLevels :  number of levels to read.\\

{\bf Returns:} \\
status code with values:\\
0: success \\
-1: bogus number of levels \\
-2: bogus number of components \\
-3: error in readlevel \\
-4: file open failed \\
This is blocking.
*/
int
ReadAnisotropicAMRHierarchyHDF5(const string& filename,
                     std::vector<DisjointBoxLayout>& a_vectGrids,
                     std::vector<LevelData<FArrayBox>* > & a_vectData,
                     Box& a_domain,
                     std::vector<IntVect>& a_refRatio,
                     int& a_numLevels);

///
/**
   Reads hierarchy of levels in HDF5 format.  Only available if the
   preprocessor macro HDF5 is defined at compilation.  Returns 0 on success.

   {\bf Arguments:}\\
handle  :  handle to input from.\\
a_vectGrids : grids at each level.\\
a_vectData :  data at each level.\\
a_domain :  domain at coarsest level.\\
a_vectRatio :  refinement ratio at all levels
(ith entry is refinement ratio between levels i and i + 1).\\
a_numLevels :  number of levels to read.\\

{\bf Returns:} \\
status code with values:\\
0: success \\
-1: bogus number of levels \\
-2: bogus number of components \\
-3: error in readlevel \\
-4: file open failed \\
This is not blocking.
*/
int
ReadAnisotropicAMRHierarchyHDF5(HDF5Handle& handle,
                     std::vector<DisjointBoxLayout>& a_vectGrids,
                     std::vector<LevelData<FArrayBox>* > & a_vectData,
                     Box& a_domain,
                     std::vector<IntVect>& a_refRatio,
                     int& a_numLevels);

#endif


///
/** Writes a plotfile using the same format as WriteAMRHierarchyHDF5, but
    for a FArrayBox.  This is useful for debugging.  *a_dataPtr is written
    to a file named fab.hdf5
*/
void writeFAB(const FArrayBox* a_dataPtr);



///
/** Writes a plotfile using the same format as WriteAMRHierarchyHDF5, but
    for a BaseFab<Real>.  This is useful for debugging.  *a_dataPtr is written
    to a file named fab.hdf5
*/
void
writeBFR(const BaseFab<Real>* a_dataPtr);


///
/** Writes a plotfile using the same format as WriteAMRHierarchyHDF5, but
    for a FArrayBox.  This is useful for debugging.  *a_dataPtr is written
    to the file given by a_filename.  If a non-empty std::vector of strings is
    given, they are used as the component names in the HDF5 header; else
    the names are set to "component_0", "component_1", ....
*/
void
writeFABname(const FArrayBox      * a_dataPtr,
             const char           * a_filename,
             const std::vector<string> & a_compNames = std::vector<string>(),
             const Real           & a_dx = 1.0);


///
/** Write a plotfile using the same format as WriteAMRHierarchyHDF5, but
    for a single LevelData<FArrayBox>.  Useful for debugging.  *a_dataPtr is
    written to a file named LDF.hdf5.
*/
void
writeLevel(const LevelData<FArrayBox>* a_dataPtr);


///
/** Write a plotfile using the same format as WriteAMRHierarchyHDF5, but
    for a single LevelData<FArrayBox>. Useful for debugging.  *a_dataPtr is
    written to the file given by a_filename.
*/
void
writeLevelname(const LevelData<FArrayBox>* a_dataPtr,
               const char*                 a_filename);

///
/** Write a plotfile using the same format as WriteAMRHierarchyHDF5, but
    for a std::vector<LevelData<FArrayBox>*>. Useful for debugging.  *a_dataPtr is
    written to the file given by a_filename.
*/
void
writestd::vectorLevelName(const std::vector<LevelData<FArrayBox>*>* a_dataPtr,
                     const std::vector<int>*                   a_refRatios,
                     const char*                          a_filename);

///
/** Write a plotfile using the same format as WriteAMRHierarchyHDF5, but
    for a single DisjointBoxLayout (do this by creating a LevelData<FArrayBox>
    and then calling writeLevel.  Useful for debugging.  *a_dataPtr is
    written to a file named DBL.hdf5.
*/
void
writeDBL(const DisjointBoxLayout* a_dataPtr);


void writeCopier(const Copier* a_copier);


///
/** Write a plotfile using the same format as WriteAMRHierarchyHDF5, but
    for a single DisjointBoxLayout (do this by creating a LevelData<FArrayBox>
    and then calling writeLevel.  Useful for debugging.  *a_dataPtr is
    written to the file given by a_filename.
*/
void
writeDBLname(const DisjointBoxLayout* a_dataPtr,
             const char*              a_filename);

///
/** Write a plotfile with levels that don't necessarily begin with
0 (will translate things so that lBase comes out as level 0, etc)
Levels in std::vector arguments should correspond with actual levels.
(and so should be defined for a_levels.begin() through a_levels.end().
   {\bf Arguments:}\\
 filename  :  file to output to.\\
 a_vectGrids : grids at each level.\\
 a_vectData :  data at each level.\\
 a_vectNames:  names of variables.\\
 a_domain :  domain at base level (given by  a_levels.begin()).\\
 a_dx     :  grid spacing at base level.\\
 a_dt     :  time step at base level.\\
 a_time     :  time.\\
 a_vectRatio :  refinement ratio at all levels
 (ith entry is refinement ratio between levels i and i + 1).\\
 a_numLevels :  number of levels to output.\\
This is blocking.

*/
void
WritePartialAMRHierarchyHDF5(const string& filename,
                             const std::vector<DisjointBoxLayout>& a_vectGrids,
                             const std::vector<LevelData<FArrayBox>* > & a_vectData,
                             const std::vector<string>& a_vectNames,
                             const Box& a_baseDomain,
                             const Real& a_baseDx,
                             const Real& a_dt,
                             const Real& a_time,
                             const std::vector<int>& a_vectRatio,
                             const Interval& a_levels);

#endif // CH_USE_HDF5

#include "Chombo_NamespaceFooter.H"

#endif
