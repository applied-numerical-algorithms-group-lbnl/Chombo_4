#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include <cstdio>
#include <fstream>
#include "Chombo_parstream.H"
#include "Chombo_SPMD.H"
#include "Chombo_BaseNamespaceHeader.H"

//
// pout() returns the standard output stream for Chombo codes.
//
// In parallel, each proc sends the stream to a separate file,
//  distinguished by having the proc num in the filename.
// In serial, the stream goes to std::cout.
//
// Since the pout() function gets used a lot, it should be as
// streamlined as possible, hence the use of shared variables.
//
// The other feature is that the application can choose the base filename
// for the parallel output files (the proc num is always appended).  The
// application is allowed to change the base name after pout() has been called
// the first time, which closes the previously open file, if any, and opens a
// new one.  The output files are _always_ truncated on open -- append is not
// supported.
//
// Special notes:
//  If running in parallel, and pout() is called before MPI_Initialize(),
//   it will return std::cout.
//  If poutFilename() is called before MPI_Initialize(), it fails because
//   the filename cannot be computed before the proc number is known.
//
// Method:
//  Serial is easy; pout() always returns std::cout and poutFilename()
//   always returns "cout".
//  Parallel is hard because the app can call setPoutBaseName() multiple
//   times.  Since nothing matters except when

////////////////////////////////////////////////////////////////

// the shared variables

static std::string       s_chombo_pout_filename ;
static std::string       s_chombo_pout_basename ;
static std::ofstream     s_chombo_pout ;

static bool              s_chombo_pout_init = false ;
static bool              s_chombo_pout_open = false ;

////////////////////////////////////////////////////////////////

// private functions: setFileName(), openFile().

#ifdef CH_MPI
// in parallel, compute the filename give the basename
//[NOTE: dont call this before MPI is initialized.]
static void setFileName()
{
  int outInterv = 1;
#ifdef CH_MPI
  char* charInterv = getenv("CH_OUTPUT_INTERVAL");
  if (charInterv != NULL)
    {
      outInterv =  atoi(charInterv);
      // If zero specified, change it to numProc() which should give pout.0 only
      if (outInterv == 0) outInterv=CH4_SPMD::numProc();
    }
#endif
  int thisProc = CH4_SPMD::procID();
  if ((thisProc % outInterv) != 0)
    {
      s_chombo_pout_filename = std::string("/dev/null");
    }
  else
    {
      static const size_t ProcnumSize = 1 + 10 + 1 ;  //'.' + 10digits + '\0'
      char procnum[ProcnumSize] ;
      snprintf( procnum ,ProcnumSize ,".%d" ,CH4_SPMD::procID() );
      s_chombo_pout_filename = s_chombo_pout_basename + procnum ;
    }
}

// in parallel, close the file if nec., open it and check for success
static void openFile()
{
  if ( s_chombo_pout_open )
  {
    s_chombo_pout.close();
  }
  s_chombo_pout.open( s_chombo_pout_filename.c_str() );
  // if open() fails, we have problems, but it's better
  // to try again later than to make believe it succeeded
  s_chombo_pout_open = (bool)s_chombo_pout ;
}

#else
// in serial, filename is always cout
static void setFileName()
{
  s_chombo_pout_filename = "cout" ;
}

// in serial, this does absolutely nothing
static void openFile()
{
}
#endif

////////////////////////////////////////////////////////////////

/// the stream that all output except error msgs should use
/** In serial this is the standard output, in parallel it is a
 *  different file on each proc (see setPoutBaseName()).
 */

std::ostream& pout()
{
#ifdef CH_MPI
  // the common case is _open == true, which just returns s_chombo_pout
  if ( ! s_chombo_pout_open )
  {
    // the uncommon cae: the file isn't opened, MPI may not be
    // initialized, and the basename may not have been set
    int flag_i, flag_f;
    MPI_Initialized(&flag_i);
    MPI_Finalized(&flag_f);
    // app hasn't set a basename yet, so set the default
    if ( ! s_chombo_pout_init )
    {
      s_chombo_pout_basename = "pout" ;
      s_chombo_pout_init = true ;
    }
    // if MPI not initialized, we cant open the file so return cout
    if ( ! flag_i || flag_f)
    {
      return std::cout; // MPI hasn't been started yet, or has ended....
    }
    // MPI is initialized, so file must not be, so open it
    setFileName() ;
    openFile() ;
    // finally, in case the open failed, return cout
    if ( ! s_chombo_pout_open )
    {
      return std::cout ;
    }
  }
  return s_chombo_pout ;
#else
  return std::cout;
#endif
}

//----------------------------------------------------------------

/// Set the base name for the parallel output files used by pout().
/**
 * If the file has already been used and this is a different name,
 * close the current file and open a new one.
 */
//[NOTE: in serial, this is irrelevant because s_chombo_pout_basename is not used.]

void setPoutBaseName( const std::string & a_Name )
{
  bool is_chombo_new_name = a_Name != s_chombo_pout_basename ;
  s_chombo_pout_basename = a_Name ;
  if ( s_chombo_pout_init && s_chombo_pout_open && is_chombo_new_name )
  {
    // open a new file
    //[NOTE: this is safe to do now because it's already been done once.]
    setFileName() ;
    openFile() ;
  }
  s_chombo_pout_init = true ;
}

//----------------------------------------------------------------

/// return the current filename as used by pout()
/** in serial, just return the string "cout";
 *  abort if MPI is not initialized.
 */
//[NOTE: to be static-initialization-safe, only app. code should call this function.]

const std::string & poutFileName()
{
#ifdef CH_MPI
  int flag;
  MPI_Initialized(&flag);
  if (flag)
  {
    if ( ! s_chombo_pout_open )
    {
      if ( ! s_chombo_pout_init )
      {
        s_chombo_pout_basename = "pout" ;
        s_chombo_pout_init = true ;
      }
      setFileName() ;
      //[NOTE: could open the file here, but we would still need
      //       code in pout() to handle the case where the file isn't
      //       open, so there's no point duplicating the code here. <dbs>]
    }
  }
  else
  {
    // There's really nothing reasonable to do in this case, since the
    // filename cannot be computed until after MPI is initialized and the
    // proc number is known.  So treat it as a programming bug.  Since MPI
    // isn't initialized, all procs must be running this code, so all procs
    // will fail.
    std::cerr << "error: poutFileName() cannot be called before MPI_Initialize()." << std::endl ;
    exit( 111 );
  }
#else
  // in serial, set the filename to "cout" and return it
  //[NOTE: yes, this resets the filename redundantly ever time. <dbs>]
  setFileName();
#endif
  return s_chombo_pout_filename ;
}

#include "Chombo_BaseNamespaceFooter.H"
