#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

//
// These are the inevitable functions that people
// can't avoid using when writing a SPMD program.
// It is a minimal set from which more elaborate
// functionality can be generated.  As always, a
// user is free to utilize the entire MPI programming
// on their own platform.  The functions are
// assured to work on all machines supported.
//

#ifndef _Chombo_SPMD_H_
#define _Chombo_SPMD_H_

#include <string>
#include "Chombo_REAL.H"
#include <vector>
#include "Chombo_MayDay.H"


#ifdef CH_MPI
#include "mpi.h"
#endif

using namespace Chombo4;
namespace CH4_SPMD
{
#ifdef CH_MPI
  struct Chombo_MPI
  {
    static MPI_Comm comm;
  };

#else
// this can be changed for debugging parallel code in serial
  extern int num_procs ;
#endif

  extern unsigned long long CH_MAX_MPI_MESSAGE_SIZE;

  extern unsigned long long CH_MaxMPISendSize;
  extern unsigned long long CH_MaxMPIRecvSize;

  int reportMPIStats();

/// local process ID
/**
   Returns the ID of the locally running process in the range 0 <=
   procID() < numProc().  This has no relation to the operating system
   pid.  There is always a procID() == 0.  */
  int procID();

  inline int CHprocID()
  {
    return procID();
  }

/// number of parallel processes
/**
   Returns the number of parallel processes running.
   Always returns at least 1. */
  unsigned int numProc();

///all ranks wait here to sync-up
/**
   All MPI ranks wait here to sync-up.  Calls MPI_Barrier(comm)  */
  void barrier(void);

  template <class T>
  size_t linearSize(const T& inputT)
  {
    return inputT.linearSize();
  }

  template <class T>
  void linearIn(T& a_outputT, const void* const inBuf)
  {
    a_outputT.linearIn(inBuf);
  }

  template <class T>
  void linearOut(void* const a_outBuf, const T& inputT)
  {
    inputT.linearOut(a_outBuf);
  }

  /// int specialization
  template < > inline
  size_t linearSize(const int& inputT)
  {
    return sizeof(int);
  }

  template < > inline
  void linearIn(int& a_outputT, const void* const a_inBuf)
  {
    const int* intbuf = (const int *) a_inBuf;
    a_outputT = *intbuf;
  }

  template < > inline
  void linearOut(void* const a_outBuf, const int& a_input)
  {
    int* intbuf = (int*)a_outBuf;
    *intbuf = a_input;
  }
  
//std::vector<T> specialization of linearSize
  template <class T>
  size_t linearListSize(const std::vector<T>& a_input, bool a_printStuff);

//std::vector<T> specialization of linearIn
  template <class T>
  void linearListIn(std::vector<T>& a_outputT,
                    size_t a_expectedSize,
                    const void* const a_inBuf,
                    bool a_printStuff);

//std::vector<T> specialization of linearOut
  template <class T>
  void linearListOut(void* const a_outBuf,
                     size_t a_expectedSize,
                     const std::vector<T>& a_inputT,
                     bool a_printStuff);
  
  /// vector<int> specialization
  template < > inline 
  size_t linearSize(const std::vector<int>& inputT)
  {
    return linearListSize(inputT, false);
  }

  template < > inline
  void linearIn(std::vector<int>& a_outputT, const void* const a_inBuf)
  {
    size_t* vecsizeptr = (size_t *) (a_inBuf);
    size_t vecsize = *vecsizeptr;
    size_t expected = (vecsize+1)*sizeof(size_t) + vecsize*sizeof(int);


    return linearListIn(a_outputT, expected, a_inBuf, true);
  }

  template < > inline
  void linearOut(void* const a_outBuf, const std::vector<int>& a_input)
  {
    size_t vecsize  = a_input.size();
    size_t expected = vecsize*sizeof(int) + (vecsize+1)*sizeof(size_t);

    return linearListOut(a_outBuf, expected, a_input, true);
  }
#ifdef CH_MPI
///gather from multiple processes
/**
   Gathers a_input of every process onto std::vector<T> a_outVec on process
   number a_dest.  It is an error if a_dest < 0 or a_dest >= numProc().
   a_outVec is numProc() long with the value of a_input on every
   process in its elements.\\

   T must have have the following functions:
   \\

   int linearSize<T>(const T& inputT); \\
   void linearIn<T>(T& a_outputT, const void* const inBuf); \\
   void linearOut<T>(void* a_outBuf, const T& inputT); \\

*/
  template <class T>
  void gather(std::vector<T>& a_outVec, const T& a_input, int a_dest);

#endif

///broadcast to every process
/**
   Broadcasts a_inAndOut to every process from process number a_src.
   It is an error if a_src < 0 or a_src >= numProc().  \\

   T must have have the following functions: \\

   int linearSize<T>(const T& inputT); \\
   void linearIn<T>(T& a_outputT, const void* const inBuf); \\
   void linearOut<T>(void* a_outBuf, const T& inputT); \\ */
  template <class T>
  void broadcast(T& a_inAndOut,  int a_src);




  class SerialTask
  {
  public:
    enum task
    {
      compute=0
    };
  };

  int GetPID(int rank);
  int GetRank(int pid);

// return id of unique processor for special serial tasks
  int
  uniqueProc(const SerialTask::task& a_task);

}

#include "Chombo_SPMDI.H"

#endif
