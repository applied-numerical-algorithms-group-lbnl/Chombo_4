#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif
#ifndef _PROTO_SPMD_H_
#define _PROTO_SPMD_H_
#include <array>
#include "Proto.H"

#ifdef PR_MPI
#include "mpi.h"
#endif

using namespace std;
namespace Proto
{
#ifdef PR_MPI
    template <typename T=void>
    struct Proto_MPI
    {
        static MPI_Comm comm;
    };
    template <typename T> MPI_Comm Proto_MPI<T>::comm = MPI_COMM_WORLD;
#else
    // this can be changed for debugging parallel code in serial
    extern int num_procs;
#endif
#define PR_MAX_MPI_MESSAGE_SIZE 30*1024*1024  
    
    /// Get Local Process ID
    /**
      Returns the ID of the locally running process in the range 0 <=
      procID() < numProc().  This has no relation to the operating system
      pid.  There is always a procID() == 0.  */
    inline int procID()
    {
#ifdef PR_MPI
        int retval;
        MPI_Comm_rank(Proto_MPI<void>::comm, &retval);
        return retval;
#else
        return 0;
#endif
    }
    
    inline int PRprocID()
    {
        return procID();
    }

    /// Get Number of Ranks
    /**
      Returns the number of parallel processes running.
      Always returns at least 1. */
    inline unsigned int numProc()
    {
#ifdef PR_MPI
        static int ret = -1;
        if (ret == -1)
        {
            MPI_Comm_size(Proto_MPI<void>::comm, &ret);
        }
        return ret;
#else
        return 1;
#endif
    }
    /// Parallel Barrier
    /**
      All MPI ranks wait here to sync-up.  Calls MPI_Barrier(comm).
      This is a no-op in the non-MPI/serial case.
    */
    inline void barrier()
    {
#ifdef PR_MPI
        MPI_Barrier(Proto_MPI<>::comm);
#endif
    }

#ifdef PR_MPI
    template<typename T>
    inline MPI_Datatype mpiDatatype()
    {
        // TODO: This is an error, but the headers wont let me use MayDay here. 
        return MPI_DATATYPE_NULL;
    }
    
    template<>
    inline MPI_Datatype mpiDatatype<int>(){return MPI_INT;}
    
    template<>
    inline MPI_Datatype mpiDatatype<double>(){return MPI_DOUBLE;}
#endif

}
#endif
