#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include <iostream>
#include <cstring>
// #extern "C" {      // The #extern "C" might have been here for a reason...
#include <unistd.h>
// }

#include "Chombo_SPMD.H"
#include "Chombo_parstream.H"
namespace CH4_SPMD
{



using std::endl;

// try a 30 Mbyte max message size and see if that helps.

unsigned long long CH_MAX_MPI_MESSAGE_SIZE = 30*1024*1024;
unsigned long long CH_MaxMPISendSize = 0;
unsigned long long CH_MaxMPIRecvSize  = 0;

int reportMPIStats()
{
  Proto::pout()<<"Chombo message size limit:"<< CH_MAX_MPI_MESSAGE_SIZE<<"\n"
        <<"Max send message size:"<<CH_MaxMPISendSize<<"\n"
        <<"Max recv message size:"<<CH_MaxMPIRecvSize<<std::endl;
  return 0;
}

std::vector<int> pids;

#ifndef CH_MPI

int procID()
{
  return 0;
}

// reset this to fool serial code into thinking its parallel
int num_procs = 1 ;

int GetRank(int pid)
{
  return 0;
}

int GetPID(int rank)
{
  CH_assert(rank == 0);
  return getpid();
}

unsigned int numProc()
{
  return num_procs;
}

#else // CH_MPI version

//int GetPID(int rank)
//{
//  if (pids.size() == 0)
//    {
//      int proc = getpid();
//      gather(pids, proc, uniqueProc(SerialTask::compute));
//      broadcast(pids, uniqueProc(SerialTask::compute));
//    }
//  if (rank<0) return -1;
//  if (rank>=pids.size()) return -2;
//  return pids[rank];
//}

int GetRank(int pid)
{
  for (int i=0; i<pids.size(); ++i)
    {
      if (pids[i]== pid) return i;
    }
  return -1;
}

// this is a global variable which indicates whether
// or not an MPI command should be used to deduce rank.
// needed for applications which switch communicators.
// set g_resetProcID=true to force next procID() call to 
// query MPI_Comm_rank
bool g_resetProcID;

int procID()
{
  static bool firstCall = true;
  static int lastProcID = 0;
  if (firstCall || g_resetProcID )
  {
    g_resetProcID = false;
    firstCall = false;

    MPI_Comm_rank(Chombo_MPI::comm, &lastProcID);
  }
  return lastProcID;
}

unsigned int numProc()
{
  static int ret = -1;
  if (ret == -1)
  {
    MPI_Comm_size(Chombo_MPI::comm, &ret);
  }
  return ret;
}

// hopefully static copy of opaque handles
MPI_Comm Chombo_MPI::comm = MPI_COMM_WORLD;

#endif // CH_MPI


// return id of unique processor for special serial tasks
int
uniqueProc(const SerialTask::task& a_task)
{
#ifdef NDEBUG
    switch (a_task)
    {
    case SerialTask::compute:
    default:
        return (0);
        //break;  // unreachable break can generate compiler warning
    }
#else
// in mpi, the debugger attaches to process 0
    return (0);
#endif
}

}
