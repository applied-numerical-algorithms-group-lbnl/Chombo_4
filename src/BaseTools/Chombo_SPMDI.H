#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _Chombo_SPMDI_H_
#define _Chombo_SPMDI_H_


#include "Chombo_parstream.H"
#include "Chombo_CH_assert.H"

namespace CH4_SPMD
{

// default implementations for linearization routines.




#ifdef CH_MPI

extern void AttachDebugger(int);
/*****************************/
//gather a_input into a a_outVec
/*****************************/
template <class T>
inline void
gather(Vector<T>& a_outVec, const T& a_input, int a_dest)
{
  CH_assert (a_dest >= 0);
  CH_assert(a_dest <  numProc());
  //now THIS size lives on THIS processor
  int isize = linearSize(a_input);

  //make stuff for linearout
  void* loclBuf = malloc(isize);
  if (loclBuf == NULL)
    Chombo4::MayDay::Error("out of memory in gather 1");

  //put linearized T into its proper buffer
  linearOut(loclBuf, a_input);

  int nProcess = numProc();
  int sendCount = 1;
  int recdCount = 1;

  //need to gather isizes onto processor a_dest
  int* vectSize = NULL;
  int* vectDisp = NULL;
  void* sendBuf = static_cast<void*>(&isize);
  //allocate received buffer
  if (procID() == a_dest)
    {
      vectSize = new int[nProcess];
      vectDisp = new int[nProcess];
    }

  int result1 = MPI_Gather(sendBuf, sendCount, MPI_INT,
                           vectSize,recdCount, MPI_INT,
                           a_dest,  Chombo_MPI::comm);

  if (result1 != MPI_SUCCESS)
    Chombo4::MayDay::Error("Gather<T> failed in MPI_Gather 1");

  //make memory for gather, linearin
  void* recdBuf = NULL;
  if (procID() == a_dest)
    {
      size_t itotsize=0;
      for (int iproc = 0; iproc < nProcess; iproc++)
        {
          vectDisp[iproc] = itotsize;
          itotsize += vectSize[iproc];
        }
      recdBuf = malloc(itotsize);
      if (recdBuf == NULL)
        {
          Chombo4::MayDay::Error("out of memory in gather 2");
        }
    }

  //gather data
  int result2 = MPI_Gatherv(loclBuf, isize, MPI_BYTE,
                            recdBuf, vectSize, vectDisp, MPI_BYTE,
                            a_dest, Chombo_MPI::comm);
  if (result2 != MPI_SUCCESS)
    Chombo4::MayDay::Error("Gather<T> failed in MPI_Gather 2");

  if (procID() == a_dest)
    {
      //calculate offset into array for current processor
      int ioffset = 0;
      a_outVec.resize(nProcess);
      //need to cast to char* to do pointer arithmetic
      char* arithPtr = (char*)recdBuf;
      for (int iproc = 0; iproc < nProcess; iproc++)
        {
          ioffset = vectDisp[iproc];
          char* thisProcBuf = arithPtr + ioffset;
          linearIn(a_outVec[iproc], thisProcBuf);
        }

      //delete memory for dest-specific arrays
      delete[] vectSize;
      delete[] vectDisp;
      free(recdBuf);
    }

  //delete memory for local buffer
  free(loclBuf);
}

/*****************************/
//broadcast T everywhere
/*****************************/
template <class T>
inline void
broadcast(T& a_inAndOut,  int a_src)
{
  CH_assert (a_src >= 0);
  CH_assert(a_src <  numProc());
  int isize;
  if (procID() == a_src)
  {
    isize = linearSize(a_inAndOut);
  }

  MPI_Bcast(&isize, 1, MPI_INT, a_src, Chombo_MPI::comm);

  void* broadBuf = malloc(isize);

  if (broadBuf == NULL)
  {
    Chombo4::MayDay::Error("out of memory in broadcast");
  }

  //take inAndOut from src and put it into broadBuf
  if (procID() == a_src)
  {
    linearOut(broadBuf, a_inAndOut);
  }

  //broadcast broadBuf to all procs
  MPI_Bcast(broadBuf, isize, MPI_BYTE, a_src, Chombo_MPI::comm);

  if (procID()==a_src)
  {
    CH_MaxMPISendSize = std::max<long long>(CH_MaxMPISendSize, isize);
  }
  else
  {
    CH_MaxMPIRecvSize = std::max<long long>(CH_MaxMPIRecvSize, isize);
  }
  //take broadBuf and put back into inAndOut if not src
  if (procID() != a_src)
  {
    linearIn(a_inAndOut, broadBuf);
  }

  //delete memory for buffer
  free(broadBuf);
}

/*****************************/
// simple Barrier
/*****************************/
inline void
barrier(void)
{
  MPI_Barrier(Chombo_MPI::comm);
}

#else
/*****************************/
//non-mpi version
/*****************************/
template <class T>
inline void
gather(Vector<T>& a_outVec, const T& a_input, int a_dest)
{
  a_outVec.resize(1);
  a_outVec[0] = a_input;
}
/*****************************/
//non-mpi version
/*****************************/
template <class T>
inline void
broadcast(T& a_inAndOut,  int a_src)
{
  //nothing to do.  in and out are the same with one proc
}
/*****************************/
//non-mpi version
/*****************************/
inline void
barrier(void)
{
  // do nothing in serial
}

#endif //the mpi thing

//*************************************
//These should work independent of MPI
//*************************************

//Vector<T> specialization of linearIn
template <class T>
void
linearListIn(Vector<T>& a_outputT, const void* const a_inBuf)
{
  //first entry is the size of the vector
  const int* const intBuf = (int*)a_inBuf;
  int vecsize = intBuf[0];
  Vector<int> vecOffset(vecsize);
  //next vecsize entries are offsets of data into buffer
  for (int ivec = 0; ivec < vecsize; ivec++)
    {
      vecOffset[ivec] = intBuf[ivec+1];
    }
  //next vecsize entries are the actual data
  //yes I could do this in one loop but that would
  // either
  // a) make it less symmetric with linearOut
  // and/or
  // b) make both of them far less readable
  a_outputT.resize(vecsize);
  const char* const charbuf = (char*)a_inBuf;
  for (int ivec = 0; ivec < vecsize; ivec++)
    {
      const char* const dataLoc = charbuf + vecOffset[ivec];
      linearIn(a_outputT[ivec], dataLoc);
    }
}

//Vector<T> specialization of linearOut
template <class T>
void
linearListOut(void* const a_outBuf, const Vector<T>& a_input)
{
  //first entry is the size of the vector
  int* const intBuf = (int*)a_outBuf;
  intBuf[0] = a_input.size();
  int vecsize = intBuf[0];
  Vector<int> vecOffset(vecsize);
  //next vecsize entries are offsets of data into buffer
  //next vecsize entries are the actual data
  int ioffset = (vecsize+1)*sizeof(int);
  for (int ivec = 0; ivec < vecsize; ivec++)
    {
      intBuf[ivec+1] = ioffset;
      vecOffset[ivec] = ioffset;
      ioffset += linearSize(a_input[ivec]);
    }
  //yes I could do this in one loop but that would
  // either
  // a) make it less symmetric with linearIn
  // and/or
  // b) make both of them far less readable
  char* const charBuf = (char*)a_outBuf;
  for (int ivec = 0; ivec < vecsize; ivec++)
    {
      char* const dataLoc = charBuf + vecOffset[ivec];
      linearOut(dataLoc, a_input[ivec]);
    }
}

//Vector<T> specialization of linearSize
template <class T>
int
linearListSize(const Vector<T>& a_input)
{
  //first entry is the size of the vector (int)
  //next vecsize entries are offsets of data into buffer (int)
  //next vecsize entries are the actual data
  int itotsize = (a_input.size() + 1)*sizeof(int);
  for (unsigned int ivec = 0; ivec < a_input.size(); ivec++)
    {
      itotsize += linearSize(a_input[ivec]);
    }
  return itotsize;
}
}


#endif
