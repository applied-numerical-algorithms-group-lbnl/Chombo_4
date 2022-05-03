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

/*****************************/
//gather a_input into a a_outVec
/*****************************/
  template <class T>
  inline void
  gather(std::vector<T>& a_outVec, const T& a_input, int a_dest, bool a_printStuff)
  {
    CH_assert (a_dest >= 0);
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
  broadcast(T& a_inAndOut,  int a_src, bool a_printStuff)
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
  gather(std::vector<T>& a_outVec, const T& a_input, int a_dest, bool a_printStuff)
  {
    a_outVec.resize(1);
    a_outVec[0] = a_input;
  }
/*****************************/
//non-mpi version
/*****************************/
  template <class T>
  inline void
  broadcast(T& a_inAndOut,  int a_src, bool a_printStuff)
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

//std::vector<T> specialization of linearIn
  template <class T>
  void
  linearListIn(std::vector<T>& a_outputT,
               size_t a_expectedSize,
               const void* const a_inBuf,
               bool a_printStuff)
  {
    //first entry is the size of the vector
    char*  charbuf = (char *) a_inBuf;
    size_t charsiz = 0;

    Chombo4::pout() << "CH4_SPMD::linearListIn: entering " << endl;

    size_t* vecsizeptr = (size_t*)(charbuf);
    size_t vecsize  = *vecsizeptr;
    Chombo4::pout() << "CH4_SPMD::linearListIn: vecsize= " << vecsize << endl;
    a_outputT.resize(vecsize);
  
    charbuf += sizeof(size_t);
    charsiz += sizeof(size_t);

    std::vector<size_t> vecOffset(vecsize);
    for (int ivec = 0; ivec < vecsize; ivec++)
    {
      size_t * sizebuf = (size_t*) charbuf;
      size_t itemsize = *sizebuf;
      vecOffset[ivec] = itemsize;
      Chombo4::pout() << "CH4_SPMD::linearListIn: "
                      << ", ivec = " << ivec
                      << ", vecOffset[ivec] = " << vecOffset[ivec] << endl;
      
      charbuf += sizeof(size_t);
      charsiz += sizeof(size_t);
    }

    for (int ivec = 0; ivec < vecsize; ivec++)
    {
      const char* const dataLoc = charbuf;
      linearIn(a_outputT, dataLoc);
      charbuf += vecOffset[ivec];
      charsiz += vecOffset[ivec];
    }

    if(charsiz != a_expectedSize)
    {
      Chombo4::pout() << "CH4_SPMD::linearListIn: WARNING size mismatch" << endl;
    }
     
    Chombo4::pout() << "CH4_SPMD::linearListIn: leaving " << endl;
  }

//std::vector<T> specialization of linearOut
  template <class T>
  void
  linearListOut(void* const a_outBuf,
                size_t a_expectedSize,
                const std::vector<T>& a_input,
                bool a_printStuff)
  {
    //first entry is the size of the vector
    char*  charbuf = (char *) a_outBuf;
    size_t charsiz = 0;  
    size_t* vecsizeptr = (size_t*)(charbuf);
    size_t vecsize  = a_input.size();
    *vecsizeptr = vecsize;
  
    charbuf += sizeof(size_t);
    charsiz += sizeof(size_t);

    std::vector<size_t> vecOffset(vecsize);
    for (int ivec = 0; ivec < vecsize; ivec++)
    {
      size_t * sizebuf = (size_t*) charbuf;
      size_t itemsize = linearSize(a_input[ivec]);
      *sizebuf        = itemsize;
      vecOffset[ivec] = itemsize;
    
      charbuf += sizeof(size_t);
      charsiz += sizeof(size_t);
    }

    for (int ivec = 0; ivec < vecsize; ivec++)
    {
      linearOut(charbuf, a_input[ivec]);
      charbuf += vecOffset[ivec];
      charsiz += vecOffset[ivec];
    }


    if(charsiz != a_expectedSize)
    {
      Chombo4::pout() << "CH4_SPMD::linearListOut: WARNING size mismatch:";
      Chombo4::pout() << " charsiz = " << charsiz ;
      Chombo4::pout() << " expected = " << a_expectedSize << endl;
    }
  }

//std::vector<T> specialization of linearSize
  template <class T>
  size_t
  linearListSize(const std::vector<T>& a_input,
                 bool a_printStuff)
  {
    //first entry is the size of the vector 
    //next vecsize entries are offsets of data into buffer 
    //next vecsize entries are the actual data
    //written pedantically but I find it comforting.
    size_t itotsize = sizeof(size_t);
    for (unsigned int ivec = 0; ivec < a_input.size(); ivec++)
    {
      itotsize += sizeof(size_t);
    }
    for (unsigned int ivec = 0; ivec < a_input.size(); ivec++)
    {
      itotsize += linearSize(a_input[ivec]);
    }
    return itotsize;
  }
}


#endif
