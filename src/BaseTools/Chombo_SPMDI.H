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
  template <class T>
  inline void  gatherSizes(vector<size_t>& a_allSizes, size_t& a_totalSize,
                           const T& a_input, int a_dest, bool a_printStuff)
  {//gather the sizes to be received from each proc
    //now THIS size lives on THIS processor
    //using unsigned long because mpi understands it
    size_t thisProcSize_t = linearSize(a_input);
    unsigned long int  thisProcSize = thisProcSize_t;
    int nProcess = CH4_SPMD::numProc();

    unsigned long int*  recvbuf = NULL;
    //allocate received buffer
    if (procID() == a_dest)
    {
      recvbuf = new unsigned long int[nProcess];
    }

    if(a_printStuff)
    {
      Chombo4::pout() << "CH4_SPMD::gatherSizes: before mpi_gather" << endl;
    }

    MPI_Gather(&thisProcSize, 1, MPI_UNSIGNED_LONG, recvbuf, 1, MPI_UNSIGNED_LONG, a_dest, Chombo_MPI::comm);
      
    if(a_printStuff)
    {
      Chombo4::pout() << "CH4_SPMD::gatherSizes: after mpi_gather" << endl;
    }
    if(procID() == a_dest)
    {
      a_allSizes.resize(nProcess);
      char* bufloc = (char*) recvbuf;
      for(int iproc = 0; iproc < nProcess; iproc++)
      {
        size_t* sizeptr = (size_t*) bufloc;
        a_allSizes[iproc] = *sizeptr;
        a_totalSize += a_allSizes[iproc];
          
        bufloc += sizeof(unsigned long int); //remember mpi does not know from size_t
      }
      if(a_printStuff)
      {
        for(int iproc = 0; iproc < nProcess; iproc++)
        {
          Chombo4::pout() << "CH4_SPMD::gatherSizes: gathered sizes : " ;
          Chombo4::pout() << "iproc = " << iproc << ", size = " << a_allSizes[iproc] << endl;
        }
      }
    }
    if(a_printStuff)
    {
      Chombo4::pout() << "CH4_SPMD::gatherSizes: cleaning up" << endl;
    }
    if (procID() == a_dest)
    {
      delete[] recvbuf;
    }
    if(a_printStuff)
    {
      Chombo4::pout() << "CH4_SPMD::gatherSizes: leaving" << endl;
    }
  }
/*****************************/
//gather a_input into a a_outVec
/*****************************/
  template <class T>
  inline void
  gather(std::vector<T>& a_outVec, const T& a_input, int a_dest, bool a_printStuff)
  {
    CH_assert (a_dest >= 0);
    CH_assert(a_dest <  numProc());
    int nProcess = CH4_SPMD::numProc();
    vector<size_t> allSizes(nProcess);
    size_t totalSize = 0;
    //gather the sizes to be received from each proc
    gatherSizes(allSizes, totalSize, a_input, a_dest, a_printStuff);

    char* recvbuf = NULL;

    //mpi wants ints. mpi gets ints
    int* vectSize = NULL;
    int* vectDisp = NULL;
    //allocate received buffer and other sundries
    if (procID() == a_dest)
    {
      vectSize = new int[nProcess];
      vectDisp = new int[nProcess];
      recvbuf = (char*) malloc(totalSize);
      
      size_t itotsize=0;
      for (int iproc = 0; iproc < nProcess; iproc++)
      {
        vectSize[iproc] = allSizes[iproc];
        vectDisp[iproc] = itotsize;
        itotsize       += allSizes[iproc];
      }
      if(itotsize != totalSize)
      {
        Chombo4::pout() << "CH4_SPMD::gather:: WARNING itotsize = "<< itotsize << ", totalSize = " << totalSize << "(they should be equal)" << endl;
      }
    }

    size_t size_t_size = linearSize(a_input);
    //mpi wants ints. mpi gets ints
    int isize = size_t_size;
    //make stuff for linearout
    char* localBuf = (char*)malloc(isize);
    linearOut(localBuf, a_input);
    
    if(a_printStuff)
    {
      Chombo4::pout() << "CH4_SPMD::gather:: going into gatherv" << endl;
    }

    MPI_Gatherv(localBuf, isize, MPI_BYTE,
                recvbuf, vectSize, vectDisp, MPI_BYTE,
                a_dest, Chombo_MPI::comm);

    if(a_printStuff)
    {
      Chombo4::pout() << "CH4_SPMD::gather:: out of gatherv " << endl;
    }
    
    if (procID() == a_dest)
    {
      a_outVec.resize(nProcess);
      char* bufloc = recvbuf;
      for (int iproc = 0; iproc < nProcess; iproc++)
      {
        linearIn(a_outVec[iproc], bufloc);
        bufloc += vectDisp[iproc];
      }

      //delete memory for dest-specific arrays
      delete[] vectSize;
      delete[] vectDisp;
      free(recvbuf);
    }

    //delete memory for local buffer
    free(localBuf);
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
    size_t size_t_size;
    if (procID() == a_src)
    {
      size_t_size = linearSize(a_inAndOut);
    }

    //mpi wants ints.  mpi gets ints.
    int isize = size_t_size;
    
    if(a_printStuff)
    {
      Chombo4::pout() << "CH4_SPMD::broadcast: before MPI_Bcast" << endl;
    }
    MPI_Bcast(&isize, 1, MPI_INT, a_src, Chombo_MPI::comm);
    
    if(a_printStuff)
    {
      Chombo4::pout() << "CH4_SPMD::broadcast: after MPI_Bcast "
                      << ", isize = " << isize << endl;      
    }
    
    char* broadBuf = (char *) malloc(isize);

    //take inAndOut from src and put it into broadBuf
    if (procID() == a_src)
    {
      if(a_printStuff)
      {
        Chombo4::pout() << "CH4_SPMD::broadcast: before linearOut " << endl;
      }
      linearOut(broadBuf, a_inAndOut);
      if(a_printStuff)
      {
        Chombo4::pout() << "CH4_SPMD::broadcast: after linearOut " << endl;
      }
    }

    //broadcast broadBuf to all procs
    if(a_printStuff)
    {
      Chombo4::pout() << "CH4_SPMD::broadcast: before second mpi_bcast " << endl;
    }

    MPI_Bcast(broadBuf, isize, MPI_BYTE, a_src, Chombo_MPI::comm);
    if(a_printStuff)
    {
      Chombo4::pout() << "CH4_SPMD::broadcast: after second mpi_bcast " << endl;
    }

    //take broadBuf and put back into inAndOut if not src
    if (procID() != a_src)
    {

      Chombo4::pout() << "CH4_SPMD::broadcast: going into linearin " << endl;
      linearIn(a_inAndOut, broadBuf);
      Chombo4::pout() << "CH4_SPMD::broadcast: coming out of linearin " << endl;
    }

    if(a_printStuff)
    {
      Chombo4::pout() << "CH4_SPMD::broadcast: just before freeing memory " << endl;
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
