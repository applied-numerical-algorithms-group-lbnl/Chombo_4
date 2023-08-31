#ifndef PROTO_MEMINFO_H
#define PROTO_MEMINFO_H

#include<map>
#include<iostream>
#include<vector>
#include <fstream>
#include<assert.h>
#include <sys/resource.h>
#include "Proto_MemType.H"

namespace Proto
{
  class describeInfo
  {
    public:
    /// describe operator+
    /**
       apply the operation + for each component of two describeInfo
       \param in a describeInfo we use to increment this describeInfo             
    */
    describeInfo& operator+=(describeInfo& in)
    {
      dataSize += in.dataSize;
      count    += in.count;
      return *this;
    }

    double dataSize;		///<  size of accumulated memory in bytes
    long unsigned int count;	///< number of malloc for a given position
  };

  // we use size_t because it's the type used in cudaMemGetInfo
  class traceInfo
  {
    public:
    typedef std::pair<size_t,size_t> pairType;

    /// insert function
    /**
       insert an element in m_info
       \param a_currentSize is the current memory size usage given with cudaGetMemInfo (device) or rusage (host)
    */
    void addElem(size_t a_currentSize)
    {
      if(m_maxSize < a_currentSize)
        m_maxSize = a_currentSize;
      pairType tracePoint = std::make_pair (a_currentSize, m_maxSize); 
      m_info.push_back(tracePoint);

    }

    /// output function
    /**
       All information are written in a txt file
       \param file is the txt file where we write the trace
    */
    void writeFile (std::ofstream &file)
    {
      for(auto &it : m_info)
      {
	file << it.first << " " << it.second << std::endl;
      }
    }

    /// accessor function
    /**
       Access to the component size 
    */
    size_t size()
    {
      return m_info.size();
    }

    private:
    size_t m_maxSize = 0;		///< store the maximum memory usage 
    std::vector<pairType> m_info;	///< store a vector of current and maximum memory usage
  };

	 
  template<MemType MEMTYPE=MEMTYPE_DEFAULT>
  class memInfo
  {
    public:
    /// constructor
    /**
       default constructor
    */
    memInfo() {}
    /// insert function
    /**
       Enter information in a 'describeInfo' element
       element are defined by their location: file + line
       if the position does not exist, we create it
       \param file is file position for a given Malloc
       \param line is line position in the file for a given Malloc
       \param dataSize is size for a given Malloc
    */
    void enterInfo(std::string file, unsigned int line, unsigned long int dataSize)
    {
      std::pair<std::string,unsigned int> position = std::make_pair(file,line);
      auto elem = m_info.find(position);
      if(elem != m_info.end())
      {
	auto& item = elem->second;
        double& mem = item.dataSize;
	long unsigned int& count = item.count;
        mem += dataSize;
	count++;
      }
      else
      {
	describeInfo item;
	item.dataSize = dataSize;
	item.count = 1;
        m_info[position] = item;
      }
    }

    /// insert function
    /**
       Enter information in a 'describeInfo' element
       element are defined by their position: file + line
       if the position does not already exist, we create it
       \param position is the position in the code for a given Malloc
       \param a_elem contains the related information for a position in an other mpi process (used with mpi)
    */
    void enterInfos(std::pair<std::string,unsigned int> a_position, describeInfo& a_elem)
    {
      auto elem = m_info.find(a_position);
      if(elem != m_info.end())
      {
	auto& item = elem->second;
        item += a_elem;
      }
      else
      {
        m_info[a_position] = a_elem;
      }
    }


    /// insert function
    /**
       add a memory trace point for the device memory
    */
    void addTracePoint()
    {
#ifdef PROTO_CUDA
	size_t tmp_free = 0;
	size_t tmp_tot = 0;
	protoMemGetInfo(&tmp_free,&tmp_tot);
	assert(tmp_tot > tmp_free);
	size_t tmp_current_size = tmp_tot - tmp_free;
	m_trace.addElem(tmp_current_size);
#endif
    }

    /// insert function
    /**
       add a memory trace point for the host memory
    */
    void addTracePointHost()
    {
      int who = RUSAGE_SELF;
      struct rusage usage;
      getrusage(who, &usage); // kilobytes
      size_t tmp_current_size = usage.ru_maxrss * 1000; // kilobytes to bytes
      m_trace.addElem(tmp_current_size);
    }

    /// check function
    /**
       check if a file exist. We avoid to rewrite a file when we print de memory trace.
    */
    bool fexists(std::string filename)
    {
      std::ifstream ifile(filename);
      return !(ifile.fail());
    }

    // display function 
    /**
       display information stored in m_info and write the memory trace
    */
    void printInfo()
    {
#ifdef PR_MPI_TMP
      mpiCollectInfo();
      int rank = procID();
      int size = numProc();
      if(rank == 0) // we display only the rank 0
      {
#endif
      std::cout << " print mem Info " << std::endl << std::endl;
      for(auto& it : m_info)
      {
        auto& position = it.first;
        std::string file = position.first;
        unsigned int line = position.second;
        double dataSize = it.second.dataSize;
        long unsigned int count = it.second.count;
        std::cout << " In " << file  << " line "   << line  << " we allocated " 
		  << dataSize * 1E-9 << " GB, in " << count << " allocation(s)"
		  << std::endl;
      }
#ifdef PR_MPI_TMP
      }
      std::string baseline = "trace-rank-" +  std::to_string(rank);
#else
      std::string baseline = "trace";
#endif
      std::string filename = baseline + ".txt";
      int i = 1;
      while(fexists(filename))
      {
        filename = baseline + "-" + std::to_string(i) + ".txt";
	i++;
      }
      std::cout << " write memory trace in " << filename  << std::endl << std::endl;
      std::ofstream file(filename);
      m_trace.writeFile(file);
    }

    /**
       trick to get only one instance of memInfo for a given MEMTYPE
    */
    static memInfo<MEMTYPE>& getMemInfo()
    {
      static memInfo global_info;
      return global_info;
    }
#ifdef PR_MPI_TMP
    /**
       collect the data of other mpi processes contained in m_info
    */
    void mpiCollectInfo()
    {
      int rank = procID();
      int size = numProc();
      if(size == 1) return; // only one MPI process
      // get m_size and m_trace size
      std::vector<int> info_sizes;
      if(rank==0)
      {
        info_sizes.resize(size);
      }
      size_t mySize = m_info.size();
      MPI_Gather(  &mySize, 
		      1, 
		      MPI_INT, 
		      info_sizes.data(), 
		      1, 
		      MPI_INT, 
		      0, 
		      Proto_MPI<void>::comm);

      // pack data
      if(rank!= 0)
      {
	std::vector<std::pair<std::pair<std::string,unsigned int>, describeInfo>> toSend;
	for(auto it:m_info)
	{
	  auto item = std::make_pair(it.first,it.second);
	  toSend.push_back(item);
	}
	size_t nBytes = sizeof(std::pair<std::pair<std::string,unsigned int>, describeInfo>)*mySize;
        MPI_Send(
			toSend.data(),
			nBytes,
			MPI_CHAR,
			0,
			666+rank, 
			Proto_MPI<void>::comm);
      }
      else
      {
        for(size_t i = 0 ; i < size ; i++) std::cout << " rank: " << rank << " size "<< info_sizes[i] << std::endl;

        for(size_t otherRank = 1 ; otherRank < size ; otherRank++)
	{
	  std::vector<std::pair<std::pair<std::string,unsigned int>, describeInfo>> toRecv;
          size_t otherSize = info_sizes[otherRank];
	  toRecv.resize(otherSize);
	  size_t nBytes = sizeof(std::pair<std::pair<std::string,unsigned int>, describeInfo>)*otherSize;
          MPI_Status status;
	  MPI_Recv(
			  toRecv.data(),
			  nBytes,
			  MPI_CHAR,
			  otherRank,
			  666+otherRank,
			  MPI_COMM_WORLD,
			  &status);

	  // unpack data
	  for(auto it : toRecv)
	  {
            enterInfos(it.first,it.second);
	  }

	}
      }
    }
#endif //PR_MPI

    private:
    std::map<std::pair<std::string,unsigned int>, describeInfo> m_info; ///< accumulate data size allocations for each malloc position in the code 
    traceInfo m_trace; 							///< store the current and maximum memory footprint for each malloc
  }; // end memInfo class

} // end PROTO namespace

// Display memory footprint per position in the code
// Extended to host and device allocations
// This function works only if the storeMemInfo is used
//turned off to kill compiler warning that it was not used
#if 0
static void PRINT_MEMORY_INFO()
{	
#ifdef memoryFootPrint
#ifdef PROTO_CUDA
  std::cout << " -------------- device ----------- " << std::endl << std::endl;
  Proto::memInfo<Proto::DEVICE>& tmpInfoDevice = Proto::memInfo<Proto::DEVICE>::getMemInfo(); 
  tmpInfoDevice.printInfo();
  std::cout << " --------------- host ------------ " << std::endl << std::endl;
#endif
  Proto::memInfo<Proto::HOST>& tmpInfoHost = Proto::memInfo<Proto::HOST>::getMemInfo(); 
  tmpInfoHost.printInfo();
#else
  std::cout << " The memory footprint option is not activated. flag: -DmemoryFootPrint (performance overhead if activated)" << std::endl;
#endif
}
#endif


// here are the macros defined to check the memory used in Proto_gpu.H and Proto_cpu.H
#ifdef memoryFootPrint
	// tmp_ to avoid doublons 
        #define countMallocDevice(cond) \
	{	\
		cond;\
	        Proto::memInfo<Proto::DEVICE>& tmpInfo = Proto::memInfo<Proto::DEVICE>::getMemInfo();\
	        tmpInfo.addTracePoint();\
	};

        #define countMallocHost(cond) \
	{	\
		cond;\
	        Proto::memInfo<Proto::HOST>& tmpInfo = Proto::memInfo<Proto::HOST>::getMemInfo();\
	        tmpInfo.addTracePointHost();\
	};

	#define storeMemInfo(MEM,B) \
	{	\
		using Proto::HOST;\
		using Proto::DEVICE;\
		Proto::memInfo<MEM>& tmpInfo = Proto::memInfo<MEM>::getMemInfo();\
		tmpInfo.enterInfo(__FILE__,__LINE__,B);\
	};
#else
	#define countMallocDevice(cond) cond
	#define countMallocHost(cond) cond
	#define storeMemInfo(MEMTYPE,BYTES) 
#endif
#endif //PROTO_MEMINFO_H
