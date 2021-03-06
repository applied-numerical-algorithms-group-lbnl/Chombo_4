
#ifndef _Proto_IrregDATA_H_
#define _Proto_IrregDATA_H_

#include <cmath>
#include <cstdlib>
#include <memory>
#include "Proto_Box.H"
#include "EBProto.H"
#include "Proto_EBIrregStruct.H"

using std::vector;
using std::shared_ptr;
namespace Proto
{

    template <CENTERING cent>
    CUDA_DECORATION
    size_t  index(EBIndex<cent>* a_indicies, unsigned int a_size, EBIndex<cent>& a_vof, int a_ivar) 
    {
      size_t ivec = 0;
      bool found = false;
      int begin = 0;
      int end = a_size;

      while (begin <= end) 
      {
        //int mid = begin + (begin - end) / 2;
        int mid = (begin + end) / 2;

        if (a_indicies[mid] == a_vof)
        {
          ivec = mid;
          found = true;
          break;
        }

        if (a_indicies[mid] < a_vof)
          begin = mid + 1;
        else
         end = mid - 1;
      }

      if(!found)
      {
        PROTO_ASSERT(found, "vof not found");
      }

      ivec += a_ivar* a_size;

      return ivec;
    }

    CUDA_DECORATION
    inline bool contains(const Point& pt, const Point& a_low, const Point& a_high);

    template<CENTERING cent>
    CUDA_DECORATION
    inline bool contains(EBIndex<cent>& a_in, const Point& a_low, const Point& a_high);

#ifndef PROTO_CUDA
#define FORCE_INLINE inline
#else  
#define FORCE_INLINE __forceinline__
#endif  
    template <CENTERING cent>
    CUDA_DECORATION
    FORCE_INLINE bool copyContains(const EBIndex<cent>& a_ind, const Point & a_srcbox_low, const Point & a_srcbox_high)
    {
      CENTERING thiscent = cent;

      bool retval = false;
      if((thiscent == CELL) || (thiscent == BOUNDARY))
      {
        retval = contains(a_ind.m_pt,a_srcbox_low,a_srcbox_high);
      }
      else if(thiscent == XFACE)
      {
        Point hipt = a_ind.m_pt;
        Point lopt = a_ind.m_pt - Point::Basis(0, 1);
        retval = (contains(lopt,a_srcbox_low,a_srcbox_high) 
               || contains(hipt,a_srcbox_low,a_srcbox_high));
      }
      else if(thiscent ==  YFACE)
      {
        Point hipt = a_ind.m_pt;
        Point lopt = a_ind.m_pt - Point::Basis(1, 1);
        retval = (contains(lopt,a_srcbox_low,a_srcbox_high) 
               || contains(hipt,a_srcbox_low,a_srcbox_high));
      }
#if DIM==3
      else if(thiscent ==  ZFACE)
      {
        Point hipt = a_ind.m_pt;
        Point lopt = a_ind.m_pt - Point::Basis(2, 1);
        retval = (contains(lopt,a_srcbox_low,a_srcbox_high) 
               || contains(hipt,a_srcbox_low,a_srcbox_high));
      }
#endif
      else
      {
        PROTO_ASSERT(false, "unaccounted-for centering");
      }
      return retval;
    }

  template <CENTERING cent, class data_t>
  struct EBBoxDataIndexer
  {
    static
    inline void
    cpu(
                unsigned int a_end, 
		EBIndex<cent>* a_indiciesSrc, unsigned int a_indiciesSizeSrc, 
		EBIndex<cent>* a_indiciesDst, unsigned int a_indiciesSizeDst,
		data_t* a_dst, const data_t* a_src, 
                unsigned int a_idst, unsigned int a_isrc,
                unsigned int a_inco,
		const Point& a_box_inf, const Point& a_box_high
    )
    {
      for(unsigned idx = 0; idx < a_end ; idx ++)
      {
        EBIndex<cent>& vof = a_indiciesDst[idx];
        if( contains(vof,a_box_inf,a_box_high) && copyContains(vof,a_box_inf,a_box_high))
        {
          for(unsigned int ico = 0; ico < a_inco; ico++)
          {
            unsigned int isrc = a_isrc + ico;
            unsigned int idst = a_idst + ico;
            size_t indDst = index(a_indiciesDst, a_indiciesSizeDst, vof, idst);
            size_t indSrc = index(a_indiciesSrc, a_indiciesSizeSrc, vof, isrc);
            a_dst[indDst] = a_src[indSrc];
          }
        }

      }
    }

#ifdef PROTO_CUDA
    __device__ static
    void
    gpu(
                unsigned int a_end, 
		EBIndex<cent>* a_indiciesSrc, unsigned int a_indiciesSizeSrc, 
		EBIndex<cent>* a_indiciesDst, unsigned int a_indiciesSizeDst,
		data_t* a_dst, const data_t* a_src, 
                unsigned int a_idst, unsigned int a_isrc,
		unsigned int a_inco,
		const Point& a_box_inf, const Point& a_box_high
    )
    {
/*  
    unsigned int idx = blockIdx.x * blockDim.x + threadIdx.x;

      if(idx < a_end)
      {
        EBIndex<cent>& vof = a_indiciesDst[idx];
        const bool do_it = copyContains(vof,a_box_inf,a_box_high);
        //if( contains(vof,a_box_inf,a_box_high) && copyContains(vof,a_box_inf,a_box_high))
        if(do_it)
        {
          size_t indDst = idx + a_indiciesSizeDst*a_idst;//index(a_indiciesDst, a_indiciesSizeDst, vof, idst);
          size_t indSrc = index(a_indiciesSrc, a_indiciesSizeSrc, vof, a_isrc);
	  __syncthreads();

#pragma unroll
          for(unsigned int ico = 0; ico < a_inco; ico++)
          {
            a_dst[indDst] = a_src[indSrc];
            indDst += a_indiciesSizeDst;
            indSrc += a_indiciesSizeSrc;
          }
        }
      }
*/
    
      for(unsigned int idx = blockIdx.x * blockDim.x + threadIdx.x ; idx < a_end ; idx += blockDim.x  * gridDim.x)
      {
        EBIndex<cent>& vof = a_indiciesDst[idx];
        const bool do_it = copyContains<cent>(vof,a_box_inf,a_box_high) && contains(vof,a_box_inf,a_box_high);
        if(do_it)
        {
          size_t indDst = idx + a_indiciesSizeDst*a_idst;
          size_t indSrc = index(a_indiciesSrc, a_indiciesSizeSrc, vof, a_isrc);

          for(unsigned int ico = 0; ico < a_inco; ico++)
          {
            a_dst[indDst] = a_src[indSrc];
            indDst += a_indiciesSizeDst;
            indSrc += a_indiciesSizeSrc;
          }
        }
      }
  }
#endif
  };
   
  template <CENTERING cent, class data_t, unsigned int ncomp>
  struct EBLinearIn
  {
    static
    inline void
    cpu(
    		EBIndex<cent>* a_indicies, unsigned int a_indiciesSize,
                void* a_buf, 
		data_t* a_dst,
                const Point& a_box_inf, const Point& a_box_high
    )
    {
      char* tmp = (char*)a_buf;
      unsigned int size = ((unsigned int*) a_buf)[0];
      tmp += sizeof(data_t); // we need to aligne double*
      data_t* data = (data_t*) tmp;
      tmp += ncomp*size*sizeof(data_t);
      EBIndex<cent>* indicies = (EBIndex<cent>*) tmp;

      for(unsigned int idx = 0; idx < size ; idx ++)
      {
        EBIndex<cent>& pt = indicies[idx];
        if( contains(pt,a_box_inf,a_box_high) && copyContains(pt,a_box_inf,a_box_high))
        {
          for(int icomp = 0; icomp < ncomp; icomp++)
          {
            size_t indxDst = index(a_indicies, a_indiciesSize, pt, icomp);
            //size_t indxSrc = index(indicies, size, pt, icomp);
            //a_dst[indxDst+icomp*a_indiciesSize] = data[indxSrc+icomp*size];
            a_dst[indxDst+icomp*a_indiciesSize] = data[idx+icomp*size];
          }
        }
      }
   }

#ifdef PROTO_CUDA
    __device__ static
    void
    gpu(
		EBIndex<cent>* a_indicies, unsigned int a_indiciesSize,
                void* a_buf, 
		data_t* a_dst,
                const Point& a_box_inf, const Point& a_box_high
    )
    {
      char* tmp = (char*)a_buf;
      unsigned int size = ((unsigned int*) a_buf)[0];
      tmp += sizeof(data_t); // we need to aligne double*
      data_t* data = (data_t*) tmp;
      tmp += ncomp*size*sizeof(data_t);
      EBIndex<cent>* indicies = (EBIndex<cent>*) tmp;

      for(unsigned int idx = blockIdx.x * blockDim.x + threadIdx.x ; idx < size ; idx += blockDim.x  * gridDim.x)
      {
        EBIndex<cent>& pt = indicies[idx];
        if( contains(pt,a_box_inf,a_box_high) && copyContains(pt,a_box_inf,a_box_high))
        //if(copyContains<cent>(pt,a_box_inf,a_box_high))
        {
          size_t indxDst = index(a_indicies, a_indiciesSize, pt, 0);
          size_t indxSrc = idx;
#pragma unroll
          for(int icomp = 0; icomp < ncomp; icomp++)
          {
            a_dst[indxDst] = data[indxSrc];
            indxDst += a_indiciesSize;
            indxSrc += size;
          }
        }
      }
   }
#endif
 };

  template <CENTERING cent, class data_t, unsigned int ncomp>
  struct EBLinearOut
  {
    static
    inline void
    cpu(
    		unsigned int a_end,
                EBIndex<cent>* a_indicies, std::size_t a_indiciesSize,
                void* a_buf, 
	        const data_t* a_src
    )
    {

      char* tmp = (char*)a_buf;
      unsigned int *size = (unsigned int*) (tmp);
      *size = a_end;
      tmp += sizeof(data_t); // we need to align double
      data_t* ptr = (data_t*) (tmp);
      tmp += ncomp*a_end*sizeof(data_t);
      EBIndex<cent>* indicies = (EBIndex<cent>*) (tmp);

      for(unsigned int idx = 0 ; idx < a_end ; idx++)
      {
        indicies[idx] = a_indicies[idx];
        for(int icomp = 0; icomp < ncomp; icomp++)
        {
          ptr[idx+icomp*a_end] = a_src[idx+icomp*a_indiciesSize];
          //ptr[idx] = a_src[idx];
        }
      }
    }

#ifdef PROTO_CUDA
    __device__ static
    void
    gpu(
    		unsigned int a_end,
                EBIndex<cent>* a_indicies, std::size_t a_indiciesSize,
                void* a_buf, 
	        const data_t* a_src
    )
    {
      unsigned int idx = blockIdx.x * blockDim.x + threadIdx.x;

      char* tmp = (char*)a_buf;

      if(threadIdx.x == 0 && blockIdx.x == 0)
      {
        unsigned int *size = (unsigned int*) (tmp);
        *size = a_end;
      }
      tmp += sizeof(data_t); // we need to align double
      data_t* ptr = (data_t*) (tmp);
      tmp += ncomp*a_end*sizeof(data_t);
      EBIndex<cent>* indicies = (EBIndex<cent>*) (tmp);

      if(idx < a_end)
      {
        indicies[idx] = a_indicies[idx];
#pragma unroll
        for(int icomp = 0; icomp < ncomp; icomp++)
        {
          ptr[idx+icomp*a_end] = a_src[idx+icomp*a_indiciesSize];
        }
      }
    }
#endif
  };


///
/**
 */
  template <CENTERING cent, class data_t,  unsigned int ncomp>
  class IrregData
  {
  public:
    typedef EBIndex<cent> index_t;
    ///
    IrregData() 
    {
      m_boxdefined = false;
    }

    /// constructor by copying
    IrregData(Box a_box, data_t* a_data, std::vector<index_t>& a_indicies) : m_boxdefined(true)
    {
      m_inputBox = a_box;
      unsigned int size = a_indicies.size();
      protoMalloc(MEMTYPE_DEFAULT,m_rawPtr,size*sizeof(data_t));
      protoMemcpy(MEMTYPE_DEFAULT,m_rawPtr, a_data, size*sizeof(data_t), protoMemcpyHostToDevice);
      m_data = ::std::shared_ptr<data_t>(m_rawPtr, [](data_t* p){protoFree(MEMTYPE_DEFAULT,p);});
      m_indicies = ::std::shared_ptr<std::vector<index_t>>(new std::vector<index_t>[1], [](std::vector<index_t>* p) { delete[] p;});
      m_indicies->resize(size);
      for(int i = 0; i < size ; i++)
        (*m_indicies)[i] = a_indicies[i];

      defineEBIrregData();
    }

    ///
    ~IrregData()
    //virtual ~IrregData()
    {
      PR_TIME("~IrregData");
      m_EBIrregDataSharedPtr.reset();
      m_indicies_shared_data.reset();
      m_data.reset();
    }

    ///needed for leveldata
    IrregData(const Box& a_box, bool a_stackAllocation=true)
    {  }

    ///
    inline std::size_t size() const 
    {
      if(m_boxdefined) return m_indicies->size()*ncomp;
      else return 0; 
    };

    ///
    inline bool defined() const {return bool(m_data);};
   

    ///
    inline void setVal(const data_t& a_val);
  

    ///
    data_t* dataPtr(unsigned int a_ivar = 0)
    {
      PR_assert(a_ivar < ncomp);
      size_t ioff = a_ivar*m_indicies->size();
      return &(m_rawPtr[ioff]);
    }

    ///
    const data_t* dataPtr(unsigned int a_ivar = 0) const
    {
      PR_assert(a_ivar < ncomp);
      size_t ioff = a_ivar*m_indicies->size();
      return &(m_rawPtr[ioff]);
    }

    data_t* data()
    {
      return m_rawPtr;
    }


    const data_t* data() const
    {
      return m_rawPtr;
    }
    
    ///
    inline size_t  index(const EBIndex<cent>& a_vof, const int& a_ivar) const
    {
      size_t ivec = 0;
      bool found = false;
      for(unsigned int i = 0; i < m_indicies->size(); i++)
      {        
        if(a_vof == (*m_indicies)[i])
        {
          ivec = i;
          found = true;
          break;
        }
      }
      if(!found)
      {
        std::cout << "index<cell> vof not found" << endl;
        abort();
      }
      PROTO_ASSERT(found, "vof not found");

      ivec += a_ivar*m_indicies->size();

      if(ivec >= size())
      {
        std::cout << "index<cell> offset too big" << endl;
        abort();
      }

      return ivec;
    }

    IrregData<cent, data_t, ncomp>& 
    operator+=(const IrregData<cent, data_t, ncomp> & a_input)
    {
      if((m_inputBox == a_input.m_inputBox))
      {
        unsigned int nsize = size();
        if(nsize == 0) return *this;
        int blockSize = 512;
        int nBlocks = (nsize + blockSize - 1) / blockSize;
        //protoLaunchKernelMemAsync((opKernel<data_t,BoxDataOp::Add>), nBlocks, blockSize, 0, protoGetCurrentStream, m_rawPtr, a_input.m_rawPtr, nsize);
        protoLaunchKernelMemAsyncT<MEMTYPE_DEFAULT, opKernel<data_t, BoxDataOp::Add>>(nBlocks, blockSize, 0, protoGetCurrentStream, m_rawPtr, a_input.m_rawPtr, nsize);
      }
      else
      {
        PROTO_ASSERT(false, "operator+= not implemented for irreg datas of different sizes");
      }
      return *this;
    }


    IrregData<cent, data_t, ncomp>& 
    operator-=(const IrregData<cent, data_t, ncomp> & a_input)
    {
      if(m_boxdefined && (m_inputBox == a_input.m_inputBox))
      {
        unsigned int nsize = size();
        if(nsize == 0) return *this;
        int blockSize = 512;
        int nBlocks = (nsize + blockSize - 1) / blockSize;
        //protoLaunchKernelMemAsync((opKernel<data_t,BoxDataOp::Subtract>), nBlocks, blockSize, 0, protoGetCurrentStream, m_rawPtr, a_input.m_rawPtr, nsize);
        protoLaunchKernelMemAsyncT<MEMTYPE_DEFAULT, opKernel<data_t,BoxDataOp::Subtract>>(nBlocks, blockSize, 0, protoGetCurrentStream, m_rawPtr, a_input.m_rawPtr, nsize);
      }
      else
      {
        PROTO_ASSERT(false, "operator+= not implemented for irreg datas of different sizes");
      }
      return *this;
    }


    IrregData<cent, data_t, ncomp>& 
    operator*=(const IrregData<cent, data_t, ncomp> & a_input)
    {
      if(m_boxdefined && (m_inputBox == a_input.m_inputBox))
      {
        unsigned int nsize = size();
        if(nsize == 0) return *this;
        int blockSize = 512;
        int nBlocks = (nsize + blockSize - 1) / blockSize;
        //protoLaunchKernelMemAsync((opKernel<data_t,BoxDataOp::Multiply>), nBlocks, blockSize, 0, protoGetCurrentStream, m_rawPtr, a_input.m_rawPtr, nsize);
        protoLaunchKernelMemAsyncT<MEMTYPE_DEFAULT, opKernel<data_t, BoxDataOp::Multiply>>(nBlocks, blockSize, 0, protoGetCurrentStream, m_rawPtr, a_input.m_rawPtr, nsize);
      }
      else
      {
        PROTO_ASSERT(false, "operator+= not implemented for irreg datas of different sizes");
      }
      return *this;
    }
    
    /// debugging only--will not work on the device
    data_t operator() (const index_t& a_vof, int a_comp) const 
    {
      size_t ind = index(a_vof, a_comp);
      return m_rawPtr[ind];
    }

    vector<EBIndex<cent> > getIndices() const
    {
      return *m_indicies;
    }

    unsigned int vecsize() const
    {
      if(m_boxdefined) return m_indicies->size();
      else return 0;
    }

    inline Box inputBox() const
    {
      return m_inputBox;
    }
    
    //aliasing define
    template<unsigned int srccomp>
    void define(IrregData<cent, data_t, srccomp>& a_input,
                unsigned int                    & a_comp)
    {
      m_boxdefined = true;
      m_indicies       = a_input.getIndicies();
      m_inputBox       = a_input.inputBox();
      m_data           = a_input.getData();
      m_rawPtr         = a_input.dataPtr(a_comp);
      EBIrregStruct<cent,data_t,ncomp> *tmp = (EBIrregStruct<cent,data_t,ncomp>*) a_input.getEBIrregDataSharedPtr().get();
      m_EBIrregDataSharedPtr = ::std::shared_ptr<EBIrregStruct<cent,data_t,ncomp>>(tmp, [](EBIrregStruct<cent,data_t,ncomp>* p){});


      //m_EBIrregDataSharedPtr  = dynamic_pointer_cast<EBIrregStruct<cent,data_t,ncomp>>(a_input.getEBIrregDataSharedPtr());
      m_indicies_shared_data = a_input.getIndiciesData();
    }           

    inline shared_ptr<EBIrregStruct<cent, data_t, ncomp>> getEBIrregDataSharedPtr() const
    {
      return m_EBIrregDataSharedPtr;
    }


    //do not use this
    shared_ptr<vector<  index_t > >  getIndicies() const
    {
      return m_indicies;
    }

    //do not use this
    shared_ptr<data_t> getData() const
    {
      return  m_data;
    }

    shared_ptr<index_t> getIndiciesData()
    {
      return m_indicies_shared_data;
    }

    inline void define(const Box& a_box,  const EBGraph& a_graph, bool a_useStack = false)
    {
      this->m_boxdefined = true;
      this->m_inputBox   = a_box;
      Box interBox = a_box & a_graph.getDomain();
      shared_ptr<vector<EBIndex<cent> > > indicies(new vector< EBIndex<cent> >() );
      *indicies = getIrregIndicies<cent>(interBox, a_graph);

      define(indicies, a_useStack);
      defineEBIrregData(a_useStack);
    }

    inline void defineEBIrregData(bool a_stackAllocation = false)
    {
      PR_TIME("defineEBIrregData");
      vector< EBIrregStruct<cent, data_t, ncomp> > hostvec = getEBIrregStruct(*m_indicies, *this);

      if(a_stackAllocation)
      {
        PR_TIME("defineEBIrregData::StackAllocation");
        StackAlloc<MEMTYPE_DEFAULT>& allocator  = Proto::DisjointBoxLayout::getStackAllocator<MEMTYPE_DEFAULT>();
	// IrregData can be empty, but we can't do zero-allocation on a stack.
	if( hostvec.size() == 0)
	{
		m_EBIrregDataPtr = nullptr;
		m_EBIrregDataSharedPtr = ::std::shared_ptr<EBIrregStruct<cent,data_t,ncomp>>(nullptr);
	}
	else
	{
		size_t memsize = hostvec.size()*sizeof(EBIrregStruct<cent,data_t,ncomp>);
        	m_EBIrregDataPtr = (EBIrregStruct<cent,data_t,ncomp>*)(allocator.alloc(memsize));
        	protoMemcpy(MEMTYPE_DEFAULT,m_EBIrregDataPtr, hostvec.data(), memsize, protoMemcpyHostToDevice);
        	// &(stack_deleter_boxdata) calls allocation.free
       		m_EBIrregDataSharedPtr = ::std::shared_ptr<EBIrregStruct<cent,data_t,ncomp>>(m_EBIrregDataPtr, &(stack_deleter_boxdata<MEMTYPE_DEFAULT>));
	}

	if(m_indicies->size() == 0)
	{
		m_indicies_data = nullptr;
		m_indicies_shared_data = ::std::shared_ptr<index_t>(nullptr);
	}
	else
	{
        	size_t memsize = sizeof(index_t) * m_indicies->size(); 
        	m_indicies_data = (index_t*)(allocator.alloc(memsize));
        	protoMemcpy(MEMTYPE_DEFAULT,m_indicies_data,m_indicies->data(),memsize,protoMemcpyHostToDevice);
        	m_indicies_shared_data = ::std::shared_ptr<index_t>(m_indicies_data, &(stack_deleter_boxdata<MEMTYPE_DEFAULT>));
	}
      }
      else 
      {
        PR_TIME("defineEBIrregData::ProtoMalloc");
        size_t memsize = hostvec.size()*sizeof(EBIrregStruct<cent,data_t,ncomp>);
        protoMalloc(MEMTYPE_DEFAULT,m_EBIrregDataPtr,memsize);
        protoMemcpy(MEMTYPE_DEFAULT,m_EBIrregDataPtr, hostvec.data(), memsize, protoMemcpyHostToDevice);
        m_EBIrregDataSharedPtr = ::std::shared_ptr<EBIrregStruct<cent,data_t,ncomp>>(m_EBIrregDataPtr, [](EBIrregStruct<cent,data_t,ncomp>* p){protoFree(MEMTYPE_DEFAULT,p);});

        memsize = sizeof(index_t) * m_indicies->size(); 
        protoMalloc(MEMTYPE_DEFAULT,m_indicies_data,memsize);
        protoMemcpy(MEMTYPE_DEFAULT,m_indicies_data,m_indicies->data(),memsize,protoMemcpyHostToDevice);
        m_indicies_shared_data = ::std::shared_ptr<index_t>(m_indicies_data, [](EBIndex<cent>* p){protoFree(MEMTYPE_DEFAULT,p);});
     } 
  }

    inline EBIrregStruct<cent,data_t, ncomp>* getEBIrregDataPtr()
    {
      return m_EBIrregDataSharedPtr.get();
    }

    inline bool hasIndex(const EBIndex<cent>& a_ind) const
    {
      bool found = false;
      //this should probably be a smarter search.   maybe sort them?
      for(int ivec = 0;(ivec < m_indicies->size() && !found); ivec++)
      {
        if((*m_indicies)[ivec] == a_ind)
        {
          found = true;
        }
      }

      return found;
    }
  
    std::size_t charsize(const Box&   a_bx,
                         unsigned int a_startcomp,
                         unsigned int a_numcomps) const
    {
      PR_TIME("IrregData::charsize");
      size_t retval = 0;
      retval += sizeof(data_t); //size is stored in a data_t type to avoid misaligned access on GPU
      retval += m_indicies->size()*(sizeof(index_t) + ncomp*sizeof(data_t));

      size_t align = sizeof(data_t);
      if(retval % align != 0) retval = (retval/ align + 1) * align;

      return retval;
    }


    inline shared_ptr<index_t> getDataIndicies()
    {
      return m_indicies_shared_data;
    }

    inline shared_ptr<index_t> getDataIndicies() const
    {
      return m_indicies_shared_data;
    }

    inline
    void copy(const IrregData<cent, data_t, ncomp>& a_src,
              const Box& a_bxFrom, unsigned int a_isrc,
              const Box& a_bxTo  , unsigned int a_idst,
              unsigned int a_inco)
    {
      PR_TIME("IrregData::copy");
      Box interBox = a_bxFrom & a_bxTo;
      auto indiciesSrc = a_src.getIndicies();
      auto indiciesDataSrc = a_src.getDataIndicies();
      unsigned int size = m_indicies->size(); 

      if(size == 0 || indiciesSrc->size()==0) return;

      int numSMs=1;
#ifdef PROTO_CUDA
      int devId = 0;
#ifdef PROTO_HIP
      hipDeviceGetAttribute(&numSMs, hipDeviceAttributeMultiprocessorCount, devId);
#else
      cudaDeviceGetAttribute(&numSMs, cudaDevAttrMultiProcessorCount, devId);
#endif
#endif
      //unsigned int threads = 256;
      //unsigned int blocks = (threads + size - 1)/threads;
      //protoLaunchKernelMemAsyncT<MEMTYPE_DEFAULT,EBBoxDataIndexer<cent,data_t>>(blocks, threads, 0, protoGetCurrentStream, 
      protoLaunchKernelMemAsyncT<MEMTYPE_DEFAULT,EBBoxDataIndexer<cent,data_t>>(32*numSMs, 256, 0, protoGetCurrentStream, 
			size, 
                        indiciesDataSrc.get(), indiciesSrc->size(),
			m_indicies_shared_data.get(), m_indicies->size(),
			m_rawPtr, 
			a_src.data(), 
			a_idst, a_isrc,
			a_inco,
			interBox.low(), interBox.high());
    }

    inline
    void linearOut(void*               a_buf,
                   const Box       &   a_bx,
                   unsigned int        a_startcomp,
                   unsigned int        a_numcomps) const
    {
      PR_TIME("IrregData::linearOut");
      unsigned int vecsize = m_indicies->size();
      const unsigned int threads = 256;
      const unsigned int blocks  = (vecsize + threads - 1) / threads; 
      assert(sizeof(data_t) >= sizeof(float));

      if(vecsize == 0) return;

      assert(vecsize == m_indicies->size());

      protoLaunchKernelMemAsyncT<MEMTYPE_DEFAULT,EBLinearOut<cent,data_t,ncomp>>(blocks, threads, 0, protoGetCurrentStream, 
                        vecsize,
                        m_indicies_shared_data.get(), m_indicies->size(),
                        a_buf,
                        data());

    }

    inline
    void  linearIn(void*                  a_buf,
                   const Box       &      a_bx,
                   unsigned int           a_startcomp,
                   unsigned int           a_numcomps)
    {
      PR_TIME("IrregData::linearIn");
      unsigned int vecsize = m_indicies->size();
      assert(sizeof(data_t) >= sizeof(float));
      if(vecsize == 0) return;

      int numSMs=1;
#ifdef PROTO_CUDA
      int devId = 0;
#ifdef PROTO_HIP
      hipDeviceGetAttribute(&numSMs, hipDeviceAttributeMultiprocessorCount, devId);
#else
      cudaDeviceGetAttribute(&numSMs, cudaDevAttrMultiProcessorCount, devId);
#endif
#endif

      protoLaunchKernelMemAsyncT<MEMTYPE_DEFAULT,EBLinearIn<cent,data_t,ncomp>>(32*numSMs, 256, 0, protoGetCurrentStream, 
                        m_indicies_shared_data.get(), vecsize,
                        a_buf,
                        m_data.get(),
                        a_bx.low(), a_bx.high());

    }
    ///
    static int preAllocatable()
    {
      return 2;
    }

    static int memTypeAllocation()
    {
      return 1;
    }

protected:
    inline void define( const shared_ptr< vector< index_t > >& a_indicies, bool a_stackAllocation = false);

    index_t*                             m_indicies_data = nullptr;
    shared_ptr<index_t>           m_indicies_shared_data = nullptr;
    bool m_boxdefined;
    Box m_inputBox; //exists if defined as part of a boxdata
    shared_ptr<vector<  index_t > >           m_indicies = nullptr;
    shared_ptr<data_t>                            m_data = nullptr;
    data_t*                                     m_rawPtr = nullptr;
    EBIrregStruct<cent,data_t,ncomp>*        m_EBIrregDataPtr = nullptr;
    shared_ptr<EBIrregStruct<cent,data_t,ncomp>> m_EBIrregDataSharedPtr = nullptr;
    IrregData(const IrregData<cent, data_t,ncomp>& a_src);
    IrregData& operator=(const IrregData<cent, data_t,ncomp>& a_src);
  };

}
#include "Proto_IrregDataImplem.H"
#include "Proto_IrregData.cpp"

#endif
