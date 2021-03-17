
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
//      for(unsigned int i = 0; i < a_size; i++)
//      {        
//        if(a_vof == a_indicies[i])
//        {
//          ivec = i;
//          found = true;
//          break;
//        }
//      }
      int begin = 0;
      int end = a_size;
      int i = (begin + a_size)/2;

      while(found == false && begin <= end)
      {        
        if(a_vof == a_indicies[i])
        {
          ivec = i;
          found = true;
          break;
        }
        if(a_vof <  a_indicies[i] )
          end = i;
	else
	  begin = i;
        i = (begin + end)/2;
      }

      if(!found)
      {
        PROTO_ASSERT(found, "vof not found");
      }

      ivec += a_ivar* a_size;

      return ivec;
    }

    CUDA_DECORATION
    inline bool contains(const Point& pt, Point& a_low, Point& a_high);

    template<CENTERING cent>
    CUDA_DECORATION
    inline bool contains(EBIndex<cent>& a_in, Point& a_low, Point& a_high);

    template <CENTERING cent>
    CUDA_DECORATION
    bool copyContains(const EBIndex<cent>& a_ind, Point & a_srcbox_low, Point & a_srcbox_high)
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

#ifdef PROTO_CUDA
    template<class data_t>
    __global__
    void printIrreg(const data_t* ptr, unsigned int size)
    {
      unsigned int idx = blockIdx.x * blockDim.x + threadIdx.x;

      if(idx < size)
	printf(" irregData[%d] = %f \n",idx,ptr[idx]);	
    }


    template <CENTERING cent, class data_t>
    __global__
    void EBBoxDataIndexer( unsigned int a_end, 
		EBIndex<cent>* a_indiciesSrc, unsigned int a_indiciesSizeSrc, 
		EBIndex<cent>* a_indiciesDst, unsigned int a_indiciesSizeDst,
		data_t* a_dst, const data_t* a_src, 
                unsigned int a_idst, unsigned int a_isrc,
		unsigned int a_inco,
		Point a_box_inf, Point a_box_high)
    {
      unsigned int idx = blockIdx.x * blockDim.x + threadIdx.x;

      if(idx < a_end)
      {
        EBIndex<cent>& vof = a_indiciesDst[idx];
        if( contains(vof,a_box_inf,a_box_high) && copyContains(vof,a_box_inf,a_box_high))
        {
#pragma unroll
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

   
    template <CENTERING cent, class data_t, unsigned int ncomp>
    __global__
    void EBLinearIn(EBIndex<cent>* a_indicies, unsigned int a_indiciesSize,
                void* a_buf, 
		data_t* a_dst,
                Point a_box_inf, Point a_box_high)
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
        {
#pragma unroll
          for(int icomp = 0; icomp < ncomp; icomp++)
          {
            size_t indxDst = index(a_indicies, a_indiciesSize, pt, icomp);
           // size_t indxSrc = index(indicies, size, pt, icomp);
           // a_dst[indxDst+icomp*a_indiciesSize] = data[indxSrc+icomp*size];
            a_dst[indxDst+icomp*a_indiciesSize] = data[idx+icomp*size];
          }
        }
      }
   }

    template <CENTERING cent, class data_t, unsigned int ncomp>
    __global__
    void EBLinearOut( unsigned int a_end,
                EBIndex<cent>* a_indicies, std::size_t a_indiciesSize,
                void* a_buf, 
	        const data_t* a_src)
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


#else

    template <CENTERING cent, class data_t>
    void EBBoxDataIndexer( unsigned int a_end, 
		EBIndex<cent>* a_indiciesSrc, unsigned int a_indiciesSizeSrc, 
		EBIndex<cent>* a_indiciesDst, unsigned int a_indiciesSizeDst,
		data_t* a_dst, const data_t* a_src, 
                unsigned int a_idst, unsigned int a_isrc,
                unsigned int a_inco,
		Point a_box_inf, Point a_box_high)

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

   
    template <CENTERING cent, class data_t, unsigned int ncomp>
    void EBLinearIn(EBIndex<cent>* a_indicies, unsigned int a_indiciesSize,
                void* a_buf, 
		data_t* a_dst,
                Point a_box_inf, Point a_box_high)
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

    template <CENTERING cent, class data_t, unsigned int ncomp>
    void EBLinearOut( unsigned int a_end,
                EBIndex<cent>* a_indicies, std::size_t a_indiciesSize,
                void* a_buf, 
	        const data_t* a_src)
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

    template<class data_t>
    void printIrreg(const data_t* ptr, unsigned int size)
    {
      for(int idx = 0 ; idx < size; idx++)
	printf(" irregData[%d] = %f \n",idx,ptr[idx]);	
    }

#endif



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
      protoMalloc(m_rawPtr, size*sizeof(data_t));
      protoMemcpy(m_rawPtr, a_data, size*sizeof(data_t), protoMemcpyHostToDevice);
      m_data = ::std::shared_ptr<data_t>(m_rawPtr, [](data_t* p){protoFree(p);});
      m_indicies = ::std::shared_ptr<std::vector<index_t>>(new std::vector<index_t>[1], [](std::vector<index_t>* p) { delete[] p;});
      m_indicies->resize(size);
      for(int i = 0; i < size ; i++)
        (*m_indicies)[i] = a_indicies[i];

      defineEBIrregData();
    }

    ///
    virtual ~IrregData()
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
        protoLaunchKernelMemAsync((opKernel<data_t,BoxDataOp::Add>), nBlocks, blockSize, 0, protoGetCurrentStream, m_rawPtr, a_input.m_rawPtr, nsize);
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
        protoLaunchKernelMemAsync((opKernel<data_t,BoxDataOp::Subtract>), nBlocks, blockSize, 0, protoGetCurrentStream, m_rawPtr, a_input.m_rawPtr, nsize);
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
        protoLaunchKernelMemAsync((opKernel<data_t,BoxDataOp::Multiply>), nBlocks, blockSize, 0, protoGetCurrentStream, m_rawPtr, a_input.m_rawPtr, nsize);
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
      m_EBIrregDataSharedPtr  = dynamic_pointer_cast<EBIrregStruct<cent,data_t,ncomp>>(a_input.getEBIrregDataSharedPtr());
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

    inline void define(const Box& a_box,  const EBGraph& a_graph)
    {
      this->m_boxdefined = true;
      this->m_inputBox   = a_box;
      Box interBox = a_box & a_graph.getDomain();
      shared_ptr<vector<EBIndex<cent> > > indicies(new vector< EBIndex<cent> >() );
      *indicies = getIrregIndicies<cent>(interBox, a_graph);

      define(indicies);
      defineEBIrregData();
    }

    inline void defineEBIrregData()
    {
      PR_TIME("defineEBIrregData");
      vector< EBIrregStruct<cent, data_t, ncomp> > hostvec = getEBIrregStruct(*m_indicies, *this);
      size_t memsize = hostvec.size()*sizeof(EBIrregStruct<cent,data_t,ncomp>);
      protoMalloc(m_EBIrregDataPtr,memsize);
      protoMemcpy(m_EBIrregDataPtr, hostvec.data(), memsize, protoMemcpyHostToDevice);
      m_EBIrregDataSharedPtr = ::std::shared_ptr<EBIrregStruct<cent,data_t,ncomp>>(m_EBIrregDataPtr, [](EBIrregStruct<cent,data_t,ncomp>* p){protoFree(p);});

      memsize = sizeof(index_t) * m_indicies->size(); 
    
      protoMalloc(m_indicies_data,memsize);
      protoMemcpy(m_indicies_data,m_indicies->data(),memsize,protoMemcpyHostToDevice);
      m_indicies_shared_data = ::std::shared_ptr<index_t>(m_indicies_data, [](EBIndex<cent>* p){protoFree(p);});

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
//      auto indicies = *m_indicies; //linearCharBufSize(m_inputBox);
      size_t retval = 0;
      retval += sizeof(data_t); //size is stored in a data_t type to avoid misaligned access on GPU
      retval += m_indicies->size()*(sizeof(index_t) + ncomp*sizeof(data_t));

      size_t align = sizeof(data_t);
      if(retval % align != 0) retval = (retval/ sizeof(align) + 1) * sizeof(align);
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
      unsigned int size = indiciesSrc->size(); 
      unsigned int threads = 256;
      unsigned int blocks = (threads + size - 1)/threads;

      protoLaunchKernel((EBBoxDataIndexer<cent,data_t>), blocks, threads, 
			size, 
                        indiciesDataSrc.get(), indiciesSrc->size(),
			m_indicies_shared_data.get(), m_indicies->size(),
			m_rawPtr, 
			a_src.data(), 
			a_idst, a_isrc,
			a_inco,
			interBox.low(), interBox.high());
    }

    // debug tool
    inline void printIrr() const
    {
	protoLaunchKernel((printIrreg<data_t>), size(), 1, data(),size());
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
//      printIrr();

      protoLaunchKernel((EBLinearOut<cent,data_t,ncomp>), blocks, threads,
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

      int numSMs;
#ifdef PROTO_CUDA
      int devId = 0;
      cudaDeviceGetAttribute(&numSMs, cudaDevAttrMultiProcessorCount, devId);
#endif

      protoLaunchKernel((EBLinearIn<cent,data_t,ncomp>), 32*numSMs, 256,
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
    inline void define( const shared_ptr< vector< index_t > >& a_indicies);

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
