#pragma once


template <CENTERING cent, class data_t>
__global__
void kernel_ebcopy
(
                unsigned int numberOfIteration,
                unsigned int* a_end,
                EBIndex<cent>** a_indiciesSrc, unsigned int* a_indiciesSizeSrc,
                EBIndex<cent>** a_indiciesDst, unsigned int* a_indiciesSizeDst,
                data_t** a_dst, const data_t** a_src,
                unsigned int* a_idst, unsigned int* a_isrc,
                unsigned int* a_inco,
                Point* a_box_inf, Point* a_box_high
)
{
    unsigned int idx = blockIdx.x * blockDim.x + threadIdx.x;

    for(int it = 0; it < numberOfIteration, it++)
      if(idx < a_end[it])
      {
        EBIndex<cent>& vof = a_indiciesDst[it][idx];
        if( contains(vof,a_box_inf,a_box_high) && copyContains(vof,a_box_inf,a_box_high))
        {
#pragma unroll
          for(unsigned int ico = 0; ico < a_inco; ico++)
          {
            unsigned int isrc = a_isrc + ico;
            unsigned int idst = a_idst + ico;
            size_t indDst = index(a_indiciesDst[it], a_indiciesSizeDst[it], vof, idst[it]);
            size_t indSrc = index(a_indiciesSrc[it], a_indiciesSizeSrc[it], vof, isrc[it]);
            a_dst[it][indDst] = a_src[it][indSrc];
          }
        }
      }

}




template <CENTERING cent, class data_t>
class EBBoxDataIndexerFused
{
public:
	EBBoxDataIndexerFused(){}


	void allocateCopy(unsigned int a_numberOfIteration)
	{
		if(isAllocated)
		{
			assert(a_numberOfIteration == m_numberOfIteration);
			return;
		}

		m_numberOfIteration = a_numberOfIteration;
		unsigned int bytes = sizeof(unsigned int);
		unsigned int nPtr  = 6;
		unsigned int size  = bytes*nPtr*a_numberOfIteration;
		protoMalloc(MEMTYPE_DEFAULT, m_end, size);
		m_indiciesSizeSrc = m_end + a_numberOfIteration;	
		m_indiciesSizeDst = m_end + 2*a_numberOfIteration;	
		m_idst            = m_end + 3*a_numberOfIteration;
		m_isrc            = m_end + 4*a_numberOfIteration;
		m_inco            = m_end + 5*a_numberOfIteration;

		bytes = sizeof(EBIndex<cent>*);
		nPtr  = 2;
		size  = bytes*nPtr*a_numberOfIteration;
		protoMalloc(MEMTYPE_DEFAULT, m_indiciesSrc, size);
		m_indiciesSizeDst = m_indiciesSrc + a_numberOfIteration;

		bytes = sizeof(data_t*)
		nPtr  = 2;
		size  = bytes*nPtr*a_numberOfIteration;
		protoMalloc(MEMTYPE_DEFAULT, m_src, size);
		m_dst = m_src + a_numberOfIteration;

		bytes = sizeof(Proto::Point*)
		nPtr  = 2;
		size  = bytes*nPtr*a_numberOfIteration;
		protoMalloc(MEMTYPE_DEFAULT, m_box_inf, size);
		m_box_high = m_box_inf + a_numberOfIteration;

		isAllocated = true;
	}

	void fillparams( copyiterator a_copier, BoxLayoutData<T>& a_dst, const Interval& a_srcComps)
	{

		int items=a_copier.size();
		assert(m_numberofiteration==items);
  		for (int n=0; n<items; n++)
  		{	
    			const motionitem& item = it[n];

    			fillparam(n, 
				a_dest[item.toindex].getIrregData(), 
            			a_src[item.fromindex].getIrregData(),
           			item.toregion.begin(),
				a_srcComps.begin(),a_srcComps.size());
  		}
	}

	void fillparam(
		unsigned int a_it, 
		irregdata<cent, data_t, ncomp>& a_dst, 
		irregdata<cent, data_t, ncomp>& a_src, 
		unsigned int a_idst, unsigned int a_src, 
		unsigned int a_inco
	)
	{
		
		fillone
		(
			a_it,
			a_src.getindicies().size(), 				// a_end
			a_src.getindicies().get(), a_src.getindicies().size(), 	// ebindex<cent>** a_indiciessrc, unsigned int* a_indiciessizesrc,
			a_dst.getindicies().get(), a_dst.getindicies().size(), 	// ebindex<cent>** a_indiciesdst, unsigned int* a_indiciessizedst,
			a_dst.data(), a_src.data(),				// data_t** a_dst, const data_t** a_src,		
			a_idst, a_isrc,						// unsigned int* a_idst, unsigned int* a_isrc,
                	a_inco, 						// unsigned int* a_inco,
                	a_box.low(), a_box.high()
		)
	}


	void fillone(   unsigned int it, unsigned int end, 
			ebindex<cent>* indiciessrc, ebindex<cent>* indiciesdst,
			unsigned int indiciessizesrc, unsigned int indiciessizedst,
			data_t* src, data_t* dst,
			unsigned int idst, unsigned int isrc, unsigned int inco,
			proto::point box_inf, proto::point box_high
		)
	{
		m_end[it] 		= end;
		m_indiciessrc[it] 	= indiciessrc;
		m_indiciesdst[it] 	= indiciesdst;
		m_indiciessizesrc[it] 	= indiciessizesrc;
		m_indiciessizedst[it] 	= indiciessizedst;
		m_src[it] 		= src;
		m_dst[it] 		= dst;
		m_idst[it] 		= idst;
		m_isrc[it] 		= isrc;
		m_inco[it] 		= inco;
		m_box_inf[it] 		= box_inf;
		m_box_high[it] 		= m_box_high;
	}

	void launch()
	{
		const unsigned int threads = 256;
		const unsigned int nBlocks = (m_numberOfIteration + threads - 1) / threads;
		protoLaunchKernel(
			kernel_ebcopy, threads, nBlocks,
			m_numberOfIteration,
			m_end,
			m_indiciesSrc, m_indiciesSizeSrc,
			m_indiciesDst, m_indiciesSizeDst,
			m_dst, m_src,
			m_idst, m_isrc,
			m_inco,
			m_box_inf, m_box_high
		);
	}


public:
	m_numberofiteration = 0;
	bool isallocated = false;
	unsigned int* m_end = nullptr;
	ebindex<cent>** m_indiciessrc = nullptr;
	ebindex<cent>** m_indiciesdst = nullptr;
	unsigned int* m_indiciessizesrc = nullptr;
	unsigned int* m_indiciessizedst = nullptr;
	data_t** m_src = nullptr;
	data_t** m_dst = nullptr;
	unsigned int* m_idst = nullptr;
	unsigned int* m_isrc = nullptr;
	unsigned int* m_inco = nullptr;
	proto::point* m_box_inf = nullptr;
	proto::point* m_box_high = nullptr;
};


template <centering cent, class data_t, unsigned int ncomp>
__global__
void kernel_eblinearout
(
	unsigned int a_numberofiteration,
	unsigned int* a_end,
	ebindex<cent>** a_indicies, std::size_t* a_indiciessize,
	void** a_buf,
	const data_t** a_src
)
{
	EBLinearOut kernel;
	for(unsigned int it = 0 ; it < a_numberofiteration ; it++)
		kernel.gpu(a_end[it], a_indicies[it], a_indiciessize[it], a_buf[it], a_src[it]);
}


template <centering cent, class data_t, unsigned int ncomp>
class EBLinearOutFused
{
	EBLinearOutFused() {}


	void allocate(unsigned int a_numberOfIaeration)
	{
		if(isAllocated)
		{
			assert(a_numberOfIteration == m_numberOfIteration);
			return;
		}

		m_numberOfIteration = a_numberOfIteration;
		unsigned int bytes = sizeof(unsigned int);
		unsigned int nPtr  = 2;
		unsigned int size  = bytes*nPtr*a_numberOfIteration;
		protoMalloc(MEMTYPE_DEFAULT, m_end, size);
		m_indiciesSizeSrc = m_end + a_numberOfIteration;	

		bytes = sizeof(EBIndex<cent>*);
		nPtr  = 1;
		size  = bytes*nPtr*a_numberOfIteration;
		protoMalloc(MEMTYPE_DEFAULT, m_indiciesSrc, size);

		bytes = sizeof(const data_t*)
		nPtr  = 1;
		size  = bytes*nPtr*a_numberOfIteration;
		protoMalloc(MEMTYPE_DEFAULT, m_src, size);

		bytes = sizeof(void*)
		nPtr  = 1;
		size  = bytes*nPtr*a_numberOfIteration;
		protoMalloc(MEMTYPE_DEFAULT, m_dst, size);

		isAllocated = true;
	}

	template<typename T>
        void fillparams( CopierBuffer* a_buffer, const BoxLayoutData<T>& a_src, const Interval& a_srcComps)
        {
		int m_numberofiteration = a_buff->m_fromMe.size();

		for (unsigned int i=0; i< m_numberofiteration; ++i)
		{
			const CopierBuffer::bufEntry& entry = m_buff->m_fromMe[i];
			char* charbuf = (char*)entry.bufPtr; // 
			auto regDat = a_src[entry.item->fromIndex].getRegData();
			::Proto::Box bx = ProtoCh::getProtoBox(entry.item->fromRegion);
			unsigned int regsize = regData.charsize(bx, a_srcComps.begin(), a_srcComps.size());
    			charbuf += regsize;
			fillparam(i, charbuf, a_src[entry.item->fromIndex].getIrregData());
		}

	}

        void fillparam(
                unsigned int a_it,
                void* a_dst,
                irregdata<cent, data_t, ncomp>& a_src
        )
        {

                fillone
                (
                        a_it,
                        a_src.getindicies().size(),                             // a_end
                        a_src.getindicies().get(), a_src.getindicies().size(),  // ebindex<cent>** a_indiciessrc, unsigned int* a_indiciessizesrc,
                        a_dst, a_src.data(),                             	// void** a_dst, const data_t** a_src,
                )
        }


        void fillOne(   unsigned int it, unsigned int end,
                        EBIndex<cent>* indiciesSrc, 
                        unsigned int indiciesSizeSrc,
                        const data_t* src, void* dst
                )
        {
                m_end[it]               = end;
                m_indiciesSrc[it]       = indiciesSrc;
                m_indiciesSizeSrc[it]   = indiciesSizeSrc;
                m_src[it]               = src;
                m_dst[it]               = dst;
        }

	void launch()
	{
		const unsigned int threads = 256;
		const unsigned int nBlocks = (a_numberofiteration + threads -1)/threads;
		protoLaunchKernel(kernel_eblinearout, nBlocks, threads,
					m_end, 
					m_indicies,
					m_indiciessize,
					m_buf, m_src
		);
	}

public:
        m_numberOfIteration = 0;
        bool isAllocated = false;
        unsigned int* m_end = nullptr;
        EBIndex<cent>** m_indiciesSrc = nullptr;
        unsigned int* m_indiciesSizeSrc = nullptr;
        const data_t** m_src = nullptr;
        void** m_dst = nullptr;
};


template <centering cent, class data_t, unsigned int ncomp>
__global__
void kernel_eblinearin
(
	unsigned int a_numberofiteration,
        EBIndex<cent>** a_indicies, unsigned int* a_indiciesSize,
        void** a_buf,
        data_t** a_dst,
        Point* a_box_inf, Point* a_box_high
)
{
	EBLinearIn kernel;
	for(unsigned int it = 0 ; it < a_numberofiteration ; it++)
		kernel.gpu(a_end[it], a_indicies[it], a_indiciessize[it], a_buf[it], a_dst[it],a_box_inf[it], a_box_high[it]);
}

template <centering cent, class data_t, unsigned int ncomp>
class EBLinearInFused
{
	EBLinearOutFused() {}


	void allocate(unsigned int a_numberOfIaeration)
	{
		if(isAllocated)
		{
			assert(a_numberOfIteration == m_numberOfIteration);
			return;
		}

		m_numberOfIteration = a_numberOfIteration;
		unsigned int bytes = sizeof(unsigned int);
		unsigned int nPtr  = 1;
		unsigned int size  = bytes*nPtr*a_numberOfIteration;
		protoMalloc(MEMTYPE_DEFAULT, m_indiciesSize, size);

		bytes = sizeof(EBIndex<cent>*);
		nPtr  = 1;
		size  = bytes*nPtr*a_numberOfIteration;
		protoMalloc(MEMTYPE_DEFAULT, m_indiciesSrc, size);

		bytes = sizeof(void*)
		nPtr  = 1;
		size  = bytes*nPtr*a_numberOfIteration;
		protoMalloc(MEMTYPE_DEFAULT, m_dst, size);

		bytes = sizeof(data_t*)
		nPtr  = 1;
		size  = bytes*nPtr*a_numberOfIteration;
		protoMalloc(MEMTYPE_DEFAULT, m_dst, size);

		bytes = sizeof(Point*)
		nPtr  = 2;
		size  = bytes*nPtr*a_numberOfIteration;
		protoMalloc(MEMTYPE_DEFAULT, m_box_inf, size);
		a_box_high = m_box_inf + m_numberOfIteration;

		isAllocated = true;
	}

//later
	template<typename T>
        void fillparams( CopierBuffer* a_buffer, dst BoxLayoutData<T>& a_dst, const Interval& a_destComps)
        {
		int m_numberofiteration = a_buff->m_toMe.size();

		for (unsigned int i=0; i< m_numberofiteration; ++i)
		{
			const CopierBuffer::bufEntry& entry = m_buff->m_toMe[i];
			char* charbuf = (char*)entry.bufPtr; // 
			auto regDat = a_src[entry.item->fromIndex].getRegData();
			::Proto::Box bx = ProtoCh::getProtoBox(entry.item->ToRegion);
			unsigned int regsize = regData.charsize(bx, a_destComps.begin(), a_destComps.size());
    			charbuf += regsize;
			fillparam(i, charbuf, a_dst[entry.item->ToIndex].getIrregData(), bx);
		}

	}

        void fillparam(
                unsigned int a_it,
                void* a_dst,
                irregdata<cent, data_t, ncomp>& a_src,
		::Proto::Box a_bx
        )
        {

                fillone
                (
                        a_it,
                        a_src.getindicies().size(),                             // a_end
                        a_src.getindicies().get(), a_src.getindicies().size(),  // ebindex<cent>** a_indiciessrc, unsigned int* a_indiciessizesrc,
                        a_dst, a_src.data(),                             	// void** a_dst, const data_t** a_src,
			a_bx.low(), a_bx.high()
                )
        }


        void fillOne(   unsigned int it, unsigned int end,
                        EBIndex<cent>* indiciesSrc, 
                        unsigned int indiciesSizeSrc,
                        void* buf, data* dst, 
			Point box_inf, Point box_high
                )
        {
                m_end[it]               = end;
                m_indiciesSrc[it]       = indiciesSrc;
                m_indiciesSizeSrc[it]   = indiciesSizeSrc;
                m_buf[it]               = buf;
                m_dst[it]               = dst;
                m_box_inf[it]           = box_inf;
                m_box_high[it]          = box_high;
        }

	void launch()
	{
		const unsigned int threads = 256;
		const unsigned int nBlocks = (a_numberofiteration + threads -1)/threads;
		protoLaunchKernel(kernel_eblinearout, nBlocks, threads,
					m_numberOfIteration,
					m_indicies,
					m_indiciessize,
					m_buf, m_dst,
					m_box_inf,
					m_box_high
		);
	}

public:
        unsigned int m_numberOfIteration = 0;
        bool isAllocated 		= false;
        unsigned int* m_end 		= nullptr;
        EBIndex<cent>** m_indiciesSrc 	= nullptr;
        unsigned int* m_indiciesSizeSrc = nullptr;
        void** m_buf 			= nullptr;
        data_t** m_dst 			= nullptr;
	Proto::Point *m_box_inf 	= nullptr;
	Proto::Point *m_box_high	= nullptr;

};
