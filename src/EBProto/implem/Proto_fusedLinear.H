#pragma once
#include "Chombo_Copier.H"
#include "Chombo_ProtoInterface.H"
#include "Proto_EBBoxData.H"

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
	EBBoxDataIndexer<cent,data_t> kernel;
    	for(int it = 0; it < numberOfIteration; it++)
		kernel.gpu(a_end[it], 
			a_indiciesSrc[it], a_indiciesSizeSrc[it], 
			a_indiciesDst[it], a_indiciesSizeDst[it], 
			a_dst[it], a_src[it], 
			a_idst[it], a_isrc[it], 
			a_inco[it], 
			a_box_inf[it], a_box_high[it]);
}

template < CENTERING cent, class data_t, unsigned int ncomp >
class EBBoxDataIndexerFused  
{
public:
	EBBoxDataIndexerFused(){}

	~EBBoxDataIndexerFused()
	{
		free();
	}

	void allocateCopy(unsigned int a_numberOfIteration)
	{
		if(m_isAllocated)
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
		m_indiciesDst = m_indiciesSrc + a_numberOfIteration;

		bytes = sizeof(data_t*);
		nPtr  = 2;
		size  = bytes*nPtr*a_numberOfIteration;
		protoMalloc(MEMTYPE_DEFAULT, m_src, size);
		m_dst = (data_t**)(m_src) + a_numberOfIteration;

		bytes = sizeof(Proto::Point*);
		nPtr  = 2;
		size  = bytes*nPtr*a_numberOfIteration;
		protoMalloc(MEMTYPE_DEFAULT, m_box_inf, size);
		m_box_high = m_box_inf + a_numberOfIteration;

		m_isAllocated = true;
	}

/*
	//void fillparams(CopyIterator a_copier, BoxLayoutData<T>& a_dst, const Interval& a_srcComps)
	void fillparams(Chombo4::CopyIterator a_copier, BoxLayoutData<EBBoxData<cent,data_t,ncomp>> a_dst, const BoxLayoutData<EBBoxData<cent,data_t,ncomp>> a_src, unsigned int a_dstComps_begin, unsigned int a_srcComps_begin, unsigned int a_srcComps_size)
	{
		if(m_isAllocated == true) return;
	
		int items=a_copier.size();
		m_numberOfIteration = items;
		allocateCopy(m_numberOfIteration);	
  		for (int n=0; n<items; n++)
  		{	
    			const Chombo4::MotionItem& item = a_copier[n];
      			::Proto::Box bxFr = ::ProtoCh::getProtoBox(item.fromRegion);
      			::Proto::Box bxTo = ::ProtoCh::getProtoBox(item.toRegion);
			::Proto::Box inter = bxFr & bxTo;
    			fillparam(n, 
				a_dst[item.toIndex].getIrregData(), 
            			a_src[item.fromIndex].getIrregData(),
           			a_dstComps_begin,
				a_srcComps_begin,a_srcComps_size,
				inter
			);
  		}
	}
*/
	void fillparam(
		unsigned int a_it, 
		IrregData<cent, data_t, ncomp>& a_dst, 
		const IrregData<cent, data_t, ncomp>& a_src, 
		unsigned int a_idst, unsigned int a_isrc, 
		unsigned int a_inco,
		Proto::Box a_box
	)
	{
		
		fillone
		(
			a_it,
			a_src.getIndicies()->size(), 				// a_end
			a_src.getDataIndicies().get(), a_src.getIndicies()->size(), 	// ebindex<cent>** a_indiciessrc, unsigned int* a_indiciessizesrc,
			a_dst.getDataIndicies().get(), a_dst.getIndicies()->size(), 	// ebindex<cent>** a_indiciesdst, unsigned int* a_indiciessizedst,
			a_src.data(), a_dst.data(),				// data_t** a_dst, const data_t** a_src,		
			a_idst, a_isrc,						// unsigned int* a_idst, unsigned int* a_isrc,
                	a_inco, 						// unsigned int* a_inco,
                	a_box.low(), a_box.high()
		);
	}


	void fillone(   unsigned int it, unsigned int end, 
			EBIndex<cent>* indiciesSrc, unsigned int indiciesSizeSrc,
			EBIndex<cent>* indiciesDst, unsigned int indiciesSizeDst,
			const data_t* src, data_t* dst,
			unsigned int idst, unsigned int isrc, unsigned int inco,
			Proto::Point box_inf, Proto::Point box_high
		)
	{
		m_end[it] 		= end;
		m_indiciesSrc[it] 	= indiciesSrc;
		m_indiciesDst[it] 	= indiciesDst;
		m_indiciesSizeSrc[it] 	= indiciesSizeSrc;
		m_indiciesSizeDst[it] 	= indiciesSizeDst;
		//m_src[it] 		= const_cast<data_t*>(src);
		m_src[it] 		= src;
		m_dst[it] 		= dst;
		m_idst[it] 		= idst;
		m_isrc[it] 		= isrc;
		m_inco[it] 		= inco;
		m_box_inf[it] 		= box_inf;
		m_box_high[it] 		= box_high;
	}

	void launch()
	{
		const unsigned int threads = 256;
		const unsigned int nBlocks = (m_numberOfIteration + threads - 1) / threads;
//		protoLaunchKernelT(MEMTYPE_DEFAULT,
//			(kernel_ebcopy<cent,data_t>), threads, nBlocks,
			kernel_ebcopy<cent,data_t><<<threads, nBlocks>>>(
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

	void free()
	{
		protoFreeGPU(m_end);
		protoFreeGPU(m_indiciesSrc);
		protoFreeGPU(m_indiciesDst);
		protoFreeGPU(m_indiciesSizeSrc);
		protoFreeGPU(m_indiciesSizeDst);
		protoFreeGPU(m_src);
		protoFreeGPU(m_dst);
		protoFreeGPU(m_dst);
		protoFreeGPU(m_idst);
		protoFreeGPU(m_isrc);
		protoFreeGPU(m_inco);
		protoFreeGPU(m_box_inf);
		protoFreeGPU(m_box_high);
	}


public:
	unsigned int m_numberOfIteration 	= 0;
	bool m_isAllocated 			= false;
	unsigned int* m_end 			= nullptr;
	EBIndex<cent>** m_indiciesSrc 		= nullptr;
	EBIndex<cent>** m_indiciesDst 		= nullptr;
	unsigned int* m_indiciesSizeSrc		= nullptr;
	unsigned int* m_indiciesSizeDst 	= nullptr;
	const data_t** m_src 			= nullptr;
	data_t** m_dst 				= nullptr;
	unsigned int* m_idst 			= nullptr;
	unsigned int* m_isrc 			= nullptr;
	unsigned int* m_inco 			= nullptr;
	Proto::Point* m_box_inf 		= nullptr;
	Proto::Point* m_box_high 		= nullptr;
};


template <CENTERING cent, class data_t, unsigned int ncomp>
__global__
void kernel_eblinearout
(
	unsigned int a_numberofiteration,
	unsigned int* a_end,
	EBIndex<cent>** a_indicies, std::size_t* a_indiciessize,
	void** a_buf,
	const data_t** a_src
)
{
	Proto::EBLinearOut<cent, data_t, ncomp> kernel;
	for(unsigned int it = 0 ; it < a_numberofiteration ; it++)
		kernel.gpu(a_end[it], a_indicies[it], a_indiciessize[it], a_buf[it], a_src[it]);
}


template <CENTERING cent, class data_t, unsigned int ncomp>
class EBLinearOutFused
{
	EBLinearOutFused() {}

	void allocate(unsigned int a_numberOfIteration)
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

		bytes = sizeof(const data_t*);
		nPtr  = 1;
		size  = bytes*nPtr*a_numberOfIteration;
		protoMalloc(MEMTYPE_DEFAULT, m_src, size);

		bytes = sizeof(void*);
		nPtr  = 1;
		size  = bytes*nPtr*a_numberOfIteration;
		protoMalloc(MEMTYPE_DEFAULT, m_dst, size);

		isAllocated = true;
	}

	template<typename T>
        void fillparams( Chombo4::CopierBuffer* a_buffer, T& a_src, unsigned int a_srcComps_begin, unsigned int a_srcComps_size)
        {
		int m_numberofiteration = a_buffer->m_fromMe.size();

		for (unsigned int i=0; i< m_numberofiteration; ++i)
		{
			const Chombo4::CopierBuffer::bufEntry& entry = a_buffer->m_fromMe[i];
			char* charbuf = (char*)entry.bufPtr; // 
			auto regData = a_src[entry.item->fromIndex].getRegData();
			Proto::Box bx = ProtoCh::getProtoBox(entry.item->fromRegion);
			unsigned int regsize = regData.charsize(bx, a_srcComps_begin, a_srcComps_size);
    			charbuf += regsize;
			fillparam(i, charbuf, a_src[entry.item->fromIndex].getIrregData());
		}

	}

        void fillparam(
                unsigned int a_it,
                void* a_dst,
                IrregData<cent, data_t, ncomp>& a_src
        )
        {

                fillone
                (
                        a_it,
                        a_src.getindicies().size(),                             // a_end
                        a_src.getindicies().get(), a_src.getindicies().size(),  // ebindex<cent>** a_indiciessrc, unsigned int* a_indiciessizesrc,
                        a_dst, a_src.data()                             	// void** a_dst, const data_t** a_src,
                );
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
		const unsigned int nBlocks = (m_numberOfIteration + threads -1)/threads;
		protoLaunchKernel(kernel_eblinearout, nBlocks, threads,
					m_end, 
					m_indiciesSrc,
					m_indiciesSizeSrc,
					m_dst, m_src
		);
	}

public:

        unsigned int m_numberOfIteration 	= 0;
        bool isAllocated 			= false;
        unsigned int* m_end 			= nullptr;
        EBIndex<cent>** m_indiciesSrc 		= nullptr;
        unsigned int* m_indiciesSizeSrc 	= nullptr;
        const data_t** m_src 			= nullptr;
        void** m_dst 				= nullptr;
};


template <CENTERING cent, class data_t, unsigned int ncomp>
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
	EBLinearIn<cent,data_t,ncomp> kernel;
	for(unsigned int it = 0 ; it < a_numberofiteration ; it++)
		kernel.gpu(a_indicies[it], a_indiciesSize[it], a_buf[it], a_dst[it],a_box_inf[it], a_box_high[it]);
}

template <CENTERING cent, class data_t, unsigned int ncomp>
class EBLinearInFused
{
	EBLinearInFused() {}


	void allocate(unsigned int a_numberOfIteration)
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
		protoMalloc(MEMTYPE_DEFAULT, m_indiciesSizeDst, size);

		bytes = sizeof(EBIndex<cent>*);
		nPtr  = 1;
		size  = bytes*nPtr*a_numberOfIteration;
		protoMalloc(MEMTYPE_DEFAULT, m_indiciesDst, size);

		bytes = sizeof(void*);
		nPtr  = 1;
		size  = bytes*nPtr*a_numberOfIteration;
		protoMalloc(MEMTYPE_DEFAULT, m_dst, size);

		bytes = sizeof(data_t*);
		nPtr  = 1;
		size  = bytes*nPtr*a_numberOfIteration;
		protoMalloc(MEMTYPE_DEFAULT, m_dst, size);

		bytes = sizeof(Proto::Point*);
		nPtr  = 2;
		size  = bytes*nPtr*a_numberOfIteration;
		protoMalloc(MEMTYPE_DEFAULT, m_box_inf, size);
		m_box_high = m_box_inf + m_numberOfIteration;

		isAllocated = true;
	}

//later
	template<typename T>
        void fillparams( Chombo4::CopierBuffer* a_buffer, T& a_dst, unsigned int a_destComps_begin, unsigned int a_destComps_size)
        {
		int m_numberofiteration = a_buffer->m_toMe.size();

		for (unsigned int i=0; i< m_numberofiteration; ++i)
		{
			const Chombo4::CopierBuffer::bufEntry& entry = a_buffer->m_toMe[i];
			char* charbuf = (char*)entry.bufPtr; // 
			auto regData = a_dst[entry.item->fromIndex].getRegData();
			::Proto::Box bx = ProtoCh::getProtoBox(entry.item->toRegion);
			unsigned int regsize = regData.charsize(bx, a_destComps_begin, a_destComps_size);
    			charbuf += regsize;
			fillparam(i, charbuf, a_dst[entry.item->toIndex].getIrregData(), bx);
		}

	}

        void fillparam(
                unsigned int a_it,
                void* a_dst,
                IrregData<cent, data_t, ncomp>& a_src,
		::Proto::Box a_bx
        )
        {
                fillone
                (
                        a_it,
                        a_src.getindicies().get(), a_src.getindicies().size(),  // ebindex<cent>** a_indiciessrc, unsigned int* a_indiciessizesrc,
                        a_dst, a_src.data(),                             	// void** a_dst, const data_t** a_src,
			a_bx.low(), a_bx.high()
                );
        }


        void fillOne(   unsigned int it, 
                        EBIndex<cent>* indiciesDst, 
                        unsigned int indiciesSizeDst,
                        void* buf, data_t* dst, 
			Point box_inf, Point box_high
                )
        {
                m_indiciesDst[it]       = indiciesDst;
                m_indiciesSizeDst[it]   = indiciesSizeDst;
                m_buf[it]               = buf;
                m_dst[it]               = dst;
                m_box_inf[it]           = box_inf;
                m_box_high[it]          = box_high;
        }

	void launch()
	{
		const unsigned int threads = 256;
		const unsigned int nBlocks = (m_numberOfIteration + threads -1)/threads;
		protoLaunchKernel(kernel_eblinearout, nBlocks, threads,
					m_numberOfIteration,
					m_indiciesDst,
					m_indiciesSizeDst,
					m_buf, m_dst,
					m_box_inf,
					m_box_high
		);
	}

public:
        unsigned int m_numberOfIteration = 0;
        bool isAllocated 		= false;
        unsigned int* m_end 		= nullptr;
        EBIndex<cent>** m_indiciesDst 	= nullptr;
        unsigned int* m_indiciesSizeDst = nullptr;
        void** m_buf 			= nullptr;
        data_t** m_dst 			= nullptr;
	Proto::Point *m_box_inf 	= nullptr;
	Proto::Point *m_box_high	= nullptr;

};
