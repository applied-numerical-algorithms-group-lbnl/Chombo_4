#ifndef __Proto_AggStencilImplem__
#define __Proto_AggStencilImplem__

namespace Proto
{

  ///
  template <CENTERING srcCenter, CENTERING dstCenter, typename data_t>
  AggStencil<srcCenter, dstCenter, data_t>::
  AggStencil(const vector< EBIndex<     dstCenter        > > & a_dstVoFs,
             const vector< LocalStencil<srcCenter, data_t> > & a_stencil,
             const EBGraph                                   & a_srcGraph,
             const EBGraph                                   & a_dstGraph,
             const Box                                       & a_srcValid,
             const Box                                       & a_dstValid,
             const Point                                     & a_srcGhost,
             const Point                                     & a_dstGhost)
  {

    PR_TIME("AggSten.constructor");
    PROTO_ASSERT(a_dstVoFs.size() == a_stencil.size(), "mismatch in aggstencil construction");
    

//    //begin debugging hooks
//    m_dstVoFs = a_dstVoFs;
//    m_stencil = a_stencil;
    //end   debugging hooks
    if(a_dstVoFs.size() > 0)
    {
      Box srcbox = a_srcValid.grow(a_srcGhost);
      Box dstbox = a_dstValid.grow(a_dstGhost);

      m_srcBox = srcbox;
      m_dstBox = dstbox;

      m_flopsThisApply = 0;        

      bool useStack = true;
      EBBoxData<srcCenter, data_t, 1> srcData(srcbox, a_srcGraph, useStack);
      EBBoxData<dstCenter, data_t, 1> dstData(dstbox, a_dstGraph, useStack);

      unsigned int realSize = 0;
      for (int idst = 0; idst < a_dstVoFs.size(); idst++)
      {
		bool noEmpty = false;
        	const LocalStencil<srcCenter, data_t>& stenLoc = a_stencil[idst];
        	for (int isten = 0; isten < stenLoc.size(); isten++)
          	{
			auto entry = stenLoc.m_entries[isten];
          		if(entry.m_weight != 0)
			{
				noEmpty = true;
			}
		}
	
		if(noEmpty)
			realSize++;
      }

      vector<uint64_t>  h_stenSizes(realSize);
      vector<uint64_t>  h_stenStart(realSize);
      m_ebstencil.resize(realSize);
      m_dstAccess.resize(realSize);

      vector<pair_t<data_t> >    h_serializedStenInfo;
      uint64_t startLoc = 0;
      unsigned int it = 0;

      for (int idst = 0; idst < a_dstVoFs.size(); idst++)
      {
        const LocalStencil<srcCenter, data_t>& stenLoc = a_stencil[idst];
	bool noEmpty = false;

        data_t sumweight = 0;
#ifdef NEWSTRUCT
        assert(stenLoc.size()< std::pow<int>(2,16)-1);
#endif
	unsigned int it2 = 0;
        for (int isten = 0; isten < stenLoc.size(); isten++)
        {
          auto entry = stenLoc.m_entries[isten];
	  if(entry.m_weight != 0)
	  {
		pair_t<data_t> tmp;
          	tmp.first  = srcData.indexInfo(entry.m_vof, 0);
          	tmp.second = entry.m_weight;
		m_ebstencil[it].push_back(tmp);
          	sumweight += entry.m_weight;
          	h_serializedStenInfo.push_back(tmp);
#ifdef COUNTER
	  	if(m_ebstencil[it][it2].first.m_dataID == 0) m_counterRegVof++;
    	  	if(m_ebstencil[it][it2].second == 0) m_counterZeroWeight++;
#endif
        	m_flopsThisApply += 3;
		noEmpty = true;
		it2++;
        	//m_flopsThisApply += 3*stenLoc.size();
	  }
        }
//      printf("agd:sumweight = %f\n", sumweight);
     
	if(noEmpty) 
	{
        	m_dstAccess[it] = dstData.indexInfo(a_dstVoFs[idst], 0);

        	h_stenSizes[it] = it2;// a_stencil[idst].size();
        	h_stenStart[it] = startLoc;
        	startLoc         += it2; // a_stencil[idst].size();
		it++;
	}
      }
#ifdef COUNTER
      m_counterTotSten = h_serializedStenInfo.size();
#endif
#ifdef PROTO_CUDA
      cudaPrep(h_serializedStenInfo, h_stenSizes, h_stenStart);
#else
      hostPrep(h_serializedStenInfo, h_stenSizes, h_stenStart);
#endif
//    printf("dstvofs size = %d, flops per apply = %d\n ",(int)a_dstVoFs.size(), m_flopsThisApply);
    }
 
  }
  ///
  template <CENTERING srcCenter, CENTERING dstCenter, typename data_t>
  void 
  AggStencil<srcCenter, dstCenter, data_t>::
  hostPrep(const vector<pair_t<data_t> >   &   h_serializedStenInfo,
           const vector<uint64_t>          &   h_stenSizes,
           const vector<uint64_t>          &   h_stenStart)

  {
    if(m_dstAccess.size() > 0)
    {
      uint64_t a_numStenEntriesTotal = h_serializedStenInfo.size();
      uint64_t numptsDst = m_dstAccess.size();
      uint64_t totStenSize = a_numStenEntriesTotal*(sizeof(pair_t<data_t>));
      uint64_t destAccSize = numptsDst*(sizeof(EBDataLoc));
      uint64_t stenLenSize = numptsDst*(sizeof(uint64_t));
      uint64_t stStartSize = numptsDst*(sizeof(uint64_t));
      ///allocate memory on the device
      d_ebstencil = (pair_t<data_t>*)malloc(totStenSize);
      d_dstaccess = (EBDataLoc*     )malloc(destAccSize);
      d_stensizes = (uint64_t*      )malloc(stenLenSize);
      d_stenstart = (uint64_t*      )malloc(stStartSize);

      ///copy from the host
      memcpy(d_dstaccess, m_dstAccess.data(),           destAccSize);
      memcpy(d_stensizes, h_stenSizes.data(),           stenLenSize);
      memcpy(d_stenstart, h_stenStart.data(),           stStartSize);
      memcpy(d_ebstencil, h_serializedStenInfo.data(),  totStenSize);
    }
  }
  /**
     lphi = L(phi)
  */
  template <CENTERING srcCenter, CENTERING dstCenter, typename data_t>
  void 
  AggStencil<srcCenter, dstCenter, data_t>::
  apply(EBBoxData<      dstCenter, data_t, 1>     & a_lph,
        const EBBoxData<srcCenter, data_t, 1>     & a_phi,
        const data_t                              & a_scale,
        bool                                        a_incrementOnly) const
  {
    if(m_dstAccess.size() > 0)
    {
      PR_assert(m_srcBox == a_phi.inputBox());
      PR_assert(m_dstBox == a_lph.inputBox());
      PR_TIME("AggSten::apply");
#ifdef PROTO_CUDA
      cudaApply(a_lph, a_phi, a_scale, a_incrementOnly);
#else
      hostApply(a_lph, a_phi, a_scale, a_incrementOnly);
#endif
      PR_FLOPS(m_flopsThisApply);
    }
  }
///
  template <typename data_t>
  void 
  hostAggStencilIndexer(int a_begin, int a_end, int idx,
                        const data_t*           a_dataPtrsSrc[2],
                        data_t*                 a_dataPtrsDst[2],
                        const pair_t<data_t>*   a_ebstencil,
                        const uint64_t*         a_stensizes,
                        const uint64_t*         a_stenstart,
                        const EBDataLoc*        a_dstaccess,
                        bool                    a_incrementOnly,
                        data_t                  a_scale)
  {
    if((idx >= a_begin) && (idx< a_end))
    {
      

      uint64_t startInd = a_stenstart[idx];
      uint64_t stenSize = a_stensizes[idx];
      const EBDataLoc& dstAcc = a_dstaccess[idx];
      data_t& lphi = a_dataPtrsDst[dstAcc.m_dataID][dstAcc.m_offset];

      if(!a_incrementOnly)
      {
        lphi = 0;
      }
      data_t sumweight = 0;
      for(uint64_t isten = 0; isten < stenSize; isten++)
      {
        const pair_t<data_t>& entry = a_ebstencil[startInd + isten];
        const EBDataLoc& srcAcc = entry.first;	
        const data_t   & weight = entry.second;
        sumweight += weight;
        const data_t& phi = a_dataPtrsSrc[srcAcc.m_dataID][srcAcc.m_offset];
        lphi += a_scale*weight*phi;
      }
      //printf("agi:sumweight = %f\n", sumweight);
    }
  }
  template <CENTERING srcCenter, CENTERING dstCenter, typename data_t>
  void 
  AggStencil<    srcCenter, dstCenter, data_t>::
  hostApply(EBBoxData<      dstCenter, data_t, 1>     & a_lph,
            const EBBoxData<srcCenter, data_t, 1>     & a_phi,
            const data_t                              & a_scale,
            bool                                        a_incrementOnly) const

  {
    int begin = 0;
    int end   = m_dstAccess.size();
    const data_t* dataPtrsSrc[2];
    data_t*       dataPtrsDst[2];
    for(int idat = 0; idat < 2; idat++)
    {
      dataPtrsSrc[idat]  = (a_phi.dataPtr(idat, 0));
      dataPtrsDst[idat]  = (a_lph.dataPtr(idat, 0));
    }
    for(int idx = 0; idx < end; idx++)
    {
    ////begin debugging hooks
//      const auto& localsten = m_stencil[idx];
//      const auto& vof       = m_dstVoFs[idx];
#if DIM==2      
//      Point debpt(2, 1);
//      int ideb = 0;
//      if(vof.m_pt == debpt)
//      {
//        ideb = 1;
//      }
#endif      
      //end debugging hooks
      hostAggStencilIndexer(begin, end, idx, dataPtrsSrc, dataPtrsDst, d_ebstencil, d_stensizes, d_stenstart, d_dstaccess, a_incrementOnly, a_scale);
    }
  }

  template <CENTERING srcCenter, CENTERING dstCenter, typename data_t>
  void 
  AggStencil<srcCenter, dstCenter, data_t>::
  hostRelease()
  {
    if(m_dstAccess.size() > 0)
    {
      free(d_ebstencil);
      free(d_dstaccess);
      free(d_stensizes);
      free(d_stenstart);
    }
  }
///
/////
template<typename T>
struct pairPtr{
  T* ptr[2];
};

#ifdef COUNTER
  template <CENTERING srcCenter, CENTERING dstCenter, typename data_t>
  unsigned int 
  AggStencil<srcCenter, dstCenter, data_t>::
  getCounterVofs()
  {
    return  m_dstAccess.size();
  }

  template <CENTERING srcCenter, CENTERING dstCenter, typename data_t>
  unsigned int 
  AggStencil<srcCenter, dstCenter, data_t>::
  getCounterRegVof()
  {
    return m_counterRegVof;
  }
 
  template <CENTERING srcCenter, CENTERING dstCenter, typename data_t>
  unsigned int 
  AggStencil<srcCenter, dstCenter, data_t>::
  getCounterZeroWeight()
  {
    return m_counterZeroWeight;

  }

  template <CENTERING srcCenter, CENTERING dstCenter, typename data_t>
  unsigned int 
  AggStencil<srcCenter, dstCenter, data_t>::
  getCounterTotSten()
  {
    return m_counterTotSten;
  }
#endif

#ifdef PROTO_CUDA
  template <CENTERING srcCenter, CENTERING dstCenter, typename data_t>
  void 
  AggStencil<srcCenter, dstCenter, data_t>::
  cudaPrep(const vector<pair_t<data_t> >   &   h_serializedStenInfo,
           const vector<uint64_t>          &   h_stenSizes,
           const vector<uint64_t>          &   h_stenStart)

  {
    uint64_t a_numStenEntriesTotal = h_serializedStenInfo.size();
    if(a_numStenEntriesTotal>0)
      {
        uint64_t numptsDst = m_dstAccess.size();
	uint64_t stStartSize = numptsDst*(sizeof(uint64_t));

#ifdef NEWSTRUCT
// new test
	uint64_t destMeshTypeSize = numptsDst*(sizeof(bool));
        uint64_t destOffsetSize   = numptsDst*(sizeof(uint64_t));
	uint64_t srcMeshTypeSize = a_numStenEntriesTotal*(sizeof(bool));
        uint64_t srcOffsetSize   = a_numStenEntriesTotal*(sizeof(uint64_t));
        uint64_t srcWeightSize   = a_numStenEntriesTotal*(sizeof(data_t));
        uint64_t stenLenSize = numptsDst*(sizeof(uint16_t));

        protoMalloc(DEVICE, d_destMeshType, destMeshTypeSize);
        protoMalloc(DEVICE, d_destOffset, destOffsetSize);
        protoMalloc(DEVICE, d_srcMeshType, srcMeshTypeSize);
        protoMalloc(DEVICE, d_srcOffset, srcOffsetSize);
        protoMalloc(DEVICE, d_srcWeight, srcWeightSize);
        protoMalloc(DEVICE, d_stenLen, stenLenSize);
        protoMalloc(DEVICE, d_stenstart, stStartSize);

        StackAlloc<HOST>& allocator  = Proto::DisjointBoxLayout::getStackAllocator<HOST>();

        bool* tmp_destMeshType    = (bool*)     (allocator.alloc(destMeshTypeSize));
        uint64_t* tmp_destOffset  = (uint64_t*) (allocator.alloc(destOffsetSize));
        bool* tmp_srcMeshType     = (bool*)     (allocator.alloc(srcMeshTypeSize));
        uint64_t* tmp_srcOffset   = (uint64_t*) (allocator.alloc(srcOffsetSize));
        data_t* tmp_srcWeight     = (data_t*)   (allocator.alloc(srcWeightSize));
        uint16_t* tmp_stenLen     = (uint16_t*) (allocator.alloc(stenLenSize));

        for(unsigned int it = 0 ; it < numptsDst ; it++)
        {
		if(m_dstAccess[it].m_dataID == 1) tmp_destMeshType[it] = true; 
		else tmp_destMeshType[it] = false;
		tmp_destOffset[it] = m_dstAccess[it].m_offset;
		tmp_stenLen[it] = h_stenSizes[it];
        }

        for(unsigned int it = 0 ; it < a_numStenEntriesTotal ; it++)
        {
		if(h_serializedStenInfo[it].first.m_dataID==1) tmp_srcMeshType[it] = true;
		else tmp_srcMeshType[it] = false;
		tmp_srcOffset[it] = h_serializedStenInfo[it].first.m_offset;
		tmp_srcWeight[it] = h_serializedStenInfo[it].second;
        }

        protoMemcpy(DEVICE, d_destMeshType, tmp_destMeshType,   destMeshTypeSize, protoMemcpyHostToDevice);
        protoMemcpy(DEVICE, d_destOffset,   tmp_destOffset,     destOffsetSize,   protoMemcpyHostToDevice);
        protoMemcpy(DEVICE, d_stenLen,      tmp_stenLen,        stenLenSize,      protoMemcpyHostToDevice);
        protoMemcpy(DEVICE, d_srcMeshType,  tmp_srcMeshType,    srcMeshTypeSize,  protoMemcpyHostToDevice);
        protoMemcpy(DEVICE, d_srcOffset,    tmp_srcOffset,      srcOffsetSize,    protoMemcpyHostToDevice);
        protoMemcpy(DEVICE, d_srcWeight,    tmp_srcWeight,      srcWeightSize,    protoMemcpyHostToDevice);
        protoMemcpy(DEVICE, d_stenstart,    h_stenStart.data(), stStartSize,      protoMemcpyHostToDevice);

	allocator.free(tmp_stenLen);
	allocator.free(tmp_srcWeight);
	allocator.free(tmp_srcOffset);
	allocator.free(tmp_srcMeshType);
	allocator.free(tmp_destOffset);
	allocator.free(tmp_destMeshType);
#else
        uint64_t stenLenSize = numptsDst*(sizeof(uint64_t));
        uint64_t totStenSize = a_numStenEntriesTotal*(sizeof(pair_t<data_t>));
        uint64_t destAccSize = numptsDst*(sizeof(EBDataLoc));

        ///allocate memory on the device
        protoMalloc(DEVICE,d_ebstencil, totStenSize);
        protoMalloc(DEVICE,d_dstaccess, destAccSize);
        protoMalloc(DEVICE,d_stensizes, stenLenSize);
        protoMalloc(DEVICE,d_stenstart, stStartSize);

        ///copy from the host
        protoMemcpy(DEVICE,d_dstaccess, m_dstAccess.data(),           destAccSize, protoMemcpyHostToDevice);
        protoMemcpy(DEVICE,d_stensizes, h_stenSizes.data(),           stenLenSize, protoMemcpyHostToDevice);
        protoMemcpy(DEVICE,d_stenstart, h_stenStart.data(),           stStartSize, protoMemcpyHostToDevice);
        protoMemcpy(DEVICE,d_ebstencil, h_serializedStenInfo.data(),  totStenSize, protoMemcpyHostToDevice);
#endif
        //figure out blocks and strides
        m_stride = 32;
        bool found = false;
        while(!found)
          {
            if((numptsDst % m_stride) == 0)
              {
                found = true;
              }
            else
              {
                m_stride /= 2;
              }
          }
        PROTO_ASSERT(found, "could not find right number of blocks and strides");
        m_blocks = numptsDst/m_stride;
        m_blocks = 1;
        m_stride =numptsDst;
      }
 
  }
    
/////

  template <CENTERING srcCenter, CENTERING dstCenter, typename data_t>
  void 
  AggStencil<srcCenter, dstCenter, data_t>::
  cudaRelease()
  {
#ifdef NEWSTRUCT
    if(d_stenstart != nullptr)
      {
	protoFree(DEVICE,d_destMeshType);
	protoFree(DEVICE,d_destOffset);
	protoFree(DEVICE,d_stenLen);
	protoFree(DEVICE,d_srcMeshType);
	protoFree(DEVICE,d_srcOffset);
	protoFree(DEVICE,d_srcWeight);
        protoFree(DEVICE,d_stenstart);
      }
#else
    if(d_ebstencil != nullptr)
      {
        protoFree(MEMTYPE_DEFAULT,d_ebstencil);
        protoFree(MEMTYPE_DEFAULT,d_dstaccess);
        protoFree(MEMTYPE_DEFAULT,d_stensizes);
        protoFree(MEMTYPE_DEFAULT,d_stenstart);
      }
#endif
  }

  template <typename data_t>
  struct aggStencilIndexer
  {
 

#ifdef NEWSTRUCT 
  static
  inline void 
  cpu(int a_begin, int a_end,
		    pairPtr<const data_t>   dataPtrsSrc,
		    pairPtr<data_t>         dataPtrsDst,
                    const uint16_t*         a_stensizes,
                    const uint64_t*         a_stenstart,
                    const bool*             a_dstMeshType,
                    const uint64_t*         a_dstOffset,
                    const bool*             a_srcMeshType,
                    const uint64_t*         a_srcOffset,
                    const data_t*           a_srcWeight,
                    bool                    a_incrementOnly,
                    data_t                  a_scale)
  {
    assert(false);
  }
   
  __device__ static
  void 
  gpu(int a_begin, int a_end,
		    pairPtr<const data_t>   dataPtrsSrc,
		    pairPtr<data_t>         dataPtrsDst,
                    const uint16_t*         a_stensizes,
                    const uint64_t*         a_stenstart,
                    const bool*             a_dstMeshType,
                    const uint64_t*         a_dstOffset,
                    const bool*             a_srcMeshType,
                    const uint64_t*         a_srcOffset,
                    const data_t*           a_srcWeight,
                    bool                    a_incrementOnly,
                    data_t                  a_scale)
  {
// test
    __shared__ const data_t* tmpSrc[2];
    __shared__ data_t* tmpDst[2];
    __shared__ const data_t* Src[512][2];
    __shared__ data_t* Dst[512][2];

    if(threadIdx.x == 0)
    {
      tmpSrc[0] = dataPtrsSrc.ptr[0];
      tmpSrc[1] = dataPtrsSrc.ptr[1];
      tmpDst[0] = dataPtrsDst.ptr[0];
      tmpDst[1] = dataPtrsDst.ptr[1];
    }

    __syncthreads();

    Src[threadIdx.x][0] = tmpSrc[0];
    Src[threadIdx.x][1] = tmpSrc[1];
    Dst[threadIdx.x][0] = tmpDst[0];
    Dst[threadIdx.x][1] = tmpDst[1];
    
    unsigned int stride = blockDim.x  * gridDim.x;

    for(unsigned int idx = a_begin + blockIdx.x * blockDim.x + threadIdx.x ; idx < a_end ; idx += stride)
    {
      const uint64_t startInd = a_stenstart[idx];
      const uint16_t stenSize = a_stensizes[idx];
      data_t acc = 0;

      #pragma loop unroll
      for(uint16_t isten = 0; isten < stenSize; isten++)
      {
	uint64_t srcIdx = startInd + isten;
        const data_t weight = a_srcWeight[srcIdx];
	const data_t* src = Src[threadIdx.x][a_srcMeshType[srcIdx]];
        const data_t phi = src[a_srcOffset[srcIdx]];
        acc += weight*phi;
      }

      data_t* dst = Dst[threadIdx.x][a_dstMeshType[idx]];
      data_t& lphi = dst[a_dstOffset[idx]];

      if(!a_incrementOnly)
      {
        lphi = a_scale*acc;
      }
      else lphi += a_scale*acc;
    }
  }
#else
  static
  inline void 
  cpu(int a_begin, int a_end,
		    pairPtr<const data_t>   dataPtrsSrc,
		    pairPtr<data_t>         dataPtrsDst,
                    const pair_t<data_t>*   a_ebstencil,
                    const uint64_t*         a_stensizes,
                    const uint64_t*         a_stenstart,
                    const EBDataLoc*        a_dstaccess,
                    bool                    a_incrementOnly,
                    data_t                  a_scale)
  {
    assert(false);
  }
   
  __device__ static
  void 
  gpu(int a_begin, int a_end,
		    pairPtr<const data_t>   dataPtrsSrc,
		    pairPtr<data_t>         dataPtrsDst,
                    const pair_t<data_t>*   a_ebstencil,
                    const uint64_t*         a_stensizes,
                    const uint64_t*         a_stenstart,
                    const EBDataLoc*        a_dstaccess,
                    bool                    a_incrementOnly,
                    data_t                  a_scale)
  {

    int idx = a_begin + threadIdx.x+blockIdx.x*blockDim.x;
    
    if(idx< a_end)
    {
      const uint64_t startInd = a_stenstart[idx];
      const uint64_t stenSize = a_stensizes[idx];
      const EBDataLoc& dstAcc = a_dstaccess[idx];

      data_t acc = 0;

      #pragma loop unroll
      for(uint64_t isten = 0; isten < stenSize; isten++)
      {
        const pair_t<data_t>& entry = a_ebstencil[startInd + isten];
        const EBDataLoc& srcAcc = entry.first;
        const data_t weight = entry.second;
        const data_t phi = dataPtrsSrc.ptr[srcAcc.m_dataID][srcAcc.m_offset];
        acc += weight*phi;
      }
      data_t& lphi = dataPtrsDst.ptr[dstAcc.m_dataID][dstAcc.m_offset];
      if(!a_incrementOnly)
      {
        lphi = a_scale*acc;
      }
      else lphi += a_scale*acc;
    }
  }
#endif
 };
/////

  template <CENTERING srcCenter, CENTERING dstCenter, typename data_t>
  void 
  AggStencil<    srcCenter, dstCenter, data_t>::
  cudaApply(EBBoxData<      dstCenter, data_t, 1>     & a_lph,
            const EBBoxData<srcCenter, data_t, 1>     & a_phi,
            const data_t                              & a_scale,
            bool                    a_incrementOnly) const

  {
    int begin = 0;
    int end   = m_dstAccess.size();
    pairPtr<const data_t> srcPtr;
    pairPtr<data_t> dstPtr;
    srcPtr.ptr[0] = a_phi.dataPtr( 0, 0);
    srcPtr.ptr[1] = a_phi.dataPtr( 1, 0);
    dstPtr.ptr[0] = a_lph.dataPtr( 0, 0);
    dstPtr.ptr[1] = a_lph.dataPtr( 1, 0);

    // m_stride or m_blocks doesn't fit. This function is temporary used.
    int stride, blocks;
//    v100tuning(end-begin, blocks, stride);
    stride = 256;
    blocks = (end + stride - 1) / stride;

    if(end - begin == 0) return;
		    
#ifdef NEWSTRUCT
      int numSMs=1;
#ifdef PROTO_CUDA
      int devId = 0;
#ifdef PROTO_HIP
      hipDeviceGetAttribute(&numSMs, hipDeviceAttributeMultiprocessorCount, devId);
#else
      cudaDeviceGetAttribute(&numSMs, cudaDevAttrMultiProcessorCount, devId);
#endif
#endif
    protoLaunchKernelT<MEMTYPE_DEFAULT,aggStencilIndexer<data_t>>(32*numSMs, 512,  
      begin, end, srcPtr, dstPtr, d_stenLen, d_stenstart, d_destMeshType, d_destOffset, d_srcMeshType, d_srcOffset, d_srcWeight, a_incrementOnly, a_scale); 
#else
    protoLaunchKernelT<MEMTYPE_DEFAULT,aggStencilIndexer<data_t>>( blocks, stride,  
      begin, end, srcPtr, dstPtr, d_ebstencil, d_stensizes, d_stenstart, d_dstaccess, a_incrementOnly, a_scale); 
#endif
  }
#endif
}
///////////
#endif
