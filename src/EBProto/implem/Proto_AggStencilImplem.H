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
    
    m_ebstencil.resize(a_dstVoFs.size());
    m_dstAccess.resize(a_dstVoFs.size());

    //begin debugging hooks
    m_dstVoFs = a_dstVoFs;
    m_stencil = a_stencil;
    //end   debugging hooks
    if(a_dstVoFs.size() > 0)
    {
      Box srcbox = a_srcValid.grow(a_srcGhost);
      Box dstbox = a_dstValid.grow(a_dstGhost);

      m_srcBox = srcbox;
      m_dstBox = dstbox;

      m_flopsThisApply = 0;        

      //turning on the stack here is a terrible idea because nothing else uses it
      bool useStack = false; 
      EBBoxData<srcCenter, data_t, 1> srcData(srcbox, a_srcGraph, useStack);
      EBBoxData<dstCenter, data_t, 1> dstData(dstbox, a_dstGraph, useStack);

      vector<uint64_t>  h_stenSizes(a_dstVoFs.size());
      vector<uint64_t>  h_stenStart(a_dstVoFs.size());
      vector<pair_t<data_t> >    h_serializedStenInfo;
      uint64_t startLoc = 0;

      for (int idst = 0; idst < a_dstVoFs.size(); idst++)
      {


        m_dstAccess[idst] = dstData.indexInfo(a_dstVoFs[idst], 0);
        m_ebstencil[idst].resize(a_stencil[idst].size());

        const LocalStencil<srcCenter, data_t>& stenLoc = a_stencil[idst];
        m_flopsThisApply += 3*stenLoc.size();

        data_t sumweight = 0;
        for (int isten = 0; isten < stenLoc.size(); isten++)
        {
          auto entry = stenLoc.m_entries[isten];
          m_ebstencil[idst][isten].first  = srcData.indexInfo(entry.m_vof, 0);
          m_ebstencil[idst][isten].second = entry.m_weight;
          sumweight += entry.m_weight;
          h_serializedStenInfo.push_back(m_ebstencil[idst][isten]);
        }
//      printf("agd:sumweight = %f\n", sumweight);
      

        h_stenSizes[idst] = a_stencil[idst].size();
        h_stenStart[idst] = startLoc;
        startLoc         += a_stencil[idst].size();

      }
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
        bool                                        a_incrementOnly,
        bool a_printStuff) const
  {
    if(m_dstAccess.size() > 0)
    {
      PR_assert(m_srcBox == a_phi.inputBox());
      PR_assert(m_dstBox == a_lph.inputBox());
      PR_TIME("AggSten::apply");
#ifdef PROTO_CUDA
      cudaApply(a_lph, a_phi, a_scale, a_incrementOnly);
#else
      hostApply(a_lph, a_phi, a_scale, a_incrementOnly, a_printStuff);
#endif
      PR_FLOPS(m_flopsThisApply);
    }
  }
///
  template <CENTERING srcCenter, CENTERING dstCenter, typename data_t>
  void 
  AggStencil<    srcCenter, dstCenter, data_t>::
  hostAggStencilIndexer(int a_begin, int a_end, int idx,
                        const data_t*           a_dataPtrsSrc[2],
                        data_t*                 a_dataPtrsDst[2],
                        const pair_t<data_t>*   a_ebstencil,
                        const uint64_t*         a_stensizes,
                        const uint64_t*         a_stenstart,
                        const EBDataLoc*        a_dstaccess,
                        bool                    a_incrementOnly,
                        data_t                  a_scale,
                        bool a_printStuff) const
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
            bool                                        a_incrementOnly,
            bool                                        a_printStuff) const

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
      bool printStuff = false;
      //begin debugging hooks
#if DIM==2      
      Point debpt(6, 16);
      auto vof = m_dstVoFs[idx];
      if(vof.m_pt == debpt)
      {
        printStuff = true;
      }
#endif      
      //end debugging hooks
      hostAggStencilIndexer(begin, end, idx, dataPtrsSrc, dataPtrsDst, d_ebstencil, d_stensizes, d_stenstart, d_dstaccess, a_incrementOnly, a_scale, printStuff);
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
        uint64_t totStenSize = a_numStenEntriesTotal*(sizeof(pair_t<data_t>));
        uint64_t destAccSize = numptsDst*(sizeof(EBDataLoc));
        uint64_t stenLenSize = numptsDst*(sizeof(uint64_t));
        uint64_t stStartSize = numptsDst*(sizeof(uint64_t));
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
    if(d_ebstencil != nullptr)
      {
        protoFree(MEMTYPE_DEFAULT,d_ebstencil);
        protoFree(MEMTYPE_DEFAULT,d_dstaccess);
        protoFree(MEMTYPE_DEFAULT,d_stensizes);
        protoFree(MEMTYPE_DEFAULT,d_stenstart);
      }
  }

  template <typename data_t>
  struct aggStencilIndexer
  {
 
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
    //protoLaunchKernel(aggStencilIndexer, blocks, stride,  
    //  begin, end, srcPtr, dstPtr, d_ebstencil, d_stensizes, d_stenstart, d_dstaccess, a_incrementOnly, a_scale);
    protoLaunchKernelT<MEMTYPE_DEFAULT,aggStencilIndexer<data_t>>( blocks, stride,  
      begin, end, srcPtr, dstPtr, d_ebstencil, d_stensizes, d_stenstart, d_dstaccess, a_incrementOnly, a_scale);
  }
#endif
}
///////////
#endif
