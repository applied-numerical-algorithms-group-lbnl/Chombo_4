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
    m_ebstencil.resize(a_dstVoFs.size());
    m_dstAccess.resize(a_dstVoFs.size());

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

      EBBoxData<srcCenter, data_t, 1> srcData(srcbox, a_srcGraph);
      EBBoxData<dstCenter, data_t, 1> dstData(dstbox, a_dstGraph);

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

        h_stenSizes[idst] = a_stencil[idst].size();
        h_stenStart[idst] = startLoc;
        startLoc         += a_stencil[idst].size();

      }
      EBProtoPrep(h_serializedStenInfo, h_stenSizes, h_stenStart);

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
      EBProtoApply(a_lph, a_phi, a_scale, a_incrementOnly);
      PR_FLOPS(m_flopsThisApply);
    }
  }

  template <CENTERING srcCenter, CENTERING dstCenter, typename data_t>
  void 
  AggStencil<srcCenter, dstCenter, data_t>::
  EBProtoPrep(const vector<pair_t<data_t> >   &   h_serializedStenInfo,
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
        protoMalloc(d_ebstencil, totStenSize);
        protoMalloc(d_dstaccess, destAccSize);
        protoMalloc(d_stensizes, stenLenSize);
        protoMalloc(d_stenstart, stStartSize);

        ///copy from the host
        protoMemcpy(d_dstaccess, m_dstAccess.data(),           destAccSize, protoMemcpyHostToDevice);
        protoMemcpy(d_stensizes, h_stenSizes.data(),           stenLenSize, protoMemcpyHostToDevice);
        protoMemcpy(d_stenstart, h_stenStart.data(),           stStartSize, protoMemcpyHostToDevice);
        protoMemcpy(d_ebstencil, h_serializedStenInfo.data(),  totStenSize, protoMemcpyHostToDevice);

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
  EBProtoRelease()
  {
    if(d_ebstencil != nullptr)
      {
        protoFree(d_ebstencil);
        protoFree(d_dstaccess);
        protoFree(d_stensizes);
        protoFree(d_stenstart);
      }
  }
/////
template<typename T>
struct pairPtr{
  T* ptr[2];
};

#ifdef PROTO_CUDA

  template <typename data_t>
  __global__
  void 
  aggStencilIndexer(int a_begin, int a_end,
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

#else

  template <typename data_t>
  void 
  aggStencilIndexer(int a_begin, int a_end,
		    pairPtr<const data_t>   dataPtrsSrc,
		    pairPtr<data_t>         dataPtrsDst,
                    const pair_t<data_t>*   a_ebstencil,
                    const uint64_t*         a_stensizes,
                    const uint64_t*         a_stenstart,
                    const EBDataLoc*        a_dstaccess,
                    bool                    a_incrementOnly,
                    data_t                  a_scale)
  {
    for(int idx = a_begin ; idx < a_end ; idx++)
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

/////

  inline void v100tuning(int nbElems, int & nbBlocks, int &blockSize)
  {
    // determine the best block size and block dim
    blockSize = 256;
    nbBlocks = ( nbElems + blockSize - 1)/ blockSize;
    if(nbBlocks < 80)
    {
      // On V100 we want at least 80 blocks;
      nbBlocks = 80;
      // figure out what is the blockSize for 80 blocks
      // 
      blockSize = ( nbElems + nbBlocks - 1) / nbBlocks;
      // as we use __syncthreads(), we want that stride modulo 32 is equal to 0
      int coeff = blockSize / 32;
      if(coeff > 0)
      {
        blockSize = coeff * 32;
        // recompute the number of blocks > 80
        nbBlocks = ( nbElems + blockSize - 1) / nbBlocks;
      }
    }
  }

  template <CENTERING srcCenter, CENTERING dstCenter, typename data_t>
  void 
  AggStencil<    srcCenter, dstCenter, data_t>::
  EBProtoApply(EBBoxData<      dstCenter, data_t, 1>     & a_lph,
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
    v100tuning(end-begin, blocks, stride);

    //protoLaunchKernel(aggStencilIndexer, m_stride, m_blocks, 
    protoLaunchKernel(aggStencilIndexer, blocks, stride,  
      begin, end, srcPtr, dstPtr, d_ebstencil, d_stensizes, d_stenstart, d_dstaccess, a_incrementOnly, a_scale);
  }
}
///////////
#endif
