#ifndef _Proto_IrregDataImplem_H_
#define _Proto_IrregDataImplem_H_

namespace Proto
{
  ////
  template <CENTERING cent, class data_t,  unsigned int ncomp>
  void 
  IrregData<cent, data_t, ncomp>::
  define(const shared_ptr<vector<index_t> >& a_indicies)
  {
    m_indicies = a_indicies;
    std::sort(m_indicies->begin(), m_indicies->end());
#ifdef PROTO_CUDA
    protoMalloc(&m_rawPtr, size()*sizeof(data_t));
    protoError err = protoGetLastError();
    if (err != protoSuccess)
    {
      fprintf(stderr, "protoGetLastError() failed at %s:%i : %s\n",
              __FILE__, __LINE__, protoGetErrorString(err));
    }

    m_data = ::std::shared_ptr<data_t>(m_rawPtr, [](data_t* p){protoFree(p);});
#else
    std::size_t count = this->size();
    m_data=::std::shared_ptr<data_t>(new data_t[count], [](data_t* p) { delete[] p;});
    m_rawPtr = m_data.get();
#endif
    makeUglyStruct();
  }
  /****/
  template <CENTERING cent, class data_t,  unsigned int ncomp>
  void 
  IrregData<cent, data_t, ncomp>::
  makeUglyStruct()
  {
    int vecsize = m_indicies->size();
    m_uglystruct =  shared_ptr< vector< EBIrregStruct<cent, data_t, ncomp> > >(new vector< EBIrregStruct<cent, data_t, ncomp> >(vecsize) );
    for(int ivec = 0; ivec < vecsize; ivec++)
    {
      EBIrregStruct<cent, data_t, ncomp>&  vecval = (*m_uglystruct)[ivec];
      auto& ugindex = (*m_indicies)[ivec];
      vecval.m_startPtr = this->data();
      vecval.m_varsize  = this->vecsize();
      vecval.m_offset   = this->index(ugindex, 0);
      vecval.m_index    = ugindex.m_pt;
    }
  }
  ////
  template <CENTERING cent, class data_t,  unsigned int ncomp>
  void 
  IrregData<cent, data_t, ncomp>::
  setVal(const data_t& a_val)
  {
#ifdef PROTO_CUDA
  unsigned int nsize = size();
  if(nsize == 0) return;
  int blockSize = 512;
  int nBlocks = (nsize + blockSize - 1) / blockSize;
  auto streamid = DisjointBoxLayout::getCurrentStream();
  protoLaunchKernelMemAsync(setValKernel, nBlocks, blockSize, 0, streamid, m_rawPtr, nsize, a_val);

#ifdef THRUST_DEBUG
  protoError err = protoGetLastError();
  if (err != protoSuccess)
  {
    fprintf(stderr, "protoGetLastError() failed at %s:%i : %s\n",
            __FILE__, __LINE__, protoGetErrorString(err));
  }
#endif

#else
    for(int ival = 0; ival < size(); ival++)
    {
      m_rawPtr[ival] = a_val;
    }
#endif
  }

}
#endif

