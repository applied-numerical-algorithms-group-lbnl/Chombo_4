#ifndef _Proto_IrregDataImplem_H_
#define _Proto_IrregDataImplem_H_

namespace Proto
{
  ////
  template <CENTERING cent, class data_t,  unsigned int ncomp>
  void 
  IrregData<cent, data_t, ncomp>::
  define(const shared_ptr<vector<index_t> >& a_indicies, bool a_stackAllocation)
  {
    m_indicies = a_indicies;
    std::sort(m_indicies->begin(), m_indicies->end());

    if(a_stackAllocation)
    {
      if(size()==0)
      {
        m_rawPtr = nullptr;
	m_data = ::std::shared_ptr<data_t>(nullptr);
      }
      else
      {
        StackAlloc<MEMTYPE_DEFAULT>& allocator  = Proto::DisjointBoxLayout::getStackAllocator<MEMTYPE_DEFAULT>();
        m_rawPtr = (data_t*)(allocator.alloc(size()*sizeof(data_t)));
        m_data = ::std::shared_ptr<data_t>(m_rawPtr, &(stack_deleter_boxdata<MEMTYPE_DEFAULT>));

      }
    }
    else
    { 
      protoMalloc(MEMTYPE_DEFAULT,m_rawPtr, size()*sizeof(data_t));
      m_data = ::std::shared_ptr<data_t>(m_rawPtr, [](data_t* p){protoFree(MEMTYPE_DEFAULT,p);});
    }
  }

  ////
  template <CENTERING cent, class data_t,  unsigned int ncomp>
  void 
  IrregData<cent, data_t, ncomp>::
  setVal(const data_t& a_val)
  {
    unsigned int nsize = size();
    if(nsize == 0) return;
    int blockSize = 512;
    int nBlocks = (nsize + blockSize - 1) / blockSize;
    //protoLaunchKernelMemAsync((opKernel<data_t,BoxDataOp::Assign>), nBlocks, blockSize, 0, protoGetCurrentStream, m_rawPtr, a_val, nsize);
    protoLaunchKernelMemAsyncT<MEMTYPE_DEFAULT,
                               opKernel<data_t,BoxDataOp::Assign>>( nBlocks, blockSize, 0, protoGetCurrentStream, m_rawPtr, a_val, nsize);
  }

}
#endif

