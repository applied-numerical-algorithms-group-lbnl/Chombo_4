
#ifndef _Proto_EBForAllImplem_H_
#define _Proto_EBForAllImplem_H_

namespace Proto
{

  template<typename T>
  inline T&
  getIrregData(T& a_s)
  {
    return a_s;
  }

  template <CENTERING cent, typename  data_t, unsigned int ncomp>
  inline IrregData<cent, data_t, ncomp>&
  getIrregData(EBBoxData<cent, data_t, ncomp>& a_s)
  {
    return a_s.getIrregData();
  }



  template <CENTERING cent, typename  data_t, unsigned int ncomp>
  inline IrregData<cent, data_t, ncomp>&
  getIrregData(const EBBoxData<cent, data_t, ncomp>& a_s)
  {
    return (IrregData<cent, data_t, ncomp>&)(a_s.getIrregData());
  }


  template<typename data_t, unsigned int ncomp>
  using uglystruct_t = EBIrregStruct<data_t, ncomp>;

  template <CENTERING cent, typename T>
  inline T
  getEBIrregStruct(const vector<EBIndex<cent> >& a_indices,
                   T& a_T)
  {
    return a_T;
  }

  template <CENTERING cent, typename  data_t, unsigned int ncomp>
  inline vector< EBIrregStruct<data_t, ncomp> >
  getEBIrregStruct(const vector<EBIndex<cent> >& a_indices,
                   IrregData<cent, data_t, ncomp>& a_s )
  {
    vector< EBIrregStruct<data_t, ncomp> > retval;
    for(int ivec = 0; ivec < a_indices.size(); ivec++)
    {
      EBIrregStruct<data_t, ncomp>  vecval;
      vecval.m_startPtr = a_s.data();
      vecval.m_varsize  = a_s.vecsize();
      vecval.m_offset   = a_s.index(a_indices[ivec], 0);
      vecval.m_index    = a_indices[ivec].m_pt;
      retval.push_back(vecval);
    }


    return retval;
  }


  template<typename T>
  inline T&
  getBoxData(T& a_s)
  {
    return a_s;
  }

  template <CENTERING cent, typename  data_t, unsigned int ncomp>
  inline BoxData<data_t, ncomp>&
  getBoxData(EBBoxData<cent, data_t, ncomp>& a_s)
  {
    return a_s.getRegData();
  }


  template <CENTERING cent, typename  data_t, unsigned int ncomp>
  inline const BoxData<data_t, ncomp>&
  getBoxData(const EBBoxData<cent, data_t, ncomp>& a_s)
  {
    return a_s.getRegData();
  }

  template<typename... Srcs>
  inline void
  emptyFunc(Srcs... a_srcs)
  {
  }

#ifdef PROTO_CUDA

  template<typename T>
  __device__ __host__
  inline T
  cudaGetVar(unsigned int ivec,  T a_s)
  {
    return a_s;
  }


  template<typename data_t, unsigned int ncomp>
  __device__ 
  inline Var<data_t, ncomp>
  cudaGetVar(unsigned int a_ivec,
             EBIrregStruct<data_t, ncomp>* a_dst)
  {
    Var<data_t, ncomp> retval;

    const EBIrregStruct<data_t, ncomp>*  rawptr = a_dst;
    const EBIrregStruct<data_t, ncomp>&  ugly   = rawptr[a_ivec];
    for(int icomp = 0; icomp < ncomp; icomp++)
    {
      retval.m_ptrs[icomp] = ugly.m_startPtr + (ugly.m_varsize*icomp);
    }
    return retval;
  }




  template< typename data_t, unsigned int ncomp, typename Func, typename... Srcs>
  __global__
  void
  vec_indexer(unsigned int a_begin, unsigned int a_end, Func a_body, 
              EBIrregStruct<data_t,ncomp>*  a_dst, Srcs... a_srcs)
  {
    int idx = threadIdx.x + blockIdx.x*blockDim.x;
    if (idx >= a_begin && idx < a_end)
    {
      a_body(cudaGetVar(idx, a_dst), cudaGetVar(idx, a_srcs)...);
    }
  }


  template< typename data_t, unsigned int ncomp, typename Func, typename... Srcs>
  __global__
  void
  vec_indexer_i(unsigned int a_begin, unsigned int a_end,Func a_body, 
                EBIrregStruct<data_t,ncomp>*  a_dst, Srcs... a_srcs)
  {
    int idx = threadIdx.x + blockIdx.x*blockDim.x;
    if (idx >= a_begin && idx < a_end)
    {
      a_body(a_dst[idx].m_index.m_tuple, cudaGetVar(idx, a_dst), cudaGetVar(idx, a_srcs)...);
    }
  }

//
#else

  template<typename T>
  inline T
  getVar(unsigned int ivec,  T a_s)
  {
    return a_s;
  }

  template<typename data_t, unsigned int ncomp>
  inline Var<data_t, ncomp>
  getVar(unsigned int a_ivec,
             EBIrregStruct<data_t, ncomp>* a_dst)
  {
    Var<data_t, ncomp> retval;

    const EBIrregStruct<data_t, ncomp>*  rawptr = a_dst;
    const EBIrregStruct<data_t, ncomp>&  ugly   = rawptr[a_ivec];

    for(int icomp = 0; icomp < ncomp; icomp++)
    {
      retval.m_ptrs[icomp] = ugly.m_startPtr + ugly.m_offset + (ugly.m_varsize*icomp);
    }
    return retval;
  }

  template<typename data_t,unsigned int ncomp,  typename Func, typename... Srcs>
  inline void
  vec_indexer(unsigned int a_begin, unsigned int a_end, const Func& a_F, 
		EBIrregStruct<data_t, ncomp>* a_dst, Srcs... a_srcs)
  {
    for(unsigned int ivec = a_begin; ivec < a_end; ivec++)
    {
      a_F(getVar(ivec, a_dst), (getVar(ivec, a_srcs))...);
    }
       
  }

  template<typename data_t, unsigned int ncomp, typename Func, typename... Srcs>
  inline void
  vec_indexer_i(unsigned int a_begin, unsigned int a_end, const Func& a_F, 
		EBIrregStruct<data_t,ncomp>* a_dst, Srcs... a_srcs)
  {
    for(unsigned int ivec = a_begin; ivec < a_end; ivec++)
    {
      Point pt = a_dst[ivec].m_index;
      a_F(pt.m_tuple, getVar(ivec, a_dst), (getVar(ivec, a_srcs))...);
    }
       
  }

#endif


//going into this srcs are EBIrregStruct<cent, data_t, ncomp>* and other stuff
  template<typename  data_t, unsigned int ncomp, typename Func, typename... Srcs>
  inline void
  protoVectorFunc(const Func& a_F, unsigned int a_Nvec, 
                 EBIrregStruct<data_t,ncomp> * a_dst,Srcs... a_srcs)
  {
    PR_TIME("vec_indexer");
    const int N = a_Nvec;
    unsigned int stride = 256;
    unsigned int blocks = (N+stride-1)/stride;

    size_t smem = 0;
    protoLaunchKernelMemAsync(vec_indexer, blocks, stride, smem, protoGetCurrentStream,
      0, N, a_F, a_dst, a_srcs...);
  }


//going into this srcs are EBIrregStruct<cent, data_t, ncomp>* and other stuff
  template<typename  data_t, unsigned int ncomp, typename Func, typename... Srcs>
  inline void
  protoVectorFunc_i(const Func& a_F, unsigned int a_Nvec, 
                   EBIrregStruct<data_t,ncomp> * a_dst,Srcs... a_srcs)
  {
    PR_TIME("vec_indexer_i");
    const int N = a_Nvec;
    unsigned int stride = 256;
    unsigned int blocks = (N+stride-1)/stride;
    
    size_t smem = 0;
    protoLaunchKernelMemAsync(vec_indexer_i, blocks, stride, smem, protoGetCurrentStream,
      0, N, a_F, a_dst, a_srcs...);
  }


  template <typename T>
  inline T protogetEBIrregStruct(T& a_T)
  {
    return a_T;
  }

  template <CENTERING cent, typename  data_t, unsigned int ncomp>
  inline  EBIrregStruct<data_t, ncomp>*
  protogetEBIrregStruct(IrregData<cent, data_t, ncomp>& a_s )
  {
    return a_s.getEBIrregDataPtr();
  }

//going into this srcs are IrregDatas and other stuff
  template<CENTERING cent, typename  data_t, unsigned int ncomp, typename Func, typename... Srcs>
  inline void
  protoEBForAllIrreg(const Func& a_F, 
                    IrregData<cent, data_t, ncomp>& a_dst,
                    Srcs&...  a_srcs)
  {

    //indicies into irreg vector that correspond to input box
    unsigned int vecsize = a_dst.vecsize();
    if(vecsize > 0)
    {
      protoVectorFunc(a_F, vecsize, protogetEBIrregStruct(a_dst), 
                     (protogetEBIrregStruct(a_srcs))...);

    }
  }


///going into this srcs are IrregDatas and other stuff
  template<CENTERING cent, typename  data_t, unsigned int ncomp, typename Func, typename... Srcs>
  inline void
  protoEBForAllIrreg_i(const Func& a_F, 
                      IrregData<cent, data_t, ncomp>& a_dst,
                      Srcs&...  a_srcs)
  {
    //indicies into irreg vector that correspond to input box
    unsigned int vecsize = a_dst.vecsize();
    if(vecsize > 0)
    {
      protoVectorFunc_i(a_F, vecsize, protogetEBIrregStruct(a_dst), 
                       (protogetEBIrregStruct(a_srcs))...);

    }
  }

  ///
  template<typename Func, typename... Srcs>
  inline void
  protoEBforall(const Func & a_F,  Box a_box, Srcs&... a_srcs)
  {
    PR_TIME("protoEBforall");	  
    forallInPlaceBase(a_F, a_box, (getBoxData(a_srcs))...);
    protoEBForAllIrreg(a_F, getIrregData(a_srcs)...);
  }


  template<typename Func, typename... Srcs>
  inline void
  protoEBforall_i(const Func & a_F,  Box a_box, Srcs&... a_srcs)
  {
    forallInPlace_i(a_F, a_box, (getBoxData(a_srcs))...);
    protoEBForAllIrreg_i(a_F, getIrregData(a_srcs)...);
  }

///version that does not send the point to the function
  template<typename Func, typename... Srcs>
  inline void ebforallInPlace(unsigned long long int a_num_flops_point,
                              const char*            a_timername,
                              const Func & a_F,  Box a_box, Srcs&... a_srcs)
  {
    PR_TIME(a_timername);

    unsigned long long int boxfloops = a_num_flops_point*a_box.size();

    protoEBforall(a_F, a_box, a_srcs...);
    PR_FLOPS(boxfloops);
  }

///version that sends the point to the function
  template<typename Func, typename... Srcs>
  inline void ebforallInPlace_i(unsigned long long int a_num_flops_point,
                                const char*            a_timername,
                                const Func & a_F,  Box a_box, Srcs&... a_srcs)
  {
    PR_TIME(a_timername);

    unsigned long long int boxfloops = a_num_flops_point*a_box.size();
    protoEBforall_i(a_F, a_box, a_srcs...);

    PR_FLOPS(boxfloops);
  }

  ///version that only touches cut cells
  template<typename Func, typename... Srcs>
  inline void ebforallIrreg(const char*            a_timername,
                            const Func & a_F,  Box a_box, Srcs&... a_srcs)

  {
    PR_TIME(a_timername);
    //the flops thing here is tricky
    protoEBForAllIrreg(a_F, getIrregData(a_srcs)...);
  }

  ///version that only touches cut cells
  template<typename Func, typename... Srcs>
  inline void ebforallIrreg_i(const char*            a_timername,
                              const Func & a_F,  Box a_box, Srcs&... a_srcs)

  {
    PR_TIME(a_timername);
    //the flops thing here is tricky
    protoEBForAllIrreg_i(a_F, getIrregData(a_srcs)...);
  }


///
  template <typename T>
  inline T&
  checkBox(Box a_inputBox, T& a_T)
  {
    return a_T;
  }
//
  template <CENTERING cent, typename  data_t, unsigned int ncomp>
  inline  EBBoxData<cent, data_t, ncomp>&
  checkBox(Box a_box, EBBoxData<cent, data_t, ncomp>& a_s )
  {
    if(a_s.inputBox() != a_box)
    {
      using std::cout;
      using std::endl;                                 // 
      cout << "input box = " << a_box << ", data box = " << a_s.inputBox() << endl;
      PROTO_ASSERT(false, "boxes do not match in ebforall");
    }
    return a_s;
  }


///simplest version
  template<typename Func, typename... Srcs>
  inline void ebforall_i(const Box& a_inputBox,
                         const Func & a_F,  Box a_box, Srcs&... a_srcs)
  {
    unsigned long long int numflops = 0;

    ebforallInPlace_i(numflops, "ebforall_i", a_F, a_box, (checkBox(a_inputBox, a_srcs))...);
  }

///simplest version
  template<typename Func, typename... Srcs>
  inline void ebforall(const Box& a_inputBox,
                       const Func & a_F,  Box a_box, Srcs&... a_srcs)
  {
    unsigned long long int numflops = 0;

    ebforallInPlace(numflops, "ebforall", a_F, a_box, (checkBox(a_inputBox, a_srcs))...);
  }

  
}

#endif

