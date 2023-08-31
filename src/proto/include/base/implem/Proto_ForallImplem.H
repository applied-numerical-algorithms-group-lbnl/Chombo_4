
// Helper API for dynamically determining MemType
template<typename T>
struct getMemType
{
    static constexpr MemType type()
    {
        return Proto::MemType::BOTH;
    }
    static MemType type_eval()
    {
        return Proto::MemType::BOTH;
    }
};

template<typename T, unsigned int C, MemType MEMTYPE,unsigned char D, unsigned char E>
struct getMemType<BoxData<T,C,MEMTYPE,D,E>>
{
  static constexpr MemType type()
  {
    return MEMTYPE;
  }
  static MemType type_eval()
  {
	  return MEMTYPE;
  }
};

template<typename T, unsigned int C, MemType MEMTYPE,unsigned char D, unsigned char E>
struct getMemType<Var<T,C,MEMTYPE,D,E>>
{
  static constexpr MemType type()
  {
    return MEMTYPE;
  }
  static MemType type_eval()
  {
	  return MEMTYPE;
  }
};

template<typename First>
constexpr MemType getMemTypeFromSrcs()
{
  return getMemType<First>::type();
 
}

template<typename First, typename Second, typename... Rest>
constexpr MemType getMemTypeFromSrcs()
{
  return (MemType)(getMemType<First>::type() & getMemTypeFromSrcs<Second, Rest...>());
}

template <typename T, unsigned int C, MemType MEMTYPE, unsigned char D, unsigned char E>
void forall_parse(Box& a_domain, const BoxData<T,C,MEMTYPE,D,E>& a_data)
{
  // Can do anything you want here later
  if(a_domain.empty())
  {
    a_domain=a_data.box();
  }
  else
  {
    a_domain = a_domain & a_data.box();
  }
}

////
template <typename T>
void forall_parse(Box& a_domain, T a_scalar)
{
  // This is a dummy function; it's meant to be empty
}

////
template <typename Last>
void forall_parse_args(Box& a_domain, Last& a_last)
{
  forall_parse(a_domain, a_last);
}

////
template <typename T, unsigned int C, unsigned char D, unsigned char E, MemType MEMTYPE, typename... Rest>
void forall_parse_args(Box& a_domain, const BoxData<T,C,MEMTYPE,D,E>& a_first, Rest&&... a_rest)
{
  forall_parse(a_domain, a_first);
  forall_parse_args(a_domain, std::forward<Rest>(a_rest)...);
}

////
template <typename T, unsigned int C, MemType MEMTYPE, unsigned char D, unsigned char E>
void forall_parse_args(Box& a_domain, const BoxData<T,C,MEMTYPE,D,E>& a_first)
{
  forall_parse(a_domain, a_first);
}

////
template <typename First , typename... Rest>
void forall_parse_args(Box& a_domain, First& a_first, Rest&&... a_rest)
{
  forall_parse(a_domain, a_first);
  forall_parse_args(a_domain, std::forward<Rest>(a_rest)...);
}

//========================================================================
// Base Forall functionality
//========================================================================

// Main version
template<typename Func, typename... Srcs>
void forallInPlaceBaseOp(unsigned long long int a_num_flops_point,
                         const char* a_timername,
                         const Func& a_F, Box a_box, Srcs&&... a_srcs)
{
  PR_TIME(a_timername);
  forallInPlaceBase(a_F, a_box, std::forward<Srcs>(a_srcs)...);
  unsigned long long int boxfloops = a_num_flops_point*a_box.size();
  PR_FLOPS(boxfloops);
}

// Main version
template<typename Func, typename... Srcs>
void forallInPlaceBase(const Func& a_F, Box a_box, Srcs&&... a_srcs)
{
  PR_TIME("forallInPlaceBase");
  protoForall(a_F, a_box, std::forward<Srcs>(a_srcs)...);
}

template<typename Func, typename... Srcs>
void forallInPlaceBaseOp_p(unsigned long long int a_num_flops_point,
                           const char* a_timername,
                           const Func& a_F, Box a_box, Srcs&&... a_srcs)
{
  PR_TIME(a_timername);
  protoForall_p(a_F, a_box, std::forward<Srcs>(a_srcs)...) ;

  unsigned long long int  boxfloops = a_num_flops_point*a_box.size();
  PR_FLOPS(boxfloops);
}


template<typename Func, typename... Srcs>
void forallInPlace_i(const Func& a_F, Box a_box, Srcs&&... a_srcs)
{
  protoForall_i(a_F, a_box, std::forward<Srcs>(a_srcs)...) ;
}

template<typename Func, typename... Srcs>
void forallInPlaceBaseOp_i(unsigned long long int a_num_flops_point,
                           const char* a_timername,
                           const Func& a_F, Box a_box, Srcs&&... a_srcs)
{
  PR_TIME(a_timername);

  forallInPlace_i(a_F, a_box, std::forward<Srcs>(a_srcs)...);

  unsigned long long int boxfloops = a_num_flops_point*a_box.size();
  PR_FLOPS(boxfloops);
}



// Point input version
template<typename Func, typename... Srcs>
inline void forallInPlaceBase_p(const Func& a_F,  Box a_box, Srcs&&... a_srcs)
{
  PR_TIME("forallInPlaceBase_p");
  return protoForall_p(a_F, a_box, std::forward<Srcs>(a_srcs)...) ;
}

//========================================================================
//  Frontend Forall Implementation
//========================================================================

template<typename T, unsigned int C, unsigned char D, unsigned char E, MemType MEMTYPE,
         typename Func, typename... Srcs>
BoxData<T,C,MEMTYPE,D,E> forallOp(unsigned long long int a_num_flops_point,   // Main ForallOp
                                  const char* a_timername,
                                  const Func& a_F, Srcs&&... a_srcs)
{
  Box intersect;
  forall_parse_args(intersect,std::forward<Srcs>(a_srcs)...);

  PROTO_ASSERT(intersect.size() > 0,
               "forall(...) inferred the domain Box to be empty.");

  BoxData<T,C,MEMTYPE,D,E> ret(intersect);
  forallInPlaceBaseOp(a_num_flops_point, a_timername, a_F, intersect, ret, std::forward<Srcs>(a_srcs)...);
  return ret;
}

// Main Forall
template<typename T, unsigned int C, unsigned char D, unsigned char E, MemType MEMTYPE,
         typename Func, typename... Srcs>
BoxData<T,C,MEMTYPE,D,E> forall(const Func& a_F, Srcs&&... a_srcs)
{
  Box intersect;
  forall_parse_args(intersect,std::forward<Srcs>(a_srcs)...);

  PROTO_ASSERT(intersect.size() > 0,
               "forall(...) inferred the domain Box to be empty.");

  BoxData<T,C,MEMTYPE,D,E> ret(intersect);
  forallInPlaceBase(a_F, intersect, ret, std::forward<Srcs>(a_srcs)...);
  return ret;
}

// Overload forallOp with Optional Box Argument
template<typename T, unsigned int C, unsigned char D, unsigned char E, MemType MEMTYPE,
         typename Func, typename... Srcs>
BoxData<T,C,MEMTYPE,D,E> forallOp(unsigned long long int a_num_flops_point,
                          const char* a_timername,
                          const Func& a_F, Box a_box, Srcs&&... a_srcs)
{
  BoxData<T,C,MEMTYPE,D,E> ret(a_box);
  forallInPlaceBaseOp(a_num_flops_point, a_timername, a_F, a_box, ret, std::forward<Srcs>(a_srcs)...);
  return ret;
}


// Overload forall with Optional Box Argument
template<typename T, unsigned int C, unsigned char D, unsigned char E, MemType MEMTYPE,
         typename Func, typename... Srcs>
BoxData<T,C,MEMTYPE,D,E> forall(const Func& a_F, Box a_box, Srcs&&... a_srcs)
{
  BoxData<T,C,MEMTYPE,D,E> ret(a_box);
  forallInPlaceBase(a_F, a_box, ret, std::forward<Srcs>(a_srcs)...);
  return ret;
}

// In Place forallOp Version
template<typename Func, typename... Srcs>
void forallInPlaceOp(unsigned long long int a_num_flops_point,
                     const char* a_timername,
                     const Func& a_F, Srcs&&... a_srcs)
{
  Box intersect;
  forall_parse_args(intersect, std::forward<Srcs>(a_srcs)...);

  PROTO_ASSERT(intersect.size() > 0,
               "forall(...) inferred the domain Box to be empty.");

  forallInPlaceBaseOp(a_num_flops_point, a_timername, a_F, intersect, std::forward<Srcs>(a_srcs)...);
}

// In Place forall Version
template<typename Func, typename... Srcs>
void forallInPlace(const Func& a_F, Srcs&&... a_srcs)
{
  Box intersect;
  forall_parse_args(intersect, std::forward<Srcs>(a_srcs)...);

  PROTO_ASSERT(intersect.size() > 0,
               "forall(...) inferred the domain Box to be empty.");

  forallInPlaceBase(a_F, intersect, std::forward<Srcs>(a_srcs)...);
}

// In Place forallop Version with Box Argument
template<typename Func, typename... Srcs>
void forallInPlaceOp(unsigned long long int a_num_flops_point,
                     const char* a_timername,
                     const Func& a_F, Box a_box, Srcs&&... a_srcs)
{
  forallInPlaceBaseOp(a_num_flops_point, a_timername, a_F, a_box, std::forward<Srcs>(a_srcs)...);
}

// In Place forall Version with Box Argument
template<typename Func, typename... Srcs>
void forallInPlace(const Func& a_F, Box a_box, Srcs&&... a_srcs)
{
  forallInPlaceBase(a_F, a_box, std::forward<Srcs>(a_srcs)...);
}

// Main forallOp version with Point Access
template<typename T, unsigned int C, unsigned char D, unsigned char E, MemType MEMTYPE,
         typename Func, typename... Srcs>
BoxData<T,C,MEMTYPE,D,E> forallOp_p(unsigned long long int a_num_flops_point,
                                    std::string a_timername,
                                    const Func& a_F, Srcs&&... a_srcs)
{
  Box intersect;
  forall_parse_args(intersect, std::forward<Srcs>(a_srcs)...);

  PROTO_ASSERT(intersect.size() > 0, "forall(...) inferred the domain Box to be empty.");

  BoxData<T,C,MEMTYPE,D,E> ret(intersect);
  forallInPlaceBaseOp_p(a_num_flops_point, a_timername, a_F, intersect, ret, std::forward<Srcs>(a_srcs)...);
  return ret;
}

// Main forall version with Point Access
template<typename T, unsigned int C, unsigned char D, unsigned char E, MemType MEMTYPE,
         typename Func, typename... Srcs>
BoxData<T,C,MEMTYPE,D,E> forall_p(const Func& a_F, Srcs&&... a_srcs)
{
  Box intersect;
  forall_parse_args(intersect, std::forward<Srcs>(a_srcs)...);

  PROTO_ASSERT(intersect.size() > 0, "forall(...) inferred the domain Box to be empty.");

  BoxData<T,C,MEMTYPE,D,E> ret(intersect);
  forallInPlaceBase_p(a_F, intersect, ret, std::forward<Srcs>(a_srcs)...);
  return ret;
}

// Main forallOp version with Point Access and Box Argument
template<typename T, unsigned int C, unsigned char D, unsigned char E, MemType MEMTYPE,
         typename Func, typename... Srcs>
BoxData<T,C,MEMTYPE,D,E> forallOp_p(unsigned long long int a_num_flops_point,
                                    const char* a_timername,
                                    const Func& a_F,  Box a_box, Srcs&&... a_srcs)
{
  BoxData<T,C,MEMTYPE,D,E> ret(a_box, true);
  forallInPlaceBaseOp_p(a_num_flops_point, a_timername, a_F, a_box, ret, std::forward<Srcs>(a_srcs)...);
  return ret;
}

// Main forall version with Point Access and Box Argument
template<typename T, unsigned int C, unsigned char D, unsigned char E, MemType MEMTYPE,
         typename Func, typename... Srcs>
BoxData<T,C,MEMTYPE,D,E> forall_p(const Func& a_F,  Box a_box, Srcs&&... a_srcs)
{
  BoxData<T,C,MEMTYPE,D,E> ret(a_box, true);
  forallInPlaceBase_p(a_F, a_box, ret, std::forward<Srcs>(a_srcs)...);
  return ret;
}

// In Place forallOp version with Point Access
template<typename Func, typename... Srcs>
inline void forallInPlaceOp_p(unsigned long long int a_num_flops_point,
                              const char* a_timername,
                              const Func& a_F, Srcs&&... a_srcs)
{
  Box intersect;
  forall_parse_args(intersect, std::forward<Srcs>(a_srcs)...);

  PROTO_ASSERT(intersect.size() > 0,
               "forall(...) inferred the domain Box to be empty.");
  forallInPlaceBaseOp_p(a_num_flops_point, a_timername, a_F, intersect, std::forward<Srcs>(a_srcs)...);
}

// In Place forall version with Point Access
template<typename Func, typename... Srcs>
inline void forallInPlace_p(const Func& a_F, Srcs&&... a_srcs)
{
  Box intersect;
  forall_parse_args(intersect, std::forward<Srcs>(a_srcs)...);

  PROTO_ASSERT(intersect.size() > 0,
               "forall(...) inferred the domain Box to be empty.");
  forallInPlaceBase_p(a_F, intersect, std::forward<Srcs>(a_srcs)...);
}

// In Place forallOp version with Point Access and Box Argument
template<typename Func, typename... Srcs>
inline void forallInPlace_p(unsigned long long int a_num_flops_point,
                            const char* a_timername,
                            const Func& a_F,  Box a_box, Srcs&&... a_srcs)
{
  forallInPlaceBase_p(a_num_flops_point, a_timername, a_F, a_box, std::forward<Srcs>(a_srcs)...);
}


// In Place forall version with Point Access and Box Argument
template<typename Func, typename... Srcs>
inline void forallInPlace_p(const Func& a_F,  Box a_box, Srcs&&... a_srcs)
{
  forallInPlaceBase_p(a_F, a_box, std::forward<Srcs>(a_srcs)...);
}


template<typename FuncStruct>
struct emptyIndexer
{
#ifdef PROTO_CUDA
  __device__ static 
  void gpu(FuncStruct a_body)
  {
  }
#endif
  static inline void cpu(FuncStruct a_body){}
};

template<typename... T>
void emptyT(T... in)
{
}


template<typename FuncStruct, typename... Srcs>
struct structIndexer
{
  #ifdef PROTO_CUDA
  __device__ static 
  void gpu(int a_begin, int a_end, FuncStruct a_body, Srcs... a_srcs)
  {
    int idx = threadIdx.x + blockIdx.x*blockDim.x;
    if (idx >= a_begin && idx < a_end)
      {
        a_body.op(a_srcs...);
      }
  }
  #endif
  
  static void cpu(int begin, int end, FuncStruct a_body, Srcs... a_srcs)
  {
    for(int i=begin;i<end; i++)
      {
        a_body.op(a_srcs...);
        emptyT(var_incr(a_srcs)...);
      }
  }
};

template<typename T>
inline T p_ref(const T& a_s, const Point& a_p) {return a_s;}

template<typename T, unsigned int C, MemType MEMTYPE, unsigned char D, unsigned char E>
inline Var<T,C,MEMTYPE,D,E> p_ref(Var<T,C,MEMTYPE,D,E>&  a_data, const Point& a_p){  a_data += a_p; return a_data;}

template<typename T, unsigned int C, MemType MEMTYPE, unsigned char D, unsigned char E>
inline Var<T,C,MEMTYPE,D,E> p_ref(const Var<T,C,MEMTYPE,D,E>& a_data, const Point& a_p){a_data += a_p; return a_data;}

template<typename T>
inline T var_incr(const T& t){return t;}

template<typename T, unsigned int C, MemType MEMTYPE, unsigned char D, unsigned char E>
Var<T, C, MEMTYPE, D, E>& var_incr(Var<T, C, MEMTYPE, D, E>& t) {return ++t;}


template<typename Func, typename... T>
inline void pencilFunc(const Func& F, int count, T... vars)
{
  for(int i=0; i < count; ++i)
  {
    F(vars...);
    emptyT(var_incr(vars)...);
  }
}

template<typename Func, typename... T>
inline void pencilFunc_p(const Func& F, int count, Point a_pt, T... vars)
{
  for(int i=0; i < count; ++i, ++a_pt[0])
  {
    F(a_pt,vars...);
    emptyT(var_incr(vars)...);
  }
}

template<typename Func, typename... T>
inline void pencilFunc_i(const Func& F, int count, Point a_pt, T... vars)
{
  for(int i=0; i < count; ++i, ++a_pt[0])
  {
    F(a_pt.m_tuple,vars...);
    emptyT(var_incr(vars)...);
  }
}

template<typename Func, typename... Srcs>
struct indexer {

  static void cpu(int a_begin, int a_end, const Box& a_box, Func a_body, Srcs... a_srcs)
  {
    Box cross = a_box.flatten(0);
    int npencil = a_box.size(0);
    
    auto last = a_box.low();
    
    for (auto iter = cross.begin(); iter != cross.end(); ++iter)
      {
        pencilFunc(a_body, npencil, p_ref(a_srcs,*iter-last)...);
        last = *iter;
      }
  }
#ifdef PROTO_CUDA
  __device__ static
  void gpu(int a_begin, int a_end, const Box& a_box, Func a_body, Srcs... a_srcs)
  {
    int idx = threadIdx.x + blockIdx.x*blockDim.x;
    if (idx >= a_begin && idx < a_end)
      {
        a_body(a_srcs...);
      }
  }
#endif
};

template<typename Func, typename... Srcs>
struct indexer_p {
  
  static void cpu(int a_begin, int a_end, const Box& a_box,
                  Func a_body, Srcs... a_srcs)
  {
    
    Box cross = a_box.flatten(0);
    int npencil = a_box.size(0);

    auto last = a_box.low();
    for (auto iter = cross.begin(); iter != cross.end(); ++iter)
      {
        pencilFunc_p(a_body, npencil, *iter, (p_ref(a_srcs,*iter-last))...);
        last = *iter;
      }
  }

  #ifdef PROTO_CUDA
  __device__ static
  void gpu(int a_begin, int a_end, const Box& a_box,
           Func a_body, Srcs... a_srcs)
  {
    int idx = threadIdx.x + blockIdx.x*blockDim.x;
    int idy = threadIdx.y + blockIdx.y*blockDim.y;
    int idz = threadIdx.z + blockIdx.z*blockDim.z;
    if (idx >= a_begin && idx < a_end)
      {
        Point p = a_box(idx,idy,idz);
        a_body(p,a_srcs...);
      }
  }
#endif
  
};

template<typename Func, typename... Srcs>
struct indexer_i {
  static void cpu(int a_begin, int a_end, const Box& a_box,
                  Func a_body, Srcs... a_srcs)
  {
    Box cross = a_box.flatten(0);
    int npencil = a_box.size(0);
    auto last = a_box.low();
    
    for (auto iter = cross.begin(); iter != cross.end(); ++iter)
      {
        pencilFunc_i(a_body, npencil, *iter, p_ref(a_srcs,*iter-last)...);
        last = *iter;
      }
  }

#ifdef PROTO_CUDA
  __device__ static
  void gpu(int a_begin, int a_end, const Box& a_box,
           Func a_body, Srcs... a_srcs)
  {
    int idx = threadIdx.x + blockIdx.x*blockDim.x;
    int idy = threadIdx.y + blockIdx.y*blockDim.y;
    int idz = threadIdx.z + blockIdx.z*blockDim.z;
    
    if (idx >= a_begin && idx < a_end)
      {
        Point p = a_box(idx,idy,idz);
        a_body(p.m_tuple,a_srcs...);
      }
  }
#endif
};



inline void v100tuningBox(Box a_box, dim3& blocks, int& blockSize)
 {
 #if DIM == 3
   if(a_box.size(0)<=512)
   {
    blocks = dim3(1,a_box.size(1),a_box.size(2));
    blockSize = a_box.size(0);
  }
  else
  {
    blockSize = 512;
    blocks = dim3((a_box.size(0)+blockSize-1)/blockSize, a_box.size(1),a_box.size(2));
  }

#else
  if(a_box.size(0)<=512)
  {
    blocks = dim3(1,a_box.size(1),1);
    blockSize = a_box.size(0);
  }
  else
  {
    blockSize = 512;
    blocks = dim3((a_box.size(0)+blockSize-1)/blockSize, a_box.size(1),1);
  }
#endif
}

inline void v100tuningBox(int a_stride, Box a_cross, dim3& blocks, int& blockSize)
{

#if DIM == 3
  if(a_stride<=512)
  {
    blocks = dim3(1,a_cross.size(1),a_cross.size(2));
    blockSize = a_stride;
  }
  else
  {
    blockSize = 512;
    blocks = dim3((a_stride+blockSize-1)/blockSize, a_cross.size(1),a_cross.size(2));
  }

#else
  if(a_stride<=512)
  {
    blocks = dim3(1,a_cross.size(1),1);
    blockSize = a_stride;
  }
  else
  {
    blockSize = 512;
    blocks = dim3((a_stride+blockSize-1)/blockSize, a_cross.size(1),1);
  }
#endif
}



template<typename Func, typename... Srcs>
inline void makeVars(const Func& a_F, Box a_box, Srcs&&... a_srcs)
{
  PR_TIME("BoxData::makevars");
  int N = a_box.size();
  int stride = a_box.size(0);

#if DIM == 3
  dim3 blocks(1,a_box.size(1),a_box.size(2));
#else
  dim3 blocks(1,a_box.size(1),1);
#endif
  {
    PR_TIME("indexer");
    size_t smem = 0;

#ifdef superDebug
    std::cout << "Try Kernel Name is: " << a_F.myname << std::endl;
#endif

    assert(stride<1024);

    protoLaunchKernelMemAsyncT<getMemTypeFromSrcs<Srcs...>(),indexer<Func, Srcs...>>( blocks, stride, smem, protoGetCurrentStream, 0, N, a_box, a_F, std::forward<Srcs>(a_srcs)...);


#ifdef superDebug
    std::cout << "Success Kernel Name is: " << a_F.myname << std::endl;
#endif
  }
}

template<typename Func, typename... Srcs>
inline void makeVarsStream(protoStream_t& a_stream, const Func& a_F, Box a_box, Srcs&&... a_srcs)
{
  PR_TIME("BoxData::makevarsstream");
  int N = a_box.size();
  int stride = a_box.size(0);
  size_t smem = 0;
#if DIM == 3
  dim3 blocks(1,a_box.size(1),a_box.size(2));
#else
  dim3 blocks(1,a_box.size(1),1);
#endif
  {
    PR_TIME("indexer");

    assert(stride<1024);
    protoLaunchKernelMemAsyncT<getMemTypeFromSrcs<Srcs...>(), indexer<Func, Srcs...>>(blocks, stride, smem, a_stream, 0, N, a_box, a_F, std::forward<Srcs>(a_srcs)...);
  }
}


template<typename FuncStruct, typename... Srcs>
inline void makeVarsStruct(protoStream_t& a_stream, const FuncStruct& a_F, Box a_box, Srcs&&... a_srcs)
{
  PR_TIME("BoxData::makevarsstruct");
  int N = a_box.size();
  int stride = a_box.size(0);
#if DIM == 3
  dim3 blocks(1,a_box.size(1),a_box.size(2));
#else
  dim3 blocks(1,a_box.size(1),1);
#endif
  size_t smem = 0;
  {
    PR_TIME("structIndexer");

    assert(stride<1024);
    protoLaunchKernelMemAsyncT<getMemTypeFromSrcs<Srcs...>,structIndexer<FuncStruct, Srcs...>>( blocks, stride, smem, a_stream, 0, N, a_F, std::forward<Srcs>(a_srcs)...);
  }
}

template<typename FuncStruct, typename... Srcs>
inline void makeVarsEmptyIndexer(protoStream_t& a_stream, const FuncStruct& a_F, Box a_box, Box a_srcBox)
{
  //a_srcs is coming into here as vars
  PR_TIME("BoxData::makevarsEmptyIndexer");

  int stride = a_box.size(0);
  int blocks = a_box.size(1);
  /*
#if DIM==3
  int Nplane    = a_box.flatten(2).size();
  int Nz        = a_box.size(2);
  int zinc      = a_srcBox.flatten(2).size();
#else
  int Nplane    = a_box.size();
  int Nz        = 1;
  int zinc      = 1; //has to be non-zero or we have an infinite loop
#endif
  */
  size_t smem = 0;
  {
    PR_TIME("zincStructIndexer");
    //Doesn't matter what memtype we pass in because this is an empyt kernel
    assert(stride<1024);
    protoLaunchKernelMemAsyncT<getMemTypeFromSrcs<Srcs...>,emptyIndexer<FuncStruct>>( blocks, stride, smem, a_stream, a_F);
  }
}

template<typename Func, typename... Srcs>
inline void makeVars_p(const Func& a_F, Box a_box, Srcs&&... a_srcs)
{
  PR_TIME("BoxData::makevars_p");
  int N = a_box.size();
  int stride = a_box.size(0);
  Box cross = a_box.flatten(0);
//  int blocks =  cross.size();
#if DIM == 3
  dim3 blocks(1,cross.size(1),cross.size(2));
#else
  dim3 blocks(1,cross.size(1),1);
#endif

  {
    PR_TIME("indexer_p");
    protoLaunchKernelT<getMemTypeFromSrcs<Srcs...>(),indexer_p<Func, Srcs...>>(blocks, stride, 0, N, a_box,
                      a_F, std::forward<Srcs>(a_srcs)...);
  }
}



template<typename Func, typename... Srcs>
inline void makeVars_i(const Func& a_F, Box a_box, Srcs&&... a_srcs)
{
  PR_TIME("BoxData::makevars_i");
  int N = a_box.size();
  int stride = a_box.size(0);
  Box cross = a_box.flatten(0);
#if DIM == 3
  dim3 blocks(1,cross.size(1),cross.size(2));
#else
  dim3 blocks(1,cross.size(1),1);
#endif
  {
    PR_TIME("indexer_i");
    protoLaunchKernelT<getMemTypeFromSrcs<Srcs...>(),indexer_i<Func, Srcs...>>(blocks, stride, 0, N, a_box,
                      a_F, std::forward<Srcs>(a_srcs)...);
  }
}


template<typename T>
inline T p_ref_cuda(const T& a_s, const Point& a_p) {return a_s;}


template<typename T, unsigned int C, MemType MEMTYPE, unsigned char D, unsigned char E>
inline Var<T,C,MEMTYPE,D,E> p_ref_cuda(BoxData<T,C,MEMTYPE,D,E>& a_data, const Point& a_p)
{
  return a_data.var(a_p);
}

template<typename T, unsigned int C, MemType MEMTYPE, unsigned char D, unsigned char E>
inline Var<T,C,MEMTYPE,D,E> p_ref_cuda(const BoxData<T,C,MEMTYPE,D,E>& a_data, const Point& a_p)
{
  return a_data.var(a_p);
}



template<typename Func, typename... Srcs>
inline void protoForall(const Func& a_F, Box a_box, Srcs&&... a_srcs)
{
  PR_TIME("protoForall");
  makeVars(a_F, a_box, p_ref_cuda(a_srcs, a_box.low())...);
}

template<typename Func, typename... Srcs>
inline void protoForallStream(protoStream_t& a_stream, const Func& a_F, Box a_box, Srcs&... a_srcs)
{
  PR_TIME("protoForallStream");
  makeVarsStream(a_stream, a_F, a_box, p_ref_cuda(a_srcs, a_box.low())...);
}

template<typename FuncStruct, typename... Srcs>
inline void protoForallStruct(protoStream_t& a_stream, const FuncStruct& a_F, Box a_box, Srcs&... a_srcs)
{
  PR_TIME("protoForallStruct");
  makeVarsStruct(a_stream, a_F, a_box, p_ref_cuda(a_srcs, a_box.low())...);
}

template<typename FuncStruct, typename... Srcs>
inline void protoForallEmptyIndexer(protoStream_t& a_stream, const FuncStruct& a_F, Box a_box, Box a_srcBox)
{
  PR_TIME("protoForallEmptyIndexer");
  makeVarsEmptyIndexer(a_stream, a_F, a_box, a_srcBox);
}


template<typename Func, typename... Srcs>
inline void protoForall_p(const Func& a_F, Box a_box, Srcs&&... a_srcs)
{
  PR_TIME("protoForall_p");
  makeVars_p(a_F, a_box, p_ref_cuda(a_srcs, a_box.low())...);
}


template<typename Func, typename... Srcs>
inline void protoForall_i(const Func& a_F, Box a_box, Srcs&&... a_srcs)
{
  PR_TIME("protoForall_i");
  makeVars_i(a_F, a_box, p_ref_cuda(a_srcs, a_box.low())...);
}
