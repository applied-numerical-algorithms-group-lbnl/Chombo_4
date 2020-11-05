#include<EBProto.H>
#include<implem/Proto_AggStencil.H>
#include<implem/Proto_AggStencilImplem.H>
using Proto::EBDataLoc;
using Proto::pairPtr;
using Proto::pair_t;

inline void test_agg_stencil_alloc(EBDataLoc*& a_ptr, unsigned int a_size, int  &compt_0, int& compt_1)
{
  a_ptr = new EBDataLoc[a_size];
  compt_1 = 0;
  compt_0 = 0;
  for(int i = 0 ; i < a_size ; i++)
  {
    a_ptr[i].m_dataID = ((unsigned int)(i*3.15)) % 2;
    if(a_ptr[i].m_dataID == 1)
    {
      a_ptr[i].m_offset = compt_1;
      compt_1++;
    }
   
    if(a_ptr[i].m_dataID == 0)
    {
      a_ptr[i].m_offset = compt_0;
      compt_0++;
    }
  }  
}

template<typename T>
inline void test_agg_stencil_alloc(pairPtr<T>& a_data, unsigned int a_size)
{
  a_data.ptr[1] = new T[a_size];
  a_data.ptr[0] = new T[a_size];
}

inline void test_agg_stencil_alloc(uint64_t*& a_data, unsigned int a_size)
{
  a_data = new uint64_t[a_size];
}

inline void test_agg_stencil_fill_sten_size_start(uint64_t* a_sizes, uint64_t* a_start, unsigned int a_size)
{
  unsigned int n = 0;
  for(int i = 0 ; i < a_size ; i++)
  {
    a_start[i] = n;
    a_sizes[i] = 3*i+2;
    n += 3*i+2;
  }
}

inline void test_agg_stencil_fill_eb_stencil(pair_t<double>*& a_eb_stencil, uint64_t* a_sizes, int a_size, int a_s0, int a_s1)
{
  int acc = 0;
  for(int i = 0; i < a_size ; i++)
    acc += a_sizes[i];
 
  a_eb_stencil = new pair_t<double>[acc];
  for(int i = 0; i < acc ; i++)
  {
    EBDataLoc tmp;
    tmp.m_dataID =  i % 2;//((unsigned int)(6.21*i)) % 2;
    if(tmp.m_dataID == 0) tmp.m_offset = acc % a_s0;
    if(tmp.m_dataID == 1) tmp.m_offset = acc % a_s1;
    
    a_eb_stencil[i].first = tmp;
    a_eb_stencil[i].second = 1;
  }
}

inline void test_agg_stencil_fill_value(pairPtr<double>& a_data, double a_val0, double a_val1, unsigned int a_size)
{
  for(int i = 0 ; i < a_size ; i++)
  {
    a_data.ptr[0][i] = a_val0;
    a_data.ptr[1][i] = a_val1;
  }
}

template<typename T>
inline void test_agg_stencil_copy_to_gpu(T* a_host, T*& a_devi, unsigned int size)
{
  unsigned int nbBytes = size * sizeof(T);
  protoMalloc(&a_devi, nbBytes);
  protoMemcpy(a_devi,a_host,nbBytes,protoMemcpyHostToDevice);
}

template<typename T>
inline void test_agg_stencil_copy_to_gpu_eb_stencil(T* a_host, T*& a_devi, uint64_t* sizes, unsigned int size)
{
  int acc = 0;
  for(int i = 0; i<size ; i++) acc += sizes[i];
  unsigned int nbBytes = acc * sizeof(T);
  protoMalloc(&a_devi, nbBytes);
  protoMemcpy(a_devi,a_host,nbBytes,protoMemcpyHostToDevice);
}

template<typename T>
inline void test_agg_stencil_copy_to_gpu(pairPtr<T>& a_host, pairPtr<T>& a_devi, unsigned int size)
{
  unsigned int nbBytes = size * sizeof(T);
  protoMalloc(&a_devi.ptr[0], nbBytes);
  protoMemcpy(a_devi.ptr[0],a_host.ptr[0],nbBytes,protoMemcpyHostToDevice);
  protoMalloc(&a_devi.ptr[1], nbBytes);
  protoMemcpy(a_devi.ptr[1],a_host.ptr[1],nbBytes,protoMemcpyHostToDevice);
}

void test_agg_stencil_print(pairPtr<double>& a_data, unsigned int a_size)
{
  for(int i =0 ; i < a_size ; i++)
    std::cout << "(" << a_data.ptr[0][i] << "," << a_data.ptr[1][i] << ") ";
  std::cout << std::endl;
}

void test_agg_stencil_get_back_data(pairPtr<double>& a_host, pairPtr<double>& a_devi, unsigned int a_size)
{
  a_host.ptr[0] = new double[a_size];
  a_host.ptr[1] = new double[a_size];

  protoMemcpy(a_host.ptr[0],a_devi.ptr[0],a_size*sizeof(double),protoMemcpyDeviceToHost);
  protoMemcpy(a_host.ptr[1],a_devi.ptr[1],a_size*sizeof(double),protoMemcpyDeviceToHost);
  protoDeviceSynchronize();
}

bool test_answer_kernel_only_using(pairPtr<double> a_res, double a_val0, double a_val1, unsigned int a_size)
{
  int compt_0 = 0;
  int compt_1 = 0;
  for(int i = 0 ; i < a_size ; i++)
  {
    int id = ((unsigned int)(i*3.15)) % 2;
    if(id == 1)
    {
      int dep = 2*i +3*(i*(i-1)/2);
      double res = (int)((3*i+2)/2) *(a_val1 + a_val0);

      if(i%2==1)
      {
        if(dep%2==1) res += a_val1;
        else res += a_val0;
      }

      if(a_res.ptr[1][compt_1] != res)
      {
	std::cout << " a_res.ptr[1]["<<compt_1<<"]: " << a_res.ptr[1][compt_1] << " diff of " << res << std::endl;
	return false;
      }
      compt_1++;
    }

    if( id == 0)
    {
      int dep = 2*i +3*(i*(i-1)/2);
      double res = 0;
      res =  (int)((3*i+2)/2) *(a_val1 + a_val0);
      if(i%2==1)
      {
        if(dep%2==1) res += a_val1;
        else res += a_val0;
      }
      if(a_res.ptr[0][compt_0] != res) 
      {
	std::cout << "i: " << i << " and a_res.ptr[0]["<<compt_0<<"]: " << a_res.ptr[0][compt_0] << " diff of " << res << std::endl;
	return false;
      }
      compt_0++;
    }
  }

  for(int i = compt_0 ; i < a_size ; i++)
    if(a_res.ptr[0][i] != 0) return false;


  for(int i = compt_1 ; i < a_size ; i++)
    if(a_res.ptr[1][i] != 0) return false;

  return true;
}

bool test_answer_increment_only(pairPtr<double> a_res, double a_val0, double a_val1, double a_val2, double a_val3, unsigned int a_size)
{
  int compt_0 = 0;
  int compt_1 = 0;
  for(int i = 0 ; i < a_size ; i++)
  {
    int id = ((unsigned int)(i*3.15)) % 2;
    if(id == 1)
    {
      int dep = 2*i +3*(i*(i-1)/2);
      double res = (int)((3*i+2)/2) *(a_val1 + a_val0) + a_val3;

      if(i%2==1)
      {
        if(dep%2==1) res += a_val1;
        else res += a_val0;
      }

      if(a_res.ptr[1][compt_1] != res)
      {
	std::cout << " a_res.ptr[1]["<<compt_1<<"]: " << a_res.ptr[1][compt_1] << " diff of " << res << std::endl;
	return false;
      }
      compt_1++;
    }

    if( id == 0)
    {
      int dep = 2*i +3*(i*(i-1)/2);
      double res = a_val2;
      res +=  (int)((3*i+2)/2) *(a_val1 + a_val0);
      if(i%2==1)
      {
        if(dep%2==1) res += a_val1;
        else res += a_val0;
      }
      if(a_res.ptr[0][compt_0] != res) 
      {
	std::cout << "i: " << i << " and a_res.ptr[0]["<<compt_0<<"]: " << a_res.ptr[0][compt_0] << " diff of " << res << std::endl;
	return false;
      }
      compt_0++;
    }
  }

  for(int i = compt_0 ; i < a_size ; i++)
    if(a_res.ptr[0][i] != a_val2) return false;


  for(int i = compt_1 ; i < a_size ; i++)
    if(a_res.ptr[1][i] != a_val3) return false;

  return true;
}

bool test_agg_stencil_is_equal(pairPtr<double> a_gpu, double* a_cpu[2], unsigned int a_size)
{
  for(int j = 0 ; j < 2 ; j++)
    for(int i = 0 ; i < a_size ; i ++)
      if(a_gpu.ptr[j][i] != a_cpu[j][i]) return false;

  return true;
}

bool run_test_agg_stencil_kernel_only_using()
{
  unsigned int size_src = 8;
  unsigned int size_dst = 8;
  unsigned int begin = 0;
  unsigned int end = size_src;
  pairPtr<double> data_ptrs_src = {nullptr,nullptr};
  pairPtr<double> data_ptrs_dst = {nullptr,nullptr};
  pair_t<double>* eb_stencil = nullptr;
  uint64_t* sten_sizes   = nullptr;
  uint64_t* sten_start   = nullptr;
  EBDataLoc* dst_access      = nullptr;
  bool increment_only = false;
  double scale = 1;

  int nb0=0;
  int nb1=0;
  test_agg_stencil_alloc(dst_access,    size_dst, nb0, nb1);
  test_agg_stencil_alloc(data_ptrs_src, size_src); 
  test_agg_stencil_alloc(data_ptrs_dst, size_dst); 
  test_agg_stencil_alloc(sten_sizes,    size_src);  
  test_agg_stencil_alloc(sten_start,    size_src);    

  test_agg_stencil_fill_sten_size_start(sten_sizes, sten_start, size_src);
  test_agg_stencil_fill_eb_stencil(eb_stencil, sten_sizes, size_src, size_src, size_dst);
  test_agg_stencil_fill_value(data_ptrs_src, 10, 100, size_src);
  test_agg_stencil_fill_value(data_ptrs_dst,  0,   0, size_dst);

  pairPtr<double> data_ptrs_src_device = {nullptr,nullptr};
  pairPtr<double> data_ptrs_dst_device = {nullptr,nullptr};
  pair_t<double>* eb_stencil_device    = nullptr;
  uint64_t* sten_sizes_device          = nullptr;
  uint64_t* sten_start_device          = nullptr;
  EBDataLoc* dst_access_device         = nullptr;

  test_agg_stencil_copy_to_gpu(sten_sizes, sten_sizes_device, size_src);
  test_agg_stencil_copy_to_gpu(sten_start, sten_start_device, size_src);
  test_agg_stencil_copy_to_gpu(dst_access, dst_access_device, size_dst);
  test_agg_stencil_copy_to_gpu(data_ptrs_src, data_ptrs_src_device, size_src);
  test_agg_stencil_copy_to_gpu(data_ptrs_dst, data_ptrs_dst_device, size_dst);
  test_agg_stencil_copy_to_gpu_eb_stencil(eb_stencil, eb_stencil_device, sten_sizes, size_src);


  pairPtr<const double> const_data_ptrs_src_device = {data_ptrs_src_device.ptr[0],data_ptrs_src_device.ptr[1]};
  const pair_t<double>* const_eb_stencil_device = (const pair_t<double>*)(eb_stencil_device);


  protoLaunchKernel(Proto::aggStencilIndexer, 1, size_src,
			(int)(begin), 
			(int)(end),
			const_data_ptrs_src_device,
			data_ptrs_dst_device,
			const_eb_stencil_device,
			(const uint64_t*)(sten_sizes_device),
			(const uint64_t*)(sten_start_device),
			(const EBDataLoc*)(dst_access_device),
			increment_only,
			scale		
		);
  pairPtr<double> result = {nullptr,nullptr};
  test_agg_stencil_get_back_data(result, data_ptrs_dst_device, size_dst);

  bool check = test_answer_kernel_only_using(result, 10, 100, size_src);
  if(!check) 
  {
    test_agg_stencil_print(data_ptrs_src, size_src);
    test_agg_stencil_print(result, size_dst);
  }
  return check;
}

bool run_test_agg_stencil_increment_only_true()
{
  unsigned int size_src = 8;
  unsigned int size_dst = 8;
  unsigned int begin = 0;
  unsigned int end = size_src;
  pairPtr<double> data_ptrs_src = {nullptr,nullptr};
  pairPtr<double> data_ptrs_dst = {nullptr,nullptr};
  pair_t<double>* eb_stencil = nullptr;
  uint64_t* sten_sizes   = nullptr;
  uint64_t* sten_start   = nullptr;
  EBDataLoc* dst_access      = nullptr;
  bool increment_only = true;
  double scale = 1;

  int nb0=0;
  int nb1=0;
  test_agg_stencil_alloc(dst_access,    size_dst, nb0, nb1);
  test_agg_stencil_alloc(data_ptrs_src, size_src); 
  test_agg_stencil_alloc(data_ptrs_dst, size_dst); 
  test_agg_stencil_alloc(sten_sizes,    size_src);  
  test_agg_stencil_alloc(sten_start,    size_src);    

  test_agg_stencil_fill_sten_size_start(sten_sizes, sten_start, size_src);
  test_agg_stencil_fill_eb_stencil(eb_stencil, sten_sizes, size_src, size_src, size_dst);
  test_agg_stencil_fill_value(data_ptrs_src, 10, 100, size_src);
  test_agg_stencil_fill_value(data_ptrs_dst,  66,   666, size_dst);

  pairPtr<double> data_ptrs_src_device = {nullptr,nullptr};
  pairPtr<double> data_ptrs_dst_device = {nullptr,nullptr};
  pair_t<double>* eb_stencil_device    = nullptr;
  uint64_t* sten_sizes_device          = nullptr;
  uint64_t* sten_start_device          = nullptr;
  EBDataLoc* dst_access_device         = nullptr;

  test_agg_stencil_copy_to_gpu(sten_sizes, sten_sizes_device, size_src);
  test_agg_stencil_copy_to_gpu(sten_start, sten_start_device, size_src);
  test_agg_stencil_copy_to_gpu(dst_access, dst_access_device, size_dst);
  test_agg_stencil_copy_to_gpu(data_ptrs_src, data_ptrs_src_device, size_src);
  test_agg_stencil_copy_to_gpu(data_ptrs_dst, data_ptrs_dst_device, size_dst);
  test_agg_stencil_copy_to_gpu_eb_stencil(eb_stencil, eb_stencil_device, sten_sizes, size_src);


  pairPtr<const double> const_data_ptrs_src_device = {data_ptrs_src_device.ptr[0],data_ptrs_src_device.ptr[1]};
  const pair_t<double>* const_eb_stencil_device = (const pair_t<double>*)(eb_stencil_device);


  protoLaunchKernel(Proto::aggStencilIndexer, 1, size_src,
			(int)(begin), 
			(int)(end),
			const_data_ptrs_src_device,
			data_ptrs_dst_device,
			const_eb_stencil_device,
			(const uint64_t*)(sten_sizes_device),
			(const uint64_t*)(sten_start_device),
			(const EBDataLoc*)(dst_access_device),
			increment_only,
			scale		
		);
  pairPtr<double> result = {nullptr,nullptr};
  test_agg_stencil_get_back_data(result, data_ptrs_dst_device, size_dst);

  bool check = test_answer_increment_only(result, 10, 100, 66, 666, size_src);
  if(!check) 
  {
    test_agg_stencil_print(data_ptrs_src, size_src);
    test_agg_stencil_print(result, size_dst);
  }
  return check;
}

bool run_test_agg_stencil_cpu_versu_gpu()
{
  unsigned int size_src = 8;
  unsigned int size_dst = 8;
  unsigned int begin = 0;
  unsigned int end = size_src;
  pairPtr<double> data_ptrs_src = {nullptr,nullptr};
  pairPtr<double> data_ptrs_dst = {nullptr,nullptr};
  pair_t<double>* eb_stencil    = nullptr;
  uint64_t* sten_sizes          = nullptr;
  uint64_t* sten_start          = nullptr;
  EBDataLoc* dst_access         = nullptr;
  bool increment_only           = false;
  double scale                  = 1;

  int nb0=0;
  int nb1=0;
  test_agg_stencil_alloc(dst_access,    size_dst, nb0, nb1);
  test_agg_stencil_alloc(data_ptrs_src, size_src); 
  test_agg_stencil_alloc(data_ptrs_dst, size_dst); 
  test_agg_stencil_alloc(sten_sizes,    size_src);  
  test_agg_stencil_alloc(sten_start,    size_src);    

  test_agg_stencil_fill_sten_size_start(sten_sizes, sten_start, size_src);
  test_agg_stencil_fill_eb_stencil(eb_stencil, sten_sizes, size_src, size_src, size_dst);
  test_agg_stencil_fill_value(data_ptrs_src, 10, 100, size_src);
  test_agg_stencil_fill_value(data_ptrs_dst,  0,   0, size_dst);

  pairPtr<double> data_ptrs_src_device = {nullptr,nullptr};
  pairPtr<double> data_ptrs_dst_device = {nullptr,nullptr};
  pair_t<double>* eb_stencil_device    = nullptr;
  uint64_t* sten_sizes_device          = nullptr;
  uint64_t* sten_start_device          = nullptr;
  EBDataLoc* dst_access_device         = nullptr;

  test_agg_stencil_copy_to_gpu(sten_sizes, sten_sizes_device, size_src);
  test_agg_stencil_copy_to_gpu(sten_start, sten_start_device, size_src);
  test_agg_stencil_copy_to_gpu(dst_access, dst_access_device, size_dst);
  test_agg_stencil_copy_to_gpu(data_ptrs_src, data_ptrs_src_device, size_src);
  test_agg_stencil_copy_to_gpu(data_ptrs_dst, data_ptrs_dst_device, size_dst);
  test_agg_stencil_copy_to_gpu_eb_stencil(eb_stencil, eb_stencil_device, sten_sizes, size_src);


  pairPtr<const double> const_data_ptrs_src_device = {data_ptrs_src_device.ptr[0],data_ptrs_src_device.ptr[1]};
  const pair_t<double>* const_eb_stencil_device = (const pair_t<double>*)(eb_stencil_device);


  protoLaunchKernel(Proto::aggStencilIndexer, 1, size_src,
			(int)(begin), 
			(int)(end),
			const_data_ptrs_src_device,
			data_ptrs_dst_device,
			const_eb_stencil_device,
			(const uint64_t*)(sten_sizes_device),
			(const uint64_t*)(sten_start_device),
			(const EBDataLoc*)(dst_access_device),
			increment_only,
			scale		
		);
  pairPtr<double> result = {nullptr,nullptr};
  test_agg_stencil_get_back_data(result, data_ptrs_dst_device, size_dst);

  const double* const_data_ptrs_src_host[2] = {data_ptrs_src.ptr[0],data_ptrs_src.ptr[1]};
  double* data_ptrs_dst_host[2]             = {data_ptrs_dst.ptr[0],data_ptrs_dst.ptr[1]};
  const pair_t<double>* const_eb_stencil     = (const pair_t<double>*) (eb_stencil);

  bool check1 = test_answer_kernel_only_using(result, 10, 100, size_src);

  for(int idx = 0 ; idx < size_dst ; idx++)
  {
	Proto::hostAggStencilIndexer<double>(
			(int)(begin),
			(int)(end),
			idx,
			const_data_ptrs_src_host,
			data_ptrs_dst_host,
			const_eb_stencil,
			(const uint64_t*)(sten_sizes),
			(const uint64_t*)(sten_start),
			(const EBDataLoc*)(dst_access),
			increment_only,
			scale
			);
  }

  bool check2 = test_agg_stencil_is_equal(result, data_ptrs_dst_host, size_dst);

  return check1 && check2;
}


bool run_test_agg_stencil_empty()
{
/*  Proto::AggStencil<Proto::CELL, Proto::CELL, double> test;

  bool check1 = test.d_ebstencil == nullptr; 
  bool check2 = test.d_dstaccess == nullptr; 
  bool check3 = test.d_stensizes == nullptr; 
  bool check4 = test.d_stenstart == nullptr; 
  bool check5 = test.m_stride == 0; 
  bool check6 = test.m_blocks == 0; 

  bool check = check1 && check2 && check3 && check4 && check5 && check6;
  return check;
*/
  return false;
}
