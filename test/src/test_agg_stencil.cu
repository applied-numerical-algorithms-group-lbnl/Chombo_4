#include<EBProto.H>
#include<implem/Proto_AggStencil.H>
#include<implem/Proto_AggStencilImplem.H>
using Proto::EBDataLoc;
//using Proto::pairPtr;
using Proto::pair_t;
#include<src/test_agg_stencil_stuff.cpp>

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


bool run_test_agg_stencil_increment_only(bool b)
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
  bool increment_only = b;
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

  bool check = test_answer_increment_only(result, 10, 100, 66, 666, size_src, increment_only);
  if(!check) 
  {
    test_agg_stencil_print(data_ptrs_src, size_src);
    test_agg_stencil_print(result, size_dst);
  }
  return check;
}

bool run_test_agg_stencil_increment_only_true()
{
  return run_test_agg_stencil_increment_only(true);
}

bool run_test_agg_stencil_increment_only_false()
{
  return run_test_agg_stencil_increment_only(false);
}
//purge
/*
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
*/

bool run_test_agg_stencil_scale(double a_scale)
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
  double scale = a_scale;

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

  bool check = test_agg_stencil_answer_scale(result, 10, 100, size_src, a_scale);
  if(!check) 
  {
    test_agg_stencil_print(data_ptrs_src, size_src);
    test_agg_stencil_print(result, size_dst);
  }
  return check;
}

bool run_test_agg_stencil_scale_0()
{
  double a = 0;
  return run_test_agg_stencil_scale(a);
}

bool run_test_agg_stencil_scale_minus10()
{
  double a = -10;
  return run_test_agg_stencil_scale(a);
}

bool run_test_agg_stencil_scale_100()
{
  double a = 100;
  return run_test_agg_stencil_scale(a);
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
