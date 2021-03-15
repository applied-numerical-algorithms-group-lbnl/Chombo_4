#include<iostream>
#include<EBProto.H>
#include<src/test_perf_dtod_kernels.cu>


int main()
{
#ifdef PROTO_CUDA
  protoSetDevice(0);
#endif
  std::cout << " perf tests " << std::endl;

  run_perf_irreg_data_copy_args(256);
  run_perf_irreg_data_copy_args(10000);
  run_perf_irreg_data_copy_args(20000);
  run_perf_irreg_data_copy_args(30000);
  run_perf_irreg_data_copy_args(40000);
  run_perf_irreg_data_copy_args(50000);
  run_perf_irreg_linear_full_args(10000);
  run_perf_irreg_linear_full_args(20000);
  run_perf_irreg_linear_full_args(30000);
  run_perf_irreg_linear_full_args(40000);
  run_perf_irreg_linear_full_args(50000);

  run_perf_irreg_linear_partial_args(1000,10000);
  run_perf_irreg_linear_partial_args(500,10000);
  run_perf_irreg_linear_partial_args(100,10000);
  run_perf_irreg_linear_partial_args(50,10000);
  run_perf_irreg_linear_partial_args(20,10000);
  run_perf_irreg_linear_partial_args(10,10000);


 
  run_perf_irreg_linear_partial_graph_args(1000,10000);
  run_perf_irreg_linear_partial_graph_args(500,10000);
  run_perf_irreg_linear_partial_graph_args(100,10000);
  run_perf_irreg_linear_partial_graph_args(50,10000);
  run_perf_irreg_linear_partial_graph_args(20,10000);
  run_perf_irreg_linear_partial_graph_args(10,10000);

/*
  run_perf_irreg_linear_partial_args(10000,100000);
  run_perf_irreg_linear_partial_args(1000,100000);
  run_perf_irreg_linear_partial_args(500,100000);
  run_perf_irreg_linear_partial_args(100,100000);
  run_perf_irreg_linear_partial_args(50,100000);
  run_perf_irreg_linear_partial_args(20,100000);
  run_perf_irreg_linear_partial_args(10,100000);
*/
  return 0;
}
