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

  run_perf_irreg_linear_full_args(500000);
  run_perf_irreg_linear_full_args(800000);
  run_perf_irreg_linear_full_args(700000);
  return 0;
}
