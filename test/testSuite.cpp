#include<iostream>
#include<src/test_IrregData.cu>
#include<src/test_ebforall.cu>
#include<src/test_ebforall_i.cu>
#include<src/test_agg_stencil.cu>

template<typename Func>
void do_test(std::string a_str, Func &fun)
{
  std::cout << "Do " << a_str << std::endl;
  bool b = fun();
  if(b) std::cout << "-> passed " << a_str << std::endl;
  else std::cout << "-> failed " << a_str << std::endl;
}
int main()
{
  cudaSetDevice(1);
  do_test("run_test_irreg_data_empty",run_test_irreg_data_empty);
//  do_test("run_test_irreg_data_has_index_empty",run_test_irreg_data_has_index_empty);
  do_test("run_test_irreg_data_use_constructor",run_test_irreg_data_use_constructor);
  do_test("run_test_irreg_data_set_val",run_test_irreg_data_set_val);
  do_test("run_test_ebforall_init",run_test_ebforall_init);
  do_test("run_test_ebforall_kernel",run_test_ebforall_kernel);
  do_test("run_test_ebforall_kernel_box_no_impact",run_test_ebforall_kernel_box_no_impact);
  do_test("run_test_ebforall_i_init",run_test_ebforall_i_init);
  do_test("run_test_ebforall_i_kernel",run_test_ebforall_i_kernel);
  do_test("run_test_ebforall_i_kernel_box_no_impact",run_test_ebforall_i_kernel_box_no_impact);
//  do_test("run_test_agg_stencil_empty",run_test_agg_stencil_empty);
  do_test("run_test_agg_stencil_kernel_only_using",run_test_agg_stencil_kernel_only_using);
  do_test("run_test_agg_stencil_cpu_versu_gpu",run_test_agg_stencil_cpu_versu_gpu);
  return 0;
}
