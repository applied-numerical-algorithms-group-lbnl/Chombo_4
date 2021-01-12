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
#ifdef PROTO_CUDA
  cudaSetDevice(1);
#endif
  do_test("run_test_irreg_data_empty",run_test_irreg_data_empty);
  do_test("run_test_irreg_data_use_constructor",run_test_irreg_data_use_constructor);
  do_test("run_test_irreg_data_set_val",run_test_irreg_data_set_val);
  do_test("run_test_irreg_copy",run_test_irreg_copy);
  do_test("run_test_irreg_copy_partial",run_test_irreg_copy_partial);
  do_test("run_test_irreg_linear_in_partial",run_test_irreg_linear_in_partial);

//  do_test("run_test_irreg_copy_partial_idst",run_test_irreg_copy_partial_idst); // test is not correct
  do_test("run_test_ebforall_init",run_test_ebforall_init);
  do_test("run_test_ebforall_kernel",run_test_ebforall_kernel);
  do_test("run_test_ebforall_vec_indexer",run_test_ebforall_vec_indexer);
  do_test("run_test_ebforall_i_init",run_test_ebforall_i_init);
  do_test("run_test_ebforall_i_kernel",run_test_ebforall_i_kernel);
  do_test("run_test_ebforall_i_vec_indexer",run_test_ebforall_i_vec_indexer);
  do_test("run_test_agg_stencil_kernel_only_using",run_test_agg_stencil_kernel_only_using);
  do_test("run_test_agg_stencil_increment_only_true",run_test_agg_stencil_increment_only_true);
  do_test("run_test_agg_stencil_increment_only_false",run_test_agg_stencil_increment_only_false);
  do_test("run_test_agg_stencil_scale_0",run_test_agg_stencil_scale_0);
  do_test("run_test_agg_stencil_scale_100",run_test_agg_stencil_scale_100);
  do_test("run_test_agg_stencil_scale_minus10",run_test_agg_stencil_scale_minus10);
  return 0;
}
