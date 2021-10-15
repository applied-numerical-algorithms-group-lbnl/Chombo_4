#include<iostream>
#include<gtest/gtest.h>
#include "src/test_IrregData.cu"
#include "src/test_ebforall.cu"
#include "src/test_ebforall_i.cu"
#include "src/test_agg_stencil.cu"

TEST(testSuite,run_test_ebforall_init) {
    EXPECT_TRUE(run_test_ebforall_init());
}
TEST(testSuite,run_test_ebforall_kernel) {
    EXPECT_TRUE(run_test_ebforall_kernel());
}
TEST(testSuite,run_test_ebforall_vec_indexer) {
    EXPECT_TRUE(run_test_ebforall_vec_indexer());
}
TEST(testSuite,run_test_ebforall_i_init) {
    EXPECT_TRUE(run_test_ebforall_i_init());
}
TEST(testSuite,run_test_ebforall_i_kernel) {
    EXPECT_TRUE(run_test_ebforall_i_kernel());
}
TEST(testSuite,run_test_ebforall_i_vec_indexer) {
    EXPECT_TRUE(run_test_ebforall_i_vec_indexer());
}
TEST(testSuite,run_test_irreg_data_empty) {
    EXPECT_TRUE(run_test_irreg_data_empty());
}
TEST(testSuite,run_test_irreg_data_use_constructor) {
    EXPECT_TRUE(run_test_irreg_data_use_constructor());
}
TEST(testSuite,run_test_irreg_data_set_val) {
    EXPECT_TRUE(run_test_irreg_data_set_val());
}
TEST(testSuite,run_test_irreg_copy) {
    EXPECT_TRUE(run_test_irreg_copy());
}
TEST(testSuite,run_test_irreg_copy_partial) {
    EXPECT_TRUE(run_test_irreg_copy_partial());
}
TEST(testSuite,run_test_irreg_linear_partial) {
    EXPECT_TRUE(run_test_irreg_linear_partial());
}
TEST(testSuite,run_test_irreg_linear_full) {
    EXPECT_TRUE(run_test_irreg_linear_full());
}
TEST(testSuite,run_test_irreg_linear_multiple) {
    EXPECT_TRUE(run_test_irreg_linear_multiple());
}
//  TEST(testSuite,run_test_irreg_linear_multiple);
TEST(testSuite,run_test_irreg_linear_partial_aliasing_define) {
    EXPECT_TRUE(run_test_irreg_linear_partial_aliasing_define());
}
//  TEST(testSuite,run_test_irreg_copy_partial_idst); // test is not correct
TEST(testSuite,run_test_agg_stencil_kernel_only_using) {
    EXPECT_TRUE(run_test_agg_stencil_kernel_only_using());
}
TEST(testSuite,run_test_agg_stencil_increment_only_true) {
    EXPECT_TRUE(run_test_agg_stencil_increment_only_true());
}
TEST(testSuite,run_test_agg_stencil_increment_only_false) {
    EXPECT_TRUE(run_test_agg_stencil_increment_only_false());
}
TEST(testSuite,run_test_agg_stencil_scale_0) {
    EXPECT_TRUE(run_test_agg_stencil_scale_0());
}
TEST(testSuite,run_test_agg_stencil_scale_100) {
    EXPECT_TRUE(run_test_agg_stencil_scale_100());
}
TEST(testSuite,run_test_agg_stencil_scale_minus10) {
    EXPECT_TRUE(run_test_agg_stencil_scale_minus10());
}
TEST(testSuite,run_test_irreg_copy_stress) {
    EXPECT_TRUE(run_test_irreg_copy_stress());
}
TEST(testSuite,run_test_irreg_linear_full_stress) {
    EXPECT_TRUE(run_test_irreg_linear_full_stress());
}
TEST(testSuite,run_test_irreg_linear_partial_stress) {
    EXPECT_TRUE(run_test_irreg_linear_partial_stress());
}
TEST(testSuite,run_test_irreg_linear_partial_aliasing_define_stress) {
    EXPECT_TRUE(run_test_irreg_linear_partial_aliasing_define_stress());
}

int main()
{
  ::testing::InitGoogleTest();
#ifdef PROTO_CUDA
  cudaSetDevice(0);
#endif
  int result = RUN_ALL_TESTS();
  return result;
}
