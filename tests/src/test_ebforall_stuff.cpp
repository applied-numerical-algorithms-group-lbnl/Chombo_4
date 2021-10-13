#pragma once

void test_ebforall_print(double* ptr, unsigned int n)
{
  for(int i = 0 ; i < n ; i++)
    std::cout << ptr[i] << " " << std::endl;
}

void test_ebforall_fill(double* a_ptr, std::vector<Proto::EBIndex<Proto::CELL>>& a_index, unsigned int a_size)
{
  for(int i = 0 ; i < a_size ; i++)
  {
    a_ptr[i] = 1;
    Proto::Point p(i,0,0);
    Proto::EBIndex<Proto::CELL> e(p,i);
    a_index.push_back(e);
  }
}

using Proto::Var; using Proto::MEMTYPE_DEFAULT; using Proto::pairPtr; 
using Proto::vecIndexer; using Proto::vecIndexer_i;
/****/
PROTO_KERNEL_START
unsigned int  kernel_test_forallF(Var<double, 1>  a_in)
{
#ifdef PROTO_CUDA  
  a_in(0) = threadIdx.x*2+3;
#else
  a_in(0) = 674;
#endif  
  return 0;
}
PROTO_KERNEL_END(kernel_test_forallF, kernel_test_forall)


PROTO_KERNEL_START
unsigned int  kernel_test_forall_valF(Var<double, 1>  a_in, double a_val)
{
  a_in(0) = a_val;
  return 0;
}
PROTO_KERNEL_END(kernel_test_forall_valF, kernel_test_forall_val)

inline unsigned int  kernel_test_forall_val_hostF(Var<double, 1>  a_in, double a_val)
{
  a_in(0) = a_val;
  return 0;
}
constexpr decltype(&kernel_test_forall_val_hostF) kernel_test_forall_val_host = kernel_test_forall_val_hostF;

PROTO_KERNEL_START
unsigned int  init_test_forallF(Var<double, 1>  a_in, double a_val)
{
  a_in(0) = a_val;
  return 0;
}
PROTO_KERNEL_END(init_test_forallF, init_test_forall)

bool test_ebforal_check_answer(double* a_ptr, unsigned int a_size)
{
  for(int i = 0 ; i < a_size ; i++)
  {
#ifdef PROTO_CUDA
    if(a_ptr[i] != 2*i+3) return false;
#else
    if(a_ptr[i] != 674) return false;
#endif
  }
  return true;
}

bool test_ebforal_check_answer_val(double* a_ptr, double a_val, unsigned int a_size)
{
  for(int i = 0 ; i < a_size ; i++)
  {
    if(a_ptr[i] != a_val) return false;
  }
  return true;
}


bool test_ebforal_check_same_result(double* a_ptr1, double* a_ptr2, unsigned int a_size)
{
  for(int i = 0 ; i < a_size ; i++)
    if(a_ptr1[i] != a_ptr2[i]) return false;

  return true;  
}
