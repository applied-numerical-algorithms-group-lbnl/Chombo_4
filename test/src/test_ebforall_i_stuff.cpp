#pragma once
void test_ebforall_i_print(double* ptr, unsigned int n)
{
  for(int i = 0 ; i < n ; i++)
    std::cout << ptr[i] << " " << std::endl;
}

void test_ebforall_i_fill(double* a_ptr, std::vector<Proto::EBIndex<Proto::CELL>>& a_index, unsigned int a_size)
{
  for(int i = 0 ; i < a_size ; i++)
  {
    a_ptr[i] = 1;
    Proto::Point p(i,0,0);
    Proto::EBIndex<Proto::CELL> e(p,i);
    a_index.push_back(e);
  }
}

using Proto::Var;
/****/
PROTO_KERNEL_START
unsigned int  kernel_test_forall_iF(int              a_pt[DIM],
					   Var<double, 1>  a_in)
{
  a_in(0) = a_pt[0] * 10;
  return 0;
}
PROTO_KERNEL_END(kernel_test_forall_iF, kernel_test_forall_i)

PROTO_KERNEL_START
unsigned int  init_test_forall_iF(int              a_pt[DIM],
					Var<double, 1>  a_in, 
						double a_val)
{
  a_in(0) = a_val;
  return 0;
}
PROTO_KERNEL_END(init_test_forall_iF, init_test_forall_i)

bool test_ebforal_i_check_answer(double* a_ptr, unsigned int a_size)
{
  for(int i = 0 ; i < a_size ; i++)
  {
    if(a_ptr[i] != i *10) return false;
  }
  return true;
}

bool test_ebforal_i_check_answer_init(double* a_ptr, unsigned int a_size)
{
  for(int i = 0 ; i < a_size ; i++)
  {
    if(a_ptr[i] != 0) return false;
  }
  return true;
}

bool test_ebforal_i_check_same_result(double* a_ptr1, double* a_ptr2, unsigned int a_size)
{
  for(int i = 0 ; i < a_size ; i++)
    if(a_ptr1[i] != a_ptr2[i]) return false;

  return true;  
}


