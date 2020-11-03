#include<iostream>

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

using Proto::Var;
/****/
PROTO_KERNEL_START
unsigned int  kernel_test_forallF(Var<double, 1>  a_in)
{
  a_in(0) = threadIdx.x*2+3;
  return 0;
}
PROTO_KERNEL_END(kernel_test_forallF, kernel_test_forall)

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
    if(a_ptr[i] != 2*i+3) return false;
  }
  return true;
}


bool test_ebforal_check_answer_init(double* a_ptr, unsigned int a_size)
{
  for(int i = 0 ; i < a_size ; i++)
  {
    if(a_ptr[i] != 0) return false;
  }
  return true;
}

bool run_test_ebforall_init()
{
  unsigned int size = 16;
  double* ptr = new double[size];
  std::vector<Proto::EBIndex<Proto::CELL>> index;
  Proto::Box bx(Proto::Point(0,0,0),Proto::Point(size-1,0,0));

  test_ebforall_fill(ptr,index,size);

  // use this constructor to initialize data on the GPU
  Proto::IrregData<Proto::CELL,double,1> fill(bx, ptr, index);

  bool check = fill.size() == size;

  double a = 0;
  cudaEBForAllIrreg(init_test_forall, bx, fill, a);

#ifdef PROTO_CUDA
  double* checkPtr = new double[size];
  double* devicPtr = fill.data();
  protoMemcpy(checkPtr,devicPtr,size*sizeof(double),protoMemcpyDeviceToHost);
#else
  double* checkPtr = fill.data();
#endif  


  bool result = test_ebforal_check_answer_init(checkPtr, size);

  if(!result) test_ebforall_print(checkPtr,size);
  assert(result);

  index.clear();
  free(ptr);  
  free(checkPtr);  
  return result;
}

bool run_test_ebforall_kernel()
{
  unsigned int size = 16;
  double* ptr = new double[size];
  std::vector<Proto::EBIndex<Proto::CELL>> index;
  Proto::Box bx(Proto::Point(0,0,0),Proto::Point(size-1,0,0));

  test_ebforall_fill(ptr,index,size);

  // use this constructor to initialize data on the GPU
  Proto::IrregData<Proto::CELL,double,1> fill(bx, ptr, index);
  cudaEBForAllIrreg(kernel_test_forall, bx, fill);

#ifdef PROTO_CUDA
  double* checkPtr = new double[size];
  double* devicPtr = fill.data();
  protoMemcpy(checkPtr,devicPtr,size*sizeof(double),protoMemcpyDeviceToHost);
#else
  double* checkPtr = fill.data();
#endif  


  bool result = test_ebforal_check_answer(checkPtr, size);
  assert(result);

  index.clear();
  free(ptr);  
  free(checkPtr);  
  return result;
}
