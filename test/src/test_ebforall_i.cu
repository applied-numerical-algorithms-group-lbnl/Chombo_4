#include<iostream>

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

bool run_test_ebforall_i_init()
{
  unsigned int size = 16;
  double* ptr = new double[size];
  std::vector<Proto::EBIndex<Proto::CELL>> index;
  Proto::Box bx(Proto::Point(0,0,0),Proto::Point(size-1,0,0));

  test_ebforall_i_fill(ptr,index,size);

  // use this constructor to initialize data on the GPU
  Proto::IrregData<Proto::CELL,double,1> fill(bx, ptr, index);

  bool check = fill.size() == size;

  double a = 0;
  cudaEBForAllIrreg_i(init_test_forall_i, bx, fill, a);

#ifdef PROTO_CUDA
  double* checkPtr = new double[size];
  double* devicPtr = fill.data();
  protoMemcpy(checkPtr,devicPtr,size*sizeof(double),protoMemcpyDeviceToHost);
#else
  double* checkPtr = fill.data();
#endif  


  bool result = test_ebforal_i_check_answer_init(checkPtr, size);

  if(!result) test_ebforall_i_print(checkPtr,size);
  assert(result);

  index.clear();
  free(ptr);  
  free(checkPtr);  
  return result;
}

bool run_test_ebforall_i_kernel()
{
  unsigned int size = 16;
  double* ptr = new double[size];
  std::vector<Proto::EBIndex<Proto::CELL>> index;
  Proto::Box bx(Proto::Point(0,0,0),Proto::Point(size-1,0,0));

  test_ebforall_i_fill(ptr,index,size);

  // use this constructor to initialize data on the GPU
  Proto::IrregData<Proto::CELL,double,1> fill(bx, ptr, index);
  cudaEBForAllIrreg_i(kernel_test_forall_i, bx, fill);

#ifdef PROTO_CUDA
  double* checkPtr = new double[size];
  double* devicPtr = fill.data();
  protoMemcpy(checkPtr,devicPtr,size*sizeof(double),protoMemcpyDeviceToHost);
#else
  double* checkPtr = fill.data();
#endif  


  bool result = test_ebforal_i_check_answer(checkPtr, size);
  assert(result);

  index.clear();
  free(ptr);  
  free(checkPtr);  
  return result;
}

bool run_test_ebforall_i_kernel_box_no_impact()
{
  unsigned int size = 16;
  double* ptr = new double[size];
  std::vector<Proto::EBIndex<Proto::CELL>> index;
  Proto::Box bx(Proto::Point(0,0,0),Proto::Point(size-1,0,0));
  Proto::Box bxminustwo(Proto::Point(2,0,0),Proto::Point(size-3,0,0));

  test_ebforall_i_fill(ptr,index,size);

  // use this constructor to initialize data on the GPU
  Proto::IrregData<Proto::CELL,double,1> fill_1(bx, ptr, index);
  Proto::IrregData<Proto::CELL,double,1> fill_2(bx, ptr, index);
  cudaEBForAllIrreg_i(kernel_test_forall_i, bx, fill_1);
  cudaEBForAllIrreg_i(kernel_test_forall_i, bxminustwo, fill_2);

#ifdef PROTO_CUDA
  double* checkPtr_2 = new double[size];
  double* checkPtr_1 = new double[size];
  double* devicPtr_2 = fill_2.data();
  double* devicPtr_1 = fill_1.data();
  protoMemcpy(checkPtr_2,devicPtr_2,size*sizeof(double),protoMemcpyDeviceToHost);
  protoMemcpy(checkPtr_1,devicPtr_1,size*sizeof(double),protoMemcpyDeviceToHost);
#else
  double* checkPtr_2 = fill_2.data();
  double* checkPtr_1 = fill_1.data();
#endif  

  bool result = test_ebforal_i_check_same_result(checkPtr_1, checkPtr_2, size);
  assert(result);

  index.clear();
  free(ptr);  
  free(checkPtr_1);  
  free(checkPtr_2);  
  return result;
}
