#include<iostream>
#include<src/test_ebforall_stuff.cpp>

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

  std::cout << " memtype =" << MEMTYPE_DEFAULT << std::endl;
  
  protoEBForAllIrreg(init_test_forall, fill, a);

  double* checkPtr = new double[size];
  double* devicPtr = fill.data();
  protoMemcpy(MEMTYPE_DEFAULT, checkPtr,devicPtr,size*sizeof(double),protoMemcpyDeviceToHost);


  bool result = test_ebforal_check_answer_val(checkPtr, double(0), size);

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
  protoEBForAllIrreg(kernel_test_forall, fill);
  double* checkPtr = new double[size];
  double* devicPtr = fill.data();
  protoMemcpy(MEMTYPE_DEFAULT, checkPtr,devicPtr,size*sizeof(double),protoMemcpyDeviceToHost);

  bool result = test_ebforal_check_answer(checkPtr, size);
  if(!result) test_ebforall_print(checkPtr,size);
  assert(result);

  index.clear();
  free(ptr);  
  free(checkPtr);  
  return result;
}

#ifdef Debug_No
bool run_test_ebforall_kernel_box_no_impact()
{
  unsigned int size = 16;
  double* ptr = new double[size];
  std::vector<Proto::EBIndex<Proto::CELL>> index;
  Proto::Box bx(Proto::Point(0,0,0),Proto::Point(size-1,0,0));
  Proto::Box bxminustwo(Proto::Point(2,0,0),Proto::Point(size-3,0,0));

  test_ebforall_fill(ptr,index,size);

  // use this constructor to initialize data on the GPU
  Proto::IrregData<Proto::CELL,double,1> fill_1(bx, ptr, index);
  Proto::IrregData<Proto::CELL,double,1> fill_2(bx, ptr, index);
  cudaEBForAllIrreg(kernel_test_forall, bx, fill_1);
  cudaEBForAllIrreg(kernel_test_forall, bxminustwo, fill_2);

#ifdef PROTO_CUDA
  double* checkPtr_2 = new double[size];
  double* checkPtr_1 = new double[size];
  double* devicPtr_2 = fill_2.data();
  double* devicPtr_1 = fill_1.data();
  protoMemcpy(MEMTYPE_DEFAULT, checkPtr_2,devicPtr_2,size*sizeof(double),protoMemcpyDeviceToHost);
  protoMemcpy(MEMTYPE_DEFAULT, checkPtr_1,devicPtr_1,size*sizeof(double),protoMemcpyDeviceToHost);
#else
  double* checkPtr_2 = fill_2.data();
  double* checkPtr_1 = fill_1.data();
#endif  

  bool result = test_ebforal_check_same_result(checkPtr_1, checkPtr_2, size);
  assert(result);

  index.clear();
  free(ptr);  
  free(checkPtr_1);  
  free(checkPtr_2);  
  return result;
}
#endif

template<typename F, typename... T>
void hackVecIndexer(unsigned int a_begin, unsigned int a_end, F func, Proto::EBIrregStruct<Proto::CELL,double, 1>* a_dst, T... args)
{
  unsigned int size = a_end - a_begin;
  protoLaunchKernelT<MEMTYPE_DEFAULT, vecIndexer<Proto::CELL,double,1, F, T...>>
			(
				1, size, //small test so nb block = 1
				0, size, func, a_dst, args...
			);
}

bool run_test_ebforall_vec_indexer()
{
  unsigned int size = 16;
  double val        = 5;
  double* ptr       = new double[size];

  std::vector<Proto::EBIndex<Proto::CELL>> index;
  Proto::Box bx(Proto::Point(0,0,0),Proto::Point(size-1,0,0));
  test_ebforall_fill(ptr,index,size);

  // use this constructor to initialize data on the GPU
  Proto::IrregData<Proto::CELL,double,1> fill(bx, ptr, index);
  Proto::EBIrregStruct<Proto::CELL,double, 1>* eb_irreg_struct_ptr = fill.getEBIrregDataPtr();

  hackVecIndexer(
			0, size,
			kernel_test_forall_val,
			eb_irreg_struct_ptr,
			val
		);

			 
  double* checkPtr = new double[size];
  double* devicPtr = fill.data();
  protoMemcpy(MEMTYPE_DEFAULT, checkPtr,devicPtr,size*sizeof(double),protoMemcpyDeviceToHost);

  bool result = test_ebforal_check_answer_val(checkPtr, val, size);
  assert(result);

  index.clear();
  free(ptr);  
  free(checkPtr);  
  return result;
}

