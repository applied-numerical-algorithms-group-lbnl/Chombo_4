#pragma once
#include "test_timer.H"
#include "test_IrregData.cu"

void run_perf_irreg_data_copy_args(unsigned int size)
{
  test_timer t_copy;

  double* ptr = new double[size];
  std::vector<Proto::EBIndex<Proto::CELL>> index;
  Proto::Box bx(Proto::Point(0,0,0),Proto::Point(size-1,0,0));

  test_irreg_data_fill(ptr,index,size);
  Proto::IrregData<Proto::CELL,double,1> in(bx, ptr, index);
  Proto::IrregData<Proto::CELL,double,1> out(bx, ptr, index);
  in.setVal(1);
  out.setVal(2);

  t_copy.begin();
  out.copy(in,bx,0,bx,0,1);
  t_copy.end();

  std::cout << " run_test_perf_copy_irreg_data_" << size << " ... " << t_copy.duration() << " ms" << std::endl;

  index.clear();
  free(ptr);
}



void run_perf_irreg_linear_full_args(unsigned int size)
{
  test_timer t_linearIn, t_linearOut;
  double* ptr = new double[size];
  double* ptr2 = new double[size];
  unsigned int inf  = 2;
  unsigned int high = 4;
  double inNumber   = 1;
  double outNumber  = 2;
  std::vector<Proto::EBIndex<Proto::CELL>> index;
  Proto::Box bx(Proto::Point(0,0,0),Proto::Point(size-1,0,0));

  test_irreg_data_fill(ptr,index,size);
  index.clear();
  test_irreg_data_fill(ptr2,index,size);
  Proto::IrregData<Proto::CELL,double,1> in(bx, ptr2, index);
  Proto::IrregData<Proto::CELL,double,1> out(bx, ptr, index);

  void* inWork;
  size_t nBytes = in.charsize(bx,0,1);
  protoMalloc(inWork, nBytes);

  in.setVal(inNumber);
  out.setVal(outNumber);

  /* test copy */
  t_linearIn.begin();
  in.linearOut(inWork,bx,0,0); 
  t_linearIn.end();
  t_linearOut.begin();
  out.linearIn(inWork,bx,0,0); 
  t_linearOut.end();

  std::cout << " run_test_perf_irreg_data_linear_in_" << size << " ... " << t_linearIn.duration() << " ms" << std::endl;
  std::cout << " run_test_perf_irreg_data_linear_out_" << size << " ... " << t_linearOut.duration() << " ms" << std::endl;
  index.clear();
  free(ptr);  
  free(ptr2);  
  protoFree(inWork); 
}
