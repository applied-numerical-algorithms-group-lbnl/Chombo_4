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

  std::cout << " run_test_perf_copy_irreg_data_" << size << " ... " << t_copy.duration() << " ms or we copy" << size*1000/t_copy.duration() << " elem/s" << std::endl;

  index.clear();
  free(ptr);
}



void run_perf_irreg_linear_full_args(unsigned int size)
{
  test_timer t_linearIn, t_linearOut;
  double* ptr = new double[size];
  double* ptr2 = new double[size];
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
  protoMalloc(DEVICE,inWork, nBytes);

  in.setVal(inNumber);
  out.setVal(outNumber);

  /* test copy */
  t_linearOut.begin();
  in.linearOut(inWork,bx,0,0); 
  t_linearOut.end();
  t_linearIn.begin();
  out.linearIn(inWork,bx,0,0); 
  t_linearIn.end();

  std::cout << " run_test_perf_irreg_data_linear_in_" << size << " ... " << t_linearIn.duration() << " ms" << std::endl;
  std::cout << " run_test_perf_irreg_data_linear_out_" << size << " ... " << t_linearOut.duration() << " ms" << std::endl;
  index.clear();
  free(ptr);  
  free(ptr2);  
  protoFree(DEVICE,inWork); 
}


void run_perf_irreg_linear_partial_args(unsigned int sizebox,unsigned int size)
{
  test_timer t_linearIn, t_linearOut;
  double* ptr = new double[size];
  double* ptr2 = new double[size];
  double inNumber   = 1;
  double outNumber  = 2;
  std::vector<Proto::EBIndex<Proto::CELL>> index;
  Proto::Box bx(Proto::Point(0,0,0),Proto::Point(size-1,0,0));

  unsigned int nbBoxes = size/sizebox;


  test_irreg_data_fill(ptr,index,size);
  index.clear();
  test_irreg_data_fill(ptr2,index,size);
  Proto::IrregData<Proto::CELL,double,1> in(bx, ptr2, index);
  Proto::IrregData<Proto::CELL,double,1> out(bx, ptr, index);

  void* inWork;

  unsigned int start = 0;
  unsigned int end = sizebox;
  size_t nBytes = 0;

  for(unsigned int i = 0; i < nbBoxes ; i++)
  {
    Proto::Box tmp(Proto::Point(start,0,0),Proto::Point(end-1,0,0));
    nBytes += in.charsize(tmp,0,1);
    start += sizebox;
    end += sizebox;
  }

  protoMalloc(DEVICE,inWork, nBytes);
  in.setVal(inNumber);
  out.setVal(outNumber);

  /* test copy */

  start = 0;
  end = sizebox;
  char* charWork = (char*)inWork;

  t_linearOut.begin();
  for(unsigned int i = 0; i < nbBoxes ; i++)
  {
    Proto::Box tmp(Proto::Point(start,0,0),Proto::Point(end-1,0,0));
    in.linearOut(charWork,tmp,0,0);
    charWork += in.charsize(tmp,0,1);
    start += sizebox;
    end += sizebox;
  }
  t_linearOut.end();

  start = 0;
  end = sizebox;
  charWork = (char*)inWork;

  t_linearIn.begin();
  for(unsigned int i = 0; i < nbBoxes ; i++)
  {
    Proto::Box tmp(Proto::Point(start,0,0),Proto::Point(end-1,0,0));
    out.linearIn(charWork,tmp,0,0);
    charWork += in.charsize(tmp,0,1);
    start += sizebox;
    end += sizebox;
  }
  t_linearIn.end();

  std::cout << " run_test_perf_irreg_data_linear_in_" << sizebox << "_" << size << " ... " << t_linearIn.duration() << " ms" << std::endl;
  std::cout << " run_test_perf_irreg_data_linear_out_" << sizebox << "_" << size << " ... " << t_linearOut.duration() << " ms" << std::endl;

  index.clear();
  free(ptr);
  free(ptr2);
  protoFree(DEVICE,inWork);
}


void run_perf_irreg_linear_partial_graph_args(unsigned int sizebox,unsigned int size)
{
  test_timer t_linearIn, t_linearOut;
  double* ptr = new double[size];
  double* ptr2 = new double[size];
  double inNumber   = 1;
  double outNumber  = 2;
  std::vector<Proto::EBIndex<Proto::CELL>> index;
  Proto::Box bx(Proto::Point(0,0,0),Proto::Point(size-1,0,0));

  unsigned int nbBoxes = size/sizebox;


  test_irreg_data_fill(ptr,index,size);
  index.clear();
  test_irreg_data_fill(ptr2,index,size);
  Proto::IrregData<Proto::CELL,double,1> in(bx, ptr2, index);
  Proto::IrregData<Proto::CELL,double,1> out(bx, ptr, index);

  void* inWork;

  unsigned int start = 0;
  unsigned int end = sizebox;
  size_t nBytes = 0;

  for(unsigned int i = 0; i < nbBoxes ; i++)
  {
    Proto::Box tmp(Proto::Point(start,0,0),Proto::Point(end-1,0,0));
    nBytes += in.charsize(tmp,0,1);
    start += sizebox;
    end += sizebox;
  }

  protoMalloc(DEVICE,inWork, nBytes);
  in.setVal(inNumber);
  out.setVal(outNumber);

  /* test copy */

  start = 0;
  end = sizebox;
  char* charWork = (char*)inWork;

  for(unsigned int i = 0; i < nbBoxes ; i++)
  {
    in.linearOut(charWork,bx,0,0);
  }

  start = 0;
  end = sizebox;
  charWork = (char*)inWork;

  cudaStream_t stream;
  cudaStreamCreate(&stream);
  cudaGraph_t graph;
  cudaGraphExec_t instance;

  cudaStreamBeginCapture(stream, cudaStreamCaptureModeGlobal);
  //{
    Proto::Box tmp(Proto::Point(start,0,0),Proto::Point(end-1,0,0));
  for(unsigned int i = 0; i < nbBoxes ; i++)
    out.linearIn(charWork,tmp,0,0);
  //  charWork += in.charsize(tmp,0,1);
    start += sizebox;
    end += sizebox;
  //}
  cudaStreamEndCapture(stream, &graph);
  cudaGraphInstantiate(&instance, graph, NULL, NULL, 0);

  cudaStreamSynchronize(stream);
  t_linearIn.begin();
  cudaGraphLaunch(instance, stream);
  cudaStreamSynchronize(stream);
  t_linearIn.end();

  std::cout << " run_test_perf_irreg_data_linear_in_" << sizebox << "_" << size << " ... " << t_linearIn.duration() << " ms" << std::endl;

  index.clear();
  free(ptr);
  free(ptr2);
  protoFree(DEVICE,inWork);
}
