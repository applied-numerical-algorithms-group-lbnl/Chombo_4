#include "EBProto.H"
#include <implem/Proto_IrregData.H>
#include "test_timer.H"

#define MAX_ELEM 500000

void test_irreg_data_fill(double* ptr, std::vector<Proto::EBIndex<Proto::CELL>>& index, unsigned int size)
{
  for(int i = 0 ; i < size ; i++)
  {
    ptr[i] = i*2;
    Proto::Point p(i,0,0);
    Proto::EBIndex<Proto::CELL> e(p,i);
    index.push_back(e);
  }
}


bool test_irreg_data_check_fill(double* ptr, std::vector<Proto::EBIndex<Proto::CELL>>& index, unsigned int size)
{
  for(int i = 0 ; i < size ; i++)
  {
    if(ptr[i] != i*2) return false;
    Proto::Point p(i,0,0);
    Proto::EBIndex<Proto::CELL> e(p,i);
    if(index[i] !=e) return false;
  }
  return true;
}

bool test_irreg_data_check_set_zero(double* ptr, unsigned int size)
{
  for(int i = 0 ; i < size ; i++)
    if(ptr[i] != 0) return false;

  return true;
}

bool run_test_irreg_data_empty()
{
  Proto::IrregData<Proto::CELL,double,1> empty;

  bool check1 = !(empty.defined());
  bool check2 = empty.vecsize() == 0;
  bool check3 = empty.size() == 0;
  //bool check3 = !(empty.hasIndex(0));
  assert(check1);
  assert(check2);
  assert(check3);
  return check1 && check2 && check3;
}

bool run_test_irreg_data_has_index_empty()
{
  Proto::IrregData<Proto::CELL,double,1> empty;
  Proto::EBIndex<Proto::CELL> index;
  bool check = !(empty.hasIndex(index));
  assert(check);
  return check;
}

bool run_test_irreg_data_use_constructor()
{
  unsigned int size = 8;
  double* ptr = new double[size];
  std::vector<Proto::EBIndex<Proto::CELL>> index;
  Proto::Box bx(Proto::Point(0,0,0),Proto::Point(size-1,0,0));

  test_irreg_data_fill(ptr,index,size);

  bool sizeOK = size == index.size();
  if(!sizeOK) std::cout << " size: " << size << "!= index.size(): " << index.size() << std::endl;
  assert(size == index.size());

  bool before = test_irreg_data_check_fill(ptr,index,size);
  assert(before);

  // use this constructor to initialize data on the GPU
  Proto::IrregData<Proto::CELL,double,1> fill(bx, ptr, index);

  bool check = fill.size() == size;

#ifdef PROTO_CUDA
  double* checkPtr = new double[size];
  double* devicPtr = fill.data();
  protoMemcpy(checkPtr,devicPtr,size*sizeof(double),protoMemcpyDeviceToHost);
#else
  double* checkPtr = fill.data();
#endif  

  bool after = test_irreg_data_check_fill(checkPtr, *(fill.getIndicies()), size);
  assert(after);

  index.clear();
  free(ptr);  

  return check && before && after;
}

bool run_test_irreg_data_set_val()
{
  unsigned int size = 8;
  double* ptr = new double[size];
  std::vector<Proto::EBIndex<Proto::CELL>> index;
  Proto::Box bx(Proto::Point(0,0,0),Proto::Point(size-1,0,0));

  test_irreg_data_fill(ptr,index,size);
  Proto::IrregData<Proto::CELL,double,1> fill(bx, ptr, index);
  fill.setVal(0);

#ifdef PROTO_CUDA
  double* checkPtr = new double[size];
  double* devicPtr = fill.data();
  protoMemcpy(checkPtr,devicPtr,size*sizeof(double),protoMemcpyDeviceToHost);
#else
  double* checkPtr = fill.data();
#endif  

  // should be false
  bool nochange = test_irreg_data_check_fill(checkPtr, *(fill.getIndicies()), size);
  assert(!nochange);

  bool change = test_irreg_data_check_set_zero(checkPtr,size);

  index.clear();
  free(ptr);  

  return change && (!nochange);
}


bool run_test_irreg_copy_args(unsigned int size)
{
  double* ptr = new double[size];
  std::vector<Proto::EBIndex<Proto::CELL>> index;
  Proto::Box bx(Proto::Point(0,0,0),Proto::Point(size-1,0,0));

  test_irreg_data_fill(ptr,index,size);
  Proto::IrregData<Proto::CELL,double,1> in(bx, ptr, index);
  Proto::IrregData<Proto::CELL,double,1> out(bx, ptr, index);
  in.setVal(1);
  out.setVal(2);

  out.copy(in,bx,0,bx,0,1); 

  double* checkPtr = new double[size];
  double* devicPtr = out.data();
  protoMemcpy(checkPtr,devicPtr,size*sizeof(double),protoMemcpyDeviceToHost);

  bool check=true;
  for(int i = 0 ; i < size ; i++)
    if(checkPtr[i] != 1) check=false;

  if(!check)
  {
    for(int i = 0 ; i < size ; i++)
      std::cout << checkPtr[i] << " ";
 
    std::cout << std::endl;
  }
  

  index.clear();
  free(ptr);  
  free(checkPtr);
 
  return check;
}

bool run_test_irreg_copy()
{
  return run_test_irreg_copy_args(8);
}


bool run_test_irreg_copy_stress()
{
	test_timer timer;
	unsigned int size = 64;
	bool b =true;
	while(size < MAX_ELEM && b)
	{
		timer.begin();
		b = run_test_irreg_copy_args(size);
		timer.end();
	 	if(b) std::cout << " run_test_irreg_copy_stress_" << size << " ok ... " << timer.duration() << " ms" << std::endl;
		else  std::cout << " error for size = " << size << std::endl;
		size *= 2;
	}

  return b;
}

bool run_test_irreg_copy_partial()
{
  unsigned int size = 8;
  double* ptr = new double[size];
  unsigned int inf  = 2;
  unsigned int high = 5;
  double inNumber   = 1;
  double outNumber  = 2;
  std::vector<Proto::EBIndex<Proto::CELL>> index;
  Proto::Box bx(Proto::Point(0,0,0),Proto::Point(size-1,0,0));
  Proto::Box bx2(Proto::Point(inf,0,0),Proto::Point(high,0,0));

  test_irreg_data_fill(ptr,index,size);
  Proto::IrregData<Proto::CELL,double,1> in(bx, ptr, index);
  Proto::IrregData<Proto::CELL,double,1> out(bx, ptr, index);
  in.setVal(inNumber);
  out.setVal(outNumber);

  /* test copy */
  out.copy(in,bx2,0,bx2,0,1); 

  double* checkPtr = new double[size];
  double* devicPtr = out.data();
  protoMemcpy(checkPtr,devicPtr,size*sizeof(double),protoMemcpyDeviceToHost);

  bool check=true;
  for(int i = 0 ; i < size ; i++)
    if(i>=inf && i <= high)
    {
      if(checkPtr[i] != inNumber) check=false;
    }
    else
    {
      if(checkPtr[i] != outNumber) check=false;
    }

  if(!check)
  {
    for(int i = 0 ; i < size ; i++)
      std::cout << checkPtr[i] << " ";
 
    std::cout << std::endl;
  }
  index.clear();
  free(ptr);  
  free(checkPtr);

  return check;
}



bool run_test_irreg_copy_partial_idst()
{
  unsigned int size = 8;
  double* ptr = new double[size];
  unsigned int inf  = 2;
  unsigned int high = 4;
  double inNumber   = 1;
  double outNumber  = 2;
  unsigned int shift = 2;
  std::vector<Proto::EBIndex<Proto::CELL>> index;
  Proto::Box bx(Proto::Point(0,0,0),Proto::Point(size-1,0,0));
  Proto::Box bx2(Proto::Point(inf,0,0),Proto::Point(high,0,0));

  test_irreg_data_fill(ptr,index,size);
  Proto::IrregData<Proto::CELL,double,1> in(bx, ptr, index);
  Proto::IrregData<Proto::CELL,double,1> out(bx, ptr, index);
  in.setVal(inNumber);
  out.setVal(outNumber);

  /* test copy */
  out.copy(in,bx2,0,bx2,shift,1); // 1 doesn't work in this case --> IrregData<Proto::CELL,double,2>  

  double* checkPtr = new double[size];
  double* devicPtr = out.data();
  protoMemcpy(checkPtr,devicPtr,size*sizeof(double),protoMemcpyDeviceToHost);

  bool check=true;
  for(int i = 0 ; i < size ; i++)
    if(i>=inf + shift && i <= high + shift)
      if(checkPtr[i] != inNumber) check=false;
    else
      if(checkPtr[i] != outNumber) check=false;

  if(!check)
  {
    for(int i = 0 ; i < size ; i++)
      std::cout << checkPtr[i] << " ";
 
    std::cout << std::endl;
  }
  index.clear();
  free(ptr);  
  free(checkPtr);  

  return check;
}


bool run_test_irreg_linear_full_args(unsigned int size)
{
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
  size_t nBytes = in.charsize(bx,0,1);//size*(sizeof(double)+sizeof(EBIndex<Proto::CELL>)) + 2*sizeof(unsigned int);
  protoMalloc(inWork, nBytes);

  in.setVal(inNumber);
  out.setVal(outNumber);

  /* test copy */
  in.linearOut(inWork,bx,0,0); 
  out.linearIn(inWork,bx,0,0); 

  double* checkPtr = new double[size];
  double* devicPtr = out.data();
  protoMemcpy(checkPtr,devicPtr,size*sizeof(double),protoMemcpyDeviceToHost);

  bool check=true;
  for(int i = 0 ; i < size ; i++)
    {
      if(checkPtr[i] != inNumber) check=false;
    }

  if(!check)
  {
    for(int i = 0 ; i < size ; i++)
      std::cout << checkPtr[i] << " ";
 
    std::cout << std::endl;
  }
  index.clear();
  free(ptr);  
  free(ptr2);  
  free(checkPtr); 
  protoFree(inWork); 

  return check;
}

bool run_test_irreg_linear_full()
{
  unsigned int size = 8;
  return run_test_irreg_linear_full_args(size);
}

bool run_test_irreg_linear_full_stress()
{
        test_timer timer;
        unsigned int size = 64;
        bool b =true;
        while(size < MAX_ELEM && b)
        {
                timer.begin();
                b = run_test_irreg_linear_full_args(size);
                timer.end();
                if(b) std::cout << " run_test_irreg_linear_full_stress_" << size << " ok ... " << timer.duration() << " ms" << std::endl;
                else  std::cout << " error for size = " << size << std::endl;
                size *= 2;
        }

  return b;
}

bool run_test_irreg_linear_partial()
{
  unsigned int size = 8;
  double* ptr = new double[size];
  double* ptr2 = new double[size];
  unsigned int inf  = 2;
  unsigned int high = 4;
  double inNumber   = 1;
  double outNumber  = 2;
  std::vector<Proto::EBIndex<Proto::CELL>> index;
  Proto::Box bx(Proto::Point(0,0,0),Proto::Point(size-1,0,0));
  Proto::Box bx2(Proto::Point(inf,0,0),Proto::Point(high,0,0));

  test_irreg_data_fill(ptr,index,size);
  index.clear();
  test_irreg_data_fill(ptr2,index,size);
  Proto::IrregData<Proto::CELL,double,1> in(bx, ptr2, index);
  Proto::IrregData<Proto::CELL,double,1> out(bx, ptr, index);

  void* inWork;
  size_t nBytes = in.charsize(bx2,0,1);//size*(sizeof(double)+sizeof(EBIndex<Proto::CELL>)) + 2*sizeof(unsigned int);
  protoMalloc(inWork, nBytes);

  in.setVal(inNumber);
  out.setVal(outNumber);

  /* test copy */
  in.linearOut(inWork,bx,0,0); 
  out.linearIn(inWork,bx2,0,0); 

  double* checkPtr = new double[size];
  double* devicPtr = out.data();
  protoMemcpy(checkPtr,devicPtr,size*sizeof(double),protoMemcpyDeviceToHost);

  bool check=true;
  for(int i = 0 ; i < size ; i++)
    if(i>=inf && i <= high )
    {
      if(checkPtr[i] != inNumber) check=false;
    }
    else
    {
      if(checkPtr[i] != outNumber) check=false;
    }

  if(!check)
  {
    for(int i = 0 ; i < size ; i++)
      std::cout << checkPtr[i] << " ";
 
    std::cout << std::endl;
  }
  index.clear();
  free(ptr);  
  free(ptr2);  
  free(checkPtr);  
  protoFree(inWork);
  return check;
}


bool run_test_irreg_linear_partial_aliasing_define_args(unsigned int size)
{
  double* ptr = new double[size];
  double* ptr2 = new double[size];
  unsigned int inf  = 2;
  unsigned int high = 4;
  double inNumber   = 1;
  double outNumber  = 2;
  std::vector<Proto::EBIndex<Proto::CELL>> index;
  Proto::Box bx(Proto::Point(0,0,0),Proto::Point(size-1,0,0));
  Proto::Box bx2(Proto::Point(inf,0,0),Proto::Point(high,0,0));

  test_irreg_data_fill(ptr,index,size);
  index.clear();
  test_irreg_data_fill(ptr2,index,size);
  Proto::IrregData<Proto::CELL,double,1> in(bx, ptr2, index);
  Proto::IrregData<Proto::CELL,double,1> out(bx, ptr, index);

  Proto::IrregData<Proto::CELL,double,1> inDefined;
  Proto::IrregData<Proto::CELL,double,1> outDefined;

  unsigned int comp = 0;
  inDefined.define(in,comp);  
  outDefined.define(out,comp);  

  void* inWork;
  size_t nBytes = in.charsize(bx2,0,1);//size*(sizeof(double)+sizeof(EBIndex<Proto::CELL>)) + 2*sizeof(unsigned int);
  protoMalloc(inWork, nBytes);

  in.setVal(inNumber);
  out.setVal(outNumber);

  /* test copy */
  inDefined.linearOut(inWork,bx,0,0); 
  outDefined.linearIn(inWork,bx2,0,0); 

  double* checkPtr = new double[size];
  double* devicPtr = out.data();
  protoMemcpy(checkPtr,devicPtr,size*sizeof(double),protoMemcpyDeviceToHost);

  bool check=true;
  for(int i = 0 ; i < size ; i++)
    if(i>=inf && i <= high )
    {
      if(checkPtr[i] != inNumber) check=false;
    }
    else
    {
      if(checkPtr[i] != outNumber) check=false;
    }

  if(!check)
  {
    for(int i = 0 ; i < size ; i++)
      std::cout << checkPtr[i] << " ";
 
    std::cout << std::endl;
  }
  index.clear();
  free(ptr);  
  free(ptr2);  
  free(checkPtr);  
  protoFree(inWork);

  return check;
}


bool run_test_irreg_linear_partial_aliasing_define()
{
  unsigned int size = 8;
  return run_test_irreg_linear_partial_aliasing_define_args(size);
}

bool run_test_irreg_linear_partial_aliasing_define_stress()
{
        test_timer timer;
        unsigned int size = 64;
        bool b =true;
        while(size < MAX_ELEM && b)
        {
                timer.begin();
                b = run_test_irreg_linear_partial_aliasing_define_args(size);
                timer.end();
                if(b) std::cout << " run_test_irreg_linear_partial_aliasing_define_stress_" << size << " ok ... " << timer.duration() << " ms" << std::endl;
                else  std::cout << " error for size = " << size << std::endl;
                size *= 2;
        }

  return b;
}

