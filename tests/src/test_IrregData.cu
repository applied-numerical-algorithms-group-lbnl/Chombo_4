#include "EBProto.H"
#include <implem/Proto_IrregData.H>
#include "test_timer.H"

#define MAX_ELEM 100000

void test_irreg_data_fill(double* ptr, std::vector<Proto::EBIndex<Proto::CELL>>& index, unsigned int size)
{
  index.resize(size);
  for(int i = 0 ; i < size ; i++)
  {
    ptr[i] = i*2;
    Proto::Point p(i,0,0);
    Proto::EBIndex<Proto::CELL> e(p,i);
    index[i]=e;
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
  protoMemcpy(MEMTYPE_DEFAULT, checkPtr,devicPtr,size*sizeof(double),protoMemcpyDeviceToHost);
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
  protoMemcpy(MEMTYPE_DEFAULT, checkPtr,devicPtr,size*sizeof(double),protoMemcpyDeviceToHost);
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
  protoMemcpy(MEMTYPE_DEFAULT, checkPtr,devicPtr,size*sizeof(double),protoMemcpyDeviceToHost);

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
  protoMemcpy(MEMTYPE_DEFAULT, checkPtr,devicPtr,size*sizeof(double),protoMemcpyDeviceToHost);

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
  protoMemcpy(MEMTYPE_DEFAULT, checkPtr,devicPtr,size*sizeof(double),protoMemcpyDeviceToHost);

  bool check=true;
  for(int i = 0 ; i < size ; i++)
    if(i>=inf + shift && i <= high + shift) {
      if(checkPtr[i] != inNumber) check=false;
    } else {
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
  protoMalloc(MEMTYPE_DEFAULT, inWork, nBytes);

  in.setVal(inNumber);
  out.setVal(outNumber);

  /* test copy */
  in.linearOut(inWork,bx,0,0); 
  out.linearIn(inWork,bx,0,0); 

  double* checkPtr = new double[size];
  double* devicPtr = out.data();
  protoMemcpy(MEMTYPE_DEFAULT, checkPtr,devicPtr,size*sizeof(double),protoMemcpyDeviceToHost);

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
  protoFree(MEMTYPE_DEFAULT, inWork); 

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


bool run_test_irreg_linear_multiple()
{
// a lot of realocation in this test
  unsigned int size = 7;

  unsigned int* sizes = new unsigned int[size];
  unsigned int base = 3;
  for(unsigned int i = 0; i < size ; i++)
  {
	sizes[i] = base;
	base *= 3;
  }

  double inNumber   = 1;
  double outNumber  = 2;
  size_t nBytes = 0;

  for(int id = 0 ; id < size ; id++)
  {
    std::vector<Proto::EBIndex<Proto::CELL>> index;
    Proto::Box bx(Proto::Point(0,0,0),Proto::Point(sizes[id]-1,0,0));
    double* ptr = new double[sizes[id]];

    test_irreg_data_fill(ptr,index,sizes[id]);
    Proto::IrregData<Proto::CELL,double,1> in(bx, ptr, index);
    nBytes += in.charsize(bx,0,1); 
    index.clear();
    free(ptr);  
  }

  void* inWork;

  protoMalloc(MEMTYPE_DEFAULT, inWork, nBytes);

  unsigned int shift = 0;


  for(int id = 0 ; id < size ; id++)
  {
    bool check=true;
    std::vector<Proto::EBIndex<Proto::CELL>> index;
    Proto::Box bx(Proto::Point(0,0,0),Proto::Point(sizes[id]-1,0,0));
    double* ptr = new double[sizes[id]];
    double* ptr2 = new double[sizes[id]];

    test_irreg_data_fill(ptr,index,sizes[id]);
    Proto::IrregData<Proto::CELL,double,1> in(bx, ptr, index);
    Proto::IrregData<Proto::CELL,double,1> out(bx, ptr, index);

    in.setVal(inNumber + id);
    out.setVal(outNumber);

    char* work = (char*) (inWork);
    work += shift;
    void *buf = (void*) work;
    /* test copy */
    in.linearOut(buf,bx,0,0); 
    out.linearIn(buf,bx,0,0); 

    double* checkPtr = new double[sizes[id]];
    double* devicPtr = out.data();
    protoMemcpy(MEMTYPE_DEFAULT, checkPtr,devicPtr,sizes[id]*sizeof(double),protoMemcpyDeviceToHost);

    for(int i = 0 ; i < sizes[id] ; i++)
    {
      if(checkPtr[i] != inNumber + id) 
      {
             std::cout << " error for size equal to " << id << " idx = " << i << " res = " << checkPtr[i] << " iNumber = " << inNumber << " and id equal to " << sizes[id] << std::endl;
	     check=false;
      }
    }

    if(!check)
    {
      for(int i = 0 ; i < sizes[id] ; i++)
        std::cout << checkPtr[i] << " ";
 
      std::cout << std::endl;
      free(checkPtr); 
      index.clear();
      free(ptr);  
      free(ptr2);  
      protoFree(MEMTYPE_DEFAULT, inWork); 
      return false;
    }
    free(checkPtr); 
    
    index.clear();
    free(ptr);  
    free(ptr2);  
    shift += in.charsize(bx,0,1); 
  }
  protoFree(MEMTYPE_DEFAULT, inWork); 


  return true;
}


bool run_test_irreg_linear_partial_args(unsigned int size)
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

  void* inWork;
  size_t nBytes = in.charsize(bx2,0,1);//size*(sizeof(double)+sizeof(EBIndex<Proto::CELL>)) + 2*sizeof(unsigned int);
  protoMalloc(MEMTYPE_DEFAULT, inWork, nBytes);

  in.setVal(inNumber);
  out.setVal(outNumber);

  /* test copy */
  in.linearOut(inWork,bx,0,0); 
  out.linearIn(inWork,bx2,0,0); 

  double* checkPtr = new double[size];
  double* devicPtr = out.data();
  protoMemcpy(MEMTYPE_DEFAULT, checkPtr,devicPtr,size*sizeof(double),protoMemcpyDeviceToHost);

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
  protoFree(MEMTYPE_DEFAULT, inWork);
  return check;
}

bool run_test_irreg_linear_partial()
{
  unsigned int size = 8;
  return run_test_irreg_linear_partial_args(size);
}
	  
bool run_test_irreg_linear_partial_stress()
{
        test_timer timer;
        unsigned int size = 64;
        bool b =true;
        while(size < MAX_ELEM && b)
        {
                timer.begin();
                b = run_test_irreg_linear_partial_args(size);
                timer.end();
                if(b) std::cout << " run_test_irreg_linear_partial_stress_" << size << " ok ... " << timer.duration() << " ms" << std::endl;
                else  std::cout << " error for size = " << size << std::endl;
                size *= 2;
        }

  return b;
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
  protoMalloc(MEMTYPE_DEFAULT, inWork, nBytes);

  in.setVal(inNumber);
  out.setVal(outNumber);

  /* test copy */
  inDefined.linearOut(inWork,bx,0,0); 
  outDefined.linearIn(inWork,bx2,0,0); 

  double* checkPtr = new double[size];
  double* devicPtr = out.data();
  protoMemcpy(MEMTYPE_DEFAULT, checkPtr,devicPtr,size*sizeof(double),protoMemcpyDeviceToHost);

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
  protoFree(MEMTYPE_DEFAULT, inWork);

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

