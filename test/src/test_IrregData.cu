#include "EBProto.H"
#include <implem/Proto_IrregData.H>

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
//  index.resize(size);
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
