#include "EBProto.H"
#include <implem/Proto_IrregData.H>


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
