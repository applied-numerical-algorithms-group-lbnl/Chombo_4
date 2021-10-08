
#include "Proto.H"
using namespace Proto;

///  after this are specific to the test
typedef Var<double,DIM> Vec;

PROTO_KERNEL_START 
void UsetUF(Vec a_U, double  a_val)
{
//  printf("in set U\n");
}
PROTO_KERNEL_END(UsetUF, UsetU)



PROTO_KERNEL_START 
void VsetVF(Vec a_V, double  a_val)
{
//  printf("in set V\n");
}
PROTO_KERNEL_END(VsetVF, VsetV)




int main(int argc, char* argv[])
{

  int nx = 4;
  Box domain(Point::Zeros(), Point::Ones(nx-1));
  BoxData<double, DIM> U(domain);
  BoxData<double, DIM> V(domain);
  BoxData<double, DIM> W(domain);
  double uval = 1;
  double vval = 2;
//  printf("going into setU\n");
  Box grid = domain;
  forallInPlace(UsetU, grid, U, uval);
#ifdef PROTO_CUDA
  protoDeviceSynchronize();
#endif
//  printf("going into setV\n");
  forallInPlace(VsetV, grid, V, vval);
#ifdef PROTO_CUDA
  protoDeviceSynchronize();
#endif
  printf("func_mess PASSED\n");
}








