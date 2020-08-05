#include "Proto.H"
#include "EBProto.H"

using namespace Proto;





///  after this are specific to the test
typedef Var<double,DIM> V;


PROTO_KERNEL_START 
void UsetUF(V a_U, double  a_val)
{
//  printf("in set U\n");
//  printf("setu: uptr[0] = %p, uptr[1] = %p\n",a_U.m_ptrs[0],a_U.m_ptrs[1]);
  for(int idir = 0; idir < DIM; idir++)
  {
    a_U(idir) = a_val;
    if(a_U(idir) != a_val)
    {
      printf("p1: values do not match \n");
      printf("setu: val = %f, uval = %f\n",a_val, a_U(idir));
    }
  }
}
PROTO_KERNEL_END(UsetUF, UsetU)


PROTO_KERNEL_START 
void setUptF(int  a_p[DIM], V a_U, double  a_val)
{
  for(int idir = 0; idir < DIM; idir++)
  {
    a_U(idir) = a_val;
    if(a_U(idir) != a_val)
    {
      printf("upt: values do not match \n");
      printf("upt: val = %f, uval = %f\n",a_val, a_U(idir));
    }
  }
}
PROTO_KERNEL_END(setUptF, setUpt)


PROTO_KERNEL_START 
void VsetVF(V a_V, double  a_val, int a_intvar)
{
//  printf("setv: vptr[0] = %p, vptr[1] = %p\n",a_V.m_ptrs[0],a_V.m_ptrs[1]);
//  printf("in set V\n");
 for(int idir = 0; idir < DIM; idir++)
 {
   a_V(idir) = a_val;
   if(a_V(idir) != a_val)
   {
     printf("setv: values do not match \n");
     //     printf("setv: val = %f, vval = %f\n",a_val, a_V(idir));
   }
 }
}
PROTO_KERNEL_END(VsetVF, VsetV)

PROTO_KERNEL_START 
void setVptF(int  a_p[DIM], V a_V, double  a_val, int a_vvar)
{
  for(int idir = 0; idir < DIM; idir++)
  {
   a_V(idir) = a_val;
   if(a_V(idir) != a_val)
   {
     printf("vpt: values do not match \n");
//     printf("setv: val = %f, vval = %f\n",a_val, a_V(idir));
   }
  }
}
PROTO_KERNEL_END(setVptF, setVpt)


PROTO_KERNEL_START 
void WsetWtoUplusVF(V a_W,
                    V a_U,
                    V a_V,
                    double  a_val)
{
  for(int idir = 0; idir < DIM; idir++)
  {
    a_W(idir) = a_U(idir) + a_V(idir);
    if(a_W(idir) != a_val)
    {
      printf("w: values do not match\n");
    }
  }

}
PROTO_KERNEL_END(WsetWtoUplusVF, WsetWtoUplusV)

PROTO_KERNEL_START 
void setWtoUplusVptF(int a_p[DIM],
                     V a_W,
                     V a_U,
                     V a_V,
                     double  a_val)
{
//  printf("setw: uptr[0] = %p, uptr[1] = %p\n" ,a_U.m_ptrs[0],a_U.m_ptrs[1]);
//  printf("setw: vptr[0] = %p, vptr[1] = %p\n" ,a_V.m_ptrs[0],a_V.m_ptrs[1]);
//  printf("setw: wptr[0] = %p, wptr[1] = %p\n" ,a_W.m_ptrs[0],a_W.m_ptrs[1]);
//  printf("in set W \n");
  for(int idir = 0; idir < DIM; idir++)
  {
    a_W(idir) = a_U(idir) + a_V(idir);
    if(a_W(idir) != a_val)
    {
      printf("wpt: values do not match\n");
    }
  }
}
PROTO_KERNEL_END(setWtoUplusVptF, setWtoUplusVpt)

int main(int argc, char* argv[])
{

  int nx = 16;
  int maxGrid = nx; //this test is for a single grid
  Box domain(Point::Zeros(), Point::Ones(nx-1));
  double dx = 1.0/domain.size(0);
  std::array<bool, DIM> periodic;
  for(int idir = 0; idir < DIM; idir++) periodic[idir]=true;
  DisjointBoxLayout grids(domain, maxGrid, periodic);

  Point dataGhost = Point::Zeros();
  Point geomGhost = Point::Zeros();

  RealVect ABC = RealVect::Unit();
  RealVect X0 = RealVect::Unit();
//  X0 *= 100.;
  X0 *= 0.5;
  RealVect origin= RealVect::Zero();
  double R = 0.25;
  shared_ptr<BaseIF>              impfunc(new SimpleEllipsoidIF(ABC, X0, R, false));
  shared_ptr<GeometryService<2> > geoserv(new GeometryService<2>(impfunc, origin, dx, domain, grids, geomGhost, 0));
  shared_ptr<LevelData<EBGraph> > graphs = geoserv->getGraphs(domain);

  for(int ibox = 0; ibox < grids.size(); ibox++)
  {
    Box grid = grids[ibox];
    EBBoxData<CELL, double, DIM> U(grid,(*graphs)[ibox]);
    EBBoxData<CELL, double, DIM> V(grid,(*graphs)[ibox]);
    EBBoxData<CELL, double, DIM> W(grid,(*graphs)[ibox]);
    unsigned long long int numFlopsPt = 0;
    double uval = 1;
    double vval = 2;
    printf("going into setU\n");
    ebforallInPlace(numFlopsPt, "setU", UsetU, grid, U, uval);
    printf("going into setV\n");
    int vvar = -1;
    ebforallInPlace(numFlopsPt, "setV", VsetV, grid, V, vval, vvar);  //tweaking signature to clarify compilers job
    double wval = 3;
    printf("going into setWtoUPlusV\n");
    ebforallInPlace(numFlopsPt, "setWtoUPlusV", WsetWtoUplusV, grid, W, U, V, wval);

    uval = 2;
    vval = 5;
    wval = 7;
    printf("going into setUpt\n");
    ebforallInPlace_i(numFlopsPt, "setU", setUpt, grid, U, uval);
    printf("going into setVpt\n");
    ebforallInPlace_i(numFlopsPt, "setV", setVpt, grid, V, vval, vvar);
    printf("going into setWpt\n");
    ebforallInPlace_i(numFlopsPt, "setWtoUPlusV", setWtoUplusVpt, grid, W, U, V, wval);

  }

}








