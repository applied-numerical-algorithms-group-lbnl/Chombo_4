#include <gtest/gtest.h>
#include "Proto.H"
#include "EBProto.H"
#include "Chombo_EBChombo.H"
#include "Chombo_REAL.H"
#include "Chombo_Box.H"
#include "Chombo_IntVect.H"
#include "Chombo_EBLevelBoxData.H"
#include "Chombo_GeometryService.H"
using Chombo4::LevelData;

///  after this are specific to the test
typedef Var<Real,DIM> V;


PROTO_KERNEL_START 
void UsetUF(V a_U, Real  a_val)
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
void setUptF(int  a_p[DIM], V a_U, Real  a_val)
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
void VsetVF(V a_V, Real  a_val, int a_intvar)
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
void setVptF(int  a_p[DIM], V a_V, Real  a_val, int a_vvar)
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
                    Real  a_val)
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
                     Real  a_val)
{
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

int checkAns(EBBoxData<CELL, Real, DIM>& calc,
             EBBoxData<CELL, Real, DIM>& ref,
             Real uval,
             Chombo4::Box grid)
{
  Bx bx = ProtoCh::getProtoBox(grid);
  ref.setVal(uval);
  Real tol = 1.0e-6;
  for(int idir = 0; idir < DIM; idir++)
  {
    int compval = EBLevelBoxData<CELL, DIM>::checkAnswer(calc, ref, bx, idir, tol);
    if(compval != 0)
    {
      return compval;
    }
  }
  return 0;
}
using Chombo4::Box;
using Chombo4::DisjointBoxLayout;
int nx = 16;

Box domain(IntVect::Zero, (nx-1)*IntVect::Unit);
Vector<Box> boxes(1, domain);
Vector<int> procs(1, 0);
Real dx = 1.0/domain.size(0);
std::array<bool, DIM> periodic;
for(int idir = 0; idir < DIM; idir++) periodic[idir]=true;
DisjointBoxLayout grids(boxes, procs);

int geomGhost = 0;

RealVect ABC = RealVect::Unit();
RealVect X0  = RealVect::Unit();
X0 *= 0.5;

RealVect origin= RealVect::Zero();
Real R = 0.25;
shared_ptr<BaseIF>              impfunc(new SimpleEllipsoidIF(ABC, X0, R, false));
shared_ptr<GeometryService<2> > geoserv(new GeometryService<2>(impfunc, origin, dx, domain, grids, geomGhost));
shared_ptr<Chombo4::LevelData<EBGraph>  > graphs = geoserv->getGraphs(domain);
Chombo4::DataIterator dit = grids.dataIterator();

TEST(testLibEBForAll,checkU) {
  int cval;
  for(int ibox = 0; ibox < dit.size(); ibox++)
  {
    Box grid = grids[dit[ibox]];
    Bx grbx = ProtoCh::getProtoBox(grid);
    EBBoxData<CELL, Real, DIM> U(grbx,(*graphs)[dit[ibox]]);
    EBBoxData<CELL, Real, DIM> reference(grbx,(*graphs)[dit[ibox]]);

    Real uval = 1;
    ebforall(grbx, UsetU, grbx, U, uval);
    cval = checkAns(U, reference, uval, grid);
    if(cval != 0)
      return;
  }
  EXPECT_EQ(cval,0);
}

TEST(testLibEBForAll,checkV) {
  int dval;
  for(int ibox = 0; ibox < dit.size(); ibox++)
  {
    Box grid = grids[dit[ibox]];
    Bx grbx = ProtoCh::getProtoBox(grid);
    EBBoxData<CELL, Real, DIM> V(grbx,(*graphs)[dit[ibox]]);
    EBBoxData<CELL, Real, DIM> reference(grbx,(*graphs)[dit[ibox]]);

    Real vval = 2;
    int vvar = -1;
    ebforall(grbx, VsetV, grbx, V, vval, vvar);  //tweaking signature to clarify compilers job
    dval = checkAns(V, reference, vval, grid);
    if(dval != 0)
      return;
  }
  EXPECT_EQ(dval,0);
}

TEST(testLibEBForAll,checkW) {
  int eval;
  for(int ibox = 0; ibox < dit.size(); ibox++)
  {
    Box grid = grids[dit[ibox]];
    Bx grbx = ProtoCh::getProtoBox(grid);
    EBBoxData<CELL, Real, DIM> U(grbx,(*graphs)[dit[ibox]]);
    EBBoxData<CELL, Real, DIM> V(grbx,(*graphs)[dit[ibox]]);
    EBBoxData<CELL, Real, DIM> W(grbx,(*graphs)[dit[ibox]]);
    EBBoxData<CELL, Real, DIM> reference(grbx,(*graphs)[dit[ibox]]);

    Real wval = 3;
    ebforall(grbx, WsetWtoUplusV, grbx, W, U, V, wval);
    eval = checkAns(W, reference, wval, grid);
    if(eval != 0)
      return;
  }
  EXPECT_EQ(eval,0);
}

TEST(testLibEBForAll,checkUpt) {
  int fval;
  for(int ibox = 0; ibox < dit.size(); ibox++)
  {
    Box grid = grids[dit[ibox]];
    Bx grbx = ProtoCh::getProtoBox(grid);
    EBBoxData<CELL, Real, DIM> U(grbx,(*graphs)[dit[ibox]]);
    EBBoxData<CELL, Real, DIM> reference(grbx,(*graphs)[dit[ibox]]);

    Real uval = 2;
    //printf("going into setUpt\n");
    ebforall_i(grbx, setUpt, grbx, U, uval);
    fval = checkAns(U, reference, uval, grid);
    if(fval != 0)
      return;
  }
  EXPECT_EQ(fval,0);
}

TEST(testLibEBForAll,checkVpt) {
  int gval;
  for(int ibox = 0; ibox < dit.size(); ibox++)
  {
    Box grid = grids[dit[ibox]];
    Bx grbx = ProtoCh::getProtoBox(grid);
    EBBoxData<CELL, Real, DIM> V(grbx,(*graphs)[dit[ibox]]);
    EBBoxData<CELL, Real, DIM> reference(grbx,(*graphs)[dit[ibox]]);

    Real vval = 5;
    int vvar = -1;
    //printf("going into setVpt\n");
    ebforall_i(grbx, setVpt, grbx, V, vval, vvar);
    gval = checkAns(V, reference, vval, grid);
    if(gval != 0)
      return;
  }
  EXPECT_EQ(gval,0);
}

TEST(testLibEBForAll,checkWpt) {
  int hval;
  for(int ibox = 0; ibox < dit.size(); ibox++)
  {
    Box grid = grids[dit[ibox]];
    Bx grbx = ProtoCh::getProtoBox(grid);
    EBBoxData<CELL, Real, DIM> U(grbx,(*graphs)[dit[ibox]]);
    EBBoxData<CELL, Real, DIM> V(grbx,(*graphs)[dit[ibox]]);
    EBBoxData<CELL, Real, DIM> W(grbx,(*graphs)[dit[ibox]]);
    EBBoxData<CELL, Real, DIM> reference(grbx,(*graphs)[dit[ibox]]);

    Real wval = 7;
    //printf("going into setWpt\n");
    ebforall_i(grbx, setWtoUplusVpt, grbx, W, U, V, wval);
    hval = checkAns(W, reference, wval, grid);
    if(hval != 0)
      return;
  }
  EXPECT_EQ(hval,0);
}
int main(int argc, char* argv[])
{
  ::testing::InitGoogleTest(&argc,argv);
  RUN_ALL_TESTS();
  return 0;
}
