#include "ShapeArray.H"
#include "crunchflow.h"

using namespace shape;

void totconc(const int ncomp, 
	     const int nspec, 
	     const int neqn,
	     const int jx, 
	     const int jy, 
	     const int jz, 
	     const int nx, 
	     const int ny, 
	     const int nz, 
	     /* const */ double *muaq_2d, // [ncomp][nspec], 
	     /* const */ double *sp10_4d, // [nz][ny][nx][ncomp + nspec], 
	     double *s_4d) // [nz][ny][nx][ncomp]
{
  ShapeArray<double, 2> muaq(muaq_2d, ncomp, nspec);
  ShapeArray<double, 4> sp10(sp10_4d, nz, ny, nx, ncomp + nspec);
  ShapeArray<double, 4> s(s_4d, nz, ny, nx, ncomp);
  
  int kk;
  double sum;

  for (int i = 1; i <= ncomp; ++i) {
    sum = 0.0;
    for (int ksp = 1; ksp <= nspec; ++ksp) {
      kk = ksp + ncomp;
      sum = sum + muaq[i-1][ksp-1] * sp10[jz-1][jy-1][jx-1][kk-1];
    }
    s[jz-1][jy-1][jx-1][i-1] = sum + sp10[jz-1][jy-1][jx-1][i-1];
  }
  return;
}
