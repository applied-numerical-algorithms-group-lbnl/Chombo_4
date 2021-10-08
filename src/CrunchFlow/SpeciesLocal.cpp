#include <stdio.h>
#include <math.h>
#include "ShapeArray.H"
#include "crunchflow.h"

using namespace shape;

void specieslocal(const int ncomp, 
		  const int nspec, 
		  const int jx, 
		  const int jy, 
		  const int jz, 
		  const int nx, 
		  const int ny, 
		  const int nz,
		  double *muaq_2d, 
		  double *sp_4d, 
		  double *gam_4d, 
		  double *sp10_4d, 
		  double *keqaq_4d)
{
  ShapeArray<double, 2> muaq(muaq_2d, ncomp, nspec);
  ShapeArray<double, 4> sp(sp_4d, nz, ny, nx, ncomp + nspec);
  ShapeArray<double, 4> gam(gam_4d, nz, ny, nx, ncomp + nspec);
  ShapeArray<double, 4> sp10(sp10_4d, nz, ny, nx, ncomp + nspec); 
  ShapeArray<double, 4> keqaq(keqaq_4d, nz, ny, nx, nspec);

  int ksp = 0;
  int i = 0;
  int nk = 0;
  double sum = 0;
  for (ksp = 1; ksp <= nspec; ksp++) 
    {
      sum = 0;
      for (i = 1; i <= ncomp; i++) 
	sum = sum + muaq[i-1][ksp-1] * ( sp[jz-1][jy-1][jx-1][i-1] + gam[jz-1][jy-1][jx-1][i-1] );

      nk = ncomp + ksp;
      sp[jz-1][jy-1][jx-1][nk-1] = keqaq[jz-1][jy-1][jx-1][ksp-1] - gam[jz-1][jy-1][jx-1][nk-1] + sum;
      sp10[jz-1][jy-1][jx-1][nk-1] = exp(sp[jz-1][jy-1][jx-1][nk-1]);
    }
  return;
}

