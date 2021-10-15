#include <math.h>
#include "ShapeArray.H"
#include "crunchflow.h"

using namespace shape;

void keqcalc2(const int ncomp, 
	      const int nrct, 
	      const int nspec, 
	      const int ngas, 
	      const int nkin,
	      const int ntemp,
	      const int nsurf_sec,
	      const int jx, 
	      const int jy, 
	      const int jz, 
	      const int nx, 
	      const int ny, 
	      const int nz,
	      /* const */ double *t_3d, // [nz][ny][nx], 
	      const bool RunIsoThermal, 
	      double *keqaq_4d, // [nz][ny][nx][nspec],
	      /* const */ int *nreactmin_1d, // [nkin], 
	      /* const */ double *keqmin_5d, // [nz][ny][nx][nkin][MAX_PATH],
	      /* const */ double *alnk_1d, // [nkin * MAX_PATH],
	      /* const */ double *eqhom_1d, // [nspec],
	      /* const */ double *as1_2d) // [NBASIS][nspec + nkin * MAX_PATH])
{
  ShapeArray<double, 3> t(t_3d, nz, ny, nx);
  ShapeArray<double, 4> keqaq(keqaq_4d, nz, ny, nx, nspec);
  ShapeArray<int, 1> nreactmin(nreactmin_1d, nkin);
  ShapeArray<double, 5> keqmin(keqmin_5d, nz, ny, nx, nkin, MAX_PATH);
  ShapeArray<double, 1> alnk(alnk_1d, nkin * MAX_PATH);
  ShapeArray<double, 1> eqhom(eqhom_1d, nspec);
  ShapeArray<double, 2> as1(as1_2d, NBASIS, nspec + nkin * MAX_PATH);

  double temp, temp2, x1, x2, x3, x4, x5;
  int msub,ksp;

  temp = t[jz-1][jy-1][jx-1] + 273.15;
  temp2 = temp * temp;
  for (int ksp = 1; ksp <= nspec; ++ksp) 
    {
      if ((ntemp == 1) || RunIsoThermal) 
	{
	  keqaq[jz-1][jy-1][jx-1][ksp-1] = -1 * M_LN10 * eqhom[ksp - 1];
	} 
      else 
	{
	  x1 = as1[1-1][ksp-1];
	  x2 = as1[2-1][ksp-1];
	  x3 = as1[3-1][ksp-1];
	  x4 = as1[4-1][ksp-1];
	  x5 = as1[5-1][ksp-1];
	  keqaq[jz-1][jy-1][jx-1][ksp-1] = -1*M_LN10*( x1*log(temp) + x2 + x3*temp + x4/temp + x5/temp2 );
	}
    }
  msub = 0;

  for (int k = 1; k <= nrct; ++k) 
    {
      for (int npIndex = 1; npIndex <= nreactmin[k-1]; ++npIndex) 
	{
	  msub = msub + 1;
	  ksp = msub + ngas + nspec;
	  if ((ntemp == 1) || RunIsoThermal) 
	    {
	      keqmin[jz-1][jy-1][jx-1][k-1][npIndex-1] = M_LN10 * alnk[msub - 1];
	    } 
	  else 
	    {
	      x1 = as1[1-1][ksp-1];
	      x2 = as1[2-1][ksp-1];
	      x3 = as1[3-1][ksp-1];
	      x4 = as1[4-1][ksp-1];
	      x5 = as1[5-1][ksp-1];
	      keqmin[jz-1][jy-1][jx-1][k-1][npIndex-1] = M_LN10*( x1*log(temp) + x2 + x3*temp + x4/temp + x5/temp2 );
	    }
	}
    }
  return;
}
