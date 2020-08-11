#include <math.h>
#include "ShapeArray.H"
#include "crunchflow.h"
 
using namespace shape;

double affinitynumerical(const int ncomp, 
			 const int jx, 
			 const int jy, 
			 const int jz, 
			 const int nx, 
			 const int ny,
			 const int nz,
			 const int np, 
			 const int nkin, 
			 const int nspec,
			 const int loopNP, 
			 int k, 
			 double *sppTMP_1d, // [ncomp + nspec]
			 double *gam_4d, // [nz][ny][nx][ncomp + nspec],
			 double *mumin_3d, // [ncomp][nkin][MAX_PATH], 
			 double *keqmin_5d, // [nz][ny][nx][nkin][MAX_PATH], 
			 double *AffinityDepend1_2d) // [nkin][MAX_PATH]
{
  ShapeArray<double, 1> sppTMP(sppTMP_1d, ncomp + nspec);
  ShapeArray<double, 4> gam(gam_4d, nz, ny, nx, ncomp + nspec);
  ShapeArray<double, 3> mumin(mumin_3d, ncomp, nkin, MAX_PATH);
  ShapeArray<double, 5> keqmin(keqmin_5d, nz, ny, nx, nkin, MAX_PATH );
  ShapeArray<double, 2> AffinityDepend1(AffinityDepend1_2d, nkin, MAX_PATH);

  double sumiap, silogTMP, silnTMP, snormTMP, siTMP, power, term1, sign;
  int i = 1;

  sumiap = 0;

  for (i = 1; i <= ncomp; i++) 
    {
      sumiap = sumiap + mumin[i-1][k-1][loopNP-1] * 
	(sppTMP[i-1] + gam[jz-1][jy-1][jx-1][i-1]);
    }

  silogTMP = (sumiap - keqmin[jz-1][jy-1][jx-1][k-1][loopNP-1])/M_LN10;
  siTMP = pow(10.0,silogTMP);
  silnTMP = silogTMP * M_LN10;

  if (siTMP > 1) {
    sign = 1.0;
  } else {
    sign = -1.0;
  }

  snormTMP = siTMP;

  if (AffinityDepend1[k-1][loopNP-1] == 1.0) {
    term1 = sign * fabs(snormTMP - 1.0);
  } else {
    term1 = sign * fabs( pow( (snormTMP - 1.0) , (AffinityDepend1[k-1][loopNP-1]) ) );
  }
  return term1;
}
