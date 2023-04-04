#include <math.h>
#include <string.h>
#include "ShapeArray.H"
#include "crunchflow.h"

using namespace shape;

void cfgamma(const int ncomp, 
	     const int nspec,
	     const int ntemp,
	     const int jx, 
	     const int jy, 
	     const int jz, 
	     const int nx, 
	     const int ny,
	     const int nz, 
	     /* const */ double *sp10_4d, // [nz][ny][nx][ncomp + nspec], 
	     /* const */ double *t_3d, // [nz][ny][nx], 
	     /* const */ double *adh_1d,
	     /* const */ double *bdh_1d,
	     /* const */ double *bdot_1d,
	     double *adhcoeff_1d, // [NBASIS],
	     double *bdhcoeff_1d, // [NBASIS],
	     double *bdtcoeff_1d, // [NBASIS],
	     double *sion_3d, // [nz][ny][nx],
	     char **ulab_1d, // ncomp+nspec
	     /* const */ double *acmp_1d,
	     /* const */ double *chg_1d,
	     double *gam_4d)
{
  ShapeArray<double, 4> sp10(sp10_4d, nz, ny, nx, ncomp + nspec);
  ShapeArray<double, 3> t(t_3d, nz, ny, nx);
  ShapeArray<double, 1> adh(adh_1d, TPOINTS);
  ShapeArray<double, 1> bdh(bdh_1d, TPOINTS);
  ShapeArray<double, 1> bdot(bdot_1d, TPOINTS); 
  ShapeArray<double, 1> adhcoeff(adhcoeff_1d, NBASIS);
  ShapeArray<double, 1> bdhcoeff(bdhcoeff_1d, NBASIS);
  ShapeArray<double, 1> bdtcoeff(bdtcoeff_1d, NBASIS);
  ShapeArray<double, 3> sion(sion_3d, nz, ny, nx);
  ShapeArray<char *, 1> ulab(ulab_1d, ncomp + nspec); 
  ShapeArray<double, 1> acmp(acmp_1d, ncomp + nspec); 
  ShapeArray<double, 1> chg(chg_1d, ncomp + nspec); 
  ShapeArray<double, 4> gam(gam_4d, nz, ny, nx, ncomp + nspec);

  double ctotal,ah,bh,bdt,sum,aa1,tempc;
  ctotal = 0.0;
  for (int ik = 0; ik < (ncomp + nspec); ++ik) {
    ctotal = ctotal + sp10[jz][jy][jx][ik];
  }
  tempc = t[jz][jy][jx];
  if (ntemp == 1) {
    ah = adh[1-1];
    bh = bdh[1-1];
    bdt = bdot[1-1];
  } else {
    ah = adhcoeff[1-1] +
      adhcoeff[2-1] * tempc +
      adhcoeff[3-1] * tempc * tempc +
      adhcoeff[4-1] * tempc * tempc * tempc +
      adhcoeff[5-1] * tempc * tempc * tempc * tempc;
    bh = bdhcoeff[1-1] +
      bdhcoeff[2-1] * tempc +
      bdhcoeff[3-1] * tempc * tempc +
      bdhcoeff[4-1] * tempc * tempc * tempc +
      bdhcoeff[5-1] * tempc * tempc * tempc * tempc;
    bdt = bdtcoeff[1-1] +
      bdtcoeff[2-1] * tempc +
      bdtcoeff[3-1] * tempc * tempc +
      bdtcoeff[4-1] * tempc * tempc * tempc +
      bdtcoeff[5-1] * tempc * tempc * tempc * tempc;
  }

  sum = 0.0;
  for (int ik = 0; ik < (ncomp + nspec); ++ik) {
    sum = sum + sp10[jz][jy][jx][ik] * chg[ik] * chg[ik];
  }

  sion[jz][jy][jx] = 0.5 * sum;
  if (sion[jz][jy][jx] > 25.0) {
    sion[jz][jy][jx] = 0.0;
  }
    
  for (int ik = 0; ik < (ncomp + nspec); ++ik) {
    if (chg[ik] == 0.0) {
      if (ulab[ik] != NULL && strncmp(ulab[ik],"H2O",3) == 0) {
	gam[jz][jy][jx][ik] = log(1.0/ctotal);
      } else if (ulab[ik] != NULL && strncmp(ulab[ik],"SiO2(aq)",8) == 0) {
	gam[jz][jy][jx][ik] = 0.08 * sion[jz][jy][jx];
      } else {
	gam[jz][jy][jx][ik] = 0.231 * sion[jz][jy][jx];
      }
    } else {
      /*aa1 = -1 * (ah * chg[ik] * chg[ik] * sqrt(array3D(sion,jx,jy,jz,nx,ny))) /
	(1.0 + array1D(acmp,ik) * bh * sqrt(array3D(sion,jx,jy,jz,nx,ny))) +
	bdt * array3D(sion,jx,jy,jz,nx,ny);*/
      aa1 = -1 * (ah * chg[ik] * chg[ik] * sqrt(sion[jz][jy][jx]) /
		  (1.0 + acmp[ik] * bh * sqrt(sion[jz][jy][jx])) +
		  bdt * sion[jz][jy][jx]);

      gam[jz][jy][jx][ik] = M_LN10 * aa1;
    }
  }
  return;
}

