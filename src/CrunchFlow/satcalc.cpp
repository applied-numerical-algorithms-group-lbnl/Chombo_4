/******************        GIMRT98     ************************
 
! Code converted using TO_F90 by Alan Miller
! Date: 2000-07-27  Time: 10:07:50
 
!************** (C) COPYRIGHT 1995,1998,1999 ******************
!*******************     C.I. Steefel      *******************
!                    All Rights Reserved

!  GIMRT98 IS PROVIDED "AS IS" AND WITHOUT ANY WARRANTY EXPRESS OR IMPLIED.
!  THE USER ASSUMES ALL RISKS OF USING GIMRT98. THERE IS NO CLAIM OF THE
!  MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE.

!  YOU MAY MODIFY THE SOURCE CODE FOR YOUR OWN USE, BUT YOU MAY NOT
!  DISTRIBUTE EITHER THE ORIGINAL OR THE MODIFIED CODE TO ANY OTHER
!  WORKSTATIONS
!**********************************************************************/

#include <math.h>
#include "ShapeArray.H"
#include "crunchflow.h"

using namespace shape;

void satcalc(const int ncomp,
	     const int nspec,
	     const int nkin,
	     const int nrct,
	     const int jx,
	     const int jy,
	     const int jz,
	     const int nx,
	     const int ny,
	     const int nz,
	     double *silog_2d, // [nkin][MAX_PATH],
	     double *mumin_3d, // [ncomp][nkin][MAX_PATH],
	     double *sp_4d, // [nz][ny][nx][ncomp + nspec],
	     double *gam_4d, // [nz][ny][nx][ncomp + nspec],
	     /* const */ double *keqmin_5d, // [nz][ny][nx][nkin][MAX_PATH],
	     /* const */ int *nreactmin_1d) // [nkin]
{
  ShapeArray<double, 2> silog(silog_2d, nkin, MAX_PATH);
  ShapeArray<double, 3> mumin(mumin_3d, ncomp, nkin, MAX_PATH);
  ShapeArray<double, 4> sp(sp_4d, nz, ny, nx, ncomp + nspec);
  ShapeArray<double, 4> gam(gam_4d, nz, ny, nx, ncomp + nspec);
  ShapeArray<double, 5> keqmin(keqmin_5d, nz, ny, nx, nkin, MAX_PATH);
  ShapeArray<int, 1> nreactmin(nreactmin_1d, nkin);
  
  double sumiap;
  for (int k = 1; k <= nrct; ++k) {
    for (int np = 1; np <= nreactmin[k-1]; ++np) {
      sumiap = 0.0;
      for (int i = 1; i <= ncomp; ++i) {
	sumiap = sumiap + mumin[i-1][k-1][np-1] *
	  ( sp[jz-1][jy-1][jx-1][i-1] + gam[jz-1][jy-1][jx-1][i-1] );
      }
      silog[k-1][np-1] = (sumiap - keqmin[jz-1][jy-1][jx-1][k-1][np-1])/M_LN10;
    }
  }
}
