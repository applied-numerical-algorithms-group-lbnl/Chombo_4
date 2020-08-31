/******************GIMRT98     ************************
      
!     Code converted using TO_F90 by Alan Miller
!     Date: 2000-07-27  Time: 10:02:15
      
!**************(C) COPYRIGHT 1995,1998,1999 ******************
!*******************C.I. Steefel      *******************
!     All Rights Reserved

!     GIMRT98 IS PROVIDED "AS IS" AND WITHOUT ANY WARRANTY EXPRESS OR IMPLIED.
!     THE USER ASSUMES ALL RISKS OF USING GIMRT98. THERE IS NO CLAIM OF THE
!     MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE.

!     YOU MAY MODIFY THE SOURCE CODE FOR YOUR OWN USE, BUT YOU MAY NOT
!     DISTRIBUTE EITHER THE ORIGINAL OR THE MODIFIED CODE TO ANY OTHER
!     WORKSTATIONS
!**********************************************************************/

#include <math.h>
#include "ShapeArray.H"
#include "crunchflow.h"

using namespace shape;

void reaction(const int ncomp,
	      const int neqn, 
	      const int nkin, 
	      const int nspec, 
	      const int jx, 
	      const int jy, 
	      const int jz,
	      const int nx,
	      const int ny, 
	      const int nz, 
	      double *snorm_2d, // [nkin][MAX_PATH], 
	      const int np, // number of pathways 
	      double *dppt_4d, // [nz][ny][nx][nkin], 
	      double *u_rate_4d, // [nz][ny][nx][nkin], 
	      /* const */ int *nreactmin_1d, // [nkin], 
	      double *actenergy_2d, // [nkin][MAX_PATH],
	      double *silog_2d, // [nkin][MAX_PATH], 
	      double *siln_2d, // [nkin][MAX_PATH],
	      double *si_2d, // [nkin][MAX_PATH], 
	      double *surf_1d, // [nkin], 
	      double *area_4d, // [nz][ny][nx][nkin], 
	      double *rmin_2d, // [nkin][MAX_PATH], 
	      double *pre_rmin_2d, // [nkin][MAX_PATH],
	      /* const */ int *ndepend_2d, // [nkin][np], 
	      const int nd, 
	      /* const */ int *idepend_3d, // [nkin][np][nd], 
	      double *depend_3d, // [nkin][MAX_PATH][MAX_DEPEND],
	      /* const */ int *itot_min_3d, // [nkin][np][nd],
	      /* const */ double *s_4d, // [nz][ny][nx][neqn], 
	      /* const */ double *gam_4d, // [nz][ny][nx][ncomp + nspec], 
	      double *sp_4d, // [nz][ny][nx][ncomp + nspec],
	      double *AffinityDepend1_2d, // [nkin][MAX_PATH], 
	      double *rate0_2d, // [nkin][MAX_PATH], 
	      double *volmol, 
	      double *mumin_3d, // [ncomp][nkin][MAX_PATH],
	      /* const */ double *keqmin_5d) // [nz][ny][nx][nkin][MAX_PATH])
{
  ShapeArray<double, 2> snorm(snorm_2d, nkin, MAX_PATH);
  ShapeArray<double, 4> dppt(dppt_4d, nz, ny, nx, nkin);
  ShapeArray<double, 4> u_rate(u_rate_4d, nz, ny, nx, nkin);
  ShapeArray<int, 1> nreactmin(nreactmin_1d, nkin);
  ShapeArray<double, 2> actenergy(actenergy_2d, nkin, MAX_PATH);
  ShapeArray<double, 2> silog(silog_2d, nkin, MAX_PATH); 
  ShapeArray<double, 2> siln(siln_2d, nkin, MAX_PATH);
  ShapeArray<double, 2> si(si_2d, nkin, MAX_PATH); 
  ShapeArray<double, 1> surf(surf_1d, nkin);
  ShapeArray<double, 4> area(area_4d, nz, ny, nx, nkin); 
  ShapeArray<double, 2> rmin(rmin_2d, nkin, MAX_PATH); 
  ShapeArray<double, 2> pre_rmin(pre_rmin_2d, nkin, MAX_PATH);
  ShapeArray<int, 2> ndepend(ndepend_2d, nkin, np); 
  ShapeArray<int, 3> idepend(idepend_3d, nkin, np, nd); 
  ShapeArray<double, 3> depend(depend_3d, nkin, MAX_PATH, MAX_DEPEND);
  ShapeArray<int, 3> itot_min(itot_min_3d, nkin, np, nd);
  ShapeArray<double, 4> s(s_4d, nz, ny, nx, neqn); 
  ShapeArray<double, 4> gam(gam_4d, nz, ny, nx, ncomp + nspec); 
  ShapeArray<double, 4> sp(sp_4d, nz, ny, nx, ncomp + nspec);
  ShapeArray<double, 2> AffinityDepend1(AffinityDepend1_2d, nkin, MAX_PATH); 
  ShapeArray<double, 2> rate0(rate0_2d, nkin, MAX_PATH); 
  ShapeArray<double, 3> mumin(mumin_3d, ncomp, nkin, MAX_PATH);
  ShapeArray<double, 5> keqmin(keqmin_5d, nz, ny, nx, nkin, MAX_PATH);
  
  double sumiap,term2,sign,term1,affinityTerm;
  int i;

  for(size_t i = 0; i < nkin; ++i)
    for(size_t j = 0; j < MAX_PATH; ++j)
      snorm[i][j] = 0.0;

  for (int k = 1; k <= nkin; ++k) 
    {
      dppt[jz-1][jy-1][jx-1][k-1] = 0.0;
      u_rate[jz-1][jy-1][jx-1][k-1] = 0.0;

      for (int npIndex = 1; npIndex <= nreactmin[k-1]; ++npIndex) 
	{
	  actenergy[k-1][npIndex-1] = 0.0;
	  sumiap = 0.0;

	  for (i = 1; i <= ncomp; ++i) {
	    sumiap = sumiap + mumin[i-1][k-1][npIndex-1] *
	      ( sp[jz-1][jy-1][jx-1][i-1] + gam[jz-1][jy-1][jx-1][i-1] );
	  }

	  silog[k-1][npIndex-1] = (sumiap - keqmin[jz-1][jy-1][jx-1][k-1][npIndex-1])/M_LN10;
	  siln[k-1][npIndex-1] = silog[k-1][npIndex-1] * M_LN10;
	  si[k-1][npIndex-1] = pow(10.0,silog[k-1][npIndex-1]);
	  surf[k-1] = area[jz-1][jy-1][jx-1][k-1];

	  if (surf[k-1] == 0.0) 
	    {
	      rmin[k-1][npIndex-1] = 0.0;
	      dppt[jz-1][jy-1][jx-1][k-1] = 0.0;
	      pre_rmin[k-1][npIndex-1] = 1.0;
	      u_rate[jz-1][jy-1][jx-1][k-1] = 0.0;
	    } 
	  else 
	    {
	      term2 = 0.0;
	      int kkMax = ndepend[k-1][npIndex-1];
	      for (int kk = 1; kk <= kkMax; ++kk) 
		{
		  i = idepend[k-1][npIndex-1][kk-1];

		  if (depend[k-1][npIndex-1][kk-1] < 0.0) 
		    {
		      if ( itot_min[k-1][npIndex-1][kk-1] == 1) 
			term2 = term2 + depend[k-1][npIndex-1][kk-1] * log(s[jz-1][jy-1][jx-1][i-1]);
		      else 
			term2 = term2 + depend[k-1][npIndex-1][kk-1] *
			  (gam[jz-1][jy-1][jx-1][i-1] + sp[jz-1][jy-1][jx-1][i-1]);
		    } 
		  else 
		    {
		      if ( itot_min[k-1][npIndex-1][kk-1] == 1) 
			term2 = term2 + depend[k-1][npIndex-1][kk-1] * log(s[jz-1][jy-1][jx-1][i-1]);
		      else 
			term2 = term2 + depend[k-1][npIndex-1][kk-1] * 
			  (gam[jz-1][jy-1][jx-1][i-1] + sp[jz-1][jy-1][jx-1][i-1]);
		    }
		} // for (int kk = 1; kk <= ndepend[k-1][npIndex-1]; ++kk)
	      if (term2 == 0.0) {
		pre_rmin[k-1][npIndex-1] = 1.0;
	      } else {
		pre_rmin[k-1][npIndex-1] = exp(term2);
	      }
	      if (si[k-1][npIndex-1] > 1.0) {
		sign = 1.0;
	      } else {
		sign = -1.0;
	      }
	      snorm[k-1][npIndex-1] = si[k-1][npIndex-1];
	      if (AffinityDepend1[k-1][npIndex-1] == 1.0) {
		term1 = sign * fabs(snorm[k-1][npIndex-1] - 1.0);
	      } else {
		term1 = sign * fabs( pow( (snorm[k-1][npIndex-1] - 1) , (AffinityDepend1[k-1][npIndex-1]) ) );
	      }
	      affinityTerm = term1;
	      rmin[k-1][npIndex-1] = surf[k-1] * rate0[k-1][npIndex-1] * pre_rmin[k-1][npIndex-1] * affinityTerm;
	      dppt[jz-1][jy-1][jx-1][k-1] = dppt[jz-1][jy-1][jx-1][k-1] + rmin[k-1][npIndex-1];
	      u_rate[jz-1][jy-1][jx-1][k-1] = u_rate[jz-1][jy-1][jx-1][k-1] + rate0[k-1][npIndex-1] *
		pre_rmin[k-1][npIndex-1] * affinityTerm * volmol[k-1];
	    }
	}
    }
  return;
}
