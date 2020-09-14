/******************GIMRT98     ************************
      
!     Code converted using TO_F90 by Alan Miller
!     Date: 2000-07-27  Time: 09:52:46
      
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

#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include "ShapeArray.H"
#include "crunchflow.h"

using namespace shape;

void assemble_local(const int ncomp, 
		    const int nspec, 
		    const int nkin, 
		    const int ikin, 
		    const int neqn, 
		    const double dt,
		    const int nx, 
		    const int ny, 
		    const int nz,
		    const int jx, 
		    const int jy, 
		    const int jz, 
		    double *satliq_3d, // [nz][ny][nx], 
		    double *por_3d, // [nz][ny][nx],
		    /* const */ double *ro_3d, // [nz][ny][nx],
		    /* const */ double *rmin_2d, // [nkin][MAX_PATH], 
		    /* const */ double *decay_correct_2d, // [nkin][ncomp],
		    double *aaa_2d, // [ncomp][ncomp],
		    /* const */ int nradmax, 
		    /* const */ int *nreactmin_1d, // [nkin], 
		    const int np,
		    double *fxx_1d, // [neqn], 
		    /* const */ double *mumin_3d, // [ncomp][nkin][MAX_PATH],
		    /* const */ double *mukin_2d, // [ncomp][ikin],
		    double *raq_tot_4d, // [nz][ny][nx][ikin],
		    /* const */ double *jac_rmin_3d, // [nkin][MAX_PATH][ncomp],
		    double *rdkin_2d, // [ncomp][ikin],
		    const int ikh2o, 
		    double *distrib,
		    /* const */ double *H2Oreacted_3d, // [nz][ny][nx], 
		    double *fjac_loc_2d, // [neqn][neqn], 
		    /* const */ double *xgram_3d) // [nz+2][ny+2][nx+3])
{
  ShapeArray<double, 3> satliq(satliq_3d, nz, ny, nx); 
  ShapeArray<double, 3> por(por_3d, nz, ny, nx);
  ShapeArray<double, 3> ro(ro_3d, nz, ny, nx);
  ShapeArray<double, 2> rmin(rmin_2d, nkin, MAX_PATH); 
  ShapeArray<double, 2> decay_correct(decay_correct_2d, nkin, ncomp);
  ShapeArray<double, 2> aaa(aaa_2d, ncomp, ncomp);
  ShapeArray<int, 1>  nreactmin(nreactmin_1d, nkin); 
  ShapeArray<double, 1> fxx(fxx_1d, neqn); 
  ShapeArray<double, 3> mumin(mumin_3d, ncomp, nkin, MAX_PATH);
  ShapeArray<double, 2> mukin(mukin_2d, ncomp, ikin);
  ShapeArray<double, 4> raq_tot(raq_tot_4d, nz, ny, nx, ikin);
  ShapeArray<double, 3> jac_rmin(jac_rmin_3d, nkin, MAX_PATH, ncomp);
  ShapeArray<double, 2> rdkin(rdkin_2d, ncomp, ikin);
  ShapeArray<double, 3> H2Oreacted(H2Oreacted_3d, nz, ny, nx); 
  ShapeArray<double, 2> fjac_loc(fjac_loc_2d, neqn, neqn); 
  ShapeArray<double, 3> xgram(xgram_3d, nz+2, ny+2, nx+3);


  double source, retardation,satl, portemp, rotemp, xgtemp, sumrct, r, sumkin,
    rxnmin, rxnaq, aq_accum,ex_accum,source_jac;
  int ind;
  double sumrd[neqn]; // allocate(sumrd(neqn),stat=ierr); sumrd = 0.0d0
  double sumjackin[neqn]; // allocate(sumjackin(neqn),stat=ierr); sumjackin = 0.0d0

  for(size_t i = 0; i < ncomp; ++i)
    for(size_t j = 0; j < ncomp; ++j)
      aaa[i][j] = 0.0;

  source = 0.0;
  retardation = 1.0;

  satl    = satliq[jz-1][jy-1][jx-1];
  portemp = por[jz-1][jy-1][jx-1];
  rotemp  = ro[jz-1][jy-1][jx-1];
  xgtemp = xgram[jz-1][jy-1][jx-1];
  r = 1.0/dt;

  for (int i = 1; i <= ncomp; ++i) {
    ind = i;

    int bkpt = 0;
    if(jz == 1 && jy == 1 && jx == 1 && ind == 1)
      bkpt = 1;

    if (nradmax > 0) {
      sumrct = 0.0;
      for (int k = 1; k <= nkin; ++k) {
	for (int npIndex = 1; npIndex <= nreactmin[k-1]; ++npIndex) {
	  if (rmin[k-1][npIndex-1] >= 0.0) {

	    sumrct = sumrct + decay_correct[k-1][i-1] * mumin[i-1][k-1][npIndex-1] * rmin[k-1][npIndex-1];

	  } else {
	    // sumrct = sumrct + array2D(decay_correct,i,k,ncomp) * array6D(mumin_decay,np,k,i,jx,jy,jz) * 
	    // array2D(rmin,np,k,np);
	    fprintf(stderr, "assemble_local_cpp.cpp: mumin_decay not implemented in crunchflow lite\n");
	    exit(1);
	  }
	}
      }
    } else {
      sumrct = 0.0;
      for (int k = 1; k <= nkin; ++k) {
	for (int npIndex = 1; npIndex <= nreactmin[k-1]; ++npIndex) {
	  sumrct = sumrct + decay_correct[k-1][i-1] * mumin[i-1][k-1][npIndex-1] * rmin[k-1][npIndex-1];
	}
      }
    }

    sumkin = 0.0;
    for (int ir = 1; ir <= ikin; ++ir) {
      sumkin = sumkin - mukin[i-1][ir-1] * raq_tot[jz-1][jy-1][jx-1][ir-1];
    }

    // Update the residual, adding reaction terms and exchange terms
    //      array1D(fxx,ind) = array1D(fxx,ind) + sumrct + satl * xgtemp * portemp * rotemp * sumkin;
    fxx[ind-1] += sumrct + satl * xgtemp * portemp * rotemp * sumkin;

    for(size_t j = 0; j < neqn; ++j)
      {
	sumrd[j] = 0.0;
	sumjackin[j] = 0.0;
      }

    if (nradmax > 0) {
      for (int k = 1; k <= nkin; ++k) {
	for (int npIndex = 1; npIndex <= nreactmin[k-1]; ++npIndex) {
	  if (mumin[i-1][k-1][npIndex-1] != 0.0) {
	    if (rmin[k-1][npIndex-1] >= 0.0) 
	      {
		for (int i2 = 1; i2 <= ncomp; ++i2) 
		  {
		    sumrd[i2-1] = sumrd[i2-1] + decay_correct[k-1][i-1] * 
		      mumin[i-1][k-1][npIndex-1] * jac_rmin[k-1][npIndex-1][i2-1];
		  }
	      } 
	    else 
	      {
		for (int i2 = 1; i2 <= ncomp; ++i2) 
		  {
		    //sumrd[i2-1] = sumrd[i2-1] + decay_correct[k-1][i-1] * 
		    //array6D(mumin_decay,npIndex,k,i,jx,jy,jz)*jac_rmin[k-1][npIndex-1][i2-1];
		    fprintf(stderr, "assemble_local_cpp.cpp: mumin_decay not implemented in crunchflow lite\n");
		    exit(1);
		  }
	      }
	  }
	}
      }
    } else {
      for (int k = 1; k <= nkin; ++k) {
	for (int npIndex = 1; npIndex <= nreactmin[k-1]; ++npIndex) {
	  if (mumin[i-1][k-1][npIndex-1] != 0.0) {
	    for (int i2 = 1; i2 <= ncomp; ++i2) 
	      {
		sumrd[i2-1] = sumrd[i2-1] + decay_correct[k-1][i-1] * 
		  mumin[i-1][k-1][npIndex-1] * jac_rmin[k-1][npIndex-1][i2-1];
	      }

	  }
	}
      }
    }

    for (int ir = 1; ir <= ikin; ++ir) {
      if (mukin[i-1][ir-1] != 0.0) {
	for (int i2 = 1; i2 <= ncomp; ++i2) {
	  sumjackin[i2-1] = sumjackin[i2-1] - mukin[i-1][ir-1] * rdkin[i2-1][ir-1];
	}
      }
    }

    for (int i2 = 1; i2 <= ncomp; ++i2) 
      {
	rxnmin = sumrd[i2-1];
	rxnaq = satl*xgtemp*portemp*rotemp * sumjackin[i2-1];

	if (i != ikh2o) 
	  {
	    aq_accum = H2Oreacted[jz-1][jy-1][jx-1] * satl * xgtemp * r * portemp * rotemp * fjac_loc[i-1][i2-1] * 
	      (1.0 + retardation * distrib[i-1]);	      
	  } 
	else 
	  {
	    aq_accum = satl *	xgtemp * r * portemp * rotemp *	fjac_loc[i-1][i2-1] *
	      (1.0 + retardation * distrib[i-1]);
	  }
	source_jac = source * fjac_loc[i-1][i2-1];
	ex_accum = 0.0;
	aaa[i-1][i2-1] = rxnmin + rxnaq + aq_accum - source_jac + ex_accum;
      } // for (int i2 = 1; i2 <= ncomp; ++i2)...
  } // for (int i = 1; i <= ncomp; ++i)...
}
