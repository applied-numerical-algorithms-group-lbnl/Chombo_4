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

#ifdef PROTO_CUDA
#  include <assert.h>
#endif

#include "ShapeArray.H"
#include "crunchflow.h"

#define get2D(A,i,j,Ni,Nj) (A[(i)*(Nj) + (j)])
#define get3D(A,i,j,k,Ni,Nj,Nk) (A[(i) * (Nj) * (Nk) + (j) * (Nk) + (k)])

using namespace shape;

#ifdef PROTO_CUDA

__global__ void assemble_local_kernel(double* dA,
			 double* dSatliq,
			 const int neqn,
			 const int ncomp,
			 const int nkin,
			 const int nx,
			 const int ny,
			 const int nz,
			 const int nradmax,
			 const double dt,
			 int *nreactmin, // [nkin]
			 double *rmin, // [nkin][MAX_PATH]
			 double *decay_correct, // [nkin][ncomp]
			 double *mumin) // [ncomp][nkin][MAX_PATH]
{
  double r = 1.0/dt;
  double sumrct;
  double sumkin;
  int ind;

  for (int i = 1; i <= ncomp; ++i)
    {
      ind = i;

    if (nradmax > 0)
      {
	sumrct = 0.0;
	for (int k = 1; k <= nkin; ++k)
	  {
	    for (int npIndex = 1; npIndex <= nreactmin[k-1]; ++npIndex)
	      {
		if (get2D(rmin,k-1,npIndex-1,nkin,MAX_PATH) /* rmin[k-1][npIndex-1] */ >= 0.0)
		  {
		    // sumrct = sumrct + decay_correct[k-1][i-1] * mumin[i-1][k-1][npIndex-1] * rmin[k-1][npIndex-1];
		    sumrct = sumrct +
		      get2D(decay_correct,k-1,i-1,nkin,ncomp) *
		      get3D(mumin,i-1,k-1,npIndex-1,ncomp,nkin,MAX_PATH) *
		      get2D(rmin,k-1,npIndex-1,nkin,MAX_PATH);
		  }
		else
		  {
		    //fprintf(stderr, "assemble_local_cpp.cpp: mumin_decay not implemented in crunchflow lite\n");
		    //exit(1);
		    return;
		  }
	      }
	  }
      }
    else
      {
	sumrct = 0.0;
	for (int k = 1; k <= nkin; ++k)
	  {
	    for (int npIndex = 1; npIndex <= nreactmin[k-1]; ++npIndex)
	      {
		// sumrct = sumrct + decay_correct[k-1][i-1] * mumin[i-1][k-1][npIndex-1] * rmin[k-1][npIndex-1];
		sumrct = sumrct +
		  get2D(decay_correct,k-1,i-1,nkin,ncomp) *
		  get3D(mumin,i-1,k-1,npIndex-1,ncomp,nkin,MAX_PATH) *
		  get2D(rmin,k-1,npIndex-1,nkin,MAX_PATH);
	      }
	  }
      }
    sumkin = 0.0;
    }
}
#endif
  
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

#ifdef PROTO_CUDA

  cudaError_t cudaStat = cudaSuccess;

  size_t pitch;
  double *dA;
  
  cudaStat = cudaMalloc ((void**)&dA, sizeof(double) * neqn * neqn );
  assert(cudaSuccess == cudaStat);

  double *dSatliq;
  cudaStat = cudaMalloc ((void**)&dSatliq, sizeof(double) * nx * ny * nz );
  assert(cudaSuccess == cudaStat);

  int *dnreactmin;
  cudaStat = cudaMalloc ((void**)&dnreactmin, sizeof(int) * nkin );
  assert(cudaSuccess == cudaStat);

  cudaStat = cudaMemcpy (dnreactmin, nreactmin_1d, sizeof(int) * nkin, cudaMemcpyHostToDevice );
  assert(cudaSuccess == cudaStat);

  double *drmin;
  cudaStat = cudaMalloc ((void**)&drmin, sizeof(double) * nkin * MAX_PATH );
  assert(cudaSuccess == cudaStat);
  cudaStat = cudaMemcpy (drmin, rmin_2d, sizeof(double) * nkin * MAX_PATH, cudaMemcpyHostToDevice );
  assert(cudaSuccess == cudaStat);

  double *ddecay_correct;
  cudaStat = cudaMalloc ((void**)&ddecay_correct, sizeof(double) * nkin * ncomp );
  assert(cudaSuccess == cudaStat);
  cudaStat = cudaMemcpy (ddecay_correct, decay_correct_2d, sizeof(double) * nkin * ncomp, cudaMemcpyHostToDevice );
  assert(cudaSuccess == cudaStat);

  double *dmumin;
  cudaStat = cudaMalloc ((void**)&dmumin, sizeof(double) * ncomp * nkin * MAX_PATH );
  assert(cudaSuccess == cudaStat);
  cudaStat = cudaMemcpy (dmumin, mumin_3d, sizeof(double) * ncomp * nkin * MAX_PATH, cudaMemcpyHostToDevice );
  assert(cudaSuccess == cudaStat);
  
  assemble_local_kernel<<<1, 1>>>(dA,
		     dSatliq,
		     neqn,
		     ncomp,
		     nkin,
		     nx,
		     ny,
		     nz,
		     nradmax,
		     dt,
		     dnreactmin,
		     drmin,
		     ddecay_correct,
		     dmumin);

#else
  
  for(size_t i = 0; i < ncomp; ++i)
    for(size_t j = 0; j < ncomp; ++j)
      aaa[i][j] = 0.0;

#endif
  
  source = 0.0;
  retardation = 1.0;

  satl    = satliq[jz-1][jy-1][jx-1];
  portemp = por[jz-1][jy-1][jx-1];
  rotemp  = ro[jz-1][jy-1][jx-1];
  xgtemp = xgram[jz-1][jy-1][jx-1];

#ifdef PROTO_CUDA

  size_t pitchSatliq;
  size_t pitchPor;
  size_t pitchRo;
  size_t pitchXGram;

  double *dPor;
  double *dRo;
  double *dXGram;
  
#endif

  r = 1.0/dt;

  for (int i = 1; i <= ncomp; ++i) {
    ind = i;

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
