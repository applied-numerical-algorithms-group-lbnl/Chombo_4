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
#include <assert.h>

#include "ShapeArray.H"
#include "crunchflow.h"

using namespace shape;

void print(double* A, const int n, const char *s)
{
  fprintf(stderr,"%s: A[%d][%d]\n", s, n, n);
  for(size_t i = 0; i < n; i++)
    {
      for(size_t j = 0; j < n; j++)
	fprintf(stderr,"%.2e\t", get2D(A,i,j,n,n));
      fprintf(stderr,"\n");
    }
}

#ifdef PROTO_CUDA

__global__ void assemble_local_kernel(double* aaa, // [neqn][neqn]
				      const int neqn,
				      const int ncomp,
				      const int nkin,
				      const int ikin,
				      const int nx,
				      const int ny,
				      const int nz,
				      const int jx, 
				      const int jy, 
				      const int jz, 
				      const int nradmax,
				      const int ikh2o,
				      const double dt,
				      int *nreactmin, // [nkin]
				      double *rmin, // [nkin][MAX_PATH]
				      double *decay_correct, // [nkin][ncomp]
				      double *mumin, // [ncomp][nkin][MAX_PATH]
				      double *raq_tot, // [nz][ny][nx][ikin]
				      double *mukin, // [ncomp][ikin]
				      double *fxx, // [neqn]
				      double *sumrd,  // [neqn]
				      double *sumjackin, // [neqn]
				      double *jac_rmin, // [nkin][MAX_PATH][ncomp]
				      double *rdkin, // [ncomp][ikin]
				      double *H2Oreacted, // [nz][ny][nx]
				      double *fjac_loc, // [neqn][neqn]
				      double *distrib, // [ncomp]
				      double* satliq, // [nz][ny][nx]
				      double* por, // [nz][ny][nx]
				      double* ro, // [nz][ny][nx]
				      double *xgram) // [nz+2][ny+2][nx+3] 
{
  double r = 1.0/dt;
  double sumrct;
  double sumkin;
  double satl = get3D(satliq,jz-1,jy-1,jx-1,nz,ny,nx);
  double portemp = get3D(por,jz-1,jy-1,jx-1,nz,ny,nx);
  double rotemp = get3D(ro,jz-1,jy-1,jx-1,nz,ny,nx);
  double xgtemp = get3D(xgram,jz-1,jy-1,jx-1,nz+2,ny+2,nx+3);
  double rxnmin;
  double rxnaq;
  double aq_accum;
  double retardation = 1.0;
  double source_jac;
  double source = 0.0;
  double ex_accum;
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

    for (int ir = 1; ir <= ikin; ++ir)
      {
	// sumkin = sumkin - mukin[i-1][ir-1] * raq_tot[jz-1][jy-1][jx-1][ir-1];
	sumkin = sumkin - get2D(mukin,i-1,ir-1,ncomp,ikin) * get4D(raq_tot,jz-1,jy-1,jx-1,ir-1,nz,ny,nx,ikin);
      }

    // Update the residual, adding reaction terms and exchange terms
    fxx[ind-1] += sumrct + satl * xgtemp * portemp * rotemp * sumkin;

    for(size_t j = 0; j < neqn; ++j)
      {
	sumrd[j] = 0.0;
	sumjackin[j] = 0.0;
      }
    if (nradmax > 0)
      {
	for (int k = 1; k <= nkin; ++k)
	  {
	    for (int npIndex = 1; npIndex <= nreactmin[k-1]; ++npIndex)
	      {
		if( get3D(mumin,i-1,k-1,npIndex-1,ncomp,nkin,MAX_PATH) != 0.0 )
		  {
		    if( get2D(rmin,k-1,npIndex-1,nkin,MAX_PATH) >= 0.0 )
		      {
			for (int i2 = 1; i2 <= ncomp; ++i2) 
			  {
			    sumrd[i2-1] = sumrd[i2-1] + 
			      get2D(decay_correct,k-1,i-1,nkin,ncomp) *
			      get3D(mumin,i-1,k-1,npIndex-1,ncomp,nkin,MAX_PATH) *
			      get3D(jac_rmin,k-1,npIndex-1,i2-1,nkin,MAX_PATH,ncomp);
			  }
		      }
		    else
		      {
			// mumin_decay not implemented in crunchflow lite
			return;
		      }
		  }
	      }
	  }
      }
    else
      {
	for (int k = 1; k <= nkin; ++k)
	  {
	    for (int npIndex = 1; npIndex <= nreactmin[k-1]; ++npIndex)
	      {
		if ( get3D(mumin,i-1,k-1,npIndex-1,ncomp,nkin,MAX_PATH) != 0.0 )
		  {
		    for (int i2 = 1; i2 <= ncomp; ++i2)
		      {
			sumrd[i2-1] = sumrd[i2-1] + 
			  get2D(decay_correct,k-1,i-1,nkin,ncomp) *
			  get3D(mumin,i-1,k-1,npIndex-1,ncomp,nkin,MAX_PATH) *
			  get3D(jac_rmin,k-1,npIndex-1,i2-1,nkin,MAX_PATH,ncomp);
		      }
		  }
	      }
	  }
      }
    for (int ir = 1; ir <= ikin; ++ir)
      {
	if ( get2D(mukin,i-1,ir-1,ncomp,ikin) != 0.0 )
	  {
	    for (int i2 = 1; i2 <= ncomp; ++i2)
	      {
		sumjackin[i2-1] = sumjackin[i2-1] -
		  get2D(mukin,i-1,ir-1,ncomp,ikin) *
		  get2D(rdkin,i2-1,ir-1,ncomp,ikin);
	      }
	  }
      }

    for (int i2 = 1; i2 <= ncomp; ++i2)
      {
	rxnmin = sumrd[i2-1];
	rxnaq = satl * xgtemp * portemp * rotemp * sumjackin[i2-1];

	if (i != ikh2o)
	  {
	    aq_accum = get3D(H2Oreacted,jz-1,jy-1,jx-1,nz,ny,nx) *
	      satl * xgtemp * r * portemp * rotemp *
	      get2D(fjac_loc,i-1,i2-1,neqn,neqn) * 
	      (1.0 + retardation * distrib[i-1]);	      
	  }
	else
	  {
	    aq_accum = satl * xgtemp * r * portemp * rotemp *
	      get2D(fjac_loc,i-1,i2-1,neqn,neqn) *
	      (1.0 + retardation * distrib[i-1]);
	  }
	source_jac = source * get2D(fjac_loc,i-1,i2-1,neqn,neqn);
	ex_accum = 0.0;
	get2D(aaa,i-1,i2-1,neqn,neqn) = rxnmin + rxnaq + aq_accum - source_jac + ex_accum;
      }
    }
}
#endif
  
void assemble_local(enum Target target,
		    const int ncomp, 
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
		    /* const */ double *distrib, // [ncomp]
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

  double source, retardation, satl, portemp, rotemp, xgtemp, sumrct, r, sumkin,
    rxnmin, rxnaq, aq_accum,ex_accum,source_jac;
  int ind;
  double sumrd[neqn]; // allocate(sumrd(neqn),stat=ierr); sumrd = 0.0d0
  double sumjackin[neqn]; // allocate(sumjackin(neqn),stat=ierr); sumjackin = 0.0d0

#ifdef PROTO_CUDA

  if(target == DEVICE)
    {
  cudaError_t cudaStat = cudaSuccess;

  double *daaa;
  cudaStat = cudaMalloc ((void**)&daaa, sizeof(double) * neqn * neqn );
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

  double *draq_tot;
  cudaStat = cudaMalloc ((void**)&draq_tot, sizeof(double) * nz * ny * nx * ikin );
  assert(cudaSuccess == cudaStat);
  cudaStat = cudaMemcpy (draq_tot, raq_tot_4d, sizeof(double) * nz * ny * nx * ikin, cudaMemcpyHostToDevice );
  assert(cudaSuccess == cudaStat);

  double *dmukin;
  cudaStat = cudaMalloc ((void**)&dmukin, sizeof(double) * ncomp * ikin );
  assert(cudaSuccess == cudaStat);
  cudaStat = cudaMemcpy (dmukin, mukin_2d, sizeof(double) * ncomp * ikin, cudaMemcpyHostToDevice );
  assert(cudaSuccess == cudaStat);

  double *dfxx;
  cudaStat = cudaMalloc ((void**)&dfxx, sizeof(double) * neqn );
  assert(cudaSuccess == cudaStat);
  cudaStat = cudaMemcpy (dfxx, fxx_1d, sizeof(double) * neqn, cudaMemcpyHostToDevice );
  assert(cudaSuccess == cudaStat);

  double *dsumrd;
  cudaStat = cudaMalloc ((void**)&dsumrd, sizeof(double) * neqn );
  assert(cudaSuccess == cudaStat);
  cudaStat = cudaMemcpy (dsumrd, sumrd, sizeof(double) * neqn, cudaMemcpyHostToDevice );
  assert(cudaSuccess == cudaStat);

  double *dsumjackin;
  cudaStat = cudaMalloc ((void**)&dsumjackin, sizeof(double) * neqn );
  assert(cudaSuccess == cudaStat);
  cudaStat = cudaMemcpy (dsumjackin, sumjackin, sizeof(double) * neqn, cudaMemcpyHostToDevice );
  assert(cudaSuccess == cudaStat);

  double *djac_rmin;
  cudaStat = cudaMalloc ((void**)&djac_rmin, sizeof(double) * nkin * MAX_PATH * ncomp );
  assert(cudaSuccess == cudaStat);
  cudaStat = cudaMemcpy (djac_rmin, jac_rmin_3d, sizeof(double) * nkin * MAX_PATH * ncomp, cudaMemcpyHostToDevice );
  assert(cudaSuccess == cudaStat);

  double *drdkin;
  cudaStat = cudaMalloc ((void**)&drdkin, sizeof(double) * ncomp * ikin );
  assert(cudaSuccess == cudaStat);
  cudaStat = cudaMemcpy (drdkin, rdkin_2d, sizeof(double) * ncomp * ikin, cudaMemcpyHostToDevice );
  assert(cudaSuccess == cudaStat);

  double *dH2Oreacted;
  cudaStat = cudaMalloc ((void**)&dH2Oreacted, sizeof(double) * nz * ny * nx );
  assert(cudaSuccess == cudaStat);
  cudaStat = cudaMemcpy (dH2Oreacted, H2Oreacted_3d, sizeof(double) * nz * ny * nx, cudaMemcpyHostToDevice );
  assert(cudaSuccess == cudaStat);

  double *dfjac_loc;
  cudaStat = cudaMalloc ((void**)&dfjac_loc, sizeof(double) * neqn * neqn );
  assert(cudaSuccess == cudaStat);
  cudaStat = cudaMemcpy (dfjac_loc, fjac_loc_2d, sizeof(double) * neqn * neqn, cudaMemcpyHostToDevice );
  assert(cudaSuccess == cudaStat);
  
  double *ddistrib;
  cudaStat = cudaMalloc ((void**)&ddistrib, sizeof(double) * ncomp );
  assert(cudaSuccess == cudaStat);
  cudaStat = cudaMemcpy (ddistrib, distrib, sizeof(double) * ncomp, cudaMemcpyHostToDevice );
  assert(cudaSuccess == cudaStat);

  double *dsatliq;
  cudaStat = cudaMalloc ((void**)&dsatliq, sizeof(double) * nz * ny * nx );
  assert(cudaSuccess == cudaStat);
  cudaStat = cudaMemcpy (dsatliq, satliq_3d, sizeof(double) * nz * ny * nx, cudaMemcpyHostToDevice );
  assert(cudaSuccess == cudaStat);

  double *dpor;
  cudaStat = cudaMalloc ((void**)&dpor, sizeof(double) * nz * ny * nx );
  assert(cudaSuccess == cudaStat);
  cudaStat = cudaMemcpy (dpor, por_3d, sizeof(double) * nz * ny * nx, cudaMemcpyHostToDevice );
  assert(cudaSuccess == cudaStat);

  double *dro;
  cudaStat = cudaMalloc ((void**)&dro, sizeof(double) * nz * ny * nx );
  assert(cudaSuccess == cudaStat);
  cudaStat = cudaMemcpy (dro, ro_3d, sizeof(double) * nz * ny * nx, cudaMemcpyHostToDevice );
  assert(cudaSuccess == cudaStat);
  
  double *dxgram;
  cudaStat = cudaMalloc ((void**)&dxgram, sizeof(double) * (nz+2) * (ny+2) * (nx+3) );
  assert(cudaSuccess == cudaStat);
  cudaStat = cudaMemcpy (dxgram, xgram_3d, sizeof(double) * (nz+2) * (ny+2) * (nx+3), cudaMemcpyHostToDevice );
  assert(cudaSuccess == cudaStat);
  
  assemble_local_kernel<<<1, 1>>>(daaa,
				  neqn,
				  ncomp,
				  nkin,
				  ikin, 
				  nx,
				  ny,
				  nz,
				  jx, 
				  jy, 
				  jz, 
				  nradmax,
				  ikh2o,
				  dt,
				  dnreactmin,
				  drmin,
				  ddecay_correct,
				  dmumin,
				  draq_tot,
				  dmukin,
				  dfxx,
				  dsumrd,
				  dsumjackin,
				  djac_rmin,
				  drdkin,
				  dH2Oreacted,
				  dfjac_loc,
				  ddistrib,
				  dsatliq,
				  dpor,
				  dro,
				  dxgram);

  //  double *A = (double *)calloc(neqn * neqn, sizeof(double));
  //cudaStat = cudaMemcpy (A, daaa, sizeof(double) * neqn * neqn, cudaMemcpyDeviceToHost );
   cudaStat = cudaMemcpy (aaa_2d, daaa, sizeof(double) * neqn * neqn, cudaMemcpyDeviceToHost );
  assert(cudaSuccess == cudaStat);

  // for(size_t i = 0; i < neqn; i++)
  //   for(size_t j = 0; j < neqn; j++)
  //     aaa[i][j] = get2D(A,i,j,neqn,neqn);
  
  //  print(&aaa[0][0], neqn, "device");
  return;
    }  
#else
  
  for(size_t i = 0; i < ncomp; ++i)
    for(size_t j = 0; j < ncomp; ++j)
      aaa[i][j] = 0.0;

#endif
  
  source = 0.0;
  retardation = 1.0;

  satl = satliq[jz-1][jy-1][jx-1];
  portemp = por[jz-1][jy-1][jx-1];
  rotemp = ro[jz-1][jy-1][jx-1];
  xgtemp = xgram[jz-1][jy-1][jx-1];
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
  //  print(&aaa[0][0], neqn, "host");
}
