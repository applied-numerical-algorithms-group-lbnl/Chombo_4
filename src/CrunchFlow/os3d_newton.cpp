#include <math.h>
#include <string.h>
#include <stdio.h>

#ifdef PROTO_CUDA
#  include <assert.h>
#  include <cuda_runtime.h>
#  include <cusolverDn.h>
#endif

#ifdef MKL
#  include <mkl_lapacke.h>
#else
#  include <lapacke.h>
#endif

#include <stdlib.h>
#include <bits/stdc++.h>
#include <chrono> 
#include "ShapeArray.H"
#include "crunchflow.h"

using namespace shape;

// SIGN(A,B) returns the value of A with the sign of B.
//
double sign(const double A, const double B)
{
  // If B >= 0 then the result is ABS(A), else it is -ABS(A).
  return( B >= 0.0 ? fabs(A) : -fabs(A) );
}

void printMatrix(int m, int n, const double*A, int lda, const char* name)
{
  for(int row = 0 ; row < m ; row++){
    for(int col = 0 ; col < n ; col++){
      double Areg = A[row + col*lda];
      printf("%s(%d,%d) = %f\n", name, row+1, col+1, Areg);
    }
  }
}

void printSolution(int n, const double* b)
{
  printf("solution: [");
  for(int col = 0 ; col < n ; col++)
    printf("%10.2e ", b[col]);
  printf("]\n");
}

void print_matrix(const int n,
		  /* const */ double *A_2d,
		  /* const */ double *b_1d,
		  const char *name)
{
  ShapeArray<double, 2> A(A_2d, n, n);
  ShapeArray<double, 1> b(b_1d, n);
  
  fprintf(stderr,"%s: C os3d_newton [A|b]: dim: %d\n", name, n);
  for(size_t i = 0; i < n; ++i)
    {
      for(size_t j = 0; j < n; ++j)
	fprintf(stderr,"%.2e\t", A[i][j]);
      fprintf(stderr,"| %.2e\n", b[i]);
    }
}

int os3d_newton(enum Target target,
		const int ncomp,
		const int nspec,
		const int nkin,
		const int nrct,
		const int ngas,
		const int ikin, 
		const int nexchange,
		const int nexch_sec,
		const int nsurf,
		const int nsurf_sec,
		const int npot,
		const int ndecay,
		const int neqn, 
		const int igamma,
		const double delt,
		const double corrmax,
		const int jx, 
		const int jy, 
		const int jz,
		int *iterat,
		int *icvg, 
		const int nx, 
		const int ny, 
		const int nz,
		const double AqueousToBulk,
		const int ntemp,
		double *sp10_4d, // [nz][ny][nx][ncomp + nspec],
		/* const */ double *t_3d, // [nz][ny][nx], 
		/* const */ double *adh_1d, // [TPOINTS]
		/* const */ double *bdh_1d, // [TPOINTS]
		/* const */ double *bdot_1d, // [TPOINTS]
		double *adhcoeff_1d, // [NBASIS],
		double *bdhcoeff_1d, // [NBASIS],
		double *bdtcoeff_1d, // [NBASIS],
		double *sion_3d, // [nz][ny][nx],
		char **ulab_1d, // [MAX_DEPEND],
		/* const */ double *acmp_1d, // [ncomp + nspec],
		/* const */ double *chg_1d, // [ncomp + nspec],
		double *gam_4d, // [nz][ny][nx][ncomp + nspec],
		double *muaq_2d, // [ncomp][nspec],
		double *sp_4d, // [nz][ny][nx][ncomp + nspec],
		double *keqaq_4d, // [nz][ny][nx][nspec],
                double *s_4d, // [nz][ny][nx][neqn],
		double *snorm_2d, // [nkin][MAX_PATH],
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
		/* const */ int *idepend_3d, // [nkin][np][nd],
		double *depend_3d, // [nkin][MAX_PATH][MAX_DEPEND],
		/* const */ int *itot_min_3d, // [nkin][np][nd],
		double *AffinityDepend1_2d, // [nkin][MAX_PATH],
		double *rate0_2d, // [nkin][MAX_PATH],
		double *volmol_1d, // [nkin]
		double *mumin_3d, // [ncomp][nkin][MAX_PATH],
		/* const */ double *keqmin_5d,		
		double *sppTMP_1d, // [ncomp + nspec]
		double *jac_rmin_3d, // [nkin][MAX_PATH][ncomp],
		/* const */ int *ivolume_1d, // [nkin],
		double *jac_sat_1d, // [ncomp],
		double *fjac_5d, // [nz][ny][nx][neqn][neqn],
		double *por_3d, // [nz][ny][nx],
		double *ro_3d, // [nz][ny][nx],
		/* const */ double *H2Oreacted_3d, // [nz][ny][nx],
		/* const */ double *xgram_3d, // [nz+2][ny+2][nx+3],
		/* const */ double *xgramOld_3d, // [nz+2][ny+2][nx+3],
		/* const */ double *sn_4d, // [nz][ny][nx][ncomp],
		double *distrib_1d, // [ncomp],
		double *fxx_1d, // [neqn],
		double *fxmax_1d, // [neqn],
		double *satliq_3d) // [nz][ny][nx]
{
  ShapeArray<double, 4> sp10(sp10_4d, nz, ny, nx, ncomp + nspec);
  ShapeArray<double, 3> t(t_3d, nz, ny, nx);
  ShapeArray<double, 1> adhcoeff(adhcoeff_1d, NBASIS);
  ShapeArray<double, 1> bdhcoeff(bdhcoeff_1d, NBASIS);
  ShapeArray<double, 1> bdtcoeff(bdtcoeff_1d, NBASIS);
  ShapeArray<double, 3> sion(sion_3d, nz, ny, nx);
  ShapeArray<char *, 1> ulab(ulab_1d, ncomp + nspec);
  ShapeArray<double, 1> acmp(acmp_1d, ncomp + nspec);
  ShapeArray<double, 1> chg(chg_1d, ncomp + nspec);
  ShapeArray<double, 4> gam(gam_4d, nz, ny, nx, ncomp + nspec);
  ShapeArray<double, 2> muaq(muaq_2d, ncomp, nspec);
  ShapeArray<double, 4> sp(sp_4d, nz, ny, nx, ncomp + nspec);
  ShapeArray<double, 4> keqaq(keqaq_4d, nz, ny, nx, nspec);
  ShapeArray<double, 4> s(s_4d, nz, ny, nx, neqn);
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
  ShapeArray<int, 2> ndepend(ndepend_2d, nkin, MAX_PATH);
  ShapeArray<int, 3> idepend(idepend_3d, nkin, MAX_PATH, MAX_DEPEND);
  ShapeArray<double, 3> depend(depend_3d, nkin, MAX_PATH, MAX_DEPEND);
  ShapeArray<int, 3> itot_min(itot_min_3d, nkin, MAX_PATH, MAX_DEPEND);
  ShapeArray<double, 2> AffinityDepend1(AffinityDepend1_2d, nkin, MAX_PATH);
  ShapeArray<double, 2> rate0(rate0_2d, nkin, MAX_PATH);
  ShapeArray<double, 1> volmol(volmol_1d, nkin);
  ShapeArray<double, 3> mumin(mumin_3d, ncomp, nkin, MAX_PATH);
  ShapeArray<double, 5> keqmin(keqmin_5d, nz, ny, nx, nkin, MAX_PATH);
  ShapeArray<double, 1> sppTMP(sppTMP_1d, ncomp + nspec);
  ShapeArray<double, 3> jac_rmin(jac_rmin_3d, nkin, MAX_PATH, ncomp);
  ShapeArray<int, 1> ivolume(ivolume_1d, nkin);
  ShapeArray<double, 1> jac_sat(jac_sat_1d, ncomp);				 
  ShapeArray<double, 5> fjac(fjac_5d, nz, ny, nx, neqn, neqn);
  ShapeArray<double, 3> por(por_3d, nz, ny, nx);
  ShapeArray<double, 3> ro(ro_3d, nz, ny, nx);
  ShapeArray<double, 3> H2Oreacted(H2Oreacted_3d, nz, ny, nx);
  ShapeArray<double, 3> xgram(xgram_3d, nz+2, ny+2, nx+3);
  ShapeArray<double, 3> xgramOld(xgramOld_3d, nz+2, ny+2, nx+3);
  ShapeArray<double, 4> sn(sn_4d, nz, ny, nx, ncomp);
  ShapeArray<double, 1> distrib(distrib_1d, ncomp);
  ShapeArray<double, 1> fxx(fxx_1d, neqn);
  ShapeArray<double, 1> fxmax(fxmax_1d, neqn);
  ShapeArray<double, 3> satliq(satliq_3d, nz, ny, nx);

#ifdef PROTO_CUDA
  cusolverDnHandle_t cusolverH = NULL;
  cudaStream_t stream = NULL;
  cusolverStatus_t status = CUSOLVER_STATUS_SUCCESS;
  cudaError_t cudaStat1 = cudaSuccess;
  cudaError_t cudaStat2 = cudaSuccess;
  cudaError_t cudaStat3 = cudaSuccess;
  cudaError_t cudaStat4 = cudaSuccess;
#endif
  
  char trans = 'N';
  const int newton = 1e3;
  const double atol = 1e-09;
  const double rtol = 1e-06;
  double aascale;
  double errmax;
  double tolmax;
  int ind;
  int nradmax = 0; // encapsulated here from initialize
  size_t i;
  double aaa[ncomp][ncomp]; // allocate(aaa(ncomp,ncomp),stat=ierr)
  double bb[ncomp]; // RHS, allocate(bb(ncomp),stat=ierr)
  double fjac_loc[neqn][neqn]; // allocate(fjac_loc(neqn,neqn),stat=ierr);, zeroed in jac_local()
  int ikh2o = 0; // species number for H2O
  double (* raq_tot)[ny][nx][IKIN] = (double (*)[ny][nx][IKIN])calloc(nz * ny * nx * IKIN, sizeof(double)); // non-contributory, consider removal
  double (* mukin)[IKIN] = (double (*)[IKIN])calloc(ncomp * IKIN, sizeof(double)); // non-contributory, consider removal
  double (* rdkin)[IKIN] = (double (*)[IKIN])calloc(ncomp * IKIN, sizeof(double)); // non-contributory, consider removal
  double (* decay_correct)[ncomp] = (double (*)[ncomp])calloc(nkin * ncomp, sizeof(double)); // non-contributory, consider removal
  double total_time_taken = 0;
  
  // nkin x ncomp
  //
  for(size_t i = 0; i < nkin; ++i)
    for(size_t j = 0; j < ncomp; ++j)
      decay_correct[i][j] = 1.0;

  if(raq_tot == NULL)
    {
      fprintf(stderr,"os3d_newton: unable to allocate dynamic memory for raq_tot, exiting.\n");
      exit(1);
    }

  if( DEBUG )
    fprintf(stderr,"os3d_newton_cpp.cpp: starting with dt = %.4e, point (%d,%d,%d)\n", delt,jx,jy,jz);

  *icvg = 0; // INTENT(OUT)
  *iterat = 0; // INTENT(OUT)

  for(unsigned int ne = 1; ne <= newton; ++ne)
    {
      *iterat = *iterat + 1;

      if(igamma == 2) 
	{
	  cfgamma(ncomp,
		  nspec,
		  ntemp,
		  jx,
		  jy,
		  jz,
		  nx,
		  ny,
		  nz,
		  sp10_4d,
		  t_3d,
		  adh_1d,
		  bdh_1d,
		  bdot_1d,
		  adhcoeff_1d,
		  bdhcoeff_1d,
		  bdtcoeff_1d,
		  sion_3d,
		  ulab_1d,
		  acmp_1d,
		  chg_1d,
		  gam_4d);
	}
  
      specieslocal(ncomp,
		   nspec,
		   jx,
		   jy,
		   jz,
		   nx,
		   ny,
		   nz,
		   muaq_2d,
		   sp_4d,
		   gam_4d,
		   sp10_4d,
		   keqaq_4d);
  
      totconc(ncomp,
	      nspec,
	      neqn,
	      jx,
	      jy,
	      jz,
	      nx,
	      ny,
	      nz,
	      muaq_2d,
	      sp10_4d,
	      s_4d);

      jac_local(ncomp,
		nspec,
		neqn,
		jx,
		jy,
		jz,
		&fjac_loc[0][0],
		nx,
		ny,
		nz,
		muaq_2d,
		sp10_4d);
  
      reaction(ncomp,
	       neqn,
	       nkin,
	       nspec,
	       jx,
	       jy,
	       jz,
	       nx,
	       ny,
	       nz,
	       snorm_2d,
	       MAX_PATH,
	       dppt_4d,
	       u_rate_4d,
	       nreactmin_1d,
	       actenergy_2d,
	       silog_2d,
	       siln_2d,
	       si_2d,
	       surf_1d,
	       area_4d,
	       rmin_2d,
	       pre_rmin_2d,
	       ndepend_2d,
	       MAX_DEPEND,
	       idepend_3d,
	       depend_3d,
	       itot_min_3d,
	       s_4d,
	       gam_4d,
	       sp_4d,
	       AffinityDepend1_2d,
	       rate0_2d,
	       volmol_1d,
	       mumin_3d,
	       keqmin_5d);
       
      jacmin(ncomp,
	     nspec,
	     nkin,
	     neqn,
	     MAX_DEPEND,
	     MAX_PATH,
	     jx,
	     jy,
	     jz,
	     nx,
	     ny,
	     nz,
	     surf_1d,
	     fjac_5d,
	     sppTMP_1d,
	     sp_4d,
	     jac_rmin_3d,
	     ivolume_1d,
	     jac_sat_1d,
	     nreactmin_1d,
	     si_2d,
	     AffinityDepend1_2d,
	     snorm_2d,
	     mumin_3d,
	     ndepend_2d,
	     idepend_3d,
	     itot_min_3d,
	     depend_3d,
	     pre_rmin_2d,
	     muaq_2d,
	     rate0_2d,
	     gam_4d,
	     keqmin_5d,
	     s_4d);
       
      fx_local(ncomp,
	       neqn,
	       delt,
	       jx,
	       jy,
	       jz,
	       nx,
	       ny,
	       nz,
	       ikh2o,
	       por_3d,
	       ro_3d,
	       H2Oreacted_3d,
	       xgram_3d,
	       s_4d,
	       xgramOld_3d,
	       sn_4d,
	       distrib_1d,
	       fxx_1d,
	       fxmax_1d,
	       satliq_3d);
      
      assemble_local(target,
		     ncomp, 
		     nspec, 
		     nkin, 
		     ikin,
		     neqn, 
		     delt,
		     nx,
		     ny,
		     nz,
		     jx,
		     jy,
		     jz,
		     satliq_3d, 
		     por_3d,
		     ro_3d, 
		     rmin_2d, 
		     &decay_correct[0][0], 
		     &aaa[0][0],
		     nradmax, 
		     nreactmin_1d, 
		     MAX_PATH,
		     fxx_1d, 
		     mumin_3d, 
		     &mukin[0][0], 
		     &raq_tot[0][0][0][0],
		     jac_rmin_3d, 
		     &rdkin[0][0], 
		     ikh2o, 
		     distrib_1d,
		     H2Oreacted_3d, 
		     &fjac_loc[0][0], 
		     xgram_3d);

      for(i = 1; i <= neqn; i++)
	{
	  aascale = 1.0;
	  for(size_t i2 = 1; i2 <= neqn; i2++)
	    aaa[i2-1][i-1] = aaa[i2-1][i-1]/aascale;
	  bb[i-1] = -fxx[i-1]/aascale;
	}

      if( DEBUG )
	{
	  fprintf(stderr,"iteration %d, solve system:\n",*iterat);
	  print_matrix(neqn, &aaa[0][0], bb, "host");
	}

#ifdef PROTO_CUDA

      if(target == DEVICE)
	{

      int info = 0;       /* host copy of error info */
      int m = neqn, lda = ncomp, ldb = ncomp, nrhs = 1;
      double *d_A = NULL; /* device copy of A */
      double *d_B = NULL; /* device copy of B */
      int *d_Ipiv = NULL; /* pivoting sequence */
      int *d_info = NULL; /* error info */
      int  lwork = 0;     /* size of workspace */
      double *d_work = NULL; /* device workspace for getrf */
      double A[ncomp][ncomp]; // Column Major Order 
      const int pivot_on = 1;

      if( DEBUG )
	if( pivot_on && DEBUG )
	  {
	    printf("pivot is on : compute P*A = L*U \n");
	  }
	else
	  {
	    printf("pivot is off: compute A = L*U (not numerically stable)\n");
	  }

      // A = aaa^T, i.e. copy Row Major Order aaa to Column Major Order A
      //
      for(size_t row = 0 ; row < ncomp; row++)
	for(size_t col = 0 ; col < ncomp; col++)
	  A[col][row] = aaa[row][col];

      if( DEBUG )
	{
	  printf("A = (matlab base-1)\n");
	  printMatrix(m, m, &A[0][0], lda, "A");
	  printf("=====\n");
	  
	  printf("B = (matlab base-1)\n");
	  printMatrix(m, 1, bb, ldb, "B");
	  printf("=====\n");
	}
      
      /* step 1: create cusolver handle, bind a stream */
      status = cusolverDnCreate(&cusolverH);
      assert(CUSOLVER_STATUS_SUCCESS == status);

      cudaStat1 = cudaStreamCreateWithFlags(&stream, cudaStreamNonBlocking);
      assert(cudaSuccess == cudaStat1);

      status = cusolverDnSetStream(cusolverH, stream);
      assert(CUSOLVER_STATUS_SUCCESS == status);

      /* step 2: copy A to device */
      cudaStat1 = cudaMalloc ((void**)&d_A, sizeof(double) * lda * m);
      cudaStat2 = cudaMalloc ((void**)&d_B, sizeof(double) * m);
      cudaStat2 = cudaMalloc ((void**)&d_Ipiv, sizeof(int) * m);
      cudaStat4 = cudaMalloc ((void**)&d_info, sizeof(int));
      assert(cudaSuccess == cudaStat1);
      assert(cudaSuccess == cudaStat2);
      assert(cudaSuccess == cudaStat3);
      assert(cudaSuccess == cudaStat4);

      cudaStat1 = cudaMemcpy(d_A, A, sizeof(double) * lda * m, cudaMemcpyHostToDevice);
      cudaStat2 = cudaMemcpy(d_B, bb, sizeof(double) * m, cudaMemcpyHostToDevice);
      assert(cudaSuccess == cudaStat1);
      assert(cudaSuccess == cudaStat2);

      /* step 3: query working space of getrf */
      status = cusolverDnDgetrf_bufferSize(cusolverH,
					   m,
					   m,
					   d_A,
					   lda,
					   &lwork);
      assert(CUSOLVER_STATUS_SUCCESS == status);

      cudaStat1 = cudaMalloc((void**)&d_work, sizeof(double)*lwork);
      assert(cudaSuccess == cudaStat1);

      /* step 4: LU factorization */
      if (pivot_on)
	{
	  status = cusolverDnDgetrf(cusolverH,
				    m,
				    m,
				    d_A,
				    lda,
				    d_work,
				    d_Ipiv,
				    d_info);
	}
      else
	{
	  status = cusolverDnDgetrf(cusolverH,
				    m,
				    m,
				    d_A,
				    lda,
				    d_work,
				    NULL,
				    d_info);
	}
      cudaStat1 = cudaDeviceSynchronize();
      assert(CUSOLVER_STATUS_SUCCESS == status);
      assert(cudaSuccess == cudaStat1);

      double LU[lda * m]; /* L and U */
      int Ipiv[m]; /* host copy of pivoting sequence */
      
      if (pivot_on)
	{
	  cudaStat1 = cudaMemcpy(Ipiv, d_Ipiv, sizeof(int) * m, cudaMemcpyDeviceToHost);
	}
      cudaStat2 = cudaMemcpy(LU, d_A, sizeof(double) * lda * m, cudaMemcpyDeviceToHost);
      cudaStat3 = cudaMemcpy(&info, d_info, sizeof(int), cudaMemcpyDeviceToHost);
      assert(cudaSuccess == cudaStat1);
      assert(cudaSuccess == cudaStat2);
      assert(cudaSuccess == cudaStat3);
      
      if ( 0 > info )
	{
	  printf("%d-th parameter is wrong \n", -info);
	  exit(1);
	}

      if( DEBUG )
	{
	  if (pivot_on)
	    {
	      printf("pivoting sequence, matlab base-1\n");
	      for(int j = 0 ; j < m ; j++)
		printf("Ipiv(%d) = %d\n", j+1, Ipiv[j]);
	    }
	  printf("L and U = (matlab base-1)\n");
	  printMatrix(m, m, LU, lda, "LU");
	  printf("=====\n");
	}
      
      /* step 5: solve A*X = B */ 
      if (pivot_on)
	{
	  status = cusolverDnDgetrs(cusolverH,
				    CUBLAS_OP_N,
				    m,
				    nrhs,
				    d_A,
				    lda,
				    d_Ipiv,
				    d_B,
				    ldb,
				    d_info);
	}
      else
	{
	  status = cusolverDnDgetrs(cusolverH,
				    CUBLAS_OP_N,
				    m,
				    nrhs,
				    d_A,
				    lda,
				    NULL,
				    d_B,
				    ldb,
				    d_info);
	}

      switch(status)
	{
	case CUSOLVER_STATUS_SUCCESS:
	  if( DEBUG )
	    printf("cusolverDnDgetrs completed successfully.\n");
	  break;
	case CUSOLVER_STATUS_NOT_INITIALIZED:
	  printf("cusolverDnDgetrs: error: the library was not initialized.\n");
	  break;
	case CUSOLVER_STATUS_INVALID_VALUE:
	  printf("cusolverDnDgetrs: error: invalid parameters were passed (n<0 or lda<max(1,n) or ldb<max(1,n)).\n");
	  break;
	case CUSOLVER_STATUS_ARCH_MISMATCH:
	  printf("cusolverDnDgetrs: error: the device only supports compute capability 2.0 and above.\n");
	  break;
	case CUSOLVER_STATUS_INTERNAL_ERROR:
	  printf("cusolverDnDgetrs: error: an internal operation failed.\n");
	  break;
	default:
	  printf("cusolverDnDgetrs: unknown error.\n");
	};

      cudaStat1 = cudaDeviceSynchronize();
      assert(CUSOLVER_STATUS_SUCCESS == status);
      
      assert(cudaSuccess == cudaStat1);

      cudaStat1 = cudaMemcpy(bb, d_B, sizeof(double) * m, cudaMemcpyDeviceToHost);
      assert(cudaSuccess == cudaStat1);

      // free resources
      if (d_A    ) cudaFree(d_A);
      if (d_B    ) cudaFree(d_B);
      if (d_Ipiv ) cudaFree(d_Ipiv);
      if (d_info ) cudaFree(d_info);
      if (d_work ) cudaFree(d_work);

      if (cusolverH   ) cusolverDnDestroy(cusolverH);
      if (stream      ) cudaStreamDestroy(stream);

	}
#endif

      if( target == HOST)
	{
      auto start = std::chrono::high_resolution_clock::now();
      
      //      CALL dgetrf(neqn,neqn,aaa,neqn,indd,info)
      //
      lapack_int info, m = neqn, n = neqn, lda = neqn, ldb = 1, nrhs = 1;
      lapack_int *ipiv = (lapack_int *) malloc(n * sizeof(lapack_int));
      
      //                    1                2 3 4          5   6
      info = LAPACKE_dgetrf(LAPACK_ROW_MAJOR,m,n,&aaa[0][0],lda,ipiv);

      if( info < 0 )
	{
	  fprintf(stderr,"os3d_newton_cpp.cpp: LAPACKE_dgetrf: the %d-th argument had an illegal value\n", info);
	  exit(1);
	}
      else if( info > 0 )
	{
	  fprintf(stderr,"os3d_newton.cpp: LAPACKE_dgetrf: U(%d,%d) is exactly zero. The factorization has been completed, but the factor U is exactly singular, and division by zero will occur if it is used to solve a system of equations.\n", info, info );
	  exit(1);
	}

      //      CALL dgetrs(trans,neqn,ione,aaa,neqn,indd,bb,neqn,info)

      info = LAPACKE_dgetrs(LAPACK_ROW_MAJOR,trans,n,nrhs,&aaa[0][0],lda,ipiv,bb,ldb);

      if( info < 0 )
	{
	  fprintf(stderr,"os3d_newton.cpp: LAPACKE_dgetrs: the %d-th argument had an illegal value\n", info);
	  exit(1);
	}

      auto end = std::chrono::high_resolution_clock::now();
      double time_taken =
	std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();
      total_time_taken += time_taken; 
	}

      if( DEBUG )
      	{
	  printSolution(ncomp, bb);
  	}
      
      errmax = 0.0;
      for(i = 1; i <= ncomp; i++)
	{
	  if (fabs(bb[i-1]) > errmax)
	    errmax = fabs(bb[i-1]);
	  if (fabs(bb[i-1]) > corrmax)
	    bb[i-1] = sign(corrmax,bb[i-1]);

	  sp[jz-1][jy-1][jx-1][i-1] += bb[i-1];
	  sp10[jz-1][jy-1][jx-1][i-1] = exp(sp[jz-1][jy-1][jx-1][i-1]);
	}

      // fprintf(stderr,"os3d_newton_cpp.cpp: solution: %2d: ", *iterat);
      // for(i = 1; i <= ncomp; i++)
      //     fprintf(stderr,"%12.4e  ", bb[i-1]);
      // fprintf(stderr,"\n");

      if (fabs(errmax) < 1e-15)
	{
	  *icvg = 0;
	  if( DEBUG )
	    printf("os3d_newton_cpp.cpp: convergence: fabs(errmax) < 1e-15\n");
	  break;
	}

      *icvg = 0;

      for(i = 1; i <= ncomp; i++)
	{
	  ind = i;
	  tolmax = atol + rtol * sp10[jz-1][jy-1][jx-1][i-1];
	  if (fabs(delt * fxx[ind-1]) > tolmax)
	    {
	      *icvg = 1;
	      if( DEBUG )
		{
		  fprintf(stderr,"os3d_newton_cpp.cpp: fabs(delt * fxx[ind-1]) > tolmax\n");
		  fprintf(stderr,"os3d_newton_cpp.cpp: delt = %e, fxx[ind-1] = %e, tolmax = %e, ind = %d, abs = %e\n",
			  delt, fxx[ind-1], tolmax,ind,fabs(delt * fxx[ind-1]));
		}
	      break;
	    }
	}

      if (*icvg == 0)
	{
	  if( DEBUG )
	    fprintf(stderr,"os3d_newton_cpp.cpp: icvg == 0\n");
	  break;
	}
    } // for(unsigned int ne = 1; ne <= newton; ++ne)

  if (*icvg == 1)
    {
      printf(" C: No convergence at local Chombo coordinate: (%d,%d,%d)", jx-1, jy-1, jz-1);
      printf(" sp10: %e", sp10[jz-1][jy-1][jx-1][i-1]);
      printf(" sp: %e", sp[jz-1][jy-1][jx-1][i-1]);
      printf(" por: %e", por[jz-1][jy-1][jx-1]);
      printf(" area: %e\n", area[jz-1][jy-1][jx-1][0]);
    }
  else
    {
      if( DEBUG )
	printf(" C: Convergence at local Chombo coordinate: (%d,%d,%d), dt = %.2e\n", jx, jy, jz, delt);
    }

  if( DEBUG )
    fprintf(stderr,"os3d_newton_cpp.cpp: completed with dt = %.4e, point (%d,%d,%d)\n", delt,jx,jy,jz);

  free(raq_tot);
  free(mukin);
  free(rdkin);
  free(decay_correct);

  if(target == HOST)
    fprintf(stderr,"os3d_newton: LU host solver: %.0f ns\n", total_time_taken);

  return *icvg;
}
