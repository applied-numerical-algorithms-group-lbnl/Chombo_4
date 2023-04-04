#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#ifdef MKL
#include <mkl_lapacke.h>
#else
#include <lapacke.h>
#endif
#include "ShapeArray.H"
#include "crunchflow.h"

using namespace shape;

int fitgamma(double *a,
	     const int ntemp,
	     double bvec[NBASIS],
	     double *vec_2d)  // [NTEMP][NBASIS]
{
  ShapeArray<double, 2> vec(vec_2d, NTEMP, NBASIS);
  char trans = 'N';

  double w[NBASIS][NBASIS];

  for (int j = 1; j <= NBASIS; ++j) {
    // array1D(bvec,j) = 0.0;
    bvec[j-1] = 0.0;
    for (int i = 1; i <= ntemp; ++i) {
      // bvec(j) = bvec(j) + adh(i)*vec(j,i)
      // array1D(bvec,j) = array1D(bvec,j) + array1D(adh,i)*array2D(vec,j,i,NBASIS);
      bvec[j-1] += a[i-1] * vec[i-1][j-1];
    }
  }

  for (int j = 1; j <= NBASIS; ++j) {
    for (int k = j; k <= NBASIS; ++k) {
      w[j-1][k-1] = 0.0;
      for (int i = 1; i <= ntemp; ++i) {
	// w(j,k) = w(j,k) + vec(j,i)*vec(k,i)
	//w[j-1][k-1] = w[j-1][k-1] + array2D(vec,j,i,NBASIS)*array2D(vec,k,i,NBASIS);
	w[j-1][k-1] = w[j-1][k-1] + vec[i-1][j-1] * vec[i-1][k-1];
      }
      if (j !=k) {
	w[k-1][j-1] = w[j-1][k-1];
      }
    }
  }

  lapack_int info,m,n,lda,ldb,nrhs;
  lapack_int *ipiv;
  m = NBASIS;
  n = NBASIS;
  lda = NBASIS;
  nrhs = 1;
  ldb = 1;
  ipiv = (lapack_int *) malloc(n*sizeof(lapack_int));

  info = LAPACKE_dgetrf(LAPACK_ROW_MAJOR,m,n,*w,lda,ipiv);

  if( info < 0 )
    {
      fprintf(stderr,"fitgamma: LAPACKE_dgetrf: the %d-th argument had an illegal value\n", info);
      exit(EXIT_FAILURE);
    }
  else if( info > 0 )
    {
      fprintf(stderr,"fitgamma: LAPACKE_dgetrf: U(%d,%d) is exactly zero. The factorization has been completed, but the factor U is exactly singular, and division by zero will occur if it is used to solve a system of equations.\n", info, info );
      exit(EXIT_FAILURE);
    }

  info = LAPACKE_dgetrs(LAPACK_ROW_MAJOR,trans,n,nrhs,*w,lda,ipiv,bvec,ldb);

  if( info < 0 )
    {
      fprintf(stderr,"fitgamma: LAPACKE_dgetrs: the %d-th argument had an illegal value\n", info);
      exit(EXIT_FAILURE);
    }

  return 0;
}
