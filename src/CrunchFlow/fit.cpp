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

void fit(const int nbasis,
	 const int ntemp,
	 double *bvec_1d,
	 /* const */ double *vec_2d,
	 double *alogk0,
	 double *bvecF,
	 double *vecF,
	 int *iflgint,
	 int *inoutint,
	 int *ntt,
	 char *nameTransfer)
{
  ShapeArray<double, 1> bvec(bvec_1d, nbasis);
  ShapeArray<double, 2> vec(vec_2d, ntemp, nbasis);
  
  char trans = 'N';
  bool NonZeroLogK;
  double vectmp[ntemp][nbasis];
  double w[nbasis][nbasis];
  double bvectmp[nbasis];
  double det;
  int nbcalc;

  *iflgint = 0;
  *ntt     = 0;
  for (int i = 1; i <= ntemp; ++i) {
    if ( alogk0[i-1] != 500.0) {
      inoutint[i-1] = 1;
      *ntt = *ntt + 1;
    } else {
      *iflgint = 1;
      inoutint[i-1] = 0;
    }
  }

  NonZeroLogK = false;

  for (int i = 1; i <= ntemp; ++i) {
    if (inoutint[i-1] !=0) {
      NonZeroLogK = true;
      break;
    }
  }

  if (!NonZeroLogK) {
    fprintf(stderr,
	    "Row of equilibrium constants all with 500\n Species or mineral: %s\n", 
	    nameTransfer);
    exit(1);
  }
 
  memset(vectmp, 0, sizeof(double) * ntemp * nbasis);
  if (*ntt == 1){
    for (int i = 1; i <= ntemp; ++i) {
      vectmp[i-1][1-1] = vec[i-1][2-1];
    }
    nbcalc = 1;
  } else if(*ntt == 2) {
    for (int i = 1; i <= ntemp; ++i) {
      vectmp[i-1][1-1] = vec[i-1][2-1];
      vectmp[i-1][2-1] = vec[i-1][3-1];
    }
    nbcalc = 2;
  } else {
    for (int i = 0; i < nbasis; ++i) 
      for (int w = 0; w < ntemp; ++w) 
	vectmp[w][i] = vec[w][i];
   nbcalc = 5;
  }

  for (int i = 0; i < nbasis; ++i) 
    bvec[i] = 0.0;
  
  memset(bvectmp, 0, sizeof(double) * nbasis);
  
  for (int j = 1; j <= nbcalc; ++j) {
    for (int i = 1; i <= ntemp; ++i) {
      if ( inoutint[i-1] == 1) {
	bvectmp[j-1] = bvectmp[j-1] + alogk0[i-1] * vectmp[i-1][j-1];
      }
    }
  }
  memset(w, 0, sizeof(double) * nbasis * nbasis);
  for (int j = 1; j <= nbcalc; ++j) {
    for (int k = j; k <= nbcalc ; ++k) {
      for (int i = 1; i <= ntemp; ++i) {
	if ( inoutint[i-1] == 1) {
	  w[k-1][j-1] = w[k-1][j-1] + vectmp[i-1][j-1] * vectmp[i-1][k-1];
	}
      }
      if (j != k ) {
	w[j-1][k-1] = w[k-1][j-1]; // make w symmetric
      }
    }
  }

  double wscale[nbcalc][nbcalc];
  double bvectmp2[nbcalc];
    
  for(int i = 0; i < nbcalc; ++i)
    for(int j = 0; j < nbcalc; ++j)
      wscale[i][j] = w[i][j];

  memcpy(bvectmp2, bvectmp, sizeof(bvectmp[0]) * nbcalc);

  lapack_int info,m,n,lda,ldb,nrhs;
  lapack_int *ipiv;
  m = nbcalc;
  n = nbcalc;
  lda = nbcalc;
  nrhs = 1;
  ldb = 1; // the leading dimension of the array specified for B
  ipiv = (lapack_int *) malloc(n*sizeof(lapack_int));

  //                    1                2 3 4       5   6
  info = LAPACKE_dgetrf(LAPACK_ROW_MAJOR,m,n,*wscale,lda,ipiv);

  if( info < 0 )
    {
      fprintf(stderr,"fit_cpp.cpp: LAPACKE_dgetrf: the %d-th argument had an illegal value\n", info);
      exit;
    }
  else if( info > 0 )
    {
      fprintf(stderr,"fit_cpp.cpp: LAPACKE_dgetrf: U(%d,%d) is exactly zero. The factorization has been completed, but the factor U is exactly singular, and division by zero will occur if it is used to solve a system of equations.\n", info, info );
      exit;
    }

  info = LAPACKE_dgetrs(LAPACK_ROW_MAJOR,trans,n,nrhs,*wscale,lda,ipiv,bvectmp,ldb);

  if( info < 0 )
    {
      fprintf(stderr,"fit_cpp.cpp: LAPACKE_dgetrs: the %d-th argument had an illegal value\n", info);
      exit;
    }

  if (nbcalc == 1) {
    bvec[2-1] = bvectmp[0];
  } else if (nbcalc == 2){
    bvec[2-1] = bvectmp[0];
    bvec[3-1] = bvectmp[1];
  } else{
    for (int i = 1; i <= 5; ++i) {
      bvec[i-1] = bvectmp[i-1];
    }
  }
  free(ipiv);
  return;
}
