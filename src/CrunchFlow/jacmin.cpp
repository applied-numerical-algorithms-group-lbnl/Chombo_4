#include <math.h>
#include "ShapeArray.H"
#include "crunchflow.h"

#define PERTURB 1e-09

using namespace shape;

void jacmin(const int ncomp, 
	    const int nspec, 
	    const int nkin, 
	    const int neqn, 
	    const int nd,
	    const int np,  
	    const int jx, 
	    const int jy, 
	    const int jz,
	    const int nx,
	    const int ny, 
	    const int nz, 
	    /* const */ double *surf_1d, // [nkin], 
	    double *fjac_5d, // [nz][ny][nx][neqn][neqn], 
	    double *sppTMP_1d, // [ncomp + nspec]
	    double *sp_4d, // [nz][ny][nx][ncomp + nspec],
	    double *jac_rmin_3d, // [nkin][MAX_PATH][ncomp],
	    /* const */ int *ivolume_1d, // [nkin], 
	    double *jac_sat_1d, // [ncomp], 
	    /* const */ int *nreactmin_1d, // [nkin], 
	    double *si_2d, // [nkin][MAX_PATH], 
	    double *AffinityDepend1_2d, // [nkin][MAX_PATH],
	    double *snorm_2d, // [nkin][MAX_PATH], 
	    double *mumin_3d, // [ncomp][nkin][MAX_PATH],
	    /* const */ int *ndepend_2d, // [nkin][np], 
	    /* const */ int *idepend_3d, // [nkin][np][nd],
	    /* const */ int *itot_min_3d, // [nkin][np][nd],
	    double *depend_3d, // [nkin][MAX_PATH][MAX_DEPEND], 
	    /* const */ double *pre_rmin_2d, // [nkin][MAX_PATH], 
	    double *muaq_2d, // [ncomp][nspec],
	    double *rate0_2d, // [nkin][MAX_PATH], 
	    double *gam_4d, // [nz][ny][nx][ncomp + nspec],
	    /* const */ double *keqmin_5d, // [nz][ny][nx][nkin][MAX_PATH], 
	    /* const */ double *s_4d) // [nz][ny][nx][neqn])
{
  ShapeArray<double, 1> surf(surf_1d, nkin); 
  ShapeArray<double, 5> fjac(fjac_5d, nz, ny, nx, neqn, neqn); 
  ShapeArray<double, 1> sppTMP(sppTMP_1d, ncomp + nspec);
  ShapeArray<double, 4> sp(sp_4d, nz, ny, nx, ncomp + nspec);
  ShapeArray<double, 3> jac_rmin(jac_rmin_3d, nkin, MAX_PATH, ncomp);
  ShapeArray<int, 1> ivolume(ivolume_1d, nkin); 
  ShapeArray<double, 1> jac_sat(jac_sat_1d, ncomp); 
  ShapeArray<int, 1> nreactmin(nreactmin_1d, nkin); 
  ShapeArray<double, 2> si(si_2d, nkin, MAX_PATH); 
  ShapeArray<double, 2> AffinityDepend1(AffinityDepend1_2d, nkin, MAX_PATH);
  ShapeArray<double, 2> snorm(snorm_2d, nkin, MAX_PATH); 
  ShapeArray<double, 3> mumin(mumin_3d, ncomp, nkin, MAX_PATH);
  ShapeArray<int, 2> ndepend(ndepend_2d, nkin, np); 
  ShapeArray<int, 3> idepend(idepend_3d, nkin, np, nd);
  ShapeArray<int, 3> itot_min(itot_min_3d, nkin, np, nd);
  ShapeArray<double, 3> depend(depend_3d, nkin, MAX_PATH, MAX_DEPEND); 
  ShapeArray<double, 2> pre_rmin(pre_rmin_2d, nkin, MAX_PATH); 
  ShapeArray<double, 2> muaq(muaq_2d, ncomp, nspec);
  ShapeArray<double, 2> rate0(rate0_2d, nkin, MAX_PATH); 
  ShapeArray<double, 4> gam(gam_4d, nz, ny, nx, ncomp + nspec);
  ShapeArray<double, 5> keqmin(keqmin_5d, nz, ny, nx, nkin, MAX_PATH); 
  ShapeArray<double, 4> s(s_4d, nz, ny, nx, neqn);

  double affinityTerm, termTMP, sign, term1;
  int ksp;
  int iValue;
  double jac_pre[MAX_PATH][ncomp];

  for (int i = 1; i <= (ncomp + nspec); i++)
    sppTMP[i-1] = sp[jz-1][jy-1][jx-1][i-1];

  for(size_t i = 0; i < nkin; ++i)
    for(size_t j = 0; j < MAX_PATH; ++j)
      for(size_t k = 0; k < ncomp; ++k)
	jac_rmin[i][j][k] = 0.0;

  for (int k = 1; k <= nkin; k++) {

    if ( ivolume[k-1] == 1) {
      continue;
    } else {

      for(size_t i = 0; i < MAX_PATH; ++i)
	for(size_t j = 0; j < ncomp; ++j)
	  jac_pre[i][j] = 0.0;

      for (int loopNP = 1; loopNP <= nreactmin[k-1]; loopNP++) {

	for(size_t j = 0; j < ncomp; ++j)
	  jac_sat[j] = 0.0;

	if (si[k-1][loopNP-1] > 1) {
	  sign = 1;
	} else {
	  sign = -1;
	}

	if (AffinityDepend1[k-1][loopNP-1] == 1) 
	  {
	    term1 = sign * fabs(snorm[k-1][loopNP-1] - 1);
	  } 
	else 
	  {
	    term1 = sign * fabs( pow(snorm[k-1][loopNP-1] - 1,
				     AffinityDepend1[k-1][loopNP-1]));
	  }
	affinityTerm = term1;

	for (int i = 1; i <= ncomp; ++i) {

	  if (mumin[i-1][k-1][loopNP-1] != 0 ) {

	    sppTMP[i-1] = sppTMP[i-1] + PERTURB;

	    termTMP = affinitynumerical(ncomp, 
					jx, 
					jy, 
					jz, 
					nx, 
					ny,
					nz,
					np, 
					nkin, 
					nspec,
					loopNP, 
					k, 
					sppTMP_1d,
					gam_4d, 
					mumin_3d, 
					keqmin_5d, 
					AffinityDepend1_2d);

	    jac_sat[i-1] = (termTMP - affinityTerm)/PERTURB;
	    sppTMP[i-1] = sppTMP[i-1] - PERTURB;
	  }
	}

	int kkMax = ndepend[k-1][loopNP-1];
	for (int kk = 1; kk <= kkMax; ++kk) 
	  {
	    iValue = idepend[k-1][loopNP-1][kk-1];
	    if ( itot_min[k-1][loopNP-1][kk-1] == 1) {
	      for (int i2 = 1; i2 <= ncomp; ++i2) {
		jac_pre[loopNP-1][i2-1] = jac_pre[loopNP-1][i2-1] +
		  fjac[jz-1][jy-1][jx-1][iValue-1][i2-1] *
		  pre_rmin[k-1][loopNP-1] /
		  s[jz-1][jy-1][jx-1][iValue-1];
	      }
	    } else {
	      if (iValue <= ncomp) {
		jac_pre[loopNP-1][iValue-1] = jac_pre[loopNP-1][iValue-1] +
		  pre_rmin[k-1][loopNP-1] *
		  depend[k-1][loopNP-1][kk-1];
	      } else {
		ksp = iValue - ncomp;
		for (int i2 = 1; i2 <= ncomp; ++i2) {
		  jac_pre[loopNP-1][i2-1] = jac_pre[loopNP-1][i2-1] +
		    pre_rmin[k-1][loopNP-1] *
		    depend[k-1][loopNP-1][kk-1] *
		    muaq[i2-1][ksp-1];
		}
	      }
	    }
	  }
	for (int i = 1; i <= ncomp; ++i) {
	  jac_rmin[k-1][loopNP-1][i-1] = surf[k-1] *
	    rate0[k-1][loopNP-1] *
	    (pre_rmin[k-1][loopNP-1] * jac_sat[i-1] + jac_pre[loopNP-1][i-1] * affinityTerm);
	}
      }
    }
  }
  return;
}



