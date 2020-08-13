#include <math.h>
#include "ShapeArray.H"
#include "crunchflow.h"

using namespace shape;

void fx_local(const int ncomp,
	      const int neqn,
	      const double dt, 
	      const int jx, 
	      const int jy, 
	      const int jz, 
	      const int nx,
	      const int ny, 
	      const int nz, 
	      const int ikh2o, 
	      double *por_3d, // [nz][ny][nx],
	      /* const */ double *ro_3d, // [nz][ny][nx],
	      /* const */ double *H2Oreacted_3d, // [nz][ny][nx], 
	      /* const */ double *xgram_3d, // [nz+2][ny+2][nx+3],
	      /* const */ double *s_4d, // [nz][ny][nx][ncomp], 
	      /* const */ double *xgramOld_3d, // [nz+2][ny+2][nx+3],
	      /* const */ double *sn_4d, // [nz][ny][nx][ncomp], 
	      double *distrib, 
	      double *fxx_1d, // [neqn], 
	      double *fxmax_1d, // [neqn],
	      double *satliq_3d) // [nz][ny][nx]
{
  ShapeArray<double, 3> por(por_3d, nz, ny, nx);
  ShapeArray<double, 3> satliq(satliq_3d, nz, ny, nx);
  ShapeArray<double, 3> ro(ro_3d, nz, ny, nx);
  ShapeArray<double, 3> H2Oreacted(H2Oreacted_3d, nz, ny, nx);
  ShapeArray<double, 3> xgram(xgram_3d, nz+2, ny+2, nx+3);
  ShapeArray<double, 3> xgramOld(xgramOld_3d, nz+2, ny+2, nx+3);
  ShapeArray<double, 4> s(s_4d, nz, ny, nx, ncomp);
  ShapeArray<double, 4> sn(sn_4d, nz, ny, nx, ncomp);
  ShapeArray<double, 1> fxx(fxx_1d, neqn);
  ShapeArray<double, 1> fxmax(fxmax_1d, neqn); 
    
  double r, retardation, satl, satgasnew, satlOld, satGasOld, source,
    gas_accum, ex_accum, aq_accum;

  int ind;
 
  r = 1.0/dt;
  retardation = 1.0;

  for(size_t i = 0; i < neqn; ++i)
    fxmax[i] = 0.0;

  satl = satliq[jz-1][jy-1][jx-1];
  satgasnew = 1.0 - satl;
  satlOld = satl;
  satGasOld = 1.0 - satlOld;

  for (int i = 1; i <= ncomp; ++i) {
    ind = i;
    source = 0.0;

    if (i != ikh2o) {
      aq_accum = r *
	por[jz-1][jy-1][jx-1] *
	ro[jz-1][jy-1][jx-1] *
	(H2Oreacted[jz-1][jy-1][jx-1] *
	 xgram[jz-1][jy-1][jx-1] *
	 satl *
	 s[jz-1][jy-1][jx-1][i-1] -
	 xgramOld[jz-1][jy-1][jx-1] *
	 satlOld *
	 sn[jz-1][jy-1][jx-1][i-1]) * (1.0 + retardation * distrib[i-1]);
    } else {
      aq_accum = r * por[jz-1][jy-1][jx-1] *
	ro[jz-1][jy-1][jx-1] *
	(xgram[jz-1][jy-1][jx-1] *
	 satl *
	 s[jz-1][jy-1][jx-1][i-1] -
	 xgramOld[jz-1][jy-1][jx-1] *
	 satlOld *
	 sn[jz-1][jy-1][jx-1][i-1]) * (1.0 + retardation * distrib[i-1]);
    }
 
    gas_accum = 0.0;
    ex_accum = 0.0;
    // array1D(fxx,ind) = aq_accum + gas_accum + ex_accum - source;
    fxx[ind-1] = aq_accum + gas_accum + ex_accum - source;
    if (fabs( fxx[ind-1] ) > fxmax[i-1]) {
      fxmax[i-1] = fabs(fxx[ind-1]);
    }
  }
  return;
}
