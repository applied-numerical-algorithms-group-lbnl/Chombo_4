#include <math.h>
#include <string.h>
#include <stdio.h>
#ifdef MKL
#include <mkl_lapacke.h>
#else
#include <lapacke.h>
#endif
#include <stdlib.h>
#include "ShapeArray.H"
#include "crunchflow.h"

using namespace shape;
//
int cfdriver(enum Target target)
{
  const int ncomp = 6;
  const int nspec = 9;
  const int nkin = 1;
  const int nrct = 1;
  const int ngas = 0;
  const int ikin = 0; 
  const int nexchange = 0;
  const int nexch_sec = 0;
  const int nsurf = 0;
  const int nsurf_sec = 0;
  const int npot = 0;
  const int ndecay = 0;
  const int neqn = 6; 
  const int igamma = 3;
  const double delt = 3.0974155783897026e-10;
  const double corrmax = 2;
  const int jx = 1; 
  const int jy = 1; 
  const int jz = 1;
  int iterat = 0;
  int icvg = 0; 
  const int nx = 1; 
  const int ny = 1; 
  const int nz = 1;
  const double AqueousToBulk = 1000;
  const int ntemp = 1;
  double sp10[nz][ny][nx][ncomp + nspec] = {
    {
      {
	{ 1.094e-05, 1.0699999999999999e-05, 1.5429999999999999e-30, 0.01, 0.01001, 1e-10, 1.1329999999999999e-09, 5.4180000000000005e-07, 1.7759999999999999e-11, 4.984e-39, 6.2359999999999994e-36, 1.6329999999999999e-38, 4.5009999999999999e-07, 2.075e-33, 1.8730000000000001e-35}
      }
    }
  };
  double t[nz][ny][nx]= {
    {
      { 25.0 }
    }
  };
  double adh[TPOINTS] = {5.4535291827807498e-312, 1.1592531020084505e-316, 0, 4.9406564584124654e-324, 0, 2.1219957909652723e-314, 1.1592531020084505e-316, 4.9406564584124654e-324};
  double bdh[TPOINTS] = {1.2404280876201673e-316, 1.2412367742692803e-316, 2.121995791459338e-314, 0, 2.1340731901051351e-314, 4.9406564584124654e-324, 0, 4.9406564584124654e-324};
  double bdot[TPOINTS] = {1.2429521701916411e-316, 1.2401190001521291e-316, 2.1219957909652723e-314, 0, 0, 1.159114763627615e-316, 0, 0};
  double adhcoeff[NBASIS] = {1.3613076615426471e-306, -1.5141925930534484e-307, 5.2314905244148958e-309, -4.0125099678955307e-311, -5.5948338571144667e-313};
  double bdhcoeff[NBASIS] = {-1.6247699247131303e-309, 1.139461441087177e-310, 1.5869502204173258e-312, -2.5815496273979835e-313, 4.653930515140055e-315};
  double bdtcoeff[NBASIS] = {-4.679928850806296e-309, 4.9692372856866895e-310, -1.5145165653969918e-311, 2.424749180310996e-314, 3.4398616844590854e-315}; 
  double sion[nz][ny][nx] = {
    {
      { 0.01001074150202 }
    }
  };
  const char *ulab[MAX_DEPEND] = {"H+", "CO2(aq)", "Ca++", "Na+", "Cl-", "Tracer", "", "", "", ""};
  double acmp[ncomp + nspec] =  {9, 3, 6, 4, 3, 1, 3.5, 4, 4.5, 3, 4, 4, 3, 4, 3};
  double chg[ncomp + nspec] = {1, 0, 2, 1, -1, 0, -1, -1, -2, 0, 1, 1, 0, 1, 0};
  double gam[nz][ny][nx][ncomp + nspec] = {
    {
      {
	{ -1.1235173163548293e-312, 0.0023124812869666201, -4.4939904043901578e-312, -1.1235173163548293e-312, -1.1235173163548293e-312, 0.0023124812869666201, -1.1235173163548293e-312, -1.1235173163548293e-312, -4.4939904043901578e-312, 0.0023124812869666201, -1.1235173163548293e-312, -1.1235173163548293e-312, 0.0023124812869666201, -1.1235173163548293e-312, 0.0023124812869666201 }
      }
    }
  };
  double muaq[ncomp][nspec] = {
    { -1, -1, -2, -2, -1, -1, 1, 0, 0},
    { 0, 1, 1, 1, 1, 0, 0, 0, 0},
    { 0, 0, 0, 1, 1, 1, 0, 1, 1},
    { 0, 0, 0, 0, 0, 0, 0, 0, 0},
    { 0, 0, 0, 0, 0, 0, 1, 1, 2},
    { 0, 0, 0, 0, 0, 0, 0, 0, 0}
  };
 
  double sp[nz][ny][nx][ncomp + nspec] = {
    {
      {
	{ -11.42308476097044, -11.445266816496414, -68.64382421644035, -4.6051701859880909, -4.6041706856550082, -23.025850929940457, -20.598396854900543, -14.428368907295416, -24.754072378364523, -88.194585845282631, -81.062724396590113, -87.007814719785415, -14.613796056647525, -75.255346915120839, -79.962936831329642 }
      }
    }
  };
  double keqaq[nz][ny][nx][nspec] = {
    {
      {
	{ -32.224908634960975, -14.601682190601872, -36.746955506942186, -30.739971518875862, -12.200385280410913, -29.58821841182974, 1.6116483852317498, -1.6127075732676093, -1.5046472548678893 }
      }
    }
  };
  double s[nz][ny][nx][ncomp] = {
    {
      {
	{ 1.084713148e-05, 1.124181776e-05, 1.5450999873139997e-30, 0.01, 0.0100104501, 1e-10 }
      }
    }
  };
  double snorm[nkin][MAX_PATH] = {
    { 0, 0, 0 }
  };
  double dppt[nz][ny][nx][nkin] = {
    {
      {
	{ 0 }
      }
    }
  };
  double u_rate[nz][ny][nx][nkin] = {
    {
      {
	{ 0 }
      }
    }
  };
  /* const */ int nreactmin[nkin] = { 3 };
  double actenergy[nkin][MAX_PATH] = {
    { 0, 0, 0 }
  };
  double silog[nkin][MAX_PATH] = {
    { 0, 0, 0 }
  };
  double siln[nkin][MAX_PATH] = {
    { 0, 0, 0 }
  };
  double si[nkin][MAX_PATH] = {
    { 0, 0, 0 }
  };
  double surf[nkin] = { 0 };
  double area[nz][ny][nx][nkin] = {
    {
      {
	{ 32000 }
      }
    }
  };
  double rmin[nkin][MAX_PATH] = {
    { 0, 0, 0 }
  };
  double pre_rmin[nkin][MAX_PATH] = {
    { 0, 0, 0 }
  };
  int ndepend[nkin][MAX_PATH] = {
    { 1, 1, 0 }
  };
  int idepend[nkin][MAX_PATH][MAX_DEPEND] = {
    {
      { 1 },
      { 2 },
      { 0 }
    }
  };
  double depend[nkin][MAX_PATH][MAX_DEPEND] = {
    {
      { 1 },
      { 1 },
      { 0 }
    }
  };
  int itot_min[nkin][MAX_PATH][MAX_DEPEND] = {
    {
      { 0 },
      { 0 },
      { 0 }
    }
  };
  double AffinityDepend1[nkin][MAX_PATH] = {
    { 1, 1, 1 }
  };
  double rate0[nkin][MAX_PATH] = {
    { 28106489.5849858, 15805.440599669666, 20.835628576367561 }
  };
  double volmol[nkin] = { 3.6933999999999996e-05 };
  double mumin[ncomp][nkin][MAX_PATH] = {
    {
      { -2, -2, -2 }
    },
    {
      { 1, 1, 1 }
    },
    {
      { 1, 1, 1 }
    },
    {
      { 0, 0, 0 }
    },
    {
      { 0, 0, 0 }
    },
    {
      { 0, 0, 0 }
    }
  };
  /* const */ double H2Oreacted[nz][ny][nx] = {
    {
      {1}
    }
  };
  /* const */ double keqmin[nz][ny][nx][nkin][MAX_PATH] = {
    {
      {
	{
	  { -17.220112876465272, -17.220112876465272, -17.220112876465272 }
	}
      }
    }
  };
  double sppTMP[ncomp + nspec] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
  double jac_rmin[nkin][MAX_PATH][ncomp] = {
    {
      { 0, 0, 0, 0, 0, 0},
      { 0, 0, 0, 0, 0, 0},
      { 0, 0, 0, 0, 0, 0}
    }
  };
  int ivolume[nkin] = { 0 };
  double jac_sat[ncomp] = { 0, 0, 0, 0, 0, 0 };
  double fjac[nz][ny][nx][neqn][neqn] = {
    {
      {
	{
	  { 0, 0, 0, 0, 0, 0},
	  { 0, 0, 0, 0, 0, 0},
	  { 0, 0, 0, 0, 0, 0},
	  { 0, 0, 0, 0, 0, 0},
	  { 0, 0, 0, 0, 0, 0},
	  { 0, 0, 0, 0, 0, 0}
	}
      }
    }
  };
  double por[nz][ny][nx] = {
    {
      { 1 }
    }
  };
  double ro[nz][ny][nx] = {
    {
      { 1000 }
    }
  };
  /* const */ double xgram[nz+2][ny+2][nx+3] = {
      {
	{1.0, 1.0, 1.0, 1.0},
	{1.0, 1.0, 1.0, 1.0},
	{1.0, 1.0, 1.0, 1.0}
      },
      {
	{1.0, 1.0, 1.0, 1.0},
	{1.0, 1.0, 1.0, 1.0},
	{1.0, 1.0, 1.0, 1.0}
      },
      {
	{1.0, 1.0, 1.0, 1.0},
	{1.0, 1.0, 1.0, 1.0},
	{1.0, 1.0, 1.0, 1.0}
      }
  };
  /* const */ double xgramOld[nz+2][ny+2][nx+3] = {
      {
	{1.0, 1.0, 1.0, 1.0},
	{1.0, 1.0, 1.0, 1.0},
	{1.0, 1.0, 1.0, 1.0}
      },
      {
	{1.0, 1.0, 1.0, 1.0},
	{1.0, 1.0, 1.0, 1.0},
	{1.0, 1.0, 1.0, 1.0}
      },
      {
	{1.0, 1.0, 1.0, 1.0},
	{1.0, 1.0, 1.0, 1.0},
	{1.0, 1.0, 1.0, 1.0}
      }
  };
  double sn[nz][ny][nx][neqn] = {
    {
      {
	{ 1.0847131479999995e-05, 1.1241817759896801e-05, 1.0319551220714718e-16, 0.0099999999999998927, 0.010010450099999893, 0.0016003316377419285 }
      }
    }
  };
  double distrib[ncomp] = { 0, 0, 0, 0, 0, 0 };
  double fxx[neqn] = { 0, 0, 0, 0, 0, 0 };
  double fxmax[neqn] = { 0, 0, 0, 0, 0, 0 };
  double satliq[nz][ny][nx] = {
    {
      { 1 }
    }
  };

  return( os3d_newton(/* enum Target */ target,
		      /* const int */ ncomp,
		      /* const int */ nspec,
		      /* const int */ nkin,
		      /* const int */ nrct,
		      /* const int */ ngas,
		      /* const int */ ikin, 
		      /* const int */ nexchange,
		      /* const int */ nexch_sec,
		      /* const int */ nsurf,
		      /* const int */ nsurf_sec, // 10
		      /* const int */ npot,
		      /* const int */ ndecay,
		      /* const int */ neqn, 
		      /* const int */ igamma,
		      /* const double */ delt,
		      /* const double */ corrmax,
		      /* const int */ jx, 
		      /* const int */ jy, 
		      /* const int */ jz,
		      /* int *i */ &iterat, // 20
		      /* int *i */ &icvg,
		      /* const int */ nx, 
		      /* const int */ ny, 
		      /* const int */ nz,
		      /* const double */ AqueousToBulk,
		      /* const int */ ntemp,
		      /* double * */ &sp10[0][0][0][0], // [nz][ny][nx][ncomp + nspec],
		      /* const */ /* double * */ &t[0][0][0], // [nz][ny][nx], 
		      /* const */ /* double * */ &adh[0], // [TPOINTS]
		      /* const */ /* double * */ &bdh[0], // [TPOINTS] 30
		      /* const */ /* double * */ &bdot[0], // [TPOINTS]
		      /* double * */ &adhcoeff[0], // [NBASIS],
		      /* double * */ &bdhcoeff[0], // [NBASIS],
		      /* double * */ &bdtcoeff[0], // [NBASIS],
		      /* double * */ &sion[0][0][0], // [nz][ny][nx],
		      /* char ** */ (char **)&ulab[0], // [MAX_DEPEND],
		      /* const */ /* double * */ &acmp[0], // [ncomp + nspec],
		      /* const */ /* double * */ &chg[0], // [ncomp + nspec],
		      /* double * */ &gam[0][0][0][0], // [nz][ny][nx][ncomp + nspec],
		      /* double * */ &muaq[0][0], // [ncomp][nspec], 40
		      /* double * */ &sp[0][0][0][0], // [nz][ny][nx][ncomp + nspec],
		      /* double * */ &keqaq[0][0][0][0], // [nz][ny][nx][nspec],
		      /* double * */ &s[0][0][0][0], // [nz][ny][nx][neqn],
		      /* double * */ &snorm[0][0], // [nkin][MAX_PATH],
		      /* double * */ &dppt[0][0][0][0], // [nz][ny][nx][nkin],
		      /* double * */ &u_rate[0][0][0][0], // [nz][ny][nx][nkin],
		      /* const */ /* int * */ &nreactmin[0], // [nkin],
		      /* double * */ &actenergy[0][0], // [nkin][MAX_PATH],
		      /* double * */ &silog[0][0], // [nkin][MAX_PATH],
		      /* double * */ &siln[0][0], // [nkin][MAX_PATH], 50
		      /* double * */ &si[0][0], // [nkin][MAX_PATH],
		      /* double * */ &surf[0], // [nkin],
		      /* double * */ &area[0][0][0][0], // [nz][ny][nx][nkin],
		      /* double * */ &rmin[0][0], // [nkin][MAX_PATH],
		      /* double * */ &pre_rmin[0][0], // [nkin][MAX_PATH],
		      /* const */ /* int * */ &ndepend[0][0], // [nkin][np],
		      /* const */ /* int * */ &idepend[0][0][0], // [nkin][np][nd],
		      /* double * */ &depend[0][0][0], // [nkin][MAX_PATH][MAX_DEPEND],
		      /* const */ /* int * */ &itot_min[0][0][0], // [nkin][np][nd],
		      /* double * */ &AffinityDepend1[0][0], // [nkin][MAX_PATH], 60
		      /* double * */ &rate0[0][0], // [nkin][MAX_PATH],
		      /* double * */ &volmol[0], // [nkin]
		      /* double * */ &mumin[0][0][0], // [ncomp][nkin][MAX_PATH],
		      /* const */ /* double * */ &keqmin[0][0][0][0][0], // [nz][ny][nx][nkin][MAX_PATH],		
		      /* double * */ &sppTMP[0], // [ncomp + nspec]
		      /* double * */ &jac_rmin[0][0][0], // [nkin][MAX_PATH][ncomp],
		      /* const */ /* int * */ &ivolume[0], // [nkin],
		      /* double * */ &jac_sat[0], // [ncomp],
		      /* double * */ &fjac[0][0][0][0][0], // [nz][ny][nx][neqn][neqn],
		      /* double * */ &por[0][0][0], // [nz][ny][nx], 70
		      /* double * */ &ro[0][0][0], // [nz][ny][nx],
		      /* const */ /* double * */ &H2Oreacted[0][0][0], // [nz][ny][nx],
		      /* const */ /* double * */ &xgram[0][0][0], // [nz+2][ny+2][nx+3],
		      /* const */ /* double * */ &xgramOld[0][0][0], // [nz+2][ny+2][nx+3],
		      /* const */ /* double * */ &sn[0][0][0][0], // [nz][ny][nx][ncomp],
		      /* double * */ &distrib[0], // [ncomp],
		      /* double * */ &fxx[0], // [neqn],
		      /* double * */ &fxmax[0], // [neqn],
		      /* double * */ &satliq[0][0][0])); // [nz][ny][nx] 79
}

int main(int argc, char *argv[])
{
  enum Target target = HOST;
  switch(argc)
    {
    case 1:
      break;
    case 2:
      if(!strcmp(argv[1],"host"))
	target = HOST;
      else if(!strcmp(argv[1],"device"))
	target = DEVICE;
      else
	{
	  fprintf(stderr,"unknown target");
	  exit(1);
	}
      break;
    default:
      fprintf(stderr,"usage: %s [target]\n",argv[0]);
      exit(1);
    }
  return(cfdriver(target));
}
