#include <stdbool.h>

#define DEBUG 0

#ifndef max
#define max(a,b)            (((a) > (b)) ? (a) : (b))
#endif

#ifndef min
#define min(a,b)            (((a) < (b)) ? (a) : (b))
#endif

//      INTEGER(I4B)                                         :: nbasis = 5
#define NBASIS 5
#define NTEMP 8

// maximum number of pathways, original np
#define MAX_PATH 10

// maximum number of dependencies, original nd, TODO: set to ncomp+nspec
#define MAX_DEPEND 15

#define MAX_SOL 10

// TODO: ikin is set to 0 in crunchflowloop
#define IKIN 25

//      integer(i4b), parameter                  :: tpoints = 8
#define TPOINTS 8

//      REAL(DP), PARAMETER                                  :: tk=273.15
#define TK 273.15

#ifdef __cplusplus
extern"C" {
#endif

  struct species
  {
    char *name;
    double nu; // stoichiometric coefficient
    struct species *next;
  };

  typedef struct species Species;

  struct component
  {
    char *name;
    double a_zero;
    double charge;
    double molecular_weight;
  };
  
  typedef struct component Component;

  struct complex
  {
    char *name;
    char *stoichiometry;
    double keq[TPOINTS];
    double a_zero;
    double charge;
    double molecular_weight;
  };

  typedef struct complex Complex;

  struct mineral
  {
    char *name;
    char *stoichiometry;
    double keq[TPOINTS];
    double molecular_weight;
    double molar_volume;
    bool read_kinetics;
  };

  typedef struct mineral Mineral;

  struct pathway
  {
    char *name;
    char *label;
    char *type;
    double rate;
    double activation_energy;
    Species *dependencyList;
    double m1;
  };

  typedef struct pathway Pathway;

  struct solution
  {
    char *name;
    Species *concentrationList;
  };

  typedef struct solution Solution;

  extern int ncomp;
  extern int nspec;
  extern int nkin;
  extern int ntemp; // number of discrete temperature values read in input_test.txt file

  extern double bvec[NBASIS]; // DIMENSION(nbasis)
  extern double vec[NTEMP][NBASIS]; // DIMENSION(nbasis,ntemp)
  extern double vecgam[NTEMP][NBASIS]; // DIMENSION(nbasis,ntemp)
  //  extern double temp;
  extern bool RunIsothermal; 

  // extern double as1[NBASIS][12]; // transposed from fortran as1(12,5)
  extern void *as1Global;

  //extern double eqhom[9]; // TODO: malloc nspec * sizeof(double) after initialize.F is translated
  extern void *eqhomGlobal;

  // extern double alnk[MAX_PATH]; // TODO: malloc nkin * max_path * sizeof(double) after initialize.F is translated, reallocate to nptot
  extern void *alnkGlobal;

  extern double tempc[TPOINTS];
  extern double a_dh[TPOINTS];
  extern double b_dh[TPOINTS];
  extern double b_dot[TPOINTS];
  extern double adh[TPOINTS]; // TODO: malloc ntemp * sizeof(double) after initialize.F is translated
  extern double bdh[TPOINTS]; // TODO: malloc ntemp * sizeof(double) after initialize.F is translated
  extern double bdot[TPOINTS]; // TODO: malloc ntemp * sizeof(double) after initialize.F is translated

  //extern char *ulab[MAX_DEPEND]; // TODO: malloc (ncomp + nspec) * sizeof(double) after initialize.F is translated
  extern void *ulabGlobal;

  // extern double chg[MAX_DEPEND]; // TODO: malloc (ncomp + nspec) * sizeof(double) after initialize.F is translated
  extern void *chgGlobal;

  // extern double acmp[MAX_DEPEND]; // TODO: malloc (ncomp + nspec) * sizeof(double) after initialize.F is translated
  extern void *acmpGlobal;

  // extern double wtaq[MAX_DEPEND]; // encapsulated into initialize.cpp as a local auto array
  extern double adhcoeff[NBASIS];
  extern double bdhcoeff[NBASIS];
  extern double bdtcoeff[NBASIS];
  extern Component *components; // calloced array of type Component

  // extern double muaq[6][9]; // TODO: malloc(ncomp * nspec * sizeof(double)
  extern void *muaqGlobal;

  // extern double keqaq[1][32][32][9]; // TODO" malloc nspec*nx*ny*nz* sizeof(double)
  // extern char *namcx[9]; // TODO: malloc(nspec * sizeof(char *)), now local to initialize

  //extern double mumin[6][1][MAX_PATH]; // transpose of (np,nkin,ncomp)
  extern void *muminGlobal;

  // extern double keqmin[1][32][32][1][MAX_PATH]; // transpose of (np,nkin,nx,ny,nz)
  // extern char *namin[1]; // (nkin), now local to initialize

  // extern int ndepend[1][MAX_PATH]; //         ! allocate(ndepend(np,nkin),stat=ierr); ndepend = 0
  extern int (* ndepend)[MAX_PATH];

  // extern int idepend[1][MAX_PATH][MAX_DEPEND]; //      ! allocate(idepend(nd,np,nkin),stat=ierr); idepend = 0    
  extern int (* idepend)[MAX_PATH][MAX_DEPEND];

  // extern double depend[1][MAX_PATH][MAX_DEPEND]; //      ! allocate(depend(nd,np,nkin),stat=ierr); depend = 0.d0   
  extern double (* depend)[MAX_PATH][MAX_DEPEND];

  // extern double rate0[1][MAX_PATH]; //          ! allocate(rate0(np,nkin),stat=ierr); rate0 = 0.d0        
  extern double (* rate0)[MAX_PATH];

  // extern int itot_min[1][MAX_PATH][MAX_DEPEND]; //     ! allocate(itot_min(nd,np,nkin),stat=ierr); itot_min = 0  
  extern int (* itot_min)[MAX_PATH][MAX_DEPEND];

  // extern int nreactmin[1]; //         ! allocate(nreactmin(nkin),stat=ierr); nreactmin = 0      
  extern int (* nreactmin);

  // extern double AffinityDepend1[1][MAX_PATH]; //  ! allocate(AffinityDepend1(np,nkin),stat=ierr) !2011-11-30
  extern double (* AffinityDepend1)[MAX_PATH];

  // extern double volmol[1]; //             ! allocate(volmol(nkin),stat=ierr)
  extern double (* volmol);

  // encapsulated into crunchflowloop
  // nkin,nx,ny,nz
  //
  // extern double dppt[1][32][32][1]; // allocate(dppt(nkin,nx,ny,nz),stat=ierr); dppt = 0.d0    
  // extern double u_rate[1][32][32][1]; // allocate(u_rate(nkin,nx,ny,nz),stat=ierr); dppt = 0.d0    
  // extern double area[1][32][32][1];  // allocate(area(nkin,nx,ny,nz),stat=ierr); area = 100.d0
  // extern double volfx[1][32][32][1]; // allocate(volfx(nkin,nx,ny,nz),stat=ierr); volfx = 1.d0

  // encapsulated into os3d_newton
  // ikin,nx,ny,nz
  //
  // extern double raq_tot[1][32][32][25]; // allocate(raq_tot(ikin,nx,ny,nz)); raq_tot = 0.d0

  // ncomp+nspec,nx,ny,nz
  //
  // extern double sp[1][32][32][MAX_DEPEND]; // allocate(sp(ncomp+nspec,nx,ny,nz),stat=ierr); sp = 0.0d0
  extern void *spGlobal;

  // extern double sp10[1][32][32][MAX_DEPEND]; // allocate(sp10(ncomp+nspec,nx,ny,nz),stat=ierr); sp10 = 0.0d0
  extern void *sp10Global;

  // extern double spold[1][32][32][MAX_DEPEND]; // allocate(spold(ncomp+nspec,nx,ny,nz),stat=ierr); spold = 0.0d0
  extern void *spoldGlobal;

  // extern double gam[1][32][32][MAX_DEPEND]; // allocate(gam(ncomp+nspec,nx,ny,nz),stat=ierr); gam = 0.0d0
  extern void *gamGlobal;

  // ncomp+nspec
  //
  // extern double sppTMP[MAX_DEPEND]; // allocate(sppTMP(ncomp+nspec),stat=ierr); sppTMP = 0.0d0, encapsulate in crunchflowloop

  // neqn
  //
  // extern double sumjackin[6]; // allocate(sumjackin(neqn),stat=ierr); sumjackin = 0.0d0, encapsulate in assemble_local
  // extern double sbdry[6]; // allocate(sbdry(neqn),stat=ierr); sbdry = 0.0d0, not actually used, disabled
  // extern double fxx[6]; // allocate(fxx(neqn),stat=ierr); fxx = 0.0d0
  // extern double fxmax[6]; // allocate(fxmax(neqn),stat=ierr); fxmax = 0.0d0
  // extern double indd[6]; // allocate(indd(neqn),stat=ierr);, not actually used, disabled

  // neqn,nx,ny,nz
  //
  // extern double s[1][32][32][6]; // allocate(s(neqn,nx,ny,nz),stat=ierr); s = 0.0d0
  extern void *sGlobal;

  // extern double sn[1][32][32][6]; // allocate(sn(neqn,nx,ny,nz),stat=ierr); sn = 0.0d0
  extern void *snGlobal;

  // neqn,neqn,nx,ny,nz
  //
  // extern double fjac[1][32][32][6][6]; // allocate(fjac(neqn,neqn,nx,ny,nz),stat=ierr); fjac = 0.0d0, encapsulate in crunchflowloop

  // neqn,neqn
  //
  // extern double fjac_loc[6][6]; // allocate(fjac_loc(neqn,neqn),stat=ierr); fjac_loc = 0.0d0, encapsulate in os3d_newton
 
  // nx,ny,nz
  //
  // extern double sion[1][32][32]; // allocate(sion(nx,ny,nz),stat=ierr); sion = 0.0d0
  extern void *sionGlobal;  

  // extern double t[1][32][32]; // allocate(t(nx,ny,nz),stat=ierr); t(:,:,:)=25.0d0
  extern void *tGlobal;  

  // extern double por[1][32][32]; // allocate(por(nx,ny,nz),stat=ierr), encapsulated in crunchflowloop
  // extern double satliq[1][32][32]; // allocate(satliq(nx,ny,nz),stat=ierr)
  // extern double ro[1][32][32]; // allocate(ro(nx,ny,nz),stat=ierr)   
  // extern double H2Oreacted[1][32][32]; // allocate(H2Oreacted(1:nx,1:ny,1:nz),stat=ierr)

  // nx+3, ny+2, nz+2
  //
  // extern double xgram[1+2][32+2][32+3]; // allocate(xgram(-1:nx+1,0:ny+1,0:nz+1),stat=ierr)
  // extern double xgramOld[1+2][32+2][32+3]; // allocate(xgramOld(-1:nx+1,0:ny+1,0:nz+1),stat=ierr)

  // np x nkin
  //
  // extern double actenergy[1][MAX_PATH]; // allocate(actenergy(np,nkin),stat=ierr); actenergy = 0.d0    
  // extern double silog[1][MAX_PATH];
  // extern double siln[1][MAX_PATH];
  // extern double si[1][MAX_PATH];

  // extern double surf[1];
  // extern double rmin[1][MAX_PATH];
  // extern double pre_rmin[1][MAX_PATH];
  // extern double snorm[1][MAX_PATH];
  // extern double jac_rmin[1][MAX_PATH][6]; // (ncomp,np,nkin)

  // ncomp
  //
  // extern double jac_sat[6];
  // extern double distrib[6];

  // ncomp x ikin
  //
  // extern double mukin[6][25]; // TODO: verify second dimension
  // extern double rdkin[6][25]; // TODO: verify second dimension

  // extern double jac_pre[MAX_PATH][6];
  // extern double decay_correct[1][6]; // allocate(decay_correct(ncomp,nkin),stat=ierr)

  extern bool read_kinetics;
  //  extern double stoichchombo[6][1][1]; // stoichchombo(1,nkinchombo,ncompchombo)

  /* void specieslocal_(int *ncomp, int *nspec, int *jx, int *jy, int *jz, */
  /* 		     double *muaq, double *sp, double *gam, double *sp10, double *keqaq, */
  /* 		     int *nx, int *ny, int *nz); */

  // SIGN(A,B) returns the value of A with the sign of B.
  //
  double sign(const double A, const double B);

  void specieslocal(const int ncomp, 
		    const int nspec, 
		    const int jx, 
		    const int jy, 
		    const int jz, 
		    const int nx, 
		    const int ny, 
		    const int nz,
		    /* const */ double *muaq_2d, // [ncomp][nspec], 
		    double *sp_4d, // [nz][ny][nx][ncomp + nspec], 
		    /* const */ double *gam_4d, // [nz][ny][nx][ncomp + nspec], 
		    double *sp10_4d, // [nz][ny][nx][ncomp + nspec], 
		    /* const */ double *keqaq_4d); // [nz][ny][nx][nspec]);

  double affinitynumerical_cpp_(const int ncomp, 
				const int jx, 
				const int jy, 
				const int jz, 
				const int nx, 
				const int ny,
				const int nz,
				const int np, 
				const int nkin, 
				const int nspec,
				const int loopNP, 
				int k, 
				double *sppTMP,
				const double *gam_4d, // [nz][ny][nx][ncomp + nspec], 
				double *mumin_3d, // [ncomp][nkin][MAX_PATH], 
				const double *keqmin_5d, // [nz][ny][nx][nkin][MAX_PATH], 
				const double *AffinityDepend1_2d); // [nkin][MAX_PATH],

  void jacmin_cpp_(const int ncomp, 
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
		   const double *surf_1d, // [nkin], 
		   double *fjac_5d, // [nz][ny][nx][neqn][neqn], 
		   double *sppTMP,
		   double *sp_4d, // [nz][ny][nx][ncomp + nspec],
		   double *jac_rmin_3d, // [nkin][MAX_PATH][ncomp],
		   const int *ivolume_1d, // [nkin], 
		   double *jac_sat, 
		   const int *nreactmin_1d, // [nkin], 
		   double *si_2d, // [nkin][MAX_PATH], 
		   double *AffinityDepend1_2d, // [nkin][MAX_PATH],
		   double *snorm_2d, // [nkin][MAX_PATH], 
		   double *mumin_3d, // [ncomp][nkin][MAX_PATH],
		   const int *ndepend_2d, // [nkin][np], 
		   const int *idepend_3d, // [nkin][np][nd],
		   const int *itot_min_3d, // [nkin][np][nd],
		   double *depend_3d, // [nkin][MAX_PATH][MAX_DEPEND], 
		   const double *pre_rmin, // [nkin][MAX_PATH], 
		   double *muaq_2d, // [ncomp][nspec],
		   double *rate0_2d, // [nkin][MAX_PATH], 
		   double *gam_4d, // [nz][ny][nx][ncomp + nspec],
		   const double *keqmin_5d, // [nz][ny][nx][nkin][MAX_PATH], 
		   const double *s_4d); // [nz][ny][nx][neqn],

  double aqueoustobulkconvert(const int jx, 
			      const int jy, 
			      const int jz,
			      const int nx,
			      const int ny, 
			      const int nz, 
			      /* const */ double *por_3d, // [nz][ny][nx], 
			      /* const */ double *satliq_3d, // [nz][ny][nx], 
			      /* const */ double *ro_3d); // [nz][ny][nx],

  void totconc(const int ncomp, 
	       const int nspec,
	       const int neqn,
	       const int jx, 
	       const int jy, 
	       const int jz, 
	       const int nx, 
	       const int ny, 
	       const int nz, 
	       /* const */ double *muaq_2d, // [ncomp][nspec], 
	       /* const */ double *sp10_4d, // [nz][ny][nx][ncomp + nspec], 
	       double *s_4d); // [nz][ny][nx][ncomp],

  void satcalc(const int ncomp,
	       const int nspec,
	       const int nkin,
	       const int nrct,
	       const int jx,
	       const int jy,
	       const int jz,
	       const int nx,
	       const int ny,
	       const int nz,
	       double *silog_2d, // [nkin][MAX_PATH],
	       double *mumin_3d, // [ncomp][nkin][MAX_PATH],
	       double *sp_4d, // [nz][ny][nx][ncomp + nspec],
	       double *gam_4d, // [nz][ny][nx][ncomp + nspec],
	       /* const */ double *keqmin_5d, // [nz][ny][nx][nkin][MAX_PATH],
	       /* const */ int *nreactmin_1d); // [nkin],

  void keqcalc2(const int ncomp, 
		const int nrct, 
		const int nspec, 
		const int ngas, 
		const int nkin,
		const int nsurf_sec,
		const int jx, 
		const int jy, 
		const int jz, 
		const int nx, 
		const int ny, 
		const int nz,
		const double *t_3d, // [nz][ny][nx], 
		const bool RunIsoThermal, 
		double *keqaq_4d, // [nz][ny][nx][nspec],
		const int *nreactmin_1d, // [nkin],
		double *keqmin_5d, // [nz][ny][nx][nkin][MAX_PATH],
		const double *alnk_1d, // [nkin * MAX_PATH],
		const double *eqhom_1d, // [nspec],
		const double *as1_2d); // [NBASIS][nspec + nkin * MAX_PATH],

  void fx_local_cpp_(const int ncomp,
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
		     const double *ro_3d, // [nz][ny][nx], 
		     const double *H2Oreacted_3d, // [nz][ny][nx], 
		     const double *xgram_3d, // [nz+2][ny+2][nx+3], 
		     const double *s_4d, // [nz][ny][nx][ncomp], 
		     const double *xgramOld_3d, // [nz+2][ny+2][nx+3],
		     const double *sn_4d, // [nz][ny][nx][ncomp], 
		     double *distrib, 
		     double *fxx_1d, // [neqn], 
		     double *fxmax, // [neqn],
      		     double *satliq_3d); // [nz][ny][nx],

  void reaction_cpp_(const int ncomp,
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
		     const int *nreactmin_1d, // [nkin], 
		     double *actenergy_2d, // [nkin][MAX_PATH],
		     double *silog_2d, // [nkin][MAX_PATH], 
		     double *siln_2d, // [nkin][MAX_PATH],
		     double *si_2d, // [nkin][MAX_PATH], 
		     double *surf_1d, // [nkin], 
		     double *area_4d, // [nz][ny][nx][nkin], 
                     double *rmin_2d, // [nkin][MAX_PATH],
                     double *pre_rmin, // [nkin][MAX_PATH],
 		     const int *ndepend_2d, // [nkin][np], 
		     const int nd, 
		     const int *idepend_3d, // [nkin][np][nd], 
		     double *depend_3d, // [nkin][MAX_PATH][MAX_DEPEND],
		     const int *itot_min_3d, // [nkin][np][nd],
		     const double *s_4d, // [nz][ny][nx][neqn], 
		     const double *gam_4d, // [nz][ny][nx][ncomp + nspec],
		     double *sp_4d, // [nz][ny][nx][ncomp + nspec],
		     double *AffinityDepend1_2d, // [nkin][MAX_PATH],
		     double *rate0_2d, // [nkin][MAX_PATH], 
		     double *volmol, 
		     double *mumin_3d, // [ncomp][nkin][MAX_PATH],
		     const double *keqmin_5d); // [nz][ny][nx][nkin][MAX_PATH],

  int fitgamma_cpp(double *);

  void jac_local_cpp_(const int ncomp, 
		      const int nspec, 
		      const int neqn,
		      const int jx, 
		      const int jy, 
		      const int jz,
		      double *fjac_loc_2d, // [neqn][neqn],
		      const int nx, 
		      const int ny, 
		      const int nz,
		      double *muaq_2d, // [ncomp][nspec], 
		      const double *sp10_4d); // [nz][ny][nx][ncomp + nspec],

  // renamed to cfgamma (crunchflow gamma), so as not to conflict with the gamma function defined in math.h
  //
  void cfgamma(const int ncomp, 
	       const int nspec, 
	       const int jx, 
	       const int jy, 
	       const int jz, 
	       const int nx, 
	       const int ny, 
	       const double *sp10_4d, // [nz][ny][nx][ncomp + nspec],
	       const double *t_3d, // [nz][ny][nx], 
	       double adhcoeff[NBASIS],
	       double bdhcoeff[NBASIS],
	       double bdtcoeff[NBASIS],
	       double *sion_3d, // [nz][ny][nx],
	       char *ulab[MAX_DEPEND], // ncomp+nspec
	       const double *acmp_1d, // [ncomp + nspec]
	       const double *chg_1d, // [ncomp + nspec]
	       const double *gam_4d); // [nz][ny][nx][ncomp + nspec],
  
  void fit_(int *nbasis, int *ntemp, double *alogk0, double *bvec, double *vec,
  	    int *iflgint, int *inoutint, int *ntt, char *nameTransfer);

  void assemble_local_cpp_(const int ncomp, 
			   const int nspec, 
			   const int nkin, 
			   const int ikin, 
			   const int neqn,
			   const double dt,
			   const int nx, 
			   const int ny, 
			   const int nz, // needed for variable-length array raq_tot 
			   const int jx, 
			   const int jy, 
			   const int jz, 
			   double *satliq_3d, // [nz][ny][nx], 
			   double *por_3d, // [nz][ny][nx],
			   const double *ro_3d, // [nz][ny][nx], 
			   const double *rmin_2d, // [nkin][MAX_PATH], 
			   const double *decay_correct_2d, // [nkin][ncomp],
			   double *aaa_2d, // [ncomp][ncomp],
			   const int nradmax, 
			   const int *nreactmin_1d, // [nkin], 
			   const int np,
			   double *fxx_1d, // [neqn], 
			   const double *mumin_3d, // [ncomp][nkin][MAX_PATH],
			   const double *mukin_2d, // [ncomp][ikin],
			   double *raq_tot_4d, // [nz][ny][nx][ikin],
			   const double *jac_rmin_3d, // [nkin][MAX_PATH][ncomp],
			   double *rdkin_2d, // [ncomp][ikin],
			   const int ikh2o, 
			   double *distrib,
			   const double *H2Oreacted_3d, // [nz][ny][nx], 
			   double *fjac_loc_2d, // [neqn][neqn], 
			   const double *xgram_3d); // [nz+2][ny+2][nx+3],

  int initialize(const int ncompchombo, 
		 const int nspecchombo, 
		 const int nkinchombo,
		 double *inflowchombo_1d, // [ncompchombo],
		 double *inflowchombosp10_1d, // [ncompchombo],
		 double *initchombo_1d, // [ncompchombo],
		 double *initchombosp_1d, // [ncompchombo + nspecchombo],
		 double *initchombosp10_1d, // [ncompchombo + nspecchombo],
		 double *stoichchombo_3d, // [ncompchombo][nkinchombo][1], 
		 const int nx, 
		 const int ny, 
		 const int nz);

  int os3d_newton_cpp_(const int ncomp,
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
		       double *t_3d, // [nz][ny][nx],
		       double adhcoeff[NBASIS],
		       double bdhcoeff[NBASIS],
		       double bdtcoeff[NBASIS], 
		       double *sion_3d, // [nz][ny][nx],
		       char *ulab[MAX_DEPEND],
		       const double *acmp_1d, // [ncomp + nspec],
		       const double *chg_1d, // [ncomp + nspec],
		       double *gam_4d, // [nz][ny][nx][ncomp + nspec],
		       double *muaq_2d, // [ncomp][nspec],
		       double *sp_4d, // [nz][ny][nx][ncomp + nspec],
		       double *keqaq_4d, // [nz][ny][nx][nspec],
		       double *s_4d, // [nz][ny][nx][ncomp],
		       double *snorm_2d, // [nkin][MAX_PATH],
		       double *dppt_4d, // [nz][ny][nx][nkin],
		       double *u_rate_4d, // [nz][ny][nx][nkin],
		       const int *nreactmin_1d, // [nkin],
		       double *actenergy_2d, // [nkin][MAX_PATH],
		       double *silog_2d, // [nkin][MAX_PATH],
		       double *siln_2d, // [nkin][MAX_PATH],
		       double *si_2d, // [nkin][MAX_PATH],
		       double *surf_1d, // [nkin],
		       double *area_4d, // [nz][ny][nx][nkin],
		       double *rmin_2d, // [nkin][MAX_PATH],
		       double *pre_rmin, // [nkin][MAX_PATH],
		       const int *ndepend_2d, // [nkin][MAX_PATH],
		       int *idepend_3d, // [nkin][MAX_PATH][MAX_DEPEND],
		       double *depend_3d, // [nkin][MAX_PATH][MAX_DEPEND],
		       int *itot_min_3d, // [nkin][MAX_PATH][MAX_DEPEND],
		       double *AffinityDepend1_2d, // [nkin][MAX_PATH],
		       double *rate0_2d, // [nkin][MAX_PATH],
		       double *volmol,
		       double *mumin_3d, // [ncomp][nkin][MAX_PATH],
		       const double *keqmin_5d, // [nz][ny][nx][nkin][MAX_PATH],
		       double *sppTMP,
		       double *jac_rmin_3d, // [nkin][MAX_PATH][ncomp],
		       const int *ivolume_1d, // [nkin],
		       double *jac_sat,
		       double *fjac_5d, // [nz][ny][nx][neqn][neqn],
		       double *por_3d, // [nz][ny][nx],
		       double *ro_3d, // [nz][ny][nx],
		       const double *H2Oreacted_3d, // [nz][ny][nx],
		       const double *xgram_3d, // [nz+2][ny+2][nx+3],
		       const double *xgramOld_3d, // [nz+2][ny+2][nx+3],
		       double *sn_4d, // [nz][ny][nx][neqn],
		       double *distrib,
		       double *fxx_1d, // [neqn],
		       double *fxmax, // [neqn],
		       double *satliq_3d); // [nz][ny][nx],

  int parseInput(char* filename);

  extern Component *components; // calloced array of type Component
  extern Complex *complexes;
  extern Mineral *minerals;
  extern Pathway *pathways;
  extern Solution *solutions;

  void printSpecies(Species *dependencyList);

#ifdef __cplusplus
}
#endif

