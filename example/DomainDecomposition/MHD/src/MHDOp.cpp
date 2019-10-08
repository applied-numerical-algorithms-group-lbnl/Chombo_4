#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "Proto.H"
#include "MHDOp.H"
#include "ProtoInterface.H"
#define PI 3.141592653589793

using     Proto::Var;
using     Proto::Stencil;
using     Proto::BoxData;
using     Proto::Point;
using     Proto::Shift;
using     Proto::forall;
using     Proto::forall_p;
typedef   Proto::Var<Real,NUMCOMPS> State;

Real MHDOp::s_gamma = 1.6666666666666666666667;
Real MHDOp::s_dx = 1.0;
Stencil<Real> MHDOp::s_laplacian;
Stencil<Real> MHDOp::s_deconvolve;
Stencil<Real> MHDOp::s_laplacian_f[DIM];
Stencil<Real> MHDOp::s_deconvolve_f[DIM];
Stencil<Real> MHDOp::s_convolve_f[DIM];
Stencil<Real> MHDOp::s_interp_H[DIM];
Stencil<Real> MHDOp::s_interp_L[DIM];
Stencil<Real> MHDOp::s_divergence[DIM];
Stencil<Real> MHDOp::s_ahead_shift[DIM];
Stencil<Real> MHDOp::s_behind_shift[DIM];
Stencil<Real> MHDOp::s_copy_f[DIM];
Copier          MHDOp::s_exchangeCopier;

typedef BoxData<Real,1,1,1> PScalar;
typedef BoxData<Real,NUMCOMPS,1,1> PVector;
using Proto::forall;
using Proto::forallOp;


PROTO_KERNEL_START
void 
consToPrimF(State&         a_W, 
            const State&   a_U,
            Real         a_gamma)
{
    Real rho = a_U(0);
    Real v2 = 0.0;
	Real B2 = 0.0;
    Real gamma = a_gamma;
    a_W(0) = rho;
    
    for (int i = 1; i <= DIM; i++)
      {
        Real v, B;
        v = a_U(i) / rho;
        B = a_U(DIM+1+i);		
        a_W(i) = v;
		a_W(DIM+1+i) = a_U(DIM+1+i);
        v2 += v*v;
		B2 += B*B;
      }
    
    a_W(DIM+1) = (a_U(DIM+1) - .5 * rho * v2 - B2/8.0/PI) * (gamma - 1.0);
}
  PROTO_KERNEL_END(consToPrimF, consToPrim)

PROTO_KERNEL_START
void 
PowellF(State&         a_P, 
        const State&   a_W)
{
#if DIM==1	
	a_P(0) = 0.;
	a_P(1) = a_W(3);
	a_P(2) = a_W(1);
	a_P(3) = a_W(1)*a_W(3);
#endif	

#if DIM==2
	a_P(0) = 0.;
	a_P(1) = a_W(4);
	a_P(2) = a_W(5);
	a_P(3) = a_W(1);
	a_P(4) = a_W(2);
	a_P(5) = a_W(1)*a_W(4) + a_W(2)*a_W(5);
#endif

#if DIM==3	
	a_P(0) = 0.;
	a_P(1) = a_W(5);
	a_P(2) = a_W(6);
	a_P(3) = a_W(7);
	a_P(4) = a_W(1);
	a_P(5) = a_W(2);
	a_P(6) = a_W(3);
	a_P(7) = a_W(1)*a_W(5) + a_W(2)*a_W(6) + a_W(3)*a_W(7);
#endif	
}
  PROTO_KERNEL_END(PowellF, Powell)
  
PROTO_KERNEL_START
void rusanovStateF(State& a_out,
                    const State& a_W_lo,
                    const State& a_W_hi,
					const State& a_F_lo,
                    const State& a_F_hi,
					const Var<Real,1>& a_lambda,
                    int   a_dir,
                    Real a_gamma)
{
    Real gamma = a_gamma;
	
#if DIM == 1
    Real rho_lo, rhou_lo, e_lo, p_lo, Bx_lo, B2_lo, v2_lo;
	Real rho_hi, rhou_hi, e_hi, p_hi, Bx_hi, B2_hi, v2_hi;
	rho_lo  = a_W_lo(0);
	rho_hi  = a_W_hi(0);
	rhou_lo = rho_lo*a_W_lo(1);
	rhou_hi = rho_hi*a_W_hi(1);
	p_lo    = a_W_lo(2);
	p_hi    = a_W_hi(2);
	Bx_lo   = a_W_lo(3);
	Bx_hi   = a_W_hi(3);
	B2_lo   = Bx_lo*Bx_lo;
	B2_hi   = Bx_hi*Bx_hi;
	v2_lo   = a_W_lo(1)*a_W_lo(1);
	v2_hi   = a_W_hi(1)*a_W_hi(1);
	e_lo    = p_lo/(gamma-1.0) + rho_lo*v2_lo/2.0 + B2_lo/8.0/PI;
	e_hi    = p_hi/(gamma-1.0) + rho_hi*v2_hi/2.0 + B2_hi/8.0/PI;
	a_out(0) = 0.5*(a_F_hi(0) + a_F_lo(0) - abs(a_lambda(0))*(rho_hi-rho_lo));
	a_out(1) = 0.5*(a_F_hi(1) + a_F_lo(1) - abs(a_lambda(0))*(rhou_hi-rhou_lo));
	a_out(2) = 0.5*(a_F_hi(2) + a_F_lo(2) - abs(a_lambda(0))*(e_hi-e_lo));
	a_out(3) = 0.5*(a_F_hi(3) + a_F_lo(3) - abs(a_lambda(0))*(Bx_hi-Bx_lo));	
#endif    
	
#if DIM == 2
    Real rho_lo, rhou_lo, rhov_lo, e_lo, p_lo, Bx_lo, By_lo, B2_lo, v2_lo;
	Real rho_hi, rhou_hi, rhov_hi, e_hi, p_hi, Bx_hi, By_hi, B2_hi, v2_hi;
	rho_lo  = a_W_lo(0);
	rho_hi  = a_W_hi(0);
	rhou_lo = rho_lo*a_W_lo(1);
	rhou_hi = rho_hi*a_W_hi(1);
	rhov_lo = rho_lo*a_W_lo(2);
	rhov_hi = rho_hi*a_W_hi(2);
	p_lo    = a_W_lo(3);
	p_hi    = a_W_hi(3);
	Bx_lo   = a_W_lo(4);
	Bx_hi   = a_W_hi(4);
	By_lo   = a_W_lo(5);
	By_hi   = a_W_hi(5);
	B2_lo   = Bx_lo*Bx_lo + By_lo*By_lo;
	B2_hi   = Bx_hi*Bx_hi + By_hi*By_hi;
	v2_lo   = a_W_lo(1)*a_W_lo(1) + a_W_lo(2)*a_W_lo(2);
	v2_hi   = a_W_hi(1)*a_W_hi(1) + a_W_hi(2)*a_W_hi(2);
	e_lo    = p_lo/(gamma-1.0) + rho_lo*v2_lo/2.0 + B2_lo/8.0/PI;
	e_hi    = p_hi/(gamma-1.0) + rho_hi*v2_hi/2.0 + B2_hi/8.0/PI;
	a_out(0) = 0.5*(a_F_hi(0) + a_F_lo(0) - abs(a_lambda(0))*(rho_hi-rho_lo));
	a_out(1) = 0.5*(a_F_hi(1) + a_F_lo(1) - abs(a_lambda(0))*(rhou_hi-rhou_lo));
	a_out(2) = 0.5*(a_F_hi(2) + a_F_lo(2) - abs(a_lambda(0))*(rhov_hi-rhov_lo));
	a_out(3) = 0.5*(a_F_hi(3) + a_F_lo(3) - abs(a_lambda(0))*(e_hi-e_lo));
	a_out(4) = 0.5*(a_F_hi(4) + a_F_lo(4) - abs(a_lambda(0))*(Bx_hi-Bx_lo));
	a_out(5) = 0.5*(a_F_hi(5) + a_F_lo(5) - abs(a_lambda(0))*(By_hi-By_lo));
#endif 

#if DIM == 3
    Real rho_lo, rhou_lo, rhov_lo, rhow_lo, e_lo, p_lo, Bx_lo, By_lo, Bz_lo, B2_lo, v2_lo;
	Real rho_hi, rhou_hi, rhov_hi, rhow_hi, e_hi, p_hi, Bx_hi, By_hi, Bz_hi, B2_hi, v2_hi;
	rho_lo  = a_W_lo(0);
	rho_hi  = a_W_hi(0);
	rhou_lo = rho_lo*a_W_lo(1);
	rhou_hi = rho_hi*a_W_hi(1);
	rhov_lo = rho_lo*a_W_lo(2);
	rhov_hi = rho_hi*a_W_hi(2);
	rhow_lo = rho_lo*a_W_lo(3);
	rhow_hi = rho_hi*a_W_hi(3);
	p_lo    = a_W_lo(4);
	p_hi    = a_W_hi(4);
	Bx_lo   = a_W_lo(5);
	Bx_hi   = a_W_hi(5);
	By_lo   = a_W_lo(6);
	By_hi   = a_W_hi(6);
	Bz_lo   = a_W_lo(7);
	Bz_hi   = a_W_hi(7);
	B2_lo   = Bx_lo*Bx_lo + By_lo*By_lo + Bz_lo*Bz_lo;
	B2_hi   = Bx_hi*Bx_hi + By_hi*By_hi + Bz_hi*Bz_hi;
	v2_lo   = a_W_lo(1)*a_W_lo(1) + a_W_lo(2)*a_W_lo(2) + a_W_lo(3)*a_W_lo(3);
	v2_hi   = a_W_hi(1)*a_W_hi(1) + a_W_hi(2)*a_W_hi(2) + a_W_hi(3)*a_W_hi(3);
	e_lo    = p_lo/(gamma-1.0) + rho_lo*v2_lo/2.0 + B2_lo/8.0/PI;
	e_hi    = p_hi/(gamma-1.0) + rho_hi*v2_hi/2.0 + B2_hi/8.0/PI;
	a_out(0) = 0.5*(a_F_hi(0) + a_F_lo(0) - abs(a_lambda(0))*(rho_hi-rho_lo));
	a_out(1) = 0.5*(a_F_hi(1) + a_F_lo(1) - abs(a_lambda(0))*(rhou_hi-rhou_lo));
	a_out(2) = 0.5*(a_F_hi(2) + a_F_lo(2) - abs(a_lambda(0))*(rhov_hi-rhov_lo));
	a_out(3) = 0.5*(a_F_hi(3) + a_F_lo(3) - abs(a_lambda(0))*(rhow_hi-rhow_lo));
	a_out(4) = 0.5*(a_F_hi(4) + a_F_lo(4) - abs(a_lambda(0))*(e_hi-e_lo));
	a_out(5) = 0.5*(a_F_hi(5) + a_F_lo(5) - abs(a_lambda(0))*(Bx_hi-Bx_lo));
	a_out(6) = 0.5*(a_F_hi(6) + a_F_lo(6) - abs(a_lambda(0))*(By_hi-By_lo));
	a_out(7) = 0.5*(a_F_hi(7) + a_F_lo(7) - abs(a_lambda(0))*(Bz_hi-Bz_lo));
#endif 	
}
PROTO_KERNEL_END(rusanovStateF, rusanovState)


PROTO_KERNEL_START
void getFluxF(State& a_F, const State& a_W, 
              int    a_d,
              Real a_gamma)
  {
	Real rho, p0, e, v2, B2, vB;  
	Real gamma = a_gamma;
	Real mult1=0.;
	Real mult2=0.;
	Real mult3=0.;	
	if (a_d == 0){mult1=1.0;}
	if (a_d == 1){mult2=1.0;}
	if (a_d == 2){mult3=1.0;}
	rho = a_W(0);
#if DIM==1	
    v2 = a_W(1)*a_W(1);
	B2 = a_W(3)*a_W(3); 
	vB = a_W(1)*a_W(3);
    p0 = a_W(2) + B2/8.0/PI;
	e  = a_W(2)/(gamma-1.0) + rho*v2/2.0 + B2/8.0/PI;
	a_F(0) = rho*a_W(1+a_d);
	a_F(1) = rho*a_W(1)*a_W(1+a_d) + mult1*p0 - a_W(3)*a_W(3+a_d)/4.0/PI;
	a_F(2) = (e+p0)*a_W(1+a_d) - a_W(3+a_d)*vB/4.0/PI;
	a_F(3) = 0.0;
#endif
	
#if DIM==2	
	v2 = a_W(1)*a_W(1) + a_W(2)*a_W(2);
	B2 = a_W(4)*a_W(4) + a_W(5)*a_W(5);
    vB = a_W(1)*a_W(4) + a_W(2)*a_W(5);
    p0 = a_W(3) + B2/8.0/PI;
    e  = a_W(3)/(gamma-1.0) + rho*v2/2.0 + B2/8.0/PI;
    a_F(0) = rho*a_W(1+a_d);
	a_F(1) = rho*a_W(1)*a_W(1+a_d) + mult1*p0 - a_W(4)*a_W(4+a_d)/4.0/PI;
	a_F(2) = rho*a_W(2)*a_W(1+a_d) + mult2*p0 - a_W(5)*a_W(4+a_d)/4.0/PI;
	a_F(3) = (e+p0)*a_W(1+a_d) - a_W(4+a_d)*vB/4.0/PI;
	a_F(4) = mult2*(a_W(1+a_d)*a_W(4) - a_W(1)*a_W(4+a_d));
	a_F(5) = mult1*(a_W(1+a_d)*a_W(5) - a_W(2)*a_W(4+a_d));	
#endif

#if DIM==3	
    v2 = a_W(1)*a_W(1) + a_W(2)*a_W(2) + a_W(3)*a_W(3);
	B2 = a_W(5)*a_W(5) + a_W(6)*a_W(6) + a_W(7)*a_W(7);
    vB = a_W(1)*a_W(5) + a_W(2)*a_W(6) + a_W(3)*a_W(7);
    p0 = a_W(4) + B2/8.0/PI;
    e  = a_W(4)/(gamma-1.0) + rho*v2/2.0 + B2/8.0/PI;
    a_F(0) = rho*a_W(1+a_d);
	a_F(1) = rho*a_W(1)*a_W(1+a_d) + mult1*p0 - a_W(5)*a_W(5+a_d)/4.0/PI;
	a_F(2) = rho*a_W(2)*a_W(1+a_d) + mult2*p0 - a_W(6)*a_W(5+a_d)/4.0/PI;
	a_F(3) = rho*a_W(3)*a_W(1+a_d) + mult3*p0 - a_W(7)*a_W(5+a_d)/4.0/PI;
	a_F(4) = (e+p0)*a_W(1+a_d) - a_W(5+a_d)*vB/4.0/PI;
	a_F(5) = (mult2+mult3)*(a_W(1+a_d)*a_W(5) - a_W(1)*a_W(5+a_d));
	a_F(6) = (mult1+mult3)*(a_W(1+a_d)*a_W(6) - a_W(2)*a_W(5+a_d));
	a_F(7) = (mult1+mult2)*(a_W(1+a_d)*a_W(7) - a_W(3)*a_W(5+a_d));
#endif	  
  }
PROTO_KERNEL_END(getFluxF, getFlux)



PROTO_KERNEL_START
void waveSpeedBoundF(Var<Real,1>& a_speed,
                     const State& a_W,
                     Real       a_gamma)
{
	Real gamma = a_gamma;
	Real rho=0., u=0., v=0., w=0., p=0., Bx=0., By=0., Bz=0., ce, af, B_mag, Bdir, udir;
	
	#if DIM == 1
	rho = a_W(0);
	u   = a_W(1);
	p   = a_W(2);
    Bx  = a_W(3);
#endif	
#if DIM == 2
	rho = a_W(0);
	u   = a_W(1);
	v   = a_W(2);
	p   = a_W(3);
    Bx  = a_W(4);
    By  = a_W(5);
#endif	
#if DIM == 3
	rho = a_W(0);
	u   = a_W(1);
	v   = a_W(2);
	w   = a_W(3);
	p   = a_W(4);
    Bx  = a_W(5);
    By  = a_W(6);
    Bz  = a_W(7);
#endif	
    a_speed(0) = 0.0;
	for (int dir = 0; dir< DIM; dir++){
		if (dir == 0){
			Bdir = Bx;
			udir = u;
		};
		if (dir == 1){
			Bdir = By;
			udir = v;
		};
		if (dir == 2){
			Bdir = Bz;
			udir = w;
		};
	
	ce = sqrt(gamma*p/rho);
	B_mag = sqrt(Bx*Bx+By*By+Bz*Bz);
	af = 0.5*(sqrt((ce*ce)+( B_mag*B_mag/(4.0*PI*rho) )+( abs(Bdir)*ce/sqrt(PI*rho) ))+
		  sqrt((ce*ce)+( B_mag*B_mag/(4.0*PI*rho) )-( abs(Bdir)*ce/sqrt(PI*rho) ))) + abs(udir);
	if (af > a_speed(0)){a_speed(0) = af;}
	}
}
PROTO_KERNEL_END(waveSpeedBoundF, waveSpeedBound)


PROTO_KERNEL_START
void lambdacalcF(Var<Real,1>& a_lambda,
                 const State& a_W_low,
                 const State& a_W_high,
				 int a_d,
                 Real a_gamma)
{
    Real gamma = a_gamma;
	Real rho=0., u=0., v=0., w=0., p=0., Bx=0., By=0., Bz=0., ce, af, B_mag, Bdir, udir;
#if DIM == 1
	rho = 0.5 * (a_W_low(0) + a_W_high(0));
	u   = 0.5 * (a_W_low(1) + a_W_high(1));
	p   = 0.5 * (a_W_low(2) + a_W_high(2));
    Bx  = 0.5 * (a_W_low(3) + a_W_high(3));
#endif	
#if DIM == 2
	rho = 0.5 * (a_W_low(0) + a_W_high(0));
	u   = 0.5 * (a_W_low(1) + a_W_high(1));
	v   = 0.5 * (a_W_low(2) + a_W_high(2));
	p   = 0.5 * (a_W_low(3) + a_W_high(3));
    Bx  = 0.5 * (a_W_low(4) + a_W_high(4));
    By  = 0.5 * (a_W_low(5) + a_W_high(5));
#endif	
#if DIM == 3
	rho = 0.5 * (a_W_low(0) + a_W_high(0));
	u   = 0.5 * (a_W_low(1) + a_W_high(1));
	v   = 0.5 * (a_W_low(2) + a_W_high(2));
	w   = 0.5 * (a_W_low(3) + a_W_high(3));
	p   = 0.5 * (a_W_low(4) + a_W_high(4));
    Bx  = 0.5 * (a_W_low(5) + a_W_high(5));
    By  = 0.5 * (a_W_low(6) + a_W_high(6));
    Bz  = 0.5 * (a_W_low(7) + a_W_high(7));
#endif	
    if (a_d == 0){
		Bdir = Bx;
		udir = u;
	};
	if (a_d == 1){
		Bdir = By;
		udir = v;
	};
	if (a_d == 2){
		Bdir = Bz;
		udir = w;
	};
	
	ce = sqrt(gamma*p/rho);
	B_mag = sqrt(Bx*Bx+By*By+Bz*Bz);
	af = 0.5*(sqrt((ce*ce)+( B_mag*B_mag/(4.0*PI*rho) )+( abs(Bdir)*ce/sqrt(PI*rho) ))+
		  sqrt((ce*ce)+( B_mag*B_mag/(4.0*PI*rho) )-( abs(Bdir)*ce/sqrt(PI*rho) )));
	a_lambda(0) = af + abs(udir);
}
PROTO_KERNEL_END(lambdacalcF, lambdacalc) 


PROTO_KERNEL_START
void BavgcalcF(State& a_Bavg,
                 const State& a_W_ave,
				 int a_d)
{
	for (int i=0; i< NUMCOMPS; i++){
		a_Bavg(i) = a_W_ave(2+DIM+a_d);
	}
}
PROTO_KERNEL_END(BavgcalcF, Bavgcalc) 


PROTO_KERNEL_START
void del_W_f_m_calcF(State& a_del_W_f_m,
                 const State& a_W_ave,
                 const State& a_W_ave_high)
{
	for (int i=0; i< NUMCOMPS; i++){
		a_del_W_f_m(i) = a_W_ave(i) - a_W_ave_high(i);
	}
}
PROTO_KERNEL_END(del_W_f_m_calcF, del_W_f_m_calc) 

PROTO_KERNEL_START
void del_W_f_p_calcF(State& a_del_W_f_p,
                 const State& a_W_ave,
                 const State& a_W_ave_low_ahead)
{
	for (int i=0; i< NUMCOMPS; i++){
		a_del_W_f_p(i) = a_W_ave_low_ahead(i) - a_W_ave(i);
	}
}
PROTO_KERNEL_END(del_W_f_p_calcF, del_W_f_p_calc)


PROTO_KERNEL_START
void del2_W_f_calcF(State& del2_W_f,
                 const State& a_W_ave,
                 const State& a_W_ave_high,
                 const State& a_W_ave_low_ahead)
{
	for (int i=0; i< NUMCOMPS; i++){
		del2_W_f(i) = 6.0*(a_W_ave_high(i) - 2.0*a_W_ave(i) + a_W_ave_low_ahead(i));
	}
}
PROTO_KERNEL_END(del2_W_f_calcF, del2_W_f_calc)


PROTO_KERNEL_START
void del2_W_c_calcF(State& del2_W_c,
                 const State& a_W_ave,
                 const State& a_W_ave_behind,
                 const State& a_W_ave_ahead)
{
	for (int i=0; i< NUMCOMPS; i++){
		del2_W_c(i) = a_W_ave_behind(i) - 2.0*a_W_ave(i) + a_W_ave_ahead(i);
	}
}
PROTO_KERNEL_END(del2_W_c_calcF, del2_W_c_calc)


PROTO_KERNEL_START
void del3_W_calcF(State& a_del3_W,
                 const State& a_del2_W_c_ahead,
                 const State& a_del2_W_c)
{
	for (int i=0; i< NUMCOMPS; i++){
		a_del3_W(i) = a_del2_W_c_ahead(i) - a_del2_W_c(i);
	}
}
PROTO_KERNEL_END(del3_W_calcF, del3_W_calc)



PROTO_KERNEL_START
void del2_W_lim_calcF(State& a_del2_W_lim,
                 const State& a_del2_W_f,
                 const State& a_del2_W_c,
                 const State& a_del2_W_c_ahead,
                 const State& a_del2_W_c_behind)
{
	for (int i=0; i< NUMCOMPS; i++){
		if ((a_del2_W_c_behind(i) >= 0.0) && (a_del2_W_c(i) >= 0.0) && (a_del2_W_c_ahead(i) >= 0.0) && (a_del2_W_f(i) >= 0.0)){
			a_del2_W_lim(i) = std::min({std::abs(a_del2_W_f(i)), 1.25*std::abs(a_del2_W_c_behind(i)),
			1.25*std::abs(a_del2_W_c(i)), 1.25*std::abs(a_del2_W_c_ahead(i))} );										
		} else if ((a_del2_W_c_behind(i) < 0.0) && (a_del2_W_c(i) < 0.0) && (a_del2_W_c_ahead(i) < 0.0) && (a_del2_W_f(i) < 0.0)){
			a_del2_W_lim(i) = -1.0*std::min({std::abs(a_del2_W_f(i)), 1.25*std::abs(a_del2_W_c_behind(i)),
			1.25*std::abs(a_del2_W_c(i)), 1.25*std::abs(a_del2_W_c_ahead(i))} );
		} else {
			a_del2_W_lim(i) = 0.0;
		}	
	}
}
PROTO_KERNEL_END(del2_W_lim_calcF, del2_W_lim_calc)


PROTO_KERNEL_START
void rho_i_calcF(State& a_rho_i,
                 const State& a_del2_W_f,
                 const State& a_del2_W_lim,
                 const State& a_W_ave,
                 const State& a_W_ave_ahead,
                 const State& a_W_ave_ahead2,
                 const State& a_W_ave_behind,
                 const State& a_W_ave_behind2)
{
	for (int i=0; i< NUMCOMPS; i++){
		double rhs = 1.0e-12*std::max({std::abs(a_W_ave_behind2(i)), std::abs(a_W_ave_behind(i)), std::abs(a_W_ave(i)), std::abs(a_W_ave_ahead(i)), std::abs(a_W_ave_ahead2(i))});
		double lhs = std::abs(a_del2_W_f(i));
		if (lhs <= rhs){	
			a_rho_i(i) = 0.0;
		} else {
			a_rho_i(i) = a_del2_W_lim(i)/a_del2_W_f(i);
		}			
	}
}
PROTO_KERNEL_END(rho_i_calcF, rho_i_calc)


PROTO_KERNEL_START
void del3_W_min_calcF(State& a_del3_W_min,
                 const State& a_del3_W_behind,
                 const State& a_del3_W,
                 const State& a_del3_W_ahead,
                 const State& a_del3_W_ahead2)
{
	for (int i=0; i< NUMCOMPS; i++){
		a_del3_W_min(i) = std::min({a_del3_W_behind(i), a_del3_W(i), a_del3_W_ahead(i), a_del3_W_ahead2(i)});			
	}
}
PROTO_KERNEL_END(del3_W_min_calcF, del3_W_min_calc)


PROTO_KERNEL_START
void del3_W_max_calcF(State& a_del3_W_max,
                 const State& a_del3_W_behind,
                 const State& a_del3_W,
                 const State& a_del3_W_ahead,
                 const State& a_del3_W_ahead2)
{
	for (int i=0; i< NUMCOMPS; i++){
		a_del3_W_max(i) = std::max({a_del3_W_behind(i), a_del3_W(i), a_del3_W_ahead(i), a_del3_W_ahead2(i)});			
	}
}
PROTO_KERNEL_END(del3_W_max_calcF, del3_W_max_calc)


PROTO_KERNEL_START
void limiter_low_calcF(State& a_W_ave_low_ahead_limited,
                 const State& a_W_ave_low_ahead,
                 const State& a_del_W_f_m,
                 const State& a_del_W_f_p,
                 const State& a_W_ave,
                 const State& a_W_ave_ahead2,
                 const State& a_W_ave_behind2,
                 const State& a_rho_i,
                 const State& a_del3_W_max,
                 const State& a_del3_W_min)
{
	double rhs = 1.0 - 1.0e-12;
	double rhs_test, lhs_test, rhs_test2, lhs_test2;
	
	for (int i=0; i< NUMCOMPS; i++){
        if (a_rho_i(i) >= rhs){
			rhs_test = a_del3_W_max(i) - a_del3_W_min(i);
	        lhs_test = 0.1*std::max({std::abs(a_del3_W_min(i)), std::abs(a_del3_W_max(i))});
			if (lhs_test <= rhs_test){
				if ((a_del_W_f_m(i) * a_del_W_f_p(i) <= 0.0) || ((a_W_ave(i) - a_W_ave_behind2(i))*(a_W_ave_ahead2(i) - a_W_ave(i)) <= 0.0)){
					rhs_test2 = 2.0*std::abs(a_del_W_f_m(i));
					lhs_test2 = std::abs(a_del_W_f_p(i));
					if (a_del_W_f_m(i) * a_del_W_f_p(i) < 0.0){
						a_W_ave_low_ahead_limited(i) = a_W_ave(i) + a_rho_i(i)*a_del_W_f_p(i);
					} else if (lhs_test2 >= rhs_test2){
						a_W_ave_low_ahead_limited(i) = a_W_ave(i) + 2.0*(1.0-a_rho_i(i))*a_del_W_f_m(i) + a_rho_i(i)*a_del_W_f_p(i);
					} else {
						a_W_ave_low_ahead_limited(i) = a_W_ave_low_ahead(i);
					}
				} else {
					rhs_test2 = 2.0*std::abs(a_del_W_f_m(i));
					lhs_test2 = std::abs(a_del_W_f_p(i));
					if (lhs_test2 >= rhs_test2){
						a_W_ave_low_ahead_limited(i) = a_W_ave(i) + 2.0*a_del_W_f_m(i);
					} else {
						a_W_ave_low_ahead_limited(i) = a_W_ave_low_ahead(i);
					}
				}		
			} else {
				a_W_ave_low_ahead_limited(i) = a_W_ave_low_ahead(i);
			}
		} else {
			a_W_ave_low_ahead_limited(i) = a_W_ave_low_ahead(i);
		}
		
	}
}
PROTO_KERNEL_END(limiter_low_calcF, limiter_low_calc)


PROTO_KERNEL_START
void limiter_high_calcF(State& a_W_ave_high_limited,
                 const State& a_W_ave_high,
                 const State& a_del_W_f_m,
                 const State& a_del_W_f_p,
                 const State& a_W_ave,
                 const State& a_W_ave_ahead2,
                 const State& a_W_ave_behind2,
                 const State& a_rho_i,
                 const State& a_del3_W_max,
                 const State& a_del3_W_min)
{
	double rhs = 1.0 - 1.0e-12;
	double rhs_test, lhs_test, rhs_test2, lhs_test2;
	
	for (int i=0; i< NUMCOMPS; i++){
        if (a_rho_i(i) >= rhs){
			rhs_test = a_del3_W_max(i) - a_del3_W_min(i);
	        lhs_test = 0.1*std::max({std::abs(a_del3_W_min(i)), std::abs(a_del3_W_max(i))});
			if (lhs_test <= rhs_test){
				if ((a_del_W_f_m(i) * a_del_W_f_p(i) <= 0.0) || ((a_W_ave(i) - a_W_ave_behind2(i))*(a_W_ave_ahead2(i) - a_W_ave(i)) <= 0.0)){
					rhs_test2 = 2.0*std::abs(a_del_W_f_p(i));
					lhs_test2 = std::abs(a_del_W_f_m(i));
					if (a_del_W_f_m(i) * a_del_W_f_p(i) < 0.0){
						a_W_ave_high_limited(i) = a_W_ave(i) - a_rho_i(i)*a_del_W_f_m(i);
					} else if (lhs_test2 >= rhs_test2){
						a_W_ave_high_limited(i) = a_W_ave(i) - 2.0*(1.0-a_rho_i(i))*a_del_W_f_p(i) - a_rho_i(i)*a_del_W_f_m(i);
					} else {
						a_W_ave_high_limited(i) = a_W_ave_high(i);
					}
				} else {
					rhs_test2 = 2.0*std::abs(a_del_W_f_p(i));
					lhs_test2 = std::abs(a_del_W_f_m(i));
					if (lhs_test2 >= rhs_test2){
						a_W_ave_high_limited(i) = a_W_ave(i) - 2.0*a_del_W_f_p(i);
					} else {
						a_W_ave_high_limited(i) = a_W_ave_high(i);
					}
				}		
			} else {
				a_W_ave_high_limited(i) = a_W_ave_high(i);
			}
		} else {
			a_W_ave_high_limited(i) = a_W_ave_high(i);
		}
		
	}
}
PROTO_KERNEL_END(limiter_high_calcF, limiter_high_calc)

void
MHDOp::
define(const DisjointBoxLayout& a_grids,
       const IntVect          & a_ghostVect) 
{
  CH_TIME("MHDOp::define");
  s_laplacian = Stencil<Real>::Laplacian();
  s_deconvolve = (-1.0/24.0)*s_laplacian + (1.0)*Shift(Point::Zeros());
  for (int dir = 0; dir < DIM; dir++)
  {
    s_laplacian_f[dir] = Stencil<Real>::LaplacianFace(dir);
    s_deconvolve_f[dir] = (-1.0/24.0)*s_laplacian_f[dir] + 1.0*Shift(Point::Zeros());
    s_convolve_f[dir] = (1.0/24.0)*s_laplacian_f[dir] + 1.0*Shift(Point::Zeros());
    s_interp_H[dir] = Stencil<Real>::CellToEdgeH(dir);
    s_interp_L[dir] = Stencil<Real>::CellToEdgeL(dir);
    s_divergence[dir] = Stencil<Real>::FluxDivergence(dir);
	s_ahead_shift[dir] = 1.0*Shift(Point::Basis(dir));
	s_behind_shift[dir] = 1.0*Shift(Point::Basis(dir)*(-1));
	s_copy_f[dir] = 1.0*Shift(Point::Zeros());
  }
  
  s_exchangeCopier.exchangeDefine(a_grids, a_ghostVect);
}

Real 
MHDOp::
proto_step(BoxData<Real,NUMCOMPS>& a_Rhs,
           const BoxData<Real,NUMCOMPS>& a_U,
           const Bx& a_rangeBox)
{
  CH_TIMERS("MHDOp::step(boxdata)");
  CH_TIMER("interp stencil eval", tint);
  CH_TIMER("riemann problem", trie);
  CH_TIMER("get_flux", tgf);
  CH_TIMER("deconvolve", tdcv);
  CH_TIMER("convolve", tcv);
  CH_TIMER("get_flux2", tgf2);
  CH_TIMER("laplacian_arith", tlap);
  CH_TIMER("divergence_eval", tdiv);
  CH_TIMER("divide_by_dx", tdx);
  a_Rhs.setVal(0.0); 
  Real gamma = s_gamma;
  Real retval;
  unsigned long long int sqrtnum = 10;  //this one is just a guess
  unsigned long long int ctoprmnum  = 4*DIM + 5;
  unsigned long long int getfluxnum = 9 + DIM*3;
  unsigned long long int wavespdnum = sqrtnum +3 + DIM;
  unsigned long long int lambdacalcnum = sqrtnum +3 + DIM; // Check
  unsigned long long int rusanovnum = sqrtnum +4 + DIM; // Check
  unsigned long long int del_W_f_m_c = sqrtnum +4 + DIM; // Check
  unsigned long long int del_W_f_p_c = sqrtnum +4 + DIM; // Check
  unsigned long long int del2_W_f_c = sqrtnum +4 + DIM; // Check
  unsigned long long int del2_W_c_c = sqrtnum +4 + DIM; // Check
  unsigned long long int del3_W_c = sqrtnum +4 + DIM; // Check
  unsigned long long int del2_W_lim_c = sqrtnum +4 + DIM; // Check
  unsigned long long int rho_i_c = sqrtnum +4 + DIM; // Check
  unsigned long long int del3_W_min_c = sqrtnum +4 + DIM; // Check
  unsigned long long int del3_W_max_c = sqrtnum +4 + DIM; // Check
  unsigned long long int limiter_low_c = sqrtnum +4 + DIM; // Check
  unsigned long long int limiter_high_c = sqrtnum +4 + DIM; // Check
  PVector W_bar = forallOp<Real,NUMCOMPS>(ctoprmnum, "consToPrim", consToPrim,a_U, gamma);
  PVector U = s_deconvolve(a_U);
  PVector W = forallOp<Real,NUMCOMPS>(ctoprmnum, "consToPrim", consToPrim,U, gamma);
  PScalar umax = forallOp<Real>(wavespdnum, "wavespeed",waveSpeedBound,a_rangeBox,W, gamma);
  retval = umax.absMax();
  PVector W_ave = s_laplacian(W_bar);
  W_ave *= (1.0/24.0);
  W_ave += W;
  for (int d = 0; d < DIM; d++)
  {
    CH_START(tint);
    PVector W_ave_low = s_interp_L[d](W_ave);
    PVector W_ave_high = s_interp_H[d](W_ave);
	
	PVector W_ave_low_ahead = s_ahead_shift[d](W_ave_low);
	PVector del_W_f_m = forallOp<Real,NUMCOMPS>(del_W_f_m_c, "del_W_f_m_calc", del_W_f_m_calc, W_ave, W_ave_high);
	PVector del_W_f_p = forallOp<Real,NUMCOMPS>(del_W_f_p_c, "del_W_f_p_calc", del_W_f_p_calc, W_ave, W_ave_low_ahead);
	PVector del2_W_f  = forallOp<Real,NUMCOMPS>(del2_W_f_c, "del2_W_f_calc", del2_W_f_calc, W_ave, W_ave_high, W_ave_low_ahead);
	PVector W_ave_ahead = s_ahead_shift[d](W_ave);
	PVector W_ave_ahead2 = s_ahead_shift[d](W_ave_ahead);
	std::cout << "Here2" << std::endl;
	PVector W_ave_behind = s_behind_shift[d](W_ave);
	std::cout << "Here3" << std::endl;
	PVector W_ave_behind2 = s_behind_shift[d](W_ave_behind);
	PVector del2_W_c  = forallOp<Real,NUMCOMPS>(del2_W_c_c, "del2_W_c_calc", del2_W_c_calc, W_ave, W_ave_behind, W_ave_ahead);
	PVector del2_W_c_ahead = s_ahead_shift[d](del2_W_c);
	PVector del2_W_c_behind = s_behind_shift[d](del2_W_c);
	PVector del3_W = forallOp<Real,NUMCOMPS>(del3_W_c, "del3_W_calc", del3_W_calc, del2_W_c_ahead, del2_W_c);
	PVector del3_W_behind = s_behind_shift[d](del3_W);
	PVector del3_W_ahead = s_ahead_shift[d](del3_W);
	PVector del3_W_ahead2 = s_ahead_shift[d](del3_W_ahead);
	PVector del2_W_lim = forallOp<Real,NUMCOMPS>(del2_W_lim_c, "del2_W_lim_calc", 
												del2_W_lim_calc, del2_W_f, del2_W_c, del2_W_c_ahead, del2_W_c_behind);
    PVector rho_i = forallOp<Real,NUMCOMPS>(rho_i_c, "rho_i_calc", rho_i_calc, del2_W_f, del2_W_lim, W_ave, 
											W_ave_ahead, W_ave_ahead2, W_ave_behind, W_ave_behind2);
	PVector del3_W_min = forallOp<Real,NUMCOMPS>(del3_W_min_c, "del3_W_min_calc", del3_W_min_calc, del3_W_behind, 
												 del3_W, del3_W_ahead, del3_W_ahead2);	
    PVector del3_W_max = forallOp<Real,NUMCOMPS>(del3_W_max_c, "del3_W_max_calc", del3_W_max_calc, del3_W_behind, 
												 del3_W, del3_W_ahead, del3_W_ahead2);												 
    PVector W_ave_low_ahead_limited = forallOp<Real,NUMCOMPS>(limiter_low_c, "limiter_low_calc", limiter_low_calc, W_ave_low_ahead, 
															  del_W_f_m, del_W_f_p, W_ave, W_ave_ahead2, W_ave_behind2, rho_i, del3_W_max,
															  del3_W_min);															  
	PVector W_ave_high_limited = forallOp<Real,NUMCOMPS>(limiter_high_c, "limiter_high_calc", limiter_high_calc, W_ave_high, 
															  del_W_f_m, del_W_f_p, W_ave, W_ave_ahead2, W_ave_behind2, rho_i, del3_W_max,
															  del3_W_min);	
														  
	W_ave_low = s_behind_shift[d](W_ave_low_ahead_limited);
	W_ave_high = s_copy_f[d](W_ave_high_limited);
	
	CH_STOP(tint);
#if DIM>1
	PVector W_low = s_deconvolve_f[d](W_ave_low);		
	PVector W_high = s_deconvolve_f[d](W_ave_high);
#else
	PVector W_low = W_ave_low;
	PVector W_high = W_ave_high;
#endif	

	PScalar Lambda_f = forallOp<Real>(lambdacalcnum, "lambdacalc", lambdacalc, W_low, W_high, d, gamma);
	PVector F_low = forallOp<Real,NUMCOMPS>(getfluxnum, "getflux", getFlux, W_low, d, gamma);
	PVector F_high = forallOp<Real,NUMCOMPS>(getfluxnum, "getflux", getFlux, W_high, d, gamma);
	PVector F_f = forallOp<Real,NUMCOMPS>(rusanovnum, "rusanov",rusanovState, W_low, W_high, F_low, F_high, Lambda_f, d, gamma);
#if DIM>1				
	PVector F_ave_f = s_convolve_f[d](F_f);
#else
	PVector F_ave_f = F_f;
#endif	

      a_Rhs += s_divergence[d](F_ave_f);
  }

  a_Rhs *= -1./s_dx;
    
  return retval;
}

void 
MHDOp::
proto_step2(BoxData<Real,NUMCOMPS>& a_Rhs,
           const BoxData<Real,NUMCOMPS>& a_U,
           const Bx& a_rangeBox)
{

  a_Rhs.setVal(0.0);
  PVector divB_term;
  divB_term *= 0.0;
  Real gamma = s_gamma;
  unsigned long long int sqrtnum = 10;  //this one is just a guess
  unsigned long long int ctoprmnum  = 4*DIM + 5;
  unsigned long long int powellnum = sqrtnum +4 + DIM; // Check
  unsigned long long int bavgnum = sqrtnum +4 + DIM; // Check
  PVector W_bar = forallOp<Real,NUMCOMPS>(ctoprmnum, "consToPrim", consToPrim,a_U, gamma);
  PVector U = s_deconvolve(a_U);
  PVector W = forallOp<Real,NUMCOMPS>(ctoprmnum, "consToPrim", consToPrim,U, gamma);
  PVector W_ave = s_laplacian(W_bar);
  W_ave *= (1.0/24.0);
  W_ave += W;
  PVector Powell_term = forallOp<Real,NUMCOMPS>(powellnum, "Powell", Powell,W_ave);
  for (int d = 0; d < DIM; d++)
  {
	PVector B_ave = forallOp<Real,NUMCOMPS>(bavgnum, "bavgcalc", Bavgcalc, W_ave, d);
	divB_term += s_divergence[d](B_ave);
	divB_term *= Powell_term;
	a_Rhs += divB_term;
  }
  a_Rhs *= -1./s_dx;  
}

Real gatherMaxWave(Real maxwaveproc)
{
  Real maxwaveall = maxwaveproc;
#ifdef CH_MPI
  Real sendBuf = maxwaveall;
  int result = MPI_Allreduce(&sendBuf, &maxwaveall, 1, MPI_CH_REAL, MPI_MAX, Chombo_MPI::comm);
  if (result != MPI_SUCCESS)
  {
    MayDay::Error("communication error in gather");
  }
#endif   
  return maxwaveall;
}


Real 
MHDOp::
step(LevelBoxData<NUMCOMPS> & a_Rhs,
     LevelBoxData<NUMCOMPS> & a_U)
{
	 //std::cout << "Reached here 3" << std::endl;
  CH_TIME("MHDOp::step(leveldata)");
  static bool initCalled =false;
  DisjointBoxLayout grids = a_U.disjointBoxLayout();
  IntVect              gv = a_U.ghostVect();
  if(!initCalled)
  {
    CH_TIME("defining stuff");
    define(grids, gv);
    initCalled = true;
  }
  a_U.exchange(s_exchangeCopier);
  Real maxwaveproc = 0;
  {
    CH_TIME("step_no_gather");
    DataIterator dit = grids.dataIterator();
#pragma omp parallel for
    for(int ibox = 0; ibox < dit.size(); ibox++)
    {
      Box grid = grids[dit[ibox]];
      Bx  pgrid = ProtoCh::getProtoBox(grid);
      BoxData<Real, NUMCOMPS>& ubd   =   a_U[dit[ibox]];
      BoxData<Real, NUMCOMPS>& a_rhsbd = a_Rhs[dit[ibox]];

      Real maxwavegrid = proto_step(a_rhsbd, ubd, pgrid);
      maxwaveproc = std::max(maxwavegrid, maxwaveproc);
    }
  }
  Real maxwaveall;
  {
    CH_TIME("gatherMaxWaveSpeed");
    maxwaveall = gatherMaxWave(maxwaveproc);
  }
  return maxwaveall;
}

void 
MHDOp::
step2(LevelBoxData<NUMCOMPS> & a_Rhs,
     LevelBoxData<NUMCOMPS> & a_U)
{
	 //std::cout << "Reached here 3" << std::endl;
  CH_TIME("MHDOp::step(leveldata)");
  static bool initCalled =false;
  DisjointBoxLayout grids = a_U.disjointBoxLayout();
  IntVect              gv = a_U.ghostVect();
  if(!initCalled)
  {
    CH_TIME("defining stuff");
    define(grids, gv);
    initCalled = true;
  }
  a_U.exchange(s_exchangeCopier);
  {
    CH_TIME("step_no_gather");
    DataIterator dit = grids.dataIterator();
#pragma omp parallel for
    for(int ibox = 0; ibox < dit.size(); ibox++)
    {
      Box grid = grids[dit[ibox]];
      Bx  pgrid = ProtoCh::getProtoBox(grid);
      BoxData<Real, NUMCOMPS>& ubd   =   a_U[dit[ibox]];
      BoxData<Real, NUMCOMPS>& a_rhsbd = a_Rhs[dit[ibox]];
      proto_step2(a_rhsbd, ubd, pgrid);
    }
  }
}

Real 
MHDOp::
maxWave(LevelBoxData<NUMCOMPS> & a_U)
{
  static bool initCalled =false;
  DisjointBoxLayout grids = a_U.disjointBoxLayout();
  IntVect              gv = a_U.ghostVect();
  if(!initCalled)
  {
    define(grids, gv);
    initCalled = true;
  }
  a_U.exchange(s_exchangeCopier);
  Real maxwaveproc = 0;
  unsigned long long int ctoprmnum  = 4*DIM + 5;
  unsigned long long int sqrtnum = 10;  //this one is just a guess
  unsigned long long int wavespdnum = sqrtnum +3 + DIM;

  Real gamma = s_gamma;
  DataIterator dit = grids.dataIterator();
#pragma omp parallel for
  for(int ibox = 0; ibox < dit.size(); ibox++)
  {
    Box grid = grids[dit[ibox]];
    Bx  pgrid = ProtoCh::getProtoBox(grid);
    BoxData<Real, NUMCOMPS>& ubd =  a_U[dit[ibox]];
    PVector U = s_deconvolve(ubd);
    PVector W = forallOp<Real,NUMCOMPS>(ctoprmnum, "consToPrim",consToPrim,ubd, gamma);
    PScalar umax = forallOp<Real>(wavespdnum, "wavespeed",waveSpeedBound,pgrid,W, gamma);
    Real maxwavegrid = umax.absMax();
    maxwaveproc = std::max(maxwavegrid, maxwaveproc);
  }
  Real maxwaveall = gatherMaxWave(maxwaveproc);

  return maxwaveall;
}