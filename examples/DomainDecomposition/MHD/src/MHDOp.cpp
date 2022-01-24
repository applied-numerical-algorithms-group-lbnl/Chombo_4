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
#include "Chombo_ProtoInterface.H"
#define PI 3.141592653589793

using     ::Proto::Var;
using     ::Proto::Stencil;
using     ::Proto::BoxData;
using     ::Proto::Point;
using     ::Proto::Shift;
using     ::Proto::forall;
using     ::Proto::forall_p;
using     ::Proto::alias;
typedef   ::Proto::Var<Real,NUMCOMPS> State;

Real MHDOp::s_gamma = 1.6666666666666666666667;
Real MHDOp::s_dx = 1.0;
Stencil<Real> MHDOp::s_laplacian;
Stencil<Real> MHDOp::s_deconvolve;
Stencil<Real> MHDOp::s_copy;
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

typedef BoxData<Real,1,::Proto::MEMTYPE_DEFAULT,1,1> PScalar;
typedef BoxData<Real,NUMCOMPS,::Proto::MEMTYPE_DEFAULT,1,1> PVector;



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
                 const State& a_del2_W_c,
                 const State& a_del2_W_c_behind)
{
	for (int i=0; i< NUMCOMPS; i++){
		a_del3_W(i) = a_del2_W_c(i) - a_del2_W_c_behind(i);
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
        if (a_rho_i(i) < rhs){
			rhs_test = a_del3_W_max(i) - a_del3_W_min(i);
	        lhs_test = 0.01*std::max({std::abs(a_del3_W_min(i)), std::abs(a_del3_W_max(i))});
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
        if (a_rho_i(i) < rhs){
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


PROTO_KERNEL_START
void limiter_check_calcF(State& a_limiter_check,
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
	double rhs_test, lhs_test, rhs_test2, lhs_test2, rhs_test3, lhs_test3;
	
	for (int i=0; i< NUMCOMPS; i++){
        if (a_rho_i(i) < rhs){
			rhs_test = a_del3_W_max(i) - a_del3_W_min(i);
	        lhs_test = 0.1*std::max({std::abs(a_del3_W_min(i)), std::abs(a_del3_W_max(i))});
			if (lhs_test <= rhs_test){
				if ((a_del_W_f_m(i) * a_del_W_f_p(i) <= 0.0) || ((a_W_ave(i) - a_W_ave_behind2(i))*(a_W_ave_ahead2(i) - a_W_ave(i)) <= 0.0)){
					rhs_test2 = 2.0*std::abs(a_del_W_f_p(i));
					lhs_test2 = std::abs(a_del_W_f_m(i));
					rhs_test3 = 2.0*std::abs(a_del_W_f_m(i));
					lhs_test3 = std::abs(a_del_W_f_p(i));
					if (a_del_W_f_m(i) * a_del_W_f_p(i) < 0.0){
						a_limiter_check(i) = 1;
					} else if ((lhs_test2 >= rhs_test2) || (lhs_test3 >= rhs_test3)){
						a_limiter_check(i) = 1;
					} else {
						a_limiter_check(i) = 0;
					}
				} else {
					rhs_test2 = 2.0*std::abs(a_del_W_f_p(i));
					lhs_test2 = std::abs(a_del_W_f_m(i));
					rhs_test3 = 2.0*std::abs(a_del_W_f_m(i));
					lhs_test3 = std::abs(a_del_W_f_p(i));
					if ((lhs_test2 >= rhs_test2)||(lhs_test2 >= rhs_test2)){
						a_limiter_check(i) = 1;
					} else {
						a_limiter_check(i) = 0;
					}
				}		
			} else {
				a_limiter_check(i) = 0;
			}
		} else {
			a_limiter_check(i) = 0;
		}
		
	}
}
PROTO_KERNEL_END(limiter_check_calcF, limiter_check_calc)


PROTO_KERNEL_START
void eta_tilde_d_calcF(Var<Real,1>& a_eta_tilde_d,
                 const State& a_W_ave,
                 const State& a_W_ave_ahead,
                 const State& a_W_ave_ahead2,
                 const State& a_W_ave_behind,
                 const State& a_W_ave_behind2,
				 int a_d)
{   
    double z1=0.85, z0=0.75, delta = 0.33, si;
	double p, p_ahead, p_ahead2, p_behind, p_behind2;
#if DIM==1	
    p = a_W_ave(2) + a_W_ave(3)*a_W_ave(3)/8.0/PI;
    p_ahead = a_W_ave_ahead(2) + a_W_ave_ahead(3)*a_W_ave_ahead(3)/8.0/PI;
    p_ahead2 = a_W_ave_ahead2(2) + a_W_ave_ahead2(3)*a_W_ave_ahead2(3)/8.0/PI;
    p_behind = a_W_ave_behind(2) + a_W_ave_behind(3)*a_W_ave_behind(3)/8.0/PI;
    p_behind2 = a_W_ave_behind2(2) + a_W_ave_behind2(3)*a_W_ave_behind2(3)/8.0/PI;
#endif
#if DIM==2	
    p = a_W_ave(3) + (a_W_ave(4)*a_W_ave(4)+a_W_ave(5)*a_W_ave(5))/8.0/PI;
    p_ahead = a_W_ave_ahead(3) + (a_W_ave_ahead(4)*a_W_ave_ahead(4)+a_W_ave_ahead(5)*a_W_ave_ahead(5))/8.0/PI;
    p_ahead2 = a_W_ave_ahead2(3) + (a_W_ave_ahead2(4)*a_W_ave_ahead2(4)+a_W_ave_ahead2(5)*a_W_ave_ahead2(5))/8.0/PI;
    p_behind = a_W_ave_behind(3) + (a_W_ave_behind(4)*a_W_ave_behind(4)+a_W_ave_behind(5)*a_W_ave_behind(5))/8.0/PI;
    p_behind2 = a_W_ave_behind2(3) + (a_W_ave_behind2(4)*a_W_ave_behind2(4)+a_W_ave_behind2(5)*a_W_ave_behind2(5))/8.0/PI;
#endif
#if DIM==3
	p = a_W_ave(4) + (a_W_ave(5)*a_W_ave(5)+a_W_ave(6)*a_W_ave(6)+a_W_ave(7)*a_W_ave(7))/8.0/PI;
	p_ahead = a_W_ave_ahead(4) + (a_W_ave_ahead(5)*a_W_ave_ahead(5)+a_W_ave_ahead(6)*a_W_ave_ahead(6)+a_W_ave_ahead(7)*a_W_ave_ahead(7))/8.0/PI;
	p_ahead2 = a_W_ave_ahead2(4) + (a_W_ave_ahead2(5)*a_W_ave_ahead2(5)+a_W_ave_ahead2(6)*a_W_ave_ahead2(6)+a_W_ave_ahead2(7)*a_W_ave_ahead2(7))/8.0/PI;
	p_behind = a_W_ave_behind(4) + (a_W_ave_behind(5)*a_W_ave_behind(5)+a_W_ave_behind(6)*a_W_ave_behind(6)+a_W_ave_behind(7)*a_W_ave_behind(7))/8.0/PI;
	p_behind2 = a_W_ave_behind2(4) + (a_W_ave_behind2(5)*a_W_ave_behind2(5)+a_W_ave_behind2(6)*a_W_ave_behind2(6)+a_W_ave_behind2(7)*a_W_ave_behind2(7))/8.0/PI;
#endif
    double arg = std::abs(p_ahead-p_behind)/std::abs(p_ahead2-p_behind2);
	if (arg > z1){si = 0.0;} 
	if (arg < z0){si = 1.0;} 
	if ((arg<=z1) && (arg>=z0)){si = 1.0 - (arg-z0)/(z1-z0);}
	
	double lhs1 = a_W_ave_behind(1+a_d)-a_W_ave_ahead(1+a_d);
	double lhs2 = std::abs(p_ahead-p_behind)/std::min({p_ahead,p_behind});
	if ((lhs1 > 0.0) && (lhs2 > delta)){
		a_eta_tilde_d(0) = si;
	} else {
		a_eta_tilde_d(0) = 0.0;
	}
}
PROTO_KERNEL_END(eta_tilde_d_calcF, eta_tilde_d_calc) 

PROTO_KERNEL_START
void eta_d_calcF(Var<Real,1>& a_eta_d,
                 const Var<Real,1>& a_eta_tilde_d,
                 const Var<Real,1>& a_eta_tilde_d_ahead,
                 const Var<Real,1>& a_eta_tilde_d_behind)
{   

	a_eta_d(0) = std::min({a_eta_tilde_d(0),a_eta_tilde_d_ahead(0),a_eta_tilde_d_behind(0)});

}
PROTO_KERNEL_END(eta_d_calcF, eta_d_calc) 

PROTO_KERNEL_START
void eta_calcF(Var<Real,1>& a_eta,
                 const Var<Real,1>& a_eta_d,
                 const Var<Real,1>& a_eta_d_old)
{   

	a_eta(0) = std::min({a_eta_d(0),a_eta_d_old(0)});

}
PROTO_KERNEL_END(eta_calcF, eta_calc) 

PROTO_KERNEL_START
void Flat_calcF(State& a_flattened,
                 const State& a_not_flattened,
                 const State& a_W_ave,
                 const Var<Real,1>& a_eta)
{   
    
	for (int i=0; i< NUMCOMPS; i++){
		a_flattened(i) = (1.0-a_eta(0))*a_not_flattened(i) + a_eta(0)*a_W_ave(i);			
	}

}
PROTO_KERNEL_END(Flat_calcF, Flat_calc) 

PROTO_KERNEL_START
void v_d_calcF(Var<Real,1>& a_v_d,
                 const State& a_W_bar,
                 const int a_d)
{
		a_v_d(0) = a_W_bar(1+a_d);			
}
PROTO_KERNEL_END(v_d_calcF, v_d_calc) 

PROTO_KERNEL_START
void viscosity1_calcF(Var<Real,1>& a_viscosity,
                 const Var<Real,1>& a_v,
                 const Var<Real,1>& a_v_behind)
{
		a_viscosity(0) = a_v(0)-a_v_behind(0);			
}
PROTO_KERNEL_END(viscosity1_calcF, viscosity1_calc) 

PROTO_KERNEL_START
void v_d2_div_calcF(Var<Real,1>& a_v_d2_div,
                 const Var<Real,1>& v_d2_ahead,
                 const Var<Real,1>& v_d2_behind,
                 const Var<Real,1>& v_d2_behind_dp,
                 const Var<Real,1>& v_d2_behind_dm)
{
		a_v_d2_div(0) = (v_d2_ahead(0)-v_d2_behind(0)+v_d2_behind_dp(0)-v_d2_behind_dm(0))/4.0;			
}
PROTO_KERNEL_END(v_d2_div_calcF, v_d2_div_calc) 


PROTO_KERNEL_START
void Fast_MS_speed_calcF(Var<Real,1>& a_Fast_MS_speed,
                 const State& a_W_bar,
				 int a_d,
                 Real a_gamma)
{
    Real gamma = a_gamma;
	Real rho=0., p=0., Bx=0., By=0., Bz=0., ce, B_mag, Bdir;
#if DIM == 1
	rho = a_W_bar(0);
	p   = a_W_bar(2);
    Bx  = a_W_bar(3);
#endif	
#if DIM == 2
	rho = a_W_bar(0);
	p   = a_W_bar(3);
    Bx  = a_W_bar(4);
    By  = a_W_bar(5);
#endif	
#if DIM == 3
	rho = a_W_bar(0);
	p   = a_W_bar(4);
    Bx  = a_W_bar(5);
    By  = a_W_bar(6);
    Bz  = a_W_bar(7);
#endif	
    if (a_d == 0){
		Bdir = Bx;
	};
	if (a_d == 1){
		Bdir = By;
	};
	if (a_d == 2){
		Bdir = Bz;
	};
	
	ce = sqrt(gamma*p/rho);
	B_mag = sqrt(Bx*Bx+By*By+Bz*Bz);
	a_Fast_MS_speed(0) = 0.5*(sqrt((ce*ce)+( B_mag*B_mag/(4.0*PI*rho) )+( abs(Bdir)*ce/sqrt(PI*rho) ))+
		  sqrt((ce*ce)+( B_mag*B_mag/(4.0*PI*rho) )-( abs(Bdir)*ce/sqrt(PI*rho) )));
}
PROTO_KERNEL_END(Fast_MS_speed_calcF, Fast_MS_speed_calc)

PROTO_KERNEL_START
void Fast_MS_speed_min_calcF(Var<Real,1>& a_Fast_MS_speed_min,
                 const Var<Real,1>& a_Fast_MS_speed,
                 const Var<Real,1>& a_Fast_MS_speed_behind)
{   
	a_Fast_MS_speed_min(0) = std::min({a_Fast_MS_speed(0),a_Fast_MS_speed_behind(0)});
}
PROTO_KERNEL_END(Fast_MS_speed_min_calcF, Fast_MS_speed_min_calc) 


PROTO_KERNEL_START
void Visc_coef_calcF(Var<Real,1>& a_Visc_coef,
                 const Var<Real,1>& a_h_lambda,
                 const Var<Real,1>& a_Fast_MS_speed_min)
{   
	if (a_h_lambda(0) <= 0){
		double temp = a_h_lambda(0)*a_h_lambda(0)/0.3/a_Fast_MS_speed_min(0)/a_Fast_MS_speed_min(0);
		a_Visc_coef(0) = a_h_lambda(0)*std::min({temp,1.0});
	} else {
		a_Visc_coef(0) = 0.0;
	}
}
PROTO_KERNEL_END(Visc_coef_calcF, Visc_coef_calc) 

PROTO_KERNEL_START
void mu_f_calcF(State& a_mu_f,
                 const Var<Real,1>& a_Visc_coef,
                 const State& a_U,
                 const State& a_U_behind)
{   
	for (int i=0; i< NUMCOMPS; i++){
		a_mu_f(i) = 0.3*a_Visc_coef(0)*(a_U(i)-a_U_behind(i));
	}
}
PROTO_KERNEL_END(mu_f_calcF, mu_f_calc) 


void
MHDOp::
define(const DisjointBoxLayout& a_grids,
       const IntVect          & a_ghostVect) 
{
  CH_TIME("MHDOp::define");
  s_laplacian = Stencil<Real>::Laplacian();
  s_deconvolve = (-1.0/24.0)*s_laplacian + (1.0)*Shift(Point::Zeros());
  s_copy = 1.0*Shift(Point::Zeros());
  for (int dir = 0; dir < DIM; dir++)
  {
    s_laplacian_f[dir] = Stencil<Real>::LaplacianFace(dir);
    s_deconvolve_f[dir] = (-1.0/24.0)*s_laplacian_f[dir] + 1.0*Shift(Point::Zeros());
    s_convolve_f[dir] = (1.0/24.0)*s_laplacian_f[dir] + 1.0*Shift(Point::Zeros());
    s_interp_H[dir] = Stencil<Real>::CellToEdgeH(dir);
    s_interp_L[dir] = Stencil<Real>::CellToEdgeL(dir);
    s_divergence[dir] = Stencil<Real>::FluxDivergence(dir);
	s_ahead_shift[dir] = 1.0*Shift(Point::Basis(dir)*(1));
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

  PVector W_bar = forall<Real,NUMCOMPS>( consToPrim,a_U, gamma);
  PVector U = s_deconvolve(a_U);
  PVector W = forall<Real,NUMCOMPS>( consToPrim,U, gamma);
  PScalar umax = forall<Real>(waveSpeedBound,a_rangeBox,W, gamma);
  retval = umax.absMax();
  PVector W_ave = s_laplacian(W_bar);
  W_ave *= (1.0/24.0);
  W_ave += W;
  PScalar eta ;
  PScalar eta_old ;
  
  for (int d = 0; d < DIM; d++)
  {

	PVector W_ave_ahead = alias(W_ave,Point::Basis(d)*(-1));
	PVector W_ave_ahead2 = alias(W_ave,Point::Basis(d)*(-2));
	PVector W_ave_ahead3 = alias(W_ave,Point::Basis(d)*(-3));
	PVector W_ave_behind = alias(W_ave, Point::Basis(d)*(1));
	PVector W_ave_behind2 = alias(W_ave,Point::Basis(d)*(2));	
	PVector W_ave_behind3 = alias(W_ave,Point::Basis(d)*(3));
	
	PScalar eta_tilde_d = forall<Real>( eta_tilde_d_calc, W_ave, W_ave_ahead, W_ave_ahead2,
	                              W_ave_behind, W_ave_behind2, d);
	PScalar eta_tilde_d_ahead = alias(eta_tilde_d,Point::Basis(d)*(-1));
	PScalar eta_tilde_d_behind = alias(eta_tilde_d,Point::Basis(d)*(1));	
	
	PScalar eta_d = forall<Real>( eta_d_calc, eta_tilde_d, eta_tilde_d_ahead, eta_tilde_d_behind);
	if (d>0){
		eta = forall<Real>(eta_calc, eta_d, eta_old);
	} else {
		eta = s_copy(eta_d);
	}	
	eta_old = s_copy(eta);
  }
  
  for (int d = 0; d < DIM; d++)
  {
    CH_START(tint);
    PVector W_ave_low = s_interp_L[d](W_ave);
    PVector W_ave_high = s_interp_H[d](W_ave);
	
	// Limiter Starts here
	PVector W_ave_low_ahead = alias(W_ave_low,Point::Basis(d)*(-1));
	PVector del_W_f_m = forall<Real,NUMCOMPS>( del_W_f_m_calc, W_ave, W_ave_high);
	PVector del_W_f_p = forall<Real,NUMCOMPS>(del_W_f_p_calc, W_ave, W_ave_low_ahead);
	PVector del2_W_f  = forall<Real,NUMCOMPS>(del2_W_f_calc, W_ave, W_ave_high, W_ave_low_ahead);
	PVector W_ave_ahead = alias(W_ave,Point::Basis(d)*(-1));
	PVector W_ave_ahead2 = alias(W_ave,Point::Basis(d)*(-2));
	PVector W_ave_behind = alias(W_ave,Point::Basis(d)*(1));
	PVector W_ave_behind2 = alias(W_ave,Point::Basis(d)*(2));
	PVector del2_W_c  = forall<Real,NUMCOMPS>( del2_W_c_calc, W_ave, W_ave_behind, W_ave_ahead);
	PVector del2_W_c_ahead = alias(del2_W_c,Point::Basis(d)*(-1));
	PVector del2_W_c_behind = alias(del2_W_c,Point::Basis(d)*(1));
	PVector del3_W = forall<Real,NUMCOMPS>( del3_W_calc, del2_W_c, del2_W_c_behind);
	PVector del3_W_behind = alias(del3_W,Point::Basis(d)*(1));
	PVector del3_W_ahead = alias(del3_W,Point::Basis(d)*(-1));
	PVector del3_W_ahead2 = alias(del3_W,Point::Basis(d)*(-2));
	PVector del2_W_lim = forall<Real,NUMCOMPS>(del2_W_lim_calc, del2_W_f, del2_W_c, del2_W_c_ahead, del2_W_c_behind);
    PVector rho_i = forall<Real,NUMCOMPS>( rho_i_calc, del2_W_f, del2_W_lim, W_ave, 
											W_ave_ahead, W_ave_ahead2, W_ave_behind, W_ave_behind2);
	PVector del3_W_min = forall<Real,NUMCOMPS>(del3_W_min_calc, del3_W_behind, 
												 del3_W, del3_W_ahead, del3_W_ahead2);	
    PVector del3_W_max = forall<Real,NUMCOMPS>( del3_W_max_calc, del3_W_behind, 
												 del3_W, del3_W_ahead, del3_W_ahead2);												 
    PVector W_ave_low_ahead_limited = forall<Real,NUMCOMPS>( limiter_low_calc, W_ave_low_ahead, del_W_f_m, del_W_f_p, W_ave, W_ave_ahead2, W_ave_behind2, rho_i, del3_W_max, del3_W_min);															  
	PVector W_ave_high_limited = forall<Real,NUMCOMPS>( limiter_high_calc, W_ave_high, del_W_f_m, del_W_f_p, W_ave, W_ave_ahead2, W_ave_behind2, rho_i, del3_W_max, del3_W_min);	
    //PVector limiter_check = forall<Real,NUMCOMPS>(limiter_check_calc, del_W_f_m, del_W_f_p, W_ave, W_ave_ahead2, W_ave_behind2, rho_i, del3_W_max, del3_W_min);		
    // Slope flattening starts here
	PVector W_ave_low_ahead_lim_flat = forall<Real,NUMCOMPS>(Flat_calc, W_ave_low_ahead_limited, W_ave, eta);
	PVector W_ave_high_lim_flat = forall<Real,NUMCOMPS>(Flat_calc, W_ave_high_limited, W_ave, eta);
	//Slope flattening ends here

	W_ave_low = alias(W_ave_low_ahead_lim_flat,Point::Basis(d)*(1));
	W_ave_high = s_copy_f[d](W_ave_high_lim_flat);
	// Limiter ends here
	
	CH_STOP(tint);
#if DIM>1
	PVector W_low = s_deconvolve_f[d](W_ave_low);		
	PVector W_high = s_deconvolve_f[d](W_ave_high);
#else
	PVector W_low = W_ave_low;
	PVector W_high = W_ave_high;
#endif	

	PScalar Lambda_f = forall<Real>(lambdacalc, W_low, W_high, d, gamma);
	PVector F_low = forall<Real,NUMCOMPS>(getFlux, W_low, d, gamma);
	PVector F_high = forall<Real,NUMCOMPS>(getFlux, W_high, d, gamma);
	PVector F_f = forall<Real,NUMCOMPS>(rusanovState, W_low, W_high, F_low, F_high, Lambda_f, d, gamma);
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

// Used to find the Powell term
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

  PVector W_bar = forall<Real,NUMCOMPS>(consToPrim,a_U, gamma);
  PVector U = s_deconvolve(a_U);
  PVector W = forall<Real,NUMCOMPS>( consToPrim,U, gamma);
  PVector W_ave = s_laplacian(W_bar);
  W_ave *= (1.0/24.0);
  W_ave += W;
  PVector Powell_term = forall<Real,NUMCOMPS>( Powell,W_ave);
  for (int d = 0; d < DIM; d++)
  {
	PVector B_ave = forall<Real,NUMCOMPS>( Bavgcalc, W_ave, d);
	divB_term += s_divergence[d](B_ave);
	divB_term *= Powell_term;
	a_Rhs += divB_term;
  }
  a_Rhs *= -1./s_dx;  
}

// Used to implement artificial viscosity
void 
MHDOp::
proto_step3(BoxData<Real,NUMCOMPS>& a_Rhs,
           const BoxData<Real,NUMCOMPS>& a_U,
           const Bx& a_rangeBox)
{

  a_Rhs.setVal(0.0);
  Real gamma = s_gamma;
  PVector W_bar = forall<Real,NUMCOMPS>( consToPrim,a_U, gamma);
  for (int d = 0; d < DIM; d++)
  { 
    PScalar v_d =  forall<Real>(v_d_calc,W_bar,d);
	PScalar v_d_behind = alias(v_d,Point::Basis(d)*(1));
	PScalar h_lambda = forall<Real>(viscosity1_calc,v_d,v_d_behind); 
	for (int d2 = 0; d2 < DIM; d2++){	
		if (d2!=d){
			PScalar v_d2 = forall<Real>(v_d_calc,W_bar,d2); 
			PScalar v_d2_ahead = alias(v_d2,Point::Basis(d2)*(-1));
			PScalar v_d2_behind = alias(v_d2,Point::Basis(d2)*(1));
			PScalar v_d2_behind_dp = alias(v_d2_ahead,Point::Basis(d)*(1));
			PScalar v_d2_behind_dm = alias(v_d2_behind,Point::Basis(d)*(1));
			PScalar v_d2_div = forall<Real>(v_d2_div_calc,v_d2_ahead,v_d2_behind,v_d2_behind_dp,v_d2_behind_dm);
			h_lambda += v_d2_div;
		}
	}
	PScalar Fast_MS_speed = forall<Real>(Fast_MS_speed_calc, W_bar, d, gamma);
	PScalar Fast_MS_speed_behind = alias(Fast_MS_speed,Point::Basis(d)*(1));
	PScalar Fast_MS_speed_min = forall<Real>(Fast_MS_speed_min_calc,Fast_MS_speed,Fast_MS_speed_behind);
	PScalar Visc_coef = forall<Real>(Visc_coef_calc,h_lambda,Fast_MS_speed_min);
	PVector a_U_behind = s_behind_shift[d](a_U);
	PVector mu_f = forall<Real,NUMCOMPS>(mu_f_calc, Visc_coef, a_U, a_U_behind);
	a_Rhs += s_divergence[d](mu_f);
  }
  a_Rhs *= -1./s_dx;  
}


// This is used to print any of the values used to find RHS (used for debugging). Should be identical to proto_step in everything except RHS.
void 
MHDOp::
proto_step_test(BoxData<Real,NUMCOMPS>& a_Rhs,
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

  PVector W_bar = forall<Real,NUMCOMPS>( consToPrim,a_U, gamma);
  PVector U = s_deconvolve(a_U);
  PVector W = forall<Real,NUMCOMPS>( consToPrim,U, gamma);
  PScalar umax = forall<Real>(waveSpeedBound,a_rangeBox,W, gamma);
  retval = umax.absMax();
  PVector W_ave = s_laplacian(W_bar);
  W_ave *= (1.0/24.0);
  W_ave += W;
  PScalar eta ;
  PScalar eta_old ;
  
  for (int d = 0; d < DIM; d++)
  {

	PVector W_ave_ahead = alias(W_ave,Point::Basis(d)*(-1));
	PVector W_ave_ahead2 = alias(W_ave,Point::Basis(d)*(-2));
	PVector W_ave_ahead3 = alias(W_ave,Point::Basis(d)*(-3));
	PVector W_ave_behind = alias(W_ave, Point::Basis(d)*(1));
	PVector W_ave_behind2 = alias(W_ave,Point::Basis(d)*(2));	
	PVector W_ave_behind3 = alias(W_ave,Point::Basis(d)*(3));
	
	PScalar eta_tilde_d = forall<Real>( eta_tilde_d_calc, W_ave, W_ave_ahead, W_ave_ahead2,
	                              W_ave_behind, W_ave_behind2, d);
	PScalar eta_tilde_d_ahead = alias(eta_tilde_d,Point::Basis(d)*(-1));
	PScalar eta_tilde_d_behind = alias(eta_tilde_d,Point::Basis(d)*(1));	
	
	PScalar eta_d = forall<Real>( eta_d_calc, eta_tilde_d, eta_tilde_d_ahead, eta_tilde_d_behind);
	if (d>0){
		eta = forall<Real>(eta_calc, eta_d, eta_old);
	} else {
		eta = s_copy(eta_d);
	}	
	eta_old = s_copy(eta);
  }
  
  for (int d = 0; d < DIM; d++)
  {
    CH_START(tint);
    PVector W_ave_low = s_interp_L[d](W_ave);
    PVector W_ave_high = s_interp_H[d](W_ave);
	
	// Limiter Starts here
	PVector W_ave_low_ahead = alias(W_ave_low,Point::Basis(d)*(-1));
	PVector del_W_f_m = forall<Real,NUMCOMPS>( del_W_f_m_calc, W_ave, W_ave_high);
	PVector del_W_f_p = forall<Real,NUMCOMPS>(del_W_f_p_calc, W_ave, W_ave_low_ahead);
	PVector del2_W_f  = forall<Real,NUMCOMPS>(del2_W_f_calc, W_ave, W_ave_high, W_ave_low_ahead);
	PVector W_ave_ahead = alias(W_ave,Point::Basis(d)*(-1));
	PVector W_ave_ahead2 = alias(W_ave,Point::Basis(d)*(-2));
	PVector W_ave_behind = alias(W_ave,Point::Basis(d)*(1));
	PVector W_ave_behind2 = alias(W_ave,Point::Basis(d)*(2));
	PVector del2_W_c  = forall<Real,NUMCOMPS>( del2_W_c_calc, W_ave, W_ave_behind, W_ave_ahead);
	PVector del2_W_c_ahead = alias(del2_W_c,Point::Basis(d)*(-1));
	PVector del2_W_c_behind = alias(del2_W_c,Point::Basis(d)*(1));
	PVector del3_W = forall<Real,NUMCOMPS>( del3_W_calc, del2_W_c, del2_W_c_behind);
	PVector del3_W_behind = alias(del3_W,Point::Basis(d)*(1));
	PVector del3_W_ahead = alias(del3_W,Point::Basis(d)*(-1));
	PVector del3_W_ahead2 = alias(del3_W,Point::Basis(d)*(-2));
	PVector del2_W_lim = forall<Real,NUMCOMPS>(del2_W_lim_calc, del2_W_f, del2_W_c, del2_W_c_ahead, del2_W_c_behind);
    PVector rho_i = forall<Real,NUMCOMPS>( rho_i_calc, del2_W_f, del2_W_lim, W_ave, 
											W_ave_ahead, W_ave_ahead2, W_ave_behind, W_ave_behind2);
	PVector del3_W_min = forall<Real,NUMCOMPS>(del3_W_min_calc, del3_W_behind, 
												 del3_W, del3_W_ahead, del3_W_ahead2);	
    PVector del3_W_max = forall<Real,NUMCOMPS>( del3_W_max_calc, del3_W_behind, 
												 del3_W, del3_W_ahead, del3_W_ahead2);												 
    PVector W_ave_low_ahead_limited = forall<Real,NUMCOMPS>( limiter_low_calc, W_ave_low_ahead, del_W_f_m, del_W_f_p, W_ave, W_ave_ahead2, W_ave_behind2, rho_i, del3_W_max, del3_W_min);															  
	PVector W_ave_high_limited = forall<Real,NUMCOMPS>( limiter_high_calc, W_ave_high, del_W_f_m, del_W_f_p, W_ave, W_ave_ahead2, W_ave_behind2, rho_i, del3_W_max, del3_W_min);	
    PVector limiter_check = forall<Real,NUMCOMPS>(limiter_check_calc, del_W_f_m, del_W_f_p, W_ave, W_ave_ahead2, W_ave_behind2, rho_i, del3_W_max, del3_W_min);		
    // Slope flattening starts here
	PVector W_ave_low_ahead_lim_flat = forall<Real,NUMCOMPS>(Flat_calc, W_ave_low_ahead_limited, W_ave, eta);
	PVector W_ave_high_lim_flat = forall<Real,NUMCOMPS>(Flat_calc, W_ave_high_limited, W_ave, eta);
	//Slope flattening ends here

	W_ave_low = alias(W_ave_low_ahead_lim_flat,Point::Basis(d)*(1));
	W_ave_high = s_copy_f[d](W_ave_high_lim_flat);
	// Limiter ends here
	
	CH_STOP(tint);
#if DIM>1
	PVector W_low = s_deconvolve_f[d](W_ave_low);		
	PVector W_high = s_deconvolve_f[d](W_ave_high);
#else
	PVector W_low = W_ave_low;
	PVector W_high = W_ave_high;
#endif	

	PScalar Lambda_f = forall<Real>(lambdacalc, W_low, W_high, d, gamma);
	PVector F_low = forall<Real,NUMCOMPS>(getFlux, W_low, d, gamma);
	PVector F_high = forall<Real,NUMCOMPS>(getFlux, W_high, d, gamma);
	PVector F_f = forall<Real,NUMCOMPS>(rusanovState, W_low, W_high, F_low, F_high, Lambda_f, d, gamma);
#if DIM>1				
	PVector F_ave_f = s_convolve_f[d](F_f);
#else
	PVector F_ave_f = F_f;
#endif	

    a_Rhs += s_copy_f[d](limiter_check);
  }

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

void 
MHDOp::
step3(LevelBoxData<NUMCOMPS> & a_Rhs,
     LevelBoxData<NUMCOMPS> & a_U)
{
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
      proto_step3(a_rhsbd, ubd, pgrid);
    }
  }
}


void 
MHDOp::
step_test(LevelBoxData<NUMCOMPS> & a_Rhs,
     LevelBoxData<NUMCOMPS> & a_U)
{
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
      proto_step_test(a_rhsbd, ubd, pgrid);
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


  Real gamma = s_gamma;
  DataIterator dit = grids.dataIterator();
#pragma omp parallel for
  for(int ibox = 0; ibox < dit.size(); ibox++)
  {
    Box grid = grids[dit[ibox]];
    Bx  pgrid = ProtoCh::getProtoBox(grid);
    BoxData<Real, NUMCOMPS>& ubd =  a_U[dit[ibox]];
    PVector U = s_deconvolve(ubd);
    PVector W = forall<Real,NUMCOMPS>(consToPrim,ubd, gamma);
    PScalar umax = forall<Real>(waveSpeedBound,pgrid,W, gamma);
    Real maxwavegrid = umax.absMax();
    maxwaveproc = std::max(maxwavegrid, maxwaveproc);
  }
  Real maxwaveall = gatherMaxWave(maxwaveproc);

  return maxwaveall;
}
