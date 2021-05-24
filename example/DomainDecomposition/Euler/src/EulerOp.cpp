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
#include "EulerOp.H"
#include "Chombo_ProtoInterface.H"
using     ::Proto::Var;
using     ::Proto::Stencil;
using     ::Proto::BoxData;
using     ::Proto::Point;
using     ::Proto::Shift;
using     ::Proto::forall;
using     ::Proto::forall_p;
typedef   ::Proto::Var<Real,NUMCOMPS> State;

Real EulerOp::s_gamma = 1.4;
Real EulerOp::s_dx = 1.0;
Stencil<Real> EulerOp::s_laplacian;
Stencil<Real> EulerOp::s_deconvolve;
Stencil<Real> EulerOp::s_laplacian_f[DIM];
Stencil<Real> EulerOp::s_deconvolve_f[DIM];
Stencil<Real> EulerOp::s_interp_H[DIM];
Stencil<Real> EulerOp::s_interp_L[DIM];
Stencil<Real> EulerOp::s_divergence[DIM];
Chombo4::Copier EulerOp::s_exchangeCopier;

typedef BoxData<Real,1,MEMTYPE_DEFAULT,1,1> PScalar;
typedef BoxData<Real,NUMCOMPS,MEMTYPE_DEFAULT,1,1> PVector;
using ::Proto::forall;
using ::Proto::forallOp;

PROTO_KERNEL_START
void 
consToPrimF(State&         a_W, 
            const State&   a_U,
            Real         a_gamma)
{
  Real rho = a_U(0);
  Real v2 = 0.0;
  Real gamma = a_gamma;
  a_W(0) = rho;
    
  for (int i = 1; i <= DIM; i++)
  {
    Real v;
    v = a_U(i) / rho;
        
    a_W(i) = v;
    v2 += v*v;
  }
    
  a_W(NUMCOMPS-1) = (a_U(NUMCOMPS-1) - .5 * rho * v2) * (gamma - 1.0);
}
PROTO_KERNEL_END(consToPrimF, consToPrim)

PROTO_KERNEL_START
void upwindStateF(State& a_out,
                  const State& a_low,
                  const State& a_high,
                  int   a_dir,
                  Real a_gamma)
{
  const Real& rhol = a_low(0);
  const Real& rhor = a_high(0);
  const Real& ul = a_low(a_dir+1);
  const Real& ur = a_high(a_dir+1);
  const Real& pl = a_low(NUMCOMPS-1);
  const Real& pr = a_high(NUMCOMPS-1);
  Real gamma = a_gamma;
  Real rhobar = (rhol + rhor)*.5;
  Real pbar = (pl + pr)*.5;
  Real ubar = (ul + ur)*.5;
  Real cbar = sqrt(gamma*pbar/rhobar);
  Real pstar = (pl + pr)*.5 + rhobar*cbar*(ul - ur)*.5;
  Real ustar = (ul + ur)*.5 + (pl - pr)/(2*rhobar*cbar);
  int sign;
  if (ustar > 0) 
  {
    sign = -1;
    for (int icomp = 0;icomp < NUMCOMPS;icomp++)
    {
      a_out(icomp) = a_low(icomp);
    }
  }
  else
  {
    sign = 1;
    for (int icomp = 0;icomp < NUMCOMPS;icomp++)
    {
      a_out(icomp) = a_high(icomp);
    }
  }
  if (cbar + sign*ubar > 0)
  {
    a_out(0) += (pstar - a_out(NUMCOMPS-1))/(cbar*cbar);
    a_out(a_dir+1) = ustar;
    a_out(NUMCOMPS-1) = pstar;
  }
}
PROTO_KERNEL_END(upwindStateF, upwindState)

PROTO_KERNEL_START
void getFluxF(State& a_F, const State& a_W, 
              int    a_d,
              Real a_gamma)
{
  Real F0 = a_W(a_d+1)*a_W(0);
  Real W2 = 0.0;
  Real gamma = a_gamma;

  a_F(0) = F0;

  for (int d = 1; d <= DIM; d++)
  {
    Real Wd = a_W(d);

    a_F(d) = Wd*F0;
    W2 += Wd*Wd;
  }

  a_F(a_d+1) += a_W(NUMCOMPS-1);
  a_F(NUMCOMPS-1) = gamma/(gamma - 1.0) * a_W(a_d+1) * a_W(NUMCOMPS-1) + 0.5 * F0 * W2;
}
PROTO_KERNEL_END(getFluxF, getFlux)

PROTO_KERNEL_START
void waveSpeedBoundF(Var<Real,1>& a_speed,
                     const State& a_W,
                     Real       a_gamma)
{
  Real gamma = a_gamma;
  a_speed(0) = DIM*sqrt(gamma*a_W(NUMCOMPS-1)/a_W(0));
  for (int dir = 1 ; dir <= DIM; dir++)
  {
    a_speed(0) += a_W(dir);
  }
}
PROTO_KERNEL_END(waveSpeedBoundF, waveSpeedBound)

#include "Chombo_NamespaceHeader.H"

void
EulerOp::
define(const DisjointBoxLayout& a_grids,
       const IntVect          & a_ghostVect) 
{
  CH_TIME("EulerOp::define");
  s_laplacian = Stencil<Real>::Laplacian();
  s_deconvolve = (-1.0/24.0)*s_laplacian + (1.0)*Shift(Point::Zeros());
  for (int dir = 0; dir < DIM; dir++)
  {
    s_laplacian_f[dir] = Stencil<Real>::LaplacianFace(dir);
    s_deconvolve_f[dir] = (-1.0/24.0)*s_laplacian_f[dir] + 1.0*Shift(Point::Zeros());
    s_interp_H[dir] = Stencil<Real>::CellToEdgeH(dir);
    s_interp_L[dir] = Stencil<Real>::CellToEdgeL(dir);
    s_divergence[dir] = Stencil<Real>::FluxDivergence(dir);
  }
  s_exchangeCopier.exchangeDefine(a_grids, a_ghostVect);
}


void
EulerOp::
proto_step(BoxData<Real,NUMCOMPS>& a_Rhs,
           const BoxData<Real,NUMCOMPS>& a_U,
           const Bx& a_rangeBox,
           Reduction<Real,Op::Abs>& a_Rxn)
{

  CH_TIMERS("EulerOp::step(boxdata)");
  CH_TIMER("interp stencil eval", tint);
  CH_TIMER("riemann problem", trie);
  CH_TIMER("get_flux", tgf);
  CH_TIMER("deconvolve", tdcv);
  CH_TIMER("get_flux2", tgf2);
  CH_TIMER("laplacian_arith", tlap);
  CH_TIMER("divergence_eval", tdiv);
  CH_TIMER("divide_by_dx", tdx);
  a_Rhs.setVal(0.0);

  Real gamma = s_gamma;
  unsigned long long int sqrtnum = 10;  //this one is just a guess
  unsigned long long int ctoprmnum  = 4*DIM + 5;
  unsigned long long int upwindnum  = sqrtnum + 25;
  unsigned long long int getfluxnum = 9 + DIM*3;
  unsigned long long int wavespdnum = sqrtnum +3 + DIM;

  PVector W_bar = forallOp<Real,NUMCOMPS>(ctoprmnum, "consToPrim", consToPrim,a_U, gamma);
  PVector U = s_deconvolve(a_U);
  PVector W = forallOp<Real,NUMCOMPS>(ctoprmnum, "consToPrim", consToPrim,U, gamma);
  PScalar umax = forallOp<Real>(wavespdnum, "wavespeed",waveSpeedBound,a_rangeBox,W, gamma);
  umax.absMax(a_Rxn);
  PVector W_ave = s_laplacian(W_bar,1.0/24.0);
  W_ave += W;
  for (int d = 0; d < DIM; d++)
  {
    CH_START(tint);
    PVector W_ave_low = s_interp_L[d](W_ave);
    PVector W_ave_high = s_interp_H[d](W_ave);
    CH_STOP(tint);
    CH_START(trie);
    PVector W_ave_f = forallOp<Real,NUMCOMPS>(upwindnum, "upwind",
      upwindState,W_ave_low, W_ave_high,d,  gamma);
    CH_STOP(trie);
    CH_START(tgf);
#if DIM>1
    PVector F_bar_f = forallOp<Real,NUMCOMPS>(getfluxnum, "getflux",getFlux, W_ave_f, d,  gamma);
#endif
    CH_STOP(tgf);
    CH_START(tdcv);
#if DIM>1
    PVector W_f = s_deconvolve_f[d](W_ave_f);
#else
    PVector W_f = W_ave_f;
#endif
    CH_STOP(tdcv);
    CH_START(tgf2);
    PVector F_ave_f = forallOp<Real,NUMCOMPS>(getfluxnum, "getflux",getFlux, W_f, d, gamma);
    CH_STOP(tgf2);
    CH_START(tlap);
#if DIM>1
    F_bar_f *= (1./24.);
    F_ave_f += s_laplacian_f[d](F_bar_f,1.0/24.0);
#endif
    CH_STOP(tlap);
    CH_START(tdiv);
    a_Rhs += s_divergence[d](F_ave_f);
    CH_STOP(tdiv);
  }
  CH_START(tdx);
  a_Rhs *= -1./s_dx;
  CH_STOP(tdx);
}

Real gatherMaxWave(Real maxwaveproc)
{
  Real maxwaveall = maxwaveproc;
#ifdef CH_MPI
  Real sendBuf = maxwaveall;
  int result = MPI_Allreduce(&sendBuf, &maxwaveall, 1, MPI_CH_REAL, MPI_MAX, CH4_SPMD::Chombo_MPI::comm);

  if (result != MPI_SUCCESS)
  {
    MayDay::Error("communication error in gather");
  }
#endif  
  
  return maxwaveall;
}

void 
EulerOp::
step(LevelBoxData<NUMCOMPS> & a_Rhs,
     LevelBoxData<NUMCOMPS> & a_U,
     Reduction<Real,Op::Abs> & a_Rxn)
{
  CH_TIME("EulerOp::step(leveldata)");
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

    for(int ibox = 0; ibox < dit.size(); ibox++)
    {
      Box grid = grids[dit[ibox]];
      Bx  pgrid = ProtoCh::getProtoBox(grid);
      BoxData<Real, NUMCOMPS>& ubd   =   a_U[dit[ibox]];
      BoxData<Real, NUMCOMPS>& rhsbd = a_Rhs[dit[ibox]];
      proto_step(rhsbd, ubd, pgrid, a_Rxn);
    }
  }
}


Real 
EulerOp::
maxWave(LevelBoxData<NUMCOMPS> & a_U,
        Reduction<Real,Op::Abs> & a_Rxn)
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
  unsigned long long int ctoprmnum  = 4*DIM + 5;
  unsigned long long int sqrtnum = 10;  //this one is just a guess
  unsigned long long int wavespdnum = sqrtnum +3 + DIM;

  Real gamma = s_gamma;
  DataIterator dit = grids.dataIterator();

  for(int ibox = 0; ibox < dit.size(); ibox++)
  {
    Box grid = grids[dit[ibox]];
    Bx  pgrid = ProtoCh::getProtoBox(grid);
    BoxData<Real, NUMCOMPS>& ubd =  a_U[dit[ibox]];
    PVector U = s_deconvolve(ubd);
    PVector W = forallOp<Real,NUMCOMPS>(ctoprmnum, "consToPrim",consToPrim,ubd, gamma);
    PScalar umax = forallOp<Real>(wavespdnum, "wavespeed",waveSpeedBound,pgrid,W, gamma);
    umax.absMax(a_Rxn);
  }
  Real maxwaveall = gatherMaxWave(a_Rxn.fetch());

  return maxwaveall;
}
#include "Chombo_NamespaceFooter.H"
