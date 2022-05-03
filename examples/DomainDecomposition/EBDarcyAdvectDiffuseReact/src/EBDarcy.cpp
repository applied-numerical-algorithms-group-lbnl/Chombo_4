
#include "EBDarcy.H"
#include "EBParabolicIntegrators.H"
#include "Chombo_ParmParse.H"
#include "Chombo_NamespaceHeader.H"
#include "DebugFunctions.H"
using Proto::Var;
/*******/
EBDarcy::
EBDarcy(shared_ptr<EBEncyclopedia<2, Real> >   & a_brit,
        shared_ptr<GeometryService<2> >        & a_geoserv,
        const DisjointBoxLayout                & a_grids,
        const Box                              & a_domain,
        const Real                             & a_dx,
        const Real                             & a_viscosity,
        const Real                             & a_permeability,
        const Real                             & a_diffusivity,
        const IntVect                          & a_nghost,
        ParabolicSolverType                      a_solver,
        EBIBC                                    a_ibc)
{
  CH_TIME("EBDarcy::define");
  m_ibc           = a_ibc;
  m_brit          = a_brit;
  m_geoserv       = a_geoserv;
  m_grids         = a_grids;
  m_domain        = a_domain;
  m_dx            = a_dx;
  m_nghost        = a_nghost;
  m_viscosity     = a_viscosity   ;
  m_permeability  = a_permeability;
  m_diffusivity   = a_diffusivity ;
  m_graphs = m_geoserv->getGraphs(m_domain);

  
  m_velo     = shared_ptr<EBLevelBoxData<CELL, DIM> >
    (new EBLevelBoxData<CELL, DIM>(m_grids, m_nghost, m_graphs));
  m_gphi     = shared_ptr<EBLevelBoxData<CELL, DIM> >
    (new EBLevelBoxData<CELL, DIM>(m_grids, m_nghost, m_graphs));
  m_sour     = shared_ptr<EBLevelBoxData<CELL, 1> >
    (new EBLevelBoxData<CELL,    1>(m_grids, m_nghost, m_graphs));
  m_rhs      = shared_ptr<EBLevelBoxData<CELL, 1> >
    (new EBLevelBoxData<CELL,    1>(m_grids, m_nghost, m_graphs));
  m_scal     = shared_ptr<EBLevelBoxData<CELL, 1  > >
    (new EBLevelBoxData<CELL,    1>(m_grids, m_nghost, m_graphs));

  string stenname = StencilNames::Poisson2;
  string bcname = StencilNames::Neumann;

  auto cell_dict = m_brit->m_cellToCell;
  Real alpha = 1; Real beta = 1; //these get reset before solve
  string helmnames[2*DIM];
  a_ibc.scalarDiffusionStencilStrings(helmnames);
  m_helmholtz = shared_ptr<EBMultigrid> 
    (new EBMultigrid(cell_dict, m_geoserv, alpha, beta, m_dx, m_grids,  
                     stenname, helmnames, bcname, m_domain, m_nghost));

  if(a_solver == BackwardEuler)
  {
    m_heatSolver = shared_ptr<BaseEBParabolic>
      (new EBBackwardEuler(m_helmholtz, m_geoserv, m_grids, m_domain, m_nghost));
  }
  else if (a_solver == CrankNicolson)
  {
    m_heatSolver = shared_ptr<BaseEBParabolic>
      (new EBCrankNicolson(m_helmholtz, m_geoserv, m_grids, m_domain, m_nghost));
  }
  else if (a_solver == TGA)
  {
    m_heatSolver = shared_ptr<BaseEBParabolic>
      (new EBTGA(m_helmholtz, m_geoserv, m_grids, m_domain, m_nghost));
  }
  else
  {
    MayDay::Error("unaccounted-for solver type");
  }
  m_advectOp = shared_ptr<EBAdvection>
    (new EBAdvection(m_brit, m_geoserv, m_velo, m_grids, m_domain, m_dx, a_ibc, m_nghost));
  
  m_ccProj  = shared_ptr<EBCCProjector>
    (new EBCCProjector(m_brit, m_geoserv, m_grids, m_domain, m_dx, m_nghost, a_ibc));
}
/*******/ 
void 
EBDarcy::
run(unsigned int a_max_step,
    Real         a_max_time,
    Real         a_cfl,
    Real         a_fixedDt,
    Real         a_tol,
    unsigned int a_maxIter,
    int          a_outputInterval,
    Real         a_coveredVal)
{
  CH_TIME("EBDarcy::run");
  bool usingFixedDt = (a_fixedDt > 0);
  bool doFileOutput = (a_outputInterval > 0);
  Real dt;

  pout() << "initializing Darcy velocity field = potential flow solution"  << endl;
  pout() << "with inflow-outflow boundary conditions (xlow to xhigh)"  << endl;
  initializeVelocity(a_tol, a_maxIter);

  //get the time step
  if(usingFixedDt)
  {
    dt = a_fixedDt;
  }
  else
  {
    dt = computeDt(a_cfl);
  }


  if(doFileOutput)
  {
    outputToFile(0, a_coveredVal, dt, 0.0);
  }

  m_time = 0;
  m_step = 0;
  while((m_step < a_max_step) && (m_time < a_max_time))
  {
    pout() << "step = " << m_step << ", time = " << m_time << " dt = " << dt << endl;

    pout() << "advecting and difussing and reacting and all that " << endl;
    advanceScalar(dt, a_tol, a_maxIter);
    m_step++;
    m_time += dt;
    if((doFileOutput) && (m_step % a_outputInterval == 0))
    {
      outputToFile(m_step, a_coveredVal, dt, m_time);
    }
  }
}
/*******/ 
void
EBDarcy::
initializeVelocity(Real         a_tol,
                   unsigned int a_maxIter)
{
  m_ccProj->project(*m_velo, *m_gphi, a_tol, a_maxIter);
}
    
/*******/ 
Real
EBDarcy::
computeDt(Real a_cfl) const
{
  CH_TIME("EBDarcy::computeDt");
  Real dtval;
  Real dtCFL = 999999999.;
  Real maxvel = 0;
  for(int idir = 0; idir < DIM; idir++)
  {
    maxvel = std::max(maxvel, m_velo->maxNorm(idir));
  }
  if(maxvel > 1.0e-16)
  {
    dtCFL = a_cfl*m_dx/maxvel;
    dtval = dtCFL;
    pout() << "maxvel = " << maxvel << ", dx = " << m_dx << ", dt = " << dtval << endl;
  }    
  else
  {
    pout() << "velocity seems to be zero--setting dt to dx" << endl;
    dtval = m_dx;
  }

  return dtval;
}
/*******/ 
PROTO_KERNEL_START 
void ParabolicRHSF(Var<Real, 1>    a_rhs,
                   Var<Real, 1>    a_divuphi,
                   Var<Real, 1>    a_source)
{
  a_rhs(0) = -a_divuphi(0) + a_source(0);
}
PROTO_KERNEL_END(ParabolicRHSF, ParabolicRHS)
/*******/ 
void
EBDarcy::
getReactionSourceTerm()
{
  m_sour->setVal(0.);
}
/*******/ 
void
EBDarcy::
advanceScalar(Real a_dt,
              Real         a_tol,    
              unsigned int a_maxIter)
{
  auto& scal = *m_scal;
  CH_TIME("EBDarcy::advanceScalar");
  scal.exchange();
  //A strange user interface, I grant you.  I did not think I would need this back door
  pout() << "get advection term divuphi" << endl;
  m_advectOp->hybridDivergence(scal, a_dt);
  //sets m_sour
  pout() << "getting reaction term R" << endl;
  getReactionSourceTerm();
  EBLevelBoxData<CELL, 1>& divuphi = m_advectOp->m_hybridDiv;

  pout() << "assembling  rhs = -divuphi + R" << endl;
  DataIterator dit = m_grids.dataIterator();
  for(unsigned int ibox = 0; ibox < dit.size(); ibox++)
  {
    unsigned long long int numflopspt = 2;
    Bx grbx = ProtoCh::getProtoBox(m_grids[dit[ibox]]);
    auto & divuphifab =     divuphi[dit[ibox]];
    auto & reactsour  =   (*m_sour)[dit[ibox]];
    auto & rhs        =   (* m_rhs)[dit[ibox]];

    ebforallInPlace(numflopspt, "ParabolicRHS", ParabolicRHS, grbx, rhs, divuphifab, reactsour);
  }
  pout() << "calling heat solver for variable "  << endl;
  //advance the parabolic equation
  m_heatSolver->advanceOneStep((*m_scal), (*m_rhs), m_diffusivity, a_dt, a_tol, a_maxIter);
}
/*******/ 
void
EBDarcy::
outputToFile(unsigned int a_step, Real a_coveredval, Real a_dt, Real a_time) const
{
  CH_TIME("EBDarcy::outputToFile");
  const EBLevelBoxData<CELL, 1> & kappa = m_advectOp->m_kappa;
  if(a_step == 0)
  {
    string filevelo = string("velo.")  + string(".hdf5");
    string filegphi = string("gphi.")  + string(".hdf5");
    writeEBLevelHDF5<DIM>(filevelo,  *m_velo, kappa, m_domain, m_graphs, a_coveredval, m_dx, a_dt, a_time);
    writeEBLevelHDF5<DIM>(filegphi,  *m_gphi, kappa, m_domain, m_graphs, a_coveredval, m_dx, a_dt, a_time);
  }
  string filescal = string("scal.")   + std::to_string(a_step) + string(".hdf5");
  writeEBLevelHDF5<1>(  filescal,  *m_scal, kappa, m_domain, m_graphs, a_coveredval, m_dx, a_dt, a_time);
}
/*******/ 
#include "Chombo_NamespaceFooter.H"

