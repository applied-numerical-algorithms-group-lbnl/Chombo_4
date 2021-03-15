
#include "EBINS.H"
#include "EBParabolicIntegrators.H"
#include "Chombo_ParmParse.H"
#include "Chombo_CH_HDF5.H"
#include "Chombo_NamespaceHeader.H"
#include "DebugFunctions.H"
using Proto::Var;
/*******/
EBINS::
EBINS(shared_ptr<EBEncyclopedia<2, Real> >   & a_brit,
      shared_ptr<GeometryService<2> >        & a_geoserv,
      const DisjointBoxLayout                & a_grids,
      const Box                              & a_domain,
      const Real                             & a_dx,
      const Real                             & a_viscosity,
      const IntVect                          & a_nghost,
      ParabolicSolverType                      a_solver,
      EBIBC                                    a_ibc,
      unsigned int                             a_num_species,  
      vector<Real> a_diffusionCoeffs)
{
  CH_TIME("EBINS::define");
  m_diffusionCoefs = a_diffusionCoeffs;
  PR_assert(m_diffusionCoefs.size() >= a_num_species);
  m_species.resize(a_num_species);
  m_brit                = a_brit;
  m_geoserv             = a_geoserv;
  m_grids               = a_grids;
  m_domain              = a_domain;
  m_dx                  = a_dx;
  m_nghost              = a_nghost;
  m_viscosity           = a_viscosity;
  defineInternals(a_ibc, a_num_species, a_viscosity, a_solver);
}  
/*******/
void
EBINS::
defineInternals(EBIBC                a_ibc,
                unsigned int         a_num_species,
                Real                 a_viscosity,
                ParabolicSolverType  a_solver)
{
  m_graphs = m_geoserv->getGraphs(m_domain);
  m_exchangeCopier.exchangeDefine(m_grids, m_nghost*IntVect::Unit);
  m_copyCopier.define(m_grids, m_grids, m_nghost*IntVect::Unit);
  
  m_velo     = shared_ptr<EBLevelBoxData<CELL, DIM> >
    (new EBLevelBoxData<CELL, DIM>(m_grids, m_nghost, m_graphs));
  m_sour     = shared_ptr<EBLevelBoxData<CELL, DIM> >
    (new EBLevelBoxData<CELL, DIM>(m_grids, m_nghost, m_graphs));
  m_divuu    = shared_ptr<EBLevelBoxData<CELL, DIM> >
    (new EBLevelBoxData<CELL, DIM>(m_grids, m_nghost, m_graphs));
  m_gphi     = shared_ptr<EBLevelBoxData<CELL, DIM> >
    (new EBLevelBoxData<CELL, DIM>(m_grids, m_nghost, m_graphs));
  m_scal     = shared_ptr<EBLevelBoxData<CELL, 1  > >
    (new EBLevelBoxData<CELL, 1  >(m_grids, m_nghost, m_graphs));
  m_sourDiff  = shared_ptr<EBLevelBoxData<CELL, 1  > >
    (new EBLevelBoxData<CELL, 1  >(m_grids, m_nghost, m_graphs));
  m_rhsDiff  = shared_ptr<EBLevelBoxData<CELL, 1  > >
    (new EBLevelBoxData<CELL, 1  >(m_grids, m_nghost, m_graphs));
  for(int ispec = 0; ispec < a_num_species; ispec++)
  {
    m_species[ispec]     = shared_ptr<EBLevelBoxData<CELL, 1  > >
      (new EBLevelBoxData<CELL, 1  >(m_grids, m_nghost, m_graphs));
  }
  
  string stenname = StencilNames::Poisson2;
  string bcname;
  if(a_viscosity == 0)
  {
    bcname = StencilNames::Neumann;
    m_eulerCalc = true;
  }
  else
  {
    bcname = StencilNames::Dirichlet;
    m_eulerCalc = false;
  }

  auto cell_dict = m_brit->m_cellToCell;
  Real alpha = 1; Real beta = 1; //these get reset before solve
  string helmnamesVelo[2*DIM];
  string helmnamesSpec[2*DIM];
  a_ibc.helmholtzStencilStrings(      helmnamesVelo);
  a_ibc.scalarDiffusionStencilStrings(helmnamesSpec);
  if(!m_eulerCalc)
  {
    m_helmholtzVelo = shared_ptr<EBMultigrid> 
      (new EBMultigrid(cell_dict, m_geoserv, alpha, beta, m_dx, m_grids,  
                       stenname, helmnamesVelo, bcname, m_domain, m_nghost));
  }
  m_helmholtzSpec = shared_ptr<EBMultigrid> 
    (new EBMultigrid(cell_dict, m_geoserv, alpha, beta, m_dx, m_grids,  
                     stenname, helmnamesSpec, StencilNames::Neumann, m_domain, m_nghost));

  if(a_solver == BackwardEuler)
  {
    if(!m_eulerCalc)
    {
      m_heatSolverVelo = shared_ptr<BaseEBParabolic>
        (new EBBackwardEuler(m_helmholtzVelo, m_geoserv, m_grids, m_domain, m_nghost));
    }
    m_heatSolverSpec = shared_ptr<BaseEBParabolic>
      (new EBBackwardEuler(m_helmholtzSpec, m_geoserv, m_grids, m_domain, m_nghost));
  }
  else if (a_solver == CrankNicolson)
  {
    if(!m_eulerCalc)
    {
      m_heatSolverVelo = shared_ptr<BaseEBParabolic>
        (new EBCrankNicolson(m_helmholtzVelo, m_geoserv, m_grids, m_domain, m_nghost));
    }
    m_heatSolverSpec = shared_ptr<BaseEBParabolic>
      (new EBCrankNicolson(m_helmholtzSpec, m_geoserv, m_grids, m_domain, m_nghost));
  }
  else if (a_solver == TGA)
  {
    if(!m_eulerCalc)
    {
      m_heatSolverVelo = shared_ptr<BaseEBParabolic>
        (new EBTGA(m_helmholtzVelo, m_geoserv, m_grids, m_domain, m_nghost));
    }
    m_heatSolverSpec = shared_ptr<BaseEBParabolic>
      (new EBTGA(m_helmholtzSpec, m_geoserv, m_grids, m_domain, m_nghost));
  }
  else
  {
    MayDay::Error("unaccounted-for solver type");
  }

  m_advectOp = shared_ptr<EBAdvection>
    (new EBAdvection(m_brit, m_geoserv, m_velo, m_grids, m_domain, m_dx, a_ibc, m_nghost));
  
  m_ccProj  = shared_ptr<EBCCProjector>
    (new EBCCProjector(m_brit, m_geoserv, m_grids, m_domain, m_dx, m_nghost, a_ibc));
  m_macProj = m_ccProj->m_macprojector;

  m_bcgAdvect = shared_ptr<BCGVelAdvect>
    (new BCGVelAdvect(m_macProj, m_helmholtzVelo, m_brit, m_geoserv, m_velo,
                      m_grids, m_domain, m_dx, m_viscosity, a_ibc, m_nghost, m_eulerCalc));


}
/*******/ 
void 
EBINS::
run(unsigned int a_max_step,
    Real         a_max_time,
    unsigned int a_startingStep,
    Real         a_startingTime,
    Real         a_cfl,
    Real         a_fixedDt,
    Real         a_tol,
    unsigned int a_numIterPres,
    unsigned int a_maxIter,
    int          a_plotfileInterval,
    int          a_checkpointInterval,
    Real         a_coveredVal)
{
  CH_TIME("EBINS::run");
  m_step = a_startingStep;
  m_time = a_startingTime;
  //Welcome to our standard interface that is at least not over specified.
  //Negative intervals turn intervalled stuff off.
  //Negative fixed time step means variable time step.
  //It is a quixotic interface at best.
  bool usingFixedDt     = (a_fixedDt            > 0);
  bool writePlotfiles   = (a_plotfileInterval   > 0);
  bool writeCheckpoints = (a_checkpointInterval > 0);

  //project the resulting field
  pout() << "projecting initial velocity"<< endl;
  auto & velo =  (*m_velo );
  auto & gphi =  (*m_gphi );

  if(m_step == 0)
  {
    pout() << "We are starting from scratch so we are projecting the initial velocity."<< endl;
    m_ccProj->project(velo, gphi, a_tol, a_maxIter);
  }
  else
  {
    pout() << "We are starting from a checkpoint so we just use the starting velocity. " << endl;
  }
    
  velo.exchange(m_exchangeCopier);

  //get the time step
  if(usingFixedDt)
  {
    m_dt = a_fixedDt;
  }
  else
  {
    m_dt = computeDt(a_cfl);
  }
  
  if(m_step == 0)
  {
    pout() << "We are starting from scratch so we are " << endl;
    pout() << "getting initial pressure using fixed point iteration " << endl;
    initializePressure(m_dt, a_tol, a_numIterPres, a_maxIter);
  }
  else
  {
    pout() << "We are starting from a checkpoint so we just use the starting pressure." << endl;
  }

  if(writePlotfiles && (m_step % a_plotfileInterval == 0))
  {
    writePlotFile(a_coveredVal);
  }
  if(writeCheckpoints && (m_step % a_checkpointInterval == 0))
  {
    writeCheckpointFile();
  }

  while((m_step < a_max_step) && (m_time < a_max_time))
  {
    pout() << "step = " << m_step << ", time = " << m_time << " dt = " << m_dt << endl;

    pout() << "advancing velocity and pressure fields " << endl;
    advanceVelocityAndPressure(m_dt, a_tol, a_maxIter);

    pout() << "advancing passive scalar" << endl;
    advanceScalar(m_dt);

    pout() << "advancing species advection/diffusion" << endl;
    advanceSpecies(m_dt, a_tol, a_maxIter);
    
    m_step++;
    m_time += m_dt;
    if(!usingFixedDt)
    {
      m_dt = computeDt(a_cfl);
    }
    if(writePlotfiles && (m_step % a_plotfileInterval == 0))
    {
      writePlotFile(a_coveredVal);
    }
    if(writeCheckpoints && (m_step % a_checkpointInterval == 0))
    {
      writeCheckpointFile();
    }
  }
}
/*******/ 
void
EBINS::
initializePressure(Real         a_dt,
                   Real         a_tol,
                   unsigned int a_numIterPres, 
                   unsigned int a_maxIter)
{
  CH_TIME("EBINS::initializePressure");
  auto & velo =  (*m_velo );
  auto & gphi =  (*m_gphi );

  EBLevelBoxData<CELL, DIM> velosave(m_grids, m_nghost, m_graphs);
  Interval interv(0, DIM-1);
  velo.copyTo(interv, velosave, interv, m_copyCopier);
  //gphi.setVal(0.);
  for(int iter = 0; iter < a_numIterPres; iter++)
  {
    advanceVelocityAndPressure(a_dt, a_tol, a_maxIter);
    velosave.copyTo(interv, velo, interv, m_copyCopier);
  }
  gphi.exchange(m_exchangeCopier);
}
/*******/ 
Real
EBINS::
computeDt(Real a_cfl) const
{
  CH_TIME("EBINS::computeDt");
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

  ParmParse pp;
  Real dtStokes = m_dx*m_dx/m_viscosity;
  bool use_stokes_dt = false;
  pp.query("use_stokes_dt", use_stokes_dt);
  if(use_stokes_dt && dtCFL > dtStokes)
  {
    dtval = dtStokes;
    pout() << "Using Stokes dt = " << dtval << endl;
  }
  return dtval;
}
/*******/ 
void
EBINS::
getAdvectiveDerivative(Real a_dt, Real a_tol, unsigned int a_maxIter)    
{
  auto& divuu = *m_divuu;
  auto& velo  = *m_velo;
  m_bcgAdvect->hybridVecDivergence(divuu, velo, a_dt, a_tol, a_maxIter);
  divuu.exchange(m_exchangeCopier);
}
/*******/ 
PROTO_KERNEL_START 
void EulerAdvanceF(Var<Real, DIM>    a_velo,
                   Var<Real, DIM>    a_divuu,
                   Var<Real, DIM>    a_gradp,
                   Real              a_dt)
{
  for(unsigned int idir = 0; idir < DIM; idir++)
  {
    a_velo(idir) -= a_dt*(a_divuu(idir) + a_gradp(idir));
  }
}
PROTO_KERNEL_END(EulerAdvanceF, EulerAdvance)

/*******/ 
void
EBINS::
advanceVelocityEuler(Real a_dt)
{
  CH_TIME("EBINS::advanceVelocityEuler");
  auto & velo =  (*m_velo );
  auto & udel =  (*m_divuu);
  auto & gphi =  (*m_gphi );

  DataIterator dit = m_grids.dataIterator();
  for(unsigned int ibox = 0; ibox < dit.size(); ibox++)
  {
    unsigned long long int numflopspt = 3*DIM;
    Bx grbx = ProtoCh::getProtoBox(m_grids[dit[ibox]]);
    auto & velofab =  velo[dit[ibox]];
    auto & udelfab =  udel[dit[ibox]];
    auto & gphifab =  gphi[dit[ibox]];  

    ebforallInPlace(numflopspt, "EulerAdvance", EulerAdvance, grbx,  
                    velofab, udelfab, gphifab, a_dt);
  }
}
/*******/ 
PROTO_KERNEL_START 
void DiffusionRHSF(Var<Real, 1>    a_rhs,
                   Var<Real, 1>    a_divuphi,
                   Var<Real, 1>    a_source)
{
  a_rhs(0) = -a_divuphi(0) + a_source(0);
}
PROTO_KERNEL_END(DiffusionRHSF, DiffusionRHS)
/*******/ 
PROTO_KERNEL_START 
void ParabolicRHSF(Var<Real, 1>    a_rhs,
                   Var<Real, 1>    a_divuu,
                   Var<Real, 1>    a_gradp)
{
  a_rhs(0) = -a_divuu(0) - a_gradp(0);
}
PROTO_KERNEL_END(ParabolicRHSF, ParabolicRHS)

/*******/ 
void
EBINS::
advanceVelocityNavierStokes(Real a_dt,                          
                            Real a_tol,                         
                            unsigned int a_maxIter)
{
  CH_TIME("EBINS::advanceVelocityNavierStokes");
  auto & sour  = *m_sour ;
  auto & velo  = *m_velo ;
  auto & divuu = *m_divuu;
  auto & gphi  = *m_gphi ;
  int ideb = 0;
  DataIterator dit = m_grids.dataIterator();
  for(unsigned int idir = 0; idir < DIM; idir++)
  {
    EBLevelBoxData<CELL, 1> scalRHS, scalVelo, scalUdelu, scalGphi;
    scalRHS.define<DIM>(  sour ,idir, m_graphs);
    scalVelo.define<DIM>( velo ,idir, m_graphs);
    scalUdelu.define<DIM>(divuu,idir, m_graphs);
    scalGphi.define<DIM>( gphi ,idir, m_graphs);

    //set source of parabolic advance to -(gphi + udelu)
    for(unsigned int ibox = 0; ibox < dit.size(); ibox++)
    {
      unsigned long long int numflopspt = 2;
      Bx grbx = ProtoCh::getProtoBox(m_grids[dit[ibox]]);
      auto & sour =    scalRHS[dit[ibox]];
      auto & udel =  scalUdelu[dit[ibox]];
      auto & gphi =  scalGphi[ dit[ibox]];  

      ebforallInPlace(numflopspt, "ParabolicRHS", ParabolicRHS, grbx, sour, udel, gphi);
    }
    pout() << "calling heat solver for variable " << idir << endl;
    //advance the parabolic equation
    m_heatSolverVelo->advanceOneStep(scalVelo, scalRHS, m_viscosity, a_dt, a_tol, a_maxIter);
    
    ideb++;
  }
 ideb++;
}
/*******/ 
PROTO_KERNEL_START 
void AddDtGphiToVeloF(Var<Real, DIM>    a_velo,
                      Var<Real, DIM>    a_gradp,
                      Real              a_dt)
{
  for(int idir = 0; idir < DIM; idir++)
  {
    a_velo(idir) += a_dt*(a_gradp(idir));
  }
}
PROTO_KERNEL_END(AddDtGphiToVeloF, AddDtGphiToVelo)
/*******/ 
PROTO_KERNEL_START 
void DivideOutDtF(Var<Real, DIM>    a_gradp,
                  Real              a_dt)
{
  for(int idir = 0; idir < DIM; idir++)
  {
    a_gradp(idir) /= a_dt;
  }
}
PROTO_KERNEL_END(DivideOutDtF, DivideOutDt)
/*******/ 
void
EBINS::
projectVelocityAndCorrectPressure(Real a_dt,
                                  Real a_tol,                         
                                  unsigned int a_maxIter)
{
  CH_TIME("EBINS::projectVelocityAndCorrectPressure");
  auto & velo =  (*m_velo);
  auto & gphi =  (*m_gphi);
  //make w = v + dt*gphi
  DataIterator dit = m_grids.dataIterator();
  for(unsigned int ibox = 0; ibox < dit.size(); ibox++)
  {
    unsigned long long int numflopspt = 2*DIM;
    Bx grbx = ProtoCh::getProtoBox(m_grids[dit[ibox]]);
    auto & velofab =  velo[dit[ibox]];
    auto & gphifab =  gphi[dit[ibox]];  

    ebforallInPlace(numflopspt, "AddDtGphiToVelo", AddDtGphiToVelo, grbx, velofab, gphifab, a_dt);
  }

  //project the resulting field
  pout() << "cc projecting vel + gphi*dt" << endl;
  m_ccProj->project(velo, gphi, a_tol, a_maxIter);

  //the resulting pressure  is = dt * gphi so we need to divide out the dt
  for(unsigned int ibox = 0; ibox < dit.size(); ibox++)
  {
    unsigned long long int numflopspt = 2*DIM;
    Bx grbx = ProtoCh::getProtoBox(m_grids[dit[ibox]]);
    auto & gphifab =  gphi[dit[ibox]];  

    ebforallInPlace(numflopspt, "DivideOutDt", DivideOutDt, grbx, gphifab, a_dt);
  }
}
/*******/ 
void
EBINS::
advanceVelocityAndPressure(Real a_dt,                          
                           Real a_tol,                         
                           unsigned int a_maxIter)
{
  CH_TIME("EBINS::advanceVelocityAndPressure");
  //get udelu
  getAdvectiveDerivative(a_dt, a_tol, a_maxIter);

  if(m_eulerCalc)
  {
    advanceVelocityEuler(a_dt);
  }
  else
  {
    advanceVelocityNavierStokes(a_dt, a_tol, a_maxIter);
  }
  projectVelocityAndCorrectPressure(a_dt, a_tol, a_maxIter);
}
/*******/ 
void
EBINS::
advanceScalar(Real a_dt)
              
{
  CH_TIME("EBINS::advanceScalar");
  
  auto& scal = *m_scal;
  scal.exchange(m_exchangeCopier);
  Real fluxval;
  ParmParse pp;
  pp.get("scalar_inflow_value",   fluxval);
  m_advectOp->advance(scal, a_dt, fluxval);
  scal.exchange(m_exchangeCopier);
}
/*******/ 
void
EBINS::
advanceSpecies(Real a_dt,
               Real         a_tol,    
               unsigned int a_maxIter)
{
  CH_TIME("EBINS::advanceSpecies");
  for(unsigned int ispec = 0; ispec < m_species.size(); ispec++)
  {
    pout() << "advancing species number "<< ispec << endl;
    auto& spec = *m_species[ispec];

    spec.exchange(m_exchangeCopier);
    //A strange user interface, I grant you.  I did not think I would need this back door
    pout() << "get advection term divuphi" << endl;
    Real fluxval;
    ParmParse pp;
    string scalname = string("species_inflow_value_") + to_string(ispec);
    pp.get(scalname.c_str(),   fluxval);
    
    m_advectOp->hybridDivergence(spec, a_dt, fluxval);
    //sets m_sour
    pout() << "getting reaction term R" << endl;
    getReactionSourceTerm(ispec);
    EBLevelBoxData<CELL, 1>& divuphi = m_advectOp->m_hybridDiv;

    pout() << "assembling  rhs = -divuphi + R" << endl;
    DataIterator dit = m_grids.dataIterator();
    for(unsigned int ibox = 0; ibox < dit.size(); ibox++)
    {
      unsigned long long int numflopspt = 2;
      Bx grbx = ProtoCh::getProtoBox(m_grids[dit[ibox]]);
      auto & divuphifab =         divuphi[dit[ibox]];
      auto & reactsour  =   (*m_sourDiff)[dit[ibox]];
      auto & rhs        =   (* m_rhsDiff)[dit[ibox]];

      ebforallInPlace(numflopspt, "DiffusionRHS", DiffusionRHS, grbx, rhs, divuphifab, reactsour);
    }
    pout() << "calling heat solver for variable "  << endl;
    //advance the parabolic equation
    Real thiscoef = m_diffusionCoefs[ispec];
    m_heatSolverSpec->advanceOneStep(spec, (*m_rhsDiff),
                                     thiscoef, a_dt, a_tol, a_maxIter);
  }
}
void
EBINS::
getReactionSourceTerm(unsigned int a_ispec)
{
  CH_TIME("EBINS::getReactionSourceTerm");
  m_sour->setVal(0.);
}
/*******/ 
void
EBINS::
writePlotFile(Real a_coveredval) const
{
  CH_TIME("EBINS::writePlotFile");
  const EBLevelBoxData<CELL, 1> & kappa = m_advectOp->m_kappa;
  string filescal = string("scal.step_")
    + std::to_string(m_step) + string(".hdf5");
  string filevelo = string("velo.step_")
    + std::to_string(m_step) + string(".hdf5");
  string filegphi = string("gphi.step_")
    + std::to_string(m_step) + string(".hdf5");
  writeEBLevelHDF5<1>(  filescal,  *m_scal, kappa, m_domain,
                        m_graphs, a_coveredval, m_dx, m_dt, m_time);
  writeEBLevelHDF5<DIM>(filevelo,  *m_velo, kappa, m_domain,
                        m_graphs, a_coveredval, m_dx, m_dt, m_time);
  writeEBLevelHDF5<DIM>(filegphi,  *m_gphi, kappa, m_domain,
                        m_graphs, a_coveredval, m_dx, m_dt, m_time);
  for(unsigned int ispec = 0; ispec < m_species.size(); ispec++)
  {
    string filespec = string("spec_") + std::to_string(ispec)
      + string(".step_") + std::to_string(m_step) + string(".hdf5");
    const auto& species = *(m_species[ispec]);
    writeEBLevelHDF5<1>(filespec,  species, kappa,
                        m_domain, m_graphs, a_coveredval,
                        m_dx, m_dt, m_time);
  }
}
/*******/
void
EBINS::
readDataFromCheckpoint(Real         & a_curTime,
                       unsigned int & a_curStep,
                       const string & a_checkpointName)
{
  CH_TIME("EBINS::readDataFromCheckpointFile");

  HDF5HeaderData header;
  HDF5Handle handle(a_checkpointName.c_str(), HDF5Handle::OPEN_RDONLY);

  header.readFromFile(handle);
  a_curTime = header.m_real["cur_time"];
  a_curStep = header.m_int ["cur_step"];

  m_velo->readFromCheckPoint(    handle, string("velo"));
  m_gphi->readFromCheckPoint(    handle, string("gphi"));
  m_scal->readFromCheckPoint(    handle, string("scal"));
  m_sour->readFromCheckPoint(    handle, string("sour"));
  m_rhsDiff->readFromCheckPoint( handle, string("rhsDiff"));
  m_sourDiff->readFromCheckPoint(handle, string("sourDiff"));
  m_divuu->readFromCheckPoint(   handle, string("divuu"));
  for(unsigned int ispec = 0; ispec < m_species.size(); ispec++)
  {
    string specname = string("spec_") + to_string(ispec);
    m_species[ispec]->readFromCheckPoint(handle, specname);
  }
  handle.close();
}
/*******/
void
EBINS::
writeCheckpointFile()
{
  CH_TIME("EBINS::writeCheckpointFile");
  string checkpointName = string("check.step_") + to_string(m_step) + string(".hdf5");
  HDF5HeaderData header;
  

  HDF5Handle handle(checkpointName.c_str(), HDF5Handle::CREATE);
  header.m_real["cur_time"]               = m_time;
  header.m_int ["cur_step"]               = m_step;  
  header.writeToFile(handle);
  
  m_velo->writeToCheckPoint(    handle, string("velo"));
  m_gphi->writeToCheckPoint(    handle, string("gphi"));
  m_scal->writeToCheckPoint(    handle, string("scal"));
  m_sour->writeToCheckPoint(    handle, string("sour"));
  m_rhsDiff->writeToCheckPoint( handle, string("rhsDiff"));
  m_sourDiff->writeToCheckPoint(handle, string("sourDiff"));
  m_divuu->writeToCheckPoint(   handle, string("divuu"));
  for(unsigned int ispec = 0; ispec < m_species.size(); ispec++)
  {
    string specname = string("spec_") + to_string(ispec);
    m_species[ispec]->writeToCheckPoint(handle, specname);
  }
  handle.close();
}
/*******/ 
#include "Chombo_NamespaceFooter.H"
