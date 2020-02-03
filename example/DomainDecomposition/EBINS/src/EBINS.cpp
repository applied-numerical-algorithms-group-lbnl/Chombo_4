
#include "EBINS.H"
#include "Chombo_NamespaceHeader.H"

/*******/
EBINS::
EBINS(shared_ptr<EBEncyclopedia<2, Real> >   & a_brit,
      shared_ptr<GeometryService<2> >        & a_geoserv,
      const DisjointBoxLayout                & a_grids,
      const Box                              & a_domain,
      const Real                             & a_dx,
      const Real                             & a_solverTolerance,
      const unsigned int                     & a_maxSolverIterations,
      const Real                             & a_coveredval,
      const IntVect                          & a_nghostdst)
{
}
/*******/ 
void 
EBINS::
run(unsigned int a_max_step,
    Real         a_max_time,
    Real         a_cfl,
    unsigned int a_outputInterval)
{
}

/*******/ 
void
EBINS::
computeDt(Real a_cfl) const
{
  Real dtval;
  Real maxvel = 0;
  for(int idir = 0; idir < DIM; idir++)
  {
    maxvel = std::max(maxvel, m_velo.maxNorm(idir));
  }
  if(maxvel > 1.0e-16)
  {
    dtval = a_cfl*m_dx/maxvel;
    pout() << "maxvel = " << maxvel << ", dx = " << a_dx << ", dt = " << a_dt << endl;
  }    
  else
  {
    pout() << "velocity seems to be zero--setting dt to dx" << endl;
    dtval = m_dx;
  }
  m_dt =  dtval;
}
/*******/ 

void
EBINS::
advanceVelocity()
{
}
/*******/ 
void
EBINS::
advanceScalar()
{
}
/*******/ 
void
EBINS::
computeVorticity()
{
}
/*******/ 
void
EBINS::
outputToFile(unsigned int a_step) const
{
  const EBLevelBoxData<CELL, 1> & kappa = m_advectOp->getKappa();
  string filescal = string("scal.") + std::to_string(a_step) + string(".hdf5");
  string filevelo = string("velo.") + std::to_string(a_step) + string(".hdf5");
  string filegphi = string("gphi.") + std::to_string(a_step) + string(".hdf5");
  string filevort = string("vort.") + std::to_string(a_step) + string(".hdf5");
  writeEBLevelHDF5<1>(  filescal,  *m_scal, kappa, m_domain, m_graphs, m_coveredval, m_dx, m_dt, m_time);
  writeEBLevelHDF5<DIM>(filevelo,  *m_velo, kappa, m_domain, m_graphs, m_coveredval, m_dx, m_dt, m_time);
  writeEBLevelHDF5<DIM>(filegphi,  *m_gphi, kappa, m_domain, m_graphs, m_coveredval, m_dx, m_dt, m_time);
#if DIM==2
  writeEBLevelHDF5<1>(  filevort,  *m_vort, kappa, m_domain, m_graphs, m_coveredval, m_dx, m_dt, m_time);
#else  
  writeEBLevelHDF5<DIM>(filevort,  *m_vort, kappa, m_domain, m_graphs, m_coveredval, m_dx, m_dt, m_time);
#endif  
}
/*******/ 
#include "Chombo_NamespaceFooter.H"

