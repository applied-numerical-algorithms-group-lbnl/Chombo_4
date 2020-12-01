
#include <cmath>
#include <memory>
#include "EBMACProjector.H"
#include "EBAdvectionFunctions.H"
#include "Chombo_ParmParse.H"
#include "Chombo_NamespaceHeader.H"
///
void
EBMACProjector::
define(shared_ptr<EBEncyclopedia<2, Real> >   & a_brit,
       shared_ptr<GeometryService<2> >        & a_geoserv,
       const DisjointBoxLayout                & a_grids,
       const Box                              & a_domain,
       const Real                             & a_dx,
       const IntVect                          & a_nghost,
       const EBIBC                            & a_ebibc)
{
  CH_TIME("EBMACProjector::define");
  m_dx      = a_dx;
  m_ebibc   = a_ebibc;
  m_grids   = a_grids;
  m_domain  = a_domain;
  m_nghost  = a_nghost;
  m_grids   = a_grids;      
  m_domain  = a_domain;     
  m_brit    = a_brit;
  m_graphs  = a_geoserv->getGraphs(m_domain);

  defineInternals(a_geoserv);
}
///
void 
EBMACProjector::
defineInternals(shared_ptr<GeometryService<2> >        & a_geoserv)
{
  CH_TIME("EBMACProjector::defineInternals");
  m_exchangeCopier.exchangeDefine(m_grids, m_nghost);
  m_rhs.define(m_grids, m_nghost, m_graphs);
  m_phi.define(m_grids, m_nghost, m_graphs);
  auto ditch = m_brit->m_cellToCell;
  Real alpha = 0; Real beta = 1; //Poisson's eqn

  string bcnames[2*DIM];
  m_ebibc.projectionStencilStrings(bcnames);

  m_solver = shared_ptr<EBMultigrid>
    (new EBMultigrid(ditch, a_geoserv, alpha, beta, m_dx, m_grids, 
                     StencilNames::Poisson2, bcnames, StencilNames::Neumann,
                     m_domain, m_nghost));

  registerStencils();
}
////
void  
EBMACProjector::
registerStencils()
{

  CH_TIME("EBMACProjector::defineInternals");
  //false is because I do not need diagonal  weights for any of these stencils
  bool needDiag = false;
  
  //interpolates from face centers to centroids
  m_brit->registerFaceStencil(StencilNames::InterpToFaceCentroid, StencilNames::NoBC,    StencilNames::NoBC, m_domain, m_domain, needDiag);

  //increment divergence by face difference
  m_brit->registerFaceToCell( StencilNames::DivergeFtoC         , StencilNames::NoBC,    StencilNames::NoBC, m_domain, m_domain, needDiag);

  //face-centered gradient of cell-centered data
  m_brit->registerCellToFace( StencilNames::MACGradient         , StencilNames::NoBC,    StencilNames::NoBC, m_domain, m_domain, needDiag, Point::Ones(2));
}
/// 
void 
EBMACProjector::
project(EBLevelFluxData<1>   & a_velo,
        EBLevelFluxData<1>   & a_gphi,
        Real a_tol, unsigned int a_maxiter)
{
  CH_TIME("EBMACProjector::project");
  // set rhs = kappa*div (vel)
  kappaDivU(m_rhs, a_velo);

//  //begin debug
//  Real rhsmax = m_rhs.maxNorm(0);
//  pout() << "rhs of mac projection = " << rhsmax << endl;
//  exit(0);
  //end debug

  //solve kappa*lapl(phi) = kappa*divu
  m_solver->solve(m_phi, m_rhs, a_tol, a_maxiter);

  //gphi = grad(phi)
  //v := v - gphi
  DataIterator dit = m_grids.dataIterator();
  int ideb = 0;
  for(int ibox = 0; ibox < dit.size(); ibox++)
  {
    
//    Bx   grid   =  ProtoCh::getProtoBox(m_grids[dit[ibox]]);
    //get face fluxes and interpolate them to centroids
    for(unsigned int idir = 0; idir < DIM; idir++)
    {
      bool initToZero = true;
      m_brit->applyCellToFace(StencilNames::MACGradient, StencilNames::NoBC, m_domain,
                              a_gphi[dit[ibox]] ,m_phi[dit[ibox]], idir, ibox, initToZero, 1.0);
    }
    applyGradBoundaryConditions(a_gphi[dit[ibox]], dit[ibox]);
    a_velo[dit[ibox]] -= a_gphi[dit[ibox]];
    ideb++;
  }
}
///
void 
EBMACProjector::
setFaceStuff(int idir, Side::LoHiSide sit, EBFluxData<Real, 1>& a_flux, Bx valbx, Real fluxval)
{  
  Bx faceBx = valbx.faceBox(idir, sit);
  //unsigned long long int numflopspt = 0;
  if(idir == 0)
  {
    //ebforallInPlace(numflopspt, "setFluxVal", setFluxVal,  faceBx,  *a_flux.m_xflux, fluxval);
    //using non-eb forall because box restriction in eb land is broken right now.   This will
    //work if there nare no cut cells near the domain boundary
    auto& regdata = a_flux.m_xflux->getRegData();
    forallInPlaceBase(setFluxVal, faceBx, regdata, fluxval);
  }
  else if(idir == 1)
  {
    //ebforallInPlace(numflopspt, "setFluxVal", setFluxVal,  faceBx,  *a_flux.m_yflux, fluxval);
    //using non-eb forall because box restriction in eb land is broken right now.   This will
    //work if there nare no cut cells near the domain boundary
    auto& regdata = a_flux.m_yflux->getRegData();
    forallInPlaceBase(setFluxVal, faceBx, regdata, fluxval);
  }
#if DIM==3          
  else if(idir == 2)
  {
    //ebforallInPlace(numflopspt, "setFluxVal", setFluxVal,  faceBx,  *a_flux.m_zflux, fluxval);
    //using non-eb forall because box restriction in eb land is broken right now.   This will
    //work if there nare no cut cells near the domain boundary
    auto& regdata = a_flux.m_zflux->getRegData();
    forallInPlaceBase(setFluxVal, faceBx, regdata, fluxval);
  }
#endif
  else
  {
    MayDay::Error("bogus idir");
  }
}
///
void 
EBMACProjector::
applyVeloBoundaryConditions(EBFluxData<Real, 1> & a_flux,
                            const DataIndex     & a_dit)
{
  Box validBox = m_grids[a_dit];

  Bx dombx = ProtoCh::getProtoBox(m_domain);
  Bx valbx = ProtoCh::getProtoBox(validBox);
  for(SideIterator sit; sit.ok(); ++sit)
  {
    Point dombnd = dombx.boundary(sit());
    Point valbnd = valbx.boundary(sit());
    for(int idir = 0; idir < DIM; idir++)
    {
      if(dombnd[idir] == valbnd[idir])
      {
        int index = ebp_index(idir, sit());
        string bcstr = m_ebibc.m_domainBC[index];
        bool setstuff = true;
        Real fluxval = 0;
        if(bcstr == string("outflow"))
        {
          //do nothing, the upwind state should already be correct
          setstuff = false;
        }
        else if(bcstr == string("inflow"))
        {
          //velocities in this context are always normal velocity
          setstuff = true;
          ParmParse pp;
          pp.get("velocity_inflow_value", fluxval);
        }
        else if(bcstr == string("slip_wall"))
        {
          setstuff = true;
          fluxval = 0;  //velocities in this context are always normal velocity
        }
        else if(bcstr == string("no_slip_wall"))
        {
          setstuff = true;
          fluxval = 0; 
        }
        else
        {
          MayDay::Error("EBMACProjector: unrecognized bc");
        }
        if(setstuff)
        {
          setFaceStuff(idir, sit(),  a_flux, valbx, fluxval);
        }
      }
    }
  }
}
///
void 
EBMACProjector::
applyGradBoundaryConditions(EBFluxData<Real, 1> & a_flux,
                            const DataIndex     & a_dit)
{
  Box validBox = m_grids[a_dit];

  Bx dombx = ProtoCh::getProtoBox(m_domain);
  Bx valbx = ProtoCh::getProtoBox(validBox);
  for(SideIterator sit; sit.ok(); ++sit)
  {
    Point dombnd = dombx.boundary(sit());
    Point valbnd = valbx.boundary(sit());
    for(int idir = 0; idir < DIM; idir++)
    {
      if(dombnd[idir] == valbnd[idir])
      {
        int index = ebp_index(idir, sit());
        string bcstr = m_ebibc.m_domainBC[index];
        Real fluxval = 0;

        Bx faceBx = valbx.faceBox(idir, sit());
        //unsigned long long int numflopspt = 0;
        if(idir == 0)
        {
          //ebforallInPlace(numflopspt, "setFluxVal", setFluxVal,  faceBx,  *a_flux.m_xflux, fluxval);
          //using non-eb forall because box restriction in eb land is broken right now.   This will
          //work if there nare no cut cells near the domain boundary
          auto& regdata = a_flux.m_xflux->getRegData();
          forallInPlaceBase(setFluxVal, faceBx, regdata, fluxval);
        }
        else if(idir == 1)
        {
          //ebforallInPlace(numflopspt, "setFluxVal", setFluxVal,  faceBx,  *a_flux.m_yflux, fluxval);
          //using non-eb forall because box restriction in eb land is broken right now.   This will
          //work if there nare no cut cells near the domain boundary
          auto& regdata = a_flux.m_yflux->getRegData();
          forallInPlaceBase(setFluxVal, faceBx, regdata, fluxval);
        }
#if DIM==3          
        else if(idir == 2)
        {
          //ebforallInPlace(numflopspt, "setFluxVal", setFluxVal,  faceBx,  *a_flux.m_zflux, fluxval);
          //using non-eb forall because box restriction in eb land is broken right now.   This will
          //work if there nare no cut cells near the domain boundary
          auto& regdata = a_flux.m_zflux->getRegData();
          forallInPlaceBase(setFluxVal, faceBx, regdata, fluxval);
        }
#endif
        else
        {
          MayDay::Error("bogus idir");
        }
      }
    }
  }
}   
///
void 
EBMACProjector::
kappaDivU(EBLevelBoxData<CELL, 1> & a_divu,
          EBLevelFluxData<1>      & a_velo)
{

  CH_TIME("EBMACProjector::kappaDivU");
  a_velo.exchange(m_exchangeCopier);
  DataIterator dit = m_grids.dataIterator();
  int ideb = 0;
  for(int ibox = 0; ibox < dit.size(); ++ibox)
  {

    Bx   grid   =  ProtoCh::getProtoBox(m_grids[dit[ibox]]);
    const EBGraph  & graph = (*m_graphs)[dit[ibox]];
    Bx  grown   =  grid.grow(ProtoCh::getPoint(m_nghost));

    //get face fluxes and interpolate them to centroids
    EBFluxData<Real, 1>  centroidFlux(grown, graph);
    EBFluxStencil<2, Real> stencils =
      m_brit->getFluxStencil(StencilNames::InterpToFaceCentroid, StencilNames::NoBC, m_domain, m_domain, ibox);
    EBFluxData<Real,1>& faceCentFlux = a_velo[dit[ibox]];
    stencils.apply(centroidFlux, faceCentFlux, true, 1.0);  //true is to initialize to zero

    applyVeloBoundaryConditions(centroidFlux, dit[ibox]);

    auto& kapdiv =  m_rhs[dit[ibox]];
    kapdiv.setVal(0.);
    for(unsigned int idir = 0; idir < DIM; idir++)
    {
      bool initToZero = false;
      m_brit->applyFaceToCell(StencilNames::DivergeFtoC, StencilNames::NoBC, m_domain, kapdiv, centroidFlux,
                              idir, ibox, initToZero, 1.0);
    }
    ideb++;
  }
}
///
#include "Chombo_NamespaceFooter.H"


