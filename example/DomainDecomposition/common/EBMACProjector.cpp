
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
       const EBIBC                            & a_ebibc,
       bool a_printStuff)
{
  CH_TIME("EBMACProjector::define");
  if(a_printStuff)
  {
    pout() << "EBMACCProjector::define begin" << endl;
  }
  m_dx      = a_dx;
  m_ebibc   = a_ebibc;
  m_grids   = a_grids;
  m_domain  = a_domain;
  m_nghost  = a_nghost;
  m_grids   = a_grids;      
  m_domain  = a_domain;     
  m_brit    = a_brit;
  m_graphs  = a_geoserv->getGraphs(m_domain);

  if(a_printStuff)
  {
    pout() << "EBMACCProjector::define going into defineInternals" << endl;
  }

  defineInternals(a_geoserv, a_printStuff);

  if(a_printStuff)
  {
    pout() << "EBMACCProjector::define end" << endl;
  }
}
///
void 
EBMACProjector::
defineInternals(shared_ptr<GeometryService<2> >        & a_geoserv, bool a_printStuff)
{
  CH_TIME("EBMACProjector::defineInternals");

  if(a_printStuff)
  {
    pout() << "EBMACCProjector::defineInternals: making data holders" << endl;
  }
  
  m_rhs.define(m_grids, m_nghost, m_graphs);
  m_phi.define(m_grids, m_nghost, m_graphs);
  auto ditch = m_brit->m_cellToCell;
  Real alpha = 0; Real beta = 1; //Poisson's eqn

  string bcnames[2*DIM];
  m_ebibc.projectionStencilStrings(bcnames);

  if(a_printStuff)
  {
    pout() << "EBMACCProjector::defineInternals: defining multigrid solver" << endl;
  }
  
  string prefix("mac_proj");
  ParmParse pp(prefix.c_str());
  bool direct_to_bottom = false;
  string bottom_solver;
  bool useWCycle;
  int numSmooth;
  pp.get("direct_to_bottom", direct_to_bottom);
  pp.get("bottom_solver", bottom_solver);
  pp.get("useWCycle", useWCycle);

  pp.get("numSmooth", numSmooth);
  m_solver = shared_ptr<EBMultigrid>
    (new EBMultigrid(ditch, a_geoserv, alpha, beta, m_dx, m_grids, 
                     StencilNames::Poisson2, bcnames, StencilNames::Neumann,
                     m_domain, m_nghost,
                     bottom_solver, direct_to_bottom, prefix, useWCycle, numSmooth,
                     a_printStuff));

  if(a_printStuff)
  {
    pout() << "EBMACCProjector::defineInternals: calling registerStencils" << endl;
  }
  
  registerStencils();
  
  if(a_printStuff)
  {
    pout() << "EBMACCProjector::defineInternals: leaving" << endl;
  }
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
        Real a_tol, unsigned int a_maxiter, bool a_printStuff)
{
  CH_TIME("EBMACProjector::project");
  // set rhs = kappa*div (vel)
  if(a_printStuff)
  {
    pout() << "ebmacproj::project: going into kappaDivU" << endl;
  }

  kappaDivU(m_rhs, a_velo, a_printStuff);

  if(a_printStuff)
  {
    pout() <<"ebmacproj::project going into solve"<< endl;
  }

  //solve kappa*lapl(phi) = kappa*divu
  m_solver->solve(m_phi, m_rhs, a_tol, a_maxiter);
  

  if(a_printStuff)
  {
    pout() <<"ebmacproj::project computing and subtracting off gradient"<< endl;
  }
  //gphi = grad(phi)
  //v := v - gphi
  DataIterator dit = m_grids.dataIterator();
  for(int ibox = 0; ibox < dit.size(); ibox++)
  {
    
    Bx   grid   =  ProtoCh::getProtoBox(m_grids[dit[ibox]]);
    bool initToZero = true;
    m_brit->applyCellToFace(StencilNames::MACGradient, StencilNames::NoBC, m_domain,
                            a_gphi[dit[ibox]] ,m_phi[dit[ibox]], ibox, initToZero, 1.0);

    applyGradBoundaryConditions(a_gphi[dit[ibox]], dit[ibox]);
    a_velo[dit[ibox]] -= a_gphi[dit[ibox]];
  }
}
///Sets fluxes at domain edges (for inflow)
///Only sets one line of data (or plane if dim==3)
void 
EBMACProjector::
setFluxAtDomainBoundary(int                  a_idir,
                        Side::LoHiSide       a_sit,
                        EBFluxData<Real, 1>& a_flux,
                        Bx                   a_valid,  
                        Real                 a_fluxval,
                        Bx                   a_domain) //this funciton is static
{  
  //this is the comparison point in setValFluxLine.
  //the discrepancy in the two sides is due to face-centered
  //data.   The lowest face has an index of  0.  The
  //highest face has an index of nx.
  Proto::Point dompt;
  
  int usei = 0; int usej = 0; int usek = 0;
  int idom = 4586; int jdom = 4586; int kdom = 4586;
  if(a_sit == Side::Lo)
  {
    dompt = a_domain.low();
  }
  else
  {
    dompt = a_domain.high() + Point::Ones();
  }
  idom = dompt[0];
  jdom = dompt[1];
#if DIM==3
  kdom = dompt[2];
#endif  
  if(a_idir == 0)
  {
    usei = 1;
    idom = dompt[0];
    Bx inputBox = a_flux.m_xflux->inputBox();
    ebforall_i(inputBox,  setFluxValLine,  inputBox, *a_flux.m_xflux,
               usei, usej, usek, idom, jdom, kdom, a_fluxval);

  }
  else if(a_idir == 1)
  {
    usej = 1;
    idom = dompt[1];
    Bx inputBox = a_flux.m_yflux->inputBox();
    ebforall_i(inputBox,  setFluxValLine,  inputBox, *a_flux.m_yflux,
               usei, usej, usek, idom, jdom, kdom, a_fluxval);
  }
#if DIM==3
  else if(a_idir == 2)
  {
    usek = 1;
    idom = dompt[2];
    Bx inputBox = a_flux.m_zflux->inputBox();
    ebforall_i(inputBox,  setFluxValLine,  inputBox, *a_flux.m_zflux,
               usei, usej, usek, idom, jdom, kdom, a_fluxval);
  }
#endif
  else
  {
    Chombo4::MayDay::Error("setFaceValue: bogus idir");
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
          setFluxAtDomainBoundary(idir, sit(),  a_flux, valbx, fluxval, m_domain);
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

        if(bcstr != string("outflow"))
        {
          setFluxAtDomainBoundary(idir,
                                  sit(),
                                  a_flux,
                                  valbx,  
                                  fluxval,
                                  dombx); //this funciton is static
        }
      }
    }
  }
}   
///
void 
EBMACProjector::
kappaDivU(EBLevelBoxData<CELL, 1> & a_divu,
          EBLevelFluxData<1>      & a_velo,
          bool a_printStuff)
{

  CH_TIME("EBMACProjector::kappaDivU");
  if(a_printStuff)
  {
    pout() << "ebmacproj::kappaDivU: exchanging velocity input" << endl;
  }
  a_velo.exchange();

  if(a_printStuff)
  {
    pout() << "applying velocity boundary conditions and applying divergence stencils" << endl;
  }
  DataIterator dit = m_grids.dataIterator();
  int ideb = 0;
  for(int ibox = 0; ibox < dit.size(); ++ibox)
  {

    Bx   grid   =  ProtoCh::getProtoBox(m_grids[dit[ibox]]);
    const EBGraph  & graph = (*m_graphs)[dit[ibox]];
    Bx  grown   =  grid.grow(ProtoCh::getPoint(m_nghost));

    //get face fluxes and interpolate them to centroids
    bool useStack = true;
    EBFluxData<Real, 1>  centroidFlux(grown, graph, useStack);
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


