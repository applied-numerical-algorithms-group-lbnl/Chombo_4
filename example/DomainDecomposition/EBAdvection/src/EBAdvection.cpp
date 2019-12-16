#include "EBAdvection.H"
#include "EBAdvectionFunctions.H"
#include "Chombo_NamespaceHeader.H"
const string EBAdvection::s_ncdivLabel          = string("Volume_Weighted_Averaging_rad_1"); //this is for the non-conservative div
const string EBAdvection::s_aveCToFLabel        = string("AverageCellToFace"); //this is to get the velocity to faces
const string EBAdvection::s_nobcsLabel          = string("no_bcs"); //none of the operators here have eb boundary conditions
const string EBAdvection::s_redistLabel         = string("Volume_Weighted_Redistribution_rad_1"); //for redistribution
const string EBAdvection::s_centInterpLabel     = string("InterpolateToFaceCentroid"); //for interpolating from face centers to face centroids
const string EBAdvection::s_slopeLowLabel       = string("Slope_Low_");   //for low  side difference
const string EBAdvection::s_slopeHighLabel      = string("Slope_High_");  //for high side difference
const string EBAdvection::s_diriLabel           = string("Dirichlet");    //for diri bcs
const string EBAdvection::s_neumLabel           = string("Neumann");      //for neum bcs
const string EBAdvection::s_divergeLabel        = string("Divergence");   //for taking the divergence of face centered stuff to cell centered result
const string EBAdvection::s_CtoFLowLabel        = string("Cell_To_Face_Low");  //for getting stuff from low  side cells to faces
const string EBAdvection::s_CtoFHighLabel       = string("Cell_To_Face_High"); //for getting stuff from high side cells to faces
using Proto::Var;
////
EBAdvection::
EBAdvection(shared_ptr<EBEncyclopedia<2, Real> >   & a_brit,
            shared_ptr<GeometryService<2> >        & a_geoserv,
            shared_ptr<EBLevelBoxData<CELL, DIM> > & a_veloCell,
            const DisjointBoxLayout                & a_grids,
            const Box                              & a_domain,
            const Real                             & a_dx,
            const IntVect                          & a_nghostsrc, 
            const IntVect                          & a_nghostdst)
{
  m_dx     = a_dx;
  m_grids  = a_grids;
  m_domain = a_domain;
  m_nghostSrc  = a_nghostsrc;
  m_nghostDst  = a_nghostdst;
  m_exchangeCopier.exchangeDefine(m_grids, m_nghostSrc);
  m_grids      = a_grids;      
  m_domain     = a_domain;     
  m_brit       = a_brit;
  m_veloCell   = a_veloCell;
  defineData(a_geoserv);
  fillKappa( a_geoserv);
  registerStencils();
}
////
void  
EBAdvection::
defineData(shared_ptr<GeometryService<2> >        & a_geoserv)
{
  m_graphs = a_geoserv->getGraphs(m_domain);
  m_kappa.define(     m_grids, m_nghostSrc, m_graphs);
  m_deltaM.define(    m_grids, m_nghostSrc, m_graphs);
  m_nonConsDiv.define(m_grids, m_nghostSrc, m_graphs);
  m_hybridDiv.define( m_grids, m_nghostSrc, m_graphs);
  m_kappaDiv.define(  m_grids, m_nghostSrc, m_graphs);

  m_exchangeCopier.exchangeDefine(m_grids, m_nghostSrc);
}
////
void  
EBAdvection::
fillKappa(shared_ptr<GeometryService<2> >        & a_geoserv)
{
  DataIterator dit = m_grids.dataIterator();
  for(int ibox = 0; ibox < dit.size(); ++ibox)
  {
    Box grid =m_grids[dit[ibox]];
    Bx  grbx = getProtoBox(grid);
    const EBGraph  & graph = (*m_graphs)[dit[ibox]];
    EBBoxData <CELL, Real, 1> & devidat = m_kappa[dit[ibox]];    
    EBHostData<CELL, Real, 1>   hostdat(devidat.box(), graph);
    //fill kappa on the host then copy to the device
    a_geoserv->fillKappa(hostdat, grid, dit[ibox], m_domain);
    // now copy to the device
    EBLevelBoxData<CELL, 1>::copyToDevice(hostdat, m_kappa[dit[ibox]]);
  }
  m_kappa.exchange(m_exchangeCopier);
}
////
void  
EBAdvection::
registerStencils()
{

  //false is because I do not need diagonal  weights for any of these stencils
  bool needDiag = false;
  m_brit->m_cellToCell->registerStencil(s_ncdivLabel , s_nobcsLabel, s_nobcsLabel, m_domain, m_domain, needDiag);
  m_brit->m_cellToCell->registerStencil(s_redistLabel, s_nobcsLabel, s_nobcsLabel, m_domain, m_domain, needDiag);

  m_brit->registerFaceStencil(s_centInterpLabel, s_nobcsLabel, s_nobcsLabel, m_domain, m_domain, needDiag);
  //no flow means dirichlet boundary conditions for normal velocities
#if 1
  //real code
  m_brit->registerCellToFace( s_aveCToFLabel , s_diriLabel , s_nobcsLabel, m_domain, m_domain, needDiag, Point::Ones());
#else
  //debug --neumann so I can use all reg
  m_brit->registerCellToFace( s_aveCToFLabel , s_neumLabel , s_nobcsLabel, m_domain, m_domain, needDiag, Point::Ones());
  //end debug
#endif

  m_brit->registerCellToFace( s_CtoFHighLabel, s_nobcsLabel, s_nobcsLabel, m_domain, m_domain, needDiag, Point::Ones());
  m_brit->registerCellToFace( s_CtoFLowLabel , s_nobcsLabel, s_nobcsLabel, m_domain, m_domain, needDiag, Point::Ones());
  m_brit->registerFaceToCell( s_divergeLabel , s_nobcsLabel, s_nobcsLabel, m_domain, m_domain, needDiag);
  for(int idir = 0; idir < DIM; idir++)
  {
    //need to set neumann bcs to set slopes to zero at domain bcs.
    string slopeLow   = s_slopeLowLabel  + std::to_string(idir);
    string slopeHigh  = s_slopeHighLabel + std::to_string(idir);
    m_brit->m_cellToCell->registerStencil(slopeLow  , s_neumLabel, s_nobcsLabel, m_domain, m_domain, needDiag, Point::Ones());
    m_brit->m_cellToCell->registerStencil(slopeHigh , s_neumLabel, s_nobcsLabel, m_domain, m_domain, needDiag, Point::Ones());
  }

}

///
void
EBAdvection::
getFaceCenteredVel(EBFluxData<Real, 1>& a_fcvel,
                   const DataIndex    & a_dit,
                   const int          & a_ibox)
{
  EBBoxData<CELL, Real, DIM>& veccell = (*m_veloCell)[a_dit];
  getFaceVelComp<XFACE>(*(a_fcvel.m_xflux),  m_brit->m_cellToXFace, veccell, 0, a_ibox);
  getFaceVelComp<YFACE>(*(a_fcvel.m_yflux),  m_brit->m_cellToYFace, veccell, 1, a_ibox);
#if DIM==3 
  getFaceVelComp<ZFACE>(*(a_fcvel.m_zflux),  m_brit->m_cellToZFace, veccell, 2, a_ibox);
#endif
}

///
void
EBAdvection::
getFaceCenteredFlux(EBFluxData<Real, 1>            & a_fcflux,
                    const EBFluxData<Real, 1>      & a_fcvel,
                    const EBBoxData<CELL, Real, 1> & a_scal,
                    const DataIndex                & a_dit,
                    const int                      & a_ibox,
                    const Real                     & a_dt)
{
  //first we compute the slopes of the data
  //then we extrapolate in space and time
  //then we solve the riemann problem to get the flux
  Bx   grid   =  ProtoCh::getProtoBox(m_grids[a_dit]);
  Bx  grown   =  grid.grow(ProtoCh::getPoint(m_nghostSrc));
  const EBGraph  & graph = (*m_graphs)[a_dit];
  EBBoxData<CELL, Real, DIM>& veccell = (*m_veloCell)[a_dit];
  bool initToZero = true;
  //compute slopes of the solution 
  //(low and high centered) in each direction

  EBBoxData<CELL, Real, DIM> slopeLo(grown, graph); 
  EBBoxData<CELL, Real, DIM> slopeHi(grown, graph);
  int ideb = 0;
  for(unsigned int idir = 0; idir < DIM; idir++)
  {
    EBBoxData<CELL, Real, 1> slopeLoComp, slopeHiComp;
    slopeLoComp.define<DIM>(slopeLo, idir);
    slopeHiComp.define<DIM>(slopeHi, idir);

    //low and high side differences
    string slopeLoLab  = s_slopeLowLabel  + std::to_string(idir);
    string slopeHiLab  = s_slopeHighLabel + std::to_string(idir);
    const auto& stenlo = m_brit->m_cellToCell->getEBStencil(slopeLoLab , s_nobcsLabel, m_domain, m_domain, a_ibox);
    const auto& stenhi = m_brit->m_cellToCell->getEBStencil(slopeHiLab , s_nobcsLabel, m_domain, m_domain, a_ibox);

    stenlo->apply(slopeLoComp, a_scal, initToZero, 1.0);
    stenhi->apply(slopeHiComp, a_scal, initToZero, 1.0);
    ideb++;
  }
//debug
//  slopeLo.setVal(0.);
//  slopeHi.setVal(0.);
//end debug

  EBFluxData<Real, 1>  scalHi(grown, graph);
  EBFluxData<Real, 1>  scalLo(grown, graph);
  for(unsigned int idir = 0; idir < DIM; idir++)
  {
    //scalar extrapolated to low side and high side face
    EBBoxData<CELL, Real, 1> scal_imh_nph(grown, graph);
    EBBoxData<CELL, Real, 1> scal_iph_nph(grown, graph);
    //extrapolate in space and time to get the inputs to the Riemann problem
    unsigned long long int numflopspt = 19 + 4*DIM;

#if 1
  //begin debug  call point version
    ebforallInPlace_i(numflopspt, "ExtrapolateScal", ExtrapolateScalPt, grown,  
                    scal_imh_nph, scal_iph_nph, a_scal, 
                    slopeLo, slopeHi, veccell, idir, a_dt, m_dx);
  //end debug;
#else
    ebforallInPlace(numflopspt, "ExtrapolateScal", ExtrapolateScal, grown,  
                    scal_imh_nph, scal_iph_nph, a_scal, 
                    slopeLo, slopeHi, veccell, idir, a_dt, m_dx);
#endif


    //we need to get the low and high states from the cell-centered holders to the face centered ones.
    //once we do that, we can solve the Rieman problem for the upwind state
    //i + 1/2 becomes the low  side of the face
    //i - 1/2 becomes the high side of the face
    m_brit->applyCellToFace(s_CtoFHighLabel, s_nobcsLabel, m_domain, scalHi, scal_imh_nph, idir, a_ibox, initToZero, 1.0);
    m_brit->applyCellToFace(s_CtoFLowLabel , s_nobcsLabel, m_domain, scalLo, scal_iph_nph, idir, a_ibox, initToZero, 1.0);
    ideb++;
  }

  //this solves the Riemann problem and sets flux = facevel*(upwind scal)
  unsigned long long int numflopspt = 2;
#if 1
  //begin debug  call point version
  ebforallInPlace_i(numflopspt, "GetUpwindFlux", GetUpwindFluxPt, a_fcflux.m_xflux->box(),
                  *a_fcflux.m_xflux, *scalLo.m_xflux, *scalHi.m_xflux, *a_fcvel.m_xflux);
  //end   debug 
  
#else
  ebforallInPlace(numflopspt, "GetUpwindFlux", GetUpwindFlux, a_fcflux.m_xflux->box(),
                  *a_fcflux.m_xflux, *scalLo.m_xflux, *scalHi.m_xflux, *a_fcvel.m_xflux);
#endif

#if 1
  //begin debug  call point version
  ebforallInPlace_i(numflopspt, "GetUpwindFlux", GetUpwindFluxPt, a_fcflux.m_yflux->box(),
                  *a_fcflux.m_yflux, *scalLo.m_yflux, *scalHi.m_yflux, *a_fcvel.m_yflux);
  //end   debug 
  
#else
  ebforallInPlace(numflopspt, "GetUpwindFlux", GetUpwindFlux, a_fcflux.m_yflux->box(),
                  *a_fcflux.m_yflux, *scalLo.m_yflux, *scalHi.m_yflux, *a_fcvel.m_yflux);
#endif

#if DIM==3
  ebforallInPlace(numflopspt, "GetUpwindFlux", GetUpwindFlux, a_fcflux.m_zflux->box(),
                  *a_fcflux.m_zflux, *scalLo.m_zflux, *scalHi.m_zflux, *a_fcvel.m_zflux);
#endif
  
}
                  
///
void
EBAdvection::
kappaConsDiv(EBLevelBoxData<CELL, 1>   & a_scal, const Real& a_dt)
{
  //coming into this we have the scalar at time = n dt
  // velocity field at cell centers. Leaving, we have filled
  // kappa* div(u scal)
  DataIterator dit = m_grids.dataIterator();
  int ideb = 0;
  for(int ibox = 0; ibox < dit.size(); ++ibox)
  {
    Bx   grid   =  ProtoCh::getProtoBox(m_grids[dit[ibox]]);
    Bx  grown   =  grid.grow(ProtoCh::getPoint(m_nghostSrc));

    const EBGraph  & graph = (*m_graphs)[dit[ibox]];
    //get face fluxes and interpolate them to centroids
    EBFluxData<Real, 1>  centroidFlux(grown, graph);
    EBFluxData<Real, 1>  faceCentFlux(grown, graph);
    EBFluxData<Real, 1>  faceCentVel( grown, graph);
    //average velocities to face centers.
    getFaceCenteredVel( faceCentVel, dit[ibox], ibox);


    getFaceCenteredFlux(faceCentFlux, faceCentVel, a_scal[dit[ibox]], dit[ibox], ibox, a_dt);
    EBFluxStencil<2, Real> stencils =   m_brit->getFluxStencil(s_centInterpLabel, s_nobcsLabel, m_domain, m_domain, ibox);
    //interpolate flux to centroids

    stencils.apply(centroidFlux, faceCentFlux, true, 1.0);  //true is to initialize to zero

    auto& kapdiv =  m_kappaDiv[dit[ibox]];
    kapdiv.setVal(0.);
//    for(unsigned int idir = DIM-1; idir >= 0; idir--)
    for(unsigned int idir = 0; idir < DIM; idir++)
    {
      bool initToZero = false;
      m_brit->applyFaceToCell(s_divergeLabel , s_nobcsLabel, m_domain, kapdiv, centroidFlux,
                              idir, ibox, initToZero, 1.0);
    }
    ideb++;
  }

}
///
void
EBAdvection::
nonConsDiv()
{
  //this makes ncdiv = divF on regular cells
  //and ncdiv = vol_weighted_ave(div) on cut cells
  DataIterator dit = m_grids.dataIterator();
  int ideb = 0;
  for(int ibox = 0; ibox < dit.size(); ++ibox)
  {
    auto& ncdiv  = m_nonConsDiv[dit[ibox]];
    auto& kapdiv =   m_kappaDiv[dit[ibox]];
    //auto& kappa =       m_kappa[dit[ibox]];
    const auto& stencil = m_brit->m_cellToCell->getEBStencil(s_ncdivLabel,s_nobcsLabel, m_domain, m_domain, ibox);
    bool initToZero = true;
    stencil->apply(ncdiv, kapdiv, initToZero, 1.0);
    ideb++;
  }
//begin debug
//  Real coveredval = -0.125;
//  Real time = 0;
//  Real dt = 1;
//  {
//    string filep("divnc.hdf5");
//    writeEBLevelHDF5<1>(  filep,  m_nonConsDiv, m_kappa, m_domain, m_graphs, coveredval, m_dx, dt, time);
//  }
//  {
//    string filep("kapdiv.hdf5");
//    writeEBLevelHDF5<1>(  filep,  m_kappaDiv, m_kappa, m_domain, m_graphs, coveredval, m_dx, dt, time);
//  }
//end debug
}
///
void
EBAdvection::
redistribute()
{
  //hybrid div comes in holding kappa*div^c + (1-kappa)div^nc
  //this redistributes delta M into the hybrid divergence.
  DataIterator dit = m_grids.dataIterator();
  for(int ibox = 0; ibox < dit.size(); ++ibox)
  {
    auto& deltaM = m_deltaM   [dit[ibox]];
    auto& hybrid = m_hybridDiv[dit[ibox]];
    const auto& stencil = m_brit->m_cellToCell->getEBStencil(s_redistLabel,s_nobcsLabel, m_domain, m_domain, ibox);
    bool initToZero = false;
    stencil->apply(hybrid, deltaM, initToZero, 1.0);
  }
}
///
void 
EBAdvection::
advance(EBLevelBoxData<CELL, 1>       & a_phi,
        const  Real                   & a_dt)
{
  a_phi.exchange(m_exchangeCopier);
  //compute kappa div^c F
  kappaConsDiv(a_phi, a_dt);

  //compute nonconservative divergence = volume weighted ave of div^c
  nonConsDiv();

  //advance solution, compute delta M
  DataIterator dit = m_grids.dataIterator();
  int ideb = 0;
  for(int ibox = 0; ibox < dit.size(); ibox++)
  {
    unsigned long long int numflopspt = 7;
    Bx grbx = ProtoCh::getProtoBox(m_grids[dit[ibox]]);
    auto & hybridDiv  =  m_hybridDiv[ dit[ibox]];
    auto & kappaDiv   =  m_kappaDiv[  dit[ibox]];
    auto & nonConsDiv =  m_nonConsDiv[dit[ibox]];  
    auto & deltaM     =  m_deltaM[    dit[ibox]]; 
    auto & kappa      =  m_kappa[     dit[ibox]];
#if 0
    ebforallInPlace(numflopspt, "HybridDivergence", HybridDivergence, grbx,  
                    hybridDiv ,
                    kappaDiv  ,
                    nonConsDiv,  
                    deltaM    , 
                    kappa     );
#else
    //debugging
    ebforallInPlace_i(numflopspt, "HybridDivergencePt", HybridDivergencePt, grbx,  
      hybridDiv ,
      kappaDiv  ,
      nonConsDiv,  
      deltaM    , 
      kappa     );
#endif 
    ideb++;
  }
  m_deltaM.exchange(m_exchangeCopier);
  //redistribute delta M
  //debug turn off redist
  //redistribute();
  //end debug
  for(int ibox = 0; ibox < dit.size(); ibox++)
  {
    auto& scalar =       a_phi[dit[ibox]];
    auto& diverg = m_hybridDiv[dit[ibox]];
    Bx grbx = ProtoCh::getProtoBox(m_grids[dit[ibox]]);
    unsigned long long int numflopspt = 2;
    ebforallInPlace(numflopspt, "AdvanceScalar", AdvanceScalar,  grbx,  
      scalar, diverg, a_dt);

    ideb++;
  }
    
}
#include "Chombo_NamespaceFooter.H"

