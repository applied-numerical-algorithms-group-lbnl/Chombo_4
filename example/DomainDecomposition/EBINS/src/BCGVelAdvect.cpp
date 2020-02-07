
#include "BCGVelAdvect.H"
#include "EBAdvectionFunctions.H"
#include "Chombo_NamespaceHeader.H"
/*******/
void 
BCGVelAdvect::
hybridDivergence(EBLevelBoxData<CELL, DIM>& a_divuu,
                 EBLevelBoxData<CELL, DIM>& a_inputVel,
                 const Real               & a_dt)
{
  //fills m_macScal with velocity extrapolated to its normal face and projected
  getAdvectionVelocity(a_inputVel, a_dt);

  //using macScal for upwinding, this gets the vector compoents of velcity to faces
  //this also gets corrected by teh gradient found in the previous step.
  getMACVectorVelocity(a_inputVel, a_dt);

  //lots of aliasing and smushing
  assembleDivergence();
}
/*******/
void  
BCGVelAdvect::
getAdvectionVelocity(EBLevelBoxData<CELL, DIM>   & a_inputVel,
                     const Real                  & a_dt)
{
  for(unsigned int idir = 0; idir < DIM; idir++)
  {
    EBLevelBoxData<CELL, 1> velcomp;
    velcomp.define(a_inputVel, idir, m_graphs);
    //source term = nu*lapl(ucomp);
    Real alpha = 0; Real beta = m_viscosity;
    m_helmholtz->resetAlphaAndBeta(alpha, beta);
    m_helmholtz->applyOp(m_source, velcomp);

    DataIterator dit = m_grids.dataIterator();
    for(int ibox = 0; ibox < dit.size(); ++ibox)
    {
      Bx   grid   =  ProtoCh::getProtoBox(m_grids[dit[ibox]]);
      Bx  grown   =  grid.grow(ProtoCh::getPoint(m_nghostSrc));
      const EBGraph  & graph = (*m_graphs)[dit[ibox]];

      //this gets the low and high side states for the riemann problem
      auto& scalfab = velcomp[dit[ibox]];
      EBFluxData<Real, 1>  scalHi(grown, graph);
      EBFluxData<Real, 1>  scalLo(grown, graph);
      auto & veccell = a_inputVel[dit[ibox]];
      auto & sourfab =   m_source[dit[ibox]];
      bcgExtrapolateScalar(scalLo, scalHi, veccell, scalfab, sourfab,
                           grown, graph, dit[ibox], ibox, a_dt);

      //get face fluxes and interpolate them to centroids
      EBFluxData<Real, 1>  faceCentVelo( grown, graph);
      EBFluxData<Real, 1>  upwindScal(   grown, graph);
      //average velocities to face centers.
      getFaceCenteredVel( faceCentVelo, dit[ibox], ibox);

      //solve for the upwind value
      unsigned long long int numflopspt = 2;
      ebforallInPlace(numflopspt, "Upwinded", Upwinded, upwindScal.m_xflux->box(),
                      *upwindScal.m_xflux, *scalLo.m_xflux, *scalHi.m_xflux,
                      *faceCentVelo.m_xflux);

      ebforallInPlace(numflopspt, "Upwinded", Upwinded, upwindScal.m_yflux->box(),
                      *upwindScal.m_yflux, *scalLo.m_yflux, *scalHi.m_yflux,
                      *faceCentVelo.m_yflux);

#if DIM==3
      ebforallInPlace(numflopspt, "Upwinded", Upwinded, upwindScal.m_zflux->box(),
                      *upwindScal.m_zflux, *scalLo.m_zflux, *scalHi.m_zflux,
                      *faceCentVelo.m_zflux);
#endif
      EBFluxData<Real, 1>& advvelfab = m_advectionVel[dit[ibox]];
      //now copy into the normal direction holder
      for(unsigned int faceDir = 0; faceDir < DIM; ++faceDir)
      {
        if(faceDir == idir)
        {
          if(idir == 0)
          {
            ebforallInPlace(numflopspt, "Copied", Copied, upwindScal.m_xflux->box(),
                            *advvelfab.m_xflux, *upwindScal.m_xflux);
          }
          else if(idir == 1)
          {
            ebforallInPlace(numflopspt, "Copied", Copied, upwindScal.m_yflux->box(),
                            *advvelfab.m_yflux, *upwindScal.m_yflux);
          }
#if DIM==3
          else if(idir == 2)
          {
            ebforallInPlace(numflopspt, "Copied", Copied, upwindScal.m_zflux->box(),
                            *advvelfab.m_zflux, *upwindScal.m_zflux);
          }
#endif
        } // end if we are on normal face
      } // end loop over facedir

    } //end loop over boxes
  }  // and loop over velocity directions
  
  // now we need to project the mac velocity
  m_macproj->project(m_advectionVel, m_macGradient, m_tol, m_maxIter);
}
/*******/
void  
BCGVelAdvect::
getMACVectorVelocity(EBLevelBoxData<CELL, DIM>   & a_inputVel,
                     const Real                  & a_dt)
{
  for(unsigned int vecDir = 0; vecDir < DIM; vecDir++)
  {
    EBLevelBoxData< CELL, 1> velcomp;
    EBLevelFluxData<1> facecomp;
    velcomp.define(a_inputVel, vecDir, m_graphs);
    facecomp.define(m_macVelocity, vecDir, m_graphs);

    //source term = nu*lapl(ucomp);
    Real alpha = 0; Real beta = m_viscosity;
    m_helmholtz->resetAlphaAndBeta(alpha, beta);
    m_helmholtz->applyOp(m_source, velcomp);

    DataIterator dit = m_grids.dataIterator();
    for(int ibox = 0; ibox < dit.size(); ++ibox)
    {
      Bx   grid   =  ProtoCh::getProtoBox(m_grids[dit[ibox]]);
      Bx  grown   =  grid.grow(ProtoCh::getPoint(m_nghostSrc));
      const EBGraph  & graph = (*m_graphs)[dit[ibox]];

      //this gets the low and high side states for the riemann problem
      auto& scalfab = velcomp[dit[ibox]];
      EBFluxData<Real, 1>  scalHi(grown, graph);
      EBFluxData<Real, 1>  scalLo(grown, graph);
      auto & veccell = a_inputVel[dit[ibox]];
      auto & sourfab =   m_source[dit[ibox]];
      bcgExtrapolateScalar(scalLo, scalHi, veccell, scalfab, sourfab,
                           grown, graph, dit[ibox], ibox, a_dt);

      //get face fluxes and interpolate them to centroids
      EBFluxData<Real, 1>&  faceCentVelo = m_advectionVel[dit[ibox]];
      EBFluxData<Real, 1>&  upwindScal   =       facecomp[dit[ibox]];

      //solve for the upwind value
      unsigned long long int numflopspt = 2;
      ebforallInPlace(numflopspt, "Upwinded", Upwinded, upwindScal.m_xflux->box(),
                      *upwindScal.m_xflux, *scalLo.m_xflux, *scalHi.m_xflux,
                      *faceCentVelo.m_xflux);

      ebforallInPlace(numflopspt, "Upwinded", Upwinded, upwindScal.m_yflux->box(),
                      *upwindScal.m_yflux, *scalLo.m_yflux, *scalHi.m_yflux,
                      *faceCentVelo.m_yflux);

#if DIM==3
      ebforallInPlace(numflopspt, "Upwinded", Upwinded, upwindScal.m_zflux->box(),
                      *upwindScal.m_zflux, *scalLo.m_zflux, *scalHi.m_zflux,
                      *faceCentVelo.m_zflux);
#endif

    } //end loop over boxes
  }  // and loop over velocity directions
}
/*******/
void 
BCGVelAdvect::
assembleDivergence()
{
}
/*******/
#include "Chombo_NamespaceFooter.H"

