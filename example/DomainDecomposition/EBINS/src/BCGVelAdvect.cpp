
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

      getUpwindState(upwindScal, faceCentVelo, scalLo, scalHi);
      EBFluxData<Real, 1>& advvelfab = m_advectionVel[dit[ibox]];
      
      //now copy into the normal direction holder
      copyComp(advvelfab, upwindScal, idir);

    } //end loop over boxes
  }  // and loop over velocity directions
  
  // now we need to project the mac velocity
  m_macproj->project(m_advectionVel, m_macGradient, m_tol, m_maxIter);
}
/*******/

//solve for the upwind value
void  
BCGVelAdvect::      
getUpwindState(EBFluxData<Real, 1>&  a_upwindScal,
               EBFluxData<Real, 1>&  a_faceCentVelo,
               EBFluxData<Real, 1>&  a_scalLo,
               EBFluxData<Real, 1>&  a_scalHi)
{
  unsigned long long int numflopspt = 0;
  ebforallInPlace(numflopspt, "Upwinded", Upwinded, a_upwindScal.m_xflux->box(),
                  *a_upwindScal.m_xflux, *a_scalLo.m_xflux, *a_scalHi.m_xflux,
                  *a_faceCentVelo.m_xflux);

  ebforallInPlace(numflopspt, "Upwinded", Upwinded, a_upwindScal.m_yflux->box(),
                  *a_upwindScal.m_yflux, *a_scalLo.m_yflux, *a_scalHi.m_yflux,
                  *a_faceCentVelo.m_yflux);

#if DIM==3
  ebforallInPlace(numflopspt, "Upwinded", Upwinded, a_upwindScal.m_zflux->box(),
                  *a_upwindScal.m_zflux, *a_scalLo.m_zflux, *a_scalHi.m_zflux,
                  *a_faceCentVelo.m_zflux);
#endif
}
/*******/

void  
BCGVelAdvect::      
copyComp(EBFluxData<Real, 1>&  a_dst,
         EBFluxData<Real, 1>&  a_src,
         int a_vecDir)
{
  unsigned numflopspt = 0;
  if(a_vecDir == 0)
  {
    ebforallInPlace(numflopspt, "Copied", Copied, a_dst.m_xflux->box(),
                    *a_dst.m_xflux, *a_src.m_xflux);
  }
  else if(a_vecDir == 1)
  {
    ebforallInPlace(numflopspt, "Copied", Copied, a_dst.m_yflux->box(),
                    *a_dst.m_yflux, *a_src.m_yflux );
  }
#if DIM==3
  else if(a_vecDir == 2)
  {
    ebforallInPlace(numflopspt, "Copied", Copied, a_dst.m_zflux->box(),
                    *a_dst.m_zflux, *a_src.m_zflux );
  }
#endif
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
      getUpwindState(upwindScal, faceCentVelo, scalLo, scalHi);

    } //end loop over boxes
  }  // and loop over velocity directions

  //subtract off pure gradient part of the velocity field
  correctVectorVelocity();
}
/*******/
//subtract off pure gradient part of the velocity field
//makes normal component of the velocity  = advection velocity at the face
//Recall the gradients only exist on the normal face.
//The tangential compoents get the average of neighboring gradients in the same direciton
//subtracted.   So the y component on the xfaces get the average of the values of the
//gradient on the neighboring y faces
void 
BCGVelAdvect::
correctVectorVelocity()
{
  for(unsigned int vecDir = 0; vecDir < DIM; vecDir++)
  {
    EBLevelFluxData<1> facecomp;
    facecomp.define(m_macVelocity, vecDir, m_graphs);

    DataIterator dit = m_grids.dataIterator();
    for(int ibox = 0; ibox < dit.size(); ++ibox)
    {
      EBFluxData<Real, 1>& advvelfab = m_advectionVel[dit[ibox]];
      EBFluxData<Real, 1>& veccomp   = m_advectionVel[dit[ibox]];
      EBFluxData<Real, 1>& macGrad   =  m_macGradient[dit[ibox]];
      //first correct facedir==vecdir direction
      //if the two directions are the same, substitute the advection velocity, which already has its gradient removed
      copyComp(veccomp, advvelfab, vecDir);

      for(unsigned int faceDir = 0; faceDir < DIM; faceDir++)
      {
        if(faceDir != vecDir)
        {
          EBCrossFaceStencil<2, Real> stencils =
            m_brit->getCrossFaceStencil(StencilNames::TanVelCorrect, StencilNames::NoBC, m_domain, m_domain, ibox);
          stencils.apply(veccomp, macGrad, faceDir, vecDir);
        }
      }
    } // end loop over facedir
  }
}

/*******/
void 
BCGVelAdvect::
assembleDivergence()
{
}
/*******/
#include "Chombo_NamespaceFooter.H"

