
#include "BCGVelAdvect.H"
#include "EBAdvectionFunctions.H"
#include "Chombo_NamespaceHeader.H"
/*******/
void 
BCGVelAdvect::
hybridDivergence(EBLevelBoxData<CELL, DIM>& a_divuu,
                 EBLevelBoxData<CELL, DIM>& a_inputVel,
                 const Real               & a_dt,
                 Real a_tolerance, unsigned int a_maxIter)    
{
  //fills m_macScal with velocity extrapolated to its normal face and projected
  getAdvectionVelocity(a_inputVel, a_dt, a_tolerance, a_maxIter);

  //using macScal for upwinding, this gets the vector compoents of velcity to faces
  //this also gets corrected by teh gradient found in the previous step.
  getMACVectorVelocity(a_inputVel, a_dt);

  //lots of aliasing and smushing
  assembleDivergence(a_divuu, a_dt);
}
/*******/
void  
BCGVelAdvect::
getAdvectionVelocity(EBLevelBoxData<CELL, DIM>   & a_inputVel,
                     const Real                  & a_dt,
                     Real a_tol, unsigned int a_maxIter)    
{
  a_inputVel.exchange(m_exchangeCopier);
  for(unsigned int idir = 0; idir < DIM; idir++)
  {
    EBLevelBoxData<CELL, 1> velcomp;
    velcomp.define<DIM>(a_inputVel, idir, m_graphs);
    //source term = nu*lapl(ucomp);
    if(m_eulerCalc)
    {
      m_source.setVal(0.);
    }
    else
    {
      Real alpha = 0; Real beta = m_viscosity;
      m_helmholtz->resetAlphaAndBeta(alpha, beta);
      m_helmholtz->applyOp(m_source, velcomp);
    }

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
  pout() << "mac projecting advection velocity" << endl;
  m_macproj->project(m_advectionVel, m_macGradient, a_tol, a_maxIter);
  m_advectionVel.exchange(m_exchangeCopier);
  pout() << "leaving getAdvectionVelocity" << endl;
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
    velcomp.define<DIM> (a_inputVel, vecDir, m_graphs);
    facecomp.define<DIM>(m_macVelocity, vecDir, m_graphs);

    //source term = nu*lapl(ucomp);
    if(m_eulerCalc)
    {
      m_source.setVal(0.); 
    }
    else
    {
      Real alpha = 0; Real beta = m_viscosity;
      m_helmholtz->resetAlphaAndBeta(alpha, beta);
      m_helmholtz->applyOp(m_source, velcomp);
    }

    m_source.exchange(m_exchangeCopier);
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
  m_macGradient.exchange(m_exchangeCopier);
  for(unsigned int vecDir = 0; vecDir < DIM; vecDir++)
  {
    EBLevelFluxData<1> facecomp;
    facecomp.define<DIM>(m_macVelocity, vecDir, m_graphs);

    DataIterator dit = m_grids.dataIterator();
    for(int ibox = 0; ibox < dit.size(); ++ibox)
    {
      auto & advvelfab = m_advectionVel[dit[ibox]];
      auto & veccomp   =       facecomp[dit[ibox]];
      auto & macGrad   =  m_macGradient[dit[ibox]];
      //first correct facedir==vecdir direction
      //if the two directions are the same, substitute the advection velocity, which already has its gradient removed
      copyComp(veccomp, advvelfab, vecDir);

      for(unsigned int faceDir = 0; faceDir < DIM; faceDir++)
      {
        if(faceDir != vecDir)
        {
          EBCrossFaceStencil<2, Real> stencils =
            m_brit->getCrossFaceStencil(StencilNames::TanVelCorrect, StencilNames::NoBC, m_domain, m_domain, ibox);
          bool initToZero = false;
          Real scale = -1.0; //subtracting off gradient
          stencils.apply(veccomp, macGrad, faceDir, vecDir, initToZero, scale);
        }
      }
    } // end loop over facedir
  }
}

/**********/
void 
BCGVelAdvect::
assembleDivergence(EBLevelBoxData<CELL, DIM>& a_divuu, 
                   const  Real              & a_dt)
{
  m_macVelocity.exchange(m_exchangeCopier);
  DataIterator dit = m_grids.dataIterator();
  for(int ivar = 0; ivar < DIM; ivar++)
  {
    EBLevelBoxData<CELL, 1> hybridDiv;
    EBLevelFluxData<1> scalarVelComp;
    hybridDiv.define<DIM>(          a_divuu, ivar, m_graphs);
    scalarVelComp.define<DIM>(m_macVelocity, ivar, m_graphs);
    //this loop fills m_kappaDiv with the conservative divergence
    for(unsigned int ibox = 0; ibox < dit.size(); ibox++)
    {
      Bx   grid   =  ProtoCh::getProtoBox(m_grids[dit[ibox]]);
      Bx  grown   =  grid.grow(ProtoCh::getPoint(m_nghostSrc));
      const EBGraph  & graph = (*m_graphs)[dit[ibox]];

      auto& faceCentVel  = m_advectionVel[dit[ibox]];
      auto& scalar       =  scalarVelComp[dit[ibox]];

      EBFluxData<Real, 1>  centroidFlux(grown, graph);
      EBFluxData<Real, 1>  faceCentFlux(grown, graph);

      assembleFlux(faceCentFlux, scalar, faceCentVel);

      EBFluxStencil<2, Real> stencils =   m_brit->getFluxStencil(s_centInterpLabel, s_nobcsLabel, m_domain, m_domain, ibox);
      //interpolate flux to centroids
      stencils.apply(centroidFlux, faceCentFlux, true, 1.0);  //true is to initialize to zero
      scalar.define(m_macVelocity[dit[ibox]], ivar);

      auto& kapdiv =  m_kappaDiv[dit[ibox]];
      getKapDivFFromCentroidFlux(kapdiv, centroidFlux, ibox);
    }
    m_kappaDiv.exchange(m_exchangeCopier);
    //this computes the non-conservative divergence and puts it into m_nonConsDiv
    nonConsDiv();
    //this forms the hybrid divergence

    //does linear combination of divnc and div c to get hybrid and  compute delta M
    kappaDivPlusOneMinKapDivNC(hybridDiv);

//    redistribute(hybridDiv);
  }
}
/*******/
#include "Chombo_NamespaceFooter.H"

