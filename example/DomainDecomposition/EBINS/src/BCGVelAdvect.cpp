
#include "BCGVelAdvect.H"
#include "EBAdvectionFunctions.H"
#include "Chombo_ParmParse.H"
#include "Chombo_NamespaceHeader.H"
///
void
BCGVelAdvect::
applyVeloFluxBCs(EBFluxData<Real, 1> & a_flux,
                 const DataIndex     & a_dit,
                 EBFluxData<Real, 1> & a_scalarLo,
                 EBFluxData<Real, 1> & a_scalarHi,
                 unsigned int a_velcomp)
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
        if(bcstr == string("outflow"))
        {
          //copy correct side of extrapolation to flux  (using extrapolated value as flux here)
          copyExtrapolatedState(a_flux, a_scalarLo, a_scalarHi, idir, sit(),  valbx);
        }
        else if(bcstr == string("inflow"))
        {
          if(idir == a_velcomp)
          {
            Real fluxval = 0;
            ParmParse pp;
            pp.get("velocity_inflow_value", fluxval);
            EBMACProjector::setFluxAtDomainBoundary(idir, sit(), a_flux, valbx, fluxval, m_domain);
          }
          else
          {
            //copy correct side of extrapolation to flux  (using extrapolated value as flux here)
            copyExtrapolatedState(a_flux, a_scalarLo, a_scalarHi, idir, sit(),  valbx);
          }
        }
        else if((bcstr == string("no_slip_wall")) || (bcstr == string("slip_wall")))
        {
          if(a_velcomp== idir)
          {
            Real fluxval = 0;
            EBMACProjector::setFluxAtDomainBoundary(idir, sit(), a_flux, valbx, fluxval, m_domain);
          }
          else
          {
            //copy correct side of extrapolation to flux  (using extrapolated value as flux here)
            copyExtrapolatedState(a_flux, a_scalarLo, a_scalarHi, idir, sit(),  valbx);
          }
        }
        else
        {
          MayDay::Error("EBAdvection: unrecognized bc");
        }
      }
    }
  }
}   
/*******/
void 
BCGVelAdvect::
copyExtrapolatedState(EBFluxData<Real, 1>& a_flux,
                      EBFluxData<Real, 1>& a_scalarLo, 
                      EBFluxData<Real, 1>& a_scalarHi, 
                      int idir, Side::LoHiSide sit, Bx valbx)
{
  Bx faceBx = valbx.faceBox(idir, sit);
  int isign = sign(sit);
  //unsigned long long int numflopspt = 0;
  if(idir == 0)
  {
    //using non-eb forall because box restriction in eb land is broken right now.   This will
    //work if there are no cut cells near the domain boundary
    auto& regflux =     a_flux.m_xflux->getRegData();
    auto& regsclo = a_scalarLo.m_xflux->getRegData();
    auto& regschi = a_scalarHi.m_xflux->getRegData();
    forallInPlaceBase(copyExtrap, faceBx, regflux, regsclo, regschi, isign);
  }
  else if(idir == 1)
  {
    //using non-eb forall because box restriction in eb land is broken right now.   This will
    //work if there are no cut cells near the domain boundary
    auto& regflux =     a_flux.m_yflux->getRegData();
    auto& regsclo = a_scalarLo.m_yflux->getRegData();
    auto& regschi = a_scalarHi.m_yflux->getRegData();
    forallInPlaceBase(copyExtrap, faceBx, regflux, regsclo, regschi, isign);
  }
#if DIM==3          
  else if(idir == 2)
  {
    //using non-eb forall because box restriction in eb land is broken right now.   This will
    //work if there are no cut cells near the domain boundary
    auto& regflux =     a_flux.m_zflux->getRegData();
    auto& regsclo = a_scalarLo.m_zflux->getRegData();
    auto& regschi = a_scalarHi.m_zflux->getRegData();
    forallInPlaceBase(copyExtrap, faceBx, regflux, regsclo, regschi, isign);
  }
#endif
  else
  {
    MayDay::Error("bogus idir");
  }
}
                      
                      
/*******/
void 
BCGVelAdvect::
hybridVecDivergence(EBLevelBoxData<CELL, DIM>& a_divuu,
                    EBLevelBoxData<CELL, DIM>& a_inputVel,
                    const Real               & a_dt,
                    Real a_tolerance, unsigned int a_maxIter)    
{
  CH_TIME("BCGAdvect::hybridDivergence");
  getMACVectorVelocity(a_inputVel, a_dt, a_tolerance, a_maxIter);    
  
  //lots of aliasing and smushing
  assembleDivergence(a_divuu, a_dt);
}
/*******/
void  
BCGVelAdvect::      
copyComp(EBFluxData<Real, 1>&  a_dst,
         EBFluxData<Real, 1>&  a_src,
         int a_vecDir)
{
  CH_TIME("BCGAdvect::copyComp");
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
                     const Real                  & a_dt,
                     Real a_tol, unsigned int a_maxIter)    
{
  CH_TIME("BCGAdvect::getMACVectorVelocity");
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

    m_source.exchange();
    DataIterator dit = m_grids.dataIterator();
    //int ideb = 0;
    for(int ibox = 0; ibox < dit.size(); ++ibox)
    {
      Bx   grid   =  ProtoCh::getProtoBox(m_grids[dit[ibox]]);
      Bx  grown   =  grid.grow(ProtoCh::getPoint(m_nghost));
      const EBGraph  & graph = (*m_graphs)[dit[ibox]];

      //this gets the low and high side states for the riemann problem
      bool useStack = true;
      auto& scalfab = velcomp[dit[ibox]];
      EBFluxData<Real, 1>  scalarHi(grown, graph, useStack);
      EBFluxData<Real, 1>  scalarLo(grown, graph, useStack);
      EBFluxData<Real, 1>  faceCentVelo( grown, graph, useStack);
      EBFluxData<Real, 1>&  upwindScal   =       facecomp[dit[ibox]];
      
      auto & veccell = a_inputVel[dit[ibox]];
      auto & sourfab =   m_source[dit[ibox]];
      bcgExtrapolateScalar(scalarLo, scalarHi, veccell, scalfab, sourfab,
                           grown, graph, dit[ibox], ibox, a_dt);

      //average velocities to face centers.
      getFaceCenteredVel( faceCentVelo, dit[ibox], ibox);

      //get face fluxes and interpolate them to centroids
      unsigned int curcomp  = vecDir;
      unsigned int doingvel = 1;
      getUpwindState(upwindScal, faceCentVelo, scalarLo, scalarHi, curcomp, doingvel);

      //enforce boundary conditions with an iron fist.
      applyVeloFluxBCs(upwindScal,  dit[ibox], 
                       scalarLo, scalarHi, curcomp);
      
      EBFluxData<Real, 1>& advvelfab = m_advectionVel[dit[ibox]];
      
      //now copy into the normal direction holder
      copyComp(advvelfab, upwindScal, vecDir);
    } //end loop over boxes
  }  // and loop over velocity directions

  // now we need to project the mac velocity
  pout() << "mac projecting advection velocity" << endl;
  bool printStuff = false;
  m_macproj->project(m_advectionVel, m_macGradient, a_tol, a_maxIter, printStuff);
  
  m_advectionVel.exchange();

  //subtract off pure gradient part of the velocity field in vector holder.
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
  CH_TIME("BCGAdvect::correctVectorVelocity");
  m_macGradient.exchange();
  DataIterator dit = m_grids.dataIterator();
  int ideb = 0;
  for(unsigned int vecDir = 0; vecDir < DIM; vecDir++)
    {
      EBLevelFluxData<1> facecomp;
      facecomp.define<DIM>(m_macVelocity, vecDir, m_graphs);

      for(int ibox = 0; ibox < dit.size(); ++ibox)
      {
        auto & advvelfab = m_advectionVel[dit[ibox]];
        auto & macGrad   =  m_macGradient[dit[ibox]];

        EBCrossFaceStencil<2, Real> stencils =
          m_brit->getCrossFaceStencil(StencilNames::TanVelCorrect, StencilNames::NoBC, m_domain, m_domain, ibox);

        auto & velcomp  = facecomp[dit[ibox]];
        for(unsigned int faceDir =0; faceDir < DIM; faceDir++)
        {
          unsigned int srcDir = vecDir;
          unsigned int dstDir = faceDir;
          if(faceDir == vecDir)
          {
            //when vecdir==facedir,
            //substitutNe the advection velocity, which already has its gradient removed
            //correct facedir==vecdir direction
            copyComp(velcomp, advvelfab, vecDir);
          }
          else
          {
            //begin debug
            //EBFluxData<Real,1> tangrad(velcomp.m_xflux->inputBox(), (*m_graphs)[dit[ibox]]);
            //stencils.apply(tangrad, macGrad, true      ,  1.0);
            //end debug

            bool setToZero=false;  Real scale = -1;  //subtracting off gradient
            stencils.apply(velcomp, macGrad, dstDir, srcDir, setToZero, scale);
          }
          ideb++;
        }
    }
    ideb++;
  }
  m_macVelocity.exchange();
}

/**********/
void 
BCGVelAdvect::
assembleDivergence(EBLevelBoxData<CELL, DIM>& a_divuu, 
                   const  Real              & a_dt)
{
  CH_TIME("BCGAdvect::assembleDivergence");
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
      Bx  grown   =  grid.grow(ProtoCh::getPoint(m_nghost));
      const EBGraph  & graph = (*m_graphs)[dit[ibox]];

      auto& faceCentVel  = m_advectionVel[dit[ibox]];
      auto& scalar       =  scalarVelComp[dit[ibox]];
      bool useStack = true;
      EBFluxData<Real, 1>  centroidFlux(grown, graph, useStack);
      EBFluxData<Real, 1>  faceCentFlux(grown, graph, useStack);

      assembleFlux(faceCentFlux, scalar, faceCentVel);

      EBFluxStencil<2, Real> stencils =   m_brit->getFluxStencil(s_centInterpLabel, s_nobcsLabel, m_domain, m_domain, ibox);
      //interpolate flux to centroids
      stencils.apply(centroidFlux, faceCentFlux, true, 1.0);  //true is to initialize to zero

      scalar.define(m_macVelocity[dit[ibox]], ivar);

      auto& kapdiv =  m_kappaDiv[dit[ibox]];
      getKapDivFFromCentroidFlux(kapdiv, centroidFlux, ibox);
    }
    m_kappaDiv.exchange();
    
    // begin debug
    //string sourcefile = string("kappaDiv.") + std::to_string(ivar) + string(".hdf5");
    //m_kappaDiv.writeToFileHDF5(sourcefile, 0.);
    //end debug
    
    //this computes the non-conservative divergence and puts it into m_nonConsDiv
    nonConsDiv();
    //this forms the hybrid divergence

    //does linear combination of divnc and div c to get hybrid and  compute delta M
    kappaDivPlusOneMinKapDivNC(hybridDiv);

    redistribute(hybridDiv);
  }
}
/*******/
#include "Chombo_NamespaceFooter.H"
