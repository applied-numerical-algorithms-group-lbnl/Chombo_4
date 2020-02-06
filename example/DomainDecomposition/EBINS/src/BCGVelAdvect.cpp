
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
#if 0  
  HERE 
  for(unsigned int idir = 0; idir < DIM; idir++)
  {
    EBLevelBoxData<CELL, 1> velcomp;
    velcomp.define(a_inputVel, idir, m_graphs);
    DataIterator dit = m_grids.dataIterator();
    int ideb = 0;
    for(int ibox = 0; ibox < dit.size(); ++ibox)
    {
      Bx   grid   =  ProtoCh::getProtoBox(m_grids[dit[ibox]]);
      Bx  grown   =  grid.grow(ProtoCh::getPoint(m_nghostSrc));
      Bx  grown   =  grid.grow(ProtoCh::getPoint(m_nghostSrc));

      const EBGraph  & graph = (*m_graphs)[dit[ibox]];
      //get face fluxes and interpolate them to centroids
      EBFluxData<Real, 1>  faceCentVel( grown, graph);
      //average velocities to face centers.
      getFaceCenteredVel( faceCentVel, dit[ibox], ibox);

    }
  }
#endif
}
/*******/
void  
BCGVelAdvect::
getMACVectorVelocity(EBLevelBoxData<CELL, DIM>   & a_inputVel,
                     const Real                  & a_dt)
{
}
/*******/
void 
BCGVelAdvect::
assembleDivergence()
{
}
/*******/
#include "Chombo_NamespaceFooter.H"

