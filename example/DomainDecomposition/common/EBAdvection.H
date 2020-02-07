#ifndef _EBAdvection_H_
#define _EBAdvection_H_
#include <cmath>
#include <memory>
#include "Proto.H"
#include "Chombo_EBDictionary.H"
#include "Chombo_EBChombo.H"
#include "Chombo_EBEncyclopedia.H"
#include "Chombo_EBLevelBoxData.H"
#include "Chombo_NamespaceHeader.H"

///class to advect scalars in an eb context (via Trebotich et al.)
class EBAdvection
{
public:
  /// 
  EBAdvection(shared_ptr<EBEncyclopedia<2, Real> >   & a_brit,
              shared_ptr<GeometryService<2> >        & a_geoserv,
              shared_ptr<EBLevelBoxData<CELL, DIM> > & a_veloCell,
              const DisjointBoxLayout                & a_grids,
              const Box                              & a_domain,
              const Real                             & a_dx,
              const IntVect                          & a_nghostsrc, 
              const IntVect                          & a_nghostdst);

  virtual ~EBAdvection()
  {
  }

  /// advance one time step (via Trebotich et al.) in  an eb context
  void 
  advance(EBLevelBoxData<CELL, 1>   & a_scal,
          const  Real               & a_dt);

  const EBLevelBoxData<CELL, 1> &
  getKappa() const
  {
    return m_kappa;
  }

  virtual void  kappaConsDiv(EBLevelBoxData<CELL, 1>   & a_scal,
                             const Real                & a_dt);
  ///fills m_hybridDiv with hybrid divergence
  virtual void 
  hybridDivergence(EBLevelBoxData<CELL, 1>   & a_scal,
                   const  Real               & a_dt);

  EBLevelBoxData<CELL, 1>                m_hybridDiv;
  EBLevelBoxData<CELL, 1>                m_kappa;

  //additional terms added to extrapolation
  EBLevelBoxData<CELL, 1>                m_source;
protected:

  //rule britannica
  shared_ptr<EBEncyclopedia<2, Real> >   m_brit;
  DisjointBoxLayout                      m_grids;
  Box                                    m_domain;
  IntVect                                m_nghostSrc;
  IntVect                                m_nghostDst;
  shared_ptr<EBLevelBoxData<CELL, DIM> > m_veloCell;
  Real m_dx;
  EBLevelBoxData<CELL, 1>                m_deltaM;
  EBLevelBoxData<CELL, 1>                m_nonConsDiv;
  EBLevelBoxData<CELL, 1>                m_kappaDiv;
  shared_ptr<LevelData<EBGraph>  >       m_graphs;
  Copier                                 m_exchangeCopier;

  static const string s_ncdivLabel;      //for non-conservative divergence
  static const string s_nobcsLabel;      //none of the operators here have eb boundary conditions
  static const string s_redistLabel;     //for redistribution
  static const string s_centInterpLabel; //for interpolation to face centroids
  static const string s_slopeLowLabel;   //for cell centered slopes
  static const string s_slopeHighLabel;  //for cell centered slopes
  static const string s_aveCToFLabel;    //get the velocity from cell centers to faces
  static const string s_diriLabel;       //label for dirichlet bcs
  static const string s_neumLabel;       //label for neumann bcs
  static const string s_divergeLabel;    //label incrementing  the divergence with face fluxes
  static const string s_CtoFLowLabel;    //for getting stuff from low  side cells to faces
  static const string s_CtoFHighLabel;   //for getting stuff from high side cells to faces

  template<CENTERING cent>
  inline void
  getFaceVelComp(EBBoxData<cent, Real, 1>                       & a_fcvel,
                 shared_ptr< EBDictionary<2, Real, CELL, cent> >& a_dict,
                 EBBoxData<CELL, Real, DIM>                     & a_velcc,
                 unsigned int                                     a_idir,
                 int                                              a_ibox)
  {
    EBBoxData<CELL, Real, 1> cellcomp;
    cellcomp.define<DIM>(a_velcc, a_idir);
    // stencils to average velocities from cells to faces
    const auto& facesten = a_dict->getEBStencil(s_aveCToFLabel, s_nobcsLabel, m_domain, m_domain, a_ibox);
    bool initToZero = true;
    facesten->apply(a_fcvel, cellcomp, initToZero, 1.0);               
  }

  void getFaceCenteredVel(EBFluxData<Real, 1>& a_fcvel,
                          const DataIndex    & a_dit,
                          const int          & a_ibox);

  void
  getFaceCenteredFlux(EBFluxData<Real, 1>            & a_fcflux,
                      const EBFluxData<Real, 1>      & a_fcvel,
                      const EBBoxData<CELL, Real, 1> & a_scal,
                      const DataIndex                & a_dit,
                      const int                      & a_ibox,
                      const Real                     & a_dt);

  void  defineData(shared_ptr<GeometryService<2> >        & a_geoserv);
  void  nonConsDiv();
  void redistribute();
  void fillKappa(shared_ptr<GeometryService<2> >        & a_geoserv);
  void registerStencils();

  void operator=(const EBAdvection& a_opin);
  void
  bcgExtrapolateScalar(EBFluxData<Real, 1>            & a_scalHi,
                       EBFluxData<Real, 1>            & a_scalLo,
                       EBBoxData<CELL, Real, DIM>     & a_veccell,
                       const EBBoxData<CELL, Real, 1> & a_scal,
                       const EBBoxData<CELL, Real, 1> & a_source,
                       const Bx                       & a_grown,
                       const EBGraph                  & a_graph,
                       const DataIndex                & a_dit,
                       const int                      & a_ibox,
                       const Real                     & a_dt);
  EBAdvection();

};
#include "Chombo_NamespaceFooter.H"

#endif
