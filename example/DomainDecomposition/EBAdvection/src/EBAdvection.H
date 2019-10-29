#ifndef _EBAdvection_H_
#define _EBAdvection_H_
#include <cmath>
#include <memory>
#include "Proto.H"
#include "EBDictionary.H"
#include "EBChombo.H"
#include "EBEncyclopedia.H"
#include "EBLevelBoxData.H"
#include "NamespaceHeader.H"

///class to advect scalars in an eb context (via Trebotich et al.)
class EBAdvection
{
public:


  static const string s_ncdivLabel;      //for non-conservative divergence
  static const string s_nobcsLabel;      //none of the operators here have eb boundary conditions
  static const string s_redistLabel;     //for redistribution
  static const string s_centInterpLabel; //for interpolation to face centroids
  static const string s_slopeLowLabel;   //for cell centered slopes
  static const string s_slopeHighLabel;  //for cell centered slopes
  static const string s_aveCToFLabel;    //get the velocity from cell centers to faces
  static const string s_diriLabel;       //label for dirichlet bcs
  static const string s_divergeLabel;    //label incrementing  the divergence with face fluxes
  
  /// 
  EBAdvection(shared_ptr<EBEncyclopedia<2, Real> >   & a_brit,
              shared_ptr<GeometryService<2> >        & a_geoserv,
              shared_ptr<EBLevelBoxData<CELL, DIM> > & a_veloCell,
              const DisjointBoxLayout                & a_grids,
              const Box                              & a_domain,
              const Real                             & a_dx,
              const IntVect                          & a_nghostsrc, 
              const IntVect                          & a_nghostdst);

  /// advance one time step (via Trebotich et al.) in  an eb context
  void 
  advance(EBLevelBoxData<CELL, 1>   & a_scal,
          const  Real               & a_dt);

protected:

  //rule britannica
  shared_ptr<EBEncyclopedia<2, Real> >   m_brit;
  DisjointBoxLayout                      m_grids;
  Box                                    m_domain;
  IntVect                                m_nghostSrc;
  IntVect                                m_nghostDst;
  shared_ptr<EBLevelBoxData<CELL, DIM> > m_veloCell;

  EBLevelBoxData<CELL, 1>                m_kappa;
  EBLevelBoxData<CELL, 1>                m_deltaM;
  EBLevelBoxData<CELL, 1>                m_nonConsDiv;
  EBLevelBoxData<CELL, 1>                m_hybridDiv;
  shared_ptr<LevelData<EBGraph>  >       m_graphs;
  Copier                                 m_exchangeCopier;

private:
  void getFaceCenteredVel(EBFluxData<Real, 1>& a_fcvel,
                          const DataIndex    & a_dit,
                          const int          & a_ibox);
  
  void  defineData(shared_ptr<GeometryService<2> >        & a_geoserv);
  void  kappaConsDiv(EBLevelBoxData<CELL, 1>   & a_scal);
  void  nonConsDiv();
  void redistribute();
  void fillKappa(shared_ptr<GeometryService<2> >        & a_geoserv);
  void registerStencils();

  void operator=(const EBAdvection& a_opin);
  EBAdvection();

};
#include "NamespaceFooter.H"

#endif
