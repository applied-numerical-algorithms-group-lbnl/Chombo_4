#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _EULEROP_H_
#define _EULEROP_H_

#pragma once

#include "Proto.H"
#include "Chombo_LevelBoxData.H"
#include "Chombo_BaseFab.H"
#include "Chombo_Copier.H"
#include "Chombo_NamespaceHeader.H"
constexpr unsigned int NUMCOMPS=DIM+2;

#define NGHOST 4
typedef ::Proto::Box Bx;
using   ::Proto::Point;
using   ::Proto::BoxData;
using     std::cout;
using     std::endl;
using     std::shared_ptr;
using   ::Proto::Stencil;
using   Op = ::Proto::Operation;
using   ::Proto::Reduction;

//global reduction of each processor's max wave
Real gatherMaxWave(Real maxwaveproc);

class EulerOp
{

public:

  static Real s_dx;
  static Real s_gamma;

  /// Compute the Right-Hand side, conforming to the RK4 interface.
  static Real step(LevelBoxData<NUMCOMPS>& a_Rhs,
                   LevelBoxData<NUMCOMPS>& a_U,
                   Reduction<Real,Op::Abs>& a_Rxn);

  //for initial time step
  static Real maxWave(LevelBoxData<NUMCOMPS>& a_U, 
                      Reduction<Real,Op::Abs>& a_Rxn); 

  EulerOp(){};
  ~EulerOp(){};

private:
  // all the messy internals.
  static Real proto_step(BoxData<Real,NUMCOMPS>& a_Rhs,
                         const BoxData<Real,NUMCOMPS>& a_U,
                         const Bx& a_dbx0,
                         Reduction<Real,Op::Abs>& a_Rxn);

  /// called internally from step
  static void define(const DisjointBoxLayout& a_grids,
                     const IntVect          & a_ghostVect) ;
  

  static Stencil<Real> s_laplacian;
  static Stencil<Real> s_deconvolve;
  static Stencil<Real> s_laplacian_f[DIM];
  static Stencil<Real> s_deconvolve_f[DIM];
  static Stencil<Real> s_interp_H[DIM];
  static Stencil<Real> s_interp_L[DIM];
  static Stencil<Real> s_divergence[DIM];
  static Copier          s_exchangeCopier;

};
#include "Chombo_NamespaceFooter.H"
#endif
