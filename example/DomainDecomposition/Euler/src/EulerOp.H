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
#include "LevelBoxData.H"
#include "BaseFab.H"
#include "Copier.H"

constexpr unsigned int NUMCOMPS=DIM+2;

#define NGHOST 4
typedef Proto::Box Bx;
using   Proto::Point;
using   Proto::BoxData;
using     std::cout;
using     std::endl;
using     std::shared_ptr;
using   Proto::Stencil;

class EulerOp
{

public:

  static Real s_dx;
  static Real s_gamma;

  /// Compute the Right-Hand side, conforming to the RK4 interface.
  static Real step(LevelBoxData<NUMCOMPS>& a_Rhs,
                   LevelBoxData<NUMCOMPS>& a_U);

  //for initial time step
  static Real maxWave(LevelBoxData<NUMCOMPS>& a_U);

  EulerOp(){};
  ~EulerOp(){};

private:
  // all the messy internals.
  static Real proto_step(BoxData<Real,NUMCOMPS>& a_Rhs,
                           const BoxData<Real,NUMCOMPS>& a_U,
                           const Bx& a_dbx0);

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

#endif
