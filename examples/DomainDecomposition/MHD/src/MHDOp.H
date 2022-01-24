#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _MHDOP_H_
#define _MHDOP_H_

#pragma once

#include "Proto.H"
#include "Chombo_LevelBoxData.H"
#include "Chombo_BaseFab.H"
#include "Chombo_Copier.H"

constexpr unsigned int NUMCOMPS=2*DIM+2;

#define NGHOST 8
typedef ::Proto::Box Bx;
using   ::Proto::Point;
using   ::Proto::BoxData;
using     std::cout;
using     std::endl;
using     std::shared_ptr;
using     std::max;
using     std::min;
using     std::abs;
using   ::Proto::Stencil;

class MHDOp
{

public:

  static Real s_dx;
  static Real s_gamma;

  /// Compute the Right-Hand side, conforming to the RK4 interface.
  static Real step(LevelBoxData<NUMCOMPS>& a_Rhs,
                   LevelBoxData<NUMCOMPS>& a_U);
						   
  /// Compute the Right-Hand side, conforming to the Euler step interface.
  static void step2(LevelBoxData<NUMCOMPS>& a_Rhs,
                   LevelBoxData<NUMCOMPS>& a_U);	

  /// Compute the Right-Hand side, conforming to the Euler step interface for artificial viscosity.
  static void step3(LevelBoxData<NUMCOMPS>& a_Rhs,
                   LevelBoxData<NUMCOMPS>& a_U);
				   
  /// Compute the testing variable Used for debugging.
  static void step_test(LevelBoxData<NUMCOMPS>& a_Rhs,
                   LevelBoxData<NUMCOMPS>& a_U);					   

  //for initial time step
  static Real maxWave(LevelBoxData<NUMCOMPS>& a_U);

  MHDOp(){};
  ~MHDOp(){};

private:
  // all the messy internals.
  static Real proto_step(BoxData<Real,NUMCOMPS>& a_Rhs,
                           const BoxData<Real,NUMCOMPS>& a_U,
                           const Bx& a_dbx0);
						   
  static void proto_step2(BoxData<Real,NUMCOMPS>& a_Rhs,
                           const BoxData<Real,NUMCOMPS>& a_U,
                           const Bx& a_dbx0);
						   
  static void proto_step3(BoxData<Real,NUMCOMPS>& a_Rhs,
                           const BoxData<Real,NUMCOMPS>& a_U,
                           const Bx& a_dbx0);					   

  static void proto_step_test(BoxData<Real,NUMCOMPS>& a_Rhs,
                           const BoxData<Real,NUMCOMPS>& a_U,
                           const Bx& a_dbx0);						   

  /// called internally from step
  static void define(const DisjointBoxLayout& a_grids,
                     const IntVect          & a_ghostVect) ;
  

  static Stencil<Real> s_laplacian;
  static Stencil<Real> s_deconvolve;
  static Stencil<Real> s_copy;
  static Stencil<Real> s_laplacian_f[DIM];
  static Stencil<Real> s_deconvolve_f[DIM];
  static Stencil<Real> s_convolve_f[DIM];
  static Stencil<Real> s_interp_H[DIM];
  static Stencil<Real> s_interp_L[DIM];
  static Stencil<Real> s_divergence[DIM];
  static Stencil<Real> s_ahead_shift[DIM];
  static Stencil<Real> s_behind_shift[DIM];
  static Stencil<Real> s_copy_f[DIM];
  
  static Copier          s_exchangeCopier;

};

#endif
