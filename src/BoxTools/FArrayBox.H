#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _FARRAYBOX_H_
#define _FARRAYBOX_H_

#ifndef WRAPPER
#include <iostream>

// #include <Pointers.H>
#include "Box.H"
#include "BaseFab.H"
#endif

#include "REAL.H"
#include "SPACE.H"
#include "NamespaceHeader.H"

typedef BaseFab<Real> FArrayBox;

#include "NamespaceFooter.H"
#endif
