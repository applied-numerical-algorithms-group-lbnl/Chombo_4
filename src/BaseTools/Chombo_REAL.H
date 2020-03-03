#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _Chombo_REAL_H_
#define _Chombo_REAL_H_

#  undef  CH_USE_FLOAT
#  undef  CH_USE_DOUBLE
#  define CH_USE_DOUBLE 1
#    define MPI_CH_REAL         MPI_DOUBLE
#    define H5T_NATIVE_REAL     H5T_NATIVE_DOUBLE
#    define BASEFAB_REAL_SETVAL 1.23456789e+300
#    define CH_BADVAL           (Real)1.0e99

typedef double Real ;

#endif
