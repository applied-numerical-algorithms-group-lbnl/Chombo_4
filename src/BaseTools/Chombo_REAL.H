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

#include <cmath>
#include <limits>
#include <algorithm>

#  undef  CH_USE_FLOAT
#  undef  CH_USE_DOUBLE
#  define CH_USE_DOUBLE 1
#    define MPI_CH_REAL         MPI_DOUBLE
#    define H5T_NATIVE_REAL     H5T_NATIVE_DOUBLE
#    define BASEFAB_REAL_SETVAL 1.23456789e+300
#    define CH_BADVAL           (Real)1.0e99

typedef double Real;

// define common floating point limits
const static Real REAL_MIN = std::numeric_limits<Real>::min();
const static Real REAL_MAX = std::numeric_limits<Real>::max();
const static Real REAL_LOW = std::numeric_limits<Real>::lowest();
const static Real REAL_EPS = std::numeric_limits<Real>::epsilon();
const static Real REAL_INF = std::numeric_limits<Real>::infinity();
// A reasonable floating point tolerance for general applications
// about 3/4 of the available digits (1e-12 for double, 1e-6 for float)
const static Real REAL_TOL = std::pow(REAL_EPS, 3./4.);

// Does a == b within floating point tolerance?
inline bool
approxEqual(Real a, Real b, Real tol=REAL_TOL)
{
    Real absA = std::abs(a);
    Real absB = std::abs(b);
    Real diff = std::abs(a - b);

    if (a == b) // shortcut for exact equality
      {
        return true;
      }
    else if (a == 0 || b == 0 || diff < REAL_LOW)
      { // handle cases where one or both numbers are close to zero or
        // the difference between the two numbers is very small
        return diff < (tol * REAL_LOW);
      }
    else // use relative error
      {
        return diff < (tol * std::min((absA + absB), REAL_MAX));
      }
}

// Does c = a + b have a meaningful contribution from both a and b?
// Default tolerance is the minimum visible value
inline bool
isSignificant(Real a, Real b, Real tol=REAL_EPS) {
    Real abs_a = std::abs(a);
    Real abs_b = std::abs(b);

    Real minVal = std::min(abs_a, abs_b);
    Real maxVal = std::max(abs_a, abs_b);

    return minVal > (tol*maxVal);
}


#endif
