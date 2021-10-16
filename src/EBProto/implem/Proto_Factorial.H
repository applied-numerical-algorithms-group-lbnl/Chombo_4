
#ifndef __Proto_FACTORIAL_H_
#define __Proto_FACTORIAL_H_

#include "Proto_RealVect.H"
#include "Proto_IndexTM.H"

namespace Proto
{
/**
 *  These helper functions calculate factorials, binomial coeff, etc.
 *  for multinomials.
 */

/// Calculates factorial for an integer
  inline double factorial(const int n)
  {
    PROTO_ASSERT(n >= 0, "negative n");
    double nfact = 1.0;
    for (int i = 2; i <= n; nfact*= i, i++);
    return nfact;
  }
/// computes x^p
  inline double POW(const double& a_x,  const int& a_p)
  {
    double retval = 1;
    for(int iexp = 0; iexp < std::abs(a_p); iexp++)
    {
      if(a_p >= 0)
      {
        retval *= a_x;
      }
      else
      {
        retval /= a_x;
      }
    }
    return retval;
  }

/// Calculates the binomial coefficient, "n choose k" 
  inline double nCk(const int n, const int k)
  {
    PROTO_ASSERT((n >= k)&&(k >= 0)&&(n >= 0), "out of range n");
    double nck = 1.0;
    for (int i = n-k+1; i <= n; nck*= i, i++);
    return nck / factorial(k);
  }

/// Calculates factorials for a multinomial
  template <int Dim>
  inline double pfactorial(const IndexTM<int, Dim> p)
  {
    double pfact = 1;
    for(int idir = 0; idir < Dim; idir++)
    {
      PROTO_ASSERT(p[idir] >= 0, "p must be positive");
      pfact *=  factorial(p[idir]); 
    }
    return pfact;
  }

/// Calculates the multinomial coefficient, "p choose k" 
  template <int Dim>
  inline double pCk(const IndexTM<int, Dim>& p, const IndexTM<int, Dim>& k)
  {
    //PROTO_ASSERT((p >= IntVect::Zero)&&(k >= IntVect::Zero)&&(p >= k));

    double pfact = 1;
    for(int idir = 0; idir < Dim; idir++)
    {
      pfact *= nCk(p[idir],k[idir]); 
    }
    return pfact;
  }


/// calculate x^p
  template <int Dim>
  inline double power(const IndexTM<double, Dim>& a_x, const IndexTM<int, Dim>& a_p)
  {
    double retval = 1;
    for(int idir = 0; idir < Dim; idir++)
    {
      if((a_p[idir] >= 0))
      {
        for(int ipow = 0; ipow < a_p[idir]; ipow++)
        {
          retval *= a_x[idir];
        }
      }
      else
      {
        for(int ipow = 0; ipow < a_p[idir]; ipow++)
        {
          retval /= a_x[idir];
        }
      }
    }
    return retval;
  }


}
#endif
