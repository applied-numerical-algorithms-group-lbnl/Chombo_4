#ifndef _CH_SMOOTHABSOLUTEVALUE_H_
#define _CH_SMOOTHABSOLUTEVALUE_H_

#include "Chombo_MayDay.H"
#include "Chombo_RealVect.H"
#include "Chombo_IntVect.H"
#include "Proto_BaseIF.H"
#include "Proto_IndexTM.H"
#include "Chombo_NamespaceHeader.H"

///
/**
   Functions to take a smooth absolute value of the difference between two functions
   A_e(a,b) = |a-b|

   max(f, g) = (a + b + |a - b|)/2
   We need to make abs(a-b) smooth.

   We take the convolution of the absolute value with cos^4(pi w/(2 delta))
 */
class SmoothAbsoluteValue
{
public:

  
  ///
  SmoothAbsoluteValue(const Proto::BaseIF*   a_f,
                      const Proto::BaseIF*   a_g,
                      const Real  &   a_delta)
  {
    m_f  = a_f;
    m_g  = a_g;
    m_d  = a_delta;
    m_pi = 4.*atan(1.0);
  }

  ///
  virtual ~SmoothAbsoluteValue()
  {; }

  ///
  virtual Real smoothAbsFMinusG(const  Proto::IndexTM<int,   DIM>& a_deriv,
                                const Proto::IndexTM<double, DIM>& a_point) const;

  ///
  /**
     returns -1 if w < -delta, 1 if w > delta, 0 otherwise
     reduces to regular |f-g| unless case == 0
  */
  void getWCase(int            & a_case,
                Real           & a_wval,
                const Proto::IndexTM<double, DIM> & a_point)const;

  // gives x^p (thought to be cheaper than calling the system function because of the int)
  static Real POW(Real a_x, int a_p)
  {
    Real retval = 1;
    for(unsigned int ipow = 0; ipow < a_p; ipow++)
    {
      retval *= a_x;
    }
    return retval;
  }
protected:

  ///
  /**
     Here is the logic of this stuff.   We have three cases that reduce 
     to two cases.    w = f(x) - g(x)
     case  1: (w > delta):  ---- whole integral is above zero
        answer = abs(w)
     case -1: (w < - delta): ---- whole integral is below zero
        answer = abs(w)
     case  0: (-delta <= w <= delta)  --- have to split integral into above and below
        answer = functionAem();
   */


  virtual Real valueAem(const Proto::IndexTM<double, DIM>& a_point) const;

  ///
  virtual Real firstDerivAem(const  Proto::IndexTM<int,   DIM>& a_deriv,
                             const Proto::IndexTM<double, DIM>& a_point) const;

  ///
  virtual Real secondDerivAem(const  Proto::IndexTM<int,   DIM>& a_deriv,
                              const Proto::IndexTM<double, DIM>& a_point) const;

  ///
  virtual Real thirdDerivAem(const  Proto::IndexTM<int,   DIM>& a_deriv,
                             const Proto::IndexTM<double, DIM>& a_point) const;

  ///
  virtual Real fourthDerivAem(const  Proto::IndexTM<int,   DIM>& a_deriv,
                              const Proto::IndexTM<double, DIM>& a_point) const;

  ///just checks nan
  bool isBogus(const Real& a_number) const;


                              
  //the two implicit functions are owned by others
  const Proto::BaseIF*   m_f;
  const Proto::BaseIF*   m_g;
  //delta = the smoothing length scale
  Real m_d;

  //pi, you know, pi = 4atan(1)
  Real m_pi;

private:

  SmoothAbsoluteValue();
  SmoothAbsoluteValue(const SmoothAbsoluteValue& a_inputIF);
  void operator=(     const SmoothAbsoluteValue& a_inputIF);

};

#include "Chombo_NamespaceFooter.H"
#endif
