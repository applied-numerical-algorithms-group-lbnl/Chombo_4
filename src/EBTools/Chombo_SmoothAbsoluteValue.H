#ifndef _CH_SMOOTHABSOLUTEVALUE_H_
#define _CH_SMOOTHABSOLUTEVALUE_H_

#include "Chombo_MayDay.H"
#include "Chombo_RealVect.H"
#include "Chombo_IntVect.H"
#include "Proto_BaseIF.H"

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
  SmoothAbsoluteValue(const BaseIF*   a_f,
                      const BaseIF*   a_g,
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
  virtual Real smoothAbsFMinusG(const  IntVect& a_deriv,
                                const RealVect& a_point) const;

  ///
  /**
     returns -1 if w < -delta, 1 if w > delta, 0 otherwise
     reduces to regular |f-g| unless case == 0
  */
  void getWCase(int            & a_case,
                Real           & a_wval,
                const RealVect & a_point)const;

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


  virtual Real valueAem(const RealVect& a_point) const;

  ///
  virtual Real firstDerivAem(const  IntVect& a_deriv,
                             const RealVect& a_point) const;

  ///
  virtual Real secondDerivAem(const  IntVect& a_deriv,
                              const RealVect& a_point) const;

  ///
  virtual Real thirdDerivAem(const  IntVect& a_deriv,
                             const RealVect& a_point) const;

  ///
  virtual Real fourthDerivAem(const  IntVect& a_deriv,
                              const RealVect& a_point) const;

  ///just checks nan
  bool isBogus(const Real& a_number) const;


                              
  //the two implicit functions are owned by others
  const BaseIF*   m_f;
  const BaseIF*   m_g;
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
