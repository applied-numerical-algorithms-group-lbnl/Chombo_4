#ifndef _ch4_SMOOTHUNION_H_
#define _ch4_SMOOTHUNION_H_

#include "Chombo_MayDay.H"
#include "Chombo_RealVect.H"
#include "Proto_BaseIF.H"
#include "Proto_IndexTM.H"

namespace Chombo4
{
///
/**
   This implicit function is the maximum of two (of more) implicit
   functions with which it is constructed.  This results in an implicit
   function whose interior is the intersection of the interiors of the given
   implicit functions.
   For the maximum, it uses the smooth absolute value Ae(a,b) \approx |a-b|
   min(a,b) = 0.5*(a+b - | a-b |)
*/
  class SmoothUnion: public Proto::BaseIF
  {
  public:

    ///
    /**
       Constructor specifying any number of implicit functions as inputs.
    */
    SmoothUnion(const std::vector<BaseIF*>& a_impFuncs,
                const Real           & a_delta);

    /// Destructor
    virtual ~SmoothUnion();

    ///
    virtual Real value(const  Proto::IndexTM<int,   DIM>& a_deriv,
                       const Proto::IndexTM<double, DIM>& a_point) const;


    virtual BaseIF* newImplicitFunction() const;

    //min(a,b) = 0.5*(a+b - | a-b |)
    virtual Real smoothMin(const  Proto::IndexTM<int,   DIM> & a_deriv,
                           const Proto::IndexTM<double, DIM> & a_point,
                           const  int     & a_closestIF,
                           const  int     & a_nextClosestIF) const;

    void findClosest(int            & a_closestIF, 
                     int            & a_nextClosestIF,
                     int            & a_numWithinDelta,
                     const Proto::IndexTM<double, DIM> & a_point) const;
  protected:
    int                  m_numFuncs; // number of implicit functions
    std::vector<BaseIF*> m_impFuncs; // implicit functions
    Real                 m_delta   ; //smoothing length scale

  private:
    SmoothUnion();
    SmoothUnion(const SmoothUnion& a_inputIF);
    void operator=(const SmoothUnion& a_inputIF);

  };

}
#endif
