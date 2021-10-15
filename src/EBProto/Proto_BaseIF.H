#ifndef __proto_BASEIF_H_
#define __proto_BASEIF_H_

#include "Proto_RealVect.H"
#include "Proto_IndexTM.H"
#include "Proto_Point.H"
#include "Proto_Box.H"


namespace Proto
{
  /// Class that specifies and implicit function
  /**
     This is the base class for an implicit function specification of geometry.
     All that is needed is a constructor/destructor, a method to give the value
     of the function at any point in space (in 2D or 3D), and a factory method
     (these will probably all be very similar).
  */
  class BaseIF
  {

  public:
    
    /// Default constructor
    BaseIF()
    {
    }

    /// Default destructor
    virtual ~BaseIF()
    {
    }

    ///
    /**
       Return the value of the function at a_point.  When delineating a domain,
       the level set value=0 represents the boundary and value<0 is inside the
       fluid.
    */
    virtual double value(const IndexTM<int   ,DIM> & a_partialDerivative,
                         const IndexTM<double,DIM>& a_point) const = 0;

    //.overwrite if you have a fast way of knowing if a whole region is regular
    virtual bool entireBoxRegular(const Proto::Box     & a_box,
                                  const RealVect& a_origin,
                                  const double  & a_dx)
    {
      return false;
    }

    ///overwrite if you have a fast way of knowing if a whole region is regular
    virtual bool entireBoxCovered(const Proto::Box&       a_box,
                                  const RealVect& a_origin,
                                  const double&   a_dx)
    {
      return false;
    }

    /// factory method
    /**
       Return a newly allocated derived class.  The responsibility
       for deleting the memory is left to the calling function.
    */
    virtual BaseIF* newImplicitFunction() const = 0;

  };

}
#endif
