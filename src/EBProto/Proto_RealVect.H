
#ifndef _Proto_REALVECT_H_
#define _Proto_REALVECT_H_

#include <cstddef>
#include <cstdlib>
#include <cstring>
#include <vector>

#include "Proto_Point.H"
#include "implem/Proto_IndexTM.H"
using std::vector;
namespace Proto
{
//template<typename T, int n>
//class IndexTM;

/// A double vector in DIM-dimensional space
/**
   The class RealVect is an implementation of a double vector in a
   DIM-dimensional space.
   RealVect values are accessed using the operator[] function, as for a normal
   C++ array.  In addition, the basic arithmetic operators have been overloaded
   to implement scaling and translation operations.
*/

  class RealVect
  {
  public:

    /**
       \name Constructors and Accessors
    */
    /*@{*/

    ///
    /**
       Construct a RealVect whose components are zero.
    */
    RealVect ()
    {
      for(int idir = 0; idir < DIM; idir++)
      {
        vect[idir] = 0;
      }
    }

    ///
    /**
       The copy constructor.
    */
    RealVect (const RealVect& rhs)
    {
      for(int idir = 0; idir < DIM; idir++)
      {
        vect[idir] = rhs.vect[idir];
      }
    }


    ///
    /**
       Construct a RealVect from an IndexTM<double,DIM>.
    */
    RealVect(const IndexTM<double,DIM>& a_tm)
    {
      for(int idir = 0; idir < DIM; idir++)
      {
        vect[idir] = a_tm[idir];
      }
    }

    ///
    /**
       Construct a RealVect from an Point by coercing each component
       from <tt>int</tt> to double.
    */
    RealVect (const Point & iv)
    {
      for (int d=0 ; d<DIM ; ++d)
      {
        vect[d] = (double)iv[d];
      }
    }


    RealVect (const double & iv)
    {
      for (int d=0 ; d<DIM ; ++d)
      {
        vect[d] = iv;
      }
    }

    ///
    /**
       The assignment operator.
    */
    RealVect& operator= (const RealVect& rhs)
    {
      for (int d=0 ; d<DIM ; ++d)
      {
        vect[d] = rhs.vect[d];
      }
      return *this;
    }

    ///
    /**
       Returns a modifiable lvalue reference to the <i>i</i>'th coordinate of the
       RealVect.
    */
    inline
    double& operator[] (int i)
    {
      return vect[i];
    }

    ///
    /**
       Returns the <i>i</i>'th coordinate of the RealVect.
    */
    inline
    const double& operator[] (int i) const
    {
      return vect[i];
    }

    /**
       \name Arithmetic Operators
    */
    /*@{*/

    ///
    /**
       Modifies this RealVect by addition of a scalar to each component.
    */
    RealVect& operator+= (double s)
    {
      for(int idir = 0; idir < DIM; idir++)
      {
        vect[idir] += s;
      }
      return *this;
    }

    ///
    /**
       Returns a RealVect that is this RealVect with a scalar s added
       to each component.
    */
    RealVect operator+ (double s) const
    {
      RealVect retval = *this;
      for(int idir = 0; idir < DIM; idir++)
      {
        retval.vect[idir] += s;
      }
      return retval;
    }

    ///
    /**
       Modifies this RealVect by component-wise addition by argument.
    */
    RealVect& operator+= (const RealVect& p)
    {
      for(int idir = 0; idir < DIM; idir++)
      {
        vect[idir] += p.vect[idir];
      }
      return *this;
    }

    ///
    /**
       Modifies this RealVect by subtraction of a scalar from each
       component.
    */
    RealVect& operator-= (double s)
    {
      for(int idir = 0; idir < DIM; idir++)
      {
        vect[idir] -= s;
      }
      return *this;
    }


    ///
    /**
       Modifies this RealVect by component-wise subtraction by argument.
    */
    RealVect& operator-= (const RealVect& p)
    {
      for(int idir = 0; idir < DIM; idir++)
      {
        vect[idir] -= p.vect[idir];
      }
      return *this;
    }

    ///
    /**
       Returns a RealVect that is this RealVect with a scalar s subtracted
       from each component.
    */
    RealVect operator- (double s) const
    {
      RealVect retval = *this;
      for(int idir = 0; idir < DIM; idir++)
      {
        retval.vect[idir] -= s;
      }
      return retval;
    }

    ///
    /**
       Modifies this RealVect by multiplying each component by a scalar.
    */
    RealVect& operator*= (double s)
    {
      for(int idir = 0; idir < DIM; idir++)
      {
        vect[idir] *= s;
      }
      return *this;
    }

    RealVect operator* (double s)
    {
      for(int idir = 0; idir < DIM; idir++)
      {
        vect[idir] *= s;
      }
      return *this;
    }


    ///
    /**

     */
    double dotProduct(const RealVect& a_rhs) const
    {
      double retval = 0;
      for(int idir = 0; idir < DIM; idir++)
      {
        retval += vect[idir]*a_rhs.vect[idir];
      }
      return retval;
    }



//XXX  ///
//XXX  /**
//XXX     Returns component-wise product of this RealVect and argument.
//XXX  */
//XXX  RealVect operator* (const RealVect& p) const;

    ///
    /**
       Returns a RealVect that is this RealVect with each component
       multiplied by a scalar.
    */
    RealVect operator* (double s) const
    {
      RealVect retval = *this;
      for(int idir = 0; idir < DIM; idir++)
      {
        retval.vect[idir] *= s;
      }
      return retval;
    }

    ///
    /**
       Modifies this RealVect by dividing each component by a scalar.
    */
    RealVect& operator/= (double s)
    {
      for(int idir = 0; idir < DIM; idir++)
      {
        vect[idir] /= s;
      }
      return *this;
    }


    ///
    /**
       Returns a RealVect that is this RealVect with each component
       divided by a scalar.
    */
    RealVect operator/ (double s) const
    {
      RealVect retval = *this;
      for(int idir = 0; idir < DIM; idir++)
      {
        retval.vect[idir] /= s;
      }
      return retval;
    }


    ///
    /**
       Unary plus -- for completeness.
    */
    RealVect operator+ () const
    {
      return RealVect(*this);
  }

    ///
    /**
       Unary minus -- negates all components of this RealVect.
    */
    RealVect operator- () const
    {
      RealVect retval = *this;
      retval *= -1.0;
      return retval;
    }

    ///
    /**
       Sum of all components of this RealVect.
    */
    double sum () const
    {
      double retval = 0;
      for(int idir = 0; idir < DIM; idir++)
      {
        retval += vect[idir];
      }

      return retval;
    }

    ///
    /**
       sqrt(sum squares)
    */
    double vectorLength() const
    {
      double len = this->radSquared();
      len = sqrt(len);

      return len;
    }

    ///
    /**
       sum squares--no square root
    */
    double radSquared() const
    {
      double len = 0;
      for (int idir = 0; idir < DIM; idir++)
      {
        len = len + vect[idir]*vect[idir];
      }

      return len;
    }


    ///
    /**
       Component with the minimum value of this RealVect (returns 0 if they are all the same).
       a_doAbs : if true then take the absolute value before comparing
    */
    int minDir(const bool& a_doAbs) const
    {
      int retval = 0;
      if(a_doAbs)
      {
        double compval = std::abs(vect[0]);
        for(int idir = 1; idir < DIM; idir++)
        {
          if(std::abs(vect[idir]) < compval)
          {
            retval = idir;
            compval = std::abs(vect[idir]);
          }
        }
      }
      else
      {
        double compval = vect[0];
        for(int idir = 1; idir < DIM; idir++)
        {
          if(vect[idir] < compval)
          {
            retval = idir;
            compval = vect[idir];
          }
        }
      }
      return retval;
    }

    ///
    /**
       Component with the maximum value of this RealVect (returns 0 if they are all the same).
       a_doAbs : if true then take the absolute value before comparing
    */
    int maxDir(const bool& a_doAbs) const
    {
      int retval = 0;
      if(a_doAbs)
      {
        double compval = std::abs(vect[0]);
        for(int idir = 1; idir < DIM; idir++)
        {
          if(std::abs(vect[idir]) > compval)
          {
            retval = idir;
            compval = std::abs(vect[idir]);
          }
        }
      }
      else
      {
        double compval = vect[0];
        for(int idir = 1; idir < DIM; idir++)
        {
          if(vect[idir] > compval)
          {
            retval = idir;
            compval = vect[idir];
          }
        }
      }
      return retval;
    }


    /*@}*/


    /**
       \name Constants
    */
    /*@{*/


    ///
    /**
       This is a RealVect all of whose components are equal to zero.
    */
    inline static RealVect Zero()
    {
      static RealVect rvzero(0.0);
      return rvzero;
    }

    ///
    /**
       This is a RealVect all of whose components are equal to one.
    */
    inline static RealVect Unit()
    {
      static RealVect rvone(1.0);
      return rvone;
    }

    /*@}*/


    /**
       The individual components of this RealVect.
    */
    double vect[DIM];

  };

  inline
  std::ostream& operator<<(std::ostream   & a_os,
                           const RealVect & a_rv)
  {
    a_os << "(";
    for(int idir = 0; idir < DIM; idir++)
    {
      a_os << a_rv[idir];
      if(idir < DIM-1)
      {
        a_os << ",";
      }
    }
    a_os << ")";
    return a_os;
  }


  inline RealVect BASISREALV(int idir)
  {
    RealVect retval = RealVect::Zero();
    retval[idir] = 1;
    return retval;
  }

  inline void convertToITM(IndexTM<double, DIM>& a_itmpoint, 
                           const RealVect      & a_rvpoint)
  {
    for(int idir = 0; idir < DIM; idir++)
    {
      a_itmpoint[idir] = a_rvpoint[idir];
    }
  }

}


#endif
