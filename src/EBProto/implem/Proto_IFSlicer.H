#ifndef _IFSLICER_H_
#define _IFSLICER_H_

#include <map>
using std::vector;

#include "Proto_IndexTM.H"
#include "Proto_BaseIF.H"
#include "Proto_Notation.H"

using std::endl;
using std::shared_ptr;
namespace Proto
{

  template <int dim> class IFSlicer
  {
  public:

    /// null constructor
    /**
       null constructor
    */
    inline IFSlicer();

    /// copy constructor
    /**
       copy constructor
    */
    inline
    IFSlicer(const IFSlicer<dim> & a_IFSlicer);
    
    /// makes a reduced IFSlicer
    /**
       makes a reduced IFSlicer
    */
    inline
    IFSlicer(const shared_ptr<IFSlicer<dim+1> > & a_slicer,     
             const int                          & a_fixedComp,  
             const double                       & a_fixedValue);

    inline
    void setParams(const shared_ptr<IFSlicer<dim+1> > & a_slicer,     
                   const int                          & a_fixedComp,  
                   const double                       & a_fixedValue);

    /// destructor
    /**
       destructor
    */
    inline
    virtual ~IFSlicer();

    /// Return the partial derivative evaluated at a_point.
    /**
       Represent the partial derivative as a list of integers.
       Return the partial derivative evaluated at a_point.
    */
    inline
    virtual double value(const IndexTM<int, dim> & a_partialDerivative,
                         const IndexTM<double,dim> & a_point)const;

  protected:
    // The higher dimensional IFSlicer being further sliced
    shared_ptr<IFSlicer<dim + 1> >  m_IFSlicer;

    // The component direction being held fixed and its value
    int  m_fixedComp;
    double m_fixedValue;

  private:
    //
    // These functions are made private to prevent use of the default
    // functions provided by the C++ compiler.
    //
    IFSlicer<dim>& operator=(const IFSlicer<dim>&);

  };

  template<> class IFSlicer<DIM>
  {
  public:
    /// null constructor
    /**
       null constructor
    */
    inline
    IFSlicer();

    /// copy constructor
    /**
       copy constructor
    */
    inline
    IFSlicer(const IFSlicer<DIM> & a_IFSlicer);

    /// if dim == DIM, construct from BaseIF
    /**
       if dim == DIM, construct from BaseIF
    */
    inline
    IFSlicer(const BaseIF & a_implicitFunction);

    /// destructor
    /**
       destructor
    */
    inline
    virtual ~IFSlicer();

    /// Return the partial derivative evaluated at a_point.
    /**
       Represent the partial derivative as a list of integers.
       Return the partial derivative evaluated at a_point.
    */
    inline
    virtual double value(const IndexTM<int, DIM> & a_partialDerivative,
                       const IndexTM<double,DIM> & a_point) const;

  protected:
    // The implicit function being sliced.  At DIM this is just a
    // container class to start the slicing process.
    BaseIF* m_implicitFunction;

  private:
    //
    // These functions are made private to prevent use of the default
    // functions provided by the C++ compiler.
    //
    IFSlicer<DIM>& operator=(const IFSlicer<DIM>&);
  };

}

#include "Proto_IFSlicerImplem.H"

#endif
