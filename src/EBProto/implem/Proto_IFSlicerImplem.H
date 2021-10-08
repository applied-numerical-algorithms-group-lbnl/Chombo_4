
#ifndef _Proto_IFSLICERIMPLEM_H_
#define _Proto_IFSLICERIMPLEM_H_

namespace Proto
{

// empty constructor
  template <int dim> IFSlicer<dim>::IFSlicer()
  {
    //set to undefined values
    m_IFSlicer     = NULL;
    m_fixedComp  = LARGEINTVAL;
    m_fixedValue = LARGEREALVAL;
  }

//copy constructor
  template <int dim >
  IFSlicer<dim>::
  IFSlicer(const IFSlicer<dim> & a_IFSlicer)
  {
    setParams(a_IFSlicer.m_IFSlicer,
              a_IFSlicer.m_fixedComp,
              a_IFSlicer.m_fixedValue);
  }

  template <int dim> 
  IFSlicer<dim>::
  IFSlicer(const shared_ptr<IFSlicer<dim+1> > & a_IFSlicer,     
           const int                          & a_fixedComp,  
           const double                       & a_fixedValue)
  {
    setParams(a_IFSlicer,
              a_fixedComp,
              a_fixedValue);
  }

  template <int dim> 
  void IFSlicer<dim>::
  setParams(const shared_ptr<IFSlicer<dim+1> > & a_IFSlicer,     
            const int                          & a_fixedComp,  
            const double                       & a_fixedValue)
  {
    // copy the inputs
    m_IFSlicer   = a_IFSlicer;
    m_fixedComp  = a_fixedComp;
    m_fixedValue = a_fixedValue;
  }

// Destructor
  template <int dim> IFSlicer<dim>::~IFSlicer()
  {
    //lower dimension doesn't own the objects referred to by any pointer.
  }

  template<int dim >double IFSlicer<dim>::value(const IndexTM<int, dim> & a_partialDerivative,
                                              const IndexTM<double,dim> & a_point) const
  {
    double retval = LARGEREALVAL;

    //use fixedComp and fixedValue to pull partialDerivative and point into one higher dimension
    IndexTM<int, dim + 1> partialDerivative;
    IndexTM<double,dim + 1> point;
    for (int idir = 0; idir < dim +1; ++idir)
    {
      if (idir < m_fixedComp)
      {
        partialDerivative[idir] = a_partialDerivative[idir];
        point            [idir] = a_point            [idir];
      }
      else if (idir > m_fixedComp)
      {
        partialDerivative[idir] = a_partialDerivative[idir - 1];
        point            [idir] = a_point            [idir - 1];
      }
      else
      {
        partialDerivative[m_fixedComp] = 0;
        point            [m_fixedComp] = m_fixedValue;
      }
    }

    //evaluate at one higher dimension
    retval = m_IFSlicer->value(partialDerivative,point);

    return retval;
  }


// empty constructor (dim == DIM)
  IFSlicer<DIM>::IFSlicer()
  {
    m_implicitFunction = NULL;
  }

//copy constructor
  IFSlicer<DIM>::IFSlicer(const IFSlicer<DIM> & a_IFSlicer)
  {
    m_implicitFunction = a_IFSlicer.m_implicitFunction->newImplicitFunction();
  }

// constructor (dim == DIM)
  IFSlicer<DIM>::IFSlicer(const BaseIF & a_implicitFunction)
  {
    m_implicitFunction = a_implicitFunction.newImplicitFunction();
  }

// Destructor (dim == DIM)
  IFSlicer<DIM>::~IFSlicer()
  {
    if (m_implicitFunction != NULL)
    {
      delete m_implicitFunction;
    }
  }

  double IFSlicer<DIM>::value(const IndexTM<int, DIM> & a_partialDerivative,
                            const IndexTM<double,DIM> & a_point) const
  {
    return m_implicitFunction->value(a_partialDerivative,a_point);
  }


}

#endif
