
#include <iostream>
using std::cerr;

#ifndef __EBEExactSolutions__
#define __EBEExactSolutions__

template <unsigned int order>
class BaseExactSolution
{
public:
  typedef       IndexTM<int, DIM>            IvDIM;
  typedef HostIrregData<CELL,IndMomDIM , 1>  VoluData;

  
  virtual ~BaseExactSolution()
  {
  }
  
  BaseExactSolution()
  {
  }

  ///
  Real operator()(const EBGraph & a_graph,
                  const Real    & a_dx,
                  const VoluData& a_voldat,
                  const   & a_cell)
  {
    Real retval = 0;
    RealVect cornerloc;
    for(int idir = 0; idir < DIM; idir++)
    {
      cornerloc[idir] = a_dx*a_cell[idir];
    }
    IndexedMoments<DIM, order> moments;
    if(a_graph.isRegular(a_cell))
    {
      moments.setToRegular(a_dx);
    }
    else if(a_graph.isCovered(a_cell))
    {
      return 0;
    }
    else
    {
      Point pt = ProtoCh::getPoint(a_cell);
      moments = a_voldata(pt, 0);
    }
      
    Real volu = moments[IvDIM::Zero]
    for(MomentIterator<DIM, order> momit; momit.ok(); ++momit)
    {
      IvDIM  power = momit();
      Real   deriv = getDerivative(power, cornerloc);
      Real   denom = (Real)(deriv.factorial());
      Real   momen = moments[power];
      Real   incre = deriv*momen/denom;
      retval += incre;
    }
    Real minvol = 1.0e-10;
    if(volu > minvol)
    {
      retval /= volu;
    }
    return retval;   
  }

  virtual Real getDerivative(IvDIM    a_power,
                             RealVect a_loc) = 0;
};

///class for always returning zero (for bcs as opposed to being for testing convergence rates)
template <unsigned int order>
class ZeroSolution: public BaseSolution<order>
{
public:
  
  typedef       IndexTM<int, DIM>            IvDIM;
  typedef HostIrregData<CELL,IndMomDIM , 1>  VoluData;
  
  virtual Real getDerivative(IvDIM    a_power,
                                    RealVect a_loc)
  {
    return 0.;
  }

  virtual ~ZeroSolution()
  {
  }
  
  ZeroSolution()
  {
  }

};

#endif
