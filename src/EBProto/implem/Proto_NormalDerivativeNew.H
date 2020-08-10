
#ifndef _Proto_NORMALDERIVATIVENEW_H_
#define _Proto_NORMALDERIVATIVENEW_H_



#include "Proto_Box.H"
#include "Proto_IFSlicer.H"

namespace Proto
{

/// This computes the derivatives of the normal of a sliced implicit function
/**
   This computes the derivatives of the normal of a sliced implicit function
*/
  template <int dim> class NormalDerivativeNew
  {
  public:
    typedef IndexTM<int,dim>  IvDim;
    typedef IndexTM<double,dim> RvDim;

    // These are maps to calculate partial derivatives, based on
    // multi-index IvDim, and of the normal, with components RvDim
    typedef map<IvDim,double,LexLT<IvDim> > ScalarPowerMap;
    typedef map<IvDim,RvDim,LexLT<IvDim> > NormalDerivativeMap;

    /// Null constructor
    /**
       Null constructor
    */
    NormalDerivativeNew()
    {;}

    /// Destructor
    /**
       Destructor
    */
    virtual ~NormalDerivativeNew()
    {;}

    /// Evaluate the normal of a BaseIF subclass
    /**
       Evaluate the normal.
       a_direction specifies which component of the normal use, 
       a_point is the point in space to evaluate the derivative, and 
       a_ifSlicer is a sliced function whose gradient is the normal.
    */
    double normal(const int           & a_direction,
                const RvDim         & a_point,
                const IFSlicer<dim> * a_implicitFunction) const
    {
      // Compute and store the magnitude of the gradient of the function.
      double lenGradIF = 0.0;
      double value = 0.0;

      for (int idir = 0; idir < dim; idir++)
      {
        double dxd = a_implicitFunction->
          value(BASISV_TM<int,dim>(idir), a_point);
        if (idir == a_direction)
          value = dxd;
        lenGradIF += dxd*dxd;
      }

      return value / sqrt(lenGradIF);
    }

    /// Evaluate all derivatives of the normal of an IFSlicer class
    /**
       Evaluate all derivatives of the normal, 
       a_maxP specifies the max sum of derivatives to 
       take in any coordinate direction, 
       a_point is the point in space to evaluate the derivatives, and 
       a_ifSlicer is a sliced function whose gradient is the normal.
    */
    NormalDerivativeMap calculateAll(const int& a_maxP,
                                     const RvDim& a_point,
                                     const IFSlicer<dim>* a_implicitFunction)
    {
      // pout() << "Normal derivatives<" << dim << "> for x=" << a_point << endl;

      ScalarPowerMap phiDerivs; // partials of \phi
      ScalarPowerMap lDerivs; // partials of l = |\grad \phi|
      NormalDerivativeMap nDerivs; // partials of the normal, n = \grad \phi / l

      // We only want the powers that would be less than
      // the power in each direction and less than the total
      Box b(Point::Zeroes(), a_maxP*Point::Unit());

      for (auto bit = b.begin(); bit != b.end(); ++bit)
      {
        IvDim pow;
        Point pt = *bit;
        for(int idir = 0; idir < dim; idir++)
        {
          pow[idir] =pt[idir];
        }
        if (pow.sum() <= a_maxP)
          nDerivs[pow] = RvDim::Zero;
      }

      // Calculate the phi derivatives we'll need
      calculatePhiDerivs(phiDerivs, a_point, a_implicitFunction, a_maxP+1);

      // In order of increasing derivatives, in slices
      typename NormalDerivativeMap::const_iterator it;
      for (it = nDerivs.begin(); it != nDerivs.end(); ++it)
      {
        IvDim p = it->first;
        addLDeriv(lDerivs, phiDerivs, nDerivs, p);
        addNDeriv(nDerivs, phiDerivs, lDerivs, p);
      }
      /*
        for (it = nDerivs.begin(); it != nDerivs.end(); ++it)
        pout() << "New algorithm Normal derivative: p=" << it->first << ", " << 
        it->second << endl;
      */

      return nDerivs;
    }

  protected:

    void calculatePhiDerivs(
      ScalarPowerMap& a_phiDerivs, 
      const RvDim& a_point,
      const IFSlicer<dim>* a_implicitFunction,
      const int& a_maxPartial)
    {
      Box b(Point::Zeroes(), a_maxPartial*Point::Unit());
      BoxIterator bit(b);
      // We only want the powers that would be less than
      // the power in each direction and less than the total
      for (auto bit  = b.begin(); bit != b.end(); ++bit)
      {
        Point pt = *bit;
        if (pt.sum() <= a_maxPartial)
        {
          IvDim deriv;
          for(int idir = 0; idir < dim; idir++)
          {
            deriv[idir]= pt[idir];
          }
          
          a_phiDerivs[deriv] = a_implicitFunction->value(deriv, a_point);
        }
      }
    }


    void addLDeriv(
      ScalarPowerMap& a_lDerivs, 
      ScalarPowerMap& a_phiDerivs, 
      NormalDerivativeMap& a_nDerivs, 
      const IvDim& a_deriv)
    {
      // If this is the first call, just calc l = |\grad phi|
      if (a_deriv == IvDim::Zero)
      {
        RvDim dphi;
        double len = 0;
        for (int d=0; d < dim; ++d)
        {
          double val = a_phiDerivs[BASISV_TM<int,dim>(d)]; 
          dphi[d] = val;
          len += val*val;
        }
        a_lDerivs[a_deriv] = sqrt(len);
        return;
      }

      // Otherwise, assume we are moving in increasing x, then y...
      int dir;
      for (dir=0; dir < dim-1; ++dir)
        if (a_deriv[dir] > 0)
          break;

      double dl = 0;
      IvDim q = a_deriv - BASISV_TM<int,dim>(dir);
      typename NormalDerivativeMap::const_iterator it;
      for (it = a_nDerivs.begin(); it != a_nDerivs.end(); ++it)
      {
        IvDim r = it->first;
        if (r <= q)
        {
          double coef = nChoosek(q,r);
          RvDim dn = a_nDerivs[r];
          for (int d=0; d < dim; ++d)
          {
            IvDim dphip = q - r + BASISV_TM<int,dim>(d) + 
              BASISV_TM<int,dim>(dir);
            dl += coef * dn[d] * a_phiDerivs[dphip];
          }
        }
      }
      a_lDerivs[a_deriv] = dl;
    }


    void addNDeriv(
      NormalDerivativeMap& a_nDerivs, 
      ScalarPowerMap& a_phiDerivs, 
      ScalarPowerMap& a_lDerivs, 
      const IvDim& a_deriv)
    {
      RvDim dn = RvDim::Zero; 
      typename NormalDerivativeMap::const_iterator it;
      for (it = a_nDerivs.begin(); it != a_nDerivs.end(); ++it)
      {
        IvDim q = it->first;
        if ((q <= a_deriv) && (q != a_deriv))
        {
          double coef = nChoosek(a_deriv,q);
          dn += (coef * a_lDerivs[a_deriv-q]) * it->second;
        }
      }

      double l = a_lDerivs[IvDim::Zero];
      RvDim dphi;
      for (int d=0; d < dim; ++d)
        dphi[d] = a_phiDerivs[a_deriv + BASISV_TM<int,dim>(d)];
      a_nDerivs[a_deriv] = (dphi - dn) / l;
    }


    // Calculates the multi-index binomial coefficient,
    // product of "n choose k" for each index
    static int nChoosek(const IvDim& a_n, const IvDim& a_k)
    {
      int numer = 1;
      int denom = 1;

      IvDim diff = a_n - a_k;
      PR_assert(diff >= IvDim::Zero);
      for (int d = 0; d < dim; ++d)
      {
        for (int i=diff[d]+1; i <= a_n[d]; ++i)
          numer *= i;
        for (int i=2; i <= a_k[d]; ++i)
          denom *= i;
      }
  
      return numer / denom;
    }

  private:
    NormalDerivativeNew(const NormalDerivativeNew& a_input);

    void operator=(const NormalDerivativeNew& a_input);

  };

}

#endif
