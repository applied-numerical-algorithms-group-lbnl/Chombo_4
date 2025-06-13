#ifndef _MOMENTS_H_
#define _MOMENTS_H_

#include "Real.hpp"

#include "paramstream.hpp"
#include "MomentIterator.hpp"

namespace hoeb
{
/** Class for holding a collection of moment information on an element
 * size is compile time known, an supports living on a GPU
 */
template <unsigned Dim, unsigned Order>
class Moments
{
public:
  /// the known size of a moment
  HOSTDEVICE
  constexpr static unsigned
  size()
  {
    return MomentIterator<Dim, Order>::size();
  }

  // Constructor
  HOSTDEVICE Moments() = default;

  /// Constructor that builds from an array
  HOSTDEVICE
  Moments(std::array<Real, Moments<Dim, Order>::size()>&& a_data) : m_moments(std::move(a_data)) {}

  /// Destructor
  HOSTDEVICE ~Moments() = default;

  /// Copy constructor
  HOSTDEVICE Moments(const Moments&) = default;

  /// Move constructor
  HOSTDEVICE Moments(Moments&&) = default;

  /// Assignment constructor
  HOSTDEVICE Moments& operator=(const Moments&) = default;

  /// Move assignment constructor
  HOSTDEVICE Moments& operator=(Moments&&) = default;

  using iterator = typename std::array<Real, MomentIterator<Dim, Order>::size()>::iterator;
  using const_iterator =
      typename std::array<Real, MomentIterator<Dim, Order>::size()>::const_iterator;

  /// standard library style iterators. Allow range-based for loops
  HOSTDEVICE
  iterator
  begin()
  {
    return m_moments.begin();
  }

  ///
  HOSTDEVICE
  iterator
  end()
  {
    return m_moments.end();
  }

  ///
  HOSTDEVICE
  const_iterator
  begin() const
  {
    return m_moments.begin();
  }

  ///
  HOSTDEVICE
  const_iterator
  end() const
  {
    return m_moments.cend();
  }

  //! Retrieve the moment from the index
  //! \params a_index The lexigraphical index to the data
  HOSTDEVICE
  Real&
  operator[](int a_i)
  {
    return m_moments[a_i];
  };

  //! Retrieve the moment from the index
  //! \params a_index The lexigraphical index to the data
  HOSTDEVICE
  const Real&
  operator[](int a_i) const
  {
    return m_moments[a_i];
  };

  //! Retrieve the moment from the index
  //! \params a_index The multi-index that's needed
  HOSTDEVICE
  Real&
  operator[](const stc::IVec<Dim>& a_index)
  {
    return m_moments[MomentIterator<Dim, Order>::indexOf(a_index)];
  };

  //! Retrieve the moment from the index
  //! \params a_index The multi-index that's needed
  HOSTDEVICE
  const Real&
  operator[](const stc::IVec<Dim>& a_index) const
  {
    return m_moments[MomentIterator<Dim, Order>::indexOf(a_index)];
  };

  /// add compoenentwise
  HOSTDEVICE Moments<Dim, Order>&
  operator+=(const Moments<Dim, Order>& a_input)
  {
    for (MomentIterator<Dim, Order> it; it.ok(); ++it)
    {
      (*this)[it()] += a_input[it()];
    }
    return *this;
  }

  /// add compoenentwise
  HOSTDEVICE Moments<Dim, Order>&
  operator-=(const Moments<Dim, Order>& a_input)
  {
    for (MomentIterator<Dim, Order> it; it.ok(); ++it)
    {
      (*this)[it()] -= a_input[it()];
    }
    return *this;
  }

  /// multiply each component by factor
  HOSTDEVICE Moments<Dim, Order>&
  operator*=(const Real& a_factor)
  {
    for (MomentIterator<Dim, Order> it; it.ok(); ++it)
    {
      (*this)[it()] *= a_factor;
    }
    return *this;
  }

  /// multiply each component by factor
  HOSTDEVICE Moments<Dim, Order>&
  operator/=(const Real& a_factor)
  {
    for (MomentIterator<Dim, Order> it; it.ok(); ++it)
    {
      (*this)[it()] /= a_factor;
    }
    return *this;
  }

  /// sum all elements
  HOSTDEVICE Real
  sum() const
  {
    Real sum = 0.;
    for (const auto& momVal : (*this))
    {
      sum += momVal;
    }
    return sum;
  }

  /// set all elements uniformly
  HOSTDEVICE void
  setVal(Real a_val = 0.)
  {
    for (auto& momVal : m_moments)
    {
      momVal = a_val;
    }
  }

  /// moment centering, expressed as an offset relative to the cell center
  HOSTDEVICE stc::RVec<Dim>
  getCenter() const
  {
    return m_center;
  }

  /// moment centering, expressed as an offset relative to the cell center
  HOSTDEVICE void
  setCenter(stc::RVec<Dim> a_cent)
  {
    m_center = a_cent;
  }

  /// make a copy of this moment into a one dimension larger moment
  HOSTDEVICE Moments<Dim + 1, Order>
  promoteDimension(const int a_newDir, const Real a_dirCenter = 0.) const
  {
    assert(a_newDir >= 0);
    assert(a_newDir < Dim + 1);

    Moments<Dim + 1, Order> newMoms;
    newMoms.setCenter(expandDimension(m_center, a_newDir));
    for (MomentIterator<Dim + 1, Order> momIt; momIt.ok(); ++momIt)
    {
      stc::IVec<Dim + 1> index = momIt.momentIndex();
      stc::IVec<Dim> indexLower = contractDimension(index, a_newDir);
      newMoms[momIt()] = std::pow(a_dirCenter, index[a_newDir]) * (*this)[indexLower];
    }
    return newMoms;
  }

  /// return the associated moment iterator, now I don't have to think about template arguments
  HOSTDEVICE
  [[nodiscard]] MomentIterator<Dim, Order>
  getMomentIterator() const
  {
    MomentIterator<Dim, Order> momIt;
    return momIt;
  }

  /** Return new moments shifted from their current center to center + a_shift
   * This works with moments of any existing size and centering
   * @param a_shift The new relative location for moments to be centered about.
   * @return        The moments about the new relative center
   */
  HOSTDEVICE
  [[nodiscard]] Moments<Dim, Order>
  getShiftedMoments(const stc::RVec<Dim> a_shift) const
  {
    Moments<Dim, Order> mdx(*this); // copy of our moment data
    mdx.m_center += a_shift;
    auto qit = getMomentIterator();     // iterator over this input
    auto pit = mdx.getMomentIterator(); // iterator over output
    for (pit.reset(); pit.ok(); ++pit)
    {
      auto p = pit.momentIndex(); // index
      Real& mp = mdx[pit()];
      for (qit.reset(); qit.ok(); ++qit)
      {
        auto q = qit.momentIndex();
        if (!(q <= p) || (q == p)) continue;
        // Otherwise, accumulate the shifted values
        const Real& mq = (*this)[qit()];
        Real pcq = 1;
        Real dxq = 1;
        for (int d = 0; d < Dim; ++d)
        {
          int pd = p[d];
          int qd = q[d];
          int nck = nCk(pd, qd);
          pcq *= nck;
          dxq *= pow(-a_shift[d], pd - qd);
        }
        mp += pcq * dxq * mq;
      }
    }
    return mdx;
  }

  /** Return new moments scaled to a specified size
   * Only works with normalized moments, ie, a unit length dx and centered about zero
   * @param a_dx  The new length scale for the moments
   * @param a_dxOrig  The original length scale for the moments
   * @return      The moments with the new dx
   */
  HOSTDEVICE
  [[nodiscard]] Moments<Dim, Order>
  getScaledMoments(const stc::RVec<Dim> a_dx, const stc::RVec<Dim> a_dxOrig = stc::make_RVec<Dim>::unit()) const
  {
    assert(m_center == stc::make_RVec<Dim>::zero());
    stc::RVec<Dim> scale = a_dx/a_dxOrig;
    Moments<Dim, Order> mdx(*this); // copy of our moment data
    for (auto mit = getMomentIterator(); mit.ok(); ++mit)
    {
      auto p = mit.momentIndex() + stc::make_IVec<Dim>::unit();
      mdx[mit()] *= hoeb::ipow<Dim>(scale, p);
    }
    return mdx;
  }

  /** Return the derivative of this moment
   * @param a_deriv  The multi-dimension derivative to evaluate
   * @return         The derivative of the moments
   */
  HOSTDEVICE
  [[nodiscard]] Moments<Dim, Order>
  getDerivative(const stc::IVec<Dim> a_deriv) const
  {
    assert(a_deriv >= 0);
    Moments<Dim, Order> mdx;
    mdx.setVal(0);
    for (auto mit = getMomentIterator(); mit.ok(); ++mit)
    {
      const stc::IVec<Dim> dvIdx = mit.momentIndex() - a_deriv;
      if ((dvIdx >= stc::make_IVec<Dim>::zero())
          && (stc::sum(dvIdx) <= Order)) // non-zero derivatives
      {
        auto polyCoef = 1;
        for (int d = 0; d != Dim; d++)
        {
          polyCoef *= hoeb::fallingFactorial(mit.momentIndex()[d], a_deriv[d]);
        }
        mdx[mit()] = (Real)polyCoef * (*this)[dvIdx];
      }
    }
    return mdx;
  }

private:
  // the moment collection
  std::array<Real, MomentIterator<Dim, Order>::size()> m_moments;
  // moment centering, relative to cell center
  stc::RVec<Dim> m_center = stc::make_RVec<Dim>::zero();
};

/// print operator
template <unsigned Dim, unsigned Order>
std::ostream&
operator<<(std::ostream& a_os, const hoeb::Moments<Dim, Order>& a_moms)
{
  MomentIterator<Dim, Order> momIt;
  a_os << "Moments: dim=" << Dim << " order=" << Order << " size=" << momIt.size() << std::endl;
  for (momIt.reset(); momIt.ok(); ++momIt)
  {
    a_os << momIt() << "\t" << momIt.momentIndex() << "\t" << a_moms[momIt()] << std::endl;
  }
  return a_os;
}

/******************************************************************************/
/** utilities for generating moments of common varieties
 */
template <unsigned Dim, unsigned Order>
struct make_moments
{
  /** generate moments with all zero values
   *  \f$ m^{\vec{p}}(\vec{x}) = 0 \f$
   */
  HOSTDEVICE static constexpr Moments<Dim, Order>
  zero() noexcept
  {
    Moments<Dim, Order> zeroMoms;
    for (auto& a : zeroMoms) a = 0;
    return zeroMoms;
  }

  /** generate moments all of Nan, should make math using this fail fast
   */
  HOSTDEVICE static constexpr Moments<Dim, Order>
  invalid() noexcept
  {
    Moments<Dim, Order> ivalMoms;
    for (auto& a : ivalMoms) a = std::numeric_limits<Real>::quiet_NaN();
    return ivalMoms;
  }

  /** generate regular moments.
   *  The moments determined by \f$ m^{\vec{p}}(\vec{x}) = \int^{h/2}_{-h/2} (h x - a)^\vec{p} d\vec{x} \f$ with defaults of a unit cell about zero
   * @param a_dx       The length of the cell the moment is about \f$ \vec{h} \f$
   * @param a_relCent  The center location the moment is evaluated about \f$ \vec{a} \f$
   */
  HOSTDEVICE static constexpr Moments<Dim, Order>
  regular(const stc::RVec<Dim> a_dx = stc::make_RVec<Dim>::unit(),
          const stc::RVec<Dim> a_relCent = stc::make_RVec<Dim>::zero()) noexcept
  {
    assert(a_dx > stc::make_RVec<Dim>::zero());
    std::array<Real, MomentIterator<Dim, Order>::size()> regMoms;
    // hard code first moment
    regMoms[0] = stc::product(a_dx);
    // all others
    MomentIterator<Dim, Order> momIt;
    for (++momIt; momIt.ok(); ++momIt)
    {
      stc::IVec<Dim> index = momIt.momentIndex();
      stc::IVec<Dim> indexP = index + 1;
      // Real moment = (std::pow(0.5 * a_dx, indexP.sum()) - std::pow(-0.5 * a_dx, indexP.sum())) /
      // indexP.product();
      Real moment = 1;
      for (int dir = 0; dir != Dim; dir++)
      {
        // moment *= std::pow(a_dx/2., indexP[dir]) * (1 - std::pow(-1, indexP[dir])) / indexP[dir];
        moment *= (std::pow(a_dx[dir] + 2. * a_relCent[dir], indexP[dir]) * std::pow(-1, index[dir])
                   + std::pow(a_dx[dir] - 2. * a_relCent[dir], indexP[dir]))
                  / (std::pow(2, indexP[dir]) * indexP[dir]);
      }
      regMoms[momIt()] = moment;
    }
    return Moments<Dim, Order>(std::move(regMoms));
  }

  /** regular flattened moments for unit size cell about cell center
   * defaults to normalized moments about the center of a of a unit cell
   */
  HOSTDEVICE static constexpr Moments<Dim, Order>
  flattenedRegular(const unsigned a_flatDir,
                   const stc::RVec<Dim> a_dx = stc::make_RVec<Dim>::unit(),
                   const stc::RVec<Dim> a_relCent = stc::make_RVec<Dim>::zero()) noexcept
  {
    assert(a_dx > stc::make_RVec<Dim>::zero());
    std::array<Real, MomentIterator<Dim, Order>::size()> regFlatMoms;
    // hard coding first moment
    regFlatMoms[0] = 1;
    // defining regular moments of lower dimension
    Moments<Dim - 1, Order> regMomsLowerDim = make_moments<Dim - 1, Order>::regular();
    // iterator over the moment
    MomentIterator<Dim, Order> momIt;
    for (++momIt; momIt.ok(); ++momIt)
    {
      stc::IVec<Dim> index = momIt.momentIndex();
      Real moment = 0;
      // a_flatDir is the direction in which the moment is to be flattened. If index[a_flatDir] is
      // non-zero, that moment =0
      if (index[a_flatDir] != 0)
      {
        moment = 0;
      }
      else
      {
        // otherwise we reduce the moment index of the higher dimension to the moment index in the
        // lower dimension by eleminating the entry at index a_flatDir
        stc::IVec<Dim - 1> indexLower = contractDimension(index, a_flatDir);
        // flattened moment at this index is now the regular moment at the truncated index
        moment = regMomsLowerDim[indexLower];
      }
      regFlatMoms[momIt()] = moment;
    }
    return Moments<Dim, Order>(std::move(regFlatMoms));
  }
};

} // namespace hoeb
#endif
