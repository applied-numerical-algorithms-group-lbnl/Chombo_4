#ifndef EXACTSOLUTIONS_H_
#define EXACTSOLUTIONS_H_

#include "MultiIndexMath.hpp"
#include "Moments.hpp"

#include <autodiff/forward/dual.hpp>

namespace hoeb
{

/** Print a type at compile-time
 *  Depends on compiler support, works for gcc and clang
 */
template <typename T>
void
print_type()
{
  std::cout << __PRETTY_FUNCTION__ << std::endl;
}

namespace exactSolution
{
/*******************************************************************************/
/** AutoDiff Real value type
 *  @tparam Order  Highest allowed derivative by AutoDiff for this type
 */
template <unsigned Order>
using ebDual = autodiff::HigherOrderDual<Order, hbr::Real>;

/** Vector of AutoDiff Real value types
 *  @tparam Order  Highest allowed derivative by AutoDiff for this type
 *  @tparam Dim    Length of vector
 */
template <unsigned Order, unsigned Dim>
using DVec = stc::Vector<ebDual<Order>, Dim>;

/*******************************************************************************/
/** Return type dispatch for autodiff::wrt
 * @tparam Diff    Number of derivatives for autodiff to evaluate with respect to
 * @tparam Order   Number of supported derivatives of the autodiff Dual type
 */
template <unsigned Diff, unsigned Order, typename... WrtArgs>
struct diffOrder;

/// no derivatives
template <unsigned Order, typename... WrtArgs>
struct diffOrder<0, Order, WrtArgs...>
{
  using wrtType = autodiff::Wrt<WrtArgs...>;
};

/// Recursively add an ebDual<Order>& type for each derivative needed
template <unsigned Diff, unsigned Order, typename... WrtArgs>
struct diffOrder
{
  using wrtType = typename diffOrder<Diff - 1, Order, WrtArgs..., ebDual<Order>&>::wrtType;
};

/*******************************************************************************/
/** Functor that turns a stc::IVec into a autodiff::wrt type object
 *  Uses recursive template meta-programming magic to work for any size types
 *  @tparam Diff   The number of actual derivatives
 *  @tparam Order  The allowable derivatives of the ebDual type
 *  @tparam Dim    The number of vector dimensions
 */
template <unsigned Diff, unsigned Order, unsigned Dim>
struct getAutoDiffWrt
{
  /** Returns a autodiff::wrt type object
   *  Uses recursive template meta-programming magic to work for any size types
   *  @param a_deriv  Derivative to evaluate
   *  @param a_loc    Location to evaluate about
   */
  HOSTDEVICE
  typename diffOrder<Diff, Order>::wrtType operator()(const stc::IVec<Dim> a_deriv,
                                                      DVec<Order, Dim>& a_loc);
};

// no derivatives
template <unsigned Order, unsigned Dim>
struct getAutoDiffWrt<0, Order, Dim>
{
  HOSTDEVICE
  typename diffOrder<0, Order>::wrtType
  operator()(const stc::IVec<3> a_deriv, DVec<Order, Dim>& a_loc)
  {
    assert(a_deriv[stc::indexMinElem(a_deriv)] >= 0);
    assert(stc::sum(a_deriv) == 0);
    typename diffOrder<0, Order>::wrtType w{};
    return w;
  }
};

// 1st derivative
template <unsigned Order, unsigned Dim>
struct getAutoDiffWrt<1, Order, Dim>
{
  HOSTDEVICE
  typename diffOrder<1, Order>::wrtType
  operator()(const stc::IVec<Dim> a_deriv, DVec<Order, Dim>& a_loc)
  {
    assert(a_deriv[stc::indexMinElem(a_deriv)] >= 0);
    assert(stc::sum(a_deriv) <= 1);
    typename diffOrder<1, Order>::wrtType w{ a_loc[stc::indexMaxElem(a_deriv)] };
    return w;
  }
};

// Nth derivatives
template <unsigned Diff, unsigned Order, unsigned Dim>
HOSTDEVICE
typename diffOrder<Diff, Order>::wrtType
getAutoDiffWrt<Diff, Order, Dim>::operator()(const stc::IVec<Dim> a_deriv, DVec<Order, Dim>& a_loc)
{
  assert(a_deriv[stc::indexMinElem(a_deriv)] >= 0);
  assert(stc::sum(a_deriv) <= Diff);
  // 1st derivative of first non-zero index (lexographical ordering)
  unsigned reduceDir = Dim;
  for (unsigned dir = 0; dir != Dim; dir++)
  {
    if (a_deriv[dir] > 0)
    {
      reduceDir = dir;
      break;
    }
  }
  stc::IVec<Dim> iv1 = stc::make_IVec<Dim>::basis(reduceDir);
  auto w1 = getAutoDiffWrt<1, Order, Dim>()(iv1, a_loc);
  // N-1 derivative of remainder - recursive
  stc::IVec<Dim> ivR = a_deriv - stc::make_IVec<Dim>::basis(reduceDir);
  auto wr = getAutoDiffWrt<Diff - 1, Order, Dim>()(ivR, a_loc);
  // combine 1st derivative with N-1 derivative
  auto wt = std::tuple_cat(w1.args, wr.args);
  typename diffOrder<Diff, Order>::wrtType w{ wt };
  return w;
}

/*******************************************************************************/
/// Helper function to generate index sequence and call autodiff::at
template <typename Vect, std::size_t... idx>
HOSTDEVICE auto
getAutoDiffAtIdx(Vect& vect, std::index_sequence<idx...>)
{
  return autodiff::at(vect[idx]...);
}

/// Get the autodiff::at object for evaluation of a vector location
template <typename Vect>
HOSTDEVICE auto
getAutoDiffAt(Vect& vect)
{
  return getAutoDiffAtIdx(vect, std::make_index_sequence<Vect::size()>{});
}

/****************************************************************************/
/** Autodiff to vector function interface
 *  this little object wraps the convenience of a vector function into a expanded list of arguments
 * that Autodiff needs to see in order to do its magic
 *  @tparam Dim   Number of dimensions/arguments of the function to evaluate
 *  @tparam Func  Vector function to evaluate, with signature F(stc::Vector<T, Dim>) -> T
 */
template <unsigned Dim, typename Func>
struct autoDiffEvalVectFunc
{
  HOSTDEVICE
  autoDiffEvalVectFunc(const Func& a_func) : m_func(&a_func) {}

  /** Autodiff to vector function interface
   *  @param x  First index value for function evaluation
   *  @param y  Additional values for evaluation
   */
  template <typename T, typename... Args>
  HOSTDEVICE T
  operator()(T x, Args... y) const
  {
    static_assert(sizeof...(y) == Dim - 1);
    return (*m_func)(stc::Vector<T, Dim>{ x, y... });
  }

  // hold a non-owning pointer to the function
  const Func* m_func;
};

/****************************************************************************/
/// Methods for function evaluation
/// These use support for arbitrary derivatives using AutoDiff
/****************************************************************************/
/** Evaluates a AutoDiff compatible scalar function at a given location
 * @tparam AutoDiffOrder The maximum supported AutoDiff derivative - must be larger than the number
 * of derivatives requested
 * @tparam Func    The function type to evaluate - should be deduced by compiler
 * @tparam Dim     Dimension of the Func - should be deduced by the compiler
 * @param a_func   The function to evaluate, expects signature a_func(stc::Vector<T, Dim>) -> Real
 * @param a_loc    The vector location to evaluate the function
 * @param a_deriv  The derivative of the function to evaluate
 */
template <unsigned AutoDiffOrder, typename Func, unsigned Dim>
HOSTDEVICE Real
evalPoint(const Func& a_func, const stc::RVec<Dim> a_loc,
          const stc::IVec<Dim> a_deriv = stc::make_IVec<Dim>::zero())
{
  static_assert(AutoDiffOrder > 0);
  static_assert(Dim > 0);
  assert(AutoDiffOrder >= stc::sum(a_deriv)); // requesting a lower autodiff order will run, but
                                              // yield garbage results
  assert(a_deriv[stc::indexMinElem(a_deriv)] >= 0);                 // derivatives must all be positive
  // make the eval location in autodiff terms
  DVec<AutoDiffOrder, Dim> dualX{ a_loc };
  auto evalPt = getAutoDiffAt(dualX);
  // print_type<decltype(evalPt)>(); // print type
  // wrap the vector valued function something autodiff can use
  autoDiffEvalVectFunc<Dim, Func> adFunc(a_func);
  // evaluate the point
  auto dv = exactSolution::getAutoDiffWrt<AutoDiffOrder, AutoDiffOrder, Dim>()(a_deriv, dualX);
  auto du = autodiff::derivatives(adFunc, dv, evalPt);
  Real deriv = du[stc::sum(a_deriv)];
  assert(std::isfinite(deriv)); // breaks on NaN/Inf
  return deriv;
}

/** Evaluates the cell average quantity of a AutoDiff compatible scalar function about a given
 * location
 * @tparam AutoDiffOrder The maximum supported AutoDiff derivative - must be larger than the number
 * of moments and derivatives requested
 * @tparam Func    The function type to evaluate - should be deduced by compiler
 * @tparam Dim     Dimension of the Func - should be deduced by the compiler
 * @param a_func   The function to evaluate, expects signature a_func(stc::Vector<T, Dim>) -> Real
 * @param a_mom    The cell moments to evaluate
 * @param a_cent   The cell center to evaluate the function about
 * @param a_dx     The cell size
 * @param a_deriv  The derivative of the function to evaluate
 */
template <unsigned AutoDiffOrder, typename Func, unsigned Dim, unsigned Order>
HOSTDEVICE Real
evalAverage(const Func& a_func, const hoeb::Moments<Dim, Order>& a_mom, const stc::RVec<Dim> a_cent,
            const stc::RVec<Dim> a_dx = stc::make_IVec<Dim>::unit(),
            const stc::IVec<Dim> a_deriv = stc::make_IVec<Dim>::zero())
{
  static_assert(AutoDiffOrder > 0);
  static_assert(Order > 0);
  static_assert(Dim > 0);
  assert(AutoDiffOrder >= Order + stc::sum(a_deriv)); // requesting a lower autodiff order will run,
                                                      // but yield garbage results
  assert(a_deriv[stc::indexMinElem(a_deriv)] >= 0);                         // derivatives must all be positive
  assert(a_deriv[stc::indexMinElem(a_deriv)] >= 0);
  // make the eval location in autodiff terms
  DVec<AutoDiffOrder, Dim> dualX{ a_cent };
  auto evalPt = getAutoDiffAt(dualX);
  // print_type<decltype(evalPt)>(); // print type
  // wrap the vector valued function something autodiff can use
  autoDiffEvalVectFunc<Dim, Func> adFunc(a_func);
  // evaluate taylor series style moment values
  // TODO - not true for 1-based moments!
  // \<f\> = \Sum_q f^{(q)}(x0) * m_q
  Real retval = 0;
  for (auto momit = a_mom.getMomentIterator(); momit.ok(); ++momit)
  {
    stc::IVec<Dim> idv = momit.momentIndex();
    auto dv
        = exactSolution::getAutoDiffWrt<AutoDiffOrder, AutoDiffOrder, Dim>()(idv + a_deriv, dualX);
    auto du = autodiff::derivatives(adFunc, dv, evalPt);
    Real deriv = du[stc::sum(idv + a_deriv)];
    assert(std::isfinite(deriv)); // breaks on NaN/Inf
    // pout() << momit.momentIndex() << a_cent << " : " << deriv << " * " << a_mom[idv] << " / " <<
    // hoeb::factorial<Dim>(idv) << std::endl;
    Real dxPow = hoeb::ipow<Dim>(a_dx, idv + stc::make_IVec<Dim>::unit()); // need factor dx^p into moments
    retval += deriv * a_mom[idv] * dxPow / hoeb::factorial<Dim>(idv);
  }
  // divide by volume for cell average
  retval /= (a_mom[0]*a_dx.product());
  // pout() << " = " <<  retval << std::endl <<  std::endl;
  //  for stability, vanishingly small volumes eval to zero
  if (a_mom[0] <= hbr::Real_eps)
  {
    retval = 0.;
  }
  return retval;
}

/****************************************************************************/
/// Common operators using autodiff functions
/****************************************************************************/

/** Evaluates the gradient of a AutoDiff compatible scalar function at a given location
 * @tparam Func    The function type to evaluate - should be deduced by compiler
 * @tparam Dim     Dimension of the Func - should be deduced by the compiler
 * @param a_func   The function to evaluate, expects signature a_func(stc::Vector<T, Dim>) -> Real
 * @param a_loc    The vector location to evaluate the function
 */
template <typename Func, unsigned Dim>
HOSTDEVICE stc::RVec<Dim>
evalGradient(const Func& a_func, const stc::RVec<Dim> a_loc)
{
  stc::RVec<Dim> grad;
  for (int d = 0; d != Dim; d++)
  {
    grad[d] = evalPoint<1>(a_func, a_loc, stc::make_IVec<Dim>::basis(d));
    assert(std::isfinite(grad[d])); // breaks on NaN/Inf
  }
  return grad;
}

/** Evaluates the divergence of a AutoDiff compatible scalar function at a given location
 * @tparam Func    The function type to evaluate - should be deduced by compiler
 * @tparam Dim     Dimension of the Func - should be deduced by the compiler
 * @param a_func   The function to evaluate, expects signature a_func(stc::Vector<T, Dim>) -> Real
 * @param a_loc    The vector location to evaluate the function
 */
template <typename Func, unsigned Dim>
HOSTDEVICE Real
evalDivergence(const Func& a_func, const stc::RVec<Dim> a_loc)
{
  Real div = 0;
  for (int d = 0; d != Dim; d++)
  {
    div += evalPoint<1>(a_func, a_loc, stc::make_IVec<Dim>::basis(d));
  }
  assert(std::isfinite(div)); // breaks on NaN/Inf
  return div;
}

/** Evaluates the laplacian of a AutoDiff compatible scalar function at a given location
 * @tparam Func    The function type to evaluate - should be deduced by compiler
 * @tparam Dim     Dimension of the Func - should be deduced by the compiler
 * @param a_func   The function to evaluate, expects signature a_func(stc::Vector<T, Dim>) -> Real
 * @param a_loc    The vector location to evaluate the function
 */
template <typename Func, unsigned Dim>
HOSTDEVICE Real
evalLaplacian(const Func& a_func, const stc::RVec<Dim> a_loc)
{
  Real lap = 0;
  for (int d = 0; d != Dim; d++)
  {
    lap += evalPoint<2>(a_func, a_loc, stc::IVec<Dim>(2 * stc::make_IVec<Dim>::basis(d)));
  }
  assert(std::isfinite(lap)); // breaks on NaN/Inf
  return lap;
}

/** Evaluates the cell averaged gradient of a AutoDiff compatible scalar function at a given
 * location
 * @tparam Func    The function type to evaluate - should be deduced by compiler
 * @tparam Dim     Dimension of the Func - should be deduced by the compiler
 * @param a_func   The function to evaluate, expects signature a_func(stc::Vector<T, Dim>) -> Real
 * @param a_mom    The cell moments to evaluate
 * @param a_cent   The cell center to evaluate the function about
 * @param a_dx     The cell size
 */
template <typename Func, unsigned Dim, unsigned Order>
HOSTDEVICE stc::RVec<Dim>
evalAvgGradient(const Func& a_func, const hoeb::Moments<Dim, Order>& a_mom,
                const stc::RVec<Dim> a_cent, const stc::RVec<Dim> a_dx = stc::make_IVec<Dim>::unit())
{
  stc::RVec<Dim> grad;
  for (int d = 0; d != Dim; d++)
  {
    grad[d] = evalAverage<Order + 1>(a_func, a_mom, a_cent, a_dx, stc::make_IVec<Dim>::basis(d));
    assert(std::isfinite(grad[d])); // breaks on NaN/Inf
  }
  return grad;
}

/** Evaluates the gradient of a AutoDiff compatible scalar function at a given location
 * @tparam Func    The function type to evaluate - should be deduced by compiler
 * @tparam Dim     Dimension of the Func - should be deduced by the compiler
 * @param a_func   The function to evaluate, expects signature a_func(stc::Vector<T, Dim>) -> Real
 * @param a_mom    The cell moments to evaluate
 * @param a_cent   The cell center to evaluate the function about
 * @param a_dx     The cell size
 */
template <typename Func, unsigned Dim, unsigned Order>
HOSTDEVICE Real
evalAvgDivergence(const Func& a_func, const hoeb::Moments<Dim, Order>& a_mom,
                  const stc::RVec<Dim> a_cent, const stc::RVec<Dim> a_dx = stc::make_IVec<Dim>::unit())
{
  Real div = 0;
  for (int d = 0; d != Dim; d++)
  {
    div += evalAverage<Order + 1>(a_func, a_mom, a_cent, a_dx, stc::make_IVec<Dim>::basis(d));
  }
  assert(std::isfinite(div)); // breaks on NaN/Inf
  return div;
}

/** Evaluates the laplacian of a AutoDiff compatible scalar function at a given location
 * @tparam Func    The function type to evaluate - should be deduced by compiler
 * @tparam Dim     Dimension of the Func - should be deduced by the compiler
 * @param a_func   The function to evaluate, expects signature a_func(stc::Vector<T, Dim>) -> Real
 * @param a_mom    The cell moments to evaluate
 * @param a_cent   The cell center to evaluate the function about
 * @param a_dx     The cell size
 */
template <typename Func, unsigned Dim, unsigned Order>
HOSTDEVICE Real
evalAvgLaplacian(const Func& a_func, const hoeb::Moments<Dim, Order>& a_mom,
                 const stc::RVec<Dim> a_cent, const stc::RVec<Dim> a_dx = stc::make_IVec<Dim>::unit())
{
  Real lap = 0;
  for (int d = 0; d != Dim; d++)
  {
    lap += evalAverage<Order + 2>(a_func, a_mom, a_cent, a_dx,
                                  stc::IVec<Dim>(2 * stc::make_IVec<Dim>::basis(d)));
  }
  assert(std::isfinite(lap)); // breaks on NaN/Inf
  return lap;
}

/****************************************************************************/
/// A collection of scalar analytical autodiff compatible functors, with vector inputs
/****************************************************************************/
/** Autodiff functor for a uniform value in space
 * @tparam Dim    The number of input variables
 */
template <unsigned Dim>
class Constant
{
public:
  /** Create function \f$ \phi(\vec{x}) = c \f$
   * @param a_val  The coefficient c
   */
  Constant(const Real a_val = 1.) : m_val(a_val) {}

  /** Evaluate function at location
   * @param a_x   Vector location to evaluate
   */
  template <typename T>
  HOSTDEVICE T
  operator()(stc::Vector<T, Dim> a_vec) const
  {
    return m_val;
  }

private:
  Real m_val;
};

/** Autodiff monomial with vector form inputs
 * @tparam Dim    The number of input variables
 */
template <unsigned Dim>
class Monomial
{
public:
  /** Create function \f$ \phi(\vec{x}) = c \vec{x}^\vec{q} = c * \Sum_d x_d ^ q_d \f$
   * @param a_powers The vector monomial powers \f$ \vec{q} \f$
   * @param a_coef   The coefficient \f$ c \f$
   */
  Monomial(const stc::IVec<Dim> a_powers = stc::make_IVec<Dim>::zero(), const Real a_coef = 1.)
      : m_powers(a_powers), m_coef(a_coef)
  {
    static_assert(Dim > 0);
  }

  /** Evaluate function at location
   * @param a_x   Vector location to evaluate
   */
  template <typename T>
  HOSTDEVICE T
  operator()(const stc::Vector<T, Dim> a_vec) const
  {
    // Careful, using ipow is required here over the normal pow
    // autodiff assumes non-integer pow, which doesn't support derivatives at 0 or powers of 0
    T val = m_coef;
    for (int dir = 0; dir != Dim; dir++)
    {
      val *= hoeb::ipow(a_vec[dir], m_powers[dir]);
    }
    return val;
  }

private:
  stc::IVec<Dim> m_powers;
  Real m_coef;
};

/** Autodiff polynomial of specified terms vector form inputs
 * @tparam Dim    The number of input variables
 * @tparam Terms  The number of monomial terms in the polynomial
 */
template <unsigned Dim, unsigned Terms>
class Polynomial
{
public:
  Polynomial() = default;

  /** Create function \f$ \phi(\vec{x}) = \Sum_i c_i \vec{x}^{\vec{q}_i} \f$
   * @param a_powers The vector monomial powers \f$ \vec{q}_i \f$
   * @param a_coef   The coefficient \f$ c_i \f$
   */
  Polynomial(std::initializer_list<std::pair<stc::IVec<Dim>, Real>> a_monomialTerms)
    {
      int i=0;
      for (auto [p, c] : a_monomialTerms)
        {
          if (i >= Terms) break;
          m_monos[i] = Monomial<Dim>(p, c);
          i++;
        }

      for (int j = a_monomialTerms.size(); j <= Terms; j++)
        {
          m_monos[j] = Monomial<Dim>(stc::make_IVec<Dim>::zero(), 0);
        }
    }

  /** Set the ith term in the function function \f$ \phi(\vec{x}) = \Sum_i c_i \vec{x}^{\vec{q}_i} \f$
   *  @param a_term  The term at index \f$ i \f$ to set
   * @param a_powers The vector monomial powers \f$ \vec{q}_i \f$
   * @param a_coef   The monomial coefficient \f$ c_i \f$
   */
  void
  setTerm(int a_term, stc::IVec<Dim> a_powers, Real a_coef = 1.)
  {
    assert(a_term >= 0);
    assert(a_term < Terms);
    m_monos[a_term] = Monomial(a_powers, a_coef);
  }

  /** Evaluate function at location
   * @param a_x   Vector location to evaluate
   */
  template <typename T>
  HOSTDEVICE T
  operator()(stc::Vector<T, Dim> a_vec) const
  {
    T val = 0;
    for (auto mono : m_monos)
    {
      val += mono(a_vec);
    }
    return val;
  }

private:
  std::array<Monomial<Dim>, Terms> m_monos;
};

/** Autodiff multi-dimensional Gaussian profile
 * @tparam Dim    The number of input variables
 */
template <unsigned Dim>
class Gaussian
{
public:
  /** Create function \f$ \phi(\vec{x}) = a \exp{\Sum_d \frac{(x_d - c_d)^2}{\sigma_d^2} } \f$
   * @param a_center     The location \f$ \vec{c} \f$ to center the Gaussian about
   * @param a_amplitude  The peak value \f$ a \f$ of the Gaussian
   * @param a_deviation  The deviation \f$ \vec{\sigma}\f$ of the Gaussian.
   *                     A deviation of zero will neglect the contribution in that component
   */
  Gaussian(const stc::RVec<Dim> a_center = stc::make_RVec<Dim>::zero(), const Real a_amplitude = 1,
           const stc::RVec<Dim> a_deviation = stc::make_RVec<Dim>::unit())
      : m_center(a_center), m_amplitude(a_amplitude), m_deviation(a_deviation)
  {
  }

  /** Evaluate function at location
   * @param a_x   Vector location to evaluate
   */
  template <typename T>
  HOSTDEVICE T
  operator()(stc::Vector<T, Dim> a_x) const
  {
    T retval = 0.;
    for (int dir = 0; dir != Dim; dir++)
    {
      if (m_deviation[dir] != 0.)
      {
        retval += hoeb::ipow((T)((a_x[dir] - m_center[dir]) / m_deviation[dir]), 2) ;
      }
    }
    return m_amplitude * exp(-retval);
  }

private:
  stc::RVec<Dim> m_center;
  Real m_amplitude;
  stc::RVec<Dim> m_deviation;
};

/** Autodiff multi-dimensional Gaussian approximation with compact support
 *  This is slightly tighter than a traditional Gaussian, with hard zeros at a cutoff
 * @tparam Dim    The number of input variables
 */
template <unsigned Dim>
class PseudoGaussian
{
public:
  /** Create function \f$ \phi(\vec{x}) = a \Prod_d (\cos((x_d - c_d)*(pi / 4 L_i))^6, 0 when |x_d| >
   * L_d \f$
   * @param a_center     The location \f$ \vec{c} \f$ to center the pseudo-Gaussian about
   * @param a_amplitude  The peak value \f$ a \f$ of the pseudo-Gaussian
   * @param a_deviation  The length \f$ \vec{L}\f$ of the pseudo-Gaussian
   *                     A deviation of zero will neglect the contribution in that component
   */
  PseudoGaussian(const stc::RVec<Dim> a_center = stc::make_RVec<Dim>::zero(),
                 const Real a_amplitude = 1,
                 const stc::RVec<Dim> a_deviation = stc::make_RVec<Dim>::unit())
      : m_center(a_center), m_amplitude(a_amplitude), m_deviation(a_deviation)
  {
  }

  /** Evaluate function at location
   * @param a_x   Vector location to evaluate
   */
  template <typename T>
  HOSTDEVICE T
  operator()(stc::Vector<T, Dim> a_x) const
  {
    T retval = m_amplitude;
    for (int dir = 0; dir != Dim; dir++)
    {
      if (m_deviation[dir] != 0.)
      {
        if (abs(a_x[dir] - m_center[dir]) > m_deviation[dir])
          {
            return 0;
          }
        retval *= hoeb::ipow((T)cos((a_x[dir] - m_center[dir]) * hoeb::Pi() / (4. * m_deviation[dir])), 6);
      }
    }
    return retval;
  }

private:
  stc::RVec<Dim> m_center;
  Real m_amplitude;
  stc::RVec<Dim> m_deviation;
};

/** Autodiff product of one dimensional trigonometric sine functions
 * @tparam Dim    The number of input variables
 */
template <unsigned Dim>
class SinProduct
{
public:
  /** Create function \f$ \phi(\vec{x}) = a \Prod_d \sin((x_d - c_d) f_d) \f$ with arguments in
   * radians
   * @param a_center     The location \f$ \vec{c} \f$ to center the function about
   * @param a_amplitude  The peak value \f$ a \f$ of the function
   * @param a_frequency  The frequency \f$ \vec{f} \f$ of the sine functions
   *                     A frequency of zero will neglect the contribution in that component
   */
  SinProduct(const stc::RVec<Dim> a_center = stc::make_RVec<Dim>::zero(),
             const Real a_amplitude = 1,
             const stc::RVec<Dim> a_frequency = stc::make_RVec<Dim>::unit())
      : m_center(a_center), m_amplitude(a_amplitude), m_frequency(a_frequency)
  {
  }

  /** Evaluate function at location
   * @param a_x   Vector location to evaluate
   */
  template <typename T>
  HOSTDEVICE T
  operator()(stc::Vector<T, Dim> a_x) const
  {
    T retval = m_amplitude;
    for (int dir = 0; dir != Dim; dir++)
    {
      if (m_frequency[dir] != 0.)
      {
        retval *= sin((a_x[dir] - m_center[dir]) * m_frequency[dir]);
      }
    }
    return retval;
  }

private:
  stc::RVec<Dim> m_center;
  Real m_amplitude;
  stc::RVec<Dim> m_frequency;
};

/** Autodiff product of one dimensional trigonometric cosine functions
 * @tparam Dim    The number of input variables
 */
template <unsigned Dim>
class CosProduct
{
public:
  /** Create function \f$ \phi(\vec{x}) = a \Prod_d \cos(x_d - c_d) f_d) \f$ with arguments in
   * radians
   * @param a_center     The location \f$ \vec{c} \f$ to center the function about
   * @param a_amplitude  The peak value \f$ a \f$ of the function
   * @param a_frequency  The frequency \f$ \vec{f} \f$ of the sine functions
   *                     A frequency of zero will neglect the contribution in that component
   */
  CosProduct(const stc::RVec<Dim> a_center = stc::make_RVec<Dim>::zero(),
             const Real a_amplitude = 1,
             const stc::RVec<Dim> a_frequency = stc::make_RVec<Dim>::unit())
      : m_center(a_center), m_amplitude(a_amplitude), m_frequency(a_frequency)
  {
  }

  /** Evaluate function at location
   * @param a_x   Vector location to evaluate
   */
  template <typename T>
  HOSTDEVICE T
  operator()(stc::Vector<T, Dim> a_x) const
  {
    T retval = m_amplitude;
    for (int dir = 0; dir != Dim; dir++)
    {
      if (m_frequency[dir] != 0.)
      {
        retval *= cos((a_x[dir] - m_center[dir]) * m_frequency[dir]);
      }
    }
    return retval;
  }

private:
  stc::RVec<Dim> m_center;
  Real m_amplitude;
  stc::RVec<Dim> m_frequency;
};

/** Autodiff product of one dimensional exponential functions
 * @tparam Dim    The number of input variables
 */
template <unsigned Dim>
class ExponentialProduct
{
public:
  /** Create function \f$ \phi(\vec{x}) = a \Prod_d \exp(x_d - c_d) f_d) \f$
   * @param a_center The location \f$ \vec{c} \f$ to center the function about
   * @param a_coeff  The scaling coefficient \f$ a \f$ of the function
   * @param a_rate   The rate \f$ \vec{f} \f$ of the sine functions
   */
  ExponentialProduct(const stc::RVec<Dim> a_center = stc::make_RVec<Dim>::zero(),
                     const Real a_coeff = 1,
                     const stc::RVec<Dim> a_rate = stc::make_RVec<Dim>::unit())
      : m_center(a_center), m_coeff(a_coeff), m_rate(a_rate)
  {
  }

  /** Evaluate function at location
   * @param a_x   Vector location to evaluate
   */
  template <typename T>
  HOSTDEVICE T
  operator()(stc::Vector<T, Dim> a_x) const
  {
    T retval = m_coeff;
    for (int dir = 0; dir != Dim; dir++)
    {
      retval *= exp((a_x[dir] - m_center[dir]) * m_rate[dir]);
    }
    return retval;
  }

private:
  stc::RVec<Dim> m_center;
  Real m_coeff;
  stc::RVec<Dim> m_rate;
};

/** Autodiff sum of functions with the same input space
 * TODO
 */
template <unsigned Dim, typename Func0, typename... FuncN>
class FunctionSum
{
public:
  /** Create function \f$ \phi(\vec{x}) = \Sum_i f(\vec{x}) \f$
   */
  FunctionSum(Func0 a_func0, FuncN... a_funcN)
      : m_funcs(a_func0, a_funcN...)
  {
  }

  /** Evaluate function at location
   * @param a_x   Vector location to evaluate
   */
  template <typename T>
  HOSTDEVICE T
  operator()(stc::Vector<T, Dim> a_x) const
  {
    T retval = std::apply([&](const auto&... f) { return (f(a_x) + ...); }, m_funcs);
    return retval;
  }

private:
  std::tuple<Func0, FuncN...> m_funcs;

};

/** Autodiff product of functions with the same input space
 * TODO
 */
template <unsigned Dim, typename Func0, typename... FuncN>
class FunctionProduct
{
public:
  /** Create function \f$ \phi(\vec{x}) = \Prod_i f(\vec{x}) \f$
   */
  FunctionProduct(Func0 a_func0, FuncN... a_funcN) : m_funcs(a_func0, a_funcN...) {}

  /** Evaluate function at location
   * @param a_x   Vector location to evaluate
   */
  template <typename T>
  HOSTDEVICE T
  operator()(stc::Vector<T, Dim> a_x) const
  {
    T retval = std::apply([&](const auto&... f) { return (f(a_x) * ...); }, m_funcs);
    return retval;
  }

private:
  std::tuple<Func0, FuncN...> m_funcs;
};

/** Autodiff sum of functions with different input spaces
 * TODO
 */
template <unsigned Dim0, unsigned Dim1, typename Func0, typename Func1>
class MultiVarFunctionSum
{
public:
  static constexpr unsigned InputDim = Dim0 + Dim1;

  /** Create function \f$ \phi(\vec{x}, \vec{u}) = f(\vec{x}) + g(\vec{u}) \f$
   */
  MultiVarFunctionSum(Func0 a_func0, Func1 a_func1) : m_func0(a_func0), m_func1(a_func1) {}

  /** Evaluate function at location
   * @param a_x   Vector location to evaluate, concatenation of \f$ \vec{x}, \vec{u} \f$
   */
  template <typename T>
  HOSTDEVICE T
  operator()(stc::Vector<T, InputDim> a_x) const
  {
    stc::Vector<T, Dim0> x0;
    for (int d = 0; d != Dim0; d++)
    {
      x0[d] = a_x[d];
    }
    stc::Vector<T, Dim1> x1;
    for (int d = 0; d != Dim1; d++)
    {
      x1[d] = a_x[d + Dim0];
    }
    return m_func0(x0) + m_func1(x1);
  }

private:
  Func0 m_func0;
  Func1 m_func1;
};

/** Autodiff product of functions with different input spaces
 * TODO
 */
template <unsigned Dim0, unsigned Dim1, typename Func0, typename Func1>
class MultiVarFunctionProduct
{
public:
  static constexpr unsigned InputDim = Dim0 + Dim1;

  /** Create function \f$ \phi(\vec{x}, \vec{u}) = f(\vec{x}) g(\vec{u})  \f$
   */
  MultiVarFunctionProduct(Func0 a_func0, Func1 a_func1) : m_func0(a_func0), m_func1(a_func1) {}

  /** Evaluate function at location
   * @param a_x   Vector location to evaluate, concatenation of \f$ \vec{x}, \vec{u} \f$
   */
  template <typename T>
  HOSTDEVICE T
  operator()(stc::Vector<T, InputDim> a_x) const
  {
    stc::Vector<T, Dim0> x0;
    for (int d = 0; d != Dim0; d++)
    {
      x0[d] = a_x[d];
    }
    stc::Vector<T, Dim1> x1;
    for (int d = 0; d != Dim1; d++)
    {
      x1[d] = a_x[d + Dim0];
    }
    return m_func0(x0) * m_func1(x1);
  }

private: Func0 m_func0;
  Func1 m_func1;
};

} // namespace exactSolution
} // namespace hoeb

#endif // EXACTSOLUTIONS_H_
