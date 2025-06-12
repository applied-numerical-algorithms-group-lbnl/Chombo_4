#ifndef _MULTIINDEXMATH_H_
#define _MULTIINDEXMATH_H_

#include <math.h>
#include "StcVector.hpp"

namespace hoeb
{
/**
 * @brief Returns the value of the mathematical constant pi.
 * @return The value of pi.
 */
HOSTDEVICE constexpr Real
Pi()
{
  return 3.1415926535897932384626433832795028841971693993751058209749445923078164062862089986280348253421170679;
}

/**
 * @brief Returns the value of the mathematical constant e.
 * @return The value of e.
 */
HOSTDEVICE constexpr Real
Exp()
{
  return 2.7182818284590452353602874713526624977572470936999595749669676277240766303535475945713821785251664274;
}

/**
 * @defgroup MultiIndexMath Multi-Index Math Functions
 * @brief A collection of functions for performing mathematical operations on vector-indices.
 * @{
 */

/**
 * @brief Computes a scalar to an integer power.
 * @tparam T The type of the input value.
 * @param a_x The input value.
 * @param a_p The power to which the input value is raised.
 * @return The result of \f$x^p\f$.
 */
template <typename T>
HOSTDEVICE constexpr T
ipow(const T a_x, const int a_p)
{
  T retval = 1;
  if (a_p > 0)
  {
    for (int iexp = 0; iexp < std::abs(a_p); iexp++)
    {
      retval *= a_x;
    }
  }
  else if (a_p < 0)
  {
    for (int iexp = 0; iexp < std::abs(a_p); iexp++)
    {
      retval /= a_x;
    }
  }
  return retval;
}

/**
 * @brief Computes the product of a scalar to integer vector power.
 * @tparam Dim The dimension of the input vector.
 * @param a_x The input value.
 * @param a_p The power to which the input value is raised.
 * @return The result of \f$\prod_{d=0}^{D} (x^{p_d})\f$.
 */
template <int Dim>
HOSTDEVICE constexpr Real
ipow(const Real a_x, const stc::IVec<Dim> a_p)
{
  Real retval = 1;
  for (int idir = 0; idir != Dim; idir++)
  {
    retval *= hoeb::ipow(a_x, a_p[idir]);
  }
  return retval;
}

/**
 * @brief Computes the product of a vector quantity to an integer
 * vector power.
 * @tparam Dim The dimension of the input vector.
 * @tparam T The type of the input vector.
 * @param a_x The input vector.
 * @param a_p The power to which the input vector is raised.
 * @return The result of \f$\prod_{d=0}^{D} (x_d^{p_d})\f$.
 */
template <int Dim, typename T>
HOSTDEVICE constexpr T
ipow(const stc::Vector<T, Dim> a_x, const stc::IVec<Dim> a_p)
{
  T retval = 1;
  for (int idir = 0; idir != Dim; idir++)
  {
    retval *= hoeb::ipow(a_x[idir], a_p[idir]);
  }
  return retval;
}

/**
 * @brief Calculates the factorial of an integer.
 * @param n The input integer.
 * @return The factorial \f$n!\f$.
 */
HOSTDEVICE
constexpr long int
factorial(const int n)
{
  assert(n >= 0); // "negative n"
  long int nfact = 1;
  for (int i = 2; i <= n; i++)
  {
    nfact *= i;
  };
  return nfact;
}

/**
 * @brief Calculates the factorials of a multinomial.
 * @tparam Dim The dimension of the input vector.
 * @param p The input vector.
 * @return The factorial \f$\prod__{d=0}^{D} {p_d}!\f$.
 */
template <int Dim>
HOSTDEVICE constexpr long int
factorial(const stc::IVec<Dim> p)
{
  long int pfact = 1;
  for (int idir = 0; idir < Dim; idir++)
  {
    pfact *= factorial(p[idir]);
  }
  return pfact;
}

/**
 * @brief Calculates the falling factorial \f$n! / (n-k)!\f$.
 * @param n The leading factorial integer upper value.
 * @param k The value at which the factorial is truncated.
 * @return The falling factorial.
 */
HOSTDEVICE
constexpr long int
fallingFactorial(const int n, const int k)
{
  assert(k >= 0);
  assert(n >= 0);
  long int val = 1;
  if (k == 0) return val;
  for (int i = 0; i != k; i++)
  {
    val *= (n - i);
  }
  return val;
}

/**
 * @brief Calculates the rising factorial \f$(n+k)! / n!\f$.
 * @param n The leading factorial integer lower value.
 * @param k The value as which the factorial is truncated.
 * @return The rising factorial.
 */
HOSTDEVICE
constexpr long int
risingFactorial(const int n, const int k)
{
  assert(k >= 0);
  assert(n >= 0);
  long int val = 1;
  if (k == 0) return val;
  for (int i = 0; i != k; i++)
  {
    val *= (n + i);
  }
  return val;
}

/**
 * @brief Calculates the binomial coefficient \f$n! / k! (n-k)!\f$.
 * @param n The number of options.
 * @param k The number of selections.
 * @return The binomial coefficient of "n choose k"
 */
HOSTDEVICE
constexpr long int
nCk(const int n, const int k)
{
  assert((n >= k) && (k >= 0) && (n >= 0)); // "out of range n"
  return fallingFactorial(n, k) / factorial(k);
}

/**
 * @brief Calculates the multi-binomial coefficient \f$\vec{p}! / \vec{k}! (\vec{n}-\vec{k})!\f$.
 * @tparam Dim The dimension of the input vectors, typically compiler deduced
 * @param p The number of vector options.
 * @param k The number of vector selections.
 * @return The multi-binomial coefficient "p choose k"
 */
template <int Dim>
HOSTDEVICE constexpr long int
pCk(const stc::IVec<Dim> p, const stc::IVec<Dim> k)
{
  long int pfact = 1;
  for (int idir = 0; idir < Dim; idir++)
  {
    pfact *= nCk(p[idir], k[idir]);
  }
  return pfact;
}

/**
 * @}
 */

/**
 * @name Vector Manipulation Functions
 * @{
 */

/**
 * @brief Reduces the dimension of a vector by removing the specified index.
 * @tparam Dim The dimension of the input vector, typically deduced
 * @tparam T The type of the input vector, typically deduced
 * @param a_vec The input vector.
 * @param a_index The index to remove.
 * @return The reduced vector.
 */
template <unsigned Dim, typename T>
HOSTDEVICE stc::Vector<T, Dim - 1>
contractDimension(const stc::Vector<T, Dim> a_vec, const int a_index)
{
  assert(a_index >= 0);
  assert(a_index < Dim);
  stc::Vector<T, Dim - 1> reducedVec;
  for (int i = 0, j = 0; i < Dim; ++i)
  {
    if (i != a_index)
    {
      reducedVec[j] = a_vec[i];
      j++;
    }
  }
  return reducedVec;
}

/**
 * @brief Reduces the dimension of a vector by removing the specified indices.
 * @tparam OutDim The dimension of the output vector.
 * @tparam InDim The dimension of the input vector, typically deduced
 * @tparam T The type of the input vector, typically deduced
 * @param a_vec The input vector.
 * @param a_index The indices to remove, marked by 0s. Non-zeros are are used to scale the input vector.
 * @return The reduced vector.
 */
template <unsigned OutDim, unsigned InDim, typename T, std::enable_if_t<(InDim >= OutDim), int> = 0>
HOSTDEVICE stc::Vector<T, OutDim>
contractDimensions(const stc::Vector<T, InDim> a_vec, const stc::Vector<int, InDim> a_index)
{
  stc::Vector<T, OutDim> reducedVec;
  for (int i = 0, j = 0; i < InDim; ++i)
  {
    if ((a_index[i] != 0) && (j < OutDim))
    {
      reducedVec[j] = a_vec[i] * a_index[i];
      j++;
    }
  }
  return reducedVec;
}

/**
 * @brief Increases the dimension of a vector by adding the specified index and value.
 * @tparam Dim The dimension of the input vector, typically deduced
 * @tparam T The type of the input vector, typically deduced
 * @param a_vec The input vector.
 * @param a_index The index to add.
 * @param a_newVal The value to add.
 * @return The expanded vector.
 */
template <unsigned Dim, typename T>
HOSTDEVICE stc::Vector<T, Dim + 1>
expandDimension(const stc::Vector<T, Dim> a_vec, const int a_index, const int a_newVal = 0)
{
  assert(a_index >= 0);
  assert(a_index < Dim + 1);
  stc::Vector<T, Dim + 1> growVec;
  growVec[a_index] = a_newVal;

  for (int i = 0; i < a_index; ++i)
  {
    growVec[i] = a_vec[i];
  }
  for (int i = a_index; i < Dim; ++i)
  {
    growVec[i + 1] = a_vec[i];
  }
  return growVec;
}

/**
 * @brief Increases the dimension of a vector by adding the specified indices and values.
 * @tparam OutDim The dimension of the output vector.
 * @tparam InDim The dimension of the input vector, typically deduced
 * @tparam T The type of the input vector, typically deduced
 * @param a_vec The input vector.
 * @param a_index The indices to add. WARNING must be in sequential order
 * @param a_newVal The values to add.
 * @return The expanded vector.
 */
template <unsigned OutDim, unsigned InDim, typename T, std::enable_if_t<(InDim <= OutDim), int> = 0>
HOSTDEVICE stc::Vector<T, OutDim>
expandDimensions(const stc::Vector<T, InDim> a_vec, const stc::Vector<int, OutDim - InDim> a_index,
                 const stc::Vector<T, OutDim - InDim> a_newVal)
{
  assert(stc::minElem(a_index) >= 0);
  assert(stc::maxElem(a_index) <= OutDim);
  //std::sort(a_index.begin(), a_index.end());

  stc::Vector<T, OutDim> growVec;

  for (int i = 0, j = 0, k = 0; i < OutDim; ++i)
  {
    if (a_index[j] == i)
    {
      growVec[i] = a_newVal[j];
      j++;
    }
    else
    {
      growVec[i] = a_vec[k];
      k++;
    }
  }
  return growVec;
}

/**
 * @}
 */

} // namespace hoeb
#endif
