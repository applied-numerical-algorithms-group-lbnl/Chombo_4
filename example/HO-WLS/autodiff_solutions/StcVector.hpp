/*
 *       __  __          __   ___      _      __ 
 *      / /_/ /__  ___  / /  /   \_________  / /__
 *     / __  / _ \/_==\/ _ \/ --</ __/ / __\/  __/
 *    /_/ /_/\__._\___/_//_/_.__/_/ /_/\___/_/\_\
 *    This file is part of HashBrick.
 *    Please refer to Copyright and License, in HashBrick's root directory.
 */

#ifndef _STCVECTOR_H_
#define _STCVECTOR_H_


/******************************************************************************/
/**
 * \file StcVector.hpp
 *
 * \brief Static vector class
 *
 *//*+*************************************************************************/

#include <iostream>
#include <utility>
#include <cassert>
#include <limits>
#include <cmath>
#include <type_traits>
#include <algorithm>
#include <bitset>
//#include <array>

#include "Parameters.hpp"
#include "Real.hpp"
#include "LinuxSupport.hpp"

/*
  Notes:
    - tag USE_GPU means enable GPU in addition to CPU
    - tag __CUDACC__ means the file is processed with nvcc.
    - tag __CUDA_ARCH__ is defined for device code.  This only works in
      functions with the __device__ qualifier.  It cannot be used in the
      class declaration.
*/

#undef HOSTDEVICE
#undef DEVICE
#ifdef __CUDACC__
#define HOSTDEVICE __host__ __device__
#define DEVICE __device__
#else
#define HOSTDEVICE
#define DEVICE
#endif

#ifdef __CUDACC__
#include "nppdefs.h"
#endif

/*
  Comments on Lambdas:
  A major goal of this class is constexpr evaluation of small-vector operations
  whenever possible.  Lambdas are not used for two reasons:
    - They cannot be constexpr until C++17
    - nvcc excessively complains about capture by reference.
  Explicit function objects are used instead

  The enable_if_t pointers need to be equated to a nullptr for constexpr use
*/

/*--------------------------------------------------------------------*
 * Spatial dimension
 *--------------------------------------------------------------------*/

#ifndef SPACEDIM
  #define SPACEDIM 3
#endif


namespace stc
{

/*--------------------------------------------------------------------*
 * Helper functions
 *--------------------------------------------------------------------*/

/*
  These are only applied to vector sizes (which ultimately derive from a
  templated parameter).  Some compilers do not yet provide a constexpr
  std::min and std::max.
*/

/// Compile-time minimum.
template <typename T>
HOSTDEVICE inline constexpr T
cmin(const T& a_x, const T& a_y) noexcept
{
  return (a_x < a_y) ? a_x : a_y;
}
/// Compile-time maximum
template <typename T>
HOSTDEVICE inline constexpr T
cmax(const T& a_x, const T& a_y) noexcept
{
  return (a_y < a_x) ? a_x : a_y;
}

/// constexpr log2
/** Note: T can be signed but if it is <= 0, this might recurse to the
 *  maximum template depth if right shift is arithmetic.  Of course, log is
 *  either infinite or complex for a_i <= 0 so avoid that!
 */
template <typename T>
HOSTDEVICE inline constexpr T
_clog2(const T a_i)
{
  static_assert(std::is_integral<T>::value, "clog2 on non-integer");
  return (a_i == (T)1) ? (T)0 : (T)1 + _clog2(a_i >> 1);
}

/// constexpr popcnt
/** Note: suitable for compile-time evaluation only.  See popcnt for the
 *  run-time version.  Will not work on a negative integer.
 */
template <typename T>
HOSTDEVICE inline constexpr T
_cpopcnt(const T a_i)
{
  CH_assert(a_i > 0);
  static_assert(std::is_integral<T>::value, "popcnt on non-integer");
  return (a_i == (T)0) ? (T)0 : (T)(a_i & (T)1) + _cpopcnt(a_i >> 1);
}

/*
  Portable popcnt (if it does not use popcnt instruction when available, get a
  decent compiler).  A test file is embedded here.  E.g., compile with
  g++ -march=native -O3 -S filename.cpp and look for popcntq in the assembly
*/
/*----------------------------------------------------------------------------*/
/*
#include <cstdio>
#include <bitset>
#include <limits>

#define HOSTDEVICE
*/
// <= 32 bit popcnt
template <typename T>
HOSTDEVICE inline auto static
popcntImpl(const T a_i, std::true_type)
{
  static_assert(sizeof(T) <= sizeof(unsigned),
                "Insufficient size of integer for type sent to 32-bit popcnt");
#ifdef __CUDA_ARCH__
  return __popc(static_cast<unsigned>(a_i));
#else
  return std::bitset<std::numeric_limits<unsigned>::digits>(a_i).count();
#endif
}
// > 32 bit popcnt
template <typename T>
HOSTDEVICE inline auto static
popcntImpl(const T a_i, std::false_type)
{
  static_assert(sizeof(T) <= sizeof(unsigned long long),
    "Insufficient size of integer for type sent to popcnt");
#ifdef __CUDA_ARCH__
  return __popcll(static_cast<unsigned long long>(a_i));
#else
  return std::bitset<std::numeric_limits<unsigned long long>::digits>(a_i)
    .count();
#endif
}
// Dispatch based on size
template <typename T>
HOSTDEVICE inline auto
popcnt(const T a_i)
{
  static_assert(std::is_integral<T>::value, "popcnt on non-integer");
  return popcntImpl(a_i, std::integral_constant<bool, (sizeof(T) <= 4)>{});
}
/*
int main()
{
  volatile unsigned u = 15;
  int i = popcnt(u);
  printf("%d\n", i);
}
*/
/*----------------------------------------------------------------------------*/

// Check for power of 2
template <typename T>
HOSTDEVICE inline constexpr bool
isPower2(const T a_i)
{
  static_assert(std::is_integral<T>::value, "isPower2 on non-integer");
  static_assert(sizeof(T) <= sizeof(unsigned long long),
    "Insufficient size of integer for type sent to isPower2");
  const unsigned long long u = a_i;
  return (a_i > 0) && !(u & (u-1));
}

//**From CH_System in Chombo

/*--------------------------------------------------------------------*/
/// Coarsen by a power of 2
/**
 *  \param[in]  a_i     Integer
 *  \param[in]  a_nref  Refinement ratio (must be power of 2)
 *  Cannot be constexpr because of popcnt
 *//*-----------------------------------------------------------------*/

template <typename T1, typename T2>
HOSTDEVICE inline T1
coarsen_sys(const T1& a_i, const T2& a_nref)
{
  static_assert(std::is_integral<T1>::value &&
                std::is_integral<T2>::value,
                "Coarsen can only be applied to integral vector types");
  CH_assert(isPower2(a_nref));
  return a_i >> popcnt(a_nref-1);
}

/*--------------------------------------------------------------------*/
/// Coarsen by a power of 2 given by 2^a_x
/**
 *  \param[in]  a_i     Integer
 *  \param[in]  a_x     Refinement ratio given by 2^a_x (a_x must be
 *                      >= 0)
 *//*-----------------------------------------------------------------*/

template <typename T1, typename T2>
HOSTDEVICE inline constexpr T1
coarsen_log2(const T1& a_i, const T2& a_x)
{
  static_assert(std::is_integral<T1>::value &&
                std::is_integral<T2>::value,
                "Coarsen can only be applied to integral vector types");
    return a_i >> a_x;
}

/*--------------------------------------------------------------------*
 * Preliminary definitions
 *--------------------------------------------------------------------*/

/*
  Normally one would use a std::array as the default implementation.  The
  reason for the custom implementation is for use with CUDA and because the
  size is changed to 32-bit by default
*/

#if true
using array_size_type = unsigned;
template <typename T, array_size_type N>
struct DefaultImpl
{
  // This is supposed to be an aggregate and use aggregate initialization so no
  // constructors are allowed
  HOSTDEVICE static constexpr array_size_type size() noexcept { return N; }
  HOSTDEVICE T* data() noexcept { return m_data; }
  HOSTDEVICE constexpr const T* data() const noexcept { return m_data; }
  HOSTDEVICE constexpr T& operator[](const array_size_type a_idx) noexcept
    { return m_data[a_idx]; }
  HOSTDEVICE constexpr const T& operator[](
    const array_size_type a_idx) const noexcept
    { return m_data[a_idx]; }
  HOSTDEVICE constexpr const T* begin() const noexcept { return m_data; }
  HOSTDEVICE T* begin() noexcept { return m_data; }
  HOSTDEVICE constexpr const T* end() const noexcept { return m_data + N; }
  HOSTDEVICE T* end() noexcept { return m_data + N; }
  T m_data[N];
};
template <typename T>
struct DefaultImpl<T, 0>
{
  HOSTDEVICE static constexpr array_size_type size() noexcept { return 0; }
  HOSTDEVICE T* data() noexcept { return static_cast<T*>(nullptr); }
  HOSTDEVICE constexpr const T* data() const noexcept
    { return static_cast<T*>(nullptr); }
  HOSTDEVICE constexpr T operator[](const array_size_type) noexcept
    { return T{}; }
  HOSTDEVICE constexpr const T operator[](
    const array_size_type a_idx) const noexcept
    { return T{}; }
  HOSTDEVICE constexpr const T* begin() const noexcept
    { return static_cast<T*>(nullptr); }
  HOSTDEVICE T* begin() noexcept { return static_cast<T*>(nullptr); }
  HOSTDEVICE constexpr const T* end() const noexcept
    { return static_cast<T*>(nullptr); }
  HOSTDEVICE T* end() noexcept { return static_cast<T*>(nullptr); }
};
#else
// Almost certainly size_t but to be formal
using array_size_type = typename std::array<int, 0>::size_type;
template <typename T, array_size_type N>
using DefaultImpl = std::array<T, N>;
#endif

// All classes that serve as reps for expressions must derive from this
class StcOp
{
};

//--Forward declarations

template <typename T>
class ScalarRep;

// Operator +
template <typename REP1, typename REP2,
          std::enable_if_t<
            std::is_base_of<stc::StcOp, std::decay_t<REP1>>::value
            &&
            std::is_base_of<stc::StcOp, std::decay_t<REP2>>::value
            &&
            std::decay_t<REP1>::is_vector::value
            &&
            std::decay_t<REP2>::is_vector::value
            &&
            !std::decay_t<REP1>::is_string::value
            &&
            !std::decay_t<REP2>::is_string::value, int> = 0
          >
HOSTDEVICE constexpr auto
operator+(REP1&& a_rep1, REP2&& a_rep2) noexcept;

template <typename REP1, typename T2,
          std::enable_if_t<
            std::is_base_of<stc::StcOp, std::decay_t<REP1>>::value
            &&
            !std::is_base_of<stc::StcOp, std::decay_t<T2>>::value
            &&
            std::decay_t<REP1>::is_vector::value
            &&
            !std::decay_t<REP1>::is_string::value
            &&
            !std::is_class<std::decay_t<T2>>::value, int> = 0
          >
HOSTDEVICE constexpr auto
operator+(REP1&& a_rep1, const T2& a_scalar) noexcept;

// Operator -
template <typename REP1, typename REP2,
          std::enable_if_t<
            std::is_base_of<stc::StcOp, std::decay_t<REP1>>::value
            &&
            std::is_base_of<stc::StcOp, std::decay_t<REP2>>::value
            &&
            std::decay_t<REP1>::is_vector::value
            &&
            std::decay_t<REP2>::is_vector::value, int> = 0
          >
HOSTDEVICE constexpr auto
operator-(REP1&& a_rep1, REP2&& a_rep2) noexcept;

template <typename REP1, typename T2,
          std::enable_if_t<
            std::is_base_of<stc::StcOp, std::decay_t<REP1>>::value
            &&
            !std::is_base_of<stc::StcOp, std::decay_t<T2>>::value
            &&
            std::decay_t<REP1>::is_vector::value
            &&
            !std::is_class<std::decay_t<T2>>::value, int> = 0
          >
HOSTDEVICE constexpr auto
operator-(REP1&& a_rep1, const T2& a_scalar) noexcept;

// Operator *
template <typename REP1, typename REP2,
          std::enable_if_t<
            std::is_base_of<stc::StcOp, std::decay_t<REP1>>::value
            &&
            std::is_base_of<stc::StcOp, std::decay_t<REP2>>::value
            &&
            std::decay_t<REP1>::is_vector::value
            &&
            std::decay_t<REP2>::is_vector::value, int> = 0
          >
HOSTDEVICE constexpr auto
operator*(REP1&& a_rep1, REP2&& a_rep2) noexcept;

template <typename REP1, typename T2,
          std::enable_if_t<
            std::is_base_of<stc::StcOp, std::decay_t<REP1>>::value
            &&
            !std::is_base_of<stc::StcOp, std::decay_t<T2>>::value
            &&
            std::decay_t<REP1>::is_vector::value
            &&
            !std::is_class<std::decay_t<T2>>::value, int> = 0
          >
HOSTDEVICE constexpr auto
operator*(REP1&& a_rep1, const T2& a_scalar) noexcept;

// Operator /
template <typename REP1, typename REP2,
          std::enable_if_t<
            std::is_base_of<stc::StcOp, std::decay_t<REP1>>::value
            &&
            std::is_base_of<stc::StcOp, std::decay_t<REP2>>::value
            &&
            std::decay_t<REP1>::is_vector::value
            &&
            std::decay_t<REP2>::is_vector::value, int> = 0
          >
HOSTDEVICE constexpr auto
operator/(REP1&& a_rep1, REP2&& a_rep2) noexcept;

template <typename REP1, typename T2,
          std::enable_if_t<
            std::is_base_of<stc::StcOp, std::decay_t<REP1>>::value
            &&
            !std::is_base_of<stc::StcOp, std::decay_t<T2>>::value
            &&
            std::decay_t<REP1>::is_vector::value
            &&
            !std::is_class<std::decay_t<T2>>::value, int> = 0
          >
HOSTDEVICE constexpr auto
operator/(REP1&& a_rep1, const T2& a_scalar) noexcept;

// Min
template <typename REP1, typename REP2,
          std::enable_if_t<
            std::is_base_of<stc::StcOp, std::decay_t<REP1>>::value
            &&
            std::is_base_of<stc::StcOp, std::decay_t<REP2>>::value
            &&
            std::decay_t<REP1>::is_vector::value
            &&
            std::decay_t<REP2>::is_vector::value, int> = 0
          >
HOSTDEVICE constexpr auto
min(REP1&& a_rep1, REP2&& a_rep2) noexcept;

// Min (right scalar)
template <typename REP1, typename T2,
          std::enable_if_t<
            std::is_base_of<stc::StcOp, std::decay_t<REP1>>::value
            &&
            !std::is_base_of<stc::StcOp, std::decay_t<T2>>::value
            &&
            std::decay_t<REP1>::is_vector::value
            &&
            !std::is_class<std::decay_t<T2>>::value, int> = 0
          >
HOSTDEVICE constexpr auto
min(REP1&& a_rep1, const T2& a_scalar) noexcept;

// Max
template <typename REP1, typename REP2,
          std::enable_if_t<
            std::is_base_of<stc::StcOp, std::decay_t<REP1>>::value
            &&
            std::is_base_of<stc::StcOp, std::decay_t<REP2>>::value
            &&
            std::decay_t<REP1>::is_vector::value
            &&
            std::decay_t<REP2>::is_vector::value, int> = 0
          >
HOSTDEVICE constexpr auto
max(REP1&& a_rep1, REP2&& a_rep2) noexcept;

// Max (right scalar)
template <typename REP1, typename T2,
          std::enable_if_t<
            std::is_base_of<stc::StcOp, std::decay_t<REP1>>::value
            &&
            !std::is_base_of<stc::StcOp, std::decay_t<T2>>::value
            &&
            std::decay_t<REP1>::is_vector::value
            &&
            !std::is_class<std::decay_t<T2>>::value, int> = 0
          >
HOSTDEVICE constexpr auto
max(REP1&& a_rep1, const T2& a_scalar) noexcept;

// Lexical less than
template <typename REP1, typename REP2,
          std::enable_if_t<
            std::decay_t<REP1>::is_vector::value
            &&
            std::decay_t<REP2>::is_vector::value, int> = 0
          >
HOSTDEVICE constexpr bool
lexLT(REP1&& a_rep1, REP2&& a_rep2) noexcept;

// Lexical greater than
template <typename REP1, typename REP2,
          std::enable_if_t<
            std::decay_t<REP1>::is_vector::value
            &&
            std::decay_t<REP2>::is_vector::value, int> = 0
          >
HOSTDEVICE constexpr bool
lexGT(REP1&& a_rep1, REP2&& a_rep2) noexcept;

// Index of the minimum element
template <typename REP,
          std::enable_if_t<std::decay_t<REP>::is_vector::value, int> = 0
          >
HOSTDEVICE constexpr array_size_type
indexMinElem(REP&& a_rep) noexcept;

// Value of the minimum element
template <typename REP,
          std::enable_if_t<std::decay_t<REP>::is_vector::value, int> = 0
          >
HOSTDEVICE constexpr auto
minElem(REP&& a_rep) noexcept;

// Index of the maximum element
template <typename REP,
          std::enable_if_t<std::decay_t<REP>::is_vector::value, int> = 0
          >
HOSTDEVICE constexpr array_size_type
indexMaxElem(REP&& a_rep) noexcept;

// Value of the maximum element
template <typename REP,
          std::enable_if_t<std::decay_t<REP>::is_vector::value, int> = 0
          >
HOSTDEVICE constexpr auto
maxElem(REP&& a_rep) noexcept;

// Index of the element with the largest most significant bit
template <typename REP,
          std::enable_if_t<std::decay_t<REP>::is_vector::value, int> = 0
          >
HOSTDEVICE constexpr array_size_type
indexMaxMSBElem(REP&& a_rep) noexcept;

// Coarsen (vector)
template <typename REP1, typename REP2,
          std::enable_if_t<
            std::decay_t<REP1>::is_vector::value
            &&
            std::decay_t<REP2>::is_vector::value, int> = 0
          >
HOSTDEVICE auto
coarsen(REP1&& a_rep1, REP2&& a_rep2) noexcept;

// Coarsen (scalar)
template <typename REP,
          std::enable_if_t<std::decay_t<REP>::is_vector::value, int> = 0
          >
HOSTDEVICE auto
coarsen(REP&& a_rep1, unsigned long long a_factor) noexcept;

// Dot product
template <typename REP1, typename REP2,
          std::enable_if_t<
            std::is_base_of<stc::StcOp, std::decay_t<REP1>>::value
            &&
            std::is_base_of<stc::StcOp, std::decay_t<REP2>>::value, int> = 0
          >
HOSTDEVICE constexpr auto
dot(REP1&& a_rep1, REP2&& a_rep2) noexcept;

// Absolute value  
template <typename REP,
          std::enable_if_t<std::decay_t<REP>::is_vector::value, int> = 0
          >
HOSTDEVICE constexpr auto
abs(REP&& a_rep) noexcept;

// Sum
template <typename REP,
          std::enable_if_t<std::decay_t<REP>::is_vector::value, int> = 0
          >
HOSTDEVICE constexpr auto
sum(REP&& a_rep) noexcept;

// Product
template <typename REP,
          std::enable_if_t<std::decay_t<REP>::is_vector::value, int> = 0
          >
HOSTDEVICE constexpr auto
product(REP&& a_rep) noexcept;

// 1-norm
template <typename REP,
          std::enable_if_t<std::decay_t<REP>::is_vector::value, int> = 0
          >
HOSTDEVICE constexpr auto
norm1(REP&& a_rep) noexcept;

// Magnitude
template <typename REP,
          std::enable_if_t<
            std::is_base_of<stc::StcOp, std::decay_t<REP>>::value, int> = 0
          >
HOSTDEVICE auto
mag(REP&& a_rep) noexcept;

/*--------------------------------------------------------------------*
 * This type is used to indicate transpose using syntax m^t
 *--------------------------------------------------------------------*/

/// Transpose type
struct TransposeType
{
};
/// Transpose marker
constexpr TransposeType t{};


/*******************************************************************************
 *
 * Makers for arrays
 *
 * A bit of wizardry to initialize the array.  The approach was taken from
 * the way make_index_sequence itself is implemented in gcc-6.1.0.
 *
 ******************************************************************************/

/*
 * The default parameter here generates the integer sequence
 */
template <typename T, array_size_type N, typename Op,
          typename ISeq = std::make_integer_sequence<array_size_type, N>>
struct Make_array;

/*
 * The sequence has type std::index_sequence<Is...> and now this specialization
 * matches
 */
template <typename T, array_size_type N, typename Op, array_size_type... Is>
struct Make_array<T, N, Op, std::integer_sequence<array_size_type, Is...>>
{
  HOSTDEVICE static constexpr DefaultImpl<T, N> with(Op&& a_op)
    {
      return { std::forward<Op>(a_op)(Is)... };  // op is applied to each
                                                 // integer in the sequence
    }
};

/*
 * Create a function so Op can be deduced
 */
template <typename T, array_size_type N, typename Op>
HOSTDEVICE constexpr auto
make_array(Op&& a_op)
{
  return Make_array<T, N, Op>::with(std::forward<Op>(a_op));
}


/*******************************************************************************
 *
 * Metaprogramming for unrolling loops
 *
 ******************************************************************************/

template <array_size_type I, typename F>
struct UnrollLoop
{
  /// Apply f(I-1) to each element
  HOSTDEVICE static constexpr void eval(F&& f)
    {
      UnrollLoop<I-1, F>::eval(std::forward<F>(f));
      f(I-1);
    }
  /// Apply f.operator()<I-1>(I-1) to each element
  HOSTDEVICE static constexpr void teval(F&& f)
    {
      UnrollLoop<I-1, F>::teval(std::forward<F>(f));
      f.template operator()<I-1>(I-1);
    }
  /// Apply f to each element in reverse order
  HOSTDEVICE static constexpr void reval(F&& f)
    {
      f(I-1);
      UnrollLoop<I-1, F>::reval(std::forward<F>(f));
    }
  /// Sum result of f
  HOSTDEVICE static constexpr auto sum(F&& f)
    {
      return UnrollLoop<I-1, F>::sum(std::forward<F>(f)) + f(I-1);
    }
  /// Product result of f
  HOSTDEVICE static constexpr auto product(F&& f)
    {
      return UnrollLoop<I-1, F>::product(std::forward<F>(f))*f(I-1);
    }
  /// And result of f
  HOSTDEVICE static constexpr bool _and(F&& f)
    {
      return UnrollLoop<I-1, F>::_and(std::forward<F>(f)) && f(I-1);
    }
  /// Or result of f
  HOSTDEVICE static constexpr bool _or(F&& f)
    {
      return UnrollLoop<I-1, F>::_or(std::forward<F>(f)) || f(I-1);
    }
  /// First true position of f
  HOSTDEVICE static constexpr array_size_type lotruepos(F&& f)
    {
      array_size_type subpos =
        UnrollLoop<I-1, F>::lotruepos(std::forward<F>(f));
      return subpos + (subpos == 0)*f(I-1)*I;
    }
  /// Last true position of f
  HOSTDEVICE static constexpr array_size_type hitruepos(F&& f)
    {
      array_size_type subpos =
        UnrollLoop<I-1, F>::hitruepos(std::forward<F>(f));
      return f(I-1) ? I : subpos;
    }
  /// Select a position from f
  HOSTDEVICE static constexpr array_size_type selectpos(F&& f)
    {
      return f(I-1, UnrollLoop<I-1, F>::selectpos(std::forward<F>(f)));
    }
};

template <typename F>
struct UnrollLoop<1, F>
{
  HOSTDEVICE static constexpr void eval(F&& f)
    {
      f(0);
    }
  HOSTDEVICE static constexpr void teval(F&& f)
    {
      f.template operator()<0>(0);
    }
  HOSTDEVICE static constexpr void reval(F&& f)
    {
      f(0);
    }
  HOSTDEVICE static constexpr auto sum(F&& f)
    {
      return f(0);
    }
  HOSTDEVICE static constexpr auto product(F&& f)
    {
      return f(0);
    }
  HOSTDEVICE static constexpr bool _and(F&& f)
    {
      return f(0);
    }
  HOSTDEVICE static constexpr bool _or(F&& f)
    {
      return f(0);
    }
  HOSTDEVICE static constexpr array_size_type lotruepos(F&& f)
    {
      return f(0);
    }
  HOSTDEVICE static constexpr array_size_type hitruepos(F&& f)
    {
      return f(0);
    }
  HOSTDEVICE static constexpr array_size_type selectpos(F&& f)
    {
      return 0;
    }
};

// In case N = 0, don't do anything
template <typename F>
struct UnrollLoop<0, F>
{
  HOSTDEVICE static constexpr auto eval(F&& f)
    { }
  HOSTDEVICE static constexpr auto teval(F&& f)
    { }
  HOSTDEVICE static constexpr auto reval(F&& f)
    { }
  HOSTDEVICE static constexpr auto sum(F&& f)
    {
      return static_cast<decltype(f(0))>(0);
    }
  HOSTDEVICE static constexpr auto product(F&& f)
    {
      return static_cast<decltype(f(0))>(1);
    }
  HOSTDEVICE static constexpr bool _and(F&& f)
    {
      return true;
    }
  HOSTDEVICE static constexpr bool _or(F&& f)
    {
      return false;
    }
  HOSTDEVICE static constexpr array_size_type lotruepos(F&& f)
    {
      return 0;
    }
  HOSTDEVICE static constexpr array_size_type hitruepos(F&& f)
    {
      return 0;
    }
  HOSTDEVICE static constexpr array_size_type selectpos(F&& f)
    {
      return 0;
    }
};

/// Apply function object 'f(I-1)' to each element of the vector
template <array_size_type I, typename F>
HOSTDEVICE constexpr void
forEachElement(F&& f)
{
  UnrollLoop<I, F>::eval(std::forward<F>(f));
}

/// Apply function object 'f.operator()<I-1>(I-1)' to each element of the vector
template <array_size_type I, typename F>
HOSTDEVICE constexpr void
forEachTElement(F&& f)
{
  UnrollLoop<I, F>::teval(std::forward<F>(f));
}

/// Apply function object 'f' to each element of the vector in reverse order
template <array_size_type I, typename F>
HOSTDEVICE constexpr void
forEachReverseElement(F&& f)
{
  UnrollLoop<I, F>::reval(std::forward<F>(f));
}

/// Apply function object 'f' to each element of the vector and sum results
template <array_size_type I, typename F>
HOSTDEVICE constexpr auto
sumEachElement(F&& f)
{
  return UnrollLoop<I, F>::sum(std::forward<F>(f));
}

/// Apply function object 'f' to each element of the vector and multiply results
template <array_size_type I, typename F>
HOSTDEVICE constexpr auto
prodEachElement(F&& f)
{
  return UnrollLoop<I, F>::product(std::forward<F>(f));
}

/// Apply function object 'f' to each element of the vector and 'and' results
template <array_size_type I, typename F>
HOSTDEVICE constexpr bool
andEachElement(F&& f)
{
  return UnrollLoop<I, F>::_and(std::forward<F>(f));
}

/// Apply function object 'f' to each element of the vector and 'or' results
template <array_size_type I, typename F>
HOSTDEVICE constexpr bool
orEachElement(F&& f)
{
  return UnrollLoop<I, F>::_or(std::forward<F>(f));
}

/// Report lowest element index where function object 'f' evaluates to true.
/** \return                 position or
 *                          std::numeric_limits<array_size_type>::max() if
 *                          all elements evaluate to false.
 *  note: for all false, lotruepos returns 0.  -1 expects wrapping and requires
 *        array_size_type to be unsigned (or alternatively choose a specific
 *        unsigned return type in hitruepos).
 */
template <array_size_type I, typename F>
HOSTDEVICE constexpr array_size_type
lowestTrueElementIndex(F&& f)
{
  static_assert(std::is_unsigned<stc::array_size_type>::value,
                "array_size_type must be unsigned");
  return UnrollLoop<I, F>::lotruepos(std::forward<F>(f)) - 1;
}

/// Report highest element index where function object 'f' evaluates to true.
/** \return                 position or
 *                          std::numeric_limits<array_size_type>::max() if
 *                          all elements evaluate to false.
 *  note: for all false, hitruepos returns 0.  -1 expects wrapping and requires
 *        array_size_type to be unsigned (or alternatively choose a specific
 *        unsigned return type in hitruepos).
 */
template <array_size_type I, typename F>
HOSTDEVICE constexpr array_size_type
highestTrueElementIndex(F&& f)
{
  static_assert(std::is_unsigned<stc::array_size_type>::value,
                "array_size_type must be unsigned");
  return UnrollLoop<I, F>::hitruepos(std::forward<F>(f)) - 1;
}

/// Select a position from sequential recursive evaluations of f(a, b)
/** \return                 position always >= 0
 */
template <array_size_type I, typename F>
HOSTDEVICE constexpr array_size_type
selectElementIndex(F&& f)
{
  return UnrollLoop<I, F>::selectpos(std::forward<F>(f));
}


/*******************************************************************************
 *
 * Class VectorImpl
 *
 * Usually built around a std::array component (or similar) but aliases are also
 * supported.  These static arrays can be any size and of any type.  Expressions
 * are also provided.
 *
 ******************************************************************************/

//--Alias implementations.  These take the same parameters as std::array and
//--discard irrelevant.

template <typename T, array_size_type>
using VectorAliasImpl = T*;

template <typename T, array_size_type>
using VectorConstAliasImpl = const T*;


/*============================================================================*/
/// Static vector
/** Normally, this class is not used directly.  There are numerous typedefs
 *  that create more accessible objects.
 *  \tparam T           Type of vector
 *  \tparam N           Size of vector
 *  \tparam Impl        Underlying implementation.  Default is DefaultImpl
 *                      (std::array or similar) but aliases are also possible
 *                      (use a type from above).
 *//*=========================================================================*/

template <typename T,
          array_size_type N,
          template<typename, array_size_type> class Impl = DefaultImpl>
class VectorImpl : StcOp
{

/*--------------------------------------------------------------------*
 * Types
 *--------------------------------------------------------------------*/

public:

  using value_type = T;
  using iterator = value_type*;
  using const_iterator = const value_type*;
  using is_scalar = std::false_type;
  using is_vector = std::true_type;
  using is_rowvec = std::false_type;
  using is_matrix = std::false_type;
  using is_string = std::false_type;
  using impl_type = Impl<T, N>;

/*--------------------------------------------------------------------*
 * Constructors and destructors
 *--------------------------------------------------------------------*/

  /// Default construct default initializes (uninitialized data)
  HOSTDEVICE VectorImpl()
    { }

  /// Construct from an implementation of same type and rank
  /** Useful when Impl is stored in shared memory
   */
  template <typename _Impl = Impl<T, N>,
            std::enable_if_t<
              std::is_same<_Impl, DefaultImpl<T, N>>::value, int> = 0
            >
  HOSTDEVICE constexpr VectorImpl(const DefaultImpl<T, N>& a_impl)
    :
    m_data(a_impl)
    { }

  /// Function object helper for array maker
  /** Returns elements from object Tarray by indexing with operator[].  Tarray
   *  has size Tsize.  If Tsize is <= a_idx, zero initialization is used
   *  instead.
   *  \tparam Tarray      Object to index with operator[]
   *  \tparam Tsize       Number of elements we can index in Tarray
   *  Function object is used instead of lambda to get constexpr.
   */
  template <typename Tarray, array_size_type Tsize>
  struct ArrayCopy_E
  {
    HOSTDEVICE constexpr ArrayCopy_E(const Tarray& a_array) noexcept
      :
      m_array(a_array)
      { }
    HOSTDEVICE constexpr T operator()(
      const array_size_type a_idx) const noexcept
      {
        return (a_idx < Tsize) ? m_array[a_idx] : (T)0;
      }
    const Tarray& m_array;
  };

  /// Construct from an implementation of same type and different rank
  /** Useful when Impl is stored in shared memory
   */
  template <array_size_type M, typename _Impl = Impl<T, N>,
            std::enable_if_t<
              std::is_same<_Impl, DefaultImpl<T, N>>::value, int> = 0
            >
  HOSTDEVICE constexpr VectorImpl(const DefaultImpl<T, M>& a_impl)
    :
    m_data(make_array<T, N>(ArrayCopy_E<DefaultImpl<T, M>, M>(a_impl)))
    { }

  /// Construct with a constant T
  struct Uniform_E
  {
    HOSTDEVICE constexpr Uniform_E(T a_x) noexcept
      :
      m_x(a_x)
      { }
    HOSTDEVICE constexpr T operator()(const array_size_type a_idx)
      const noexcept
      {
        return m_x;
      }
    T m_x;
  };
  // This cannot be used in implicit conversions
  HOSTDEVICE constexpr explicit VectorImpl(const T a_x)
    :
    m_data(make_array<T, N>(Uniform_E(a_x)))
    { }

  /// Constructor that only builds from rvalues.
  /** Use with make_vector helper functions
   */
  HOSTDEVICE constexpr VectorImpl(DefaultImpl<T, N>&& a_data) noexcept
    :
    m_data(std::move(a_data))
    { }

  // Builder for construction from arguments
  template <array_size_type M>
  struct ArrayArgs_E
  {
    template <typename... Args>
    HOSTDEVICE constexpr ArrayArgs_E(const Args&... a_args) noexcept
      :
      m_args({ (T)a_args... })
      { }
    HOSTDEVICE constexpr T operator()(const array_size_type a_idx)
      const noexcept
      {
        return (a_idx < M) ? m_args[a_idx] : (T)0;
      }
    const DefaultImpl<T, M> m_args;
  };

  /// Construct from arguments
  /** Generally, this expects N arguments of some types that are convertible to
   *  T, but it tolerates more or less (fills remainder with T(0).  This insists
   *  on at least two args to avoid conflict with the constructor that takes a
   *  single arg (constant T above).  In practice, we should check that all Args
   *  are convertible to T but that is more difficult without std::conjunction
   *  provided by c++17.  For now, only the first two arguments are checked.
   */
  template <typename Arg0, typename Arg1, typename... Args,
            std::enable_if_t<
              std::is_convertible<Arg0, T>::value
              &&
              std::is_convertible<Arg1, T>::value
              &&
              std::is_same<Impl<T, N>, DefaultImpl<T, N>>::value, int> = 0
            >
  HOSTDEVICE constexpr VectorImpl(const Arg0& a_arg0,
                                  const Arg1& a_arg1,
                                  const Args&... a_args)
    :
    m_data(make_array<T, N>(ArrayArgs_E<2 + sizeof...(Args)>(
                              a_arg0, a_arg1, a_args...)))
    { }

//--Construct: from pointer (use at your own risk)

private:  // For targets of dispatch

  // If destination is constant alias
  HOSTDEVICE explicit VectorImpl(const T* a_ptr, std::false_type) noexcept
    :
    m_data(a_ptr)
    { }

  // If destination is alias
  HOSTDEVICE explicit VectorImpl(T* a_ptr, std::false_type) noexcept
    :
    m_data(a_ptr)
    { }

  // If destination has storage
  /* This does not make an alias.  Data is taken from the array at a_ptr.
   */
  HOSTDEVICE explicit VectorImpl(const T* a_ptr, std::true_type) noexcept
    :
    m_data(make_array<T, N>(
             [a_ptr]
             (const array_size_type a_idx)
             { return a_ptr[a_idx]; }))
    { }

public:

  /// Dispatch for construction from const pointer (use at your own risk)
  HOSTDEVICE constexpr VectorImpl(const T* a_ptr) noexcept
    :
    VectorImpl(a_ptr, std::is_same<Impl<T, N>, DefaultImpl<T, N>>{})
    { }

  /// Dispatch for construction from pointer (use at your own risk)
  HOSTDEVICE constexpr VectorImpl(T* a_ptr) noexcept
    :
    VectorImpl(a_ptr, std::is_same<Impl<T, N>, DefaultImpl<T, N>>{})
    { }

//--Construct: from initializer list

  // // Helper function object for initializer_list constructor.  Used instead of
  // // lambda to allow constexpr.
  // struct IL_E
  // {
  //   HOSTDEVICE constexpr IL_E(const std::initializer_list<T> a_il) noexcept
  //     :
  //     m_il(a_il)
  //     { }
  //   HOSTDEVICE constexpr T operator()(const array_size_type a_idx)
  //     const noexcept
  //     {
  //       return (a_idx < m_il.size()) ? m_il.begin()[a_idx] : (T)0;
  //     }
  //   const std::initializer_list<T> m_il;
  // };

  // /// Constructor taking initializer... and forwarded to Make_array
  // HOSTDEVICE constexpr VectorImpl(std::initializer_list<T> a_il) noexcept
  //   :
  //   m_data(make_array<T, N>(IL_E(a_il)))
  //   { }

  // Helper function object for reversing initializer_list.  Used instead of
  // lambda to allow constexpr.
  struct ILReverse_E
  {
    HOSTDEVICE constexpr ILReverse_E(
      Impl<T, N>&                    a_array,
      const std::initializer_list<T> a_il) noexcept
      :
      m_array(a_array),
      m_il(a_il)
      { }
    HOSTDEVICE constexpr void operator()(const array_size_type a_idx)
      const noexcept
      {
        m_array[a_idx] = m_il.begin()[N - 1 - a_idx];
      }
    Impl<T, N>& m_array;
    const std::initializer_list<T> m_il;
  };

  /// Define that reverses initializer list
  /** Be careful with this because all expressions in the initializer list
   *  are evaluated in order.
   */
  HOSTDEVICE VectorImpl& reverse(std::initializer_list<T> a_il) noexcept
    {
      assert(N <= a_il.size());
      ILReverse_E f(m_data, a_il);
      forEachElement<N>(f);
      return *this;
    }

  /// Copy constructor
  VectorImpl(const VectorImpl&) = default;
  /// Move constructor
  VectorImpl(VectorImpl&&) = default;
  /// Assignment constructor
  VectorImpl& operator=(const VectorImpl&) = default;
  /// Move assignment constructor
  VectorImpl& operator=(VectorImpl&&) = default;
  /// Destructor
  ~VectorImpl() = default;

//--Construct from smaller vector and additional terms

  /// Function object helper for array maker
  template <array_size_type N1, typename Tarray, typename... OtherVals>
  struct ArrayVecAndVals_E
  {
    static constexpr array_size_type N2 = sizeof...(OtherVals);
    HOSTDEVICE constexpr ArrayVecAndVals_E(
      const Tarray&       a_array,
      const OtherVals&... a_otherVals) noexcept
      :
      m_array(a_array),
      m_otherVals({ a_otherVals... })
      { }
    HOSTDEVICE constexpr T operator()(
      const array_size_type a_idx) noexcept
      {
        return (a_idx < N1) ? m_array[a_idx] : m_otherVals[a_idx - N1];
      }
    const Tarray& m_array;
    DefaultImpl<T, N2> m_otherVals;
  };

  /// Construct a large vector from a small one plus additional terms
  /** Not recommended for large objects in runtime operations as T is stored
   *  in the builder.  On the other hand, this supports constexpr construction.
   *  For POD, this is recommended.
   */
  template <array_size_type M,
            template<typename, array_size_type> class Impl2,
            typename... OtherVals,
            std::enable_if_t<
              // Impl2 has smaller size and there are other values
              (M < N) && (sizeof...(OtherVals) > 0) &&
              (M + sizeof...(OtherVals) >= N) &&
              // Impl is not an alias
              std::is_same<Impl<T, N>, DefaultImpl<T, N>>::value, int> = 0
            >
  HOSTDEVICE constexpr VectorImpl(
    const VectorImpl<T, M, Impl2>& a_smallVec,
    const OtherVals&...            a_otherVals) noexcept
    :
    m_data(make_array<T, N>(ArrayVecAndVals_E<M, Impl2<T, M>, OtherVals...>(
                              a_smallVec.m_data,
                              a_otherVals...)))
    { }

  /// Construct from a vector-like type and other additional terms
  /** For this constructor, Vec is anything with operator[] and it's size is
   *  assumed to be at least N-sizeof...(OtherVals).
   *  Not recommended for large objects in runtime operations as T is stored
   *  in the builder.  On the other hand, this supports constexpr construction.
   *  For POD, this is recommended.
   */
  //**FIXME suggests that need for a general construction from something
  //**      using operator[] but without other vals.
  template <typename Vec,
            typename... OtherVals,
            std::enable_if_t<
              // Must be object with operator[]
              (std::is_pointer<Vec>::value ||
               std::is_array<Vec>::value ||
               std::is_class<Vec>::value)
              &&
              // But not part of Stc::Vector
              !std::is_base_of<StcOp, std::decay_t<Vec>>::value
              &&
              // Need at least one other value
              (sizeof...(OtherVals) > 0)
              &&
              // Avoid dispatch from constructors taking a pointer
              (std::conjunction<std::bool_constant<
                 std::is_same<T, OtherVals>::value>...>::value),
              int> = 0
            >
  HOSTDEVICE constexpr VectorImpl(
    const Vec&          a_smallVec,
    const OtherVals&... a_otherVals) noexcept
    :
    m_data(make_array<T, N>(
             ArrayVecAndVals_E<N - sizeof...(OtherVals), Vec, OtherVals...>(
               a_smallVec,
               a_otherVals...)))
    { }

//--Construct from leading terms and a smaller vector.  Only 1 and 2 leading
//--terms are supported due to the fact the template argument deduction does
//--not work for a leading parameter pack.

  /// Function object helper for array maker
  template <array_size_type M, typename Tarray>
  struct ArrayValsAndVec_E
  {
    template <array_size_type M1 = M, std::enable_if_t<(M1 == 1), int> = 0>
    HOSTDEVICE constexpr ArrayValsAndVec_E(
      const T&      a_val0,
      const Tarray& a_array) noexcept
      :
      m_leadingVals({ a_val0 }),
      m_array(a_array)
      { }
    template <array_size_type M1 = M, std::enable_if_t<(M1 == 2), int> = 0>
    HOSTDEVICE constexpr ArrayValsAndVec_E(
      const T&      a_val0,
      const T&      a_val1,
      const Tarray& a_array) noexcept
      :
      m_leadingVals({ a_val0, a_val1 }),
      m_array(a_array)
      { }
    template <array_size_type M1 = M, std::enable_if_t<(M1 == 3), int> = 0>
    HOSTDEVICE constexpr ArrayValsAndVec_E(
      const T&      a_val0,
      const T&      a_val1,
      const T&      a_val2,
      const Tarray& a_array) noexcept
      :
      m_leadingVals({ a_val0, a_val1, a_val2 }),
      m_array(a_array)
      { }
    HOSTDEVICE constexpr T operator()(
      const array_size_type a_idx) const noexcept
      {
        return (a_idx < M) ? m_leadingVals[a_idx] : m_array[a_idx - M];
      }
    DefaultImpl<T, M> m_leadingVals;
    const Tarray& m_array;
  };

  /// Construct a large vector from one leading term and a smaller vector
  /** Not recommended for large objects in runtime operations as T is stored
   *  in the builder.  On the other hand, this supports constexpr construction.
   *  For POD, this is recommended.
   */
  template <array_size_type M,
            template<typename, array_size_type> class Impl2,
            std::enable_if_t<
              // Impl2 has sufficient size
              M + 1 >= N &&
              // Impl is not an alias
              std::is_same<Impl<T, N>, DefaultImpl<T, N>>::value, int> = 0
            >
  HOSTDEVICE constexpr VectorImpl(
    const T&                       a_val0,
    const VectorImpl<T, M, Impl2>& a_smallVec) noexcept
    :
    m_data(make_array<T, N>(ArrayValsAndVec_E<1, Impl2<T, M>>(
                              a_val0, a_smallVec.m_data)))
    { }

  /// Construct a large vector from two leading terms and a smaller vector
  /** Not recommended for large objects in runtime operations as T is stored
   *  in the builder.  On the other hand, this supports constexpr construction.
   *  For POD, this is recommended.
   */
  template <array_size_type M,
            template<typename, array_size_type> class Impl2,
            std::enable_if_t<
              // Impl2 has sufficient size
              M + 2 >= N &&
              // Impl is not an alias
              std::is_same<Impl<T, N>, DefaultImpl<T, N>>::value, int> = 0
            >
  HOSTDEVICE constexpr VectorImpl(
    const T&                       a_val0,
    const T&                       a_val1,
    const VectorImpl<T, M, Impl2>& a_smallVec) noexcept
    :
    m_data(make_array<T, N>(ArrayValsAndVec_E<2, Impl2<T, M>>(
                              a_val0, a_val1, a_smallVec.m_data)))
    { }

  /// Construct a large vector from three leading terms and a smaller vector
  /** Not recommended for large objects in runtime operations as T is stored
   *  in the builder.  On the other hand, this supports constexpr construction.
   *  For POD, this is recommended.
   */
  template <array_size_type M,
            template<typename, array_size_type> class Impl2,
            std::enable_if_t<
              // Impl2 has sufficient size
              M + 3 >= N &&
              // Impl is not an alias
              std::is_same<Impl<T, N>, DefaultImpl<T, N>>::value, int> = 0
            >
  HOSTDEVICE constexpr VectorImpl(
    const T&                       a_val0,
    const T&                       a_val1,
    const T&                       a_val2,
    const VectorImpl<T, M, Impl2>& a_smallVec) noexcept
    :
    m_data(make_array<T, N>(ArrayValsAndVec_E<3, Impl2<T, M>>(
                              a_val0, a_val1, a_val2, a_smallVec.m_data)))
    { }

//--Alias (via constructor): if source is vector and destination is alias

  /// Alias a const vector
  //  Note: In GCC, you do not need S and can use T everywhere to force the
  //        same type.  However, clang++-8 and others insist on having the type
  //        for the std::enable_if_t.  And then we need the static assert.
  template <typename S,
            array_size_type M,
            template<typename, array_size_type> class Impl2,
            std::enable_if_t<
              (M >= N)
              &&
              std::is_same<Impl<S, M>, VectorConstAliasImpl<S, M>>::value,
              int> = 0
            >
  HOSTDEVICE constexpr VectorImpl(
    const VectorImpl<S, M, Impl2>& a_vec,
    const array_size_type          a_idxBeg = 0) noexcept
    :
    m_data(a_vec.data() + a_idxBeg)
    {
      static_assert(std::is_same<S, T>::value, "Alias must have same type");
      assert(a_idxBeg + N <= M);
    }

  // Alias a non-const vector
  template <typename S,
            array_size_type M,
            template<typename, array_size_type> class Impl2,
            std::enable_if_t<
              (M >= N)
              &&
              std::is_same<Impl<S, M>, VectorAliasImpl<S, M>>::value, int> = 0
            >
  HOSTDEVICE VectorImpl(VectorImpl<S, M, Impl2>& a_vec,
    const array_size_type                        a_idxBeg = 0) noexcept
    :
    m_data(a_vec.data() + a_idxBeg)
    {
      static_assert(std::is_same<S, T>::value, "Alias must have same type");
      assert(a_idxBeg + N <= M);
    }

//--Copy: offset copy by a_idxBeg from source vector or alias, requires M >= N

  /// Helper function object for array maker
  /** Returns elements from object Tarray by indexing with operator[] but
   *  starting at offset m_idxBeg.  Tarray has size Tsize.  If Tsize is <=
   *  a_idx, zero initialization is used instead.
   *  \tparam Tarray      Object to index with operator[]
   *  \tparam Tsize       Number of elements we can index in Tarray
   *  Function object is used instead of lambda to get constexpr.
   */
  template <typename Tarray, array_size_type Tsize>
  struct ArrayCopyAtOffset_E
  {
    HOSTDEVICE constexpr ArrayCopyAtOffset_E(
      const Tarray&         a_array,
      const array_size_type a_idxBeg) noexcept
      :
      m_idxBeg(a_idxBeg),
      m_array(a_array)
      { }
    HOSTDEVICE constexpr T operator()(
      const array_size_type a_idx) const noexcept
      {
        return (a_idx + m_idxBeg < Tsize) ? m_array[a_idx + m_idxBeg] : (T)0;
      }
    const array_size_type m_idxBeg;
    const Tarray& m_array;
  };

  // Offset copy
  template <typename S,
            array_size_type M,
            template<typename, array_size_type> class Impl2,
            std::enable_if_t<
              (M >= N)
              &&
              std::is_same<Impl<S, M>, DefaultImpl<S, M>>::value, int> = 0
            >
  HOSTDEVICE constexpr VectorImpl(
    const VectorImpl<S, M, Impl2>& a_vec,
    const array_size_type          a_idxBeg) noexcept
    :
    m_data(make_array<T, N>(ArrayCopyAtOffset_E<Impl2<S, M>, M>(a_vec.m_data,
                                                                a_idxBeg)))
    {
      static_assert(std::is_same<S, T>::value, "Source must have same type");
    }

//--Copy: if source is alias and destination is vector

  // For a const alias
  template <typename S,
            std::enable_if_t<
              std::is_same<Impl<S, N>, DefaultImpl<S, N>>::value, int> = 0
            >
  HOSTDEVICE constexpr VectorImpl(
    const VectorImpl<S, N, VectorConstAliasImpl>& a_alias) noexcept
    :
    m_data(make_array<T, N>(ArrayCopy_E<VectorConstAliasImpl<S, N>, N>(
                              a_alias.m_data)))
    {
      static_assert(std::is_same<S, T>::value, "Vector alias must have same "
                    "type");
    }

  // For a non-const alias
  template <typename S,
            std::enable_if_t<
              std::is_same<Impl<S, N>, DefaultImpl<S, N>>::value, int> = 0
            >
  HOSTDEVICE constexpr VectorImpl(
    const VectorImpl<S, N, VectorAliasImpl>& a_alias) noexcept
    :
    m_data(make_array<T, N>(ArrayCopy_E<VectorAliasImpl<S, N>, N>(
                              a_alias.m_data)))
    {
      static_assert(std::is_same<S, T>::value, "Vector alias must have same "
                    "type");
    }

//--Construct from an expression template given by Rep

  template <typename REP,
            std::enable_if_t<
              // These first two make sure the destination is not an alias
              !std::is_same<Impl<typename std::decay_t<REP>::value_type, N>,
                            VectorAliasImpl<
                              typename std::decay_t<REP>
                              ::value_type, N>>::value
              &&
              !std::is_same<Impl<typename std::decay_t<REP>::value_type, N>,
                            VectorConstAliasImpl<
                              typename std::decay_t<REP>
                              ::value_type, N>>::value
              &&
              // These make sure the source is not a vector of the same size
              // (where we prefer to use other constructors... for sanity)
              !std::is_base_of<VectorImpl<T, N, DefaultImpl>,
                               std::decay_t<REP>>::value
              &&
              !std::is_base_of<VectorImpl<T, N, VectorConstAliasImpl>,
                               std::decay_t<REP>>::value
              &&
              !std::is_base_of<VectorImpl<T, N, VectorAliasImpl>,
                               std::decay_t<REP>>::value
              &&
              std::is_base_of<StcOp, std::decay_t<REP>>::value, int> = 0
            >
  HOSTDEVICE constexpr VectorImpl(REP&& a_rhs) noexcept
    :
    m_data(make_array<T, N>(
             ArrayCopy_E<
               std::decay_t<REP>,
               cmin(VectorImpl::size(), std::decay_t<REP>::size())
             >(a_rhs)))
    { }

/*--------------------------------------------------------------------*
 * Member functions
 *--------------------------------------------------------------------*/

  /// Size of the vector
  HOSTDEVICE static constexpr array_size_type size() noexcept
    {
      return N;
    }

//--If this is an alias, change the pointer

  /// Alias another const vector
  //  Note: In GCC, you do not need S and can use T everywhere to force the
  //        same type.  However, clang++-8 and others insist on having the type
  //        for the std::enable_if_t.  And then we need the static assert.
  template <typename S,
            array_size_type M,
            template<typename, array_size_type> class Impl2,
            std::enable_if_t<
              std::is_same<Impl<S, M>, VectorConstAliasImpl<S, M>>::value,
              int> = 0
            >
  HOSTDEVICE VectorImpl& alias(
    const VectorImpl<S, M, Impl2>& a_vec,
    const array_size_type          a_idxBeg = 0) noexcept
    {
      static_assert(std::is_same<S, T>::value, "Alias must have same type");
      m_data = &a_vec[a_idxBeg];
      return *this;
    }

  /// Alias another non-const vector
  template <typename S,
            array_size_type M,
            template<typename, array_size_type> class Impl2,
            std::enable_if_t<
              std::is_same<Impl<S, M>, VectorAliasImpl<S, M>>::value, int> = 0
            >
  HOSTDEVICE VectorImpl& alias(VectorImpl<S, M, Impl2>& a_vec,
                               const array_size_type    a_idxBeg = 0) noexcept
    {
      static_assert(std::is_same<S, T>::value, "Alias must have same type");
      m_data = &a_vec[a_idxBeg];
      return *this;
    }

  /// Alias a pointer
  HOSTDEVICE VectorImpl& alias(T*                    a_ptr,
                               const array_size_type a_idxBeg = 0) noexcept
    {
      m_data = a_ptr + a_idxBeg;
      return *this;
    }

  /// Alias a const pointer
  HOSTDEVICE VectorImpl& alias(const T*              a_ptr,
                               const array_size_type a_idxBeg = 0) noexcept
    {
      m_data = a_ptr + a_idxBeg;
      return *this;
    }

//--Operations that evaluate expressions

  template <typename REP2>
  struct ArrayEq_E
  {
    HOSTDEVICE constexpr ArrayEq_E(Impl<T, N>& a_array,
                                   const REP2& a_rhs) noexcept
      :
      m_array(a_array),
      m_rhs(a_rhs)
      { }
    HOSTDEVICE constexpr void operator()(
      const array_size_type a_idx) const noexcept
      {
        m_array[a_idx] = m_rhs[a_idx];
      }
    Impl<T, N>& m_array;
    const REP2&  m_rhs;
  };

  /// Evaluate an expression template given by Rep (=)
  //  Was updated to also support assignment between different representations
  /** This is never called if RHS is any vector implementation
   *  \tparam     REP     Representation of expression
   *  \param[in]  a_rhs   Right-hand-side expression
   */
  template <typename REP,
            std::enable_if_t<
              !std::is_same<VectorImpl, std::decay_t<REP>>::value &&
              // I'm not sure why we had this before...
              // !std::is_base_of<VectorImpl<T, N, DefaultImpl>,
              //                  std::decay_t<REP>>::value
              // &&
              // !std::is_base_of<VectorImpl<T, N, VectorConstAliasImpl>,
              //                  std::decay_t<REP>>::value
              // &&
              // !std::is_base_of<VectorImpl<T, N, VectorAliasImpl>,
              //                  std::decay_t<REP>>::value
              // &&
              std::is_base_of<StcOp, std::decay_t<REP>>::value, int> = 0
            >
  HOSTDEVICE VectorImpl& operator=(REP&& a_rhs) noexcept
    {
      ArrayEq_E<std::decay_t<REP>> f(m_data, a_rhs);
      constexpr array_size_type I = cmin(VectorImpl::size(),
                                         std::decay_t<REP>::size());
      // If not for nvcc, here we could use a lambda since LHS cannot be
      // constexpr, e.g.:
      // auto f = [&](const array_size_type a_idx)
      //   {
      //     m_data[a_idx] = a_rhs[a_idx];
      //   };
      forEachElement<I>(f);
      return *this;
    }

  /// Evaluate an expression template given by Rep (+=)
  template <typename REP,
            std::enable_if_t<
              std::is_base_of<StcOp, std::decay_t<REP>>::value, int> = 0
            >
  HOSTDEVICE VectorImpl& operator+=(REP&& a_rhs) noexcept
    {
      this->operator=(*this + static_cast<const REP&>(a_rhs));
      return *this;
    }

  /// Evaluate an expression template given by Rep (-=)
  template <typename REP,
            std::enable_if_t<
              std::is_base_of<StcOp, std::decay_t<REP>>::value, int> = 0
            >
  HOSTDEVICE VectorImpl& operator-=(REP&& a_rhs) noexcept
    {
      this->operator=(*this - static_cast<const REP&>(a_rhs));
      return *this;
    }

  /// Evaluate an expression template given by Rep (*=)
  template <typename REP,
            std::enable_if_t<
              std::is_base_of<StcOp, std::decay_t<REP>>::value, int> = 0
            >
  HOSTDEVICE VectorImpl& operator*=(REP&& a_rhs) noexcept
    {
      this->operator=(*this * static_cast<const REP&>(a_rhs));
      return *this;
    }

  /// Evaluate an expression template given by Rep (/=)
  template <typename REP,
            std::enable_if_t<
              std::is_base_of<StcOp, std::decay_t<REP>>::value, int> = 0
            >
  HOSTDEVICE VectorImpl& operator/=(REP&& a_rhs) noexcept
    {
      this->operator=(*this / static_cast<const REP&>(a_rhs));
      return *this;
    }

//--x= scalars

  /// Assign constant
  HOSTDEVICE VectorImpl& operator=(const T& a_val) noexcept
    {
      this->operator=(stc::ScalarRep<T>(a_val));
      return *this;
    }

  /// Add scalar to each element
  HOSTDEVICE VectorImpl& operator+=(const T& a_val) noexcept
    {
      this->operator=(*this + a_val);
      return *this;
    }

  /// Subtract scalar from each element
  HOSTDEVICE VectorImpl& operator-=(const T& a_val) noexcept
    {
      this->operator=(*this - a_val);
      return *this;
    }

  /// Multiply scalar with each element
  HOSTDEVICE VectorImpl& operator*=(const T& a_val) noexcept
    {
      this->operator=(*this * a_val);
      return *this;
    }

  /// Divide each element by scalar
  HOSTDEVICE VectorImpl& operator/=(const T& a_val) noexcept
    {
      this->operator=(*this / a_val);
      return *this;
    }

//--Other modifiers

  /* Since these involve evaluation, all REP operands to expressions can be cast
     to constant reference going forward.  In other words, none of these
     results are used to build up further expressions.  '*this' must be an
     L-value and will naturally proceed as a reference. */

  /// Component-wise minimum with other Vector (modifies this)
  template <typename REP,
            std::enable_if_t<
              std::is_base_of<StcOp, std::decay_t<REP>>::value, int> = 0
            >
  HOSTDEVICE VectorImpl& min(REP&& a_rhs) noexcept
    {
      this->operator=(stc::min(*this, static_cast<const REP&>(a_rhs)));
      return *this;
    }

  /// Component-wise minimum with scalar (modifies this)
  HOSTDEVICE VectorImpl& min(const T& a_val) noexcept
    {
      this->operator=(stc::min(*this, a_val));
      return *this;
    }

  /// Component-wise maximum with other Vector (modifies this)
  template <typename REP,
            std::enable_if_t<
              std::is_base_of<StcOp, std::decay_t<REP>>::value, int> = 0
            >
  HOSTDEVICE VectorImpl& max(REP&& a_rhs) noexcept
    {
      this->operator=(stc::max(*this, static_cast<const REP&>(a_rhs)));
      return *this;
    }

  /// Component-wise maximum with scalar (modifies this)
  HOSTDEVICE VectorImpl& max(const T& a_val) noexcept
    {
      this->operator=(stc::max(*this, a_val));
      return *this;
    }

  /* The coarsen implemented here requires two complement, that >> is an
     arithmetic right shift for signed integers, and that the right hand side
     is a power of 2.  The former seems always true on modern architectures and
     compilers.
   */

  /// Modify vector by component-wise integer projection (modifies this)
  template <typename REP,
            std::enable_if_t<
              std::is_base_of<StcOp, std::decay_t<REP>>::value, int> = 0
            >
  HOSTDEVICE VectorImpl& coarsen(REP&& a_rhs) noexcept
    {
      this->operator=(stc::coarsen(*this, static_cast<const REP&>(a_rhs)));
      return *this;
    }

  /// Modify vector by integer projection (modifies this)
  HOSTDEVICE VectorImpl& coarsen(const unsigned long long a_factor) noexcept
    {
      this->operator=(stc::coarsen(*this, a_factor));
      return *this;
    }

  /// Scale (scalar multiplication)
  HOSTDEVICE VectorImpl& scale(const T& a_val) noexcept
    {
      this->operator*=(a_val);
      return *this;
    }

  /// Reflection
  /**
     Modifies vector by reflecting it in the plane defined by a_refIdx
     and with normal in the direction of a_dir.
     Directions are based at zero.
  */
  HOSTDEVICE VectorImpl& reflect(const T&        a_refIdx,
                                 array_size_type a_dir) noexcept
    {
      m_data[a_dir] = -m_data[a_dir] + 2*a_refIdx;
      return *this;
    }

  /// Shift by a_val in direction a_idx
  HOSTDEVICE VectorImpl& shift(array_size_type a_dir, const T& a_val) noexcept
    {
      m_data[a_dir] += a_val;
      return *this;
    }

  /// Shift by another vector (vector addition)
  template <typename REP,
            std::enable_if_t<
              std::is_base_of<StcOp, std::decay_t<REP>>::value, int> = 0
            >
  HOSTDEVICE VectorImpl& shift(REP&& a_rhs) noexcept
    {
      this->operator+=(static_cast<const REP&>(a_rhs));
      return *this;
    }

  /// Diagonal shift of vector (scalar addition)
  HOSTDEVICE VectorImpl& diagShift(const T& a_val) noexcept
    {
      this->operator+=(a_val);
      return *this;
    }

  /// Reverse elements
  HOSTDEVICE VectorImpl& reverse() noexcept
    {
      forEachElement<N/2>(
        [this]
        (const array_size_type a_idx)
          {
#ifdef __CUDA_ARCH__
            T tmp = m_data[a_idx];
            m_data[a_idx] = m_data[N - 1 - a_idx];
            m_data[N - 1 - a_idx] = tmp;
#else
            std::swap(m_data[a_idx], m_data[N - 1 - a_idx]);
#endif
          });
      return *this;
    }

//--Access

  /// Constant indexing
  HOSTDEVICE constexpr const T& operator[](
    const array_size_type a_idx) const noexcept
    {
      return m_data[a_idx];
    }

  /// Indexing
  template <typename _Impl = Impl<T, N>,
            std::enable_if_t<
              !std::is_same<_Impl,
                            VectorConstAliasImpl<T, N>>::value, int> = 0
            >
  HOSTDEVICE constexpr T& operator[](const array_size_type a_idx) noexcept
    {
      return m_data[a_idx];
    }

  /// Indexing
  /** Returns const T&, which is required for VectorConstAliasImpl
   */
  template <typename _Impl = Impl<T, N>,
            std::enable_if_t<
              std::is_same<_Impl,
                           VectorConstAliasImpl<T, N>>::value, int> = 0
            >
  HOSTDEVICE constexpr const T& operator[](const array_size_type a_idx) noexcept
    {
      return m_data[a_idx];
    }

  /// Constant indexing (function/"lambda-like")
  HOSTDEVICE constexpr const T& operator()(
    const array_size_type a_idx) const noexcept
    {
      return m_data[a_idx];
    }

  /// Indexing (function/"lambda-like")
  template <typename _Impl = Impl<T, N>,
            std::enable_if_t<
              !std::is_same<_Impl,
                            VectorConstAliasImpl<T, N>>::value, int> = 0
            >
  HOSTDEVICE constexpr T& operator()(const array_size_type a_idx) noexcept
    {
      return m_data[a_idx];
    }

  /// Indexing (function/"lambda-like")
  /** Returns const T&, which is required for VectorConstAliasImpl
   */
  template <typename _Impl = Impl<T, N>,
            std::enable_if_t<
              std::is_same<_Impl,
                           VectorConstAliasImpl<T, N>>::value, int> = 0
            >
  HOSTDEVICE constexpr const T& operator()(const array_size_type a_idx) noexcept
    {
      return m_data[a_idx];
    }

  /// Constant data pointer
  HOSTDEVICE const T* data() const noexcept
    {
      return &m_data[0];
    }

  /// Data pointer
  template <typename _Impl = Impl<T, N>,
            std::enable_if_t<
              !std::is_same<_Impl,
                            VectorConstAliasImpl<T, N>>::value, int> = 0
            >
  HOSTDEVICE T* data() noexcept
    {
      return &m_data[0];
    }

  /// Data pointer
  /** Returns const T*, which is required for VectorConstAliasImpl
   */
  template <typename _Impl = Impl<T, N>,
            std::enable_if_t<
              std::is_same<_Impl,
                           VectorConstAliasImpl<T, N>>::value, int> = 0
            >
  HOSTDEVICE const T* data() noexcept
    {
      return &m_data[0];
    }

  /// Constant data pointer
  HOSTDEVICE const T* dataPtr() const noexcept
    {
      return &m_data[0];
    }

  /// Data pointer
  template <typename _Impl = Impl<T, N>,
            std::enable_if_t<
              !std::is_same<_Impl,
                            VectorConstAliasImpl<T, N>>::value, int> = 0
            >
  HOSTDEVICE T* dataPtr() noexcept
    {
      return &m_data[0];
    }

  /// Data pointer
  /** Returns const T*, which is required for VectorConstAliasImpl
   */
  template <typename _Impl = Impl<T, N>,
            std::enable_if_t<
              std::is_same<_Impl,
                           VectorConstAliasImpl<T, N>>::value, int> = 0
            >
  HOSTDEVICE const T* dataPtr() noexcept
    {
      return &m_data[0];
    }

  /// Return the implementation (for copying to another aggregate)
  HOSTDEVICE constexpr const Impl<T, N>& getImpl() const
    {
      return m_data;
    }

  /// Constant begin iterator (DefaultImpl)
  template <typename _Impl = Impl<T, N>,
            std::enable_if_t<
              std::is_same<_Impl,
                           DefaultImpl<T, N>>::value, int> = 0
            >
  HOSTDEVICE constexpr const_iterator cbegin() const noexcept
    {
      return m_data.begin();
    }

  /// Constant begin iterator (VectorAliasImpl or VectorConstAliasImpl)
  template <typename _Impl = Impl<T, N>,
            std::enable_if_t<
              !std::is_same<_Impl,
                            DefaultImpl<T, N>>::value, int> = 0
            >
  HOSTDEVICE constexpr const_iterator cbegin() const noexcept
    {
      return m_data;
    }

  /// Constant begin iterator (DefaultImpl)
  template <typename _Impl = Impl<T, N>,
            std::enable_if_t<
              std::is_same<_Impl,
                           DefaultImpl<T, N>>::value, int> = 0
            >
  HOSTDEVICE constexpr const_iterator begin() const noexcept
    {
      return m_data.begin();
    }

  /// Constant begin iterator (VectorAliasImpl or VectorConstAliasImpl)
  template <typename _Impl = Impl<T, N>,
            std::enable_if_t<
              !std::is_same<_Impl,
                            DefaultImpl<T, N>>::value, int> = 0
            >
  HOSTDEVICE constexpr const_iterator begin() const noexcept
    {
      return m_data;
    }

  /// Begin iterator (DefaultImpl)
  template <typename _Impl = Impl<T, N>,
            std::enable_if_t<
              std::is_same<_Impl,
                           DefaultImpl<T, N>>::value, int> = 0
            >
  HOSTDEVICE iterator begin() noexcept
    {
      return m_data.begin();
    }

  /// Begin iterator (VectorAliasImpl)
  template <typename _Impl = Impl<T, N>,
            std::enable_if_t<
              std::is_same<_Impl,
                           VectorAliasImpl<T, N>>::value, int> = 0
            >
  HOSTDEVICE iterator begin() noexcept
    {
      return m_data;
    }

  /// Begin iterator (VectorConstAliasImpl)
  template <typename _Impl = Impl<T, N>,
            std::enable_if_t<
              std::is_same<_Impl,
                           VectorConstAliasImpl<T, N>>::value, int> = 0
            >
  HOSTDEVICE const_iterator begin() noexcept
    {
      return m_data;
    }

  /// Constant end iterator (DefaultImpl)
  template <typename _Impl = Impl<T, N>,
            std::enable_if_t<
              std::is_same<_Impl,
                           DefaultImpl<T, N>>::value, int> = 0
            >
  HOSTDEVICE constexpr const_iterator cend() const noexcept
    {
      return m_data.end();
    }

  /// Constant end iterator (VectorAliasImpl or VectorConstAliasImpl)
  template <typename _Impl = Impl<T, N>,
            std::enable_if_t<
              !std::is_same<_Impl,
                            DefaultImpl<T, N>>::value, int> = 0
            >
  HOSTDEVICE constexpr const_iterator cend() const noexcept
    {
      return m_data + N;
    }

  /// Constant end iterator (DefaultImpl)
  template <typename _Impl = Impl<T, N>,
            std::enable_if_t<
              std::is_same<_Impl,
                           DefaultImpl<T, N>>::value, int> = 0
            >
  HOSTDEVICE constexpr const_iterator end() const noexcept
    {
      return m_data.end();
    }

  /// Constant end iterator (VectorAliasImpl or VectorConstAliasImpl)
  template <typename _Impl = Impl<T, N>,
            std::enable_if_t<
              !std::is_same<_Impl,
                            DefaultImpl<T, N>>::value, int> = 0
            >
  HOSTDEVICE constexpr const_iterator end() const noexcept
    {
      return m_data + N;
    }

  /// End iterator (DefaultImpl)
  template <typename _Impl = Impl<T, N>,
            std::enable_if_t<
              std::is_same<_Impl,
                           DefaultImpl<T, N>>::value, int> = 0
            >
  HOSTDEVICE iterator end() noexcept
    {
      return m_data.end();
    }

  /// End iterator (VectorAliasImpl)
  template <typename _Impl = Impl<T, N>,
            std::enable_if_t<
              std::is_same<_Impl,
                           VectorAliasImpl<T, N>>::value, int> = 0
            >
  HOSTDEVICE iterator end() noexcept
    {
      return m_data + N;
    }

  /// End iterator (VectorConstAliasImpl)
  template <typename _Impl = Impl<T, N>,
            std::enable_if_t<
              std::is_same<_Impl,
                           VectorConstAliasImpl<T, N>>::value, int> = 0
            >
  HOSTDEVICE const_iterator end() noexcept
    {
      return m_data + N;
    }

//--Operations (but it is preferred to use the external variants)

  /// Sum of all components
  HOSTDEVICE constexpr T sum() const noexcept
    {
      return stc::sum(*this);
    }

  /// Product of all components
  HOSTDEVICE constexpr T product() const noexcept
    {
      return stc::product(*this);
    }

  /// 1-norm of all components
  HOSTDEVICE constexpr T norm1() const noexcept
    {
      return stc::norm1(*this);
    }

  /// Lexical less-than
  /** Returns true if this vector is lexically less than the argument.
      A vector MUST BE either lexically less than, lexically greater
      than, or equal to another IntVect.

      v1 is lexically less than v2 if:

      in 2-D:<br>
      v1[0] < v2[0] || (v1[0] == v2[0] && v1[1] < v2[1]);

      in 3-D:<br>
      v1[0] < v2[0] || (v1[0] == v2[0] && (v1[1] < v2[1] ||
        (v1[1] == v2[1] && v1[2] < v2[2])));
   */
  template <typename REP,
            std::enable_if_t<
              std::is_base_of<StcOp, std::decay_t<REP>>::value, int> = 0
            >
  HOSTDEVICE constexpr bool lexLT(REP&& a_rhs) const noexcept
    {
      return stc::lexLT(*this, static_cast<const REP&>(a_rhs));
    }

  /// Lexical greater-than
  /** Returns true if this IntVect is lexically greater than the
      argument.  An IntVect MUST BE either lexically less than,
      lexically greater than, or equal to another IntVect.

      v1 is lexically greater than v2 if:

      in 2-D:<br>
      v1[0] > v2[0] || (v1[0] == v2[0] && v1[1] > v2[1]);

      in 3-D:<br>
      v1[0] > v2[0] || (v1[0] == v2[0] && (v1[1] > v2[1] ||
        (v1[1] == v2[1] && v1[2] > v2[2])));
  */
  template <typename REP,
            std::enable_if_t<
              std::is_base_of<StcOp, std::decay_t<REP>>::value, int> = 0
            >
  HOSTDEVICE constexpr bool lexGT(REP&& a_rhs) const noexcept
    {
      return stc::lexGT(*this, static_cast<const REP&>(a_rhs));
    }

//--Depreciated

#ifndef __CUDACC__
  // Use copy or assignment constructor instead
  VectorImpl copy() const noexcept
    {
      return *this;
    }

  // Use operator[] instead
  void setVal(const array_size_type a_idx, const T& val) noexcept
    {
      m_data[a_idx] = val;
    }

  // set all values
  void setAll(const T& val) noexcept
    {
      *this = val;
    }
  
  // Use dataPtr functions instead
  const T* getVect() const noexcept
    {
      return dataPtr();
    }

  // Use mag instead
  constexpr auto vectorLength() const noexcept
    {
      return stc::mag(*this);
    }

  // Use dot product instead
  constexpr auto radSquared() const noexcept
    {
      return stc::dot(*this, *this);
    }
  
  // The minimum value of the vector
  constexpr auto min() const noexcept
    {
      return m_data[stc::indexMinElem(*this)];
    }

  // The maximum value of the vector
  constexpr auto max() const noexcept
    {
      return m_data[stc::indexMaxElem(*this)];
    }
  
  // Use indexMinElem instead (and abs if required)
  constexpr array_size_type minDir(const bool& a_doAbs) const noexcept
    {
      return (a_doAbs) ?
        stc::indexMinElem(stc::abs(*this)) :
        stc::indexMinElem(*this);
    }

  // Use indexMaxElem instead (and abs if required)
  constexpr array_size_type maxDir(const bool& a_doAbs) const noexcept
    {
      return (a_doAbs) ?
        stc::indexMaxElem(stc::abs(*this)) :
        stc::indexMaxElem(*this);
    }

  // Use dot instead
  template <typename REP,
            std::enable_if_t<
              std::is_base_of<StcOp, std::decay_t<REP>>::value, int> = 0
            >
  constexpr T dotProduct(REP&& a_rhs) const noexcept
    {
      return stc::dot(*this, static_cast<const REP&>(a_rhs));
    }
#endif

/*--------------------------------------------------------------------*
 * Member data
 *--------------------------------------------------------------------*/

//--Friend declarations so different implementations can access data

  template <typename T1, array_size_type R1,
            template <typename, array_size_type> class Impl1>
  friend class VectorImpl;

protected:

  Impl<T, N> m_data;                  ///< Underlying implementation giving data
                                      ///< for the vector
};


/*==============================================================================
 *
 * Types of vectors that are more accessible
 *
 *============================================================================*/

//--Vector of any type and size

template <typename T, array_size_type N>
using Vector = VectorImpl<T, N, DefaultImpl>;
template <typename T, array_size_type N>
using VectorAlias = VectorImpl<T, N, VectorAliasImpl>;
template <typename T, array_size_type N>
using VectorConstAlias = VectorImpl<T, N, VectorConstAliasImpl>;

//--Integer vector of any size

template <array_size_type N>
using IVec = VectorImpl<int, N, DefaultImpl>;
template <array_size_type N>
using IVecAlias = VectorImpl<int, N, VectorAliasImpl>;
template <array_size_type N>
using IVecConstAlias = VectorImpl<int, N, VectorConstAliasImpl>;

//--Real vector of any size

template <array_size_type N>
using RVec = VectorImpl<Real, N, DefaultImpl>;
template <array_size_type N>
using RVecAlias = VectorImpl<Real, N, VectorAliasImpl>;
template <array_size_type N>
using RVecConstAlias = VectorImpl<Real, N, VectorConstAliasImpl>;

//--Spatial vector of type Real and size spatial dimension

using SVec = VectorImpl<Real, SPACEDIM, DefaultImpl>;
using SVecAlias = VectorImpl<Real, SPACEDIM, VectorAliasImpl>;
using SVecConstAlias = VectorImpl<Real, SPACEDIM, VectorConstAliasImpl>;

//--Quaternion of type Real and size spatial dimension

using QVec = VectorImpl<Real, SPACEDIM+1, DefaultImpl>;
using QVecAlias = VectorImpl<Real, SPACEDIM+1, VectorAliasImpl>;
using QVecConstAlias = VectorImpl<Real, SPACEDIM+1, VectorConstAliasImpl>;

//--Helper to make alias types of a vector

// Modifiable alias
template <typename REP,
          std::enable_if_t<std::is_same_v<
                             typename std::decay_t<REP>::impl_type,
                             DefaultImpl<typename std::decay_t<REP>::value_type,
                                         std::decay_t<REP>::size()>>,
                           int> = 0>
struct make_alias
{
  using type = VectorImpl<typename std::decay_t<REP>::value_type,
                          std::decay_t<REP>::size(),
                          VectorAliasImpl>;
};

template <typename REP>
using make_alias_t = typename make_alias<REP>::type;

// Const alias
template <typename REP,
          std::enable_if_t<std::is_same_v<
                             typename std::decay_t<REP>::impl_type,
                             DefaultImpl<typename std::decay_t<REP>::value_type,
                                         std::decay_t<REP>::size()>>,
                           int> = 0>
struct make_constalias
{
  using type = VectorImpl<typename std::decay_t<REP>::value_type,
                          std::decay_t<REP>::size(),
                          VectorConstAliasImpl>;
};

template <typename REP>
using make_constalias_t = typename make_constalias<REP>::type;


/*******************************************************************************
 *
 * Makers for vectors
 *
 * Each maker constructs a special type of Vector, e.g, unit vectors.  There are
 * various ways to do this.  Here, we define static member functions in a class.
 * Each maker function uses a function object.  Again, lambdas are avoided since
 * they cannot be made constexpr until C++17.
 *
 *  \tparam     T       Type of array
 *  \tparam     N       Size of array
 *
 *//*+*************************************************************************/

/// Maker class supporting Vector Impl
template <typename T, array_size_type N>
struct make_vector
{
  /*
   * Using a lambda is ideal but we cannot make the constexpr until C++17
   */
  // static auto zero()
  //   {
  //     return VectorImpl<T, N, DefaultImpl>(
  //       make_array<T, N>(
  //         [](const array_size_type a_idx)
  //         {
  //           return (T)0;
  //         }));
  //   }
  /*
   * So instead, our function objects are declared before each function
   */

  struct Zero_E
  {
    HOSTDEVICE constexpr T operator()(
      const array_size_type a_idx) const noexcept
      {
        return (T)0;
      }
  };
  /// Set all elements to 0
  HOSTDEVICE static constexpr auto zero() noexcept
    {
      return Vector<T, N>(make_array<T, N>(Zero_E()));
    }

  struct Unit_E
  {
    HOSTDEVICE constexpr T operator()(
      const array_size_type a_idx) const noexcept
      {
        return (T)1;
      }
  };
  /// Set all elements to 1
  HOSTDEVICE static constexpr auto unit() noexcept
    {
      return Vector<T, N>(make_array<T, N>(Unit_E()));
    }

  struct Unit_x_E
  {
    HOSTDEVICE constexpr T operator()(
      const array_size_type a_idx) const noexcept
      {
        return (T)(a_idx == 0);
      }
  };
  /// Set first element to 1, others to 0
  HOSTDEVICE static constexpr auto unit_x() noexcept
    {
      return Vector<T, N>(make_array<T, N>(Unit_x_E()));
    }

  struct Unit_y_E
  {
    HOSTDEVICE constexpr T operator()(
      const array_size_type a_idx) const noexcept
      {
        return (T)(a_idx == 1);
      }
  };
  /// Set second element to 1, others to 0
  HOSTDEVICE static constexpr auto unit_y() noexcept
    {
      return Vector<T, N>(make_array<T, N>(Unit_y_E()));
    }

  struct Unit_z_E
  {
    HOSTDEVICE constexpr T operator()(
      const array_size_type a_idx) const noexcept
      {
        return (T)(a_idx == 2);
      }
  };
  /// Set third element to 1, others to 0
  HOSTDEVICE static constexpr auto unit_z() noexcept
    {
      return Vector<T, N>(make_array<T, N>(Unit_z_E()));
    }

  struct Unit_w_E
  {
    HOSTDEVICE constexpr T operator()(
      const array_size_type a_idx) const noexcept
      {
        return (T)(a_idx == 3);
      }
  };
  /// Set fourth element to 1, others to 0
  HOSTDEVICE static constexpr auto unit_w() noexcept
    {
      return Vector<T, N>(make_array<T, N>(Unit_w_E()));
    }

  struct Basis_E
  {
    HOSTDEVICE constexpr Basis_E(const array_size_type a_eidx) noexcept
      :
      m_eidx(a_eidx)
      { }
    HOSTDEVICE constexpr T operator()(
      const array_size_type a_idx) const noexcept
      {
        return (T)(a_idx == m_eidx);
      }
    const array_size_type m_eidx;
  };
  /// Set argument index to 1, others to 0
  HOSTDEVICE static constexpr auto basis(const array_size_type a_eidx) noexcept
    {
      return Vector<T, N>(make_array<T, N>(Basis_E(a_eidx)));
    }

  struct Fill_E
  {
    HOSTDEVICE constexpr Fill_E(const T& a_val) noexcept
      :
      m_val(a_val)
      { }
    HOSTDEVICE constexpr T operator()(
      const array_size_type a_idx) const noexcept
      {
        return m_val;
      }
    const T m_val;
  };
  /// Fill all values
  HOSTDEVICE static constexpr auto fill(const T& a_val) noexcept
    {
      return Vector<T, N>(make_array<T, N>(Fill_E(a_val)));
    }

  template <int SzVec, typename Vec, typename... Idxs>
  struct VecAndIdxs_E
  {
    HOSTDEVICE constexpr VecAndIdxs_E(
      const Vec&     a_vec,
      const Idxs&... a_idxs) noexcept
      :
      m_vec(a_vec),
      m_idxs({ a_idxs... })
      { }
    HOSTDEVICE constexpr int operator()(
      const array_size_type a_idx) const noexcept
      {
        return (a_idx < SzVec) ? m_vec[a_idx] : m_idxs[a_idx - SzVec];
      }
    const Vec& m_vec;
    DefaultImpl<T, sizeof...(Idxs)> m_idxs;
  };
  /// Compose from a vector and standalone indices
  /** The number of entries taken from a_vec is N - sizeof...(Idxs)
   */
  template <typename Vec, typename... Idxs>
  HOSTDEVICE static constexpr auto vecAndIdxs(const Vec&     a_vec,
                                              const Idxs&... a_idxs) noexcept
    {
      return Vector<T, N>(
        make_array<T, N>(
          VecAndIdxs_E<N - sizeof...(Idxs), Vec, Idxs...>(a_vec, a_idxs...)));
    }

  template <typename Vec>
  struct MDIdxStride_E
  {
    HOSTDEVICE constexpr MDIdxStride_E(
      stc::Vector<T, N>& a_mdidx,
      const Vec&         a_stride,
      array_size_type&   a_lin) noexcept
      :
      m_mdidx(a_mdidx),
      m_stride(a_stride),
      m_lin(a_lin)
      { }
    HOSTDEVICE constexpr void operator()(
      const array_size_type a_idx) const noexcept
      {
        const int idx = a_idx + 1;
        m_mdidx[idx] = m_lin/m_stride[idx];
        m_lin -= m_mdidx[idx]*m_stride[idx];
      }
    stc::Vector<T, N>& m_mdidx;
    const Vec& m_stride;
    array_size_type& m_lin;
  };
  /// Compose a multidimensional vector index from a linear one with strides
  /** See also @c fromLinearLog2Dim
   */
  template <typename Vec>
  HOSTDEVICE static constexpr auto mdidx(
    array_size_type a_lin,
    const Vec&      a_stride) noexcept
    {
      stc::Vector<T, N> mdidxRet(0);
      forEachReverseElement<N-1>(MDIdxStride_E<Vec>(mdidxRet, a_stride, a_lin));
      mdidxRet[0] = a_lin;
      return mdidxRet;
    }
};

//--More accessible definitions

template <array_size_type N>
using make_IVec = make_vector<int, N>;
template <array_size_type N>
using make_RVec = make_vector<Real, N>;
using make_SVec = make_vector<Real, SPACEDIM>;
using make_QVec = make_vector<Real, SPACEDIM+1>;

// Handy definitions
constexpr auto SVec_zero = make_SVec::zero();
constexpr auto SVec_unit = make_SVec::unit();
constexpr auto SVec_unit_x = make_SVec::unit_x();
constexpr auto SVec_unit_y = make_SVec::unit_y();
constexpr auto SVec_unit_z = make_SVec::unit_z();
constexpr auto QVec_zero = make_QVec::zero();
constexpr auto QVec_unit = make_QVec::unit();


/*******************************************************************************
 */
/// Vector (or matrix) representation of a scalar
/** Scalars are promoted to this class for use in expressions
 *
 *  \tparam T           Type of scalar
 *
 *//*+*************************************************************************/

template <typename T>
class ScalarRep : StcOp
{

//--Types

public:

  using value_type = T;
  using is_scalar = std::true_type;
  using is_vector = std::false_type;
  using is_rowvec = std::false_type;
  using is_matrix = std::false_type;
  using is_string = std::false_type;

//--Constructors and destructors

  /// Constructor
  HOSTDEVICE constexpr ScalarRep(const T& a_scalar) noexcept
    :
    m_scalar(a_scalar)
    { }
  /// Copy constructor
  constexpr ScalarRep(const ScalarRep&) noexcept = default;
  /// Move constructor
  constexpr ScalarRep(ScalarRep&&) noexcept = default;
  /// Assignment constructor
  ScalarRep& operator=(const ScalarRep&) = delete;
  /// Move assignment constructor
  ScalarRep& operator=(ScalarRep&&) noexcept = delete;
  /// Destructor
  ~ScalarRep() = default;

//--Member functions

  /// Size of the vector representation of the scalar
  /** Since scalars don't have a size, it is set to a maximum
   */
  HOSTDEVICE static constexpr array_size_type size() noexcept
    {
#ifdef __CUDA_ARCH__
      // array_size_type is expected to be unsigned
      return NPP_MAX_32U;
#else
      return std::numeric_limits<array_size_type>::max();
#endif
    }
  HOSTDEVICE static constexpr array_size_type rank() noexcept
    {
#ifdef __CUDA_ARCH__
      // array_size_type is expected to be unsigned
      return NPP_MAX_32U;
#else
      return std::numeric_limits<array_size_type>::max();
#endif
    }

  /// Access
  HOSTDEVICE constexpr auto operator[](const int a_idx) const noexcept
    {
      return m_scalar;
    }

  /// Access (function/"lambda-like")
  HOSTDEVICE constexpr auto operator()(const int a_idx) const noexcept
    {
      return m_scalar;
    }

  const T m_scalar;                   ///< Stored scalar
};


/*******************************************************************************
 *
 * Operator Representations
 *
 * There is a single representation for operator types, e.g., binary and unary.
 * Each operator is defined by a function that creates the appropriate
 * representations with a function object.  Lambdas could be used except we make
 * use of constexpr and that is not supported until C++17.  You will find the
 * function objects before each operator function.
 *
 * We are quite careful about moving rvalues (temporary objects) which are
 * stored in class, and copying lvalue references, which are stored as
 * references.  In otherwords, anything that started as an lvalue is referred to
 * throughout the chain of representations; anything that started as an rvalue
 * is stored in the representation class and moved when required (the original
 * ceases to exist).
 *
 * Enable_if is used in the operators.  The operands must derive from StcOp.
 * Otherwise, the operators bind to everything under the sun...
 *
 ******************************************************************************/


/*============================================================================*/
/// Representation of a binary operation between vectors
/** \tparam REP1        Left operand
 *  \tparam REP2        Right operand
 *  \tparam BinaryOp    Function object performing a desired binary operation
 *
 *  \note
 *  <ul>
 *    <li> The function object operates on the underyling types which must be
 *         available through the representations
 *  </ul>
 *  
 *//*=========================================================================*/

//--For basic Ops, the representations are evaluated at the index before
//--invoking the BinaryOp operator.

template <typename REP1, typename REP2, typename BinaryOp>
class VectorBinaryOpRep : public BinaryOp, private StcOp
{
public:
  using value_type = typename BinaryOp::value_type;
  using is_scalar = std::false_type;
  using is_vector = std::true_type;
  using is_rowvec = std::false_type;
  using is_matrix = std::false_type;
  using is_string = std::false_type;
private:
  /* Store REP instead of REP&&.  The former stores lvalues as a reference and
   * rvalues directly.  The latter stores lvalue and rvalue references.  We
   * decide to use the former because it allows for temporaries to persist
   * across a scope.  Rvalues are automatically moved, lvalues can be explicitly
   * moved.
   */
  REP1 m_rep1;
  REP2 m_rep2;
public:
  template <typename... OpArgs>
  HOSTDEVICE constexpr VectorBinaryOpRep(REP1&&      a_rep1,
                                         REP2&&      a_rep2,
                                         OpArgs&&... a_opArgs) noexcept
    :
    BinaryOp(std::forward<OpArgs>(a_opArgs)...),
    m_rep1(std::forward<REP1>(a_rep1)),
    m_rep2(std::forward<REP2>(a_rep2))
    { }
  HOSTDEVICE constexpr VectorBinaryOpRep(
    const VectorBinaryOpRep& a_other) noexcept
    :
    BinaryOp(a_other),
    m_rep1(a_other.m_rep1),
    m_rep2(a_other.m_rep2)
    { }
  HOSTDEVICE constexpr VectorBinaryOpRep(VectorBinaryOpRep&& a_other) noexcept
    :
    BinaryOp(std::move(a_other)),
    m_rep1(std::forward<REP1>(a_other.m_rep1)),
    m_rep2(std::forward<REP2>(a_other.m_rep2))
    { }
  /// Assignment constructor
  VectorBinaryOpRep& operator=(const VectorBinaryOpRep&) = delete;
  /// Move assignment constructor
  VectorBinaryOpRep& operator=(VectorBinaryOpRep&&) noexcept = delete;
  /// Size of the vector representation (minimum of Reps)
  HOSTDEVICE static constexpr array_size_type size() noexcept
    {
      return cmin(std::decay_t<REP1>::size(), std::decay_t<REP2>::size());
    }

//--Access

  /// Array indexing
  HOSTDEVICE constexpr auto operator[](
    const array_size_type a_idx) const noexcept
    {
      return this->eval(m_rep1[a_idx], m_rep2[a_idx]);
    }
  /// Function/"lambda-like"
  HOSTDEVICE constexpr auto operator()(
    const array_size_type a_idx) const noexcept
    {
      return this->eval(m_rep1(a_idx), m_rep2(a_idx));
    }

//--Operations

  /// Sum of all components
  HOSTDEVICE constexpr value_type sum() const noexcept
    {
      return stc::sum(*this);
    }
  /// Product of all components
  HOSTDEVICE constexpr value_type product() const noexcept
    {
      return stc::product(*this);
    }
  /// 1-norm of all components
  HOSTDEVICE constexpr value_type norm1() const noexcept
    {
      return stc::norm1(*this);
    }
};

//--For ROps, the representations are not evaluated at the index before
//--passing to the BinaryROp operator.  This is often used for more complicated
//--operations where multiple indices of the Reps must be accessed to compute
//--a single value (e.g., matrix-vector multiplication), but we still want to
//--compose an expression.

template <typename REP1, typename REP2, typename BinaryROp>
class VectorBinaryROpRep : public BinaryROp, private StcOp
{
public:
  using value_type = typename BinaryROp::value_type;
  using is_scalar = std::false_type;
  using is_vector = std::true_type;
  using is_rowvec = std::false_type;
  using is_matrix = std::false_type;
  using is_string = std::false_type;
private:
  REP1 m_rep1;
  REP2 m_rep2;
public:
  template <typename... OpArgs>
  HOSTDEVICE constexpr VectorBinaryROpRep(REP1&&      a_rep1,
                                          REP2&&      a_rep2,
                                          OpArgs&&... a_opArgs) noexcept
    :
    BinaryROp(std::forward<OpArgs>(a_opArgs)...),
    m_rep1(std::forward<REP1>(a_rep1)),
    m_rep2(std::forward<REP2>(a_rep2))
    { }
  HOSTDEVICE constexpr VectorBinaryROpRep(
    const VectorBinaryROpRep& a_other) noexcept
    :
    BinaryROp(a_other),
    m_rep1(a_other.m_rep1),
    m_rep2(a_other.m_rep2)
    { }
  HOSTDEVICE constexpr VectorBinaryROpRep(
    VectorBinaryROpRep&& a_other) noexcept
    :
    BinaryROp(std::move(a_other)),
    m_rep1(std::forward<REP1>(a_other.m_rep1)),
    m_rep2(std::forward<REP2>(a_other.m_rep2))
    { }
  /// Assignment constructor
  VectorBinaryROpRep& operator=(const VectorBinaryROpRep&) = delete;
  /// Move assignment constructor
  VectorBinaryROpRep& operator=(VectorBinaryROpRep&&) noexcept = delete;
  /// Size of the vector representation (minimum of Reps)
  HOSTDEVICE static constexpr array_size_type size() noexcept
    {
      return cmin(std::decay_t<REP1>::size(), std::decay_t<REP2>::size());
    }

//--Access

  /// Array indexing
  HOSTDEVICE constexpr auto operator[](
    const array_size_type a_idx) const noexcept
    {
      return this->eval(m_rep1, m_rep2, a_idx);
    }
  /// Function/"lambda-like"
  HOSTDEVICE constexpr auto operator()(
    const array_size_type a_idx) const noexcept
    {
      return this->eval(m_rep1, m_rep2, a_idx);
    }

//--Operations

  /// Sum of all components
  HOSTDEVICE constexpr value_type sum() const noexcept
    {
      return stc::sum(*this);
    }
  /// Product of all components
  HOSTDEVICE constexpr value_type product() const noexcept
    {
      return stc::product(*this);
    }
  /// 1-norm of all components
  HOSTDEVICE constexpr value_type norm1() const noexcept
    {
      return stc::norm1(*this);
    }
};


/*============================================================================*/
/// Representation of a unary operation on a vector
/** \tparam REP         Operand
 *  \tparam UnaryOp     Function object performing a desired unary operation
 *
 *  \note
 *  <ul>
 *    <li> The function object operates on the underyling type which must be
 *         available through the representation
 *  </ul>
 *  
 *//*=========================================================================*/

//--For basic Ops, the representations are evaluated at the index before
//--invoking the UnaryOp operator.

template <typename REP, typename UnaryOp>
class VectorUnaryOpRep : public UnaryOp, private StcOp
{
public:
  using value_type = typename UnaryOp::value_type;
  using is_scalar = std::false_type;
  using is_vector = std::true_type;
  using is_rowvec = std::false_type;
  using is_matrix = std::false_type;
  using is_string = std::false_type;
private:
  REP m_rep;
public:
  template <typename... OpArgs>
  HOSTDEVICE constexpr VectorUnaryOpRep(REP&&       a_rep,
                                        OpArgs&&... a_opArgs) noexcept
    :
    UnaryOp(std::forward<OpArgs>(a_opArgs)...),
    m_rep(std::forward<REP>(a_rep))
    { }
  HOSTDEVICE constexpr VectorUnaryOpRep(
    const VectorUnaryOpRep& a_other) noexcept
    :
    UnaryOp(a_other),
    m_rep(a_other.m_rep)
    { }
  HOSTDEVICE constexpr VectorUnaryOpRep(VectorUnaryOpRep&& a_other) noexcept
    :
    UnaryOp(std::move(a_other)),
    m_rep(std::forward<REP>(a_other.m_rep))
    { }
  /// Assignment constructor
  VectorUnaryOpRep& operator=(const VectorUnaryOpRep&) = delete;
  /// Move assignment constructor
  VectorUnaryOpRep& operator=(VectorUnaryOpRep&&) noexcept = delete;
  /// Size of the vector representation
  HOSTDEVICE static constexpr array_size_type size() noexcept
    {
      return std::decay_t<REP>::size();
    }

//--Access

  /// Array indexing
  HOSTDEVICE constexpr auto operator[](
    const array_size_type a_idx) const noexcept
    {
      return this->eval(m_rep[a_idx]);
    }
  /// Function/"lambda-like"
  HOSTDEVICE constexpr auto operator()(
    const array_size_type a_idx) const noexcept
    {
      return this->eval(m_rep(a_idx));
    }

//--Operations

  /// Sum of all components
  HOSTDEVICE constexpr value_type sum() const noexcept
    {
      return stc::sum(*this);
    }
  /// Product of all components
  HOSTDEVICE constexpr value_type product() const noexcept
    {
      return stc::product(*this);
    }
  /// 1-norm of all components
  HOSTDEVICE constexpr value_type norm1() const noexcept
    {
      return stc::norm1(*this);
    }
};

//--For ROps, the representation is not evaluated at the index before
//--passing to the UnaryROp operator.  This is often used for more complicated
//--operations where multiple indices of the Rep must be accessed to compute
//--a single value (e.g., matrix-vector multiplication), but we still want to
//--compose an expression.

template <typename REP, typename UnaryROp>
class VectorUnaryROpRep : public UnaryROp, private StcOp
{
public:
  using value_type = typename UnaryROp::value_type;
  using is_scalar = std::false_type;
  using is_vector = std::true_type;
  using is_rowvec = std::false_type;
  using is_matrix = std::false_type;
  using is_string = std::false_type;
private:
  REP m_rep;
public:
  template <typename... OpArgs>
  HOSTDEVICE constexpr VectorUnaryROpRep(REP&&       a_rep,
                                         OpArgs&&... a_opArgs) noexcept
    :
    UnaryROp(std::forward<OpArgs>(a_opArgs)...),
    m_rep(std::forward<REP>(a_rep))
    { }
  HOSTDEVICE constexpr VectorUnaryROpRep(
    const VectorUnaryROpRep& a_other) noexcept
    :
    UnaryROp(a_other),
    m_rep(a_other.m_rep)
    { }
  HOSTDEVICE constexpr VectorUnaryROpRep(VectorUnaryROpRep&& a_other) noexcept
    :
    UnaryROp(std::move(a_other)),
    m_rep(std::forward<REP>(a_other.m_rep))
    { }
  /// Assignment constructor
  VectorUnaryROpRep& operator=(const VectorUnaryROpRep&) = delete;
  /// Move assignment constructor
  VectorUnaryROpRep& operator=(VectorUnaryROpRep&&) noexcept = delete;
  /// Size of the vector representation
  HOSTDEVICE static constexpr array_size_type size() noexcept
    {
      return std::decay_t<REP>::size();
    }

//--Access

  /// Array indexing
  HOSTDEVICE constexpr auto operator[](
    const array_size_type a_idx) const noexcept
    {
      return this->eval(m_rep, a_idx);
    }
  /// Function/"lambda-like"
  HOSTDEVICE constexpr auto operator()(
    const array_size_type a_idx) const noexcept
    {
      return this->eval(m_rep, a_idx);
    }

//--Operations

  /// Sum of all components
  HOSTDEVICE constexpr value_type sum() const noexcept
    {
      return stc::sum(*this);
    }
  /// Product of all components
  HOSTDEVICE constexpr value_type product() const noexcept
    {
      return stc::product(*this);
    }
  /// 1-norm of all components
  HOSTDEVICE constexpr value_type norm1() const noexcept
    {
      return stc::norm1(*this);
    }
};


/*******************************************************************************
 *
 * Operators
 *
 ******************************************************************************/

/*--------------------------------------------------------------------*
 * Addition
 *--------------------------------------------------------------------*/

template <typename T1, typename T2>
struct BinaryOpAdd
{
  using value_type = decltype(T1() + T2());
  HOSTDEVICE constexpr auto eval(const T1& a_x1, const T2& a_x2) const noexcept
    {
      return a_x1 + a_x2;
    }
};

//--Operator between two vector representations

template <typename REP1, typename REP2,
          std::enable_if_t<
            std::is_base_of<stc::StcOp, std::decay_t<REP1>>::value
            &&
            std::is_base_of<stc::StcOp, std::decay_t<REP2>>::value
            &&
            std::decay_t<REP1>::is_vector::value
            &&
            std::decay_t<REP2>::is_vector::value
            &&
            !std::decay_t<REP1>::is_string::value
            &&
            !std::decay_t<REP2>::is_string::value, int>
          >
HOSTDEVICE constexpr auto
operator+(REP1&& a_rep1, REP2&& a_rep2) noexcept
{
  return stc::VectorBinaryOpRep<REP1, REP2,
                                stc::BinaryOpAdd<
                                  typename std::decay_t<REP1>::value_type,
                                  typename std::decay_t<REP2>::value_type>>(
    std::forward<REP1>(a_rep1),
    std::forward<REP2>(a_rep2));
}

//--Operator between left scalar and right vector representation

template <typename T1, typename REP2,
          std::enable_if_t<
            !std::is_base_of<stc::StcOp, std::decay_t<T1>>::value
            &&
            std::is_base_of<stc::StcOp, std::decay_t<REP2>>::value
            &&
            !std::is_class<std::decay_t<T1>>::value
            &&
            std::decay_t<REP2>::is_vector::value
            &&
            !std::decay_t<REP2>::is_string::value, int> = 0
          >
HOSTDEVICE constexpr auto
operator+(const T1& a_scalar, REP2&& a_rep2) noexcept
{
  return stc::VectorBinaryOpRep<stc::ScalarRep<T1>, REP2,
                                stc::BinaryOpAdd<
                                  T1,
                                  typename std::decay_t<REP2>::value_type>>(
    stc::ScalarRep<T1>(a_scalar),
    std::forward<REP2>(a_rep2));
}

//--Operator between left vector representation and right scalar

template <typename REP1, typename T2,
          std::enable_if_t<
            std::is_base_of<stc::StcOp, std::decay_t<REP1>>::value
            &&
            !std::is_base_of<stc::StcOp, std::decay_t<T2>>::value
            &&
            std::decay_t<REP1>::is_vector::value
            &&
            !std::decay_t<REP1>::is_string::value
            &&
            !std::is_class<std::decay_t<T2>>::value, int>
          >
HOSTDEVICE constexpr auto
operator+(REP1&& a_rep1, const T2& a_scalar) noexcept
{
  return stc::VectorBinaryOpRep<REP1, stc::ScalarRep<T2>,
                                stc::BinaryOpAdd<
                                  typename std::decay_t<REP1>::value_type,
                                  T2>>(
    std::forward<REP1>(a_rep1),
    stc::ScalarRep<T2>(a_scalar));
}

/*--------------------------------------------------------------------*
 * Subtraction
 *--------------------------------------------------------------------*/
  
template <typename T1, typename T2>
struct BinaryOpSubtract
{
  using value_type = decltype(T1() - T2());
  HOSTDEVICE constexpr auto eval(const T1& a_x1, const T2& a_x2) const noexcept
    {
      return a_x1 - a_x2;
    }
};

//--Operator between two vector representations

template <typename REP1, typename REP2,
          std::enable_if_t<
            std::is_base_of<stc::StcOp, std::decay_t<REP1>>::value
            &&
            std::is_base_of<stc::StcOp, std::decay_t<REP2>>::value
            &&
            std::decay_t<REP1>::is_vector::value
            &&
            std::decay_t<REP2>::is_vector::value, int>
          >
HOSTDEVICE constexpr auto
operator-(REP1&& a_rep1, REP2&& a_rep2) noexcept
{
  return stc::VectorBinaryOpRep<REP1, REP2,
                                stc::BinaryOpSubtract<
                                  typename std::decay_t<REP1>::value_type,
                                  typename std::decay_t<REP2>::value_type>>(
    std::forward<REP1>(a_rep1),
    std::forward<REP2>(a_rep2));
}

//--Operator between left scalar and right vector representation

template <typename T1, typename REP2,
          std::enable_if_t<
            !std::is_base_of<stc::StcOp, std::decay_t<T1>>::value
            &&
            std::is_base_of<stc::StcOp, std::decay_t<REP2>>::value
            &&
            !std::is_class<std::decay_t<T1>>::value
            &&
            std::decay_t<REP2>::is_vector::value, int> = 0
          >
HOSTDEVICE constexpr auto
operator-(const T1& a_scalar, REP2&& a_rep2) noexcept
{
  return stc::VectorBinaryOpRep<stc::ScalarRep<T1>, REP2,
                                stc::BinaryOpSubtract<
                                  T1,
                                  typename std::decay_t<REP2>::value_type>>(
    stc::ScalarRep<T1>(a_scalar),
    std::forward<REP2>(a_rep2));
}

//--Operator between left vector representation and right scalar

template <typename REP1, typename T2,
          std::enable_if_t<
            std::is_base_of<stc::StcOp, std::decay_t<REP1>>::value
            &&
            !std::is_base_of<stc::StcOp, std::decay_t<T2>>::value
            &&
            std::decay_t<REP1>::is_vector::value
            &&
            !std::is_class<std::decay_t<T2>>::value, int>
          >
HOSTDEVICE constexpr auto
operator-(REP1&& a_rep1, const T2& a_scalar) noexcept
{
  return stc::VectorBinaryOpRep<REP1, stc::ScalarRep<T2>,
                                stc::BinaryOpSubtract<
                                  typename std::decay_t<REP1>::value_type,
                                  T2>>(
    std::forward<REP1>(a_rep1),
    stc::ScalarRep<T2>(a_scalar));
}

/*--------------------------------------------------------------------*
 * Multiplication
 *--------------------------------------------------------------------*/

template <typename T1, typename T2>
struct BinaryOpMultiply
{
  using value_type = decltype(T1()*T2());
  HOSTDEVICE constexpr auto eval(const T1& a_x1, const T2& a_x2) const noexcept
    {
      return a_x1*a_x2;
    }
};

//--Operator between two vector representations

template <typename REP1, typename REP2,
          std::enable_if_t<
            std::is_base_of<stc::StcOp, std::decay_t<REP1>>::value
            &&
            std::is_base_of<stc::StcOp, std::decay_t<REP2>>::value
            &&
            std::decay_t<REP1>::is_vector::value
            &&
            std::decay_t<REP2>::is_vector::value, int>
          >
HOSTDEVICE constexpr auto
operator*(REP1&& a_rep1, REP2&& a_rep2) noexcept
{
  return stc::VectorBinaryOpRep<REP1, REP2,
                                stc::BinaryOpMultiply<
                                  typename std::decay_t<REP1>::value_type,
                                  typename std::decay_t<REP2>::value_type>>(
    std::forward<REP1>(a_rep1),
    std::forward<REP2>(a_rep2));
}

//--Operator between left scalar and right vector representation

template <typename T1, typename REP2,
          std::enable_if_t<
            !std::is_base_of<stc::StcOp, std::decay_t<T1>>::value
            &&
            std::is_base_of<stc::StcOp, std::decay_t<REP2>>::value
            &&
            !std::is_class<std::decay_t<T1>>::value
            &&
            std::decay_t<REP2>::is_vector::value, int> = 0
          >
HOSTDEVICE constexpr auto
operator*(const T1& a_scalar, REP2&& a_rep2) noexcept
{
  return stc::VectorBinaryOpRep<stc::ScalarRep<T1>, REP2,
                                stc::BinaryOpMultiply<
                                  T1,
                                  typename std::decay_t<REP2>::value_type>>(
    stc::ScalarRep<T1>(a_scalar),
    std::forward<REP2>(a_rep2));
}

//--Operator between left vector representation and right scalar

template <typename REP1, typename T2,
          std::enable_if_t<
            std::is_base_of<stc::StcOp, std::decay_t<REP1>>::value
            &&
            !std::is_base_of<stc::StcOp, std::decay_t<T2>>::value
            &&
            std::decay_t<REP1>::is_vector::value
            &&
            !std::is_class<std::decay_t<T2>>::value, int>
          >
HOSTDEVICE constexpr auto
operator*(REP1&& a_rep1, const T2& a_scalar) noexcept
{
  return stc::VectorBinaryOpRep<REP1, stc::ScalarRep<T2>,
                                stc::BinaryOpMultiply<
                                  typename std::decay_t<REP1>::value_type,
                                  T2>>(
    std::forward<REP1>(a_rep1),
    stc::ScalarRep<T2>(a_scalar));
}

/*--------------------------------------------------------------------*
 * Division
 *--------------------------------------------------------------------*/

template <typename T1, typename T2>
struct BinaryOpDivision
{
  using value_type = decltype(T1() / T2());
  HOSTDEVICE constexpr auto eval(const T1& a_x1, const T2& a_x2) const noexcept
    {
      return a_x1/a_x2;
    }
};

//--Operator between two vector representations

template <typename REP1, typename REP2,
          std::enable_if_t<
            std::is_base_of<stc::StcOp, std::decay_t<REP1>>::value
            &&
            std::is_base_of<stc::StcOp, std::decay_t<REP2>>::value
            &&
            std::decay_t<REP1>::is_vector::value
            &&
            std::decay_t<REP2>::is_vector::value, int>
          >
HOSTDEVICE constexpr auto
operator/(REP1&& a_rep1, REP2&& a_rep2) noexcept
{
  return stc::VectorBinaryOpRep<REP1, REP2,
                                stc::BinaryOpDivision<
                                  typename std::decay_t<REP1>::value_type,
                                  typename std::decay_t<REP2>::value_type>>(
    std::forward<REP1>(a_rep1),
    std::forward<REP2>(a_rep2));
}

//--Operator between left scalar and right vector representation

template <typename T1, typename REP2,
          std::enable_if_t<
            !std::is_base_of<stc::StcOp, std::decay_t<T1>>::value
            &&
            std::is_base_of<stc::StcOp, std::decay_t<REP2>>::value
            &&
            !std::is_class<std::decay_t<T1>>::value
            &&
            std::decay_t<REP2>::is_vector::value, int> = 0
          >
HOSTDEVICE constexpr auto
operator/(const T1& a_scalar, REP2&& a_rep2) noexcept
{
  return stc::VectorBinaryOpRep<stc::ScalarRep<T1>, REP2,
                                stc::BinaryOpDivision<
                                  T1,
                                  typename std::decay_t<REP2>::value_type>>(
    stc::ScalarRep<T1>(a_scalar),
    std::forward<REP2>(a_rep2));
}

//--Operator between left vector representation and right scalar

template <typename REP1, typename T2,
          std::enable_if_t<
            std::is_base_of<stc::StcOp, std::decay_t<REP1>>::value
            &&
            !std::is_base_of<stc::StcOp, std::decay_t<T2>>::value
            &&
            std::decay_t<REP1>::is_vector::value
            &&
            !std::is_class<std::decay_t<T2>>::value, int>
          >
HOSTDEVICE constexpr auto
operator/(REP1&& a_rep1, const T2& a_scalar) noexcept
{
  return stc::VectorBinaryOpRep<REP1, stc::ScalarRep<T2>,
                                stc::BinaryOpDivision<
                                  typename std::decay_t<REP1>::value_type,
                                  T2>>(
    std::forward<REP1>(a_rep1),
    stc::ScalarRep<T2>(a_scalar));
}

/*--------------------------------------------------------------------*
 * Modulus
 *--------------------------------------------------------------------*/

template <typename T1, typename T2>
struct BinaryOpModulus
{
  using value_type = decltype(T1() % T2());
  HOSTDEVICE constexpr auto eval(const T1& a_x1, const T2& a_x2) const noexcept
    {
      return a_x1%a_x2;
    }
};

//--Operator between two vector representations

template <typename REP1, typename REP2,
          std::enable_if_t<
            std::is_base_of<stc::StcOp, std::decay_t<REP1>>::value
            &&
            std::is_base_of<stc::StcOp, std::decay_t<REP2>>::value
            &&
            std::decay_t<REP1>::is_vector::value
            &&
            std::decay_t<REP2>::is_vector::value, int> = 0
          >
HOSTDEVICE constexpr auto
operator%(REP1&& a_rep1, REP2&& a_rep2) noexcept
{
  return stc::VectorBinaryOpRep<REP1, REP2,
                                stc::BinaryOpModulus<
                                  typename std::decay_t<REP1>::value_type,
                                  typename std::decay_t<REP2>::value_type>>(
    std::forward<REP1>(a_rep1),
    std::forward<REP2>(a_rep2));
}

//--Operator between left scalar and right vector representation

template <typename T1, typename REP2,
          std::enable_if_t<
            !std::is_base_of<stc::StcOp, std::decay_t<T1>>::value
            &&
            std::is_base_of<stc::StcOp, std::decay_t<REP2>>::value
            &&
            !std::is_class<std::decay_t<T1>>::value
            &&
            std::decay_t<REP2>::is_vector::value, int> = 0
          >
HOSTDEVICE constexpr auto
operator%(const T1& a_scalar, REP2&& a_rep2) noexcept
{
  return stc::VectorBinaryOpRep<stc::ScalarRep<T1>, REP2,
                                stc::BinaryOpModulus<
                                  T1,
                                  typename std::decay_t<REP2>::value_type>>(
    stc::ScalarRep<T1>(a_scalar),
    std::forward<REP2>(a_rep2));
}

//--Operator between left vector representation and right scalar

template <typename REP1, typename T2,
          std::enable_if_t<
            std::is_base_of<stc::StcOp, std::decay_t<REP1>>::value
            &&
            !std::is_base_of<stc::StcOp, std::decay_t<T2>>::value
            &&
            std::decay_t<REP1>::is_vector::value
            &&
            !std::is_class<std::decay_t<T2>>::value, int> = 0
          >
HOSTDEVICE constexpr auto
operator%(REP1&& a_rep1, const T2& a_scalar) noexcept
{
  return stc::VectorBinaryOpRep<REP1, stc::ScalarRep<T2>,
                                stc::BinaryOpModulus<
                                  typename std::decay_t<REP1>::value_type,
                                  T2>>(
    std::forward<REP1>(a_rep1),
    stc::ScalarRep<T2>(a_scalar));
}

/*--------------------------------------------------------------------*
 * And
 *--------------------------------------------------------------------*/

template <typename T1, typename T2>
struct BinaryOpAnd
{
  using value_type = decltype(T1() & T2());
  HOSTDEVICE constexpr auto eval(const T1& a_x1, const T2& a_x2) const noexcept
    {
      return a_x1 & a_x2;
    }
};

//--Operator between two vector representations

template <typename REP1, typename REP2,
          std::enable_if_t<
            std::is_base_of<stc::StcOp, std::decay_t<REP1>>::value
            &&
            std::is_base_of<stc::StcOp, std::decay_t<REP2>>::value
            &&
            std::decay_t<REP1>::is_vector::value
            &&
            std::decay_t<REP2>::is_vector::value, int> = 0
          >
HOSTDEVICE constexpr auto
operator&(REP1&& a_rep1, REP2&& a_rep2) noexcept
{
  return stc::VectorBinaryOpRep<REP1, REP2,
                                stc::BinaryOpAnd<
                                  typename std::decay_t<REP1>::value_type,
                                  typename std::decay_t<REP2>::value_type>>(
    std::forward<REP1>(a_rep1),
    std::forward<REP2>(a_rep2));
}

//--Operator between left scalar and right vector representation

template <typename T1, typename REP2,
          std::enable_if_t<
            !std::is_base_of<stc::StcOp, std::decay_t<T1>>::value
            &&
            std::is_base_of<stc::StcOp, std::decay_t<REP2>>::value
            &&
            !std::is_class<std::decay_t<T1>>::value
            &&
            std::decay_t<REP2>::is_vector::value, int> = 0
          >
HOSTDEVICE constexpr auto
operator&(const T1& a_scalar, REP2&& a_rep2) noexcept
{
  return stc::VectorBinaryOpRep<stc::ScalarRep<T1>, REP2,
                                stc::BinaryOpAnd<
                                  T1,
                                  typename std::decay_t<REP2>::value_type>>(
    stc::ScalarRep<T1>(a_scalar),
    std::forward<REP2>(a_rep2));
}

//--Operator between left vector representation and right scalar

template <typename REP1, typename T2,
          std::enable_if_t<
            std::is_base_of<stc::StcOp, std::decay_t<REP1>>::value
            &&
            !std::is_base_of<stc::StcOp, std::decay_t<T2>>::value
            &&
            std::decay_t<REP1>::is_vector::value
            &&
            !std::is_class<std::decay_t<T2>>::value, int> = 0
          >
HOSTDEVICE constexpr auto
operator&(REP1&& a_rep1, const T2& a_scalar) noexcept
{
  return stc::VectorBinaryOpRep<REP1, stc::ScalarRep<T2>,
                                stc::BinaryOpAnd<
                                  typename std::decay_t<REP1>::value_type,
                                  T2>>(
    std::forward<REP1>(a_rep1),
    stc::ScalarRep<T2>(a_scalar));
}

/*--------------------------------------------------------------------*
 * Or
 *--------------------------------------------------------------------*/

template <typename T1, typename T2>
struct BinaryOpOr
{
  using value_type = decltype(T1() | T2());
  HOSTDEVICE constexpr auto eval(const T1& a_x1, const T2& a_x2) const noexcept
    {
      return a_x1 | a_x2;
    }
};

//--Operator between two vector representations

template <typename REP1, typename REP2,
          std::enable_if_t<
            std::is_base_of<stc::StcOp, std::decay_t<REP1>>::value
            &&
            std::is_base_of<stc::StcOp, std::decay_t<REP2>>::value
            &&
            std::decay_t<REP1>::is_vector::value
            &&
            std::decay_t<REP2>::is_vector::value, int> = 0
          >
HOSTDEVICE constexpr auto
operator|(REP1&& a_rep1, REP2&& a_rep2) noexcept
{
  return stc::VectorBinaryOpRep<REP1, REP2,
                                stc::BinaryOpOr<
                                  typename std::decay_t<REP1>::value_type,
                                  typename std::decay_t<REP2>::value_type>>(
    std::forward<REP1>(a_rep1),
    std::forward<REP2>(a_rep2));
}

//--Operator between left scalar and right vector representation

template <typename T1, typename REP2,
          std::enable_if_t<
            !std::is_base_of<stc::StcOp, std::decay_t<T1>>::value
            &&
            std::is_base_of<stc::StcOp, std::decay_t<REP2>>::value
            &&
            !std::is_class<std::decay_t<T1>>::value
            &&
            std::decay_t<REP2>::is_vector::value, int> = 0
          >
HOSTDEVICE constexpr auto
operator|(const T1& a_scalar, REP2&& a_rep2) noexcept
{
  return stc::VectorBinaryOpRep<stc::ScalarRep<T1>, REP2,
                                stc::BinaryOpOr<
                                  T1,
                                  typename std::decay_t<REP2>::value_type>>(
    stc::ScalarRep<T1>(a_scalar),
    std::forward<REP2>(a_rep2));
}

//--Operator between left vector representation and right scalar

template <typename REP1, typename T2,
          std::enable_if_t<
            std::is_base_of<stc::StcOp, std::decay_t<REP1>>::value
            &&
            !std::is_base_of<stc::StcOp, std::decay_t<T2>>::value
            &&
            std::decay_t<REP1>::is_vector::value
            &&
            !std::is_class<std::decay_t<T2>>::value, int> = 0
          >
HOSTDEVICE constexpr auto
operator|(REP1&& a_rep1, const T2& a_scalar) noexcept
{
  return stc::VectorBinaryOpRep<REP1, stc::ScalarRep<T2>,
                                stc::BinaryOpOr<
                                  typename std::decay_t<REP1>::value_type,
                                  T2>>(
    std::forward<REP1>(a_rep1),
    stc::ScalarRep<T2>(a_scalar));
}

/*--------------------------------------------------------------------*
 * XOr
 *--------------------------------------------------------------------*/

template <typename T1, typename T2>
struct BinaryOpXOr
{
  using value_type = decltype(T1() | T2());
  HOSTDEVICE constexpr auto eval(const T1& a_x1, const T2& a_x2) const noexcept
    {
      return a_x1 ^ a_x2;
    }
};

//--Operator between two vector representations

template <typename REP1, typename REP2,
          std::enable_if_t<
            std::is_base_of<stc::StcOp, std::decay_t<REP1>>::value
            &&
            std::is_base_of<stc::StcOp, std::decay_t<REP2>>::value
            &&
            std::decay_t<REP1>::is_vector::value
            &&
            std::decay_t<REP2>::is_vector::value, int> = 0
          >
HOSTDEVICE constexpr auto
operator^(REP1&& a_rep1, REP2&& a_rep2) noexcept
{
  return stc::VectorBinaryOpRep<REP1, REP2,
                                stc::BinaryOpXOr<
                                  typename std::decay_t<REP1>::value_type,
                                  typename std::decay_t<REP2>::value_type>>(
    std::forward<REP1>(a_rep1),
    std::forward<REP2>(a_rep2));
}

//--Operator between left scalar and right vector representation

template <typename T1, typename REP2,
          std::enable_if_t<
            !std::is_base_of<stc::StcOp, std::decay_t<T1>>::value
            &&
            std::is_base_of<stc::StcOp, std::decay_t<REP2>>::value
            &&
            !std::is_class<std::decay_t<T1>>::value
            &&
            std::decay_t<REP2>::is_vector::value, int> = 0
          >
HOSTDEVICE constexpr auto
operator^(const T1& a_scalar, REP2&& a_rep2) noexcept
{
  return stc::VectorBinaryOpRep<stc::ScalarRep<T1>, REP2,
                                stc::BinaryOpXOr<
                                  T1,
                                  typename std::decay_t<REP2>::value_type>>(
    stc::ScalarRep<T1>(a_scalar),
    std::forward<REP2>(a_rep2));
}

//--Operator between left vector representation and right scalar

template <typename REP1, typename T2,
          std::enable_if_t<
            std::is_base_of<stc::StcOp, std::decay_t<REP1>>::value
            &&
            !std::is_base_of<stc::StcOp, std::decay_t<T2>>::value
            &&
            std::decay_t<REP1>::is_vector::value
            &&
            !std::is_class<std::decay_t<T2>>::value, int> = 0
          >
HOSTDEVICE constexpr auto
operator^(REP1&& a_rep1, const T2& a_scalar) noexcept
{
  return stc::VectorBinaryOpRep<REP1, stc::ScalarRep<T2>,
                                stc::BinaryOpXOr<
                                  typename std::decay_t<REP1>::value_type,
                                  T2>>(
    std::forward<REP1>(a_rep1),
    stc::ScalarRep<T2>(a_scalar));
}

/*--------------------------------------------------------------------*
 * Complement
 *--------------------------------------------------------------------*/

template <typename T>
struct UnaryOpComplement
{
  using value_type = T;
  HOSTDEVICE constexpr auto eval(const T& a_x) const noexcept
    {
      return ~a_x;
    }
};

//--Operator

template <typename REP,
          std::enable_if_t<std::decay_t<REP>::is_vector::value, int> = 0
          >
HOSTDEVICE constexpr auto
operator~(REP&& a_rep) noexcept
{
  return stc::VectorUnaryOpRep<REP,
                               stc::UnaryOpComplement<
                                 typename std::decay_t<REP>::value_type>>(
    std::forward<REP>(a_rep));
}

/*--------------------------------------------------------------------*
 * Min
 *--------------------------------------------------------------------*/

template <typename T1, typename T2>
struct BinaryOpMin
{
  using value_type = decltype(T1() + T2());
  HOSTDEVICE constexpr auto eval(const T1& a_x1, const T2& a_x2) const noexcept
    {
#ifdef __CUDA_ARCH__
      return ::min(a_x1, a_x2);
#else
      return std::min(a_x1, a_x2);
#endif
    }
};

//--Operator between two vector representations

template <typename REP1, typename REP2,
          std::enable_if_t<
            std::is_base_of<stc::StcOp, std::decay_t<REP1>>::value
            &&
            std::is_base_of<stc::StcOp, std::decay_t<REP2>>::value
            &&
            std::decay_t<REP1>::is_vector::value
            &&
            std::decay_t<REP2>::is_vector::value, int>
          >
HOSTDEVICE constexpr auto
min(REP1&& a_rep1, REP2&& a_rep2) noexcept
{
  return stc::VectorBinaryOpRep<REP1, REP2,
                                stc::BinaryOpMin<
                                  typename std::decay_t<REP1>::value_type,
                                  typename std::decay_t<REP2>::value_type>>(
    std::forward<REP1>(a_rep1),
    std::forward<REP2>(a_rep2));
}

//--Operator between left scalar and right vector representation

template <typename T1, typename REP2,
          std::enable_if_t<
            !std::is_base_of<stc::StcOp, std::decay_t<T1>>::value
            &&
            std::is_base_of<stc::StcOp, std::decay_t<REP2>>::value
            &&
            !std::is_class<std::decay_t<T1>>::value
            &&
            std::decay_t<REP2>::is_vector::value, int> = 0
          >
HOSTDEVICE constexpr auto
min(const T1& a_scalar, REP2&& a_rep2) noexcept
{
  return stc::VectorBinaryOpRep<stc::ScalarRep<T1>, REP2,
                                stc::BinaryOpMin<
                                  T1,
                                  typename std::decay_t<REP2>::value_type>>(
    stc::ScalarRep<T1>(a_scalar),
    std::forward<REP2>(a_rep2));
}

//--Operator between left vector representation and right scalar

template <typename REP1, typename T2,
          std::enable_if_t<
            std::is_base_of<stc::StcOp, std::decay_t<REP1>>::value
            &&
            !std::is_base_of<stc::StcOp, std::decay_t<T2>>::value
            &&
            std::decay_t<REP1>::is_vector::value
            &&
            !std::is_class<std::decay_t<T2>>::value, int>
          >
HOSTDEVICE constexpr auto
min(REP1&& a_rep1, const T2& a_scalar) noexcept
{
  return stc::VectorBinaryOpRep<REP1, stc::ScalarRep<T2>,
                                stc::BinaryOpMin<
                                  typename std::decay_t<REP1>::value_type,
                                  T2>>(
    std::forward<REP1>(a_rep1),
    stc::ScalarRep<T2>(a_scalar));
}

/*--------------------------------------------------------------------*
 * Max
 *--------------------------------------------------------------------*/

template <typename T1, typename T2>
struct BinaryOpMax
{
  using value_type = decltype(T1() + T2());
  HOSTDEVICE constexpr auto eval(const T1& a_x1, const T2& a_x2) const noexcept
    {
#ifdef __CUDA_ARCH__
      return ::max(a_x1, a_x2);
#else
      return std::max(a_x1, a_x2);
#endif
    }
};

//--Operator between two vector representations

template <typename REP1, typename REP2,
          std::enable_if_t<
            std::is_base_of<stc::StcOp, std::decay_t<REP1>>::value
            &&
            std::is_base_of<stc::StcOp, std::decay_t<REP2>>::value
            &&
            std::decay_t<REP1>::is_vector::value
            &&
            std::decay_t<REP2>::is_vector::value, int>
          >
HOSTDEVICE constexpr auto
max(REP1&& a_rep1, REP2&& a_rep2) noexcept
{
  return stc::VectorBinaryOpRep<REP1, REP2,
                                stc::BinaryOpMax<
                                  typename std::decay_t<REP1>::value_type,
                                  typename std::decay_t<REP2>::value_type>>(
    std::forward<REP1>(a_rep1),
    std::forward<REP2>(a_rep2));
}

//--Operator between left scalar and right vector representation

template <typename T1, typename REP2,
          std::enable_if_t<
            !std::is_base_of<stc::StcOp, std::decay_t<T1>>::value
            &&
            std::is_base_of<stc::StcOp, std::decay_t<REP2>>::value
            &&
            !std::is_class<std::decay_t<T1>>::value
            &&
            std::decay_t<REP2>::is_vector::value, int> = 0
          >
HOSTDEVICE constexpr auto
max(const T1& a_scalar, REP2&& a_rep2) noexcept
{
  return stc::VectorBinaryOpRep<stc::ScalarRep<T1>, REP2,
                                stc::BinaryOpMax<
                                  T1,
                                  typename std::decay_t<REP2>::value_type>>(
    stc::ScalarRep<T1>(a_scalar),
    std::forward<REP2>(a_rep2));
}

//--Operator between left vector representation and right scalar

template <typename REP1, typename T2,
          std::enable_if_t<
            std::is_base_of<stc::StcOp, std::decay_t<REP1>>::value
            &&
            !std::is_base_of<stc::StcOp, std::decay_t<T2>>::value
            &&
            std::decay_t<REP1>::is_vector::value
            &&
            !std::is_class<std::decay_t<T2>>::value, int>
          >
HOSTDEVICE constexpr auto
max(REP1&& a_rep1, const T2& a_scalar) noexcept
{
  return stc::VectorBinaryOpRep<REP1, stc::ScalarRep<T2>,
                                stc::BinaryOpMax<
                                  typename std::decay_t<REP1>::value_type,
                                  T2>>(
    std::forward<REP1>(a_rep1),
    stc::ScalarRep<T2>(a_scalar));
}

/*
 * Beyond this point, we do not test that operands derive from StcOp.  If
 * any other class define is_vector, those test will have to be added to
 * enable_if
 */

/*
 * Comment on reduction operations.  These operations return a single value
 * and therefore must evaluate immediately.  They can be identified by some
 * application of *eachElement<N>(f).  To support lambdas, 'f' is accessed
 * with f(idx).  For unary reductions and some binary reductions, we either
 * directly pass the rep, or construct a new rep and pass this as 'f'.  See
 * 'sum', 'norm1', and 'dot' for examples.  If constructing a new rep, operands
 * are cast to const references since the operator must evalute immediately.
 * For some binary reductions, such as '==' below, there is no rep to reuse.
 * Instead, we use custom function objects that just take references.
 */

/*--------------------------------------------------------------------*
 * Equivalence
 *--------------------------------------------------------------------*/

template <typename T1, typename T2>
struct BinaryOpEquivalence
{
  using value_type = bool;
  HOSTDEVICE constexpr auto eval(const T1& a_x1, const T2& a_x2) const noexcept
    {
      return a_x1 == a_x2;
    }
};

//--Equivalence between two vector representations that immediately evaluates
//--and returns

template <typename REP1, typename REP2,
          std::enable_if_t<
            std::is_base_of<stc::StcOp, std::decay_t<REP1>>::value
            &&
            std::is_base_of<stc::StcOp, std::decay_t<REP2>>::value
            &&
            std::decay_t<REP1>::is_vector::value
            &&
            std::decay_t<REP2>::is_vector::value
            &&
            !std::decay_t<REP1>::is_string::value
            &&
            !std::decay_t<REP2>::is_string::value, int> = 0
          >
HOSTDEVICE constexpr bool
operator==(REP1&& a_rep1, REP2&& a_rep2) noexcept
{
  using T1 = typename std::decay_t<REP1>::value_type;
  using T2 = typename std::decay_t<REP2>::value_type;
  constexpr stc::array_size_type N = std::decay_t<REP1>::size();
  static_assert(std::is_same<T1, T2>::value && N == std::decay_t<REP2>::size(),
                "Invalid comparison between different vector types");
  auto f = stc::VectorBinaryOpRep<REP1, REP2,
                                  stc::BinaryOpEquivalence<T1, T2>>(
    std::forward<REP1>(a_rep1),
    std::forward<REP2>(a_rep2));
  return stc::andEachElement<N>(f);
}

//--Equivalence between left scalar and right vector representations that
//--immediately evaluates and returns

template <typename T1, typename REP2,
          std::enable_if_t<
            !std::is_base_of<stc::StcOp, std::decay_t<T1>>::value
            &&
            std::is_base_of<stc::StcOp, std::decay_t<REP2>>::value
            &&
            !std::is_class<std::decay_t<T1>>::value
            &&
            std::decay_t<REP2>::is_vector::value
            &&
            !std::decay_t<REP2>::is_string::value, int> = 0
          >
HOSTDEVICE constexpr bool
operator==(const T1& a_scalar, REP2&& a_rep2) noexcept
{
  using T2 = typename std::decay_t<REP2>::value_type;
  constexpr stc::array_size_type N = std::decay_t<REP2>::size();
  static_assert(std::is_same<T1, T2>::value,
                "Invalid comparison between different vector types");
  auto f = stc::VectorBinaryOpRep<stc::ScalarRep<T1>, REP2,
                                  stc::BinaryOpEquivalence<T1, T2>>(
    stc::ScalarRep<T1>(a_scalar),
    std::forward<REP2>(a_rep2));
  return stc::andEachElement<N>(f);
}

//--Equivalence between left vector and right scalar representations that
//--immediately evaluates and returns

template <typename REP1, typename T2,
          std::enable_if_t<
            std::is_base_of<stc::StcOp, std::decay_t<REP1>>::value
            &&
            !std::is_base_of<stc::StcOp, std::decay_t<T2>>::value
            &&
            std::decay_t<REP1>::is_vector::value
            &&
            !std::decay_t<REP1>::is_string::value
            &&
            !std::is_class<std::decay_t<T2>>::value, int> = 0
          >
HOSTDEVICE constexpr bool
operator==(REP1&& a_rep1, const T2& a_scalar) noexcept
{
  using T1 = typename std::decay_t<REP1>::value_type;
  constexpr stc::array_size_type N = std::decay_t<REP1>::size();
  static_assert(std::is_same<T1, T2>::value,
                "Invalid comparison between different vector types");
  auto f = stc::VectorBinaryOpRep<REP1, stc::ScalarRep<T2>,
                                  stc::BinaryOpEquivalence<T1, T2>>(
    std::forward<REP1>(a_rep1),
    stc::ScalarRep<T2>(a_scalar));
  return stc::andEachElement<N>(f);
}

//--Equivalence between two vector representations that returns a representation

template <typename REP1, typename REP2,
          std::enable_if_t<
            std::is_base_of<stc::StcOp, std::decay_t<REP1>>::value
            &&
            std::is_base_of<stc::StcOp, std::decay_t<REP2>>::value
            &&
            std::decay_t<REP1>::is_vector::value
            &&
            std::decay_t<REP2>::is_vector::value, int> = 0
          >
HOSTDEVICE constexpr auto
eq(REP1&& a_rep1, REP2&& a_rep2) noexcept
{
  return stc::VectorBinaryOpRep<REP1, REP2,
                                stc::BinaryOpEquivalence<
                                  typename std::decay_t<REP1>::value_type,
                                  typename std::decay_t<REP2>::value_type>>(
    std::forward<REP1>(a_rep1),
    std::forward<REP2>(a_rep2));
}

//--Equivalence between left scalar and right vector representations that
//--returns a representation

template <typename T1, typename REP2,
          std::enable_if_t<
            !std::is_base_of<stc::StcOp, std::decay_t<T1>>::value
            &&
            std::is_base_of<stc::StcOp, std::decay_t<REP2>>::value
            &&
            !std::is_class<std::decay_t<T1>>::value
            &&
            std::decay_t<REP2>::is_vector::value, int> = 0
          >
HOSTDEVICE constexpr auto
eq(const T1& a_scalar, REP2&& a_rep2) noexcept
{
  return stc::VectorBinaryOpRep<stc::ScalarRep<T1>, REP2,
                                stc::BinaryOpEquivalence<
                                  T1,
                                  typename std::decay_t<REP2>::value_type>>(
    stc::ScalarRep<T1>(a_scalar),
    std::forward<REP2>(a_rep2));
}

//--Equivalence between left vector and right scalar representations that
//--returns a representation

template <typename REP1, typename T2,
          std::enable_if_t<
            std::is_base_of<stc::StcOp, std::decay_t<REP1>>::value
            &&
            !std::is_base_of<stc::StcOp, std::decay_t<T2>>::value
            &&
            std::decay_t<REP1>::is_vector::value
            &&
            !std::is_class<std::decay_t<T2>>::value, int> = 0
          >
HOSTDEVICE constexpr auto
eq(REP1&& a_rep1, const T2& a_scalar) noexcept
{
  return stc::VectorBinaryOpRep<REP1, stc::ScalarRep<T2>,
                                stc::BinaryOpEquivalence<
                                  typename std::decay_t<REP1>::value_type,
                                  T2>>(
    std::forward<REP1>(a_rep1),
    stc::ScalarRep<T2>(a_scalar));
}

/*--------------------------------------------------------------------*
 * Non-equivalence
 *--------------------------------------------------------------------*/

//--Equivalence between two vector representations that immediately evaluates
//--and returns

// The operands must have the same type, size, and values for true equivalence

template <typename REP1, typename REP2,
          std::enable_if_t<
            std::is_base_of<stc::StcOp, std::decay_t<REP1>>::value
            &&
            std::is_base_of<stc::StcOp, std::decay_t<REP2>>::value
            &&
            std::decay_t<REP1>::is_vector::value
            &&
            std::decay_t<REP2>::is_vector::value
            &&
            !std::decay_t<REP1>::is_string::value
            &&
            !std::decay_t<REP2>::is_string::value, int> = 0
          >
HOSTDEVICE constexpr bool
operator!=(REP1&& a_rep1, REP2&& a_rep2) noexcept
{
  using T1 = typename std::decay_t<REP1>::value_type;
  using T2 = typename std::decay_t<REP2>::value_type;
  constexpr stc::array_size_type N = std::decay_t<REP1>::size();
  if (!std::is_same<T1, T2>::value ||
      std::decay_t<REP1>::size() != std::decay_t<REP2>::size())
    {
      return true;
    }
  const auto f = stc::VectorBinaryOpRep<REP1, REP2,
                                        stc::BinaryOpEquivalence<T1, T2>>(
    std::forward<REP1>(a_rep1),
    std::forward<REP2>(a_rep2));
  return !stc::andEachElement<N>(f);
}

//--Non-equivalence between left scalar and right vector representations that
//--immediately evaluates and returns

// The operands must have the same type

template <typename T1, typename REP2,
          std::enable_if_t<
            !std::is_base_of<stc::StcOp, std::decay_t<T1>>::value
            &&
            std::is_base_of<stc::StcOp, std::decay_t<REP2>>::value
            &&
            !std::is_class<std::decay_t<T1>>::value
            &&
            std::decay_t<REP2>::is_vector::value
            &&
            !std::decay_t<REP2>::is_string::value, int> = 0
          >
HOSTDEVICE constexpr bool
operator!=(const T1& a_scalar, REP2&& a_rep2) noexcept
{
  using T2 = typename std::decay_t<REP2>::value_type;
  constexpr stc::array_size_type N = std::decay_t<REP2>::size();
  if (!std::is_same<T1, T2>::value)
    {
      return true;
    }
  auto f = stc::VectorBinaryOpRep<stc::ScalarRep<T1>, REP2,
                                  stc::BinaryOpEquivalence<T1, T2>>(
    stc::ScalarRep<T1>(a_scalar),
    std::forward<REP2>(a_rep2));
  return !stc::andEachElement<N>(f);
}

//--Non-equivalence between left vector and right scalar representations that
//--immediately evaluates and returns

// The operands must have the same type

template <typename REP1, typename T2,
          std::enable_if_t<
            std::is_base_of<stc::StcOp, std::decay_t<REP1>>::value
            &&
            !std::is_base_of<stc::StcOp, std::decay_t<T2>>::value
            &&
            std::decay_t<REP1>::is_vector::value
            &&
            !std::decay_t<REP1>::is_string::value
            &&
            !std::is_class<std::decay_t<T2>>::value, int> = 0
          >
HOSTDEVICE constexpr bool
operator!=(REP1&& a_rep1, const T2& a_scalar) noexcept
{
  using T1 = typename std::decay_t<REP1>::value_type;
  constexpr stc::array_size_type N = std::decay_t<REP1>::size();
  if (!std::is_same<T1, T2>::value)
    {
      return true;
    }
  auto f = stc::VectorBinaryOpRep<REP1, stc::ScalarRep<T2>,
                                  stc::BinaryOpEquivalence<T1, T2>>(
    std::forward<REP1>(a_rep1),
    stc::ScalarRep<T2>(a_scalar));
  return !stc::andEachElement<N>(f);
}

/*--------------------------------------------------------------------*
 * Less than
 *--------------------------------------------------------------------*/

template <typename T1, typename T2>
struct BinaryOpLessThan
{
  using value_type = bool;
  HOSTDEVICE constexpr auto eval(const T1& a_x1, const T2& a_x2) const noexcept
    {
      return a_x1 < a_x2;
    }
};

//--Less than between two vector representations that immediately evaluates and
//--returns

template <typename REP1, typename REP2,
          std::enable_if_t<
            std::is_base_of<stc::StcOp, std::decay_t<REP1>>::value
            &&
            std::is_base_of<stc::StcOp, std::decay_t<REP2>>::value
            &&
            std::decay_t<REP1>::is_vector::value
            &&
            std::decay_t<REP2>::is_vector::value, int> = 0
          >
HOSTDEVICE constexpr bool
operator<(REP1&& a_rep1, REP2&& a_rep2) noexcept
{
  using T1 = typename std::decay_t<REP1>::value_type;
  using T2 = typename std::decay_t<REP2>::value_type;
  constexpr stc::array_size_type N = std::decay_t<REP1>::size();
  static_assert(std::is_same<T1, T2>::value && N == std::decay_t<REP2>::size(),
                "Invalid comparison between different vector types");
  auto f = stc::VectorBinaryOpRep<REP1, REP2,
                                  stc::BinaryOpLessThan<T1, T2>>(
    std::forward<REP1>(a_rep1),
    std::forward<REP2>(a_rep2));
  return stc::andEachElement<N>(f);
}

//--Less than between left scalar and right vector representations that
//--immediately evaluates and returns

template <typename T1, typename REP2,
          std::enable_if_t<
            !std::is_base_of<stc::StcOp, std::decay_t<T1>>::value
            &&
            std::is_base_of<stc::StcOp, std::decay_t<REP2>>::value
            &&
            !std::is_class<std::decay_t<T1>>::value
            &&
            std::decay_t<REP2>::is_vector::value, int> = 0
          >
HOSTDEVICE constexpr bool
operator<(const T1& a_scalar, REP2&& a_rep2) noexcept
{
  using T2 = typename std::decay_t<REP2>::value_type;
  constexpr stc::array_size_type N = std::decay_t<REP2>::size();
  static_assert(std::is_same<T1, T2>::value,
                "Invalid comparison between different vector types");
  auto f = stc::VectorBinaryOpRep<stc::ScalarRep<T1>, REP2,
                                  stc::BinaryOpLessThan<T1, T2>>(
    stc::ScalarRep<T1>(a_scalar),
    std::forward<REP2>(a_rep2));
  return stc::andEachElement<N>(f);
}

//--Less than between left vector and right scalar representations that
//--immediately evaluates and returns

template <typename REP1, typename T2,
          std::enable_if_t<
            std::is_base_of<stc::StcOp, std::decay_t<REP1>>::value
            &&
            !std::is_base_of<stc::StcOp, std::decay_t<T2>>::value
            &&
            std::decay_t<REP1>::is_vector::value
            &&
            !std::is_class<std::decay_t<T2>>::value, int> = 0
          >
HOSTDEVICE constexpr bool
operator<(REP1&& a_rep1, const T2& a_scalar) noexcept
{
  using T1 = typename std::decay_t<REP1>::value_type;
  constexpr stc::array_size_type N = std::decay_t<REP1>::size();
  static_assert(std::is_same<T1, T2>::value,
                "Invalid comparison between different vector types");
  auto f = stc::VectorBinaryOpRep<REP1, stc::ScalarRep<T2>,
                                  stc::BinaryOpLessThan<T1, T2>>(
    std::forward<REP1>(a_rep1),
    stc::ScalarRep<T2>(a_scalar));
  return stc::andEachElement<N>(f);
}

//--Less than between two vector representations that returns a representation

template <typename REP1, typename REP2,
          std::enable_if_t<
            std::is_base_of<stc::StcOp, std::decay_t<REP1>>::value
            &&
            std::is_base_of<stc::StcOp, std::decay_t<REP2>>::value
            &&
            std::decay_t<REP1>::is_vector::value
            &&
            std::decay_t<REP2>::is_vector::value, int> = 0
          >
HOSTDEVICE constexpr auto
lt(REP1&& a_rep1, REP2&& a_rep2) noexcept
{
  return stc::VectorBinaryOpRep<REP1, REP2,
                                stc::BinaryOpLessThan<
                                  typename std::decay_t<REP1>::value_type,
                                  typename std::decay_t<REP2>::value_type>>(
    std::forward<REP1>(a_rep1),
    std::forward<REP2>(a_rep2));
}

//--Less than between left scalar and right vector representations that returns
//--a representation

template <typename T1, typename REP2,
          std::enable_if_t<
            !std::is_base_of<stc::StcOp, std::decay_t<T1>>::value
            &&
            std::is_base_of<stc::StcOp, std::decay_t<REP2>>::value
            &&
            !std::is_class<std::decay_t<T1>>::value
            &&
            std::decay_t<REP2>::is_vector::value, int> = 0
          >
HOSTDEVICE constexpr auto
lt(const T1& a_scalar, REP2&& a_rep2) noexcept
{
  return stc::VectorBinaryOpRep<stc::ScalarRep<T1>, REP2,
                                stc::BinaryOpLessThan<
                                  T1,
                                  typename std::decay_t<REP2>::value_type>>(
    stc::ScalarRep<T1>(a_scalar),
    std::forward<REP2>(a_rep2));
}

//--Less than between left vector and right scalar representations that returns
//--a representation

template <typename REP1, typename T2,
          std::enable_if_t<
            std::is_base_of<stc::StcOp, std::decay_t<REP1>>::value
            &&
            !std::is_base_of<stc::StcOp, std::decay_t<T2>>::value
            &&
            std::decay_t<REP1>::is_vector::value
            &&
            !std::is_class<std::decay_t<T2>>::value, int> = 0
          >
HOSTDEVICE constexpr auto
lt(REP1&& a_rep1, const T2& a_scalar) noexcept
{
  return stc::VectorBinaryOpRep<REP1, stc::ScalarRep<T2>,
                                stc::BinaryOpLessThan<
                                  typename std::decay_t<REP1>::value_type,
                                  T2>>(
    std::forward<REP1>(a_rep1),
    stc::ScalarRep<T2>(a_scalar));
}

/*--------------------------------------------------------------------*
 * Less than or equal
 * Note: you cannot use !> because for vectors, > does not imply !<.
 *       In other words, both can be false.
 *--------------------------------------------------------------------*/

template <typename T1, typename T2>
struct BinaryOpLessThanEqual
{
  using value_type = bool;
  HOSTDEVICE constexpr auto eval(const T1& a_x1, const T2& a_x2) const noexcept
    {
      return a_x1 <= a_x2;
    }
};

//--Less than equal between two vector representations that immediately
//--evaluates and returns

template <typename REP1, typename REP2,
          std::enable_if_t<
            std::is_base_of<stc::StcOp, std::decay_t<REP1>>::value
            &&
            std::is_base_of<stc::StcOp, std::decay_t<REP2>>::value
            &&
            std::decay_t<REP1>::is_vector::value
            &&
            std::decay_t<REP2>::is_vector::value, int> = 0
          >
HOSTDEVICE constexpr bool
operator<=(REP1&& a_rep1, REP2&& a_rep2) noexcept
{
  using T1 = typename std::decay_t<REP1>::value_type;
  using T2 = typename std::decay_t<REP2>::value_type;
  constexpr stc::array_size_type N = std::decay_t<REP1>::size();
  static_assert(std::is_same<T1, T2>::value && N == std::decay_t<REP2>::size(),
                "Invalid comparison between different vector types");
  auto f = stc::VectorBinaryOpRep<REP1, REP2,
                                  stc::BinaryOpLessThanEqual<T1, T2>>(
    std::forward<REP1>(a_rep1),
    std::forward<REP2>(a_rep2));
  return stc::andEachElement<N>(f);
}

//--Less than equal between left scalar and right vector representations that
//--immediately evaluates and returns

template <typename T1, typename REP2,
          std::enable_if_t<
            !std::is_base_of<stc::StcOp, std::decay_t<T1>>::value
            &&
            std::is_base_of<stc::StcOp, std::decay_t<REP2>>::value
            &&
            !std::is_class<std::decay_t<T1>>::value
            &&
            std::decay_t<REP2>::is_vector::value, int> = 0
          >
HOSTDEVICE constexpr bool
operator<=(const T1& a_scalar, REP2&& a_rep2) noexcept
{
  using T2 = typename std::decay_t<REP2>::value_type;
  constexpr stc::array_size_type N = std::decay_t<REP2>::size();
  static_assert(std::is_same<T1, T2>::value,
                "Invalid comparison between different vector types");
  auto f = stc::VectorBinaryOpRep<stc::ScalarRep<T1>, REP2,
                                  stc::BinaryOpLessThanEqual<T1, T2>>(
    stc::ScalarRep<T1>(a_scalar),
    std::forward<REP2>(a_rep2));
  return stc::andEachElement<N>(f);
}

//--Less than equal between left vector and right scalar representations that
//--immediately evaluates and returns

template <typename REP1, typename T2,
          std::enable_if_t<
            std::is_base_of<stc::StcOp, std::decay_t<REP1>>::value
            &&
            !std::is_base_of<stc::StcOp, std::decay_t<T2>>::value
            &&
            std::decay_t<REP1>::is_vector::value
            &&
            !std::is_class<std::decay_t<T2>>::value, int> = 0
          >
HOSTDEVICE constexpr bool
operator<=(REP1&& a_rep1, const T2& a_scalar) noexcept
{
  using T1 = typename std::decay_t<REP1>::value_type;
  constexpr stc::array_size_type N = std::decay_t<REP1>::size();
  static_assert(std::is_same<T1, T2>::value,
                "Invalid comparison between different vector types");
  auto f = stc::VectorBinaryOpRep<REP1, stc::ScalarRep<T2>,
                                  stc::BinaryOpLessThanEqual<T1, T2>>(
    std::forward<REP1>(a_rep1),
    stc::ScalarRep<T2>(a_scalar));
  return stc::andEachElement<N>(f);
}

//--Less than equal between two vector representations that returns a
//--representation

template <typename REP1, typename REP2,
          std::enable_if_t<
            std::is_base_of<stc::StcOp, std::decay_t<REP1>>::value
            &&
            std::is_base_of<stc::StcOp, std::decay_t<REP2>>::value
            &&
            std::decay_t<REP1>::is_vector::value
            &&
            std::decay_t<REP2>::is_vector::value, int> = 0
          >
HOSTDEVICE constexpr auto
le(REP1&& a_rep1, REP2&& a_rep2) noexcept
{
  return stc::VectorBinaryOpRep<REP1, REP2,
                                stc::BinaryOpLessThanEqual<
                                  typename std::decay_t<REP1>::value_type,
                                  typename std::decay_t<REP2>::value_type>>(
    std::forward<REP1>(a_rep1),
    std::forward<REP2>(a_rep2));
}

//--Less than equal between left scalar and right vector representations that
//--returns a representation

template <typename T1, typename REP2,
          std::enable_if_t<
            !std::is_base_of<stc::StcOp, std::decay_t<T1>>::value
            &&
            std::is_base_of<stc::StcOp, std::decay_t<REP2>>::value
            &&
            !std::is_class<std::decay_t<T1>>::value
            &&
            std::decay_t<REP2>::is_vector::value, int> = 0
          >
HOSTDEVICE constexpr auto
le(const T1& a_scalar, REP2&& a_rep2) noexcept
{
  return stc::VectorBinaryOpRep<stc::ScalarRep<T1>, REP2,
                                stc::BinaryOpLessThanEqual<
                                  T1,
                                  typename std::decay_t<REP2>::value_type>>(
    stc::ScalarRep<T1>(a_scalar),
    std::forward<REP2>(a_rep2));
}

//--Less than equal between left vector and right scalar representations that
//--returns a representation

template <typename REP1, typename T2,
          std::enable_if_t<
            std::is_base_of<stc::StcOp, std::decay_t<REP1>>::value
            &&
            !std::is_base_of<stc::StcOp, std::decay_t<T2>>::value
            &&
            std::decay_t<REP1>::is_vector::value
            &&
            !std::is_class<std::decay_t<T2>>::value, int> = 0
          >
HOSTDEVICE constexpr auto
le(REP1&& a_rep1, const T2& a_scalar) noexcept
{
  return stc::VectorBinaryOpRep<REP1, stc::ScalarRep<T2>,
                                stc::BinaryOpLessThanEqual<
                                  typename std::decay_t<REP1>::value_type,
                                  T2>>(
    std::forward<REP1>(a_rep1),
    stc::ScalarRep<T2>(a_scalar));
}

/*--------------------------------------------------------------------*
 * Greater than
 *--------------------------------------------------------------------*/

//--Greater than between two vector representations that immediately evaluates
//--and returns

template <typename REP1, typename REP2,
          std::enable_if_t<
            std::is_base_of<stc::StcOp, std::decay_t<REP1>>::value
            &&
            std::is_base_of<stc::StcOp, std::decay_t<REP2>>::value
            &&
            std::decay_t<REP1>::is_vector::value
            &&
            std::decay_t<REP2>::is_vector::value, int> = 0
          >
HOSTDEVICE constexpr bool
operator>(REP1&& a_rep1, REP2&& a_rep2) noexcept
{
  using T1 = typename std::decay_t<REP1>::value_type;
  using T2 = typename std::decay_t<REP2>::value_type;
  constexpr stc::array_size_type N = std::decay_t<REP1>::size();
  static_assert(std::is_same<T1, T2>::value && N == std::decay_t<REP2>::size(),
                "Invalid comparison between different vector types");
  auto f = stc::VectorBinaryOpRep<REP2, REP1,
                                  stc::BinaryOpLessThan<T2, T1>>(
    std::forward<REP2>(a_rep2),
    std::forward<REP1>(a_rep1));
  return stc::andEachElement<N>(f);
}

//--Greater than between left scalar and right vector representations that
//--immediately evaluates and returns

template <typename T1, typename REP2,
          std::enable_if_t<
            !std::is_base_of<stc::StcOp, std::decay_t<T1>>::value
            &&
            std::is_base_of<stc::StcOp, std::decay_t<REP2>>::value
            &&
            !std::is_class<std::decay_t<T1>>::value
            &&
            std::decay_t<REP2>::is_vector::value, int> = 0
          >
HOSTDEVICE constexpr bool
operator>(const T1& a_scalar, REP2&& a_rep2) noexcept
{
  using T2 = typename std::decay_t<REP2>::value_type;
  constexpr stc::array_size_type N = std::decay_t<REP2>::size();
  static_assert(std::is_same<T1, T2>::value,
                "Invalid comparison between different vector types");
  auto f = stc::VectorBinaryOpRep<REP2, stc::ScalarRep<T1>,
                                  stc::BinaryOpLessThan<T2, T1>>(
    std::forward<REP2>(a_rep2),
    stc::ScalarRep<T1>(a_scalar));
  return stc::andEachElement<N>(f);
}

//--Greater than between left vector and right scalar representations that
//--immediately evaluates and returns

template <typename REP1, typename T2,
          std::enable_if_t<
            std::is_base_of<stc::StcOp, std::decay_t<REP1>>::value
            &&
            !std::is_base_of<stc::StcOp, std::decay_t<T2>>::value
            &&
            std::decay_t<REP1>::is_vector::value
            &&
            !std::is_class<std::decay_t<T2>>::value, int> = 0
          >
HOSTDEVICE constexpr bool
operator>(REP1&& a_rep1, const T2& a_scalar) noexcept
{
  using T1 = typename std::decay_t<REP1>::value_type;
  constexpr stc::array_size_type N = std::decay_t<REP1>::size();
  static_assert(std::is_same<T1, T2>::value,
                "Invalid comparison between different vector types");
  auto f = stc::VectorBinaryOpRep<stc::ScalarRep<T2>, REP1,
                                  stc::BinaryOpLessThan<T2, T1>>(
    stc::ScalarRep<T2>(a_scalar),
    std::forward<REP1>(a_rep1));
  return stc::andEachElement<N>(f);
}

/*--------------------------------------------------------------------*
 * Greater than or equal
 * Note: you cannot use !< because for vectors, < does not imply !>.
 *       In other words, both can be false.
 *--------------------------------------------------------------------*/

//--Greater than equal immediately evaluates and returns

template <typename REP1, typename REP2,
          std::enable_if_t<
            std::is_base_of<stc::StcOp, std::decay_t<REP1>>::value
            &&
            std::is_base_of<stc::StcOp, std::decay_t<REP2>>::value
            &&
            std::decay_t<REP1>::is_vector::value
            &&
            std::decay_t<REP2>::is_vector::value, int> = 0
          >
HOSTDEVICE constexpr bool
operator>=(REP1&& a_rep1, REP2&& a_rep2) noexcept
{
  using T1 = typename std::decay_t<REP1>::value_type;
  using T2 = typename std::decay_t<REP2>::value_type;
  constexpr stc::array_size_type N = std::decay_t<REP1>::size();
  static_assert(std::is_same<T1, T2>::value && N == std::decay_t<REP2>::size(),
                "Invalid comparison between different vector types");
  auto f = stc::VectorBinaryOpRep<REP2, REP1,
                                  stc::BinaryOpLessThanEqual<T2, T1>>(
    std::forward<REP2>(a_rep2),
    std::forward<REP1>(a_rep1));
  return stc::andEachElement<N>(f);
}

//--Greater than equal between left scalar and right vector representations that
//--immediately evaluates and returns

template <typename T1, typename REP2,
          std::enable_if_t<
            !std::is_base_of<stc::StcOp, std::decay_t<T1>>::value
            &&
            std::is_base_of<stc::StcOp, std::decay_t<REP2>>::value
            &&
            !std::is_class<std::decay_t<T1>>::value
            &&
            std::decay_t<REP2>::is_vector::value, int> = 0
          >
HOSTDEVICE constexpr bool
operator>=(const T1& a_scalar, REP2&& a_rep2) noexcept
{
  using T2 = typename std::decay_t<REP2>::value_type;
  constexpr stc::array_size_type N = std::decay_t<REP2>::size();
  static_assert(std::is_same<T1, T2>::value,
                "Invalid comparison between different vector types");
  auto f = stc::VectorBinaryOpRep<REP2, stc::ScalarRep<T1>,
                                  stc::BinaryOpLessThanEqual<T2, T1>>(
    std::forward<REP2>(a_rep2),
    stc::ScalarRep<T1>(a_scalar));
  return stc::andEachElement<N>(f);
}

//--Greater than equal between left vector and right scalar representations that
//--immediately evaluates and returns

template <typename REP1, typename T2,
          std::enable_if_t<
            std::is_base_of<stc::StcOp, std::decay_t<REP1>>::value
            &&
            !std::is_base_of<stc::StcOp, std::decay_t<T2>>::value
            &&
            std::decay_t<REP1>::is_vector::value
            &&
            !std::is_class<std::decay_t<T2>>::value, int> = 0
          >
HOSTDEVICE constexpr bool
operator>=(REP1&& a_rep1, const T2& a_scalar) noexcept
{
  using T1 = typename std::decay_t<REP1>::value_type;
  constexpr stc::array_size_type N = std::decay_t<REP1>::size();
  static_assert(std::is_same<T1, T2>::value,
                "Invalid comparison between different vector types");
  auto f = stc::VectorBinaryOpRep<stc::ScalarRep<T2>, REP1,
                                  stc::BinaryOpLessThanEqual<T2, T1>>(
    stc::ScalarRep<T2>(a_scalar),
    std::forward<REP1>(a_rep1));
  return stc::andEachElement<N>(f);
}

/*--------------------------------------------------------------------*
 * Lexical less than.  First < has lower index than first >
 *--------------------------------------------------------------------*/

// Use immediately
template <typename REP1, typename REP2>
struct BinaryOpLessThan_E
{
  using value_type = bool;
  HOSTDEVICE constexpr BinaryOpLessThan_E(const REP1& a_rep1,
                                          const REP2& a_rep2) noexcept
    :
    m_rep1(a_rep1),
    m_rep2(a_rep2)
    { }
  HOSTDEVICE constexpr bool operator()(const array_size_type a_idx)
    const noexcept
    {
      return m_rep1[a_idx] < m_rep2[a_idx];
    }
  const REP1& m_rep1;
  const REP2& m_rep2;
};

//--Lexical less than immediately evaluates and returns

template <typename REP1, typename REP2,
          std::enable_if_t<
            std::decay_t<REP1>::is_vector::value
            &&
            std::decay_t<REP2>::is_vector::value, int>
          >
HOSTDEVICE constexpr bool
lexLT(REP1&& a_rep1, REP2&& a_rep2) noexcept
{
  using T1 = typename std::decay_t<REP1>::value_type;
  using T2 = typename std::decay_t<REP2>::value_type;
  static_assert(std::is_same<T1, T2>::value &&
                std::decay_t<REP1>::size() == std::decay_t<REP2>::size(),
                "Invalid lexLT between different vector types");
  constexpr array_size_type N = std::decay_t<REP1>::size();
  using FLT = stc::BinaryOpLessThan_E<std::decay_t<REP1>, std::decay_t<REP2>>;
  auto flt = FLT(a_rep1, a_rep2);
  stc::array_size_type idxlt = stc::lowestTrueElementIndex<N>(flt);
  using FGT = stc::BinaryOpLessThan_E<std::decay_t<REP2>, std::decay_t<REP1>>;
  auto fgt = FGT(a_rep2, a_rep1);
  stc::array_size_type idxgt = stc::lowestTrueElementIndex<N>(fgt);
  return idxlt < idxgt;
}

/*--------------------------------------------------------------------*
 * Lexical greater than.  First > has lower index than first <
 *--------------------------------------------------------------------*/

//--Lexical greater than immediately evaluates and returns

template <typename REP1, typename REP2,
          std::enable_if_t<
            std::decay_t<REP1>::is_vector::value
            &&
            std::decay_t<REP2>::is_vector::value, int>
          >
HOSTDEVICE constexpr bool
lexGT(REP1&& a_rep1, REP2&& a_rep2) noexcept
{
  using T1 = typename std::decay_t<REP1>::value_type;
  using T2 = typename std::decay_t<REP2>::value_type;
  static_assert(std::is_same<T1, T2>::value &&
                std::decay_t<REP1>::size() == std::decay_t<REP2>::size(),
                "Invalid lexLT between different vector types");
  constexpr array_size_type N = std::decay_t<REP1>::size();
  using FLT = stc::BinaryOpLessThan_E<std::decay_t<REP1>, std::decay_t<REP2>>;
  auto flt = FLT(a_rep1, a_rep2);
  stc::array_size_type idxlt = stc::lowestTrueElementIndex<N>(flt);
  using FGT = stc::BinaryOpLessThan_E<std::decay_t<REP2>, std::decay_t<REP1>>;
  auto fgt = FGT(a_rep2, a_rep1);
  stc::array_size_type idxgt = stc::lowestTrueElementIndex<N>(fgt);
  return idxgt < idxlt;
}

/*--------------------------------------------------------------------*
 * Z-curve ordering less than (Morton ordering)
 *--------------------------------------------------------------------*/

//--Morton immediately evaluates and returns

template <typename REP1, typename REP2,
          std::enable_if_t<
            std::decay_t<REP1>::is_vector::value
            &&
            std::decay_t<REP2>::is_vector::value, int> = 0
          >
HOSTDEVICE constexpr bool
morton(REP1&& a_rep1, REP2&& a_rep2) noexcept
{
  constexpr array_size_type I = std::decay_t<REP1>::size();
  static_assert(I == std::decay_t<REP2>::size(),
                "Vectors must have same size for Morton comparison");
  using T = typename std::decay_t<REP1>::value_type;
  static_assert(std::is_same<T, typename std::decay_t<REP2>::value_type>::value,
                "Vectors must have same type for Morton comparison");
  // Evaluate here since elements are repeatedly accessed in indexMaxMSBElem
  const Vector<T, I> xorVal = a_rep1 ^ a_rep2;
  const array_size_type idxMaxMSB =
    indexMaxMSBElem<const Vector<T, I>&>(xorVal);
  return a_rep1[idxMaxMSB] < a_rep2[idxMaxMSB];
}

/*--------------------------------------------------------------------*
 * Position of minimum element (lowest position of equal elements)
 *--------------------------------------------------------------------*/

//--indexMinElem immediately evaluates and returns

// Use immediately
template <typename REP>
struct MinElem_E
{
  using T = typename std::decay_t<REP>::value_type;
  HOSTDEVICE constexpr MinElem_E(const REP& a_rep) noexcept
    :
    m_rep(a_rep)
    { }
  HOSTDEVICE constexpr array_size_type operator()(
    const array_size_type a_idx,
    const array_size_type a_idxOther) const noexcept
    {
      return (m_rep(a_idx) < m_rep(a_idxOther)) ? a_idx : a_idxOther;
    }
  const REP& m_rep;
};

template <typename REP,
          std::enable_if_t<std::decay_t<REP>::is_vector::value, int>
          >
HOSTDEVICE constexpr array_size_type
indexMinElem(REP&& a_rep) noexcept
{
  constexpr array_size_type N = std::decay_t<REP>::size();
  MinElem_E<std::decay_t<REP>> f(a_rep);
  return selectElementIndex<N>(f);
}

/*--------------------------------------------------------------------*
 * Value of minimum element
 *--------------------------------------------------------------------*/

template <typename REP,
          std::enable_if_t<std::decay_t<REP>::is_vector::value, int>
          >
HOSTDEVICE constexpr auto
minElem(REP&& a_rep) noexcept
{
  return a_rep[indexMinElem(a_rep)];
}

/*--------------------------------------------------------------------*
 * Position of maximum element (lowest position of equal elements)
 *--------------------------------------------------------------------*/

//--indexMaxElem immediately evaluates and returns

// Use immediately
template <typename REP>
struct MaxElem_E
{
  using T = typename std::decay_t<REP>::value_type;
  HOSTDEVICE constexpr MaxElem_E(const REP& a_rep) noexcept
    :
    m_rep(a_rep)
    { }
  HOSTDEVICE constexpr array_size_type operator()(
    const array_size_type a_idx,
    const array_size_type a_idxOther) const noexcept
    {
      return (m_rep(a_idx) > m_rep(a_idxOther)) ? a_idx : a_idxOther;
    }
  const REP& m_rep;
};

template <typename REP,
          std::enable_if_t<std::decay_t<REP>::is_vector::value, int>
          >
HOSTDEVICE constexpr array_size_type
indexMaxElem(REP&& a_rep) noexcept
{
  constexpr array_size_type N = std::decay_t<REP>::size();
  MaxElem_E<std::decay_t<REP>> f(a_rep);
  return selectElementIndex<N>(f);
}

/*--------------------------------------------------------------------*
 * Value of maximum element
 *--------------------------------------------------------------------*/

template <typename REP,
          std::enable_if_t<std::decay_t<REP>::is_vector::value, int>
          >
HOSTDEVICE constexpr auto
maxElem(REP&& a_rep) noexcept
{
  return a_rep[indexMaxElem(a_rep)];
}

/*--------------------------------------------------------------------*
 * Position of element with largest MSB (highest position of equal
 * elements)
 *--------------------------------------------------------------------*/

//--indexMaxMSBElem immediately evaluates and returns

// Use immediately
template <typename REP>
struct MaxMSBElem_E
{
  using T = typename std::decay_t<REP>::value_type;
  HOSTDEVICE constexpr MaxMSBElem_E(const REP& a_rep) noexcept
    :
    m_rep(a_rep)
    { }
  HOSTDEVICE constexpr array_size_type operator()(
    const array_size_type a_idx,
    const array_size_type a_idxOther) const noexcept
    {
      const T x = m_rep(a_idx);
      const T y = m_rep(a_idxOther);
      return (x < y && x < (x^y)) ? a_idxOther : a_idx;
    }
  const REP& m_rep;
};

template <typename REP,
          std::enable_if_t<std::decay_t<REP>::is_vector::value, int>
          >
HOSTDEVICE constexpr array_size_type
indexMaxMSBElem(REP&& a_rep) noexcept
{
  constexpr array_size_type N = std::decay_t<REP>::size();
  MaxMSBElem_E<std::decay_t<REP>> f(a_rep);
  return selectElementIndex<N>(f);
}

/*--------------------------------------------------------------------*
 * Log 2 of an integer.  This is optimized for constexpr output and
 * assumes the input is a power of 2 (otherwise this gives the floor
 * of log2).  It is required that a_x > 0 but T can be signed or
 * unsigned.  It is mostly intended for use with coarsen_log2.  The
 * function coarsen is better for run-time values since the log2 is
 * achieved by a popcnt which probably maps to a single instruction on
 * most architectures.
 *--------------------------------------------------------------------*/

template <typename T>
struct UnaryOpCLog2
{
  using value_type = T;
  HOSTDEVICE constexpr auto eval(const T& a_x) const noexcept
    {
      return _clog2(a_x);
    }
};

//--Operator for log2

template <typename REP,
          std::enable_if_t<std::decay_t<REP>::is_vector::value, int> = 0
          >
HOSTDEVICE constexpr auto
clog2(REP&& a_rep1) noexcept
{
  static_assert(std::is_integral<typename std::decay_t<REP>::value_type>
                ::value,
                "clog2 can only be applied to integral vector types");
  return stc::VectorUnaryOpRep<REP,
                               stc::UnaryOpCLog2<
                                 typename std::decay_t<REP>::value_type>>(
    std::forward<REP>(a_rep1));
}

/*--------------------------------------------------------------------*
 * Log 2 of an integer.  This is optimized for runtime evaluation and
 * assumes the input is a power of 2 (otherwise this gives the floor
 * of log2).  It is required that a_x > 0 but T can be signed or
 * unsigned.
 *--------------------------------------------------------------------*/

template <typename T>
struct UnaryOpLog2
{
  using value_type = T;
  HOSTDEVICE constexpr auto eval(const T& a_x) const noexcept
    {
      return popcnt(a_x - (T)1);
    }
};

//--Operator for log2

template <typename REP,
          std::enable_if_t<std::decay_t<REP>::is_vector::value, int> = 0
          >
HOSTDEVICE constexpr auto
log2(REP&& a_rep1) noexcept
{
  static_assert(std::is_integral<typename std::decay_t<REP>::value_type>
                ::value,
                "log2 can only be applied to integral vector types");
  return stc::VectorUnaryOpRep<REP,
                               stc::UnaryOpLog2<
                                 typename std::decay_t<REP>::value_type>>(
    std::forward<REP>(a_rep1));
}

/*--------------------------------------------------------------------*
 * Compile-time coarsen (only for integral vector types and factor
 * must be a power of 2).
 *--------------------------------------------------------------------*/

//--Used for two vectors

template <typename T1, typename T2>
struct BinaryOpCoarsenLog2
{
  using value_type = T1;
  HOSTDEVICE constexpr auto eval(const T1& a_x1, const T2& a_x2) const noexcept
    {
      return coarsen_log2(a_x1, a_x2);
    }
};

//--Operator between two vector representations (coarsen x by y)

template <typename REP1, typename REP2,
          std::enable_if_t<
            std::decay_t<REP1>::is_vector::value
            &&
            std::decay_t<REP2>::is_vector::value, int> = 0
          >
HOSTDEVICE constexpr auto
ccoarsen(REP1&& a_rep1, REP2&& a_rep2) noexcept
{
  static_assert(std::is_integral<typename std::decay_t<REP1>::value_type>
                ::value &&
                std::is_integral<typename std::decay_t<REP2>::value_type>
                ::value,
                "Coarsen can only be applied to integral vector types");
  return stc::VectorBinaryOpRep<REP1,
                                stc::VectorUnaryOpRep<
                                  REP2,
                                  stc::UnaryOpCLog2<
                                    typename std::decay_t<REP2>::value_type>>,
                                stc::BinaryOpCoarsenLog2<
                                  typename std::decay_t<REP1>::value_type,
                                  typename std::decay_t<REP2>::value_type>>(
    std::forward<REP1>(a_rep1),
    clog2(std::forward<REP2>(a_rep2)));
}

//--Operator for coarsening by a scalar factor (coarsen x by y)

template <typename REP,
          std::enable_if_t<std::decay_t<REP>::is_vector::value, int> = 0
          >
HOSTDEVICE constexpr auto
ccoarsen(REP&& a_rep1, unsigned long long a_factor) noexcept
{
  static_assert(std::is_integral<typename std::decay_t<REP>::value_type>
                ::value,
                "Coarsen can only be applied to integral vector types");
  return stc::VectorBinaryOpRep<REP, stc::ScalarRep<unsigned>,
                                stc::BinaryOpCoarsenLog2<
                                  typename std::decay_t<REP>::value_type,
                                  unsigned>>(
    std::forward<REP>(a_rep1),
    static_cast<unsigned>(_clog2(a_factor)));
}

//--Operator between two vector representations (coarsen x by 2^y)

template <typename REP1, typename REP2,
          std::enable_if_t<
            std::decay_t<REP1>::is_vector::value
            &&
            std::decay_t<REP2>::is_vector::value, int> = 0
          >
HOSTDEVICE constexpr auto
ccoarsen_log2(REP1&& a_rep1, REP2&& a_rep2) noexcept
{
  static_assert(std::is_integral<typename std::decay_t<REP1>::value_type>
                ::value &&
                std::is_integral<typename std::decay_t<REP2>::value_type>
                ::value,
                "Coarsen can only be applied to integral vector types");
  return stc::VectorBinaryOpRep<REP1,
                                REP2,
                                stc::BinaryOpCoarsenLog2<
                                  typename std::decay_t<REP1>::value_type,
                                  typename std::decay_t<REP2>::value_type>>(
    std::forward<REP1>(a_rep1),
    std::forward<REP2>(a_rep2));
}

//--Operator for coarsening by a scalar factor (coarsen x by 2^y)

template <typename REP,
          std::enable_if_t<std::decay_t<REP>::is_vector::value, int> = 0
          >
HOSTDEVICE constexpr auto
ccoarsen_log2(REP&& a_rep1, unsigned a_factor) noexcept
{
  static_assert(std::is_integral<typename std::decay_t<REP>::value_type>
                ::value,
                "Coarsen can only be applied to integral vector types");
  return stc::VectorBinaryOpRep<REP, stc::ScalarRep<unsigned>,
                                stc::BinaryOpCoarsenLog2<
                                  typename std::decay_t<REP>::value_type,
                                  unsigned>>(
    std::forward<REP>(a_rep1),
    a_factor);
}

/*--------------------------------------------------------------------*
 * Coarsen (only for integral vector types and factor must be a power
 * of 2).  Requires twos complement and that >> is arithmetic.
 * These are not constexpr because of popcnt.
 *--------------------------------------------------------------------*/

template <typename T1, typename T2>
struct BinaryOpCoarsen
{
  using value_type = T1;
  HOSTDEVICE auto eval(const T1& a_x1, const T2& a_x2) const noexcept
    {
      return coarsen_sys(a_x1, a_x2);
      // assert(isPower2(a_x2));
      // return a_x1 >> stc::popcnt(a_x2 - 1);
    }
};

//--Operator between two vector representations

template <typename REP1, typename REP2,
          std::enable_if_t<
            std::decay_t<REP1>::is_vector::value
            &&
            std::decay_t<REP2>::is_vector::value, int>
          >
HOSTDEVICE auto
coarsen(REP1&& a_rep1, REP2&& a_rep2) noexcept
{
  static_assert(std::is_integral<typename std::decay_t<REP1>::value_type>
                ::value &&
                std::is_integral<typename std::decay_t<REP2>::value_type>
                ::value,
                "Coarsen can only be applied to integral vector types");
  return stc::VectorBinaryOpRep<REP1, REP2,
                                stc::BinaryOpCoarsen<
                                  typename std::decay_t<REP1>::value_type,
                                  typename std::decay_t<REP2>::value_type>>(
    std::forward<REP1>(a_rep1),
    std::forward<REP2>(a_rep2));
}

template <typename T>
struct UnaryOpCoarsenLog2Cache
{
  using value_type = T;
  HOSTDEVICE UnaryOpCoarsenLog2Cache(const unsigned a_pc)
    :
    m_pc(a_pc)
    { }
  const unsigned m_pc;
  HOSTDEVICE auto eval(const T& a_x) const noexcept
    {
      return coarsen_log2(a_x, m_pc);
    }
};

//--Operator for coarsening by a scalar factor

template <typename REP,
          std::enable_if_t<std::decay_t<REP>::is_vector::value, int>
          >
HOSTDEVICE auto
coarsen(REP&& a_rep1, unsigned long long a_factor) noexcept
{
  static_assert(std::is_integral<typename std::decay_t<REP>::value_type>
                ::value,
                "Coarsen can only be applied to integral vector types");
  assert(stc::isPower2(a_factor));  // Do we want this?
  return stc::VectorUnaryOpRep<REP,
                               stc::UnaryOpCoarsenLog2Cache<
                                 typename std::decay_t<REP>::value_type>>(
    std::forward<REP>(a_rep1), stc::popcnt(a_factor - 1));
}

/*--------------------------------------------------------------------*
 * Scale (same as operator* but type restricted)
 *--------------------------------------------------------------------*/

template <typename REP,
          std::enable_if_t<std::decay_t<REP>::is_vector::value, int> = 0
          >
HOSTDEVICE constexpr auto
scale(REP&&                                        a_rep,
      const typename std::decay_t<REP>::value_type a_scalar) noexcept
{
  return std::forward<REP>(a_rep)*a_scalar;
}

/*--------------------------------------------------------------------*
 * Reflection
 *--------------------------------------------------------------------*/

template <typename REP>
struct UnaryROpReflect
{
  using value_type = typename REP::value_type;
  constexpr UnaryROpReflect(const value_type      a_refIdx,
                            const array_size_type a_dir)
    :
    m_refIdx(a_refIdx),
    m_dir(a_dir)
    { }
  const value_type m_refIdx;
  const array_size_type m_dir;
  HOSTDEVICE constexpr auto eval(const REP&            a_rep,
                                 const array_size_type a_idx) const noexcept
    {
      return (a_idx == m_dir) ? -a_rep[a_idx] + 2*m_refIdx : a_rep[a_idx];
    }
};

//--Operator for reflecting about plane defined by a_refIdx with normal in the
//--direction of a_dir

template <typename REP,
          std::enable_if_t<std::decay_t<REP>::is_vector::value, int> = 0
          >
HOSTDEVICE constexpr auto
reflect(REP&&                                        a_rep,
        const typename std::decay_t<REP>::value_type a_refIdx,
        const array_size_type                        a_dir) noexcept
{
  using F = stc::UnaryROpReflect<std::decay_t<REP>>;
  return stc::VectorUnaryROpRep<REP, F>(std::forward<REP>(a_rep),
                                        a_refIdx,
                                        a_dir);
}

/*--------------------------------------------------------------------*
 * Shift (same as operator+ but type restricted)
 *--------------------------------------------------------------------*/

template <typename REP>
struct UnaryROpShiftDir
{
  using value_type = typename REP::value_type;
  HOSTDEVICE constexpr UnaryROpShiftDir(const array_size_type a_dir,
                                        const value_type      a_shift)
    :
    m_shift(a_shift),
    m_dir(a_dir)
    { }
  const value_type m_shift;
  const array_size_type m_dir;
  HOSTDEVICE constexpr auto eval(const REP&            a_rep,
                                 const array_size_type a_idx) const noexcept
    {
      return (a_idx == m_dir) ? a_rep[a_idx] + m_shift : a_rep[a_idx];
    }
};

//--Operator for shifting in a single direction

template <typename REP1, typename T2,
          std::enable_if_t<
            std::decay_t<REP1>::is_vector::value
            &&
            !std::is_class<std::decay_t<T2>>::value, int> = 0
          >
HOSTDEVICE constexpr auto
shift(REP1&& a_rep, const array_size_type a_dir, const T2 a_shift) noexcept
{
  using F = stc::UnaryROpShiftDir<std::decay_t<REP1>>;
  return stc::VectorUnaryROpRep<REP1, F>(
    std::forward<REP1>(a_rep),
    a_dir,
    static_cast<typename std::decay_t<REP1>::value_type>(a_shift));
}

// Same as BinaryOpAdd but using T1 as return type
template <typename T1, typename T2>
struct BinaryOpShiftVector
{
  using value_type = T1;
  HOSTDEVICE constexpr T1 eval(const T1& a_x1, const T2& a_x2) const noexcept
    {
      return a_x1 + a_x2;
    }
};

//--Operator for shifting by a vector

template <typename REP1, typename REP2,
          std::enable_if_t<
            std::decay_t<REP2>::is_vector::value
            &&
            std::decay_t<REP2>::is_vector::value, int> = 0
          >
HOSTDEVICE constexpr auto
shift(REP1&& a_rep1, REP2&& a_rep2shiftBy) noexcept
{
  return stc::VectorBinaryOpRep<REP1, REP2,
                                stc::BinaryOpShiftVector<
                                  typename std::decay_t<REP1>::value_type,
                                  typename std::decay_t<REP2>::value_type>>(
    std::forward<REP1>(a_rep1),
    std::forward<REP2>(a_rep2shiftBy));
}

//--Operator for shifting by a scalar

template <typename REP1, typename T2,
          std::enable_if_t<
            std::decay_t<REP1>::is_vector::value
            &&
            !std::is_class<std::decay_t<T2>>::value, int> = 0
          >
HOSTDEVICE constexpr auto
diagShift(REP1&& a_rep1, const T2& a_shiftBy) noexcept
{
  return stc::VectorBinaryOpRep<REP1, stc::ScalarRep<T2>,
                                stc::BinaryOpShiftVector<
                                  typename std::decay_t<REP1>::value_type,
                                  T2>>(
    std::forward<REP1>(a_rep1),
    stc::ScalarRep<T2>(a_shiftBy));
}

/*--------------------------------------------------------------------*
 * Dot product
 *--------------------------------------------------------------------*/

//--Dot immediately evaluates and returns.  Sizes must be the same.

template <typename REP1, typename REP2,
          std::enable_if_t<
            std::is_base_of<stc::StcOp, std::decay_t<REP1>>::value
            &&
            std::is_base_of<stc::StcOp, std::decay_t<REP2>>::value, int>
          >
HOSTDEVICE constexpr auto
dot(REP1&& a_rep1, REP2&& a_rep2) noexcept
{
  constexpr array_size_type I = std::decay_t<REP1>::size();
  static_assert(I == std::decay_t<REP2>::size(),
                "Vectors must have same size for dot product");
  // Whatever REPs were, all we want going forward are constant references.
  // The cast is required to avoid possible storage in operator*
  return stc::sumEachElement<I>(static_cast<const REP1&>(a_rep1)*
                                static_cast<const REP2&>(a_rep2));
}

/*--------------------------------------------------------------------*
 * Cross product
 *--------------------------------------------------------------------*/

  /* We try to avoid evaluating an expression twice.  Each operand would
     have to be evalualted twice at each index to fully determine the cross
     product.  This is why temporaries are used instead of composing a new
     expression using ROp, which would be an alternative approach.
  */

template <typename REP1, typename REP2>
struct BinaryOpCross
{
  using value_type = decltype(typename REP1::value_type()*
                              typename REP2::value_type());
  HOSTDEVICE constexpr BinaryOpCross(const REP1& a_rep1,
                                     const REP2& a_rep2) noexcept
    :
    m_rep1(a_rep1),
    m_rep2(a_rep2)
    { }
  HOSTDEVICE constexpr auto operator()(
    const array_size_type a_idx) const noexcept
    {
      const array_size_type a_i1 = (a_idx + 1) % 3;
      const array_size_type a_i2 = (a_idx + 2) % 3;
      return m_rep1[a_i1]*m_rep2[a_i2] - m_rep1[a_i2]*m_rep2[a_i1];
    }
  const REP1& m_rep1;
  const REP2& m_rep2;
};

//--Reps are both VectorImpl

template <typename REP1, typename REP2>
HOSTDEVICE constexpr auto
crossImpl(REP1&& a_rep1, REP2&& a_rep2, std::true_type) noexcept
{
  constexpr array_size_type I = std::decay_t<REP1>::size();
  using T = typename std::decay_t<REP1>::value_type;
  using F = stc::BinaryOpCross<std::decay_t<REP1>, std::decay_t<REP2>>;
  auto f = F(a_rep1, a_rep2);
  return Vector<T, I>(make_array<T, I>(f));
}

//--At least one Rep is an expression

template <typename REP1, typename REP2>
HOSTDEVICE constexpr auto
crossImpl(REP1&& a_rep1, REP2&& a_rep2, std::false_type) noexcept
{
  constexpr array_size_type I = std::decay_t<REP1>::size();
  using T = typename std::decay_t<REP1>::value_type;
  Vector<T, I> vec1 = a_rep1;  // Evaluates a_rep1 if expression
  Vector<T, I> vec2 = a_rep2;  // Evaluates a_rep2 if expression
  using F = stc::BinaryOpCross<decltype(vec1), decltype(vec2)>;
  auto f = F(vec1, vec2);
  return Vector<T, I>(make_array<T, I>(f));
}

//--Cross immediately evaluates and returns.  Size must be 3.
//--The entry to cross simply does a dispatch to the implementations

template <typename REP1, typename REP2,
          std::enable_if_t<
            std::decay_t<REP1>::is_vector::value
            &&
            std::decay_t<REP2>::is_vector::value, int> = 0
          >
HOSTDEVICE constexpr auto
cross(REP1&& a_rep1, REP2&& a_rep2) noexcept
{
  constexpr array_size_type I = std::decay_t<REP1>::size();
  static_assert(I == 3,
                "Vectors must have size 3 for cross product");
  static_assert(I == std::decay_t<REP2>::size(),
                "Vectors must have same size for cross product");
  using T = typename std::decay_t<REP1>::value_type;
  static_assert(std::is_same<T, typename std::decay_t<REP2>::value_type>::value,
                "Vectors must have same type for cross product");
  return stc::crossImpl(
    std::forward<REP1>(a_rep1),
    std::forward<REP2>(a_rep2),
    std::integral_constant<
      bool,
      (
        std::is_base_of<stc::VectorImpl<T, I, DefaultImpl>,
                        std::decay_t<REP1>>::value
        ||
        std::is_base_of<stc::VectorImpl<T, I, stc::VectorConstAliasImpl>,
                        std::decay_t<REP1>>::value
        ||
        std::is_base_of<stc::VectorImpl<T, I, stc::VectorAliasImpl>,
                        std::decay_t<REP1>>::value
      )
      &&
      (
        std::is_base_of<stc::VectorImpl<T, I, DefaultImpl>,
                        std::decay_t<REP2>>::value
        ||
        std::is_base_of<stc::VectorImpl<T, I, stc::VectorConstAliasImpl>,
                        std::decay_t<REP2>>::value
        ||
        std::is_base_of<stc::VectorImpl<T, I, stc::VectorAliasImpl>,
                        std::decay_t<REP2>>::value
      )
    >{});
}

/*--------------------------------------------------------------------*
 * Negation
 *--------------------------------------------------------------------*/

template <typename T>
struct UnaryOpNegation
{
  using value_type = T;
  HOSTDEVICE constexpr auto eval(const T& a_x) const noexcept
    {
      return -a_x;
    }
};

//--Operator

template <typename REP,
          std::enable_if_t<std::decay_t<REP>::is_vector::value, int> = 0
          >
HOSTDEVICE constexpr auto
operator-(REP&& a_rep) noexcept
{
  return stc::VectorUnaryOpRep<REP,
                               stc::UnaryOpNegation<
                                 typename std::decay_t<REP>::value_type>>(
    std::forward<REP>(a_rep));
}

/*--------------------------------------------------------------------*
 * Component-wise absolute value
 *--------------------------------------------------------------------*/

template <typename T>
struct UnaryOpAbs
{
  using value_type = T;
  HOSTDEVICE constexpr auto eval(const T& a_x) const noexcept
    {
      // Note: clang++-6.0 does not support constexpr abs
      // return std::abs(a_x);
      return (a_x < static_cast<T>(0)) ? -a_x : a_x;
    }
};

//--Operator

template <typename REP,
          std::enable_if_t<std::decay_t<REP>::is_vector::value, int>
          >
HOSTDEVICE constexpr auto
abs(REP&& a_rep) noexcept
{
  return stc::VectorUnaryOpRep<REP,
                               stc::UnaryOpAbs<
                                 typename std::decay_t<REP>::value_type>>(
    std::forward<REP>(a_rep));
}

/*--------------------------------------------------------------------*
 * Sum of all components
 *--------------------------------------------------------------------*/

template <typename REP,
          std::enable_if_t<std::decay_t<REP>::is_vector::value, int>
          >
HOSTDEVICE constexpr auto
sum(REP&& a_rep) noexcept
{
  constexpr array_size_type I = std::decay_t<REP>::size();
  return stc::sumEachElement<I>(std::forward<REP>(a_rep));
}

/*--------------------------------------------------------------------*
 * Product of all components
 *--------------------------------------------------------------------*/

template <typename REP,
          std::enable_if_t<std::decay_t<REP>::is_vector::value, int>
          >
HOSTDEVICE constexpr auto
product(REP&& a_rep) noexcept
{
  constexpr array_size_type I = std::decay_t<REP>::size();
  return stc::prodEachElement<I>(std::forward<REP>(a_rep));
}

/*--------------------------------------------------------------------*
 * 1-norm of all components
 *--------------------------------------------------------------------*/

template <typename REP,
          std::enable_if_t<std::decay_t<REP>::is_vector::value, int>
          >
HOSTDEVICE constexpr auto
norm1(REP&& a_rep) noexcept
{
  constexpr array_size_type I = std::decay_t<REP>::size();
  // Whatever REP was, all we want going forward is a constant reference
  // The cast is required to avoid possible storage in abs
  return stc::sumEachElement<I>(stc::abs(static_cast<const REP&>(a_rep)));
}

/*--------------------------------------------------------------------*
 * Magnitude (Euclidean norm of the vector)
 * This is not constexpr because std::sqrt is not
 *--------------------------------------------------------------------*/

//--Mag immediately evaluates and returns.  Type of result is floating point.

template <typename REP,
          std::enable_if_t<
            std::is_base_of<stc::StcOp, std::decay_t<REP>>::value, int>
          >
HOSTDEVICE inline auto
mag(REP&& a_rep) noexcept
{
  constexpr array_size_type I = std::decay_t<REP>::size();
  using T = typename std::decay_t<REP>::value_type;
  using Tfp = std::conditional_t<std::is_floating_point<T>::value, T, Real>;
  auto f = [&](const array_size_type a_idx)
    {
      return std::pow((Tfp)a_rep[a_idx], 2);
    };
  return std::sqrt(stc::sumEachElement<I>(f));
}

//--This high-precision version first sorts the vector from small to high.
//--It is quite useless for small vectors but I am keeping it around to remember
//--the technique.

template <typename REP,
          std::enable_if_t<
            std::is_base_of<stc::StcOp, std::decay_t<REP>>::value, int> = 0
          >
inline auto
mag_mp(REP&& a_rep) noexcept
{
  constexpr array_size_type I = std::decay_t<REP>::size();
  using T = typename std::decay_t<REP>::value_type;
  using Tfp = std::conditional_t<std::is_floating_point<T>::value, T, Real>;
  stc::Vector<T, I> vec = a_rep;  // Evaluates a_rep if expression

  // Build an indirection array (0, 1, ..., I)
  DefaultImpl<array_size_type, I> indirect(
    stc::make_array<array_size_type, I>([](const array_size_type a_idx)
                                        { return a_idx; }));
  // Sort the indirection array
  std::sort(indirect.begin(), indirect.end(),
            [&](const array_size_type i, const array_size_type j)
            {
              return std::abs(vec[i]) < std::abs(vec[j]);
            });

  auto f = [&](const array_size_type a_idx)
    {
      return std::pow((Tfp)vec[indirect[a_idx]], 2);
    };
  return std::sqrt(stc::sumEachElement<I>(f));
}

/*--------------------------------------------------------------------*
 * Unit vector
 * This is not constexpr because std::sqrt (from mag) is not.
 * This is optimized by tag dispatch for Reps of type VectorImpl and
 * those that are expressions.  The latter are evaluated into a
 * temporary.
 *--------------------------------------------------------------------*/

//--Standard unit

//--Rep is a VectorImpl

template <typename REP>
HOSTDEVICE inline auto
unitImpl(REP&& a_rep, std::true_type)
{
  using T = typename std::decay_t<REP>::value_type;
  // mag returns a floating type
  T norm = static_cast<T>(mag(std::forward<REP>(a_rep)));
  return a_rep/norm;
}

//--Rep is an expression.

  /* We try to avoid evaluating an expression twice.  Since the expression
     must be evaluated to get the norm and the unit, we first store into a
     temporary
  */

template <typename REP>
HOSTDEVICE inline auto
unitImpl(REP&& a_rep, std::false_type)
{
  constexpr array_size_type I = std::decay_t<REP>::size();
  using T = typename std::decay_t<REP>::value_type;
  stc::Vector<T, I> vec = a_rep;       // Evaluates a_rep if expression
  T norm = static_cast<T>(mag(vec));    // mag returns a floating type
  return std::move(vec)/norm;           // Move into new rep
}

//--The entry to unit simply does a dispatch to the implementations

template <typename REP,
          std::enable_if_t<
            std::is_base_of<stc::StcOp, std::decay_t<REP>>::value, int> = 0
          >
HOSTDEVICE inline auto
unit(REP&& a_rep)
{
  constexpr array_size_type I = std::decay_t<REP>::size();
  using T = typename std::decay_t<REP>::value_type;
  return stc::unitImpl(
    std::forward<REP>(a_rep),
    std::integral_constant<
      bool,
      std::is_base_of<stc::VectorImpl<T, I, DefaultImpl>,
                      std::decay_t<REP>>::value
      ||
      std::is_base_of<stc::VectorImpl<T, I, stc::VectorConstAliasImpl>,
                      std::decay_t<REP>>::value
      ||
      std::is_base_of<stc::VectorImpl<T, I, stc::VectorAliasImpl>,
                      std::decay_t<REP>>::value
    >{});
}

/*--------------------------------------------------------------------*
 * Compile-time stride.  This works but is inefficient for run-time
 * values.
 *--------------------------------------------------------------------*/

template <typename REP>
struct UnaryROpStride
{
  using value_type = typename REP::value_type;
  struct CStrideAtIdx
  {
    HOSTDEVICE constexpr CStrideAtIdx(const REP&            a_rep,
                                      const array_size_type a_idx) noexcept
      :
      m_rep(a_rep),
      m_idx(static_cast<int>(a_idx) - 1)
      { }
    const REP& m_rep;
    const int m_idx;
    HOSTDEVICE constexpr value_type operator()(const array_size_type a_idx)
      const noexcept
      {
        return (static_cast<int>(a_idx) > m_idx) ?
          static_cast<value_type>(1) : m_rep[a_idx];
      }
  };
  HOSTDEVICE constexpr auto eval(const REP&            a_rep,
                                 const array_size_type a_idx) const noexcept
    {
      constexpr array_size_type I = REP::size();
      CStrideAtIdx f(a_rep, a_idx);
      return prodEachElement<I>(f);
    }
};

//--Operator for stride

template <typename REP,
          std::enable_if_t<std::decay_t<REP>::is_vector::value, int> = 0
          >
HOSTDEVICE constexpr auto
cstride(REP&& a_rep) noexcept
{
  using F = stc::UnaryROpStride<std::decay_t<REP>>;
  return stc::VectorUnaryROpRep<REP, F>(std::forward<REP>(a_rep));
}

/*--------------------------------------------------------------------*
 * Run-time stride
 *--------------------------------------------------------------------*/

template <typename REP>
struct Stride_E
{
  static constexpr array_size_type N = REP::size();
  using T = typename REP::value_type;
  HOSTDEVICE Stride_E(Vector<T, N>& a_array, const REP& a_rep) noexcept
    :
    m_array(a_array),
    m_rep(a_rep)
    { }
  HOSTDEVICE void operator()(const array_size_type a_idx) noexcept
    {
      const int idx = a_idx + 1;
      m_array[idx] = m_array[idx-1]*m_rep[idx-1];
    }
  Vector<T, N>& m_array;
  const REP&    m_rep;
};

//--Stride immediately evaluates and creates a new vector

template <typename REP,
          std::enable_if_t<
            std::is_base_of<stc::StcOp, std::decay_t<REP>>::value
            &&
            std::decay_t<REP>::is_vector::value, int> = 0
          >
HOSTDEVICE inline auto
stride(REP&& a_rep)
{
  constexpr array_size_type I = std::decay_t<REP>::size();
  using T = typename std::decay_t<REP>::value_type;
  stc::Vector<T, I> array;
  array[0] = static_cast<T>(1);
  Stride_E<std::decay_t<REP>> f(array, a_rep);
  forEachElement<I-1>(f);
  // return std::move(array);  // Move into new destination
  return array;  // Above would prohibit copy elision.  Compilers should either
                 // use elision or automatically select move here
}

/*--------------------------------------------------------------------*
 * Construct from linear index and log2 of dimensions
 *--------------------------------------------------------------------*/

template <typename REP>
struct FromLinearLog2_E
{
  static constexpr array_size_type N = REP::size();
  using T = typename REP::value_type;
  HOSTDEVICE FromLinearLog2_E(stc::Vector<T, N>& a_array,
                              const REP&         a_log2DimRep,
                              array_size_type&   a_lin) noexcept
    :
    m_array(a_array),
    m_log2DimRep(a_log2DimRep),
    m_lin(a_lin),
    m_shift(stc::sumEachElement<N-1>(a_log2DimRep))
    { }
  HOSTDEVICE void operator()(const array_size_type a_idx) noexcept
    {
      const int idx = a_idx + 1;
      m_array[idx] = m_lin >> m_shift;
      m_lin -= m_array[idx] << m_shift;
      m_shift -= m_log2DimRep[a_idx];
    }
  stc::Vector<T, N>& m_array;
  const REP& m_log2DimRep;
  array_size_type& m_lin;
  array_size_type m_shift;
};

//--From linear immediately evaluates and creates a new vector

template <typename REP,
          std::enable_if_t<
            std::is_base_of<stc::StcOp, std::decay_t<REP>>::value
            &&
            std::decay_t<REP>::is_vector::value, int> = 0
          >
HOSTDEVICE inline auto
fromLinearLog2Dim(array_size_type a_lin, const REP& a_log2DimRep)
{
  constexpr array_size_type I = std::decay_t<REP>::size();
  using T = typename std::decay_t<REP>::value_type;
  stc::Vector<T, I> array;
  FromLinearLog2_E<std::decay_t<REP>> f(array, a_log2DimRep, a_lin);
  forEachReverseElement<I-1>(f);
  array[0] = a_lin;
  return array;
}


/*--------------------------------------------------------------------*
 * Transpose
 *--------------------------------------------------------------------*/

/// Transpose represenation operator
/** This is called a wrap because its sole purpose is to deliver m_rep to
 *  some other operation.  It should never be saved on the LHS.
 */
template <typename REP>
class VectorTransposeOpWrap : private StcOp
{
public:
  using REP_type = REP;
  using T = typename std::decay_t<REP>::value_type;
  using is_scalar = std::false_type;
  using is_vector = std::false_type;
  using is_rowvec = std::true_type;
  using is_matrix = std::false_type;
  using is_string = std::false_type;
private:
  REP&& m_rep;
public:
  HOSTDEVICE constexpr VectorTransposeOpWrap(REP&& a_rep) noexcept
    :
    m_rep(std::forward<REP>(a_rep))
    { }
  VectorTransposeOpWrap(const VectorTransposeOpWrap&) = default;
  VectorTransposeOpWrap(VectorTransposeOpWrap&&) = default;
  /// Assignment constructor
  VectorTransposeOpWrap& operator=(const VectorTransposeOpWrap&) = delete;
  /// Move assignment constructor
  VectorTransposeOpWrap& operator=(VectorTransposeOpWrap&&) = delete;
  /// Size of the vector representation
  HOSTDEVICE static constexpr array_size_type size() noexcept
    {
      return std::decay_t<REP>::size();
    }
  /// Access
  HOSTDEVICE constexpr REP&& get() const noexcept
    {
      return std::forward<REP>(m_rep);
    }
};

/// Transpose Operator
template <typename REP,
          std::enable_if_t<std::decay_t<REP>::is_vector::value, int> = 0
          >
constexpr auto
operator^(REP&& a_rep, stc::TransposeType) noexcept
{
  return stc::VectorTransposeOpWrap<REP>(std::forward<REP>(a_rep));
}


/*******************************************************************************
 *
 * Class VectorImpl: external related functions
 *
 ******************************************************************************/

/*--------------------------------------------------------------------*/
/// Output the vector
/** \param[in] a_os     The output stream
 *  \param[in] a_vec    Vector to output
 *  \return             The output stream
 *//*-----------------------------------------------------------------*/

template <typename T,
          stc::array_size_type N,
          template<typename, stc::array_size_type> class Impl>
inline std::ostream&
operator<<(std::ostream&                      a_os,
           const stc::VectorImpl<T, N, Impl>& a_vec)
{
  a_os << '(';
  if (N > 0)
    {
      a_os << a_vec[0];
    }
  for (int i = 1; i != N; ++i)
    {
      a_os << ',' << a_vec[i];
    }
  a_os << ')';
  return a_os;
}

/*--------------------------------------------------------------------*/
//  Input the vector
/** \param[in] a_os     The input stream
 *  \param[in] a_vec    Vector to load
 *  \return             The input stream
 *//*-----------------------------------------------------------------*/

template <typename T,
          stc::array_size_type N,
          template<typename, stc::array_size_type> class Impl>
inline std::istream&
operator>>(std::istream&                a_is,
           stc::VectorImpl<T, N, Impl>& a_vec)
{
  constexpr size_t maxIgnore = std::numeric_limits<std::streamsize>::max();
  a_is >> std::ws;
  char c;
  a_is >> c;
  a_is.putback(c);
  if (c == '(')
    {
      a_is.ignore(maxIgnore, '(') >> a_vec[0];
      for (int i = 1; i != N; ++i)
        {
          a_is.ignore(maxIgnore, ',') >> a_vec[i];
        }
      a_is.ignore(maxIgnore, ')');
    }
  else if (c == '<')
    {
      a_is.ignore(maxIgnore, '<') >> a_vec[0];
      for (int i = 1; i != N; ++i)
        {
          a_is.ignore(maxIgnore, ',') >> a_vec[i];
        }
      a_is.ignore(maxIgnore, '>');
    }
  else
    {
      // "operator>>(istream&,stc::VectorImpl&): expected \'(\' or \'<\'"
      abort();
    }
  if (a_is.fail())
    abort();  // "operator>>(istream&,stc::VectorImpl&) failed";
  return a_is;
}

/*==============================================================================
 * Comparison function objects
 *============================================================================*/

/*--------------------------------------------------------------------*/
//  Lexicographic less-than function object
/** 
 *//*-----------------------------------------------------------------*/

template <typename T,
          stc::array_size_type N,
          template<typename, stc::array_size_type> class Impl = DefaultImpl>
struct CompareLT
{
  using Vect = stc::VectorImpl<T, N, Impl>;
  constexpr bool operator() (const Vect& a_x, const Vect& a_y) const noexcept
    {
      return stc::lexLT(a_x, a_y);
    }
};

/*--------------------------------------------------------------------*/
//  Equal function object
/** 
 *//*-----------------------------------------------------------------*/

template <typename T,
          stc::array_size_type N,
          template<typename, stc::array_size_type> class Impl = DefaultImpl>
struct CompareEQ
{
  using Vect = stc::VectorImpl<T, N, Impl>;
  constexpr bool operator() (const Vect& a_x, const Vect& a_y) const noexcept
    {
      return a_x == a_y;
    }
};

/*--------------------------------------------------------------------*/
//  Z-curve (Morton ordering) less-than function object
/** 
 *//*-----------------------------------------------------------------*/

template <typename T,
          stc::array_size_type N,
          template<typename, stc::array_size_type> class Impl = DefaultImpl>
struct MortonLT
{
  using Vect = stc::VectorImpl<T, N, Impl>;
  constexpr bool operator() (const Vect& a_x, const Vect& a_y) const noexcept
    {
      return stc::morton(a_x, a_y);
    }
};

/*==============================================================================
 * Meta-programming using a vector
 *============================================================================*/

/// A nested loop over some integer vector of rank N
template <array_size_type I, array_size_type N, typename F>
struct NestedLoop
{
  // Unit stride
  HOSTDEVICE static constexpr void loop(const IVec<N>& lo,
                                        const IVec<N>& hi,
                                        IVec<N>& iv,
                                        F&& f)
    {
      int& i = iv[I-1];
      const int i_end = hi[I-1] + 1;
      for (i = lo[I-1]; i != i_end; ++i)
        {
          NestedLoop<I-1, N, F>::loop(lo, hi, iv, std::forward<F>(f));
        }
    }
  // Unit stride
  static constexpr void loop_host(const IVec<N>& lo,
                                  const IVec<N>& hi,
                                  IVec<N>& iv,
                                  F&& f)
    {
      int& i = iv[I-1];
      const int i_end = hi[I-1] + 1;
      for (i = lo[I-1]; i != i_end; ++i)
        {
          NestedLoop<I-1, N, F>::loop_host(lo, hi, iv, std::forward<F>(f));
        }
    }
  // Given stride
  HOSTDEVICE static constexpr void loopStride(const IVec<N>& lo,
                                              const IVec<N>& hi,
                                              const IVec<N>& stride,
                                              IVec<N>& iv,
                                              F&& f)
    {
      int& i = iv[I-1];
      const int i_last = hi[I-1];
      const int i_stride = stride[I-1];
      for (i = lo[I-1]; i <= i_last; i += i_stride)
        {
          NestedLoop<I-1, N, F>::loopStride(lo, hi, stride, iv,
                                            std::forward<F>(f));
        }
    }
  // Given stride
  static constexpr void loopStride_host(const IVec<N>& lo,
                                        const IVec<N>& hi,
                                        const IVec<N>& stride,
                                        IVec<N>& iv,
                                        F&& f)
    {
      int& i = iv[I-1];
      const int i_last = hi[I-1];
      const int i_stride = stride[I-1];
      for (i = lo[I-1]; i <= i_last; i += i_stride)
        {
          NestedLoop<I-1, N, F>::loopStride_host(lo, hi, stride, iv,
                                                 std::forward<F>(f));
        }
    }
};

// Terminate when N = 1 (execute the body of the loop)
template <array_size_type N, typename F>
struct NestedLoop<1, N, F>
{
  // Unit stride
  HOSTDEVICE static constexpr void loop(const IVec<N>& lo,
                                        const IVec<N>& hi,
                                        IVec<N>& iv,
                                        F&& f)
    {
      int& i = iv[0];
      const int i_end = hi[0] + 1;
      for (i = lo[0]; i != i_end; ++i)
        {
          f(iv);
        }
    }
  // Unit stride
  static constexpr void loop_host(const IVec<N>& lo,
                                  const IVec<N>& hi,
                                  IVec<N>& iv,
                                  F&& f)
    {
      int& i = iv[0];
      const int i_end = hi[0] + 1;
      for (i = lo[0]; i != i_end; ++i)
        {
          f(iv);
        }
    }
  // Given stride
  HOSTDEVICE static constexpr void loopStride(const IVec<N>& lo,
                                              const IVec<N>& hi,
                                              const IVec<N>& stride,
                                              IVec<N>& iv,
                                              F&& f)
    {
      int& i = iv[0];
      const int i_last = hi[0];
      const int i_stride = stride[0];
      for (i = lo[0]; i <= i_last; i += i_stride)
        {
          f(iv);
        }
    }
  // Given stride
  static constexpr void loopStride_host(const IVec<N>& lo,
                                        const IVec<N>& hi,
                                        const IVec<N>& stride,
                                        IVec<N>& iv,
                                        F&& f)
    {
      int& i = iv[0];
      const int i_last = hi[0];
      const int i_stride = stride[0];
      for (i = lo[0]; i <= i_last; i += i_stride)
        {
          f(iv);
        }
    }
};          

// In case N = 0, don't do anything
template <array_size_type N, typename F>
struct NestedLoop<0, N, F>
{
  HOSTDEVICE static constexpr void loop(const IVec<N>& lo,
                                        const IVec<N>& hi,
                                        IVec<N>& iv,
                                        F&& f)
    { }
  static constexpr void loop_host(const IVec<N>& lo,
                                  const IVec<N>& hi,
                                  IVec<N>& iv,
                                  F&& f)
    { }
  HOSTDEVICE static constexpr void loopStride(const IVec<N>& lo,
                                              const IVec<N>& hi,
                                              const IVec<N>& stride,
                                              IVec<N>& iv,
                                              F&& f)
    { }
  static constexpr void loopStride_host(const IVec<N>& lo,
                                        const IVec<N>& hi,
                                        const IVec<N>& stride,
                                        IVec<N>& iv,
                                        F&& f)
    { }
};

/// Execute function object 'f(iv)' within nested loop over bounds
/** Note: the bounds are inclusive, lo <= iv <= hi
 */
template <typename REP1, typename REP2, typename F,
          std::enable_if_t<
            std::is_base_of<stc::StcOp, std::decay_t<REP1>>::value
            &&
            std::is_base_of<stc::StcOp, std::decay_t<REP2>>::value, int> = 0
          >
HOSTDEVICE constexpr void
nestedLoop(REP1&& loRep, REP2&& hiRep, F&& f)
{
  constexpr array_size_type N = std::decay_t<REP1>::size();
  static_assert(N == std::decay_t<REP2>::size(),
                "Nested loop must have same size for upper and lower bounds");
  const IVec<N> lo = loRep;  // Evaluate if have expression
  const IVec<N> hi = hiRep;
  IVec<N> iv;                // Current index in the loop
  NestedLoop<N, N, F>::loop(lo, hi, iv, std::forward<F>(f));
}

/// Execute function object 'f(iv)' within nested loop over bounds
/** Note: the bounds are inclusive, lo <= iv <= hi
 */
template <typename REP1, typename REP2, typename F,
          std::enable_if_t<
            std::is_base_of<stc::StcOp, std::decay_t<REP1>>::value
            &&
            std::is_base_of<stc::StcOp, std::decay_t<REP2>>::value, int> = 0
          >
constexpr void
nestedLoop_host(REP1&& loRep, REP2&& hiRep, F&& f)
{
  constexpr array_size_type N = std::decay_t<REP1>::size();
  static_assert(N == std::decay_t<REP2>::size(),
                "Nested loop must have same size for upper and lower bounds");
  const IVec<N> lo = loRep;  // Evaluate if have expression
  const IVec<N> hi = hiRep;
  IVec<N> iv;                // Current index in the loop
  NestedLoop<N, N, F>::loop_host(lo, hi, iv, std::forward<F>(f));
}

/// Execute function object 'f(iv)' within nested loop over bounds given stride
/** Note: the bounds are inclusive, lo <= iv <= hi
 */
template <typename REP1, typename REP2, typename REP3, typename F,
          std::enable_if_t<
            std::is_base_of<stc::StcOp, std::decay_t<REP1>>::value
            &&
            std::is_base_of<stc::StcOp, std::decay_t<REP2>>::value
            &&
            std::is_base_of<stc::StcOp, std::decay_t<REP3>>::value, int> = 0
          >
HOSTDEVICE constexpr void
nestedLoop(REP1&& loRep, REP2&& hiRep, REP3&& strideRep, F&& f)
{
  constexpr array_size_type N = std::decay_t<REP1>::size();
  static_assert(N == std::decay_t<REP2>::size(),
                "Nested loop must have same size for upper and lower bounds");
  static_assert(N == std::decay_t<REP3>::size(),
                "Nested loop must have same size for bounds and stride");
  const IVec<N> lo     = loRep;  // Evaluate if have expression
  const IVec<N> hi     = hiRep;
  const IVec<N> stride = strideRep;
  IVec<N> iv;                    // Current index in the loop
  NestedLoop<N, N, F>::loopStride(lo, hi, stride, iv, std::forward<F>(f));
}

/// Execute function object 'f(iv)' within nested loop over bounds given stride
/** Note: the bounds are inclusive, lo <= iv <= hi
 */
template <typename REP1, typename REP2, typename REP3, typename F,
          std::enable_if_t<
            std::is_base_of<stc::StcOp, std::decay_t<REP1>>::value
            &&
            std::is_base_of<stc::StcOp, std::decay_t<REP2>>::value
            &&
            std::is_base_of<stc::StcOp, std::decay_t<REP3>>::value, int> = 0
          >
constexpr void
nestedLoop_host(REP1&& loRep, REP2&& hiRep, REP3&& strideRep, F&& f)
{
  constexpr array_size_type N = std::decay_t<REP1>::size();
  static_assert(N == std::decay_t<REP2>::size(),
                "Nested loop must have same size for upper and lower bounds");
  static_assert(N == std::decay_t<REP3>::size(),
                "Nested loop must have same size for bounds and stride");
  const IVec<N> lo     = loRep;  // Evaluate if have expression
  const IVec<N> hi     = hiRep;
  const IVec<N> stride = strideRep;
  IVec<N> iv;                    // Current index in the loop
  NestedLoop<N, N, F>::loopStride_host(lo, hi, stride, iv, std::forward<F>(f));
}

}  // namespace stc

#endif  /* ! defined _STCVECTOR_H_ */
