#ifndef _MOMENTITERATOR_H_
#define _MOMENTITERATOR_H_

#include <array>

#include "StcVector.hpp"
#include "MultiIndexMath.hpp"

namespace hoeb
{
/** This provides a way to iterate of the vector indices of moments
 * Its size is compile time known
 */
template <unsigned Dim, unsigned Order>
class MomentIterator
{
public:
  // generate the size
  HOSTDEVICE constexpr static unsigned
  size()
  {
    return nCk(Dim + Order, Order);
  }

  ///
  HOSTDEVICE
  constexpr MomentIterator() : m_current(0) {}

  /// Destructor
  ~MomentIterator() = default;

  /// Copy constructor
  MomentIterator(const MomentIterator&) = default;

  /// Move constructor
  MomentIterator(MomentIterator&&) = default;

  /// Assignment constructor
  MomentIterator& operator=(const MomentIterator&) = default;

  /// Move assignment constructor
  MomentIterator& operator=(MomentIterator&&) = default;

  using iterator = typename std::array<stc::IVec<Dim>, size()>::iterator;
  using const_iterator = typename std::array<stc::IVec<Dim>, size()>::const_iterator;

  /// standard library style iterators. Allow range-based for loops
  HOSTDEVICE
  iterator
  begin()
  {
    return s_indicies.begin();
  }

  ///
  HOSTDEVICE
  iterator
  end()
  {
    return s_indicies.end();
  }

  ///
  HOSTDEVICE const_iterator
  begin() const
  {
    return s_indicies.begin();
  }

  ///
  HOSTDEVICE
  const_iterator
  end() const
  {
    return s_indicies.cend();
  }

  ///
  HOSTDEVICE
  void
  reset()
  {
    m_current = 0;
  }

  /// increment
  HOSTDEVICE
  void
  operator++()
  {
    m_current++;
  }

  /// Returns the moment of the current multi-index of this MomentIterator.
  HOSTDEVICE
  unsigned
  operator()() const
  {
    return m_current;
  }

  ///
  HOSTDEVICE
  stc::IVec<Dim>
  momentIndex() const
  {
    return s_indicies[m_current];
  }

  /// Returns true if this MomentIterator's location is within its IndexedMoment.
  HOSTDEVICE
  bool
  ok() const
  {
    return (m_current < size());
  }

  /// Returns the moment of the current multi-index of this MomentIterator.
  HOSTDEVICE constexpr static unsigned
  indexOf(const stc::IVec<Dim>& a_index)
  {
    assert(a_index[stc::indexMinElem(a_index)] >= 0);
    assert(stc::sum(a_index) <= Order);
    int idx = a_index[0];
    int rem = Order;
    for (int d = Dim - 1; d != 0; d--)
    {
      idx += nCk(d + 1 + rem, rem) - nCk(d + 1 + rem - a_index[d], rem - a_index[d]);
      rem -= a_index[d];
    }
    return idx;
  }

private:
  // constexpr static //FIXME would be great if this was constexpr
  std::array<stc::IVec<Dim>, size()> s_indicies = MomentIterator<Dim, Order>::setMultiIndicies();
  /// iterator index
  unsigned m_current;

  // generate the moment power/index list
  HOSTDEVICE constexpr static std::array<stc::IVec<Dim>, MomentIterator<Dim, Order>::size()>
  setMultiIndicies()
  {
    constexpr unsigned size = MomentIterator<Dim, Order>::size();
    std::array<stc::IVec<Dim>, size> indicies;
    auto index = stc::make_IVec<Dim>::zero();
    for (unsigned ix = 0; ix < size; ++ix)
    {
      // If the sum is too large, shift extras to the right
      for (int d = 0; (stc::sum(index) > (int)(Order)) && (d < Dim - 1); ++d)
      {
        index[d] = 0;
        ++index[d + 1];
      }
      indicies[ix] = index;
      ++index[0];
    }
    return indicies;
  }
};

} // namespace hoeb
#endif
