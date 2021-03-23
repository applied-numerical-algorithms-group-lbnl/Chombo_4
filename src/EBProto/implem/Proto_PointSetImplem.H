#ifndef __Proto_PointSetImplem__
#define __Proto_PointSetImplem__
namespace Proto
{

///
PointSet::
PointSet(const Box& a_box)
{
  for(BoxIterator bit = a_box.begin(); bit != a_box.end();  ++bit)
  {
    *this |= *bit;
  }
}
///
  PointSet::
  PointSet(const PointSet& a_sivs)
  {
    m_stdSet = a_sivs.m_stdSet; 
  }
///
  void 
  PointSet::
  define(const Box& a_box)
  {
    *this = PointSet(a_box);
  }
///
  void 
  PointSet::
  define(const PointSet& a_sivs)
  {
    m_stdSet = a_sivs.m_stdSet; 
  }
///
  PointSet& 
  PointSet::
  operator=(const PointSet& a_sivs)
  {
    m_stdSet = a_sivs.m_stdSet; 
    return *this;
  }
///
  PointSet& 
  PointSet::
  operator|=(const PointSet& a_sivs)
  {
    const std::set<Point> inputset = a_sivs.m_stdSet;
    std::set<Point>::iterator it;
    for(it = inputset.begin(); it!=  inputset.end(); ++it)
    {
      m_stdSet.insert(*it);
    }
    return *this;
  }
  
///

  PointSet & 
  PointSet::
  operator|=(const Point& a_iv)
  {
    m_stdSet.insert(a_iv);
    return *this;
  }
///
  PointSet& 
  PointSet::
  operator|=(const Box& a_box)
  {
    for(BoxIterator bit = a_box.begin(); bit != a_box.end();  ++bit)
    {
      m_stdSet.insert(*bit);
    }
    return *this;
  }
///
  PointSet& 
  PointSet::
  operator&=(const PointSet& a_sivs)
  {
    if(&a_sivs != this)
    {
      std::set<Point> newSet;
      std::set<Point>::iterator it;
      for(it = m_stdSet.begin(); it != m_stdSet.end(); ++it)
      {
        const Point& iv = *it;
        if(contains(iv) && a_sivs.contains(iv))
        {
          newSet.insert(iv);
        }
      }
      m_stdSet = newSet;
    }
    return *this;
  }

///and
  PointSet& 
  PointSet::
  operator&=(const Box& a_box)
  {
    std::set<Point> newSet;
    std::set<Point>::iterator it;

    for(it = m_stdSet.begin(); it!=  m_stdSet.end(); ++it)
    {
      const Point& iv = *it;
      if(a_box.contains(iv))
      {
        newSet.insert(iv);
      }
    }
    m_stdSet = newSet;
    return *this;
  }
///not
  PointSet& 
  PointSet::
  operator-=(const PointSet& a_sivs)
  {
    std::set<Point>::iterator it;
    //leaving out the ++it because  erase 
    for(it = m_stdSet.begin(); it!=  m_stdSet.end(); )
    {
      if(a_sivs.contains(*it))
      {
        m_stdSet.erase(it++);
      }
      else
      {
        ++it;
      }
    }
    return *this;
  }
///not
  PointSet& 
  PointSet::
  operator-=(const Point& a_iv)
  {
    if(contains(a_iv))
    {
      m_stdSet.erase(m_stdSet.find(a_iv));
    }
    return *this;
  }
///not
  PointSet& 
  PointSet::
  operator-=(const Box& a_box)
  {
    for(BoxIterator bit = a_box.begin(); bit != a_box.end();  ++bit)
    {
      *this -= *bit;
    }
    return *this;
  }
///
  bool 
  PointSet::
  operator==(const PointSet& a_lhs) const
  {
    if(a_lhs.m_stdSet.size() != m_stdSet.size())
    {
      return false;
    }

    bool retval = true;
    std::set<Point>::iterator it;
    for(it = m_stdSet.begin(); it!=  m_stdSet.end(); ++it)
    {
      if((!contains(*it)) || (!a_lhs.contains(*it)))
      {
        retval = false;
        break;
      }
    }
    return retval;
  }

///
  bool 
  PointSet::
  contains(const Point& a_iv) const
  {
    std::set<Point>::iterator it = m_stdSet.find(a_iv);
    return (it != m_stdSet.end());
  }

///
  bool 
  PointSet::
  contains(const Box& a_box) const
  {
    bool retval = true;
    for(BoxIterator bit = a_box.begin(); bit != a_box.end();  ++bit)
    {
      if(!contains(*bit))
      {
        retval = false;
        break;
      }
    }
    return retval;
  }

///
  void 
  PointSet::
  grow(int igrow)
  {
    PointSet newSet;
    std::set<Point>::iterator it;
    for(it = m_stdSet.begin(); it != m_stdSet.end(); ++it)
    {
      const Point& iv = *it;
      Box grid(iv, iv);
      Box bgrow = grid.grow(igrow);
      newSet |= bgrow;
    }
    *this = newSet;
  }

///
  void 
  PointSet::
  grow(int idir, int igrow)
  {
    PointSet newSet;
    std::set<Point>::iterator it;
    for(it = m_stdSet.begin(); it != m_stdSet.end(); ++it)
    {
      const Point& iv = *it;
      Box grid(iv, iv);
      Box grown = grid.grow(idir, igrow);
      newSet |= grown;
    }
    *this = newSet;
  }


///
  void 
  PointSet::
  refine(int iref)
  {
    PointSet newSet;
    std::set<Point>::iterator it;
    for(it = m_stdSet.begin(); it != m_stdSet.end(); ++it)
    {
      const Point& iv = *it;
      Box grid(iv, iv);
      Box refgrid = grid.refine(iref);
      newSet |= refgrid;
    }
    *this = newSet;
  }

///
  void 
  PointSet::
  coarsen(int iref)
  {
    std::set<Point> newSet;
    std::set<Point>::iterator it;
    for(it = m_stdSet.begin(); it != m_stdSet.end(); ++it)
    {
      const Point& iv = *it;
      Point coariv = iv/iref;
      newSet.insert(coariv);
    }
    m_stdSet = newSet;
  }

///
  void 
  PointSet::
  shift(const Point& a_iv)
  {
    std::set<Point> newSet;
    std::set<Point>::iterator it;
    for(it = m_stdSet.begin(); it != m_stdSet.end(); ++it)
    {
      Point iv = *it;
      iv += a_iv;
      newSet.insert(iv);
    }
    m_stdSet = newSet;
  }

///
  void 
  PointSet::
  clear()
  {
    std::set<Point> newSet;
    m_stdSet = newSet;
  }

///
  Box 
  PointSet::
  minBox() const
  {
    int bignum = 100000;
    Point lo = bignum*Point::Unit();
    Point hi =-bignum*Point::Unit();
    std::set<Point>::iterator it;
    for(it = m_stdSet.begin(); it != m_stdSet.end(); ++it)
    {
      const Point& iv = *it;
      for(int idir = 0; idir < DIM; idir++)
      {
        lo[idir] = std::min(lo[idir], iv[idir]);
        hi[idir] = std::max(hi[idir], iv[idir]);
      }
    }
  
    Box retval(lo, hi);
    return retval;
  }

///
  bool 
  PointSet::
  isEmpty() const
  {
    return (m_stdSet.size() == 0);
  }
///
  void 
  PointSet::
  makeEmpty() 
  {
    clear();
  }

///
  int
  PointSet::
  numPts() const
  {
    return m_stdSet.size();
  }

///
  PointSetIterator::
  PointSetIterator()
  {
    m_ivs = NULL;
  }

///
  PointSetIterator::
  PointSetIterator(const PointSet& ivs)
  {
    m_ivs = &ivs;
    m_iter = m_ivs->m_stdSet.begin();
  }

///
  void 
  PointSetIterator::
  define(const PointSet& a_ivs)
  {
    m_ivs = &a_ivs;
    m_iter = m_ivs->m_stdSet.begin();
  }

///
  const Point& 
  PointSetIterator::
  operator()() const 
  {
    return *m_iter;
  }

///
  bool 
  PointSetIterator::
  ok() const
  {
    PR_assert(m_ivs != NULL);
    return (m_iter != m_ivs->m_stdSet.end());
  }

///
  void 
  PointSetIterator::
  operator++()
  {
    m_iter++;
  }

///
  void 
  PointSetIterator::
  begin()
  {
    PR_assert(m_ivs != NULL);
    m_iter = m_ivs->m_stdSet.begin();
  }

///
  void 
  PointSetIterator::
  end()
  {
    PR_assert(m_ivs != NULL);
    m_iter = m_ivs->m_stdSet.end();
  }

///
  void 
  PointSetIterator::
  clear()
  {
    m_ivs = NULL;
  }
}
#endif
