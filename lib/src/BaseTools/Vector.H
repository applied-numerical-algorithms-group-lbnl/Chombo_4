#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _VECTOR_H_
#define _VECTOR_H_

#include <vector>
#include <iostream>
#include <algorithm>

#include "BaseNamespaceHeader.H"

template <class T>
class Vector: public std::vector<T>
{
public:
  Vector() : std::vector<T>()
  {
  }


  Vector(unsigned int num) : std::vector<T>(num)
  {
  }

  const std::vector<T>& constStdVector() const
  {
    return *this;
  }

  std::vector<T>& stdVector()
  {
    return *this;
  }

  Vector(unsigned int num, const T& val) : std::vector<T>(num, val)
  {
  }

  virtual ~Vector()
  {
  }

  void sort()
  {
    std::sort(this->begin(), this->end());
  }
};

//class ostream;
using std::ostream;
template <class T>
ostream& operator<<(ostream& os, const Vector<T>& vec)
{
  for (unsigned int i=0; i<vec.size(); i++) os<<vec[i]<<" ";
  return os;
}

#include "BaseNamespaceFooter.H"
#endif
