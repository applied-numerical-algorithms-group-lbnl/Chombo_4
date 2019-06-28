#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _REFCOUNTEDPTR___H_
#define _REFCOUNTEDPTR___H_
#include "BaseNamespaceHeader.H"

template<class T>
class RefCountedPtr: public std::shared_ptr<T>
{
public:
  RefCountedPtr():std::shared_ptr<T>()
  {
  }
  
  virtual ~RefCountedPtr()
  {
  }

  RefCountedPtr(const T* a_ptr):std::shared_ptr<T>(a_ptr)
  {
  }


  RefCountedPtr(T* a_ptr):std::shared_ptr<T>(a_ptr)
  {
  }

  bool operator==(const T* const a_ptr) const
  {
    return (&(*this) == a_ptr);
  }
  
  unsigned int refCount() const
  {
    return this->use_count();
  }

};

#include "BaseNamespaceFooter.H"
#endif
