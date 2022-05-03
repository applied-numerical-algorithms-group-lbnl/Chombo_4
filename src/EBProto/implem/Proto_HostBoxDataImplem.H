#ifndef __Proto_HostBOXDATAIMPLEM__
#define __Proto_HostBOXDATAIMPLEM__


#pragma once

namespace Proto
{

// Default Constructor
  template <class T, unsigned int C> 
  HostBoxData<T,C>::HostBoxData()
    : m_rawPtr(NULL) {}

// Box Constructor
  template <class T, unsigned int C> 
  HostBoxData<T,C>::HostBoxData(const Box& a_box, bool a_stackAlloc) 
  {
    define(a_box, a_stackAlloc);
  }

// Define
  template <class T, unsigned int C> 
  void HostBoxData<T,C>::define(const Box& a_box, bool a_stackAlloc)
  {
    PR_TIME("HostBoxData::define(Box) (memory allocation)");
    m_box=a_box;

    std::size_t count = this->size();
    m_data=::std::shared_ptr<T>(new T [count], [](T* p) { delete[] p;});
    m_rawPtr = m_data.get();

  }

// Destructor
  template <class T, unsigned int C> 
  HostBoxData<T,C>::~HostBoxData()
  {}


  // Copy From
  template <class T,unsigned int C>
  template<unsigned int Csrc>
  void HostBoxData<T,C>::copy(const HostBoxData<T,Csrc> &  a_src,
                              const Box                 &  a_srcBox,
                              unsigned int                 a_srcComp,
                              const Box                 &  a_destBox,
                              unsigned int                 a_destComp,
                              unsigned int                 a_numcomps)
  {
    PR_TIME("HostBoxData::copy");
    PR_assert(a_srcBox.size() == a_destBox.size());
    Point offset = a_destBox.low() - a_srcBox.low();
    Box interbox = a_srcBox & m_box;
    for(BoxIterator bit = interbox.begin(); bit != interbox.end();  ++bit)
    {
      Point srcpt = *bit;
      Point dstpt = srcpt + offset;
      for (unsigned int cc = 0; cc < a_numcomps; cc++)
      {
        int scomp = a_srcComp  + cc;
        int dcomp = a_destComp + cc;
        (*this)(dstpt, dcomp) = a_src(srcpt, scomp);
      }
    }
  }
  ///
  template <class T, unsigned int C>
  void HostBoxData<T,C>::setVal(const T& a_val)
  {
    PR_TIME("HostBoxData::setval");
    for(BoxIterator bit = m_box.begin(); bit != m_box.end();  ++bit)
    {
      Point pt = *bit;
      for (unsigned int cc = 0; cc < C; cc++)
      {
        (*this)(pt, cc) = a_val;
      }

    }
  }

  template <class T, unsigned int C>
  T& HostBoxData<T,C>::operator()(const Point& a_pt, unsigned int a_c)
  {
    int index = m_box.index(a_pt);
    int vsize = m_box.size();
    return m_rawPtr[a_c*vsize + index];
  }

  ///
  template <class T, unsigned int C>
  const T& HostBoxData<T,C>::operator()(const Point& a_pt, unsigned int a_c) const
  {
    int index = m_box.index(a_pt);
    int vsize = m_box.size();
    return m_rawPtr[a_c*vsize + index];
  }
}
#endif
