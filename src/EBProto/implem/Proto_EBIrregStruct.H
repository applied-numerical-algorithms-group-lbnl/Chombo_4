#ifndef _Proto_EBIrregStruct_H_
#define _Proto_EBIrregStruct_H_

#include <cmath>
#include <cstdlib>
#include <memory>
#include "base/Proto_Box.H"
#include "EBProto.H"

 namespace Proto
 {
  template <CENTERING cent, class data_t,  unsigned int ncomp>
  class
  EBIrregStruct
  {
public:
    data_t*        m_startPtr;
    unsigned int   m_varsize;
    unsigned int   m_offset;
    Point          m_index;

    EBIrregStruct() : m_index()
    {
      m_startPtr = nullptr;
      m_varsize  = 0;
      m_offset   = 0;
    }
    ~EBIrregStruct(){};
   };
}

#endif
