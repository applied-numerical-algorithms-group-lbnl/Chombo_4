
#ifndef _Proto_EBIrregStruct_H_
#define _Proto_EBIrregStruct_H_

#include <cmath>
#include <cstdlib>
#include <memory>
#include "Proto_Box.H"
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

    template<unsigned int ncomp2>
    EBIrregStruct(EBIrregStruct<cent, data_t,ncomp2> a_in)
    {
      m_startPtr = a_in.m_startPtr;
      m_varsize  = a_in.m_varsize;
      m_offset   = a_in.m_offset;
      m_index    = a_in.m_index;
    }

    EBIrregStruct() : m_index()
    {
      m_startPtr = nullptr;
      m_varsize  = 0;
      m_offset   = 0;
    }
    virtual ~EBIrregStruct(){}
   };

}





#endif
