
#ifndef _Proto_EBIrregStruct_H_
#define _Proto_EBIrregStruct_H_

#include <cmath>
#include <cstdlib>
#include <memory>
#include "Proto_Box.H"
#include "EBProto.H"

namespace Proto
{
  template <class data_t, unsigned int ncomp>
  struct
  EBIrregStruct
  {
    data_t*        m_startPtr;
    unsigned int   m_varsize;
    unsigned int   m_offset;
    Point          m_index;
  };

}

#endif
