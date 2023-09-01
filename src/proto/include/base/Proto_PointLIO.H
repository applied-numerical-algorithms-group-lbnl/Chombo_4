#ifndef __PROTO_SPMD_H__
#define __PROTO_SPMD_H__

#include "Proto_Point.H"
#include <vector>

#ifdef PROTO_CUDA
#include "implem/Proto_DisjointBoxLayout.H"
#endif

namespace Proto
{

  template <class T>
  int linearSize(const T& inputT);

  template <class T>
  void linearIn(T& a_outputT, const void* const inBuf);

  template <class T>
  void linearOut(void* const a_outBuf, const T& inputT);



  //===========  specializations =============================

  template < >
  void linearIn<std::vector<Point>>(std::vector<Point>& a_outputT, const void* const a_inBuf)
  {
    std::size_t* buffer = (std::size_t*)a_inBuf;
    a_outputT.resize(*buffer);
    char* bufferC = (char*)a_inBuf;
    bufferC+=sizeof(std::size_t);
    memcpy(a_output_T.data(), bufferC, sizeof(Point)*(*buffer));
  }

  template < >
  void linearOut<std::vector<Point>>(void* const a_outBuf, const std::vector<Point>& a_inputT)
  {
    std::size_t length = &((std::size_t*)a_outBuf);
    length= a_inputT.size();
    char* buffer = (char*)a_outBuf;
    buffer += sizeof(std::size_t);
    memcpy(buffer, a_inputT.data(), sizeof(Point)*length);
  }

  template < >
  int linearSize<std::vector<Point>>(const std::vector<Point>& a_input)
  {
    return sizeof(std::size_t)+a_input.size()*sizeof(Point);
  }
  //=======================default implementations

  template <class T>
  int linearSize(const T& inputT)
  {
    return inputT.linearSize();
  }

  template <class T>
  void linearIn(T& a_outputT, const void* const inBuf)
  {
    a_outputT.linearIn(inBuf);
  }
  
  template <class T>
  void linearOut(void* const a_outBuf, const T& inputT)
  {
    inputT.linearOut(a_outBuf);
  }
  
  
} // end namespace Proto

#endif  // __PROTO_SPMD_H__
