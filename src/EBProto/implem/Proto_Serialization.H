
#ifndef _Proto_Serialization_H_
#define _Proto_Serialization_H_



namespace ProtoSer
{
  ///return serialization length (in number of characters)
  template <typename data_t> inline size_t
  linearCharBufSize(const data_t & a_dat)
  {
    return a_dat.linearSize();
  }


  ///define from string of bits
  template <typename data_t> inline  void
  linearInFromCharBuf(data_t & a_dat, const char* a_buf)
  {
    a_dat.linearIn(a_buf); 
  }

  ///fill string of bits from internals
  template <typename data_t> inline void
  linearOutToCharBuf(char* a_buf, const data_t & a_dat)
  {
    a_dat.linearOut(a_buf); 
  }

///and now for all the template specializations needed for intrinsic data types.
  template < >
  inline size_t
  linearCharBufSize<double>(const double & a_dat)
  {
    return sizeof(double);
  }

  template < >
  inline void
  linearInFromCharBuf<double>(double & a_dat, const char* a_buf)
  {
    const double* dbuf = (const double *)a_buf;
    a_dat = *dbuf;
  }


  template < >
  inline void
  linearOutToCharBuf<double>(char* a_buf, const double & a_dat)
  {
    double* dbuf = (double *)a_buf;
    *dbuf = a_dat;
  }



  template < >
  inline size_t
  linearCharBufSize<int>(const int & a_dat)
  {
    return sizeof(int);
  }

  template < >
  inline void
  linearInFromCharBuf<int>(int & a_dat, const char* a_buf)
  {
    const int* dbuf = (const int *)a_buf;
    a_dat = *dbuf;
  }


  template < >
  inline void
  linearOutToCharBuf<int>(char* a_buf, const int & a_dat)
  {
    int* dbuf = (int *)a_buf;
    *dbuf = a_dat;
  }

  template < >
  inline size_t
  linearCharBufSize<float>(const float & a_dat)
  {
    return sizeof(float);
  }

  template < >
  inline void
  linearInFromCharBuf<float>(float & a_dat, const char* a_buf)
  {
    const float* dbuf = (const float *)a_buf;
    a_dat = *dbuf;
  }


  template < >
  inline void
  linearOutToCharBuf<float>(char* a_buf, const float & a_dat)
  {
    float* dbuf = (float *)a_buf;
    *dbuf = a_dat;
  }

}


#endif
