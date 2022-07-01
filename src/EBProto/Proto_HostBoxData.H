#ifndef _PROTO_HostBoxData_H_
#define _PROTO_HostBoxData_H_

#include "base/Proto_Box.H"
#include "Chombo_parstream.H"
#include "implem/Proto_Serialization.H"

namespace Proto
{
/// Multidimensional Rectangular Array that lives on the host
/**
   this is mainly for graph-like data thave lives on the host
*/
  template <class data_t=double, unsigned int ncomp=1>
  class HostBoxData 
  {
  public:

    ///
    HostBoxData();

    int nComp() const
    {
      return ncomp;
    }
    
    ///
    HostBoxData(const Box& a_box, bool a_useStackAlloc = false);

    /// aliasing for when the device is a lie
    HostBoxData(const BoxData<data_t, ncomp>& a_input);
    
    /// aliasing for when the device is a lie
    void define(const BoxData<data_t, ncomp>& a_input);
    
    /// 
    void define(const Box& a_box, bool a_useStackAlloc = false);
   
    /// Destructor.
    virtual ~HostBoxData(); 

    ///
    template<unsigned int Csrc>
    void copy(const HostBoxData<data_t,Csrc>  & a_dsrc,
              const Box                       & a_srcBox,
              unsigned int                      a_srcComp,
              const Box                       & a_destBox,
              unsigned int                      a_destComp,
              unsigned int                      a_numcomp);
    
    ///this DOES include sizeof(data_t)  
    virtual std::size_t charsize(const Box&   a_bx,
                                 unsigned int a_startcomp,
                                 unsigned int a_numcomps,
                                 bool a_printStuff = false) const
    {
      PR_TIME("hostboxdata::charsize");
      using ProtoSer::linearCharBufSize;

      Box binter = a_bx & m_box;
      
      size_t regsize = 0;
      for(auto bit = binter.begin(); bit != binter.end(); ++bit)
      {
        Point pt = *bit;
        for(int icomp = 0; icomp < ncomp; icomp++)
        {
          regsize +=  linearCharBufSize<data_t>((*this)(pt, icomp));
        }
      }

      if(a_printStuff)
      {
        Chombo4::pout() << "HostBoxData::charsize for box "<< a_bx << " returning " << regsize << endl;
      }
      return regsize;
    }

    virtual void  linearIn(void*                  a_buf,
                           size_t                 a_expectedCharSize,
                           const Box       &      a_bx,
                           unsigned int           a_startcomp,
                           unsigned int           a_numcomps,
                           bool a_printStuff = false) 
    {
      PR_TIME("hostboxdata::linearin");
      PR_assert((a_startcomp + a_numcomps-1) < ncomp)
      using ProtoSer::linearCharBufSize;
      using ProtoSer::linearInFromCharBuf;

      if(a_printStuff)
      {
        Chombo4::pout() << "entering HostBoxData::linearIn for box "<< a_bx << endl;
      }
      char* charbuf = (char*)a_buf;
      std::size_t  charsiz = 0;
      Box binter = a_bx & m_box;
      for(auto bit = binter.begin(); bit != binter.end(); ++bit)
      {
        Point pt = *bit;
        for(int icomp = 0; icomp < a_numcomps; icomp++)
        {
          data_t& dataval = (*this)(pt, icomp + a_startcomp);
          linearInFromCharBuf<data_t>(dataval , charbuf);
          std::size_t datasiz =  linearCharBufSize<data_t>(dataval);
          charbuf  += datasiz;
          charsiz  += datasiz;
        }
      }
      if(a_printStuff)
      {
        Chombo4::pout() << "HostBoxData::linearIn leaving for box "<< a_bx << ", with charsiz = " << charsiz  << endl;
        if(charsiz != a_expectedCharSize)
        {
          Chombo4::pout() << "HostBoxData::linearIn: WARNING: "
                          << "expected charsiz = " << a_expectedCharSize
                          << ", actual size = "<< charsiz << endl;
        }
      }
    }
    virtual void linearOut(void*               a_buf,
                           size_t              a_expectedCharSize,
                           const Box       &   a_bx,
                           unsigned int        a_startcomp,
                           unsigned int        a_numcomps,
                           bool a_printStuff = false) const
    {
      PR_TIME("hostboxdata::linearout");
      PR_assert((a_startcomp + a_numcomps-1) < ncomp);
      using ProtoSer::linearCharBufSize;
      using ProtoSer::linearOutToCharBuf;
      if(a_printStuff)
      {
        Chombo4::pout() << "HostBoxData::linearOut for box "<< a_bx << endl;
      }
      char* charbuf = (char*)a_buf;
      std::size_t  charsiz = 0;
      Box binter = a_bx & m_box;
      for(auto bit = binter.begin(); bit != binter.end(); ++bit)
      {
        Point pt = *bit;
        for(int icomp = 0; icomp < ncomp; icomp++)
        {
          const data_t& dataval = (*this)(pt, icomp);  // 
          linearOutToCharBuf<data_t>(charbuf, dataval);
          std::size_t datasiz = linearCharBufSize<data_t>(dataval);
          charbuf +=  datasiz;
          charsiz +=  datasiz;
        }
      }
      if(a_printStuff)
      {
        Chombo4::pout() << "HostBoxData::linearOut for box "<< a_bx << ", with charsiz = " << charsiz  << endl;
        if(charsiz != a_expectedCharSize)
        {
          Chombo4::pout() << "HostBoxData::linearOut: WARNING expected charsiz ="
                          << a_expectedCharSize << endl;
        }
      }
    }
    static int preAllocatable()
    {
      return 1; // dyanmic allocatable.
    }

    ///
    inline const data_t& operator()(const Point& a_pt, 
                                    unsigned int  a_c) const;

    ///
    inline data_t& operator()(const Point& a_pt, 
                              unsigned int  a_c);
  
    ///
    inline Box box() const {return m_box;};
    
    ///
    inline std::size_t size() const {return m_box.size()*ncomp;};

    ///
    inline bool defined() const {return bool(m_data);};
   

    ///
    void setVal(const data_t& a_val);
  

    ///
    data_t* data()
    {
      return m_rawPtr;
    }

    ///
    const data_t* data() const
    {
      return m_rawPtr;
    }
    inline size_t index(const Point a_pt,
                        unsigned int a_c) const
    {
      PR_assert(m_box.contains(a_pt));

      PR_assert(a_c < ncomp);

      size_t m = m_box.size();
      size_t k = m_box.index(a_pt);
      return k+m*a_c;
    }

  private:
    /// Copy constructor/assignment forbidden for all the usual reasons.
    HostBoxData(const HostBoxData<data_t,ncomp>& a_src);
    HostBoxData& operator=(const HostBoxData<data_t,ncomp>& a_src);
    bool m_stackAlloc;
    Box                       m_box;    ///< Box defining the domain of *this
    ::std::shared_ptr<data_t> m_data;   ///< Data array
    data_t*                   m_rawPtr; ///< Raw pointer to the data
  };

} //end Proto namespace

#include "implem/Proto_HostBoxDataImplem.H"
#endif //end include guard
