
#ifndef _Proto_EBHostFluxData_H_
#define _Proto_EBHostFluxData_H_

#include <cmath>
#include <cstdlib>
#include <memory>
#include "Proto.H"
#include "EBProto.H"

using std::vector;
using std::shared_ptr;
namespace Proto
{
/// EBFluxData that lives on the host
/**
 */
  ///
  template <typename  data_t, unsigned int ncomp>
  class EBHostFluxData
  {
  public:
    ///
    EBHostFluxData()
    {
    }

    virtual ~EBHostFluxData()
    {
    }


    /// sets all face data to the input value
    inline void setVal(data_t a_val)
    {
      m_xflux.setVal(a_val);
      m_yflux.setVal(a_val);
#if DIM==3      
      m_zflux.setVal(a_val);
#endif      
    }

    /// 
   /**
      The input box is cell centered.  Each face data holder has a box = surrounding nodes of the box in its direction.  Calls define function.
    */
    EBHostFluxData(const Box& a_box, const EBGraph& a_graph, bool a_useStack=false)
    {
      define(a_box, a_graph);
    }
    
    /// define all data
    /**
       The input box is cell centered.  Each face data holder has a box = surrounding nodes of the box in its direction.
    */
    inline void   define(const Box& a_box, const EBGraph& a_graph)
    {
      m_box = a_box;
      auto* ptrx = new EBHostData<XFACE, data_t, ncomp>(a_box, a_graph);
      auto* ptry = new EBHostData<YFACE, data_t, ncomp>(a_box, a_graph);
      m_xflux = shared_ptr< EBHostData<XFACE, data_t, ncomp> >(ptrx); 
      m_yflux = shared_ptr< EBHostData<YFACE, data_t, ncomp> >(ptry);
#if DIM==3      
      auto* ptrz = new EBHostData<ZFACE, data_t, ncomp>(a_box, a_graph);
      m_zflux = shared_ptr< EBHostData<ZFACE, data_t, ncomp> >(ptrz);
#endif      

    }
    ///data  on x coordinate faces
    shared_ptr< EBHostData<XFACE, data_t, ncomp> >    m_xflux;
    ///data  on y coordinate faces
    shared_ptr< EBHostData<YFACE, data_t, ncomp> >    m_yflux;
    ///data  on z coordinate faces
    shared_ptr< EBHostData<ZFACE, data_t, ncomp> >    m_zflux;
    //graph stuff needed for norm

    ///
    data_t maxNorm(const Box& a_grid, unsigned int a_icomp) const
    {
      data_t xval = Proto::maxNormFace<XFACE>(*m_xflux, a_grid, a_icomp);
      data_t yval = Proto::maxNormFace<YFACE>(*m_yflux, a_grid, a_icomp);
      data_t maxval = std::max(xval, yval);
#if DIM==3
      data_t zval = Proto::maxNormFace<ZFACE>(*m_zflux, a_grid, a_icomp);
      maxval = std::max(maxval, zval);
#endif
      return maxval;
    }

    int nComp() const
    {
      return ncomp;
    }
    ///for data transfers.  datasize does NOT include sizeof(data_t) (to conform with boxdata::size)
    void dataInfo(data_t*& a_dataPtr, size_t& a_datasize, unsigned int a_whichfab, unsigned int a_idir) const
    {
      if(a_idir == 0)
      {
        m_xflux->dataInfo(a_dataPtr, a_datasize, a_whichfab);
      }
      else if(a_idir == 1)
      {
        m_yflux->dataInfo(a_dataPtr, a_datasize, a_whichfab);
      }
#if DIM==3
      else if (a_idir == 2)
      {
        m_zflux->dataInfo(a_dataPtr, a_datasize, a_whichfab);
      }
#endif
      else
      {
        PROTO_ASSERT(false, "ebhosfluxdata::dataInfo: bogus idir");
      }
    }


    virtual
    void copy(const EBHostFluxData< data_t, ncomp>      &  a_src,
              const Box                                 &  a_srcbox,
              int                                          a_srccomp,
              const Box                                 &  a_dstbox,
              int                                          a_dstcomp,
              int                                          a_numcomp)
    {
      m_xflux->copy(*a_src.m_xflux, a_srcbox, a_srccomp, a_dstbox, a_dstcomp, a_numcomp);
      m_yflux->copy(*a_src.m_yflux, a_srcbox, a_srccomp, a_dstbox, a_dstcomp, a_numcomp);
#if DIM==3      
      m_zflux->copy(*a_src.m_zflux, a_srcbox, a_srccomp, a_dstbox, a_dstcomp, a_numcomp);
#endif      
    }

    ///this DOES include sizeof(T)  
    virtual std::size_t charsize(const Box&   a_bx,
                                 unsigned int a_startcomp,
                                 unsigned int a_numcomps,
                                 bool a_printStuff = false) const
    {
      size_t sumsize = 3*sizeof(size_t);
      sumsize += m_xflux->charsize(a_bx, a_startcomp, a_numcomps);
      sumsize += m_yflux->charsize(a_bx, a_startcomp, a_numcomps);
#if DIM==3      
      sumsize += m_zflux->charsize(a_bx, a_startcomp, a_numcomps);
#endif      
      return sumsize;
      
    }
    virtual void  linearIn(void*                  a_buf,
                           size_t                 a_expectedCharSize,
                           const Box       &      a_bx,
                           unsigned int           a_startcomp,
                           unsigned int           a_numcomps,
                           bool a_printStuff = false)
    {
      PR_assert((a_startcomp + a_numcomps-1) < ncomp);
      char*  charbuf = (char*)a_buf;
      size_t charsiz = 0;
      size_t* longbuf = (size_t *)(charbuf);

      size_t xSize = *longbuf;    longbuf++;
      charbuf += sizeof(size_t);
      charsiz += sizeof(size_t);
      
      size_t ySize = *longbuf;    longbuf++;
      charbuf += sizeof(size_t);
      charsiz += sizeof(size_t);
#if DIM==3      
      size_t zSize = *longbuf;    longbuf++;
      charbuf += sizeof(size_t);
      charsiz += sizeof(size_t);
#endif      

      
      m_xflux->linearIn(charbuf,   xSize, a_bx, a_startcomp, a_numcomps);
      charbuf += xSize;
      charsiz += xSize;
      
      m_yflux->linearIn(charbuf,   ySize, a_bx, a_startcomp, a_numcomps);
      charbuf += ySize;
      charsiz += ySize;
#if DIM==3      
      m_zflux->linearIn(charbuf,   zSize, a_bx, a_startcomp, a_numcomps);
      charbuf += zSize;
      charsiz += zSize;
#endif
    }
    virtual void linearOut(void*               a_buf,
                           size_t              a_expectedCharSize,
                           const Box       &   a_bx,
                           unsigned int        a_startcomp,
                           unsigned int        a_numcomps,
                           bool a_printStuff = false) const
    {
      PR_assert((a_startcomp + a_numcomps-1) < ncomp);
      char* charbuf = (char*)a_buf;
      size_t charsiz = 0;
      
      size_t xSize = m_xflux->charsize(a_bx, a_startcomp, a_numcomps);
      size_t ySize = m_yflux->charsize(a_bx, a_startcomp, a_numcomps);

#if DIM==3      
      size_t zSize = 0;
      zSize = m_zflux->charsize(a_bx, a_startcomp, a_numcomps);
#endif      
      
      size_t* longbuf = (size_t*)(charbuf);
      *longbuf = xSize;      longbuf++;
      charbuf += sizeof(size_t);
      charsiz += sizeof(size_t);
      
      *longbuf = ySize;      longbuf++;
      charbuf += sizeof(size_t);
      charsiz += sizeof(size_t);
      
#if DIM==3      
      *longbuf = zSize;      longbuf++;
      charbuf += sizeof(size_t);
      charsiz += sizeof(size_t);
#endif
        
      m_xflux->linearOut(charbuf,  xSize, a_bx, a_startcomp, a_numcomps);
      charbuf += xSize;
      charsiz += xSize;
      
      m_yflux->linearOut(charbuf,  ySize, a_bx, a_startcomp, a_numcomps);
      charbuf += ySize;
      charsiz += ySize;
#if DIM==3      
      m_zflux->linearOut(charbuf,  zSize, a_bx, a_startcomp, a_numcomps);
      charbuf += zSize;
      charsiz += zSize;
#endif      
    }
    ///
    static int preAllocatable()
    {
      return 2; 
    }

    static int memTypeAllocation()
    {
      return 0;
    }

    Box box() const
    {
      return m_box;
    }
  protected:
    Box m_box;
    EBHostFluxData(const            EBHostFluxData<data_t, ncomp>& a_src);
    EBHostFluxData& operator=(const EBHostFluxData<data_t, ncomp>& a_src);

    

  };

}

#endif
