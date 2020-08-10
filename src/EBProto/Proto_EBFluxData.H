
#ifndef _Proto_EBFluxDATA_H_
#define _Proto_EBFluxDATA_H_

#include <cmath>
#include <cstdlib>
#include <memory>
#include "Proto.H"
#include "EBProto.H"

using std::vector;
using std::shared_ptr;
namespace Proto
{
/// Convenience structure that holds coordinate face data in all directions.
/**
   This is just so I can just declare one thing and have all the fluxes.
 */
  ///
  template <typename  data_t, unsigned int ncomp>
  class EBFluxData
  {
  public:
    ///
    EBFluxData()
    {
    }

    virtual ~EBFluxData()
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
    EBFluxData<data_t, ncomp>&
    operator-=(const EBFluxData<data_t, ncomp>& a_input)
    {
      (*m_xflux) -= (*a_input.m_xflux);
      (*m_yflux) -= (*a_input.m_yflux);
#if DIM==3      
      (*m_zflux) -= (*a_input.m_zflux);
#endif      
      return *this;
    }

    template<unsigned int Csrc>
    void copy(const EBFluxData<data_t, Csrc> & a_src,
              const Box&                  a_srcBox,
              unsigned int               a_srcComp,
              const Box&                  a_destBox,
              unsigned int               a_destComp,
              unsigned int               a_numcomp)
    {
      PROTO_ASSERT(false, "how did I get here");
    }

    ///
    EBFluxData<data_t, ncomp>&
    operator+=(const EBFluxData<data_t, ncomp>& a_input)
    {
      (*m_xflux) += (*a_input.m_xflux);
      (*m_yflux) += (*a_input.m_yflux);
#if DIM==3      
      (*m_zflux) += (*a_input.m_zflux);
#endif      
      return *this;
    }
    
    ///
    EBFluxData<data_t, ncomp>&
    operator*=(const EBFluxData<data_t, ncomp>& a_input)
    {
      (*m_xflux) *= (*a_input.m_xflux);
      (*m_yflux) *= (*a_input.m_yflux);
#if DIM==3      
      (*m_zflux) *= (*a_input.m_zflux);
#endif      
      return *this;
    }


    ///
    EBFluxData<data_t, ncomp>&
    operator*=(data_t a_input)
    {
      (*m_xflux) *= (a_input);
      (*m_yflux) *= (a_input);
#if DIM==3      
      (*m_zflux) *= (a_input);
#endif      
      return *this;
    }

    ///aliasing define
    template <unsigned int srccomp>
    void define(EBFluxData<data_t, srccomp> & a_input,
                unsigned int                  a_comp)
    {
      auto* ptrx = new EBBoxData<XFACE, data_t, ncomp>();
      auto* ptry = new EBBoxData<YFACE, data_t, ncomp>();
      m_xflux = shared_ptr< EBBoxData<XFACE, data_t, ncomp> >(ptrx); 
      m_yflux = shared_ptr< EBBoxData<YFACE, data_t, ncomp> >(ptry);
      ptrx->define(*(a_input.m_xflux), a_comp);
      ptry->define(*(a_input.m_yflux), a_comp);
#if DIM==3      
      auto* ptrz = new EBBoxData<ZFACE, data_t, ncomp>();
      m_zflux = shared_ptr< EBBoxData<ZFACE, data_t, ncomp> >(ptrz);
      ptrz->define(*(a_input.m_zflux), a_comp);
#endif      
      
    }
    
    /// 
   /**
      The input box is cell centered.  Each face data holder has a box = surrounding nodes of the box in its direction.  Calls define function.
    */
    EBFluxData(const Box& a_box, const EBGraph& a_graph, bool a_useStack=false)
    {
      define(a_box, a_graph);
    }
    
    /// define all data
    /**
       The input box is cell centered.  Each face data holder has a box = surrounding nodes of the box in its direction.
    */
    inline void   define(const Box& a_box, const EBGraph& a_graph, bool a_useStack=false)
    {
      auto* ptrx = new EBBoxData<XFACE, data_t, ncomp>(a_box, a_graph, a_useStack);
      auto* ptry = new EBBoxData<YFACE, data_t, ncomp>(a_box, a_graph, a_useStack);
      m_xflux = shared_ptr< EBBoxData<XFACE, data_t, ncomp> >(ptrx); 
      m_yflux = shared_ptr< EBBoxData<YFACE, data_t, ncomp> >(ptry);
#if DIM==3      
      auto* ptrz = new EBBoxData<ZFACE, data_t, ncomp>(a_box, a_graph, a_useStack);
      m_zflux = shared_ptr< EBBoxData<ZFACE, data_t, ncomp> >(ptrz);
#endif      

    }
    ///data  on x coordinate faces
    shared_ptr< EBBoxData<XFACE, data_t, ncomp> >    m_xflux;
    ///data  on y coordinate faces
    shared_ptr< EBBoxData<YFACE, data_t, ncomp> >    m_yflux;
    ///data  on z coordinate faces
    shared_ptr< EBBoxData<ZFACE, data_t, ncomp> >    m_zflux;

    ///for data transfers.  datasize does NOT include sizeof(data_t) (to conform with boxdata::size)
    void dataInfo(data_t*& a_dataPtr, size_t& a_datasize, unsigned int a_whichfab, unsigned int a_idir)
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
      else if(a_idir == 2)
      {
        m_zflux->dataInfo(a_dataPtr, a_datasize, a_whichfab);
      }
#endif
      else
      {
        PROTO_ASSERT(false, "ebhosfluxdata::dataInfo: bogus idir");
      }
    }
    

  protected:
    EBFluxData(const            EBFluxData<data_t, ncomp>& a_src);
    EBFluxData& operator=(const EBFluxData<data_t, ncomp>& a_src);

    

  };

}

#endif
