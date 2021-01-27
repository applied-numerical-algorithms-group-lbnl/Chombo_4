
#ifndef _Proto_EBBoxDATA_H_
#define _Proto_EBBoxDATA_H_

#include <cmath>
#include <cstdlib>
#include <memory>
#include <iomanip>
#include "Proto.H"
#include "EBProto.H"
#include "implem/Proto_IrregData.H"
using std::vector;
using std::shared_ptr;
namespace Proto
{
  struct EBDataLoc
  {
    int    m_dataID;   //which pointer to index
    size_t m_offset;   //distance from start of pointer 
  };

  ///Holder for EB data over a box
  /**
     This is the general holder for all EB.
     if cent == CELL, this is over cells.
     if cent == BOUNDARY, the data only exists at cut cells
     if cent == XFACE, YFACE, ZFACE, the data is at coordinate face.
     All box inputs are assumed to be cll centered
   */
  template <CENTERING cent, typename  data_t, unsigned int ncomp>
  class EBBoxData 
  {
  public:
    static const bool s_stackDefault = false;

    ///
    EBBoxData()
    {
    }

    virtual ~EBBoxData()
    {
    }

    ///
    EBBoxData<cent, data_t, ncomp>&
    operator-=(const EBBoxData<cent, data_t, ncomp>& a_input)
    {
      PROTO_ASSERT((a_input.m_inputBox == m_inputBox), "operator-= box mismatch");
      m_regData -= a_input.m_regData;
      m_irrData -= a_input.m_irrData;

      return *this;
    }

    template<unsigned int Csrc>
    void copy(const EBBoxData<cent, data_t, Csrc> & a_src,
              const Box&                  a_srcBox,
              unsigned int               a_srcComp,
              const Box&                  a_destBox,
              unsigned int               a_destComp,
              unsigned int               a_numcomp)
    {
      PR_TIME("EBBoxData::copy");
      {
        PR_TIME("copy_reg_data");
        m_regData.copy(a_src.m_regData, a_srcBox, a_srcComp, a_destBox, a_destComp, a_numcomp);
      }
      {
        PR_TIME("copy_irr_data");
        m_irrData.copy(a_src.m_irrData, a_srcBox, a_srcComp, a_destBox, a_destComp, a_numcomp);
      }
    }
    ///
    EBBoxData<cent, data_t, ncomp>&
    operator+=(const EBBoxData<cent, data_t, ncomp>& a_input)
    {
      PROTO_ASSERT((a_input.m_inputBox == m_inputBox), "operator+= box mismatch");

      m_regData += a_input.m_regData;
      m_irrData += a_input.m_irrData;

      return *this;
    }
    
    ///
    EBBoxData<cent, data_t, ncomp>&
    operator*=(const EBBoxData<cent, data_t, ncomp>& a_input)
    {
      PROTO_ASSERT((a_input.m_inputBox == m_inputBox), "operator*= box mismatch");
      m_regData *= a_input.m_regData;
      m_irrData *= a_input.m_irrData;

      return *this;
    }


    inline
    void scale(data_t                                a_scale)
    {
      EBBoxData<cent, data_t, ncomp> scalefab(m_inputBox, m_graph);
      scalefab.setVal(a_scale);
      (*this) *= scalefab;
    }

    inline
    void incr(const EBBoxData<cent, data_t, ncomp>& a_input,
              data_t                                a_scale)
    {
      EBBoxData<cent, data_t, ncomp> incramnt(m_inputBox, m_graph);
      incramnt.setVal(0);
      incramnt.scale(a_scale);
      (*this) += incramnt;
    }
    

    ///
    virtual  void setVal(data_t a_val)
    {
      m_regData.setVal(a_val);
      m_irrData.setVal(a_val);
    }

    
    ///will not work on device--debugging use only
    inline  data_t operator()(const EBIndex<cent>& a_vof, 
                              unsigned int  a_ivar) const
    {
      EBDataLoc index = indexInfo(a_vof,  a_ivar);
      data_t retval;
      if(index.m_dataID == 0)
      {
        retval = m_regData(a_vof.m_pt, a_ivar);
      }
      else
      {
        retval = m_irrData(a_vof, a_ivar);
      }
      return retval;

    }

    ///for AggStencil
    inline const data_t* dataPtr(int a_dataType, int a_ivar) const
    {
      const data_t*  retval = nullptr;
      if(a_dataType == 0)
      {
        retval = m_regData.dataPtr(a_ivar);
      }
      else
      {
        retval = m_irrData.dataPtr(a_ivar);
      }
      return retval;
    }


    ///for AggStencil
    inline data_t* dataPtr(int a_dataType, int a_ivar) 
    {
      data_t*  retval = nullptr;
      if(a_dataType == 0)
      {
        retval = m_regData.dataPtr(a_ivar);
      }
      else
      {
        retval = m_irrData.dataPtr(a_ivar);
      }
      return retval;
    }

    ///needed to apply stencils
    BoxData<data_t,ncomp>& getRegData() 
    {
      return m_regData;
    }

    ///needed to apply stencils
    const BoxData<data_t,ncomp>& getRegData() const
    {
      return m_regData;
    }


    ///needed to apply stencils
    IrregData<cent, data_t,ncomp>& getIrregData() 
    {
      return m_irrData;
    }

    ///needed to apply stencils
    const IrregData<cent, data_t, ncomp>& getIrregData() const
    {
      return m_irrData;
    }

    ///needed for Proto::LevelData.  meaningless here becasue we need the graph.
    EBBoxData(const Box& a_box, bool a_stackAlloc)
    {
    }


    //aliasing define
    template<unsigned int srccomp>
    void define(EBBoxData<cent, data_t, srccomp>& a_input,
                unsigned int                      a_comp)
    {
      m_inputBox = a_input.getRegData().box();
      m_graph = a_input.ebgraph();
      m_irrData.define(a_input.getIrregData(), a_comp);
      m_regData.define(a_input.getRegData()  , a_comp);
    }           


    
    ///for data transfers.  datasize does NOT include sizeof(data_t) (to conform with boxdata::size)
    void dataInfo(data_t*& a_dataPtr, size_t& a_datasize, unsigned int a_whichfab)
    {
      PR_assert((a_whichfab == 0) || (a_whichfab == 1));
      if(a_whichfab == 0)
      {
        a_dataPtr  =  m_regData.data();
        a_datasize =  m_regData.size();
      }
      else 
      {
        a_dataPtr  =  m_irrData.data();
        a_datasize =  m_irrData.size();
      }
    }

    /// input box is cell centered---semantic is somewhat different depending on centering
//    inline virtual  void define(const Box& a_box, const EBGraph& a_graph) = 0;

    Box box() const
    {
      return m_regData.box();
    }

    EBGraph ebgraph() const
    {
      return m_graph;
    }

    Box inputBox() const
    {
      return m_inputBox;
    }
    /// 
    /**
         Input box (a_box)  is cell centered. the box is adjusted for centering if it is cent==XFACE YFACE or ZFACE.
         All irregular data (data at cut cells is held in  a sparse holder (IrregData).    Regular data is held in a BoxData.
    */
    EBBoxData(const Box& a_box, const EBGraph& a_graph, bool a_useStack = false)
    {
      define(a_box, a_graph, a_useStack);
    }

    inline void define(const Box& a_box, const EBGraph& a_graph, bool a_useStack= false)
    {
      this->m_inputBox = a_box;
      this->m_graph = a_graph;

      this->m_irrData.define(a_box, this->m_graph);  
      this->m_regData.define(a_box, a_useStack);  //stack allocation turned off for now
    }

    /// for AggStencil -- not really part of the public interface.
    inline EBDataLoc indexInfo(const EBIndex<cent>& a_index, int a_ivar) const
    {
      //this ought to work even if the index is outside the domain
      //should just send it to reg fab
      EBDataLoc retval;
      if(m_irrData.hasIndex(a_index))
      {
        retval.m_dataID   = 1;
        retval.m_offset   = this->m_irrData.index(a_index, a_ivar);
      }
      else
      {
        retval.m_dataID   = 0;
        retval.m_offset   = this->m_regData.index(a_index.m_pt, a_ivar);
      }
      return retval;
    }


// for exchange
  virtual std::size_t charsize(const Box&   a_bx,
                               unsigned int a_startcomp,
                               unsigned int a_numcomps) const
  {
    PR_TIME("EBBoxData::charsize");
    PR_assert((a_startcomp + a_numcomps-1) < ncomp);
    size_t regsize = m_regData.charsize(a_bx, a_startcomp, a_numcomps);
    size_t irrsize = m_irrData.charsize(a_bx, a_startcomp, a_numcomps);
    size_t sumsize = regsize + irrsize;

    return sumsize;
  }


  virtual void  linearIn(void*                  a_buf,
                         const Box       &      a_bx,
                         unsigned int           a_startcomp,
                         unsigned int           a_numcomps)
  {
    PR_TIME("EBBoxData::linearIn");
    PR_assert((a_startcomp + a_numcomps-1) < ncomp);
    char* charbuf = (char*)a_buf;
    m_regData.linearIn(charbuf,   a_bx, a_startcomp, a_numcomps);
    charbuf += m_regData.charsize(a_bx, a_startcomp, a_numcomps);
    m_irrData.linearIn(charbuf,   a_bx, a_startcomp, a_numcomps);
  }

  virtual void linearOut(void*               a_buf,
                         const Box       &   a_bx,
                         unsigned int        a_startcomp,
                         unsigned int        a_numcomps) const
  {
    PR_TIME("EBBoxData::linearOut");
    PR_assert((a_startcomp + a_numcomps-1) < ncomp);
    char* charbuf = (char*)a_buf;
    m_regData.linearOut(charbuf,  a_bx, a_startcomp, a_numcomps);
    charbuf += m_regData.charsize(a_bx, a_startcomp, a_numcomps);
    m_irrData.linearOut(charbuf,  a_bx, a_startcomp, a_numcomps);
  }
  ///
  static int preAllocatable()
  {
    return 1;
  }

  static int memTypeAllocation()
  {
    return 1;
  }
  protected:
    EBBoxData(const            EBBoxData<cent, data_t, ncomp>& a_src);
    EBBoxData& operator=(const EBBoxData<cent, data_t, ncomp>& a_src);

    Box m_inputBox;
    BoxData<        data_t, ncomp>     m_regData;
    IrregData<cent, data_t, ncomp>     m_irrData;
    EBGraph                            m_graph;
  };


  //for debugging --- only works on host
  template< typename data_t>
  inline void 
  maxLocCell(const EBBoxData<CELL, data_t, 1>& data, const Box& a_grid, string a_prefix)
  {
    using std::cout; 
    using std::endl;
    using std::ios;
    cout    << std::setprecision(8)
            << std::setiosflags(ios::showpoint)
            << std::setiosflags(ios::scientific);

    EBGraph graph = data.ebgraph();
    data_t maxval = 0;
    EBIndex<CELL> bigshot;
    for(auto bit = a_grid.begin(); bit != a_grid.end(); ++bit)
    {
      auto pt = *bit;
      vector<EBIndex<CELL> > vofs = graph.getVoFs(pt);
      for(int ivof = 0; ivof < vofs.size(); ivof++)
      {
        auto vof =vofs[ivof];
        data_t thisval = std::abs(data(vof, 0));
        if(thisval > maxval)
        {
          maxval = thisval;
          bigshot = vof;
        }
      }
    }
    cout << a_prefix << "maxval = " << maxval << ", loc = "<< bigshot.m_pt << endl;
  }

  //for debugging --- only works on host
  template<CENTERING cent, typename data_t>
  inline void 
  maxLocFace(const EBBoxData<cent, data_t, 1>& data, const Box& a_grid, string a_prefix)
  {
    using std::cout; 
    using std::endl;
    cout    << std::setprecision(8)
            << std::setiosflags(std::ios::showpoint)
            << std::setiosflags(std::ios::scientific);

    EBGraph graph = data.ebgraph();
    data_t maxval = 0;
    EBIndex<cent> bigshot;
    for(auto bit = a_grid.begin(); bit != a_grid.end(); ++bit)
    {
      auto pt = *bit;
      vector<EBIndex<CELL> > vofs = graph.getVoFs(pt);
      for(int ivof = 0; ivof < vofs.size(); ivof++)
      {
        auto vof =vofs[ivof];
        for(SideIterator sit; sit.ok(); ++sit)
        {
          vector<EBIndex<cent> > faces = getFaces<cent>(vof, sit(), graph);
          for(int ifac = 0; ifac < faces.size(); ifac++)
          {
            auto fac = faces[ifac];
            data_t thisval = std::abs(data(fac, 0));
            if(thisval > maxval)
            {
              maxval = thisval;
              bigshot = fac;
            }
          }
        }
      }
    }
    cout << a_prefix << "maxval = " << maxval << ", loc = "<< bigshot.m_pt << endl;
  }





///version that does not send the point to the function
  template<typename Func, typename... Srcs>
  inline void ebforallInPlace(unsigned long long int a_num_flops_point,
                              const char*            a_timername,
                              const Func & a_F,  Box a_box, Srcs&... a_srcs);

///version that sends the point to the function
  template<typename Func, typename... Srcs>
  inline void ebforallInPlace_i(unsigned long long int a_num_flops_point,
                                const char*            a_timername,
                                const Func & a_F,  Box a_box, Srcs&... a_srcs);

  ///version that only touches cut cells
  template<typename Func, typename... Srcs>
  inline void ebforallIrreg(const char*            a_timername,
                            const Func & a_F,  Box a_box, Srcs&... a_srcs);


  ///version that only touches cut cells
  template<typename Func, typename... Srcs>
  inline void ebforallIrreg_i(const char*            a_timername,
                              const Func & a_F,  Box a_box, Srcs&... a_srcs);

  ///simplest version
  template<typename Func, typename... Srcs>
  inline void ebforall_i(const Box& a_inputBox,
                         const Func & a_F,  Box a_box, Srcs&... a_srcs);

///simplest version
  template<typename Func, typename... Srcs>
  inline void ebforall(const Box& a_inputBox,
                       const Func & a_F,  Box a_box, Srcs&... a_srcs);

}



#include "implem/Proto_EBForAllImplem.H"
#endif

