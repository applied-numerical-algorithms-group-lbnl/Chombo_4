
#ifndef _Proto_EBHostBoxDATA_H_
#define _Proto_EBHostBoxDATA_H_

#include <cmath>
#include <cstdlib>
#include <memory>
#include "Proto.H"
#include "EBProto.H"
#include "Proto_EBBoxData.H" //for aliasing
using std::vector;
using std::shared_ptr;
namespace Proto
{
  ///Holder for EB data over a box
  /**
     This is the general holder for all EB.
     if cent == CELL, this is over cells.
     if cent == BOUNDARY, the data only exists at cut cells
     if cent == XFACE, YFACE, ZFACE, the data is at coordinate face.
     All box inputs are assumed to be cll centered
   */
  template <CENTERING cent, typename  data_t, unsigned int ncomp>
  class EBHostData 
  {
  public:
    static const bool s_stackDefault = false;

    ///
    EBHostData()
    {
    }

    int nComp() const
    {
      return ncomp;
    }
    virtual ~EBHostData()
    {
    }

    inline EBDataLoc indexInfo(const EBIndex<cent>& a_index, int a_ivar) const
    {
      //this ought to work even if the index is outside the domain
      //should just send it to reg fab
      EBDataLoc retval;
      if(m_irrData.contains(a_index))
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

    inline  data_t& operator()(const EBIndex<cent>& a_vof, 
                               unsigned int  a_ivar) 
    {
      EBDataLoc index = indexInfo(a_vof,  a_ivar);
      data_t* retval;
      if(index.m_dataID == 0)
      {
        retval = &(m_regData(a_vof.m_pt, a_ivar));
      }
      else
      {
        retval = &(m_irrData(a_vof, a_ivar));
      }
      return *retval;

    }

    inline  const data_t& operator()(const EBIndex<cent>& a_vof, 
                                     unsigned int  a_ivar)  const
    {
      EBDataLoc index = indexInfo(a_vof,  a_ivar);
      const data_t* retval;
      if(index.m_dataID == 0)
      {
        retval = &(m_regData(a_vof.m_pt, a_ivar));
      }
      else
      {
        retval = &(m_irrData(a_vof, a_ivar));
      }
      return *retval;

    }

    ///needed for Proto::LevelData.  meaningless here becasue we need the graph.
    EBHostData(const Box& a_box, bool a_stackAlloc)
    {
    }


    //aliasing define
    template<unsigned int srccomp>
    void define(EBHostData<cent, data_t, srccomp>& a_input,
                unsigned int                      a_comp)
    {
      m_inputBox = a_input.getRegData().box();
      m_graph = a_input.ebgraph();
      m_irrData.define(a_input.getIrregData(), a_comp);
      m_regData.define(a_input.getRegData()  , a_comp);
    }           


    
    EBHostData(EBBoxData<cent, data_t, ncomp>& a_input)
    {
      define(a_input);
    }
    
    //aliasing define from device data for when the device is just the cpu
    void define(EBBoxData<cent, data_t, ncomp>& a_input)
    {
#ifdef PROTO_CUDA
      PROTO_ASSERT(false, "aliasing constructor called between host and device data holders");
#endif      
      m_inputBox = a_input.getRegData().box();
      m_graph    = a_input.ebgraph();
      m_irrData.define(a_input.getIrregData());
      m_regData.define(a_input.getRegData()  );
    }           

    ///for aggstencil mostly.  datasize does NOT include sizeof(data_t) (to conform with boxdata::size)
    void dataInfo(data_t*& a_dataPtr, size_t& a_datasize, unsigned int a_whichfab) const
    {
      PR_assert((a_whichfab == 0) || (a_whichfab == 1));
      if(a_whichfab == 0)
      {
        a_dataPtr  =  (data_t*)(m_regData.data());
        a_datasize =  m_regData.size();
      }
      else 
      {
        a_dataPtr  =  (data_t*)(m_irrData.data());
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
    EBHostData(const Box& a_box, const EBGraph& a_graph, bool a_useStack = false)
    {

      define(a_box, a_graph, a_useStack);
    }

    inline void define(const Box& a_box, const EBGraph& a_graph, bool a_useStack= false)
    {
      this->m_inputBox = a_box;
      this->m_graph = a_graph;

      this->m_irrData.define(a_box, this->m_graph);
//begin debug
//      Chombo4::pout() << "in box define for EBHostData: " ;
//      Chombo4::pout() << "irrdata.vecsize = " << this->m_irrData.vecsize() << std::endl;
//end debug      
      this->m_regData.define(a_box);  
    }

    HostBoxData<data_t, ncomp>& getRegData()
    {
      return m_regData;
    }
    
    HostIrregData<cent, data_t, ncomp>&  getIrrData()
    {
      return m_irrData;
    }

    inline
    void setVal(data_t val)
    {
      m_regData.setVal(val);
      m_irrData.setVal(val);
    }

    inline
    void copy(const EBHostData<cent, data_t, ncomp>& a_src,
              const Box& a_bxFrom, unsigned int isrc,
              const Box& a_bxTo  , unsigned int idst,
              unsigned int inco)
    {
      {
        m_regData.copy(a_src.m_regData, a_bxFrom, isrc, a_bxTo, idst, inco);
      }
      {
        m_irrData.copy(a_src.m_irrData, a_bxFrom, isrc, a_bxTo, idst, inco);
      }
    }
              
    ///this DOES include sizeof(T)  
    virtual std::size_t charsize(const Box&   a_bx,
                                 unsigned int a_startcomp,
                                 unsigned int a_numcomps,
                                 bool a_printStuff = false) const
    {
      PR_assert((a_startcomp + a_numcomps-1) < ncomp);
      size_t retval = 2*sizeof(size_t); //regsize and irregsize
      size_t regsize = m_regData.charsize(a_bx, a_startcomp, a_numcomps, a_printStuff);
      size_t irrsize = m_irrData.charsize(a_bx, a_startcomp, a_numcomps, a_printStuff);
      size_t sumsize = regsize + irrsize;
      retval += sumsize;
      return retval;
      
    }
    virtual void  linearIn(void*                  a_buf,
                           const size_t a_expectedCharSize,
                           const Box       &      a_bx,
                           unsigned int           a_startcomp,
                           unsigned int           a_numcomps,
                           bool a_printStuff = false)
    {
      PR_assert((a_startcomp + a_numcomps-1) < ncomp);
      if(a_printStuff)
      {
        Chombo4::pout() << "EBHostData::linearIn: entering for box "<< a_bx << endl;
      }
      size_t charsiz = 0;
      char*  charbuf = (char*)a_buf;
      size_t* regSizePtr = (size_t*)(charbuf);
      size_t  regSize = *regSizePtr;

      charbuf  += sizeof(size_t);
      charsiz  += sizeof(size_t);

      size_t* irrSizePtr = (size_t*)charbuf;
      size_t  irrSize = *irrSizePtr;

      charbuf += sizeof(size_t);
      charsiz += sizeof(size_t);
      
      if(a_printStuff)
      {
        size_t totSize = regSize + irrSize + 2*sizeof(size_t);
        Chombo4::pout() << "EBHostData::linearIn: totSize = " << totSize << endl;
        if(totSize != a_expectedCharSize)
        {
          Chombo4::pout() << "EBHostData::linearIn: WARNING: "
                          << "expectedCharSize = " << a_expectedCharSize
                          << ", regSize = " << regSize 
                          << ", irrSize = " << irrSize 
                          << ", totSize = " << totSize << endl;
        }
      }
      
      size_t expectedCharSize = regSize; 
      m_regData.linearIn(charbuf,   regSize, a_bx, a_startcomp, a_numcomps, a_printStuff);
      charbuf += regSize;
      charsiz += regSize;
      if(a_printStuff)
      {
        Chombo4::pout() << "EBHostData::linearIn:going into IrregData::linearIn" << endl;
      }

      expectedCharSize = irrSize; 
      m_irrData.linearIn(charbuf, expectedCharSize,  a_bx, a_startcomp, a_numcomps, a_printStuff);
      charbuf += irrSize;
      charsiz += irrSize;

      if(a_printStuff)
      {
        if(charsiz != a_expectedCharSize)
        {
          Chombo4::pout() << "EBHostData::linearIn:WARNING:"
            << "  charsiz  = " << charsiz 
            << ", expected = " << a_expectedCharSize << endl;
        }     
        Chombo4::pout() << "EBHostData::linearIn:leaving" << endl;
      }

    }
    virtual void linearOut(void*               a_buf,
                           const size_t a_expectedCharSize,
                           const Box       &   a_bx,
                           unsigned int        a_startcomp,
                           unsigned int        a_numcomps,
                           bool a_printStuff = false) const
    {
      PR_assert((a_startcomp + a_numcomps-1) < ncomp);
      if(a_printStuff)
      {
        Chombo4::pout() << "EBHostData::linearOut:entering for box " << a_bx << endl;
      }
      size_t regSize = m_regData.charsize(a_bx, a_startcomp, a_numcomps, a_printStuff);
      size_t irrSize = m_irrData.charsize(a_bx, a_startcomp, a_numcomps, a_printStuff);

      if(a_printStuff)
      {
        size_t totSize = regSize + irrSize + 2*sizeof(size_t);
        Chombo4::pout() << "EBHostData::linearOut: regSize = " << regSize << endl;
        Chombo4::pout() << "EBHostData::linearOut: irrSize = " << irrSize << endl;
        Chombo4::pout() << "EBHostData::linearOut: totSize = " << totSize << endl;
        if(totSize != a_expectedCharSize)
        {
          Chombo4::pout() << "EBHostData::linearOut: WARNING: "
            << "expectedCharSize = "  << a_expectedCharSize 
            << ", actual tot size  = " << totSize << endl;
        }
      }
      char*  charbuf = (char*)(a_buf);
      size_t charsiz = 0;
      
      size_t* regSizePtr = (size_t*)(charbuf);
      *regSizePtr = regSize;
      charbuf += sizeof(size_t);
      charsiz += sizeof(size_t);
      
      size_t* irrSizePtr = (size_t*)(charbuf);
      *irrSizePtr = irrSize;
      charbuf += sizeof(size_t);
      charsiz += sizeof(size_t);
      
      size_t expectedCharSize = regSize; 
      m_regData.linearOut(charbuf,  expectedCharSize, a_bx, a_startcomp, a_numcomps, a_printStuff);
      charbuf += regSize;
      charsiz += regSize;
      
      expectedCharSize = irrSize; 
      m_irrData.linearOut(charbuf, expectedCharSize,  a_bx, a_startcomp, a_numcomps, a_printStuff);
      charbuf += irrSize;
      charsiz += irrSize;
      if(a_printStuff)
      {
        if(charsiz != a_expectedCharSize)
        {
          Chombo4::pout() << "EBHostData::linearOut: WARNING: "
            << "expectedCharSize = "  << a_expectedCharSize 
            << ", actual char size = "  << charsiz << endl;
        }
        Chombo4::pout() << "EBHostData::linearOut: leaving" << endl;
      }
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

  protected:
    EBHostData(const            EBHostData<cent, data_t, ncomp>& a_src);
    EBHostData& operator=(const EBHostData<cent, data_t, ncomp>& a_src);

    Box m_inputBox;
    HostBoxData<        data_t, ncomp>     m_regData;
    HostIrregData<cent, data_t, ncomp>     m_irrData;
    EBGraph                                m_graph;
  };

  ///template specialization for cell
  template <typename  data_t, unsigned int ncomp>
  data_t maxNormCell(EBIndex<CELL>& a_vofmax,
                     const EBHostData<CELL, data_t, ncomp>& a_src,
                     const Box& a_grid, unsigned int a_icomp)
                     
    {
      data_t maxval = 0; 
      using Proto::BoxIterator;
      for(auto bit = a_grid.begin(); bit != a_grid.end(); ++bit)
      {
        Point pt = *bit;
        vector<EBIndex<CELL> > vofspt = a_src.ebgraph().getVoFs(pt);
        for(int ivec = 0 ; ivec < vofspt.size(); ivec++)
        {
          const EBIndex<CELL>& vof = vofspt[ivec];
          data_t datval = a_src(vof, a_icomp);
          data_t absval = std::abs(datval);
          if(absval > maxval)
          {
            maxval = absval;
            a_vofmax = vof;
          }
        }
      }
      return maxval;
    }

  ///template specialization for cell
  template <typename  data_t, unsigned int ncomp>
  data_t maxCell(const EBHostData<CELL, data_t, ncomp>& a_src,
                     const Box& a_grid, unsigned int a_icomp)
  {
    data_t maxval = -1e300;
    using Proto::BoxIterator;
    for(auto bit = a_grid.begin(); bit != a_grid.end(); ++bit)
      {
        Point pt = *bit;
        vector<EBIndex<CELL> > vofspt = a_src.ebgraph().getVoFs(pt);
        for(int ivec = 0 ; ivec < vofspt.size(); ivec++)
          {
            const EBIndex<CELL>& vof = vofspt[ivec];
            data_t datval = a_src(vof, a_icomp);
            if(datval > maxval)
              {
                maxval = datval;
              }
          }
      }
    return maxval;
  }

    ///template specialization for cell
  template <typename  data_t, unsigned int ncomp>
  data_t minCell(const EBHostData<CELL, data_t, ncomp>& a_src,
                     const Box& a_grid, unsigned int a_icomp)
  {
    data_t minval = 1e300;
    using Proto::BoxIterator;
    for(auto bit = a_grid.begin(); bit != a_grid.end(); ++bit)
      {
        Point pt = *bit;
        vector<EBIndex<CELL> > vofspt = a_src.ebgraph().getVoFs(pt);
        for(int ivec = 0 ; ivec < vofspt.size(); ivec++)
          {
            const EBIndex<CELL>& vof = vofspt[ivec];
            data_t datval = a_src(vof, a_icomp);
            if(datval < minval)
              {
                minval = datval;
              }
          }
      }
    return minval;
  }

  ///
  template <CENTERING cent, typename  data_t, unsigned int ncomp>
  data_t maxNormFace(EBIndex<cent>& a_facemax,
                     const EBHostData<cent, data_t, ncomp>& a_src,
                     const Box& a_grid, unsigned int a_icomp)
  {
    data_t maxval = 0; 
    using Proto::BoxIterator;
    for(auto bit = a_grid.begin(); bit != a_grid.end(); ++bit)
    {
      Point pt = *bit;
      vector<EBIndex<CELL> > vofspt = a_src.ebgraph().getVoFs(pt);
      for(int ivec = 0 ; ivec < vofspt.size(); ivec++)
      {
        const EBIndex<CELL>& vof = vofspt[ivec];
        for(SideIterator sit; sit.ok(); ++sit)
        {
          vector<EBIndex<cent> > faces = getFaces<cent>(vof, sit(),  a_src.ebgraph());
          for(unsigned int iface = 0; iface < faces.size(); iface++)
          {
            data_t datval = a_src(faces[iface], a_icomp);
            data_t absval = std::abs(datval);
            if(absval > maxval)
            {
              maxval = absval;
              a_facemax = faces[iface];
            }
          }
        }
      }
    }
    return maxval;
  }


}
#endif
