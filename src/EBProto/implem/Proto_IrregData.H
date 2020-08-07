
#ifndef _Proto_IrregDATA_H_
#define _Proto_IrregDATA_H_

#include <cmath>
#include <cstdlib>
#include <memory>
#include "Proto_Box.H"
#include "EBProto.H"
#include "Proto_EBIrregStruct.H"

using std::vector;
using std::shared_ptr;
namespace Proto
{

///
/**
 */
  template <CENTERING cent, class data_t,  unsigned int ncomp>
  class IrregData
  {
  public:
    typedef EBIndex<cent> index_t;
    ///
    IrregData()
    {
      m_boxdefined = false;
    }

    ///
    virtual ~IrregData()
    {
#ifdef PROTO_CUDA
      protoFree(m_uglystructDevPtr);      
#endif
    }

    ///needed for leveldata
    IrregData(const Box& a_box, bool a_stackAllocation=true)
    {  }

    ///
    inline std::size_t size() const {return m_indicies->size()*ncomp;};

    ///
    inline bool defined() const {return bool(m_data);};
   

    ///
    inline void setVal(const data_t& a_val);
  

    ///
    data_t* dataPtr(unsigned int a_ivar = 0)
    {
      PR_assert(a_ivar < ncomp);
      size_t ioff = a_ivar*m_indicies->size();
      return &(m_rawPtr[ioff]);
    }

    ///
    const data_t* dataPtr(unsigned int a_ivar = 0) const
    {
      PR_assert(a_ivar < ncomp);
      size_t ioff = a_ivar*m_indicies->size();
      return &(m_rawPtr[ioff]);
    }

    data_t* data()
    {
      return m_rawPtr;
    }


    const data_t* data() const
    {
      return m_rawPtr;
    }
    
    ///
    inline size_t  index(const EBIndex<cent>& a_vof, const int& a_ivar) const
    {
      size_t ivec = 0;
      bool found = false;
      for(unsigned int i = 0; i < m_indicies->size(); i++)
      {        
        if(a_vof == (*m_indicies)[i])
        {
          ivec = i;
          found = true;
          break;
        }
      }
      if(!found)
      {
        std::cout << "index<cell> vof not found" << endl;
        abort();
      }
      PROTO_ASSERT(found, "vof not found");

      ivec += a_ivar*m_indicies->size();

      if(ivec >= size())
      {
        std::cout << "index<cell> offset too big" << endl;
        abort();
      }

      return ivec;
    }

    IrregData<cent, data_t, ncomp>& 
    operator+=(const IrregData<cent, data_t, ncomp> & a_input)
    {
#ifdef PROTO_CUDA
      if(m_boxdefined && (m_inputBox == a_input.m_inputBox))
      {
        data_t* castPtrInp = (data_t*)(a_input.m_rawPtr);
        thrust::device_ptr<data_t> devptrLoc = thrust::device_pointer_cast(m_rawPtr);
        thrust::device_ptr<data_t> devptrInp = thrust::device_pointer_cast(castPtrInp);
        unsigned int nsize = size();
        thrust::transform(thrust_proto::par.on(DisjointBoxLayout::getCurrentStream()),
                          devptrLoc, devptrLoc+nsize, devptrInp, devptrLoc, thrust::plus<data_t>());
        
      }
      else
      {
        PROTO_ASSERT(false, "operator+= not implemented for irreg datas of different sizes");
      }
#else
      for(int i = 0; i < m_indicies->size(); i++)
      {
        for(int icomp = 0; icomp < ncomp; icomp++)
        {
          size_t dstind =             index((*m_indicies)[i], icomp);
          size_t srcind =     a_input.index((*m_indicies)[i], icomp);
          m_rawPtr[dstind] += a_input.m_rawPtr[srcind];
        }
      }

#endif
      return *this;
    }


    IrregData<cent, data_t, ncomp>& 
    operator-=(const IrregData<cent, data_t, ncomp> & a_input)
    {
#ifdef PROTO_CUDA
      if(m_boxdefined && (m_inputBox == a_input.m_inputBox))
      {
        data_t* castPtrInp = (data_t*)(a_input.m_rawPtr);
        thrust::device_ptr<data_t> devptrLoc = thrust::device_pointer_cast(m_rawPtr);
        thrust::device_ptr<data_t> devptrInp = thrust::device_pointer_cast(castPtrInp);
        unsigned int nsize = size();
        thrust::transform(thrust_proto::par.on(DisjointBoxLayout::getCurrentStream()),
                          devptrLoc, devptrLoc+nsize, devptrInp, devptrLoc, thrust::plus<data_t>());
        
      }
      else
      {
        PROTO_ASSERT(false, "operator+= not implemented for irreg datas of different sizes");
      }
#else
      for(int i = 0; i < m_indicies->size(); i++)
      {
        for(int icomp = 0; icomp < ncomp; icomp++)
        {
          size_t dstind =             index((*m_indicies)[i], icomp);
          size_t srcind =     a_input.index((*m_indicies)[i], icomp);
          m_rawPtr[dstind] -= a_input.m_rawPtr[srcind];
        }
      }

#endif
      return *this;
    }


    IrregData<cent, data_t, ncomp>& 
    operator*=(const IrregData<cent, data_t, ncomp> & a_input)
    {
#ifdef PROTO_CUDA
      if(m_boxdefined && (m_inputBox == a_input.m_inputBox))
      {
        data_t* castPtrInp = (data_t*)(a_input.m_rawPtr);
        thrust::device_ptr<data_t> devptrLoc = thrust::device_pointer_cast(m_rawPtr);
        thrust::device_ptr<data_t> devptrInp = thrust::device_pointer_cast(castPtrInp);
        unsigned int nsize = size();
        thrust::transform(thrust_proto::par.on(DisjointBoxLayout::getCurrentStream()),
                          devptrLoc, devptrLoc+nsize, devptrInp, devptrLoc, thrust::multiplies<data_t>());
        
      }
      else
      {
        PROTO_ASSERT(false, "operator+= not implemented for irreg datas of different sizes");
      }
#else
      for(int i = 0; i < m_indicies->size(); i++)
      {
        for(int icomp = 0; icomp < ncomp; icomp++)
        {
          size_t dstind =             index((*m_indicies)[i], icomp);
          size_t srcind =     a_input.index((*m_indicies)[i], icomp);
          m_rawPtr[dstind] *= a_input.m_rawPtr[srcind];
        }
      }

#endif
      return *this;
    }
    
    /// debugging only--will not work on the device
    data_t operator() (const index_t& a_vof, int a_comp) const
    {
      size_t ind = index(a_vof, a_comp);
      return m_rawPtr[ind];
    }

    vector<EBIndex<cent> > getIndices(const Box& a_box) const
    {
      vector<EBIndex<cent> > retval;
      CENTERING tcent = cent;
      Box compBox = a_box;
      if(tcent == XFACE)
      {
        compBox = a_box.grow(0, 1);
      }
      else if(tcent == YFACE)
      {
        compBox = a_box.grow(1, 1);
      }
      else if(tcent == ZFACE)
      {
        compBox = a_box.grow(1, 1);
      }
      
      for(int i = 0; i < m_indicies->size(); i++)
      {
        if(a_box.contains((*m_indicies)[i].m_pt))
        {
          retval.push_back((*m_indicies)[i]);
        }
      }
      return retval;
    }

    unsigned int vecsize() const
    {
      return m_indicies->size();
    }

    inline Box inputBox() const
    {
      return m_inputBox;
    }
    
    //aliasing define
    template<unsigned int srccomp>
    void define(IrregData<cent, data_t, srccomp>& a_input,
                unsigned int                    & a_comp)
    {
      m_indicies = a_input.getIndicies();
      m_inputBox = a_input.inputBox();
      m_data     = a_input.getData();
      m_rawPtr   = a_input.dataPtr(a_comp);
    }           

    //do not use this
    shared_ptr<vector<  index_t > >  getIndicies() const
    {
      return m_indicies;
    }

    //do not use this
    shared_ptr<data_t> getData() const
    {
      return  m_data;
    }

    inline void define(const Box& a_box,  const EBGraph& a_graph)
    {
      this->m_boxdefined = true;
      this->m_inputBox   = a_box;
      Box interBox = a_box & a_graph.getDomain();
      shared_ptr<vector<EBIndex<cent> > > indicies(new vector< EBIndex<cent> >() );
      *indicies = getIrregIndicies<cent>(interBox, a_graph);

      define(indicies);
    }

    inline bool hasIndex(const EBIndex<cent>& a_ind) const
    {
      bool found = false;
      //this should probably be a smarter search.   maybe sort them?
      for(int ivec = 0;(ivec < m_indicies->size() && !found); ivec++)
      {
        if((*m_indicies)[ivec] == a_ind)
        {
          found = true;
        }
      }

      return found;
    }

    ///
    shared_ptr< vector< EBIrregStruct<cent, data_t, ncomp> > > getUglyStructPtr()
    {
      return m_uglystruct;
    }
#ifdef PROTO_CUDA
    EBIrregStruct<cent, data_t, ncomp>*  getUglyStructDevicePtr()
    {
      return m_uglystructDevPtr;
    }
#endif
  protected:
    inline void define( const shared_ptr< vector< index_t > >& a_indicies);

    /// for fast forall
    shared_ptr< vector< EBIrregStruct<cent, data_t, ncomp> > > m_uglystruct;
    void makeUglyStruct();
#ifdef PROTO_CUDA
    EBIrregStruct<cent, data_t, ncomp>*  m_uglystructDevPtr;
#endif
    bool m_boxdefined;
    Box m_inputBox; //exists if defined as part of a boxdata
    shared_ptr<vector<  index_t > >        m_indicies;
    shared_ptr<data_t>         m_data;
    data_t*                    m_rawPtr;
    IrregData(const IrregData<cent, data_t,ncomp>& a_src);
    IrregData& operator=(const IrregData<cent, data_t,ncomp>& a_src);


  };

}
#include "Proto_IrregDataImplem.H"

#endif