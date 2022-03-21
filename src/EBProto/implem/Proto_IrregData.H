
#ifndef _Proto_IrregDATA_H_
#define _Proto_IrregDATA_H_

#include <cmath>
#include <cstdlib>
#include <memory>
#include "base/Proto_Box.H"
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

    ///to conform with mpi interface
    int nComp() const
    {
      return ncomp;
    }
    
    /// constructor by copying
    IrregData(Box a_box, data_t* a_data, std::vector<index_t>& a_indicies) : m_boxdefined(true)
    {
      m_inputBox = a_box;
      unsigned int size = a_indicies.size();
      protoMalloc(MEMTYPE_DEFAULT, m_rawPtr,size*sizeof(data_t));
      protoMemcpy(MEMTYPE_DEFAULT, m_rawPtr, a_data, size*sizeof(data_t), protoMemcpyHostToDevice);
      m_data = ::std::shared_ptr<data_t>(m_rawPtr, [](data_t* p){protoFree(MEMTYPE_DEFAULT, p);});
      m_indicies = ::std::shared_ptr<std::vector<index_t>>(new std::vector<index_t>[1], [](std::vector<index_t>* p) { delete[] p;});
      m_indicies->resize(size);
      for(int i = 0; i < size ; i++)
        (*m_indicies)[i] = a_indicies[i];
      defineEBIrregData();
    }

    ///
    virtual ~IrregData()
    {
    }

    ///needed for leveldata
    IrregData(const Box& a_box, bool a_stackAllocation=true)
    {  }

    ///
    inline std::size_t size() const 
    {
      if(m_boxdefined) return m_indicies->size()*ncomp;
      else return 0; 
    };

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
      if((m_inputBox == a_input.m_inputBox))
      {
        unsigned int nsize = size();
        if(nsize == 0) return *this;
        int blockSize = 512;
        int nBlocks = (nsize + blockSize - 1) / blockSize;
        //protoLaunchKernelMemAsync((opKernel<data_t,BoxDataOp::Add>), nBlocks, blockSize, 0, protoGetCurrentStream, m_rawPtr, a_input.m_rawPtr, nsize);
        protoLaunchKernelMemAsyncT<MEMTYPE_DEFAULT, opKernel<data_t, BoxDataOp::Add>>(nBlocks, blockSize, 0, protoGetCurrentStream, m_rawPtr, a_input.m_rawPtr, nsize);
      }
      else
      {
        PROTO_ASSERT(false, "operator+= not implemented for irreg datas of different sizes");
      }
      return *this;
    }


    IrregData<cent, data_t, ncomp>& 
    operator-=(const IrregData<cent, data_t, ncomp> & a_input)
    {
      if(m_boxdefined && (m_inputBox == a_input.m_inputBox))
      {
        unsigned int nsize = size();
        if(nsize == 0) return *this;
        int blockSize = 512;
        int nBlocks = (nsize + blockSize - 1) / blockSize;
        //protoLaunchKernelMemAsync((opKernel<data_t,BoxDataOp::Subtract>), nBlocks, blockSize, 0, protoGetCurrentStream, m_rawPtr, a_input.m_rawPtr, nsize);
        protoLaunchKernelMemAsyncT<MEMTYPE_DEFAULT, opKernel<data_t,BoxDataOp::Subtract>>(nBlocks, blockSize, 0, protoGetCurrentStream, m_rawPtr, a_input.m_rawPtr, nsize);
      }
      else
      {
        PROTO_ASSERT(false, "operator+= not implemented for irreg datas of different sizes");
      }
      return *this;
    }


    IrregData<cent, data_t, ncomp>& 
    operator*=(const IrregData<cent, data_t, ncomp> & a_input)
    {
      if(m_boxdefined && (m_inputBox == a_input.m_inputBox))
      {
        unsigned int nsize = size();
        if(nsize == 0) return *this;
        int blockSize = 512;
        int nBlocks = (nsize + blockSize - 1) / blockSize;
        //protoLaunchKernelMemAsync((opKernel<data_t,BoxDataOp::Multiply>), nBlocks, blockSize, 0, protoGetCurrentStream, m_rawPtr, a_input.m_rawPtr, nsize);
        protoLaunchKernelMemAsyncT<MEMTYPE_DEFAULT, opKernel<data_t, BoxDataOp::Multiply>>(nBlocks, blockSize, 0, protoGetCurrentStream, m_rawPtr, a_input.m_rawPtr, nsize);
      }
      else
      {
        PROTO_ASSERT(false, "operator+= not implemented for irreg datas of different sizes");
      }
      return *this;
    }
    
    /// debugging only--will not work on the device
    data_t operator() (const index_t& a_vof, int a_comp) const
    {
      size_t ind = index(a_vof, a_comp);
      return m_rawPtr[ind];
    }

    vector<EBIndex<cent> > getIndices() const
    {
      return *m_indicies;
    }

    unsigned int vecsize() const
    {
      if(m_boxdefined) return m_indicies->size();
      else return 0;
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
      m_boxdefined = true;
      m_indicies       = a_input.getIndicies();
      m_inputBox       = a_input.inputBox();
      m_data           = a_input.getData();
      m_rawPtr         = a_input.dataPtr(a_comp);
      defineEBIrregData();
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
      defineEBIrregData();
    }

    inline void defineEBIrregData()
    {
      vector< EBIrregStruct<cent, data_t, ncomp> > hostvec = getEBIrregStruct(*m_indicies, *this);
      size_t memsize = hostvec.size()*sizeof(EBIrregStruct<cent, data_t, ncomp>);
      protoMalloc(MEMTYPE_DEFAULT, m_EBIrregDataPtr,memsize);
      protoMemcpy(MEMTYPE_DEFAULT, m_EBIrregDataPtr, hostvec.data(), memsize, protoMemcpyHostToDevice);
    }

    inline EBIrregStruct<cent, data_t, ncomp>* getEBIrregDataPtr()
    {
      return m_EBIrregDataPtr;
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
  
  protected:
    inline void define( const shared_ptr< vector< index_t > >& a_indicies);

    /// for fast forall
    bool m_boxdefined;
    Box m_inputBox; //exists if defined as part of a boxdata
    shared_ptr<vector<  index_t > >           m_indicies = nullptr;
    shared_ptr<data_t>                            m_data = nullptr;
    data_t*                                     m_rawPtr = nullptr;
    EBIrregStruct<cent, data_t, ncomp>* m_EBIrregDataPtr = nullptr;
    IrregData(const IrregData<cent, data_t,ncomp>& a_src);
    IrregData& operator=(const IrregData<cent, data_t,ncomp>& a_src);
  };

}
#include "Proto_IrregDataImplem.H"

#endif
