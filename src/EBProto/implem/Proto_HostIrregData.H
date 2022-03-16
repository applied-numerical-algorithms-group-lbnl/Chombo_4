
#ifndef _HostIrregData_H_
#define _HostIrregData_H_


#include <vector>
#include <set>
#include "Proto_EBGraph.H"
#include "EBProto.H"
#include "Proto_EBIrregStruct.H"
#include "Proto_Serialization.H"
namespace Proto
{
///
/**
*/
///
/**
 */
  template <CENTERING cent, class data_t,  unsigned int ncomp>
  class HostIrregData
  {
  public:

    ///
    HostIrregData()
    {
    }

    ///
    virtual ~HostIrregData()
    {
    }


    ///needed for leveldata
    HostIrregData(const Box& a_box, bool a_stackAllocation=true)
    {  }

    ///
    inline std::size_t size() const {return m_indexVec->size()*ncomp;};

    ///
    inline bool defined() const {return bool(m_data);};
   

    ///
    data_t* dataPtr(unsigned int a_ivar = 0)
    {
      PR_assert(a_ivar < ncomp);
      size_t ioff = a_ivar*m_indexVec->size();
      return &(m_rawPtr[ioff]);
    }

    ///
    const data_t* dataPtr(unsigned int a_ivar = 0) const
    {
      PR_assert(a_ivar < ncomp);
      size_t ioff = a_ivar*m_indexVec->size();
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
      for(unsigned int i = 0; i < m_indexVec->size(); i++)
      {        
        if(a_vof == (*m_indexVec)[i])
        {
          ivec = i;
          found = true;
          break;
        }
      }
      if(!found)
      {
        std::cout << "index<cell> vof not found" << std::endl;
        abort();
      }
      PROTO_ASSERT(found, "vof not found");

      ivec += a_ivar*m_indexVec->size();

      if(ivec >= size())
      {
        std::cout << "index<cell> offset too big" << std::endl;
        abort();
      }

      return ivec;
    }

    data_t& operator() (const EBIndex<cent>& a_vof, int a_comp)
    {
      size_t ind = index(a_vof, a_comp);
      return m_rawPtr[ind];
    }
    
    /// 
    const data_t& operator() (const EBIndex<cent>& a_vof, int a_comp) const
    {
      size_t ind = index(a_vof, a_comp);
      return m_rawPtr[ind];
    }
    


    unsigned int vecsize() const
    {
      return m_indexVec->size();
    }

    //aliasing define
    template<unsigned int srccomp>
    void define(HostIrregData<cent, data_t, srccomp>& a_input,
                unsigned int                        & a_comp)
    {
      m_indexVec = a_input.getIndicies();
      m_indexSet = a_input.getIndexSet();
      m_inputBox = a_input.inputBox();
      m_data     = a_input.getData();
      m_rawPtr   = a_input.dataPtr(a_comp);
    }           

    Box inputBox() const
    {
      return m_inputBox;
    }
    //do not use this
    shared_ptr<vector<  EBIndex<cent> > >  getIndicies() const
    {
      return m_indexVec;
    }

    //do not use this
    shared_ptr<data_t> getData() const
    {
      return  m_data;
    }

    inline void define(const Box& a_box,  const EBGraph& a_graph)
    {
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
      for(int ivec = 0;(ivec < m_indexVec->size() && !found); ivec++)
      {
        if((*m_indexVec)[ivec] == a_ind)
        {
          found = true;
        }
      }

      return found;
    }
    
    bool contains(const EBIndex<cent>& a_index) const
    {
      typename std::set<EBIndex<cent> >::iterator it = m_indexSet->find(a_index);
      return (it != m_indexSet->end());
    }

    void setVal(data_t val)
    {
      for(int ivec = 0;ivec < m_indexVec->size() ; ivec++)
      {
        for(int icomp = 0; icomp < ncomp; icomp++)
        {
          (*this)((*m_indexVec)[ivec], icomp) = val;
        }
      }
    }

    bool copyContains(const EBIndex<cent>& a_ind, const Box & a_srcbox) const
    {
      CENTERING thiscent = cent;

      bool in_box = false;
      if((thiscent == CELL) || (thiscent == BOUNDARY))
      {        
        in_box = a_srcbox.contains(a_ind.m_pt);
      }
      else if(thiscent == XFACE)
      {
        Point hipt = a_ind.m_pt;
        Point lopt = a_ind.m_pt - Point::Basis(0, 1);
        in_box = (a_srcbox.contains(hipt) || a_srcbox.contains(lopt));
      }
      else if(thiscent ==  YFACE)
      {
        Point hipt = a_ind.m_pt;
        Point lopt = a_ind.m_pt - Point::Basis(1, 1);
        in_box = (a_srcbox.contains(lopt) || a_srcbox.contains(hipt));
      }
#if DIM==3      
      else if(thiscent ==  ZFACE)
      {
        Point hipt = a_ind.m_pt;
        Point lopt = a_ind.m_pt - Point::Basis(2, 1);
        in_box = (a_srcbox.contains(lopt) || a_srcbox.contains(hipt));
      }
#endif
      else
      {
        PROTO_ASSERT(false, "unaccounted-for centering");
      }
      //now check it is in the actual list
      bool found = this->contains(a_ind);


      return (in_box && found);
    }
    ///
    vector<EBIndex<cent>> getIndiciesInBox(const Box & a_bx) const
    {
      vector<EBIndex<cent>> indicies;

 
      for (unsigned int i=0; i<m_indexVec->size(); ++i)
      {
        const EBIndex<cent> & ind = (*m_indexVec)[i];
        if(this->copyContains(ind, a_bx))
          {
            indicies.push_back(ind);
          } 
      }
      return indicies;
    }

    ///this DOES include sizeof(T)  
    std::size_t charsize(const Box&   a_bx,
                         unsigned int a_startcomp,
                         unsigned int a_numcomps,
                         bool a_printStuff = false) const
    {
      PR_TIME("HostIrregData::charsize");
      using ProtoSer::linearCharBufSize;
      size_t retval = linearCharBufSize(m_inputBox);
      retval += sizeof(int); //size of the vector;
      const auto& data     = *this;
      const auto& indicies = *m_indexVec;
      for(int ipt = 0; ipt < indicies.size(); ipt++)
      {
        retval += linearCharBufSize(indicies[ipt]);
      }
      
      for(int icomp = 0; icomp < ncomp; icomp++)
      {
        for(int ipt = 0; ipt < indicies.size(); ipt++)
        {
          const auto& pt = indicies[ipt];
          retval += linearCharBufSize(data(pt, icomp));
        }
      }
      if(a_printStuff)
      {
        Chombo4::pout() << "hostirregdata::charsize for box " << a_bx << " returning " << retval << endl;
      }
      return retval;
    }

    /// recreates this object into a stream of bytes
    inline
    void energize(void*  a_buf,
                  bool   a_printStuff = false) const
    {
      PR_TIME("HostIrregData::energize");
      using ProtoSer::linearCharBufSize;
      using ProtoSer::linearOutToCharBuf;
      size_t charsiz = 0;
      char*  charbuf = (char *) a_buf;
      //input box
      linearOutToCharBuf(charbuf, m_inputBox);
      charbuf += linearCharBufSize(m_inputBox);
      charsiz += linearCharBufSize(m_inputBox);

      //size of the vector;
      int* intbuf  = (int*) charbuf;
      int vecsize = m_indexVec->size();
      *intbuf = vecsize;
      charsiz += sizeof(int); 
      charbuf += sizeof(int);
      auto& indicies = *m_indexVec;

      //output the indicies
      for(int ipt = 0; ipt < indicies.size(); ipt++)
      {
        linearOutToCharBuf(charbuf, indicies[ipt]);
        charsiz += linearCharBufSize(indicies[ipt]);
        charbuf += linearCharBufSize(indicies[ipt]);
      }
      
      const auto& data     = *this;
      //now output the actual data
      for(int icomp = 0; icomp < ncomp; icomp++)
      {
        for(int ipt = 0; ipt < indicies.size(); ipt++)
        {
          auto pt = indicies[ipt];
          const auto & dataval = data(pt, icomp);
          linearOutToCharBuf(charbuf, dataval);
          charsiz += linearCharBufSize(data(pt, icomp));
          charbuf += linearCharBufSize(data(pt, icomp));
        }
      }
    }

   ///creates this object from the aforementioned stream of bytes
    HostIrregData(void *a_buf,
                  bool a_printStuff = false)
    {
      PR_TIME("HostIrregData::char_constructor");
      if(a_printStuff)
      {
        Chombo4::pout() << "entering hostirregdata buffer constructor"   << endl;
      }
      using ProtoSer::linearInFromCharBuf;
      using ProtoSer::linearCharBufSize;
      size_t charsiz = 0;
      char*  charbuf = (char *) a_buf;
      //input box
      linearInFromCharBuf(m_inputBox, charbuf);
      charbuf += linearCharBufSize(m_inputBox);
      charsiz += linearCharBufSize(m_inputBox);

      //size of the vector;
      int* intbuf  = (int*) charbuf;
      int vecsize = *intbuf;
      charsiz += sizeof(int); 
      charbuf += sizeof(int);
      shared_ptr<vector<EBIndex<cent> > > indicies_sp(new vector<EBIndex<cent> >(vecsize));
      auto& indicies = *indicies_sp;

      //input the indicies
      for(int ipt = 0; ipt < indicies.size(); ipt++)
      {
        linearInFromCharBuf(indicies[ipt], charbuf);
        charsiz += linearCharBufSize(indicies[ipt]);
        charbuf += linearCharBufSize(indicies[ipt]);
      }

      //define the data pointer and all that.
      //set is filled below
      define(indicies_sp);
      
      const auto& data     = *this;
      //now input the actual data
      for(int icomp = 0; icomp < ncomp; icomp++)
      {
        for(int ipt = 0; ipt < indicies.size(); ipt++)
        {
          auto pt = indicies[ipt];
          data_t dataval;
          linearInFromCharBuf(dataval, charbuf);
          size_t ind = index(pt,  icomp);
          m_rawPtr[ind] = dataval;
          charsiz += linearCharBufSize(data(pt, icomp));
          charbuf += linearCharBufSize(data(pt, icomp));
        }
      }
      if(a_printStuff)
      {
        Chombo4::pout() << "leaving hostirregdata buffer constructor charsiz=" << charsiz   << endl;
      }

    }
    
    ///
    inline
    void linearOut(void*               a_buf,
                   const Box       &   a_bx,
                   unsigned int        a_startcomp,
                   unsigned int        a_numcomps,
                   bool a_printStuff = false) const
    {
      PR_TIME("HostIrregData::linearOut");
      energize(a_buf, a_printStuff);
    }
    
    ///
    inline
    void  linearIn(void*                  a_buf,
                   const Box       &      a_bx,
                   unsigned int           a_startcomp,
                   unsigned int           a_numcomps,
                   bool a_printStuff = false)
    {
      PR_TIME("HostIrregData::linearIn");
      
      HostIrregData<cent, data_t, ncomp> src(a_buf, a_printStuff);
      copy(src, a_bx, a_startcomp, a_bx, a_startcomp, a_numcomps, a_printStuff);
      
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
    
    inline
    void copy(const HostIrregData<cent, data_t, ncomp>& a_src,
              const Box& a_bxFrom, unsigned int a_isrc,
              const Box& a_bxTo  , unsigned int a_idst,
              unsigned int a_inco, bool a_printStuff = false)
    {
      PR_TIME("HostIrregData::copy");
      Box interBox = a_bxFrom & a_bxTo;

      for(int iind = 0; iind < m_indexVec->size(); iind++)
      {
        auto& vof = (*m_indexVec)[iind];
        if((interBox.contains(vof.m_pt)) && (a_src.copyContains(vof, interBox)))
        {
          for(unsigned int ico = 0; ico < a_inco; ico++)
          {
            unsigned int isrc = a_isrc + ico;
            unsigned int idst = a_idst + ico;
            (*this)(vof, idst) = a_src(vof, isrc);
          }
        }
      }
    }


    //just used for debugging
    inline void define(const Box& a_box,
                       const shared_ptr< vector< EBIndex<cent> > >& a_indicies)
    {
      this->m_inputBox   = a_box;
      define(a_indicies);
    }
  protected:
    ///
    inline void define( const shared_ptr< vector< EBIndex<cent> > >& a_indexVec)
    {
      m_indexVec = a_indexVec;
      std::size_t count = this->size();
      m_data=::std::shared_ptr<data_t>(new data_t[count], [](data_t* p) { delete[] p;});
      m_rawPtr = m_data.get();
      fillSetFromVector();
    }

    inline void fillSetFromVector() 
    {
      std::shared_ptr<std::set< EBIndex<cent> > > blerg( new std::set< EBIndex<cent> >());
      m_indexSet =  blerg;
      for(int ivec = 0; ivec < m_indexVec->size(); ivec++)
      {
        m_indexSet->insert((*m_indexVec)[ivec]);
      }
    }


    Box m_inputBox; //exists if defined as part of a boxdata
    std::shared_ptr<std::vector<  EBIndex<cent> > >       m_indexVec;
    std::shared_ptr<std::set<     EBIndex<cent> > >       m_indexSet;
    std::shared_ptr<data_t>                               m_data;
    data_t*                                               m_rawPtr;

    ///insert usual reasons here
    HostIrregData(const HostIrregData<cent, data_t,ncomp>& a_src);
    HostIrregData& operator=(const HostIrregData  <cent, data_t,ncomp>& a_src);


  };

}

#endif
