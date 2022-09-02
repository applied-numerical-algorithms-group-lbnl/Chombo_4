
#ifndef _HostIrregData_H_
#define _HostIrregData_H_

using std::vector;
using std::shared_ptr;
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
    typedef EBIndex<cent> index_t;
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
    inline std::size_t size() const {return m_indicies->size()*ncomp;};

    ///
    inline bool defined() const {return bool(m_data);};
   

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
        std::cout << "index<cell> vof not found" << std::endl;
        abort();
      }
      PROTO_ASSERT(found, "vof not found");

      ivec += a_ivar*m_indicies->size();

      if(ivec >= size())
      {
        std::cout << "index<cell> offset too big" << std::endl;
        abort();
      }

      return ivec;
    }

    data_t& operator() (const index_t& a_vof, int a_comp)
    {
      size_t ind = index(a_vof, a_comp);
      return m_rawPtr[ind];
    }
    
    /// 
    const data_t& operator() (const index_t& a_vof, int a_comp) const
    {
      size_t ind = index(a_vof, a_comp);
      return m_rawPtr[ind];
    }
    


    unsigned int vecsize() const
    {
      return m_indicies->size();
    }

    //aliasing define
    template<unsigned int srccomp>
    void define(HostIrregData<cent, data_t, srccomp>& a_input,
                unsigned int                        & a_comp)
    {
      m_indicies = a_input.getIndicies();
      m_inputBox = a_input.inputBox();
      m_data     = a_input.getData();
      m_rawPtr   = a_input.dataPtr(a_comp);
    }           

    Box inputBox() const
    {
      return m_inputBox;
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

    inline void printDebugStuff(const EBIndex<cent>& a_index, const EBGraph& a_graph) const
    {
#if DIM==3
      
      using std::cout;
      using std::endl;
      if(this->hasIndex(a_index))
      {
        cout << "this fab DOES have the index" << endl;
      }
      else 
      {
        cout << "this fab DOES NOT have the index" << endl;
      }
      EBIndex<CELL>  voflo = a_index.getVoF(Side::Lo);
      EBIndex<CELL>  vofhi = a_index.getVoF(Side::Hi);

      if(a_graph.isIrregular(voflo.m_pt))
      {
        cout << "graph thinks " << voflo.m_pt << " is IRREGULAR" << endl;
      }                                            
      else if(a_graph.isCovered(voflo.m_pt))       
      {                                            
        cout << "graph thinks " << voflo.m_pt << " is COVERED  " << endl;
      }                                            
      else if(a_graph.isRegular(voflo.m_pt))       
      {                                            
        cout << "graph thinks " << voflo.m_pt << " is REGULAR  " << endl;
      }                                            
                                                   
      if(a_graph.isIrregular(vofhi.m_pt))          
      {                                            
        cout << "graph thinks " << vofhi.m_pt << " is IRREGULAR" << endl;
      }                                            
      else if(a_graph.isCovered(vofhi.m_pt))       
      {                                            
        cout << "graph thinks " << vofhi.m_pt << " is COVERED  " << endl;
      }                                            
      else if(a_graph.isRegular(vofhi.m_pt))       
      {                                            
        cout << "graph thinks " << vofhi.m_pt << " is REGULAR  " << endl;
      }
      vector<EBIndex<YFACE> > faces = a_graph.getYFaces(vofhi, Side::Lo);
      if(faces.size() > 0)
      {
        cout << "graph thinks there IS     a face between the two vofs" << endl;
      }
      else
      {
        cout << "graph thinks there IS NOT a face between the two vofs" << endl;
      }
#endif 
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
      for(int ivec = 0;(ivec < m_indicies->size() && !found); ivec++)
      {
        if((*m_indicies)[ivec] == a_ind)
        {
          found = true;
        }
      }

      return found;
    }
    
    bool contains(const EBIndex<cent>& a_index) const
    {
      return hasIndex(a_index);
    }

    void setVal(data_t val)
    {
      for(int ivec = 0;ivec < m_indicies->size() ; ivec++)
      {
        for(int icomp = 0; icomp < ncomp; icomp++)
        {
          (*this)((*m_indicies)[ivec], icomp) = val;
        }
      }
    }

    bool copyContains(const index_t& a_ind, const Box & a_srcbox) const
    {
      CENTERING thiscent = cent;

      bool retval = false;
      if((thiscent == CELL) || (thiscent == BOUNDARY))
      {        
        retval = a_srcbox.contains(a_ind.m_pt);
      }
      else if(thiscent == XFACE)
      {
        Point hipt = a_ind.m_pt;
        Point lopt = a_ind.m_pt - Point::Basis(0, 1);
        retval = (a_srcbox.contains(hipt) || a_srcbox.contains(lopt));
      }
      else if(thiscent ==  YFACE)
      {
        Point hipt = a_ind.m_pt;
        Point lopt = a_ind.m_pt - Point::Basis(1, 1);
        retval = (a_srcbox.contains(lopt) || a_srcbox.contains(hipt));
      }
#if DIM==3      
      else if(thiscent ==  ZFACE)
      {
        Point hipt = a_ind.m_pt;
        Point lopt = a_ind.m_pt - Point::Basis(2, 1);
        retval = (a_srcbox.contains(lopt) || a_srcbox.contains(hipt));
      }
#endif
      else
      {
        PROTO_ASSERT(false, "unaccounted-for centering");
      }
      return retval;
    }
    ///
    vector<index_t> getIndiciesInBox(const Box & a_bx) const
    {
      vector<index_t> indicies;

 
      for (unsigned int i=0; i<m_indicies->size(); ++i)
      {
        const index_t & ind = (*m_indicies)[i];
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
                         unsigned int a_numcomps) const
    {
      PR_TIME("HostIrregData::charsize");
      using ProtoSer::linearCharBufSize;
      size_t retval = linearCharBufSize(m_inputBox);
      retval += sizeof(int); //size of the vector;
      const auto& data     = *this;
      const auto& indicies = *m_indicies;
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
      return retval;
    }

    /// recreates this object into a stream of bytes
    inline
    void energize(void*  a_buf) const
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
      int vecsize = m_indicies->size();
      *intbuf = vecsize;
      charsiz += sizeof(int); 
      charbuf += sizeof(int);
      auto& indicies = *m_indicies;

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
    HostIrregData(void *a_buf)
    {
      PR_TIME("HostIrregData::char_constructor");
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
    }
    
    ///
    inline
    void linearOut(void*               a_buf,
                   const Box       &   a_bx,
                   unsigned int        a_startcomp,
                   unsigned int        a_numcomps) const
    {
      PR_TIME("HostIrregData::linearOut");
      energize(a_buf);
    }

    inline
    void  linearIn(void*                  a_buf,
                   const Box       &      a_bx,
                   unsigned int           a_startcomp,
                   unsigned int           a_numcomps)
    {
      PR_TIME("HostIrregData::linearIn");
      HostIrregData<cent, data_t, ncomp> src(a_buf);
      copy(src, a_bx, a_startcomp, a_bx, a_startcomp, a_numcomps);
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
              unsigned int a_inco)
    {
      PR_TIME("HostIrregData::copy");
      Box interBox = a_bxFrom & a_bxTo;

      for(int iind = 0; iind < m_indicies->size(); iind++)
      {
        auto& vof = (*m_indicies)[iind];
        if((interBox.contains(vof.m_pt)) && (a_src.copyContains(vof, interBox)))
        {
          for(unsigned int ico = 0; ico < a_inco; ico++)
          {
            unsigned int isrc = a_isrc + ico;
            unsigned int idst = a_idst + ico;
            if (this->contains(vof) && a_src.contains(vof))
              {
                (*this)(vof, idst) = a_src(vof, isrc);
              }
          }
        }
      }
    }


    //just used for debugging
    inline void define(const Box& a_box,
                       const shared_ptr< vector< index_t > >& a_indicies)
    {
      this->m_inputBox   = a_box;
      define(a_indicies);
    }
  protected:
    ///
    inline void define( const shared_ptr< vector< index_t > >& a_indicies)
    {
      m_indicies = a_indicies;
      std::size_t count = this->size();
      m_data=::std::shared_ptr<data_t>(new data_t[count], [](data_t* p) { delete[] p;});
      m_rawPtr = m_data.get();
    }



    Box m_inputBox; //exists if defined as part of a boxdata
    shared_ptr<vector<  index_t > >        m_indicies;
    shared_ptr<data_t>         m_data;
    data_t*                    m_rawPtr;
    HostIrregData(const HostIrregData<cent, data_t,ncomp>& a_src);
    HostIrregData& operator=(const HostIrregData  <cent, data_t,ncomp>& a_src);


  };

}

#endif
