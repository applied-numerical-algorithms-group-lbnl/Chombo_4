#ifndef _AGGSTENCIL_H_
#define _AGGSTENCIL_H_

#include "EBProto.H"

namespace Proto
{
//  using std::pair;
  template<typename data_t>
  struct pair_t
  {
    EBDataLoc first;
    data_t    second;
  };

  using std::vector;
  ///replacement for Chombo's vofstencil and facestencil
  template<CENTERING cent, typename data_t>
  class LocalStencil
  {
  public:
    struct entry_t
    {
      EBIndex<cent> m_vof;
      data_t        m_weight;
      entry_t(const EBIndex<cent>& a_vof,
              const data_t       & a_weight)
      {
        m_vof = a_vof;
        m_weight = a_weight;
      }

      entry_t()
      {
      }

    };
    vector<entry_t>  m_entries;
    //diagonal weight
    data_t           m_betaWeight;    
    //whether diagonal weight has been set
    //most stencils do not need it
    bool             m_betaWeightSet; 

    void
    setBetaWeight(data_t a_value)
    {
      m_betaWeight    = a_value;
      m_betaWeightSet = true;
    }
    ///
    LocalStencil()
    {
      m_betaWeightSet = false;
    }
    
    ///
    ~LocalStencil()
    {}
    
    ///
    inline void clear()
    {
      m_entries.resize(0);
    }
    
    ///
    int size() const
    {
      return m_entries.size();
    }
    
    ///this api should be stable
    inline void add(const entry_t& a_entry)
    {
      add(a_entry.m_vof, a_entry.m_weight);
    }

    ///this might change (if a variable number argument is added, for example)
    inline void add(const EBIndex<cent> a_vof, const data_t& a_weight)
    {
      bool alreadyhere = false;
      for (int ivof = 0; ivof < m_entries.size(); ivof++)
      {
        if(m_entries[ivof].m_vof == a_vof)
        {
          alreadyhere = true;
          m_entries[ivof].m_weight += a_weight;
        }
      }
      if (!alreadyhere)
      {
        m_entries.push_back(entry_t(a_vof, a_weight));
      }
    }
   
    ///add full stencil
    LocalStencil<cent, data_t>& operator+=(const LocalStencil<cent, data_t>& a_input)
    {
      if(&a_input != this)
      {
        for(int ient = 0; ient < a_input.size(); ient++)
        {
          add(a_input.m_entries[ient]);
        }
      }
      return *this;
    }

    //return true if vof found
    bool getValueAtVoF(data_t& a_value, const EBIndex<cent> a_vof)
    {
      bool retval = false;
      for(int ient = 0; ient < m_entries.size(); ient++)
      {
        if(m_entries[ient].m_vof == a_vof)
        {
          retval = true;
          a_value = m_entries[ient].m_weight;
        }
      }
      return retval;
    }
    ///scaling operator
    void operator*=  (const data_t & a_scaling)
    {
      for(int ient = 0; ient < m_entries.size(); ient++)
      {
        m_entries[ient].m_weight *= a_scaling;
      }
    }
    
    void print() const
    {
      using std::cout;
      using std::endl;
      using std::setprecision;
      using std::ios;
      cout << "local stencil contains" << endl;
      cout    << std::setprecision(8)
              << std::setiosflags(ios::showpoint)
              << std::setiosflags(ios::scientific);
      for(int ivec =  0; ivec < m_entries.size(); ivec++)
      {
        cout << "(" << m_entries[ivec].m_vof.m_pt << "," << m_entries[ivec].m_weight << endl;
      }
    }
  };


  /// 
  /**
  */
  template <CENTERING srcCenter, CENTERING dstCenter, typename data_t>
  class AggStencil
  {
  public:
    ///
    /**
    */
    AggStencil(const vector< EBIndex<     dstCenter        > > & a_dstVoFs,
               const vector< LocalStencil<srcCenter, data_t> > & a_stencil,
               const EBGraph                                   & a_srcGraph,
               const EBGraph                                   & a_dstGraph,
               const Box                                       & a_srcValid,
               const Box                                       & a_dstValid,
               const Point                                     & a_srcGhost,
               const Point                                     & a_dstGhost);


    ///
    ~AggStencil()
    {
#ifdef PROTO_CUDA   
      //just contains a bunch of cuda frees
      cudaRelease();
#else
      hostRelease();
#endif 
    }

     ///
    /**
       lphi = L(phi)
     */
    void apply(EBBoxData<      dstCenter, data_t, 1>     & a_lph,
               const EBBoxData<srcCenter, data_t, 1>     & a_phi,
               const data_t                              & a_scale,
               bool                                        a_incrementOnly) const;

#ifdef COUNTER
    unsigned int getCounterVofs();
    unsigned int getCounterRegVof();
    unsigned int getCounterZeroWeight();
    unsigned int getCounterTotSten();
#endif

  protected:
    //  /this stuff lives on the host
    typedef vector< pair_t<data_t>  > stencil_t;
    vector<stencil_t>    m_ebstencil;
    vector<EBDataLoc>    m_dstAccess;

//    //begin debugging hooks
//    vector< EBIndex<     dstCenter        > > m_dstVoFs;
//    vector< LocalStencil<srcCenter, data_t> > m_stencil;
//    //end debugging hooks
    void hostApply(EBBoxData<      dstCenter, data_t, 1>     & a_lph,
                   const EBBoxData<srcCenter, data_t, 1>     & a_phi,
                   const data_t                              & a_scale,
                   bool                                      a_incrementOnly) const;

    void hostPrep(const vector<pair_t<data_t> >   &   h_serializedStenInfo,
                  const vector<uint64_t>          &   h_stenSizes,
                  const vector<uint64_t>          &   h_stenStart);

    void hostRelease();

#ifdef COUNTER
    unsigned int m_counterTotSten = 0;
    unsigned int m_counterRegVof = 0;
    unsigned int m_counterZeroWeight = 0;
#endif

#ifdef PROTO_CUDA
    void cudaPrep(const vector<pair_t<data_t> >   &   h_serializedStenInfo,
                  const vector<uint64_t> &   h_stenSizes,
                  const vector<uint64_t> &   h_stenStart);


    void cudaRelease();
    void cudaApply(EBBoxData<      dstCenter, data_t, 1>     & a_lph,
                   const EBBoxData<srcCenter, data_t, 1>     & a_phi,
                   const data_t                              & a_scale,
                   bool                                a_incrementOnly) const;


    int m_stride = 0;
    int m_blocks = 0;
#endif
    //this is for flop counts
    int  m_flopsThisApply = 0;        

#ifdef NEWSTRUCT
    /// this stuff lives on the pinned memory
    uint16_t*         d_stenLen = nullptr;
    bool*             d_destMeshType = nullptr;
    uint64_t*         d_destOffset = nullptr;
    bool*             d_srcMeshType = nullptr;
    uint64_t*         d_srcOffset = nullptr;
    data_t*           d_srcWeight = nullptr;

// for compilation
    pair_t<data_t>*      d_ebstencil = nullptr;
    EBDataLoc*           d_dstaccess = nullptr;
    uint64_t*            d_stensizes = nullptr;
    uint64_t*            d_stenstart = nullptr;
#else

    /// this stuff lives on the device
    pair_t<data_t>*      d_ebstencil = nullptr;
    EBDataLoc*           d_dstaccess = nullptr;
    uint64_t*            d_stensizes = nullptr;
    uint64_t*            d_stenstart = nullptr;
#endif

    Box m_srcBox, m_dstBox;
  private:
    /// disallowed operators.  
    AggStencil();
    AggStencil& operator=(const AggStencil& stenin);
    AggStencil(const AggStencil& stenin);

  };
}
#include "Proto_AggStencilImplem.H"

#endif
