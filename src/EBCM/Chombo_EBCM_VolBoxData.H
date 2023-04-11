
#ifndef _Chombo_EBCM_VolBoxData_H_
#define _Chombo_EBCM_VolBoxData_H_

#include <cmath>
#include <cstdlib>
#include <memory>
#include "Proto.H"
#include "EBProto.H"
#include "Chombo_EBCM_Graph.H"

using std::vector;
using std::shared_ptr;

namespace EBCM
{
  ///Holder for host EBCM data over a box.   volume only.
  template <typename  data_t,  int ncomp, int order>
  class VolBoxData
  {
  public:
    typedef EBCM::EBCM_Graph< order >         ebcm_graph_t;
    typedef Proto::HostBoxData<data_t, ncomp> pr_hostdat_t;
    typedef Chombo4::MayDay                   ch_mayday;
    ///
    VolBoxData()
    {
    }


    int nComp() const
    {
      return ncomp;
    }

    void define(const ebcm_graph_t& a_graph)
    {
      m_graph = a_graph;
      pr_box prgid = ProtoCh::getProtoBox(a_graph.m_grid);
      m_regdat.define(prgid);
      if(a_graph.hasIrregular())
      {
        int numIrreg = a_graph.m_volumes->size();
        int nvar = numIrreg*ncomp;
        m_irrdat.resize(nvar);
      }
    }

    

    data_t* dataPtr(const pr_pt& a_pt, unsigned int a_icomp)
    {
      if(a_icomp >= ncomp)
      {
        ch_mayday::Error("VolBoxData access icomp too big");
      }
      data_t* retvalptr = NULL;
      
      bool useRegular = true;
      int mapval = -4586;
      if(m_graph.hasIrregular())
      {
        mapval = (*m_graph.m_map)(a_pt, 0);
        useRegular = (mapval < 0);
      }
      if(useRegular)
      {
        retvalptr = &(m_regdat(a_pt, a_icomp));
      }
      else
      {
        int numVols = m_graph.m_volumes->size();
        int vecvalcomp0 = mapval;
        int idatavec  = vecvalcomp0 + a_icomp*numVols;
        retvalptr = &(m_irrdat[idatavec]);
      }
      return retvalptr;
    }
    ///
    /**
       I will have to do something smarter than this for any performance  but for now, 
       since I am working alone on this and there is nobody here to mock me, I will try 
       to get this paper done with just point accessor functions.    
       It should not be a big deal to make this work with AggStencil.
    */
    data_t&       operator() (const pr_pt& a_pt, unsigned int a_icomp)
    {
      data_t* retvalptr = dataPtr(a_pt, a_icomp);
      return *retvalptr;
    }
    const data_t& operator() (const pr_pt& a_pt, unsigned int a_icomp) const
    {
      
      VolBoxData<data_t, ncomp, order>* castme = (VolBoxData<data_t, ncomp, order>*) (this);
      const data_t* retvalptr = castme->dataPtr(a_pt, a_icomp);
      return *retvalptr;
    }

    void setVal(data_t a_value)
    {
      m_regdat.setVal(a_value);
      for(int idat = 0; idat < m_irrdat.size(); idat++)
      {
        m_irrdat[idat] = a_value;
      }
    }
  private:
    //private because their use is not obvious.    see fancy operator() functions
    ebcm_graph_t   m_graph;
    pr_hostdat_t   m_regdat;
    vector<data_t> m_irrdat;
  };


}
#endif
