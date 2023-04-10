#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _EBCMHostLevelData_H____
#define _EBCMHostLevelData_H____

#include "Proto.H"
#include "EBProto.H"
#include "Chombo_ProtoInterface.H"
#include "Chombo_EBChombo.H"
#include "Chombo_ProtoFactories.H"

#include "Chombo_SPMD.H"
#include "Chombo_EBCM_VolBoxData.H"
#include "Chombo_EBDataChoreography.H"
#include "Chombo_GeometryService.H"
#include "Chombo_AMRIO.H"
#include <memory>
namespace EBCM
{
  /// Simple structure to hold data on a level in an ebcm context.  volume data only
  /**
   */
  template<typename data_t, int  NUMCOMPS, int order>
  class HostLevelData
  {
  public:
    typedef EBCM::VolBoxData<data_t, NUMCOMPS, order>            ebcm_fab_t;
    typedef Chombo4::NullConstructorDataFactory<    ebcm_fab_t>  null_factory_t;
    typedef CH4_Data_Choreography::DistributedData< ebcm_fab_t > mpi_frame_t;
    typedef    EBCM::MetaDataLevel<        order >               ebcm_meta_t;
    typedef    EBCM::SubVolumeIterator<    order >               ebcm_subvol_it;
    typedef Chombo4::GeometryService<      order >               ch_geoserv;
    typedef Chombo4::DisjointBoxLayout                           ch_dbl;
    typedef Chombo4::ProblemDomain                               ch_probdom;
    typedef Chombo4::DataIterator                                ch_dit;
    typedef Chombo4::IntVect                                     ch_iv;
    typedef Chombo4::FArrayBox                                   ch_fab;
    typedef Chombo4::BoxIterator                                 ch_bit;
    typedef Chombo4::LevelData<ch_fab>                           ch_ldfab;
    
    HostLevelData(const shared_ptr<ebcm_meta_t>& a_meta,
                  const ch_iv                  & a_ghost)
    {
      m_meta = a_meta;

      const auto& graphs = (m_meta->m_graphs);
      const auto& grids  = (m_meta->m_grids);
      m_data = shared_ptr<mpi_frame_t>(new mpi_frame_t(grids, a_ghost, null_factory_t()));
      auto dit = grids.dataIterator();
      for(int ibox =0; ibox < dit.size(); ibox++)
      {
        auto      &   fab = (*m_data)[dit[ibox]];
        const auto& graph = (*graphs)[dit[ibox]];
        fab.define(graph);
      }
    }

    void    writeToHDF5(string a_filename) const;
    {
      ch_ldfab fabdata(m_meta->m_grids, ncomp, ch_iv::Zero);
      ch_dit dit =     m_meta->m_grids.dataIterator();
      for(int ibox = 0; ibox < dit.size(); ibox++)
      {
        auto& fabdata = fabdata[dit[ibox]];
        const auto& ebcmdata = (*m_data)(dit[ibox]);
        const auto& grid = (m_meta->m_grids)[dit[ibox]];
        for(ch_bit bit(grid); bit.ok(); ++bit)
        {
          auto  iv = bit();
          ch_pt pt = ProtoCh::getPoint(iv);
          for(int icomp = 0; icomp < ncomp; icomp++)
          {
            fabdata(iv, icomp) = ebcmdata(pt, icomp);
          }
        }
        writeLevelname(&fabdata, a_filename.c_str());
      }
    }
    shared_ptr< ebcm_meta_t >  m_meta;
    shared_ptr< mpi_frame_t >  m_data;
    
  private:    //only disallowing weak construction 
    HostLevelData();
  }; //end class hostleveldata
  
} // end namespace ebcm


#endif
