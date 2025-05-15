#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _EBCMFluxLevelData_H____
#define _EBCMFluxLevelData_H____

#include "Proto.H"
#include "EBProto.H"
#include "Chombo_ProtoInterface.H"
#include "Chombo_EBChombo.H"
#include "Chombo_ProtoFactories.H"

#include "Chombo_SPMD.H"
#include "Proto_HostBoxData.H"
#include "Chombo_EBDataChoreography.H"
#include "Chombo_GeometryService.H"
#include "Chombo_AMRIO.H"
#include <memory>
namespace EBCM
{
  /// 
  /**
     Simple structure to hold fluxes on a level in an ebcm context.  
   */
  template<typename data_t, int  ncomp, int order>
  class FluxLevelData
  {
  public:
    typedef Proto::HostBoxData<data_t, ncomp>                    pr_fab_t;
    typedef Chombo4::NullConstructorDataFactory<    pr_fab_t>    null_factory_t;
    typedef CH4_Data_Choreography::DistributedData< pr_fab_t >   mpi_frame_t;
    typedef    EBCM::MetaDataLevel<        order >               ebcm_meta_t;
    typedef Chombo4::GeometryService<      order >               ch_geoserv;
    typedef Chombo4::DisjointBoxLayout                           ch_dbl;
    typedef Chombo4::ProblemDomain                               ch_probdom;
    typedef Chombo4::DataIterator                                ch_dit;
    typedef Chombo4::DataIndex                                   ch_index;
    typedef Chombo4::IntVect                                     ch_iv;
    typedef Chombo4::FArrayBox                                   ch_fab;
    typedef Chombo4::BoxIterator                                 ch_bit;
    typedef Chombo4::LevelData<ch_fab>                           ch_ldfab;
    typedef Chombo4::MayDay                                      ch_mayday;


    
    FluxLevelData(const shared_ptr<ebcm_meta_t>& a_meta,
                  const int                    & a_ghost)
    {
      m_meta = a_meta;
      //const auto& graphs = (m_meta->m_graphs);
      const auto& grids  = (m_meta->m_grids);
      for(int idir  = 0; idir < DIM; idir++)
      {
        m_data[idir] = shared_ptr<mpi_frame_t>(new mpi_frame_t(grids, a_ghost*ch_iv::Unit, null_factory_t()));
        ch_box face_dom = m_meta->m_domain;
        face_dom.grow(idir, 1);
        auto dit = grids.dataIterator();
        for(int ibox =0; ibox < dit.size(); ibox++)
        {
          auto valid = grids[dit[ibox]];
          auto grown = valid.grow(a_ghost-1);
          grown.grow(idir, 1);
          grown &= face_dom;
          (*(m_data[idir]))[dit[ibox]].define(grown);
        }
      }
    }

    void    writeToHDF5(string a_prefix) const
    {
#ifdef CH_USE_HDF5
      for(int idir = 0; idir < DIM; idir++)
      {
        string filename = a_prefix + to_string(idir) + string(".hdf5");
        const auto& data_ptr =  m_data[idir];
        ch_ldfab levelfabs(m_meta->m_grids, ncomp, ch_iv::Zero);
        ch_dit dit =     m_meta->m_grids.dataIterator();
        const auto& mpibit = (*data_ptr);
        for(int ibox = 0; ibox < dit.size(); ibox++)
        {
          auto& fabdata            = levelfabs[dit[ibox]];
          const auto& ebcmdata        = mpibit[dit[ibox]];
          const auto& grid = (m_meta->m_grids)[dit[ibox]];
          for(ch_bit bit(grid); bit.ok(); ++bit)
          {
            auto  iv = bit();
            pr_pt pt = ProtoCh::getPoint(iv);
            for(int icomp = 0; icomp < ncomp; icomp++)
            {
              fabdata(iv, icomp) = ebcmdata(pt, icomp);
            } //end loop over components
          }   //end loop over cells in box
        }     //end loop over boxes
        Chombo4::pout() <<
          "fluxdata: writeToHDF5 is writing to "
                        << filename
                        << std::endl;
        writeLevelname(&levelfabs, filename.c_str());
      } //end loop over directions
#else
      Chombo4::pout() <<
        "Confusion reigns. \n The function writeToHDF5 was called but CH_USE_HDF5 is false.  \n This message is to let you know that happened."
                      << std::endl;
#endif      
    }
    
    shared_ptr< ebcm_meta_t >  m_meta;
    shared_ptr< mpi_frame_t >  m_data[DIM];

    void exchange(bool a_printStuff)
    {
      for(int idir = 0; idir < DIM; idir++)
      {
        m_data[idir]->exchange(a_printStuff);
      }
    }
  private:    //only disallowing weak construction 
    FluxLevelData();
  }; //end class fluxleveldata
  
} // end namespace ebcm


#endif
