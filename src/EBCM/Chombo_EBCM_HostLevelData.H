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
#include "Proto_HostBoxData.H"
#include "Chombo_EBDataChoreography.H"
#include "Chombo_GeometryService.H"
#include "Chombo_AMRIO.H"
#include <memory>
namespace EBCM
{
  /// Simple structure to hold data on a level in an ebcm context.  volume data only
  /**
   */
  template<typename data_t, int  ncomp, int order>
  class HostLevelData
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
    typedef Chombo4::IntVect                                     ch_iv;
    typedef Chombo4::FArrayBox                                   ch_fab;
    typedef Chombo4::BoxIterator                                 ch_bit;
    typedef Chombo4::LevelData<ch_fab>                           ch_ldfab;
    typedef Chombo4::MayDay                                      ch_mayday;

    void operator+=(const HostLevelData<data_t, ncomp, order>& a_input)
    {
      ch_dit dit =     m_meta->m_grids.dataIterator();
      const auto& mpibit = (*m_data);
      for(int ibox = 0; ibox < dit.size(); ibox++)
      {
        const auto& inpdata         = a_input[dit[ibox]];
        auto&       outdata         = (*this)[dit[ibox]];
        const auto& grid  = (m_meta->m_grids)[dit[ibox]];
        for(ch_bit bit(grid); bit.ok(); ++bit)
        {
          auto  iv = bit();
          pr_pt pt = ProtoCh::getPoint(iv);
          for(int icomp = 0; icomp < ncomp; icomp++)
          {
            outdata(pt, icomp) += indata(pt, icomp);
          } //end loop over components
        }   //end loop over cells in box
      }     //end loop over boxes
    }

    void operator-=(const HostLevelData<data_t, ncomp, order>& a_input)
    {
      ch_dit dit =     m_meta->m_grids.dataIterator();
      const auto& mpibit = (*m_data);
      for(int ibox = 0; ibox < dit.size(); ibox++)
      {
        const auto& inpdata         = a_input[dit[ibox]];
        auto&       outdata         = (*this)[dit[ibox]];
        const auto& grid  = (m_meta->m_grids)[dit[ibox]];
        for(ch_bit bit(grid); bit.ok(); ++bit)
        {
          auto  iv = bit();
          pr_pt pt = ProtoCh::getPoint(iv);
          for(int icomp = 0; icomp < ncomp; icomp++)
          {
            outdata(pt, icomp) -= indata(pt, icomp);
          } //end loop over components
        }   //end loop over cells in box
      }     //end loop over boxes

    }

    void setVal(data_t a_value)
    {
      ch_dit dit =     m_meta->m_grids.dataIterator();
      const auto& mpibit = (*m_data);
      for(int ibox = 0; ibox < dit.size(); ibox++)
      {
        const auto& data  = (*this)[dit[ibox]];
        const auto& grid  = (m_meta->m_grids)[dit[ibox]];
        for(ch_bit bit(grid); bit.ok(); ++bit)
        {
          auto  iv = bit();
          pr_pt pt = ProtoCh::getPoint(iv);
          for(int icomp = 0; icomp < ncomp; icomp++)
          {
            data(pt, icomp) = a_value;
          } //end loop over components
        }   //end loop over cells in box
      }     //end loop over boxes
      
    }
    
    double maxNorm() const
    {
      double maxval = -4586;
      ch_dit dit =     m_meta->m_grids.dataIterator();
      const auto& mpibit = (*m_data);
      for(int ibox = 0; ibox < dit.size(); ibox++)
      {
        const auto& data  = (*this)[dit[ibox]];
        const auto& grid  = (m_meta->m_grids)[dit[ibox]];
        for(ch_bit bit(grid); bit.ok(); ++bit)
        {
          auto  iv = bit();
          pr_pt pt = ProtoCh::getPoint(iv);
          for(int icomp = 0; icomp < ncomp; icomp++)
          {
            double locvalue = data(pt, icomp);
            maxval = std::max(maxval, std::abs(locvalue));
          } //end loop over components
        }   //end loop over cells in box
      }     //end loop over boxes
      
#ifdef CH_MPI    
      double sendBuf = maxval;
      int result = MPI_Allreduce(&sendBuf, &maxval, 1, MPI_CH_REAL, MPI_MAX, Chombo_MPI::comm);

      if (result != MPI_SUCCESS)
      {
        MayDay::Error("Communication error in EBCM_HostLevelData::maxNorm");
      }
    
#endif
      return maxval;
    }
    HostLevelData(const shared_ptr<ebcm_meta_t>& a_meta,
                  const int                    & a_ghost)
    {
      m_meta = a_meta;
      const auto& graphs = (m_meta->m_graphs);
      const auto& grids  = (m_meta->m_grids);
      m_data = shared_ptr<mpi_frame_t>(new mpi_frame_t(grids, a_ghost*ch_iv::Unit, null_factory_t()));
      auto dit = grids.dataIterator();
      for(int ibox =0; ibox < dit.size(); ibox++)
      {
        auto valid = grids[dit[ibox]];
        auto grown = valid.grow(a_ghost);
        grown &= a_meta->m_domain;
        (*m_data)[dit[ibox]].define(grown);
      }
    }

    void    writeToHDF5(string a_filename) const
    {
      ch_ldfab levelfabs(m_meta->m_grids, ncomp, ch_iv::Zero);
      ch_dit dit =     m_meta->m_grids.dataIterator();
      const auto& mpibit = (*m_data);
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
      writeLevelname(&levelfabs, a_filename.c_str());
    }
    
    shared_ptr< ebcm_meta_t >  m_meta;
    shared_ptr< mpi_frame_t >  m_data;

    void exchange(bool a_printStuff)
    {
      m_data->exchange(a_printStuff);
    }
  private:    //only disallowing weak construction 
    HostLevelData();
  }; //end class hostleveldata
  
} // end namespace ebcm


#endif
