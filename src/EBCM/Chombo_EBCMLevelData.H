#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _EBCMLevelData_H____
#define _EBCMLevelData_H____

#include "Proto.H"
#include "EBProto.H"
#include "Chombo_LevelData.H"
#include "Chombo_AMRIO.H"
#include "Chombo_ProtoInterface.H"
#include "Chombo_EBChombo.H"
#include "Chombo_ProtoFactories.H"
#include "Chombo_LevelBoxData.H"
#include "Chombo_SPMD.H"
#include "Proto_EBBoxData.H"
#include "Proto_EBHostData.H"
#include "Chombo_EBDataChoreography.H"
#include "Chombo_GeometryService.H"
#include <memory>
namespace EBCM
{
  ///
  /**
   */
  template<typename data_t, int  NUMCOMPS, int order>
  class EBCMLevelData
    : public CH4_Data_Choreography::DistributedData<EBCellMergeData<data_t, NUMCOMPS> >
  {
  public:
    typedef Chombo4::GeometryService<      order >   ch_geoserv;
    typedef    EBCM::MetaDataLevel<        order >   ebcm_meta;
    typedef    EBCM::SubVolumeIterator<    order >   ebcm_subvol_it;
    typedef Chombo4::DisjointBoxLayout ch_dbl;
    typedef Chombo4::ProblemDomain ch_probdom;
    typedef Chombo4::DataIterator  ch_dit;

    EBCMLevelData(const ebcm_meta& a_meta)
    {
    }

  private:
    //
    EBCMLevelData();
  }; //end class ebcmleveldata
  
} // end namespace ebcm


#endif
