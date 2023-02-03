#ifndef __PROTO_MomentHolder__
#define __PROTO_MomentHolder__

#include "implem/Proto_IrregGraph.H"
#include "EBProto.H"

// This consists of a
namespace Proto
{
  // Used in a hack to allow indexing dimensional components. This save a lot of code duplication higher up
  template<CENTERING>
  struct typenameTag {};

  // Structs to make sending data around sane.
  // These exist to hide deep design choices that are uncomfortable to work with

  template <int order>
  struct MomentHolder
  {
    // Yes are a lot of aliases here
    // auto is your friend if you just need to access them, which is usually all that is required
    typedef IndexedMoments<DIM  , order> IndMomDIM;
    typedef IndexedMoments<DIM-1, order> IndMomSDMinOne;
    typedef MomentIterator<DIM  , order> MomItDIM;
    typedef MomentIterator<DIM-1, order> MomItSDMinOne;
    typedef HostIrregData<CELL    ,  IndMomDIM , 1>  VoluData;
    typedef HostIrregData<BOUNDARY,  IndMomDIM , 1>  EBFaData;
    typedef HostIrregData<XFACE, IndMomSDMinOne, 1>  XFacData;
    typedef HostIrregData<YFACE, IndMomSDMinOne, 1>  YFacData;
    typedef HostIrregData<ZFACE, IndMomSDMinOne, 1>  ZFacData;
    typedef HostIrregData<CELL    ,  IndMomDIM , 1>  EBNormalData;

    // Don't call me, look in GeometryService
    MomentHolder(const VoluData     & a_VoluData,
                 const EBFaData     & a_EBFaceData,
                 const XFacData     & a_XFacData,
                 const YFacData     & a_YFacData,
                 const ZFacData     & a_ZFacData,
                 const EBNormalData & a_XEBNormData,
                 const EBNormalData & a_YEBNormData,
                 const EBNormalData & a_ZEBNormData)
      {
        // I promise not to change things. Is casting away consts good? No. But this is a bandaid for bigger problems
        m_VoluData = const_cast<VoluData*>(&a_VoluData);
        m_EBFaceData = const_cast<EBFaData*>(&a_EBFaceData);
#if DIM >= 1
        m_XFaceData = &a_XFacData;
        m_XEBNormData = &a_XEBNormData;
#endif
#if DIM >= 2
        m_YFaceData = &a_YFacData;
        m_YEBNormData = &a_YEBNormData;
#endif
#if DIM >= 3
        m_ZFaceData = &a_ZFacData;
        m_ZEBNormData = &a_ZEBNormData;
#endif
      }

    inline const VoluData&
    volumeMoments() const
      {
        return *m_VoluData;
      }

    inline const EBFaData&
    ebFaceMoments() const
      {
        return *m_EBFaceData;
      }

  private://dispatched access to hide template types
    const HostIrregData<XFACE, IndMomSDMinOne, 1>&
    faceData_dispatch(typenameTag<XFACE>) const;

    const HostIrregData<YFACE, IndMomSDMinOne, 1>&
    faceData_dispatch(typenameTag<YFACE>) const;

    const HostIrregData<ZFACE, IndMomSDMinOne, 1>&
    faceData_dispatch(typenameTag<ZFACE>) const;

  public:
    /// ONE ACCESSOR TO RULE THEM ALL
    /// call this with odd syntax "a_momentData.template faceMoments<XFACE>();"
    template<CENTERING cent>
    const HostIrregData<cent, IndMomSDMinOne, 1>&
    faceMoments() const
      {
        return faceData_dispatch(typenameTag<cent>{});
        //MayDay::Error("invalid face centering");
      }

    inline const EBNormalData&
    ebNormal(const int a_dir) const
      {
        CH_assert((a_dir >= 0) && (a_dir < SpaceDim));
        if (a_dir == 0)
          {
            return *m_XEBNormData;
          }
        else if (a_dir == 1)
          {
            return *m_YEBNormData;
          }
        else if (a_dir == 2)
          {
            return *m_ZEBNormData;
          }
        Chombo4::MayDay::Error("EBnormal access error");
        return *m_XEBNormData; // shut up the compiler warnings
      }

  private:
    // Raw pointers to all kinds of information
    // I don't own these, I just make them accessible elsewhere
    Box m_ValidBox;
    EBGraph* m_Graph = nullptr;
    VoluData* m_VoluData = nullptr;
    EBFaData* m_EBFaceData = nullptr;
//#if DIM >= 1
    XFacData* m_XFaceData = nullptr;
    EBNormalData* m_XEBNormData = nullptr;
//#elif DIM >= 2
    YFacData* m_YFaceData = nullptr;
    EBNormalData* m_YEBNormData = nullptr;
//#elif DIM >= 3
    ZFacData* m_ZFaceData = nullptr;
    EBNormalData* m_ZEBNormData = nullptr;
//#endif

  };


  template<int order>
  const HostIrregData<XFACE, IndexedMoments<DIM-1, order>, 1>&
  MomentHolder<order>::faceData_dispatch(typenameTag<XFACE>) const
  {
    CH_assert(DIM >= 0);
    return *(this->m_XFaceData);
  }

  template<int order>
  const HostIrregData<YFACE, IndexedMoments<DIM-1, order>, 1>&
  MomentHolder<order>::faceData_dispatch(typenameTag<YFACE>) const
  {
    CH_assert(DIM >= 1);
    return *(this->m_YFaceData);
  }

  template<int order>
  const HostIrregData<ZFACE, IndexedMoments<DIM-1, order>, 1>&
  MomentHolder<order>::faceData_dispatch(typenameTag<ZFACE>) const
  {
    CH_assert(DIM >= 2);
    return *(this->m_ZFaceData);
  }

}

#endif
