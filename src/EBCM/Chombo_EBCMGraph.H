#ifndef __PROTO_EBCMGraph__
#define __PROTO_EBCMGraph__

#include "implem/Proto_IrregGraph.H"
#include "EBProto.H"
#include "Proto_RealVect.H"
using std::vector;
using std::set;

namespace EBCM
{
  class BlobVolIndex;
  
  ///BlobFaceIndex are made by collecting up EBIndex<cent> faces along a merged cell.  can be an EB face
  template<CENTERING cent>
  class BlobFaceIndex
  {
  public:
    //centroid location
    RealVect m_cent;
    BlobVolIndex m_lo;
    BlobVolIndex m_hi;
    IndexedMoments<DIM> m_facmom;
    defineAsRegular();

  };

  
  class BlobVoluIndex
  {
  public:
    vector<BlobFaceIndex<XFACE> > m_xfaces;
    vector<BlobFaceIndex<YFACE> > m_yfaces;
    vector<BlobFaceIndex<ZFACE> > m_zfaces;
    BlobFaceIndex<      BOUNDARY> m_ebface;
    
    bool                          m_regular;
    //centroid location
    RealVect m_cent;
    IndexedMoments<DIM> m_volmom;
    IndexedMoments<DIM> m_ebfmom;
    Point    m_pt; //just used as a unique identifier
    defineAsRegular();
  };

  class BlobGraph
  {
  public:
    //tells which cells eb thinks are regular, irregular covered
    shared_ptr<EBGraph> m_graph;

    //going to be different than EBGraph info because of coarsening
    PointSet            m_regular;
    PointSet            m_irregular;
    PointSet            m_covered;
    
    
    bool isRegular(Point a_pt) const
    {
      return m_regular.contains(a_pt);
    }
    
    bool isIrregular(Point a_pt) const
    {
      return m_irregular.contains(a_pt);
    }

    bool isCovered(Point a_pt) const
    {
      return m_covered.contains(a_pt);
    }
    
  };
}


#endif
