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

  class boxmap_t
  {
  public:
  
    //going to be different than EBGraph info because of coarsening
    PointSet            m_regular;
    PointSet            m_irregular;
    PointSet            m_covered;
    boxmap_t()
    {
      m_regular.clear();
      m_irregular.clear();
      m_covered.clear();
    }
  }
  class BlobGraph
  {
  public:

    shared_ptr<std::set<EBIndex<CELL> > >
    fillVolumesFromEB(shared_ptr<EBGraph> a_graph,
                      shared_ptr<EBData>  a_data,
                      Real                a_dx)

    {
      shared_ptr<std::set<EBIndex<CELL> > > retval(new std::set<EBIndex<CELL> >() );
      auto valid= a_graph.validBox();
      for(auto bit = valid.begin(); bit != valid.end(); ++bit)
      {
        if(a_graph->isCovered(*bit))
        {
          //covered bits are not part of the solution 
        }
        else if( a_graph->isRegular())
        {
          
        }
        auto vofs = a_graph->getVoFs(*bit);
        for(int ivof = 0; ivof < vofs.size(); ivof++)
        {
          const auto& vof = vofs[ivof];
        }
      }

      return retval;
    }

    BlobGraph(shared_ptr<EBGraph> a_graph,
              shared_ptr<EBData>  a_data,
              Real                a_dx,
              bool                a_mergeSmallCells) const
    {

      shared_ptr<std::set<EBIndex<CELL> > > allEBVols = fillVolumesFromEB(a_graph, a_data, a_dx);

      if(a_mergeSmallCells)
      {
        m_volumes = mergeSmallVolumes(allEBVols);
      }
      else
      {
        m_volumes = allEBVols;
      }
      return retval;
                         
    }
  private:
              
   vector<BlobVoluIndex> m_volumes;
              
    //this reflects the information in the EBIS.   in this unstructured format, none of this stuff matters.
    boxmap_t   m_boxmap;

      
    //weak construction is weak
    BlobGraph();
  };
}


#endif
