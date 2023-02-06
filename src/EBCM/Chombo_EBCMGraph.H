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
  
  ///
  template<CENTERING cent>
  class BlobFaceIndex
  {
  public:
    //centroid location
    RealVect m_cent;
    EBIndex<CELL>  m_lo;
    EBIndex<CELL>  m_hi;
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
    
    vector<EBIndex<CELL> > m_subblobs;
    void defineAsRegular()
    {
      MayDay::Error("not implemented");
    }
  };


  ///
  /**
     unpacks EB description into a fully unstructured set of volumes 
     Handy set calculus is provided to assist in keeping the final data structure sparse.
     return value is me anticipating the need for an error code
  */
  inline  int
  unpackEBDescription(std::set<EBIndex<CELL> >              & a_soluVolumes,
                      const EBGraph                         & a_graph)

  {
    a_regIrregCov.clear();
    a_soluVolumes.clear();
    shared_ptr<std::set<EBIndex<CELL> > > retval(new std::set<EBIndex<CELL> >() );
    auto valid= a_graph.validBox();
    for(auto bit = valid.begin(); bit != valid.end(); ++bit)
    {
      if(a_graph->isCovered(*bit))
      {
        //covered bits are not part of the solution so do not go into the set
      }
      else if( a_graph->isRegular())
      {
          
        EBIndex<CELL> regvol;
        regvol.defineAsRegular(*bit);
        a_soluVolumes.insert(regvol);
      }
      else
      {
        auto vofs = a_graph->getVoFs(*bit);
        for(int ivof = 0; ivof < vofs.size(); ivof++)
        {
          a_soluVolumes.insert(vofs[ivof]);
        }
      }
      
      return 0;
    }
  }
  ///
  /**
     Takes fully unstructured EB data and makes it sparse again by using the 
     Blob structures.   
     The set that results is a fully unstructured (and merged, if you ask for it) graph.
     a_ebVolumesLeft this gets emptied as the merger goes on.    This side effect 
     on a temporary data structure allows me to not copy it.   
  */
  inline  int
  makeTheBlobs(std::set<BlobIndex     >              & a_soluVolumes,
               std::set<EBIndex<CELL> >              & a_ebVolumesLeft,//should be null at exit
               const EBGraph                         & a_graph,
               const EBData                          & a_data,
               const double                          & a_dx,
               bool                                    a_mergeSmallCells)
  {

    
    if(!a_ebVolumesLeft.isEmpty())
    {
      int isize = a_ebVolumesLeft.size();
      pout() << "EBCMGraph.H makeTheBlobs: "
             << isize << " eb volumes are still left over" << endl;
      MayDay::Error();
    }
    return 0;
  }

  ///Translate irregular graph into the sparse  representation that we actually store.
  /**
     Basically, we just do not add regular volumes and set the integer in the map to the correct thing.
     a_volumes: irregular volumes   we do not store unmerged regular cell volumes              
     a_map:     map to box stuff regular, -2 covered, >= 0 is the index into m_volumes.  merged cells both hold the same index.
     a_soluVolumes: full irregular representation
     
  **/
  inline int
  sparsifyTheBlobs(vector<BlobVoluIndex>                 & a_volumes,
                   HostBoxData<int  , 1>                 & a_map,
                   std::set<BlobIndex  >                 & a_soluVolumes)
  {
  }


  ///Class to make merged sparse embedded boundary description
  /**
   */
  class BlobGraph
  {
  public:


    BlobGraph(shared_ptr<EBGraph> a_graph,
              shared_ptr<EBData>  a_data,
              Real                a_dx,
              bool                a_mergeSmallCells) const
    {

      std::set< EBIndex<CELL> >  allEBInd;
      std::set< BlobVoluIndex >  allBlobs;
      
      upackEBDescription(allEBInd, a_graph);
      makeTheBlobs(allBlobs,  allEBInd, a_graph, a_data, a_dx, a_mergeSmallCells);

      //translate irregulkar graph into sparse  representation that we actually store
      sparsifyTheBlobs(m_volumes, m_map, allBlobs)
      return retval;
                         
    }
  private:

    //irregular volumes   we do not store unmerged regular cell volumes              
    vector<BlobVoluIndex> m_volumes;
    //-1 regular, -2 covered, >= 0 is the index into m_volumes.  merged cells both hold the same index.
    HostBoxData<int, 1>   m_map;

      
    //weak construction is weak
    BlobGraph();
  };
}


#endif
