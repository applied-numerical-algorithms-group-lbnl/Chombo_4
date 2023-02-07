#ifndef __Chombo_EBCMGraph__
#define __Chombo_EBCMGraph__

#include "implem/Proto_IrregGraph.H"
#include "EBProto.H"
#include "Proto_RealVect.H"
using std::vector;
using std::set;

namespace EBCM
{
  using std::set;
  using std::shared_ptr;

  typedef Proto::EBGraph       pr_ebgrap;
  typedef Proto::EBData        pr_ebdata;
  typedef Proto::pr_vof pr_vof;
  
  ///
  template<CENTERING cent>
  class EBCMFace
  {
  public:
    //centroid location
    RealVect m_cent;
    pr_vof  m_lo;
    pr_vof  m_hi;
    IndexedMoments<DIM> m_facmom;
    defineAsRegular();

  };

  
  class EBCMVolu
  {
  public:
    vector<EBCMFace<XFACE> > m_xfaces;
    vector<EBCMFace<YFACE> > m_yfaces;
    vector<EBCMFace<ZFACE> > m_zfaces;
    EBCMFace<      BOUNDARY> m_ebface;

    
    bool                          m_regular;
    //centroid location
    RealVect m_cent;
    IndexedMoments<DIM> m_volmom;
    IndexedMoments<DIM> m_ebfmom;
    Point    m_pt; //just used as a unique identifier
    
    vector<pr_vof > m_subblobs;
    void defineAsRegular(double a_dx)
    {
      MayDay::Error("not implemented");
    }
  };


  ///create single-index ebcmvolume from EB description
  inline EBCMVolu
  createEBCMVolume(const pr_vof                          & a_vof,
                   const EBGraph                         & a_graph,
                   const EBData                          & a_data,
                   const double                          & a_dx)
  {
    EBCMVolu retval;
    if(a_graph.isCovered(a_vof.m_pt))
    {
      MayDay::Error("cannot create EBCMVolume from a covered cell");
    }
    else if(a_graph.isRegular(a_vof.m_pt))
    {
      retval.defineAsRegular(a_dx);
    }
    else
    {
      HERE
    }

    return retval;
  }
  
  ///
  /**
     unpacks EB description into a fully unstructured set of volumes 
     Handy set calculus is provided to assist in keeping the final data structure sparse.
     return value is me anticipating the need for an error code
  */
  inline  int
  unpackEBDescription(shared_ptr <set<pr_vof > > & a_soluVolumes,
                      const pr_ebgrap                   & a_graph)

  {

    a_soluVolumes = shared_ptr<set<pr_vof > >(new set<pr_vof >() );
    auto valid= a_graph.validBox();
    for(auto bit = valid.begin(); bit != valid.end(); ++bit)
    {
      if(a_graph->isCovered(*bit))
      {
        //covered bits are not part of the solution so do not go into the set
      }
      else if( a_graph->isRegular())
      {
          
        pr_vof regvol;
        regvol.defineAsRegular(*bit);
        soluVolumes->insert(regvol);
      }
      else
      {
        auto vofs = a_graph->getVoFs(*bit);
        for(int ivof = 0; ivof < vofs.size(); ivof++)
        {
          a_soluVolumes->insert(vofs[ivof]);
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
  */
  inline  int
  makeTheBlobs(shared_ptr<set<EBCMVolu > >           & a_soluVolumes,
               shared_ptr<set<pr_vof   > >           & a_ebVolumes, 
               const EBGraph                         & a_graph,
               const EBData                          & a_data,
               const double                          & a_dx,
               bool                                    a_mergeSmallCells)
  {
    a_soluVolumes = shared_ptr< set< EBCMVolu > >(new set<EBCMVolu>());
    for(auto vof = a_ebVolumes.begin(); vof != a_ebVolumes.end(); ++vof)
    {
      EBCMVolu volume = createEBCMVolume(*vof, a_graph, a_data, a_dx);
      a_soluVolumes.insert(volume);
    }
    
    if(a_mergeSmallCells)
    {
      //return != 0 means we still have small volumes left
      int errval = mergeVolumesOnBox(a_soluVolumes);
      if(errval != 0)
      {
        MayDay::Error("maketheblobs: mergeVolumes on box says we have a small vof somewhere");
      }
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
  sparsifyTheBlobs(vector<EBCMVolu>                 & a_volumes,
                   HostBoxData<int  , 1>            & a_map,
                   std::set<EBCMVolu>               & a_soluVolumes)
  {
  }


  ///Class to make merged sparse embedded boundary description
  /**
    */
  class EBCMGraph
  {
  public:

    EBCMGraph(shared_ptr<pr_ebgrap>  a_graph,
              shared_ptr<pr_ebdata>  a_data,
              Real                   a_dx,
              bool                   a_mergeSmallCells) const:
      m_graph(a_graph), m_data(a_data), m_dx(a_dx)
    {

      shared_ptr< set< pr_vof > >        allEBInd;
      shared_ptr< set< EBCMVolu > > allBlobs;

      ///take all EB stuff and make it into an irregular grid without prejudice about regular or irregular
      //(those distinctions will change here because of merger making regulars into irregulars)
      upackEBDescription(allEBInd, a_graph);
      //put all moment info into a (still fully irregular) a set of volumes
      makeTheBlobs(allBlobs,  allEBInd, a_graph, a_data, a_dx, a_mergeSmallCells);

      //translate irregular graph and datat into sparse  representation that we actually store
      sparsifyTheBlobs(m_volumes, m_map, allBlobs)
      return retval;
                         
    }
  private:
    ///input EB stuff will be necessary if we want to go between descriptions
    //need  dx because we compute regular moment info on the fly.
    shared_ptr<pr_ebgrap>  m_graph;
    shared_ptr<pr_ebdata>  m_data;
    Real                   m_dx;

    //irregular volumes   we do not store unmerged regular cell volumes              
    vector<EBCMVolu> m_volumes;
    //-1 regular, -2 covered, >= 0 is the index into m_volumes.  merged cells both hold the same index.
    HostBoxData<int, 1>   m_map;

    double m_dx;
    Box
    //weak construction is weak
    EBCMGraph();
  };
}


#endif
