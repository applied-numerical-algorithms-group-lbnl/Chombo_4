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
  typedef Chombo::GeometryService
  typedef Proto::pr_vof pr_vof;

  ///
  template<CENTERING cent, int order>
  class EBCMFace
  {
  public:
    //centroid location
    RealVect m_cent;
    pr_vof  m_lo;
    pr_vof  m_hi;
    IndexedMoments<DIM-1, order> m_facmom;
    void defineAsRegular(pr_vof a_lo, pr_vof a_hi, Real a_dx)
    {
      m_lo = a_lo;
      m_hi = a_hi;
      m_facmom.setRegular(a_dx);
    }

  };


  template<int order>
  class EBCMVolu
  {
  public:
    vector<EBCMFace<XFACE> > m_xfaceslo;
    vector<EBCMFace<YFACE> > m_yfaceslo;
    vector<EBCMFace<ZFACE> > m_zfaceslo;
    vector<EBCMFace<XFACE> > m_xfaceshi;
    vector<EBCMFace<YFACE> > m_yfaceshi;
    vector<EBCMFace<ZFACE> > m_zfaceshi;
    EBCMFace<      BOUNDARY> m_ebface;
    
    
    bool                          m_regular;
    //centroid location in global space (not relative to cell center)
    RealVect m_cent;
    IndexedMoments<DIM, order> m_volmom;
    IndexedMoments<DIM, order> m_ebfmom;
    Point    m_pt; //just used as a unique identifier
    

    template <CENTERING cent>
    void fillFaceData(vector< EBCMFace< cent> >                       & a_facedata,
                      const vector< EBIndex<  cent> >                 & a_ebfaces,
                      const HostIrregData<   cent, IndMomSDMinOne, 1> & a_facedat)
    {

      a_facedata.resize(a_ebfaces.size());
      for(int ifac = 0; ifac < a_ebfaces.size(); ifac)
      {
        a_facedata[ifac].m_lo = a_ebfaces[iface].getVoF(Side::Lo);
        a_facedata[ifac].m_hi = a_ebfaces[iface].getVoF(Side::Hi);
        a_facedata[ifac].m_facmom = a_facedat(a_ebfaces[iface], 0);
      }
    }

    void defineAsIrregular(const pr_vof                          & a_vof,
                           const shared_ptr<ch_geoserv>          & a_geoserv,
                           const Box                             & a_domain,
                           const ch_dit                          & a_dit,
                           const double                          & a_dx)
    {
      auto graph   = (*(a_geoserv->getGraphs(   a_domain)))[a_dit];
      auto voludat = (*(a_geoserv->getVoluData( a_domain)))[a_dit];
      auto xfacdat = (*(a_geoserv->getXFaceData(a_domain)))[a_dit];
      auto yfacdat = (*(a_geoserv->getYFaceData(a_domain)))[a_dit];

      auto xebfaclo= graph.getXFaces(Side::Lo);
      auto xebfachi= graph.getXFaces(Side::Hi);
      auto yebfaclo= graph.getYFaces(Side::Lo);
      auto yebfachi= graph.getYFaces(Side::Hi);

      fillFaceData<XFACE>(m_xfaceslo, xebfaclo, xfacdat);
      fillFaceData<YFACE>(m_yfaceslo, yebfaclo, yfacdat);
#if DIM=3                     
      auto zfacdat = (*(a_geoserv->getZFaceData(a_domain)))[a_dit];
      auto zebfaclo= graph.getZFaces(Side::Lo);
      auto zebfachi= graph.getZFaces(Side::Hi);
      fillFaceData<ZFACE>(m_zfaceslo, zebfaclo, zfacdat);
#endif
      m_regular = false;
      m_pt = a_vof.m_pt;
      EBIndex<BOUNDARY> ebf = a_vof.getCutFace();
      auto levvol_ptr = a_geoserv->getVoluData(a_domain);
      auto levebf_ptr = a_geoserv->getEBFaData(a_domain);
      const auto& voldat = (*levvol_ptr)[a_dit];
      const auto& ebfdat = (*levebf_ptr)[a_dit];

      m_volmom = voldat(a_vof, 0);
      m_ebfmom = ebfdat(ebf, 0);;
      defineCentroid(a_dx);
    }

    void defineCentroid(Real a_dx)
    {
      //relies upon volmom  and dx being set
      for(int idir = 0; idir < DIM; idir++)
      {
        Real centerloc = (m_pt[idir] + 0.5)*a_dx;
        Real scaledCent = m_volmom[1+idir];
        m_centroid[idir]= centerloc + (scaledCent/a_dx);
      }
    }
    void defineAsRegular(pr_pt a_pt, double a_dx)
    {
      m_regular = true;
      m_pt = a_pt;
      m_ebfmom.setToZero();
      m_volmom.setRegular(a_dx);
      {
        m_xfaceslo.resize(1);
        m_xfaceshi.resize(1);
        pr_pt lopt = m_pt - Point::Basis(0);
        pr_pt hipt = m_pt + Point::Basis(0);
        m_xfaceslo[0].setToRegular(lopt, m_pt, a_dx);
        m_xfaceshi[0].setToRegular(m_pt, hipt, a_dx);
      }
      {
        m_yfaceslo.resize(1);
        m_yfaceshi.resize(1);
        pr_pt lopt = m_pt - Point::Basis(1);
        pr_pt hipt = m_pt + Point::Basis(1);
        m_xfaceslo[0].setToRegular(lopt, m_pt, a_dx);
        m_xfaceshi[0].setToRegular(m_pt, hipt, a_dx);
      }
#if DIM==3
      {
        m_zfaceslo.resize(1);
        m_zfaceshi.resize(1);
        pr_pt lopt = m_pt - Point::Basis(2);
        pr_pt hipt = m_pt + Point::Basis(2);
        m_zfaceslo[0].setToRegular(lopt, m_pt, a_dx);
        m_zfaceshi[0].setToRegular(m_pt, hipt, a_dx);
      }
#endif      
      defineCentroid(a_dx);
    }
  };


  ///create single-index ebcmvolume from EB description
  template <int order>
  inline EBCMVolu
  createEBCMVolume(const pr_vof                          & a_vof,
                   const shared_ptr<ch_geoserv>          & a_geoserv,
                   const Box                             & a_domain,
                   const ch_dit                          & a_dit,
                   const double                          & a_dx)
  {
    const EBGraph   & a_graph = (*a_geoserv)
    HERE
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
      retval.defineAsIrregular(a_vof, a_geoserv, a_doomain, a_dit, a_dx)
    }

    return retval;
  }
  
  ///
  /**
     unpacks EB description into a fully unstructured set of volumes 
     return value is me anticipating the need for an error code
  */
  inline  int
  unpackEBDescription(shared_ptr <set<pr_vof > >        & a_soluVolumes,
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

  inline int
  mergeVolumesOnBox(shared_ptr<set<EBCMVolu> > & a_soluVolumes)
  {
    //error! not written yet! 
    return -1;
  }
  ///
  /**
     Takes fully unstructured EB data and puts it into EBCM volumes that 
     Blob structures.   
     The set that results is a fully unstructured (and merged, if you ask for it) graph.
     return value = 0 means no errors.
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
        return errval;
      }
    }
    return 0;
  }

  ///Translate irregular graph  and data into the sparse  representation that we actually store.
  /**
     Basically, we just do not add regular volumes and set the integer in the map to the correct thing.
     a_volumes: irregular volumes   we do not store unmerged regular cell volumes              
     a_map: -1= regular -2 covered, >= 0 is the index into m_volumes.  
     Merged cells both hold the same vector index so iterating over this can double count.
     a_soluVolumes: full irregular representation (discard after use).
     
  **/
  inline int
  sparsifyTheBlobs(vector<EBCMVolu>                 & a_volumes,
                   HostBoxData<int  , 1>            & a_map,
                   std::set<EBCMVolu>               & a_soluVolumes)
  {
    for(int ivol = 0; ivol
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
      unpackEBDescription(allEBInd, a_graph);
      //put all moment info into a (still fully irregular) a set of volumes
      makeTheBlobs(allBlobs,  allEBInd, a_graph, a_data, a_dx, a_mergeSmallCells);

      //translate irregular graph and data into sparse  representation that we actually store
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
