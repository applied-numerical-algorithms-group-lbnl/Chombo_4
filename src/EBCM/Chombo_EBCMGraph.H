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

    pr_vof  m_lo;
    pr_vof  m_hi;
    IndexedMoments<DIM-1, order> m_facmom;
    void defineAsRegular(pr_vof a_lo, pr_vof a_hi, Real a_dx)
    {                                                  // 
      m_lo = a_lo;
      m_hi = a_hi;
      m_facmom.setRegular(a_dx);
    }

  };

  template<int order>
  class IrregFace
  {
  public:

    pr_pt m_pt;
    IndexedMoments<DIM, order> m_ebfmom;
    void defineAsRegular(pr_vof a_vof);
    {
      m_pt =  a_vof.m_pt;

      m_ebfmom.setToZero();
    }

  };
  

  template<int order>
  class EBCMVolu
  {
  public:
    //
    vector<EBCMFace<XFACE> >  m_xfaceslo;
    vector<EBCMFace<YFACE> >  m_yfaceslo;
    vector<EBCMFace<ZFACE> >  m_zfaceslo;
    vector<EBCMFace<XFACE> >  m_xfaceshi;
    vector<EBCMFace<YFACE> >  m_yfaceshi;
    vector<EBCMFace<ZFACE> >  m_zfaceshi;
    vector<IrregFace<order> > m_irregFaces;
    
    bool                          m_regular;
    //centroid location in global space (not relative to cell center)
    RealVect m_centroid;
    Real     m_kappa; //volume fraction
    IndexedMoments<DIM, order> m_volmom;

    //used as a unique identifier.  
    //Also volume moments are centered around this cell center.  Hence dx.
    Point    m_pt; 
    Real     m_dx;

    vector<Point> m_cells; //all cells marked for this entry.   

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
      m_regular = false;
      m_pt = a_vof.m_pt;
      m_dx = a_dx;
      m_cells.resize(1, m_pt);
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
      m_cells.resize(1, m_pt);
      EBIndex<BOUNDARY> ebf = a_vof.getCutFace();
      auto levvol_ptr = a_geoserv->getVoluData(a_domain);
      auto levebf_ptr = a_geoserv->getEBFaData(a_domain);
      const auto& voldat = (*levvol_ptr)[a_dit];
      const auto& ebfdat = (*levebf_ptr)[a_dit];

      m_volmom = voldat(a_vof, 0);
      m_irregfaces.resize(1);
      m_irregfaces[0].m_ebfmom = ebfdat(ebf, 0);
      m_irregfaces[0].m_pt = m_pt;
      
      defineDerivedMeasures(a_dx);
    }

    void defineDerivedMeasures(Real a_dx)
    {
      Real cellVol = D_TERM(a_dx; *a_dx; *a_dx);
      m_kappa = m_volmom[0]/cellVol;
      //relies upon volmom  and dx being set
      for(int idir = 0; idir < DIM; idir++)
      {
        cellVol *= a_dx;
        Real centerloc = (m_pt[idir] + 0.5)*a_dx;
        Real scaledCent = m_volmom[1+idir];
        m_centroid[idir]= centerloc + scaledCent;
      }
    }
    void defineAsRegular(pr_pt a_pt, double a_dx)
    {
      m_regular = true;
      m_pt = a_pt;
      m_dx = a_dx;
      m_cells.resize(1, m_pt);
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
      defineDerivedMeasures(a_dx);
    }
  };


  ///create single-index ebcmvolume from EB description
  template <int order>
  inline EBCMVolu<order>
  createEBCMVolume(const pr_vof                          & a_vof,
                   const shared_ptr<ch_geoserv>          & a_geoserv,
                   const Box                             & a_domain,
                   const ch_dit                          & a_dit,
                   const double                          & a_dx)
  {
    const EBGraph   & a_graph = (*a_geoserv);

    EBCMVolu<order> retval;
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
      retval.defineAsIrregular(a_vof, a_geoserv, a_doomain, a_dit, a_dx);
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
      } //end irregular
    } //end loop over box
      return 0;
  } //end function


  template<int order>
  inline int
  findBigNeighbor(EBCMVolu<order>                      & a_bigNeigh,
                  const EBCMVolu<order>                & a_start, 
                  shared_ptr< set<EBCMVolu<order> > >    a_remaining)
  {
    bool found = false;
    for(int idir = 0; idir < DIM; idir++)
    {
      for(SideIterator sit; sit.ok(); ++sit)
      {
        int isign  = sign(sit());
        pr_pt newpt = a_start.m_pt + pr_pt::Basis(idir*isign);
        if(!found)
        {
          for(auto volit = a_remaining->begin(); volit != a_remaining->end(); ++volit)
          {
            const auto& volu = *it;
            Real sumkappa = a_start.m_kappa + volu.m_kappa;
            if( (volu.m_pt == newpt) && (volu.m_regular || (sumkappa >= 1.)))
            {
              found = true;
              a_bigNeigh = volu;
            }
          }
        }
      }
    }
    if(!found)
    {
      pout() << "Function findBigNeighbor could not find a big neighbor; returning -1." << endl;
      return -1;
    }
    return 0;
  }

  

  template<int order, CENTERING cent>
  inline void
  mergeFaces(vector<EBCMFace<cent> >       &  a_mergedFacesLo,
             const vector<EBCMFace<cent> > &  a_startFacesLo,
             const vector<EBCMFace<cent> > &  a_bigNeighFacesLo)
  {
    a_mergedFacesLo = a_startFaces;
    a_mergedFacesLo.insert(a_mergedFacesLo.end(), a_bigNeighFaces.start(), a_bigNeighFaces.end());
  }

  template<int order, CENTERING cent>
  inline void
  pickFaces(vector<EBCMFace<cent> >       &  a_mergedFacesLo,
             vector<EBCMFace<cent> >       &  a_mergedFacesHi,
             const vector<EBCMFace<cent> > &  a_hiVoFFacesLo,
             const vector<EBCMFace<cent> > &  a_hiVoFFacesHi,
             const vector<EBCMFace<cent> > &  a_loVoFFacesLo,
             const vector<EBCMFace<cent> > &  a_loVoFFacesHi)
  {
    a_mergedFacesLo = a_loVoFFacesLo;
    a_mergedFacesHi = a_hiVoFFacesHi;
  }
  
  ///code that fills EBCMVolu objects
  /**
     Face objects contain moment informaiton
     Data that needs to be filled:

     vector<EBCMFace<XFACE> >  m_xfaceslo;
     vector<EBCMFace<YFACE> >  m_yfaceslo;
     vector<EBCMFace<ZFACE> >  m_zfaceslo;
     vector<EBCMFace<XFACE> >  m_xfaceshi;
     vector<EBCMFace<YFACE> >  m_yfaceshi;
     vector<EBCMFace<ZFACE> >  m_zfaceshi;
     vector<IrregFace<order> > m_irregFaces;

     //only true if it started as a regular cell and never got merged
     bool                          m_regular;

     //centroid location in global space (not relative to cell center)
     RealVect m_centroid;
     Real     m_kappa; //volume fraction
     IndexedMoments<DIM, order> m_volmom;

    //used as a unique identifier.  
     //Also volume moments are centered around this cell center.  Hence dx.
     Point    m_pt; 
     Real     m_dx
     vector<Point> m_cells; //all cells marked for this entry.   
  */


  ///  get how far to shift moments when merging two volumes
  template<int order>
  inline void
  getShifts(RealVect                             & a_shiftOne,
            RealVect                             & a_shiftTwo,
            const EBCMVolu<order>                & a_volMer,
            const EBCMVolu<order>                & a_volOne,
            const EBCMVolu<order>                & a_volTwo)
  {
    RealVect locMer, locOne, locTwo;
    for(int idir = 0; idir < DIM; idir++)
    {
      locMer[idir] = a_volMer.m_dx*(a_volMer.m_pt[idir] + 0.5);
      locOne[idir] = a_volOne.m_dx*(a_volOne.m_pt[idir] + 0.5);
      locTwo[idir] = a_volTwo.m_dx*(a_volTwo.m_pt[idir] + 0.5);
    }
    a_shiftOne = locMer - locOne;
    a_shiftTwo = locMer - locTwo;
  }

  ///
  template<int order>
  inline void
  shiftAndMergeVolumeData(EBCMVolu<order>                      & a_volMer,
                          const EBCMVolu<order>                & a_volOne,
                          const EBCMVolu<order>                & a_volTwo)
  {
    RealVect shiftOne, shiftTwo;
    
    getShifts(shiftOne, shiftTwo, a_volMer, a_volOne, a_volTwo)
    auto momone = a_volOne.m_volmom;
    auto momtwo = a_volTwo.m_volmom;
    momone.shift(shiftOne);
    momtwo.shift(shiftTwo);
    a_volMer.m_volmom = momone + momtwo;
  }
  

  ///
  template<int order>
  inline void
  shiftAndMergeIrregData(EBCMVolu<order>                      & a_volMer,
                         const EBCMVolu<order>                & a_volOne,
                         const EBCMVolu<order>                & a_volTwo)
  {
    ///this does the same thing as mergeFaces (just slap the two vectors together into the third),
    ///the moments are all wrong, though, and need to be shifted
    a_volMer.m_irregFaces = a_volOne.m_irregFaces;
    a_volMer.m_irregFaces.insert(a_volMer.m_irrregFaces.end(), a_volTwo.m_irregFaces.start(), a_volTwo.m_irregFaces.end());
    RealVect shiftOne, shiftTwo;
    
    getShifts(shiftOne, shiftTwo, a_volMer, a_volOne, a_volTwo)
    auto momone = a_volOne.m_volmom;
    auto momtwo = a_volTwo.m_volmom;
    //remember that we put the volones in first
    int volOneSize = a_volOne.size(); //usually one or zero
    int curVol = 0;
    for(auto it = a_volMer.m_irregFaces.begin();  it != a_volMer.m_irregFaces.end(); ++it)
    {
      auto& curmom = *it;
      if(curVol < volOneSize)
      {
        curmom.shift(shiftOne);
      }
      else
      {
        curmom.shift(shiftTwo);
      }
    }
  }
  
  template<int order>
  inline int
  mergeTwoVolumessIntoThird(EBCMVolu<order>                      & a_volMer,
                            const EBCMVolu<order>                & a_volOne,
                            const EBCMVolu<order>                & a_volTwo)
  {
    a_volMer.m_regular = false;
    a_volMer.m_pt = a_volTwo.m_pt; //just having   the bigger one have the point
    a_volMer.m_cells = a_volOne.m_cells;
    a_volMer.m_cells.insert(a_volMer.m_cells.end(), a_volTwo.m_cells.begin(), a_volTwo.m_cells.end());
    shiftAndMergeVolumeData(a_volMer, a_volOne, a_volTwo);
    shiftAndMergeIrregData( a_volMer, a_volOne, a_volTwo);
    
    pr_pt diff = a_volOne.m_pt - a_volTwo.m_pt;
    if(diff[0] == 0) 
    {
      //No difference in x.  The points are adjacent
      mergeFaces<order, XFACE>(a_volMer.m_xfaceslo, a_volOne.m_xfaceslo, a_volTwo.m_xfaceslo);
      mergeFaces<order, XFACE>(a_volMer.m_xfaceshi, a_volOne.m_xfaceshi, a_volTwo.m_xfaceshi);
    }
    if(diff[1] == 0) 
    {
      //No difference in x.  The points are adjacent
      mergeFaces<order, YFACE>(a_volMer.m_yfaceslo, a_volOne.m_yfaceslo, a_volTwo.m_yfaceslo);
      mergeFaces<order, YFACE>(a_volMer.m_yfaceshi, a_volOne.m_yfaceshi, a_volTwo.m_yfaceshi);
    }
#if DIM==3
    if(diff[2] == 0) 
    {
      //No difference in x.  The points are adjacent
      mergeFaces<order, ZFACE>(a_volMer.m_zfaceslo, a_volOne.m_zfaceslo, a_volTwo.m_zfaceslo);
      mergeFaces<order, ZFACE>(a_volMer.m_zfaceshi, a_volOne.m_zfaceshi, a_volTwo.m_zfaceshi);
    }
#endif
    if(std::abs(diff[0] == 1))
    {
      bool volOneHigher = diff[0] == 1;

      if(volOneHigher)
      {
        pickFaces<order, XFACE>(a_volMer.m_xfaceslo, a_volMer.m_xfaceshi,
                                a_volOne.m_xfaceslo, a_volOne.m_xfaceshi, 
                                a_volTwo.m_xfaceslo, a_volTwo.m_xfaceshi);
      }
      else
      {
        pickFaces<order, XFACE>(a_volMer.m_xfaceslo, a_volMer.m_xfaceshi,
                                a_volTwo.m_xfaceslo, a_volTwo.m_xfaceshi, 
                                a_volOne.m_xfaceslo, a_volOne.m_xfaceshi);
      }
    }
    if(std::abs(diff[1] == 1))
    {
      bool volOneHigher = diff[1] == 1;
      if(volOneHigher)
      {
        pickFaces<order, YFACE>(a_volMer.m_yfaceslo, a_volMer.m_yfaceshi,
                                a_volOne.m_yfaceslo, a_volOne.m_yfaceshi, 
                                a_volTwo.m_yfaceslo, a_volTwo.m_yfaceshi);
      }
      else
      {
        pickFaces<order, YFACE>(a_volMer.m_yfaceslo, a_volMer.m_yfaceshi,
                                a_volTwo.m_yfaceslo, a_volTwo.m_yfaceshi, 
                                a_volOne.m_yfaceslo, a_volOne.m_yfaceshi);
      }
    }
#if DIM==3
    if(std::abs(diff[2] == 1))
    {
      bool volOneHigher = diff[2] == 1;
      if(volOneHigher)
      {
        pickFaces<order, ZFACE>(a_volMer.m_zfaceslo, a_volMer.m_zfaceshi,
                                a_volOne.m_zfaceslo, a_volOne.m_zfaceshi, 
                                a_volTwo.m_zfaceslo, a_volTwo.m_zfaceshi);
      }
      else
      {
        pickFaces<order, ZFACE>(a_volMer.m_zfaceslo, a_volMer.m_zfaceshi,
                                a_volTwo.m_zfaceslo, a_volTwo.m_zfaceshi, 
                                a_volOne.m_zfaceslo, a_volOne.m_zfaceshi);
      }
    }

#else    

    return 0;
 }

  template<int order>
  inline int
  mergeWithBigNeighbor(EBCMVolu<order>                      & a_volMer,
                       EBCMVolu<order>                      & a_bigNeigh,
                       const EBCMVolu<order>                & a_start, 
                       shared_ptr< set<EBCMVolu<order> > >    a_remaining)
  {
    {
      int neival = findBigNeighbor(a_bigNeigh, a_start, a_remaining);
      if(!neival)
      {
        pout() << "mergeWithBigNeighbor: findBigNeighbor returned with error code " << neival << endl;
        return neival;
      }
    }
    {
      int mergeval = mergeTwoVolumesIntoThird(a_volMer, a_bigNeigh, a_start);
      if(!mergeval)
      {
        pout() << "Function mergeTwoVolumesIntoThird failed; returning error code " << mergeval << endl;
        return mergeval;
      }
  }
    
  ///modify a fully irregular description to merge neighboring volumes to eliminate smaller cells
  template <int order>
  inline int
  mergeVolumesOnBox(shared_ptr<set<EBCMVolu<order> > > & a_soluVolumes)
  {
    shared_ptr<set<EBCMVolu<order> > >     returnVolumes(new set<EBCMVolu<order> >());
    shared_ptr<set<EBCMVolu<order> > >  volumesRemaining(new set<EBCMVolu<order> >());
    *volumesRemaining = *a_soluVolumes;
    for(auto volit = a_soluVolumes->begin(); volit != a_soluVolumes->end(); ++volit)
    {
      const auto& volu = *it;
      if(volumesRemaining->find(volu) != volumeRemaining.end())
      {
        volumesRemaining->erase(volu);
        if(volu.m_regular)
        {
          returnVolumes->insert(volu);
        }
        else
        {
          EBCMVolu<order> merged;
          EBCMVolu<order> bigNeigh;
          int neival = mergeWithBigNeighbor(merged, bigNeigh, volu,  volumesRemaining);
          if(!neival)
          {
            volumesRemaining.erase(bigNeigh);
            returnVolumes.insert(merged);
          }
          else
          {
            pout() << "mergeVolumes: mergeWithBigNeighbor failed with error code " << retval << endl;
            return neival;
          }
        }
      }
    }

    //replace input argument with new set.
    *a_soluVolumes = *returnVolumes;
    
    //check to see if we have been successful
    bool stillHaveSmallVolumes = false;
    for(auto volit = a_soluVolumes.begin(); volit != a_soluVolumes.end(); ++volit)
    {
      const auto volu = *it;
      if(volu.kappa < 1.0)
      {
        stillHaveSmallVolumes = true;
      }
    }
    if(stillHaveSmallVolumes)
    {
      pout() << "mergeVolumesOnBox:  Leaving with a small volume still there; returning -1." << endl;
      return -1;
    }
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
                   std::set<EBCMVolu>               & a_soluVolumes,
                   const pr_ebgrap                  & a_graph))
  {
    a_volumes.resize(0);
    a_map.define(a_graph.box());
    a_map.setVal(-2);
    
    //define the map and set it all to covered.   the input set contains only regular and irregular
    for(auto it = a_soluVolumes.begin(); it != soluVolumes.end(); ++it)
    {
      const auto& volu = *it;
      if(volu.m_regular)
      {
        a_map(volu.m_pt, 0) = -1;
      }
      else
      {
        for(int ipt = 0; ipt < volu.m_cells.size(); ipt++)
        {
          a_map(volu.m_cells[ipt], 0) = a_volumes.size();  //set all merged cells to the vector value.
          a_volumes.push_back(volu);
        } //end loop over cell merger
      } //end if irregular
    } //end loop over solution volumes
  }//end function


///Class to make merged sparse embedded boundary description
/**
 */
template< int order>
class EBCMGraph
{
  public:

    EBCMGraph(shared_ptr<pr_ebgrap>  a_graph,
              shared_ptr<pr_ebdata>  a_data,
              Real                   a_dx,
              bool                   a_mergeSmallCells) const:
      m_graph(a_graph), m_data(a_data), m_dx(a_dx)
    {

      shared_ptr< set< pr_vof > >   allEBInd;
      shared_ptr< set< EBCMVolu > > allBlobs;

      ///take all EB stuff and make it into an irregular grid without prejudice about regular or irregular
      //(those distinctions will change here because of merger making regulars into irregulars)
      unpackEBDescription(allEBInd, a_graph);
      //put all moment info into a (still fully irregular) a set of volumes
      makeTheBlobs(allBlobs,  allEBInd, a_graph, a_data, a_dx, a_mergeSmallCells);

      //translate irregular graph and data into sparse  representation that we actually store
      sparsifyTheBlobs(m_volumes, m_map, allBlobs, a_graph)
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
