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
  using Chombo4::pout;
  using Chombo4::MayDay;
  
  typedef Proto::EBGraph       pr_ebgrap;
  
  typedef Proto::EBIndex<CELL> pr_vof;
  typedef Proto::Point         pr_pt;
  typedef Chombo4::Box         ch_box;
  typedef Chombo4::DataIndex   ch_dit;
  
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
    void defineAsRegular(pr_vof a_vof)
    {
      m_pt =  a_vof.m_pt;

      m_ebfmom.setToZero();
    }

  };

  template<int order>
  class EBCMVolu
  {
  public:
    typedef Chombo4::GeometryService<order> ch_geoserv;
    typedef Proto::IndexedMoments<DIM-1, order> IndMomSDMinOne;
    typedef Proto::IndexedMoments<DIM  , order> IndMomSpaceDim;
    //
    vector<EBCMFace<XFACE, order> >  m_xfaceslo;
    vector<EBCMFace<YFACE, order> >  m_yfaceslo;
    vector<EBCMFace<ZFACE, order> >  m_zfaceslo;
    vector<EBCMFace<XFACE, order> >  m_xfaceshi;
    vector<EBCMFace<YFACE, order> >  m_yfaceshi;
    vector<EBCMFace<ZFACE, order> >  m_zfaceshi;
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

    bool operator < (const EBCMVolu<order>& a_input) const
    {
      return (m_pt < a_input.m_pt);
    }
    vector<Point> m_cells; //all cells marked for this entry.   

    template <CENTERING cent>
    void fillFaceData(vector< EBCMFace< cent, order> >                & a_facedata,
                      const vector< EBIndex<  cent > >                & a_ebfaces,
                      const HostIrregData<   cent, IndMomSDMinOne, 1> & a_facedat)
    {

      a_facedata.resize(a_ebfaces.size());
      for(int iface = 0; iface < a_ebfaces.size(); ++iface)
      {
        a_facedata[iface].m_lo = a_ebfaces[iface].getVoF(Side::Lo);
        a_facedata[iface].m_hi = a_ebfaces[iface].getVoF(Side::Hi);
        a_facedata[iface].m_facmom = a_facedat(a_ebfaces[iface], 0);
      }
    }

    void defineAsIrregular(const pr_vof                          & a_vof,
                           const shared_ptr<ch_geoserv>          & a_geoserv,
                           const ch_box                          & a_domain,
                           const ch_dit                          & a_dit,
                           const double                          & a_dx)
    {
      m_regular = false;
      m_pt = a_vof.m_pt;
      m_dx = a_dx;
      m_cells.resize(1, m_pt);
      const auto & graph   = (*(a_geoserv->getGraphs(   a_domain)))[a_dit];
      const auto & xfacdat = (*(a_geoserv->getXFaceData(a_domain)))[a_dit];
      const auto & yfacdat = (*(a_geoserv->getYFaceData(a_domain)))[a_dit];

      auto xebfaclo= graph.getXFaces(a_vof, Side::Lo);
      auto xebfachi= graph.getXFaces(a_vof, Side::Hi);
      auto yebfaclo= graph.getYFaces(a_vof, Side::Lo);
      auto yebfachi= graph.getYFaces(a_vof, Side::Hi);

      fillFaceData<XFACE>(m_xfaceslo, xebfaclo, xfacdat);
      fillFaceData<XFACE>(m_xfaceshi, xebfachi, xfacdat);
      fillFaceData<YFACE>(m_yfaceslo, yebfaclo, yfacdat);
      fillFaceData<YFACE>(m_yfaceshi, yebfachi, yfacdat);
#if DIM==3                     
      const auto & zfacdat = (*(a_geoserv->getZFaceData(a_domain)))[a_dit];
      auto         zebfaclo= graph.getZFaces(a_vof, Side::Lo);
      auto         zebfachi= graph.getZFaces(a_vof, Side::Hi);
      fillFaceData<ZFACE>(m_zfaceslo, zebfaclo, zfacdat);
      fillFaceData<ZFACE>(m_zfaceshi, zebfachi, zfacdat);
#endif
      m_regular = false;
      m_pt = a_vof.m_pt;
      m_cells.resize(1, m_pt);
      EBIndex<BOUNDARY> ebf = a_vof.getCutFace();
      auto levvol_ptr = a_geoserv->getVoluData(  a_domain);
      auto levebf_ptr = a_geoserv->getEBFaceData(a_domain);
      const auto& voldat = (*levvol_ptr)[a_dit];
      const auto& ebfdat = (*levebf_ptr)[a_dit];

      m_volmom = voldat(a_vof, 0);
      m_irregFaces.resize(1);
      m_irregFaces[0].m_ebfmom = ebfdat(ebf, 0);
      m_irregFaces[0].m_pt = m_pt;
      
      defineDerivedMeasures(a_dx);
    }

    void defineDerivedMeasures(Real a_dx)
    {
      Real cellVol = D_TERM(a_dx, *a_dx, *a_dx);
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
      m_irregFaces.resize(0);
      m_volmom.setRegular(a_dx);
      {
        m_xfaceslo.resize(1);
        m_xfaceshi.resize(1);
        pr_pt lopt = m_pt - Point::Basis(0);
        pr_pt hipt = m_pt + Point::Basis(0);
        pr_vof lovof(lopt, 0);
        pr_vof hivof(hipt, 0);
        pr_vof ptvof(m_pt, 0);
        m_xfaceslo[0].defineAsRegular(lovof, ptvof, a_dx);
        m_xfaceshi[0].defineAsRegular(ptvof, hivof, a_dx);
      }
      {
        m_yfaceslo.resize(1);
        m_yfaceshi.resize(1);
        pr_pt lopt = m_pt - Point::Basis(1);
        pr_pt hipt = m_pt + Point::Basis(1);
        pr_vof lovof(lopt, 0);
        pr_vof hivof(hipt, 0);
        pr_vof ptvof(m_pt, 0);

        m_yfaceslo[0].defineAsRegular(lovof, ptvof, a_dx);
        m_yfaceshi[0].defineAsRegular(ptvof, hivof, a_dx);
      }
#if DIM==3
      {
        m_zfaceslo.resize(1);
        m_zfaceshi.resize(1);
        pr_pt lopt = m_pt - Point::Basis(2);
        pr_pt hipt = m_pt + Point::Basis(2);
        pr_vof hivof(hipt, 0);
        pr_vof ptvof(m_pt, 0);
        m_zfaceslo[0].defineAsRegular(lovof, ptvof, a_dx);
        m_zfaceshi[0].defineAsRegular(ptvof, hivof, a_dx);
      }
#endif      
      defineDerivedMeasures(a_dx);
    }
  };


  ///Class to make merged sparse embedded boundary description
  /**
   */
  template< int order>
  class EBCMGraph
  {
  public:
    /**
       First a bunch of static functions that will probably be useful in other contexts
       They are not standalone because I want to use templated typedefs in declarations.
    **/
    ///Speaking of which,
    typedef Chombo4::GeometryService<order> ch_geoserv;
    
    ///create single-index ebcmvolume from EB description
    static
    inline EBCMVolu<order>
    createEBCMVolume(const pr_vof                          & a_vof,
                     const shared_ptr<ch_geoserv>          & a_geoserv,
                     const ch_box                          & a_domain,
                     const ch_dit                          & a_dit,
                     const double                          & a_dx)
    {
      const auto& graphs = a_geoserv->getGraphs(a_domain);
      const auto& ebgrap = (*graphs)[a_dit];

      EBCMVolu<order> retval;
      if(ebgrap.isCovered(a_vof.m_pt))
      {
        Chombo4::MayDay::Error("cannot create EBCMVolume from a covered cell");
      }
      else if(ebgrap.isRegular(a_vof.m_pt))
      {
        retval.defineAsRegular(a_vof.m_pt, a_dx);
      }
      else
      {
        retval.defineAsIrregular(a_vof, a_geoserv, a_domain, a_dit, a_dx);
      }

      return retval;
    }
  
    ///
    /**
       unpacks EB description into a fully unstructured set of volumes 
       return value is me anticipating the need for an error code
    */
    static
    inline int
    unpackEBDescription(shared_ptr <set<pr_vof > >        & a_soluVolumes,
                        const pr_ebgrap                   & a_graph)

    {

      a_soluVolumes = shared_ptr<set<pr_vof > >(new set<pr_vof >() );
      auto valid= a_graph.validBox();
      for(auto bit = valid.begin(); bit != valid.end(); ++bit)
      {
        auto pt = (*bit);
        if(a_graph.isCovered(pt))
        {
          //covered bits are not part of the solution so do not go into the set
        }
        else if( a_graph.isRegular(pt))
        {
          
          pr_vof regvol;
          regvol.defineAsRegular(pt);
          a_soluVolumes->insert(regvol);
        }
        else
        {
          auto vofs = a_graph.getVoFs(*bit);
          for(int ivof = 0; ivof < vofs.size(); ivof++)
          {
            a_soluVolumes->insert(vofs[ivof]);
          }
        } //end irregular
      } //end loop over box
      return 0;
    } //end function


    static
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
              const auto& volu = *volit;
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
        Chombo4::pout() << "Function findBigNeighbor could not find a big neighbor; returning -1." << endl;
        return -1;
      }
      return 0;
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
    static
    inline void
    getShifts(IndexTM<double, DIM>                 & a_shiftOne,
              IndexTM<double, DIM>                 & a_shiftTwo,
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
      RealVect shiftOne = locMer - locOne;
      RealVect shiftTwo = locMer - locTwo;
      for(int idir = 0; idir < DIM; idir++)
      {
        a_shiftOne[idir] = shiftOne[idir];
        a_shiftTwo[idir] = shiftTwo[idir];
      }
    }
  
    ///
    static
    inline void
    shiftAndMergeVolumeData(EBCMVolu<order>                      & a_volMer,
                            const EBCMVolu<order>                & a_volOne,
                            const EBCMVolu<order>                & a_volTwo)
    {
      IndexTM<double, DIM> shiftOne, shiftTwo;
    
      getShifts(shiftOne, shiftTwo, a_volMer, a_volOne, a_volTwo);
      auto momone = a_volOne.m_volmom;
      auto momtwo = a_volTwo.m_volmom;
      momone.shift(shiftOne);
      momtwo.shift(shiftTwo);
      a_volMer.m_volmom  = momone; 
      a_volMer.m_volmom += momtwo;
    }
  

    ///make two vectors into one long one when merged faces align
    template< CENTERING cent>
    static inline void
    mergeFaces(vector<EBCMFace<      cent, order> >       &  a_mergedFaces,
               const vector<EBCMFace<cent, order> >       &  a_startFaces,
               const vector<EBCMFace<cent, order> >       &  a_bigNeighFaces)
    {
      a_mergedFaces = a_startFaces;
      a_mergedFaces.insert(a_mergedFaces.end(), a_bigNeighFaces.begin(), a_bigNeighFaces.end());
    }

    ///pick the correct faces for low and high if we are side by side
    template< CENTERING cent>
    static inline void
    pickFaces(vector<EBCMFace<      cent, order> >       &  a_mergedFacesLo,
              vector<EBCMFace<      cent, order> >       &  a_mergedFacesHi,
              const vector<EBCMFace<cent, order> >       &  a_hiVoFFacesLo,
              const vector<EBCMFace<cent, order> >       &  a_hiVoFFacesHi,
              const vector<EBCMFace<cent, order> >       &  a_loVoFFacesLo,
              const vector<EBCMFace<cent, order> >       &  a_loVoFFacesHi)
    {
      a_mergedFacesLo = a_loVoFFacesLo;
      a_mergedFacesHi = a_hiVoFFacesHi;
    }

    ///
    static inline void
    shiftAndMergeIrregData(EBCMVolu<order>                      & a_volMer,
                           const EBCMVolu<order>                & a_volOne,
                           const EBCMVolu<order>                & a_volTwo)
    {
      ///This does the same thing as mergeFaces (just slap the two vectors together into the third).
      ///The moments are all wrong, though, and need to be shifted
      a_volMer.m_irregFaces = a_volOne.m_irregFaces;
      a_volMer.m_irregFaces.insert(a_volMer.m_irregFaces.end(), a_volTwo.m_irregFaces.begin(), a_volTwo.m_irregFaces.end());
      IndexTM<double, DIM> shiftOne, shiftTwo;

      int curmer = 0;
      for(int iface = 0; iface < a_volOne.m_irregFaces.size(); iface++)
      {
        a_volMer.m_irregFaces[curmer].m_ebfmom.shift(shiftOne);
        curmer++;
      }
      for(int iface = 0; iface < a_volTwo.m_irregFaces.size(); iface++)
      {
        a_volMer.m_irregFaces[curmer].m_ebfmom.shift(shiftTwo);
        curmer++;
      }
    }
  
    static inline int
    mergeTwoVolumesIntoThird(EBCMVolu<order>                      & a_volMer,
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
        mergeFaces<XFACE>(a_volMer.m_xfaceslo, a_volOne.m_xfaceslo, a_volTwo.m_xfaceslo);
        mergeFaces<XFACE>(a_volMer.m_xfaceshi, a_volOne.m_xfaceshi, a_volTwo.m_xfaceshi);
      }
      if(diff[1] == 0) 
      {
        //No difference in y.  The points are adjacent
        mergeFaces<YFACE>(a_volMer.m_yfaceslo, a_volOne.m_yfaceslo, a_volTwo.m_yfaceslo);
        mergeFaces<YFACE>(a_volMer.m_yfaceshi, a_volOne.m_yfaceshi, a_volTwo.m_yfaceshi);
      }
#if DIM==3
      if(diff[2] == 0) 
      {
        //No difference in z.  The points are adjacent
        mergeFaces<ZFACE>(a_volMer.m_zfaceslo, a_volOne.m_zfaceslo, a_volTwo.m_zfaceslo);
        mergeFaces<ZFACE>(a_volMer.m_zfaceshi, a_volOne.m_zfaceshi, a_volTwo.m_zfaceshi);
      }
#endif
      if(std::abs(diff[0] == 1))
      {
        bool volOneHigher = diff[0] == 1;

        if(volOneHigher)
        {
          pickFaces<XFACE>(a_volMer.m_xfaceslo, a_volMer.m_xfaceshi,
                           a_volOne.m_xfaceslo, a_volOne.m_xfaceshi, 
                           a_volTwo.m_xfaceslo, a_volTwo.m_xfaceshi);
        }
        else
        {
          pickFaces<XFACE>(a_volMer.m_xfaceslo, a_volMer.m_xfaceshi,
                           a_volTwo.m_xfaceslo, a_volTwo.m_xfaceshi, 
                           a_volOne.m_xfaceslo, a_volOne.m_xfaceshi);
        }
      }
      if(std::abs(diff[1] == 1))
      {
        bool volOneHigher = diff[1] == 1;
        if(volOneHigher)
        {
          pickFaces<YFACE>(a_volMer.m_yfaceslo, a_volMer.m_yfaceshi,
                           a_volOne.m_yfaceslo, a_volOne.m_yfaceshi, 
                           a_volTwo.m_yfaceslo, a_volTwo.m_yfaceshi);
        }
        else
        {
          pickFaces<YFACE>(a_volMer.m_yfaceslo, a_volMer.m_yfaceshi,
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
          pickFaces<ZFACE>(a_volMer.m_zfaceslo, a_volMer.m_zfaceshi,
                           a_volOne.m_zfaceslo, a_volOne.m_zfaceshi, 
                           a_volTwo.m_zfaceslo, a_volTwo.m_zfaceshi);
        }
        else
        {
          pickFaces<ZFACE>(a_volMer.m_zfaceslo, a_volMer.m_zfaceshi,
                           a_volTwo.m_zfaceslo, a_volTwo.m_zfaceshi, 
                           a_volOne.m_zfaceslo, a_volOne.m_zfaceshi);
        }
      }

#endif

      return 0;
    }

    static inline int
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
    }    
  
    ///modify a fully irregular description to merge neighboring volumes to eliminate smaller cells
    static inline int
    mergeVolumesOnBox(shared_ptr<set<EBCMVolu<order> > > & a_soluVolumes)
    {
      shared_ptr<set<EBCMVolu<order> > >  returningVolumes(new set<EBCMVolu<order> >());
      shared_ptr<set<EBCMVolu<order> > >  volumesRemaining(new set<EBCMVolu<order> >());
      *volumesRemaining = *a_soluVolumes;
      for(auto volit = a_soluVolumes->begin(); volit != a_soluVolumes->end(); ++volit)
      {
        const auto& volu = *volit;
        if(volumesRemaining->find(volu) != volumesRemaining->end())
        {
          volumesRemaining->erase(volu);
          if(volu.m_regular)
          {
            returningVolumes->insert(volu); //if regular, no need to merge with anything.
          }
          else
          {
            EBCMVolu<order> merged;
            EBCMVolu<order> bigNeigh;
            int neival = mergeWithBigNeighbor(merged, bigNeigh, volu,  volumesRemaining);
            if(!neival)
            {
              volumesRemaining->erase(bigNeigh);
              returningVolumes->insert(merged);
            }
            else
            {
              pout() << "mergeVolumes: mergeWithBigNeighbor failed with error code " << neival << endl;
              return neival;
            }
          }
        }
      }

      //replace input argument with new set.
      *a_soluVolumes = *returningVolumes;
    
      //check to see if we have been successful
      bool stillHaveSmallVolumes = false;
      for(auto volit = a_soluVolumes->begin(); volit != a_soluVolumes->end(); ++volit)
      {
        const auto volu = *volit;
        if(volu.m_kappa < 1.0)
        {
          stillHaveSmallVolumes = true;
        }
      }
      if(stillHaveSmallVolumes)
      {
        pout() << "mergeVolumesOnBox:  Leaving with a small volume still there; returning -1." << endl;
        return -1;
      }
      return 0;
    }

    ///
    /**
       Takes fully unstructured EB data and puts it into EBCM volumes that 
       Blob structures.   
       The set that results is a fully unstructured (and merged, if you ask for it) graph.
       return value = 0 means no errors.
    */
    inline  int
    makeTheBlobs(shared_ptr<set<EBCMVolu <order> > >   & a_soluVolumes,
                 shared_ptr<set<pr_vof   > >           & a_ebVolumes, 
                 const shared_ptr<ch_geoserv>          & a_geoserv,
                 const ch_box                          & a_domain,
                 const ch_dit                          & a_dit,
                 const double                          & a_dx,
                 bool                                    a_mergeSmallCells)
    {
      a_soluVolumes = shared_ptr< set< EBCMVolu<order> > >(new set<EBCMVolu<order> >());
      for(auto vof = a_ebVolumes->begin(); vof != a_ebVolumes->end(); ++vof)
      {
        EBCMVolu<order> volume = createEBCMVolume(*vof, a_geoserv, a_domain, a_dit, a_dx);
        a_soluVolumes->insert(volume);
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
    static inline int
    sparsifyTheBlobs(shared_ptr<vector<EBCMVolu<order> > >     & a_volumes,
                     shared_ptr<HostBoxData<int  , 1>  >       & a_map,
                     shared_ptr< set<EBCMVolu<order> > >       & a_soluVolumes,
                     const pr_ebgrap                           & a_graph)
    {
      a_volumes = shared_ptr<vector<EBCMVolu<order> > >(new vector<EBCMVolu<order>());
      a_map    =  shared_ptr<HostBoxData<    int, 1 > >(new HostBoxData<int, 1>());
      a_map->setVal(-2);
    
      //define the map and set it all to covered.   the input set contains only regular and irregular
      for(auto it = a_soluVolumes->begin(); it != a_soluVolumes->end(); ++it)
      {
        const auto& volu = *it;
        if(volu.m_regular)
        {
          (*a_map)(volu.m_pt, 0) = -1;
        }
        else
        {
          for(int ipt = 0; ipt < volu.m_cells.size(); ipt++)
          {
            //set all merged cells to the vector value.
            (*a_map)(volu.m_cells[ipt], 0) = a_volumes->size();  
            a_volumes0->push_back(volu);
          } //end loop over cell merger
        } //end if irregular
      }
      return 0;
    }





    //strong construction would be better but way less convenient
    EBCMGraph()
    {
      m_defined = false;
    }

    inline void 
    define(const shared_ptr<ch_geoserv>          & a_geoserv,
           const ch_box                          & a_domain,
           const ch_dit                          & a_dit,
           const double                          & a_dx,
           bool                   a_mergeSmallCells) 
    {
      m_defined = true;
      m_dx         = a_dx;
      m_domain     = a_domain;
      m_volumes.resize(0);

      const auto& graphs = a_geoserv->getGraphs(a_domain);
      const auto& ebgrap = (*graphs)[a_dit];
      //keeping the allregular
      if(ebgrap.isAllRegular())
      {
        m_allRegular  = true;
        m_allCovered  = false;
      }
      else if(ebgrap.isAllCovered())
      {
        m_allRegular = false;
        m_allCovered = true;
        
      }
      else
      {
        m_allRegular = false;
        m_allCovered = false;
        shared_ptr< set< pr_vof           > > allEBInd;
        shared_ptr< set< EBCMVolu <order> > > allBlobs;

        ///take all EB stuff and make it into an irregular grid without prejudice about regular or irregular
        //(those distinctions will change here because of merger making regulars into irregulars)
        unpackEBDescription(allEBInd, ebgrap);
        //put all moment info into a (still fully irregular) a set of volumes
        makeTheBlobs(allBlobs,  allEBInd, a_geoserv, a_domain, a_dit, a_dx, a_mergeSmallCells);

        //translate irregular graph and data into sparse  representation that we actually store
        sparsifyTheBlobs(m_volumes, m_map, allBlobs, ebgrap);

      }
    }
    
    //whether define function has been called on this object
    /**
       If this == false, everything else is not defined.
    **/
    bool m_defined;

    ///
    /**
       I am not bothering with fancy accessor functions.   Member data is left in the clear.  Use it wisely.
       Remember the map and the vector are not defined if allRegular or allCovered are true.   This should save a 
       ton of memory where there is not much EB.
    **/
    bool m_allRegular;
    bool m_allCovered;
    //need  dx because we compute regular moment info on the fly.
    Real                   m_dx;
    ch_box                 m_domain;
    /**
       Irregular volume information.   We do not store unmerged regular cell volumes.
       The stuff below is not defined if m_allRegular or m_allCovered are true
    **/
    shared_ptr< vector<EBCMVolu<order> > > m_volumes;
    //-1 regular, -2 covered, >= 0 is the index into m_volumes.  merged cells both hold the same index.
    shared_ptr< HostBoxData<int, 1> >      m_map;


    bool isRegular(ch_pt a_pt) const
    {
      return (m_map(a_pt, 0) == -1);
    }
    bool isCovered(ch_pt a_pt) const
    {
      return (m_map(a_pt, 0) == -2);
    }

 
  }; //end ebcmgraph class

  template <int order>
  class MetaDataLevel
  {
  public:
    typedef Chombo4::GeometryService<order>      ch_geoserv;
    typedef Chombo4::Box                         ch_box;
    typedef Chombo4::DisjointBoxLayout           ch_dbl;
    typedef Chombo4::LevelGeometryService<order> ch_levelgeo;
    typedef Chombo4::NullConstructorDataFactory<   EBCMGraph<order> >      cmgraph_factory_t;    
    typedef CH4_Data_Choreography::DistributedData<EBCMGraph<order> >      cmgraph_distrib_t;
    
    MetaDataLevel(const shared_ptr<ch_geoserv>& a_geoserv,
                  const ch_box                & a_domain,
                  const Real                  & a_dx,
                  bool                          a_mergeSmallCells)
    {
      shared_ptr<ch_levelgeo> eb_description = (*a_geoserv)[a_domain];
      m_grids  = eb_description->getDBL();
      m_domain = a_domain;
      m_dx     = a_dx;
      auto ghost = eb_description->ghostVect();
      m_graphs   = shared_ptr<cmgraph_distrib_t >
        ( new cmgraph_distrib_t(m_grids,  ghost, cmgraph_factory_t()));
      auto dit = m_grids.dataIterator();
      for(int ibox = 0; ibox < dit.size(); ibox++)
      {
        (*m_graphs)[dit[ibox]].define(a_geoserv, a_domain, dit[ibox], a_dx, a_mergeSmallCells);
      }
    }

    shared_ptr<cmgraph_distrib_t>  m_graphs;
    ch_dbl                         m_grids;
    ch_box                         m_domain;
    Real                           m_dx;
  }; //end class EBCMMetaDataLevel


  ///way to get through the volumes on a graph
  template <int order>
  class SubVolumeIterator
  {
  public:
    SubVolumeIter(const EBCMGraph<order> & a_graph, const pr_box& a_region)
    {
      //I cannot just iterate through box points or I double count if there is merging.
      //First we get the  regular volumes that fit the region; we then check the irregular volumes.
      for (auto it = a_region.begin(); it != a_region.end(); ++it)
      {
        if(a_graph.isRegular(*it))
        {
          m_regularPoints.push_back(*it);
        }
      }
      for(int ivec = 0; ivec < a_graph.m_volumes.size(); ivec++)
      {
        const auto& volu = a_graph.m_volumes[ivec];
        if(a_region.contains(volu.m_pt))
        {
          m_irregularIndicies.push_back(ivec);
        }
      }
      m_current = 0;
      m_totalSize = m_irregIndicies.size() + m_regularPoints.size();
      m_dx    = a_graph.m_dx;
      m_graph = a_graph;
    }

    void begin()
    {
      m_current = 0;
    }
    bool ok() const
    {
      return (m_current < m_totalSize);
    }

    EBCMVolu<order> operator*()
    {
      EBCMVolu<order> retval;
      if(m_current < m_regularPoints.size())
      {
        retval.defineAsRegular(m_regularPoints[m_current], m_dx);
        reurn retval;
      }
      else if(m_current < m_totalSize())
      {
        unsigned int irregInd = m_current - m_regularPoints.size();
        unsigned int graphind = m_irregIndicies[irregInd];
        return (*m_graph.m_volumes)[graphind];
      }
      else
      {
        ch_mayday("bogus dude");
      }
      //should never get here but this will shut up compilers
      return retval;
    }
    
  private:
    SubVolumeIterator();
    
    EBCMGraph    m_graph;
    unsigned int m_current, m_totalSize;
    vector<ch_pt>         m_regularPoints;
    vector<unsigned int>  m_irregIndicies;
    double                m_dx;
  }
} //end namespace EBCM


#endif
