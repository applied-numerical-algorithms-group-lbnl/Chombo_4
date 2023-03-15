#ifndef __Chombo_EBCMGraph__
#define __Chombo_EBCMGraph__

#include "implem/Proto_IrregGraph.H"
#include "EBProto.H"
#include "Proto_RealVect.H"
#include "Chombo_ProtoInterface.H"
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
      auto startpt    = a_start.m_pt;
      auto startkappa = a_start.m_kappa;
      bool found = false;
      for(int idir = 0; idir < DIM; idir++)
      {
        for(SideIterator sit; sit.ok(); ++sit)
        {
          int isign  = sign(sit());
          pr_pt newpt = a_start.m_pt + pr_pt::Basis(idir, isign);
          if(!found)
          {
            for(auto volit = a_remaining->begin(); volit != a_remaining->end(); ++volit)
            {
              const auto& volu = *volit;
              if(volu.m_pt == newpt)
              {
                if(volu.m_regular)
                {
                  Chombo4::pout() << "found one!" << endl;
                  found = true;
                  a_bigNeigh = volu;
                }
                else
                {
                  auto neikappa =  volu.m_kappa;
                  auto sumkappa = startkappa + neikappa;
                  if(sumkappa >= 1.0)
                  {
                    Chombo4::pout() << "found one!" << endl;
                    found = true;
                    a_bigNeigh = volu;
                  }
                }
              }
            }
          }
        }
      }
      if(!found)
      {
        
        Chombo4::pout() << "Function findBigNeighbor:a_start.m_pt = " << a_start.m_pt << endl;
        Chombo4::pout() << "Function findBigNeighbor:could not find a big neighbor; returning -1." << endl;
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
    
    ///
    static
    inline EBCMVolu<order>
    mergeManyIntoOne(vector<EBCMVolu<order>& a_many)
    {
      EBCMVolu<order> retval;
      if(a_many.size() == 1)
      {
        retval = a_many[0];
      }
      else
      {
        EBCMVolu<order> third = mergeTwoVolumesIntoThird(a_volsLeft[0], a_volsLeft[1]);
        //we reduce the size of the vector by one with every call.
        vector<EBCMVolu<order> > volsLeft(a_many.size() - 1);
        volsLeft[0] = third;
        //-2 below because we just replaced the first two
        for(int ivol = 0; ivol < a_many.size() - 2; ivol++)
        {
          //+1 because  used one slot for the merged one
          volsLeft[ivol] = a_many[ivol+1];
        }
        //the rare instance where recusion is just way easier.
        retval = mergeManyIntoOne(volsLeft);
      }

      return retval;
    }
    ///
    static
    inline EBCMVolu<order>
    getGhostedVolume(vector<pr_vof>                        & a_allvofs,
                     double                                & a_sumkappa,
                     const  pr_vof                         & a_startvof,
                     const shared_ptr<ch_geoserv>          & a_geoserv,
                     const ch_box                          & a_domain,
                     const ch_dit                          & a_dit,
                     const double                          & a_dx,
                     const shared_ptr<set<pr_vof   > >     & a_vofsLeft,
                     const pr_ebgrap                       & a_graph,
                     int numGhost)
    {
      EBCMVolu retval;
      pr_box ptbox(a_startvol.m_pt, a_startvol.m_pt);
      pr_box grownbox = ptbox.grow(1);
      grownbox &= a_graph.getRegion();
      a_allvofs.resize(0);
      for(auto bit = grownbox.begin(); bit != grownbox.end(); ++bit)
      {
        auto pt = *bit;
        vector<pr_vof> vofs = a_graph.getVoFs(pt);
        for(int ivof = 0; ivof < vofs.size(); ivof++)
        {
          //have to be certain to not double count here.
          if(a_vofsLeft.find(vofs[ivof]) != a_vofsLeft.end())
          {
            a_allvofs.push_back(vofs[ivof;]
          }
        }
      }
      
      vector<EBCMVolu<order> > allvols(a_allvofs.size());
      for(int ivol = 0; ivol <         a_allvofs.size(); ivol++)
      {
        auto curpt = a_allvofs[ivol].m_pt;
        if(a_graph.isCovered(curpt))
        {
          Chombo4::MayDay::Error("cannot create EBCMVolume from a covered cell");
        }
        else if(a_graph.isRegular(curpt))
        {
          allvols[ivol].defineAsRegular(curpt, a_dx);
        }
        else
        {
          allvols[ivol].defineAsIrregular(allvofs[ivof], a_geoserv, a_domain, a_dit, a_dx);
        }
      }

      //make two into one until there are no more
      retval = mergeManyIntoOne(allvols[ivol], 0);
      Chombo4::MayDay::Error("not implemented");

      return retval;
    }
    
    ///
    static
    inline EBCMVolu<order>
    getMergedVolume(vector<pr_vof>                        & a_allvofs,
                    const  pr_vof                         & a_startvof,
                    const shared_ptr<ch_geoserv>          & a_geoserv,
                    const ch_box                          & a_domain,
                    const ch_dit                          & a_dit,
                    const double                          & a_dx,
                    const shared_ptr<set<pr_vof   > >     & a_vofsLeft,
                    const pr_ebgrap                       & a_graph)
    {
      double sumkappa = 0;
      int numghost = 1;
      EBCMVolu<order> retval;
      {
        int numghost = 1;
        retval = getGhostedVolume(a_allvofs, sumkappa, a_startvof, a_geoserv, a_domain, a_dit, a_dx, a_vofsLeft, a_graph,numghost);
      }
      if(sumkappa < 1)
      {
        int numghost = 2;
        retval = getGhostedVolume(a_allvofs, sumkappa, a_startvof, a_geoserv, a_domain, a_dit, a_dx, a_vofsLeft, a_graph,numghost);
      }
      if(sumkappa < 1)
      {
        MayDay::Error("could not find a full cell even with ghost rad = 2");
      }
      return retval;
    }
    ///
    static
    inline EBCMVolu<order>
    createEBCMVolume(vector<pr_vof>                        & a_mergevec,
                     const pr_vof                          & a_vof,
                     const shared_ptr<ch_geoserv>          & a_geoserv,
                     const ch_box                          & a_domain,
                     const ch_dit                          & a_dit,
                     const double                          & a_dx,
                     const shared_ptr<set<pr_vof   > >     & a_vofsLeft,
                     const pr_ebgrap                       & a_graph,
                     const bool                            & a_mergeSmallCells)
    {
      EBCMVolu<order> startvol, retval;
      a_mergevec.resize(1, a_vof);
      
      if(a_graph.isCovered(a_vof.m_pt))
      {
        Chombo4::MayDay::Error("cannot create EBCMVolume from a covered cell");
      }
      else if(a_graph.isRegular(a_vof.m_pt))
      {
        startvol.defineAsRegular(a_vof.m_pt, a_dx);
      }
      else
      {
        startvol.defineAsIrregular(a_vof, a_geoserv, a_domain, a_dit, a_dx);
      }
      
      if( (startvol.m_kappa == 1.) || (!a_mergeSmallCells))
      {
        retval = startvol;
      }
      else
      {
        retval = getMergedVolume(a_mergevec,
                                 a_vof, 
                                 a_geoserv,        
                                 a_domain,         
                                 a_dit,            
                                 a_dx,             
                                 a_vofsLeft,       
                                 a_graph);
                                 
      }
      return retval;
    }
    ///
    /**
       Takes fully unstructured EB data and puts it into EBCM volumes that 
       Blob structures.   
       The set that results is a fully unstructured (and merged, if you ask for it) graph.
       return value = 0 means no errors.
    */
    static inline  int
    makeTheBlobs(shared_ptr<vector<EBCMVolu<order> > >     & a_volumes,
                 shared_ptr<HostBoxData<int  , 1>  >       & a_map,
                 shared_ptr<set<pr_vof   > >               & a_inputVoFs, 
                 const shared_ptr<ch_geoserv>              & a_geoserv,
                 const ch_box                              & a_domain,
                 const ch_dit                              & a_dit,
                 const double                              & a_dx,
                 const pr_ebgrap                           & a_graph,
                 bool                                        a_mergeSmallCells)
    {
      a_volumes = shared_ptr<vector<EBCMVolu<order> > >(new vector<EBCMVolu<order> > ());
      a_map    =  shared_ptr<HostBoxData<    int, 1 > >(new HostBoxData<int, 1>      ());
      //define the map and set it all to covered.   the input set contains only regular and irregular
      a_map->define(a_graph.getRegion());
      a_map->setVal(-2);
      shared_ptr<set<pr_vof   > >  vofsLeft(new set<pr_vof>() );
      *vofLeft = a_inputVoFs;
      
      for(auto vofit = a_ebVolumes->begin(); vofit != a_ebVolumes->end(); ++vofit)
      {
        auto vof = *vofit;
        if(a_graph->isRegular(vof.m_pt))
        {
          //only reset it if it is set to covered.   It can also be reset to the vector value if merged
          int mapval = (*a_map)(vof.m_pt, 0);
          if(mapval == -2)
          {
            //again, can be reset in the irreg section
            (*a_map)(vof.m_pt, 0) = -1;
          }
        }
        else if(a_graph->isIrregular(vof.m_pt))
        {
          if(vofsLeft->find(vof) != vofsLeft->end())
          {
            vector<pr_vof> mergevec;
            EBCMVolu<order> volume = createEBCMVolume(mergevec, *vof, , a_geoserv, a_domain, a_dit, a_dx,
                                                      vofsLeft, a_graph, a_mergeSmallCells);
            int vecInd = a_volumes.size();
            for(int ivec = 0; ivec < mergevec.size(); ivec++)
            {
              (*a_map)(mergevec[ivec], 0) = vecInd;
            }
            a_volumes->push_back(volume);
          }
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
           bool                                    a_mergeSmallCells) 
    {
      m_defined = true;
      m_dx         = a_dx;
      m_domain     = a_domain;

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
        shared_ptr< set< pr_vof > > allEBInd;

        ///take all EB stuff and make it into an irregular grid without prejudice about regular or irregular
        //(those distinctions will change here because of merger making regulars into irregulars)
        unpackEBDescription(allEBInd, ebgrap);
        
        //put all moment info into a (still fully irregular) a set of volumes
        //translate all this mixed stuff into sparse  representation that we actually store
        makeTheBlobs(m_volumes, m_map, allBlobs, ebgrap, allEBInd, a_geoserv, a_domain, a_dit, a_dx, a_mergeSmallCells);

      }
    }
    
    //whether define function has been called on this object
    /**
       If this == false, everything else is not defined.
    **/
    bool m_defined;

    ///
    /**
       I am not bothering with fancy accessor functions.   Member data is left in the clear.  
       Use it wisely.
       Remember the map and the vector are not defined if allRegular or allCovered are true.   
       This should save a ton of memory where there is not much EB.
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

    ///
    bool isRegular(pr_pt a_pt) const
    {
      bool retval;
      if(m_allRegular)
      {
        retval = true;
      }
      else if (m_allCovered)
      {
        retval = false;
      }
      else
      {
        
        retval = ((*m_map)(a_pt, 0) == -1);
      }
      return retval;
    }

    ///
    bool isCovered(pr_pt a_pt) const
    {
      bool retval;
      if(m_allRegular)
      {
        retval = false;
      }
      else if(m_allCovered)
      {
        retval = true;
      }
      else
      {
        retval = ((*m_map)(a_pt, 0) == -2);
      }
      return  retval;
    }

    bool hasIrregular() const
    {
      bool     regularOrCovered = m_allRegular  || m_allCovered;
      return (!regularOrCovered);
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


  ///
  /**
     Given a specfied subbox covered by EBCMGraph, this iterator will move through the volumes 
     within the subbox.   They are divide by m_pt so this could get weird when merging happens.
     This can be used to find the volumes within the neighborhood of a particular cell.
   */
  template <int order>
  class SubVolumeIterator
  {
  public:
    SubVolumeIterator(const EBCMGraph<order> & a_graph, const ch_box& a_region)
    {
      m_regularPoints.resize(0);
      m_irregIndicies.resize(0);
      Proto::Box prregion = ProtoCh::getProtoBox(a_region);
      for (auto it = prregion.begin(); it != prregion.end(); ++it)
      {
        if(a_graph.isRegular(*it))
        {
          m_regularPoints.push_back(*it);
        }
      }
      if(a_graph.hasIrregular())
      {
        for(int ivec = 0; ivec < a_graph.m_volumes->size(); ivec++)
        {
          const auto& volu = (*(a_graph.m_volumes))[ivec];
          if(a_region.contains(volu.m_pt))
          {
            m_irregIndicies.push_back(ivec);
          }
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
        return retval;
      }
      else if(m_current < m_totalSize)
      {
        unsigned int irregInd = m_current - m_regularPoints.size();
        unsigned int graphind = m_irregIndicies[irregInd];
        return (*m_graph.m_volumes)[graphind];
      }
      else
      {
        Chombo4::MayDay::Error("subvolumeiterator: bogus ebcm iteration");
      }
      //should never get here but this will shut up compilers
      return retval;

    }

    void operator++()
    {
      m_current++;
    }
    
  private:
    SubVolumeIterator();
    
    EBCMGraph<order>    m_graph;
    unsigned int m_current, m_totalSize;
    vector<pr_pt>         m_regularPoints;
    vector<unsigned int>  m_irregIndicies;
    double                m_dx;
  }; //end class subvolumeiterator
} //end namespace EBCM


#endif
