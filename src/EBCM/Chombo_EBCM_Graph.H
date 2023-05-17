#ifndef __Chombo_EBCM_Graph__
#define __Chombo_EBCM_Graph__

#include "implem/Proto_IrregGraph.H"
#include "EBProto.H"
#include "Proto_RealVect.H"
#include "Chombo_ProtoInterface.H"
using std::vector;
using std::set;

namespace EBCM
{
  using std::set;
  using std::endl;
  using std::shared_ptr;
  using Chombo4::pout;
  using Chombo4::MayDay;
//  using CH4_SPMD::linearIn;  
//  using CH4_SPMD::linearOut;  
//  using CH4_SPMD::linearSize;
  
  typedef Proto::EBGraph                           pr_ebgrap;
  typedef Proto::EBIndex<CELL>                     pr_vof;
  typedef Proto::Point                             pr_pt;
  typedef Proto::Box                               pr_box;
  typedef Proto::HostBoxData<int, 1>               pr_hbd_int;
  typedef Chombo4::Box                             ch_box;
  typedef Chombo4::IntVect                         ch_iv;
  typedef Chombo4::DataIterator                    ch_dit;
  typedef Chombo4::DataIndex                       ch_datind;
  typedef Chombo4::MayDay                          ch_mayday;

  ///get the serialization size of this vector
  template<typename entry_t>
  size_t
  vector_charsize(const vector<entry_t>& a_vec)
  {
    //start with the size of the vector as an int
    size_t charsiz = sizeof(int);
    for(int ivec = 0; ivec < a_vec.size(); ivec++)
    {
      charsiz += sizeof(size_t); //size of entry
      charsiz += a_vec[ivec].charsize(); //data for entry
    }
    return charsiz;
  }

  ///output this vector to a bytestream
  template<typename entry_t>
  size_t
  vector_energize(char*                   a_outbuf,
                  const vector<entry_t> & a_vec,
                  size_t                  a_expected_size)
  {
    char*  charbuf = a_outbuf;
    size_t charsiz = 0;
    {
      int* intbuf = (int*) charbuf;
      *intbuf = a_vec.size();
      charbuf += sizeof(int);
      charsiz += sizeof(int);
    }
    for(int ivec = 0; ivec < a_vec.size(); ivec++)
    {
      const auto& entry = a_vec[ivec];
      size_t entrysize  = entry.charsize();
        
      size_t* sizbuf = (size_t*) charbuf;
      *sizbuf = entrysize;
      charbuf += sizeof(size_t);
      charsiz += sizeof(size_t);
        
        
      entry.energize(charbuf, entrysize);
      charbuf += entrysize;
      charsiz += entrysize;
    } 
    if(charsiz != a_expected_size)
    {
      ch_mayday::Error("EBCMVolu::energize: size check failed");
    }
    return charsiz;
  }

  ///define this vector from a bytestream
  template<typename entry_t>
  size_t
  vector_linearIn(char*             a_inbuf,
                  vector<entry_t> & a_vec,
                  size_t            a_expected_size)
  {
    char*  charbuf = a_inbuf;
    size_t charsiz = 0;
    int* intbuf = (int*) charbuf;
    int vecsize = *intbuf;
    charbuf += sizeof(int);
    charsiz += sizeof(int);
      
    a_vec.resize(vecsize);
    for(int ivec = 0; ivec < a_vec.size(); ivec++)
    {
      auto& entry = a_vec[ivec];
        
      size_t* sizbuf = (size_t*) charbuf;
      size_t entrysize  = *sizbuf;
      charbuf += sizeof(size_t);
      charsiz += sizeof(size_t);
        
      entry.linearIn(charbuf, entrysize);
      charbuf += entrysize;
      charsiz += entrysize;
    } 
    if(charsiz != a_expected_size)
    {
      ch_mayday::Error("EBCMVolu::energize: size check failed");
    }
    return charsiz;
  }
  
  ///
  template<CENTERING cent, int order>
  class EBCM_Face
  {
  public:

    pr_vof  m_lo;
    pr_vof  m_hi;
    IndexedMoments<DIM-1, order> m_facmom;
    void defineAsRegular(pr_vof a_lo, pr_vof a_hi, double a_dx)
    {                                                  
      m_lo = a_lo;
      m_hi = a_hi;
      m_facmom.setRegular(a_dx);
    }

    bool operator==(const EBCM_Face& a_input) const
    {
      return ((m_lo == a_input.m_lo)  && (m_hi == a_input.m_hi));
    }

    //for debugging.   proper IO would involve the moments.
    inline void print() const
    {
      int icent = cent;
      pout() << icent << "-side face:( " << m_lo << " , " << m_hi << " )" << endl;
    }

    
    size_t charsize() const
    {
      size_t charsiz = 0;
      charsiz += m_lo.charsize();
      charsiz += m_hi.charsize();
      charsiz += m_facmom.charsize();
      return charsiz;
    }
    size_t energize(char* a_outbuf, size_t a_expected_size) const
    {
      size_t charsiz = 0;
      char*  charbuf = (char*) a_outbuf;

      m_lo.energize(charbuf);
      charbuf += m_lo.charsize();
      charsiz += m_lo.charsize();

      m_hi.energize(charbuf);
      charbuf += m_hi.charsize();
      charsiz += m_hi.charsize();

      m_facmom.energize(charbuf);
      charbuf += m_facmom.charsize();
      charsiz += m_facmom.charsize();


      if(charsiz != a_expected_size)
      {
        ch_mayday::Error("EBCMFace::energize: size check error");
      }
      return charsiz;
    }

    ///
    size_t linearIn(char* a_inbuf, size_t a_expected_size) 
    {

      size_t charsiz = 0;
      char*  charbuf = (char*) a_inbuf;

      m_lo.linearIn(charbuf);
      charbuf += m_lo.charsize();
      charsiz += m_lo.charsize();

      m_hi.linearIn(charbuf);
      charbuf += m_hi.charsize();
      charsiz += m_hi.charsize();

      m_facmom.linearIn(charbuf);
      charbuf += m_facmom.charsize();
      charsiz += m_facmom.charsize();


      if(charsiz != a_expected_size)
      {
        ch_mayday::Error("EBCMFace::linearIn: size check error");
      }
      return charsiz;

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

    inline void print() const
    {
      pout() <<  "cut face:(" << m_pt << ")" << endl;
    }

    size_t charsize() const
    {
      size_t charsiz = 0;
      charsiz +=     m_pt.charsize();
      charsiz += m_ebfmom.charsize();
      return charsiz;
    }
    void energize(char* a_outbuf, size_t a_expected_size) const
    {
      size_t charsiz = 0;
      char*  charbuf = (char*) a_outbuf;

      size_t ptsize = m_pt.charsize();
      m_pt.energize(charbuf, ptsize);
      charbuf += m_pt.charsize();
      charsiz += m_pt.charsize();


      m_ebfmom.energize(charbuf);
      charbuf += m_ebfmom.charsize();
      charsiz += m_ebfmom.charsize();

      if(charsiz != a_expected_size)
      {
        ch_mayday::Error("EBCM::IrregFace::energize: size check error");
      }
      return;

    }
    void linearIn(char* a_inbuf, size_t a_expected_size) 
    {
      size_t charsiz = 0;
      char*  charbuf = (char*) a_inbuf;
      size_t ptsize = m_pt.charsize();

      m_pt.linearIn(charbuf);
      charbuf += m_pt.charsize();
      charsiz += m_pt.charsize();

      m_ebfmom.linearIn(charbuf);
      charbuf += m_ebfmom.charsize();
      charsiz += m_ebfmom.charsize();

      if(charsiz != a_expected_size)
      {
        ch_mayday::Error("EBCMFace::linearIn: size check error");
      }
      return;

    }
  };

  template<int order>
  class EBCM_Volu
  {
  public:

    ///Data:
    ///this is another one of those classes with its data just hanging out in the breeze
    typedef Chombo4::GeometryService<order> ch_geoserv;
    typedef Proto::IndexedMoments<DIM-1, order> IndMomSDMinOne;
    typedef Proto::IndexedMoments<DIM  , order> IndMomSpaceDim;
    typedef EBCM::EBCM_Face<XFACE, order> ebcm_xface;
    typedef EBCM::EBCM_Face<YFACE, order> ebcm_yface;
    typedef EBCM::EBCM_Face<ZFACE, order> ebcm_zface;
    typedef EBCM::IrregFace<       order> ebcm_iface;
    
    bool                          m_regular;
    //centroid location in global space (not relative to cell center)
    RealVect   m_centroid;

    //volume fraction
    double     m_kappa; 
    IndexedMoments<DIM, order> m_volmom;

    //used as a unique identifier.  
    //Also volume moments are centered around this cell center.  Hence dx.
    Point      m_pt; 
    double     m_dx;
    
    //
    vector<ebcm_xface>  m_xfaces;
    vector<ebcm_yface>  m_yfaces;
    vector<ebcm_zface>  m_zfaces;
    vector<ebcm_iface>  m_ifaces;

    vector<Point> m_cells; //all cells marked for this entry.   

    ///The usual set of functions involved in serialization

    ///fill a bytestream with the necessary data to reproduce this class
    void energize(char* a_outbuf, size_t a_expected_size) const
    {
      char*  charbuf = a_outbuf;
      size_t charsiz = 0;
      //flag for regular or not
      int iflag = -486;
      if(m_regular)
      {
        iflag = 1;
      }
      {
        ///flag for m_regular
        int* intbuf = (int *)charbuf;
        *intbuf = iflag;
        charbuf += sizeof(int);
        charsiz += sizeof(int);

        ///m_dx
        double* doubbuf = (double*) charbuf;
        *doubbuf = m_dx;
        charbuf += sizeof(double);
        charsiz += sizeof(double);
      }
      for(int idir = 0; idir < DIM; idir++)
      {
        //m_pt 
        int* intbuf = (int*)charbuf;
        *intbuf = m_pt[idir];
        charbuf += sizeof(int);
        charsiz += sizeof(int);
        
        //m_centroid
        double* doubbuf = (double*) charbuf;
        *doubbuf = m_centroid[idir];
        charbuf += sizeof(double);
        charsiz += sizeof(double);

      }
      
      //5 vectors
      size_t xsize = vector_charsize(m_xfaces);
      size_t ysize = vector_charsize(m_xfaces);
      size_t zsize = vector_charsize(m_xfaces);
      size_t isize = vector_charsize(m_ifaces);
      size_t csize = vector_charsize(m_cells);
      {
        size_t* sizbuf = (size_t*) charbuf;
        *sizbuf = xsize; sizbuf++;
        *sizbuf = ysize; sizbuf++;
        *sizbuf = zsize; sizbuf++;
        *sizbuf = isize; sizbuf++;
        *sizbuf = csize;
      
        charbuf += 5*sizeof(size_t);
        charsiz += 5*sizeof(size_t);
      }

      vector_energize(charbuf, m_xfaces, xsize); charbuf += xsize; charsiz += xsize;
      vector_energize(charbuf, m_yfaces, ysize); charbuf += ysize; charsiz += ysize;
      vector_energize(charbuf, m_zfaces, zsize); charbuf += zsize; charsiz += zsize;
      vector_energize(charbuf, m_ifaces, isize); charbuf += isize; charsiz += isize;
      vector_energize(charbuf, m_cells , csize); charbuf += csize; charsiz += csize;

      //
      if(charsiz != a_expected_size)
      {
        ch_mayday::Error("EBCMVolu::energize: size check failed");
      }
    }
    
    ///define this object via bytestream
    void linearIn(const char* a_inbuf, size_t a_expected_size)
    {
      char*  charbuf = (char*)a_inbuf;
      size_t charsiz = 0;
      //flag for regular or not
      int iflag = -486;
      if(m_regular)
      {
        iflag = 1;
      }
      {
        ///flag for m_regular
        int* intbuf = (int *)charbuf;
        *intbuf = iflag;
        charbuf += sizeof(int);
        charsiz += sizeof(int);

        ///m_dx
        double* doubbuf = (double*) charbuf;
        *doubbuf = m_dx;
        charbuf += sizeof(double);
        charsiz += sizeof(double);
      }
      for(int idir = 0; idir < DIM; idir++)
      {
        //m_pt 
        int* intbuf = (int*)charbuf;
        *intbuf = m_pt[idir];
        charbuf += sizeof(int);
        charsiz += sizeof(int);
        
        //m_centroid
        double* doubbuf = (double*) charbuf;
        *doubbuf = m_centroid[idir];
        charbuf += sizeof(double);
        charsiz += sizeof(double);
      }
      
      //5 vectors
      size_t xsize,ysize, zsize, isize, csize;
      {
        size_t* sizbuf = (size_t*) charbuf;
        xsize = *sizbuf; sizbuf++;
        ysize = *sizbuf; sizbuf++;
        zsize = *sizbuf; sizbuf++;
        isize = *sizbuf; sizbuf++;
        csize = *sizbuf;
      
        charbuf += 5*sizeof(size_t);
        charsiz += 5*sizeof(size_t);
      }
        
      vector_linearIn(charbuf, m_xfaces, xsize); charbuf += xsize; charsiz += xsize;
      vector_linearIn(charbuf, m_yfaces, ysize); charbuf += ysize; charsiz += ysize;
      vector_linearIn(charbuf, m_zfaces, zsize); charbuf += zsize; charsiz += zsize;
      vector_linearIn(charbuf, m_ifaces, isize); charbuf += isize; charsiz += isize;
      vector_linearIn(charbuf, m_cells , csize); charbuf += csize; charsiz += csize;

      //
      if(charsiz != a_expected_size)
      {
        ch_mayday::Error("EBCMVolu::linearIn: size check failed");
      }
      
    }
    
    ///Get the bytestream size required to hold this object
    size_t charsize() const
    {
      size_t charsiz = 0;
      //flag for regular or not
      charsiz += sizeof(int);
      ///m_dx
      charsiz += sizeof(double);

      for(int idir = 0; idir < DIM; idir++)
      {
        //m_pt , m_centroid
        charsiz += sizeof(int);
        charsiz += sizeof(double);
      }
      
      //5 vectors each need a length
      charsiz += 5*sizeof(size_t);

      size_t xsize = vector_charsize(m_xfaces);
      size_t ysize = vector_charsize(m_xfaces);
      size_t zsize = vector_charsize(m_xfaces);
      size_t isize = vector_charsize(m_ifaces);
      size_t csize = vector_charsize(m_cells);
        

      charsiz += xsize;
      charsiz += ysize;
      charsiz += zsize;
      charsiz += isize;
      charsiz += csize;

      return charsiz;
    }

    ///older name for energize (bit more gravitas, I guess)
    void linearOut(char* a_outbuf, size_t a_expected_size) const
    {
      energize(a_outbuf);
    }
    
    inline void print() const
    {
      pout() << "volume :(" << m_pt << "), ";
      pout() << "volmom = " << endl;
      m_volmom.print();
    }
    
    bool operator < (const EBCM_Volu<order>& a_input) const
    {
      return (m_pt < a_input.m_pt);
    }
    //keeping it simple---just looking for matching faces
    template<CENTERING cent>
    static bool foundMatch(const vector< EBCM_Face<cent, order>  > & a_vecOne,
                           const vector< EBCM_Face<cent, order>  > & a_vecTwo)
    {
      bool matchFound = false;
      for(int ione= 0; ione < a_vecOne.size(); ione++)
      {
        for(int itwo= 0; itwo < a_vecTwo.size(); itwo++)
        {
          if(a_vecOne[ione] == a_vecTwo[itwo])
          {
            matchFound = true;
          }
        }
      }

      return matchFound;
    }

    //see if ones low faces matches the others high faces.
    static bool connected(const EBCM_Volu<order>& a_volOne,
                          const EBCM_Volu<order>& a_volTwo)
    {
      bool connected = false;
      if(!connected)
      {        
        connected = foundMatch<XFACE>(a_volOne.m_xfaces, a_volTwo.m_xfaces);
      }
      if(!connected)
      {
        connected = foundMatch<YFACE>(a_volOne.m_yfaces, a_volTwo.m_yfaces);
      }

#if DIM==3
      if(!connected)
      {
        connected = foundMatch<ZFACE>(a_volOne.m_zfaces, a_volTwo.m_zfaces);
      }
#endif      

      return connected;
    }
    template <CENTERING cent>
    void fillFaceData(      vector< EBCM_Face< cent, order> >           & a_facedata,
                      const vector< EBIndex<   cent > >                 & a_ebfaces,
                      const HostIrregData<     cent, IndMomSDMinOne, 1> & a_facedat)
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
                           const ch_datind                       & a_dit,
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

      auto xebfac = xebfaclo;
      auto yebfac = yebfaclo;
      xebfac.insert(xebfac.end(), xebfachi.begin(), xebfachi.end());
      yebfac.insert(yebfac.end(), yebfachi.begin(), yebfachi.end());
      
      fillFaceData<XFACE>(m_xfaces, xebfac, xfacdat);
      fillFaceData<YFACE>(m_yfaces, yebfac, yfacdat);

#if DIM==3                     
      auto zebfaclo= graph.getZFaces(a_vof, Side::Lo);
      auto zebfachi= graph.getZFaces(a_vof, Side::Hi);
      
      auto zebfac = zebfaclo;
      zebfac.insert(zebfac.end(), zebfachi.begin(), zebfachi.end());
      const auto & zfacdat = (*(a_geoserv->getZFaceData(a_domain)))[a_dit];
      fillFaceData<ZFACE>(m_zfaces, zebfac, zfacdat);
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
      m_ifaces.resize(1);
      m_ifaces[0].m_ebfmom = ebfdat(ebf, 0);
      m_ifaces[0].m_pt = m_pt;
      
      defineDerivedMeasures(a_dx);
    }

    void defineDerivedMeasures(double a_dx)
    {
      double cellVol = D_TERM(a_dx, *a_dx, *a_dx);
      m_kappa = m_volmom[0]/cellVol;
      //relies upon volmom  and dx being set
      for(int idir = 0; idir < DIM; idir++)
      {
        cellVol *= a_dx;
        double centerloc = (m_pt[idir] + 0.5)*a_dx;
        double scaledCent = m_volmom[1+idir];
        m_centroid[idir]= centerloc + scaledCent;
      }
    }
    void defineAsRegular(pr_pt a_pt, double a_dx)
    {
      m_regular = true;
      m_pt = a_pt;
      m_dx = a_dx;
      m_cells.resize(1, m_pt);
      m_ifaces.resize(0);
      m_volmom.setRegular(a_dx);
      {
        m_xfaces.resize(2);
        pr_pt lopt = m_pt - Point::Basis(0);
        pr_pt hipt = m_pt + Point::Basis(0);
        pr_vof lovof(lopt, 0);
        pr_vof hivof(hipt, 0);
        pr_vof ptvof(m_pt, 0);
        m_xfaces[0].defineAsRegular(lovof, ptvof, a_dx);
        m_xfaces[1].defineAsRegular(ptvof, hivof, a_dx);
      }
      {
        m_yfaces.resize(2);
        pr_pt lopt = m_pt - Point::Basis(1);
        pr_pt hipt = m_pt + Point::Basis(1);
        pr_vof lovof(lopt, 0);
        pr_vof hivof(hipt, 0);
        pr_vof ptvof(m_pt, 0);

        m_yfaces[0].defineAsRegular(lovof, ptvof, a_dx);
        m_yfaces[1].defineAsRegular(ptvof, hivof, a_dx);
      }
#if DIM==3
      {
        m_zfaces.resize(2);
        pr_pt lopt = m_pt - Point::Basis(2);
        pr_pt hipt = m_pt + Point::Basis(2);
        pr_vof lovof(lopt, 0);
        pr_vof hivof(hipt, 0);
        pr_vof ptvof(m_pt, 0);
        m_zfaces[0].defineAsRegular(lovof, ptvof, a_dx);
        m_zfaces[1].defineAsRegular(ptvof, hivof, a_dx);
      }
#endif      
      defineDerivedMeasures(a_dx);
    }
  };// end class EBCM_volu


  ///
  /**
     Class holds a optionally-merged sparse embedded boundary description. 
  */
  template< int order>
  class EBCM_Graph
  {
  public:
    /**
       First a bunch of static functions that will probably be useful in other contexts
       They are not standalone because I want to use templated typedefs in declarations.
    **/
    ///Speaking of which,
    typedef Chombo4::GeometryService<order> ch_geoserv;
    typedef EBCM_Volu<order>                ebcm_volu;
  
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



  
    ///code that fills EBCM_Volu objects
    /**
       Face objects contain moment informaiton
       Data that needs to be filled:

       vector<ebcm_xface >  m_xfaceslo;
       vector<EBCM_Face<YFACE> >  m_yfaceslo;
       vector<EBCM_Face<ZFACE> >  m_zfaceslo;
       vector<ebcm_xface   >  m_xfaceshi;
       vector<EBCM_Face<YFACE> >  m_yfaceshi;
       vector<EBCM_Face<ZFACE> >  m_zfaceshi;
       vector<IrregFace<order> > m_ifaces;

       //only true if it started as a regular cell and never got merged
       bool                          m_regular;

       //centroid location in global space (not relative to cell center)
       RealVect m_centroid;
       double     m_kappa; //volume fraction
       IndexedMoments<DIM, order> m_volmom;

       //used as a unique identifier.  
       //Also volume moments are centered around this cell center.  Hence dx.
       Point    m_pt; 
       double     m_dx
       vector<Point> m_cells; //all cells marked for this entry.   
    */


    ///  get how far to shift moments when merging two volumes
    static
    inline void
    getShifts(IndexTM<double, DIM>                 & a_shiftOne,
              IndexTM<double, DIM>                 & a_shiftTwo,
              const ebcm_volu                      & a_volMer,
              const ebcm_volu                      & a_volOne,
              const ebcm_volu                      & a_volTwo)
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
    shiftAndMergeVolumeData(      ebcm_volu                      & a_volMer,
                            const ebcm_volu                      & a_volOne,
                            const ebcm_volu                      & a_volTwo)
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
  

    ///
    static inline void
    shiftAndMergeIrregData(      ebcm_volu                      & a_volMer,
                           const ebcm_volu                      & a_volOne,
                           const ebcm_volu                      & a_volTwo)
    {
      ///This does the same thing as mergeFaces (just slap the two vectors together into the third).
      ///The moments are all wrong, though, and need to be shifted
      a_volMer.m_ifaces = a_volOne.m_ifaces;
      a_volMer.m_ifaces.insert(a_volMer.m_ifaces.end(), a_volTwo.m_ifaces.begin(), a_volTwo.m_ifaces.end());
      IndexTM<double, DIM> shiftOne, shiftTwo;

      int curmer = 0;
      for(int iface = 0; iface < a_volOne.m_ifaces.size(); iface++)
      {
        a_volMer.m_ifaces[curmer].m_ebfmom.shift(shiftOne);
        curmer++;
      }
      for(int iface = 0; iface < a_volTwo.m_ifaces.size(); iface++)
      {
        a_volMer.m_ifaces[curmer].m_ebfmom.shift(shiftTwo);
        curmer++;
      }
    }


    //the vectors should be very short 
    template <CENTERING cent>
    static inline void
    addAllOfFirstNotInSecondToThird(vector<      EBCM_Face<cent, order> >  & a_merged,
                                    const vector<EBCM_Face<cent, order> >  & a_vecOne,
                                    const vector<EBCM_Face<cent, order> >  & a_vecTwo)
    {
      for(int ione= 0; ione < a_vecOne.size(); ione++)
      {
        const auto& oneval = a_vecOne[ione];

        bool matchFound = false;
        for(int itwo= 0; itwo < a_vecTwo.size(); itwo++)
        {
          const auto& twoval = a_vecTwo[itwo];
          if(oneval == twoval) //the old one -two
          {
            matchFound = true;
          }
        }
        if(!matchFound)
        {
          a_merged.push_back(oneval);
        }
      }
    }

      
    ///
    /**
       makes the face centered structures of merged volume volMer from face structures
       in volOne and volTwo.   If we detect matching faces, that is the line over which 
       we are merging and those faces go away.  All other faces get thrown into the mix.
    */
    template <CENTERING cent>
    static inline void
    consolidateFaces(vector<      EBCM_Face<cent, order> >  & a_merged,
                     const vector<EBCM_Face<cent, order> >  & a_vecOne,
                     const vector<EBCM_Face<cent, order> >  & a_vecTwo,
                     bool a_printStuff)
    {
      a_merged.resize(0);
      //add all faces from vecOne that have no match in vecTwo
      addAllOfFirstNotInSecondToThird(a_merged, a_vecOne, a_vecTwo);
      
      //add all faces from vecTwo that have no match in vecOne
      addAllOfFirstNotInSecondToThird(a_merged, a_vecTwo, a_vecOne);

      if(a_printStuff)
      {
        int icent = cent;
        pout() << icent << "-side face sizes" 
               << ": merged = "  << a_merged.size() 
               << ", vecone = "  << a_vecOne.size() 
               << ", vectwo = "  << a_vecTwo.size() << endl;
      }
    }
    
    ///Consolidate in each direction    static inline void
    static inline void
    consolidateFacesAllDirections(      ebcm_volu        & a_volMer,
                                  const ebcm_volu        & a_volOne,
                                  const ebcm_volu        & a_volTwo,
                                  bool a_printStuff)
    {

      consolidateFaces(a_volMer.m_xfaces, a_volOne.m_xfaces, a_volTwo.m_xfaces,a_printStuff); 
      consolidateFaces(a_volMer.m_yfaces, a_volOne.m_yfaces, a_volTwo.m_yfaces,a_printStuff);
#if DIM==3                                                                    
      consolidateFaces(a_volMer.m_zfaces, a_volOne.m_zfaces, a_volTwo.m_zfaces,a_printStuff);
#endif      
    }


    ///make EBCM_Volu volMer that has eaten both volOne and volTwo
    static inline int
    mergeTwoVolumesIntoThird(      ebcm_volu                      & a_volMer,
                             const ebcm_volu                      & a_volOne,
                             const ebcm_volu                      & a_volTwo,
                             double a_dx, bool a_printStuff)
    {
      a_volMer.m_regular = false;
      a_volMer.m_pt = a_volTwo.m_pt; //just having   the bigger one have the point
      a_volMer.m_cells = a_volOne.m_cells;
      a_volMer.m_cells.insert(a_volMer.m_cells.end(), a_volTwo.m_cells.begin(), a_volTwo.m_cells.end());
      shiftAndMergeVolumeData(a_volMer, a_volOne, a_volTwo);
      shiftAndMergeIrregData( a_volMer, a_volOne, a_volTwo);

      consolidateFacesAllDirections( a_volMer, a_volOne, a_volTwo, a_printStuff);
      
      a_volMer.defineDerivedMeasures(a_dx);
      return 0;
    } //end funciton mergeTwoIntoThird

    static inline bool
    volumesConnected(const ebcm_volu      & a_volOne,
                     const ebcm_volu      & a_volTwo)
    {
      return EBCM_Volu<order>::connected(a_volOne, a_volTwo); // 
    }
    /// 
    /**
       Recursive function that merges a vector of volumes into a single volume.
       Fails with MayDay::Error if the vector is not has disconnected components.
    **/ 
    static
    inline ebcm_volu      
    mergeManyIntoOne(const vector<ebcm_volu>& a_many,
                     double a_dx, bool a_printStuff)
    {
      ebcm_volu retval;
      if(a_many.size() == 1)
      {
        retval = a_many[0];
      }
      else
      {
        //find a vof connected to first vof
        int vofind0 = 0;
        const auto & volu0 = a_many[0];
        int vofind1 = -1;
        bool found = false;
        //start at 1 because vofind0= 0
        for(int ivec = 1; (ivec < a_many.size()) && (!found); ivec++)
        {
          if(volumesConnected(volu0, a_many[ivec]))
          {
            vofind1 = ivec;
            found = true;
          }
        }
        if(!found)
        {
          Chombo4::MayDay::Error("mergeManyIntoOne error: could not find connected volume");
        }
        ebcm_volu   third; 
        int errflag =
          mergeTwoVolumesIntoThird(third, a_many[vofind0], a_many[vofind1],
                                   a_dx, a_printStuff);
        if(errflag != 0)
        {
          MayDay::Error("mergeManyIntoOne error: function mergeTwoVolumesIntoOne failed");
        }

        //we reduce the size of the vector by one with every call.
        vector<ebcm_volu> volsLeft(a_many.size() - 1);
        volsLeft[0] = third;
        int ileft = 1;
        for(int imany = 0; imany < a_many.size(); imany++)
        {
          if((imany != vofind0) && (imany != vofind1))
          {
            volsLeft[ileft] = a_many[imany];
            ileft++;
          }
        }
        //the rare instance where recusion is just way easier.
        retval = mergeManyIntoOne(volsLeft, a_dx, a_printStuff);
      
      }
      return retval;
    } //end function mergeManyIntoOne
    

    ///
    /**
       This coarsens, then refines the input cell to get a 2x2x2
       (at least in cells, some may be covered cells) region
       merge it all together.    This algorithm should be independent 
       of box division decisions because of the blocking factor thing.
       Absolutely not guaranteed to get kappa >= 1.
       --dtg
    **/
    static
    inline ebcm_volu
    getRefinedCoarseVolume(vector<pr_vof>                        & a_vofsEaten,
                           const  pr_vof                         & a_startvof,
                           const shared_ptr<ch_geoserv>          & a_geoserv,
                           const ch_box                          & a_domain,
                           const ch_datind                       & a_dit,
                           const double                          & a_dx,
                           const shared_ptr<set<pr_vof   > >     & a_vofsLeft,
                           const pr_ebgrap                       & a_graph,
                           bool a_printStuff)
    {
      ebcm_volu retval;
      pr_pt finept = a_startvof.m_pt;
      pr_pt coarpt = finept.coarsen(2);
      pr_box coarbox(coarpt, coarpt);
      pr_box finebox = coarbox.refine(2);
      pr_box grownbox = finebox;
      grownbox &= a_graph.getRegion();
      a_vofsEaten.resize(0);
      if(a_printStuff)
      {
        pout() << "grownbox = " << grownbox << endl;
      }
      for(auto bit = grownbox.begin(); bit != grownbox.end(); ++bit)
      {
        auto pt = *bit;
        vector<pr_vof> vofs = a_graph.getVoFs(pt);
        for(int ivof = 0; ivof < vofs.size(); ivof++)
        {
          //have to be certain to not double count here.
          //vofsLeft gets updated after createEBCM_Volume returns
          bool stillhere = (a_vofsLeft->find(vofs[ivof]) != a_vofsLeft->end());
          if(stillhere)
          {
            a_vofsEaten.push_back(vofs[ivof]);
            if(a_printStuff)
            {
              pout() << "adding   " << vofs[ivof].m_pt << endl;
            }
          }
          else
          {
            if(a_printStuff)
            {
              pout() << "skipping " << vofs[ivof].m_pt << endl;
            }
          }
          
        }
      } // end loop over cells in grown box      

      if(a_printStuff)
      {
        pout() << "getGhostedVolume: a_vofsEaten.size() == " << a_vofsEaten.size() << endl;
      }
      vector<ebcm_volu > allvols(a_vofsEaten.size());
      for(int ivol = 0; ivol <         a_vofsEaten.size(); ivol++)
      {
        auto curpt = a_vofsEaten[ivol].m_pt;
        if(a_graph.isCovered(curpt))
        {
          Chombo4::MayDay::Error("cannot create EBCM_Volume from a covered cell");
        }
        else if(a_graph.isRegular(curpt))
        {
          allvols[ivol].defineAsRegular(curpt, a_dx);
        }
        else
        {
          allvols[ivol].defineAsIrregular(a_vofsEaten[ivol], a_geoserv, a_domain, a_dit, a_dx);
        }

      } // end loop over input vofs in the ghosted region

      //make two into one until there are no more
      retval = mergeManyIntoOne(allvols, a_dx, a_printStuff);

      return retval;
    }
    
    
    ///
    /**
       Creates volume from EB description.   If merger is necessary and called for,
       it calls for merger.   VoFsLeft gets modified by calling function.
    */
    static
    inline ebcm_volu
    createEBCM_Volume(vector<pr_vof>                        & a_vofsEaten,
                      const pr_vof                          & a_startvof,
                      const shared_ptr<ch_geoserv>          & a_geoserv,
                      const ch_box                          & a_domain,
                      const ch_datind                       & a_dit,
                      const double                          & a_dx,
                      const shared_ptr<set<pr_vof   > >     & a_vofsLeft,
                      const pr_ebgrap                       & a_graph,
                      bool                                    a_mergeSmallCells,
                      bool                                    a_printStuff)
    {
      ebcm_volu startvol, retval;
      if(a_graph.isCovered(a_startvof.m_pt))
      {
        Chombo4::MayDay::Error("cannot create EBCM_Volume from a covered cell");
      }
      else if(a_graph.isRegular(a_startvof.m_pt))
      {
        startvol.defineAsRegular(a_startvof.m_pt, a_dx);
      }
      else
      {
        startvol.defineAsIrregular(a_startvof, a_geoserv, a_domain, a_dit, a_dx);
      }
      
      if( (startvol.m_kappa == 1.) || (!a_mergeSmallCells))
      {
        retval = startvol;
        a_vofsEaten.resize(1, a_startvof);
      }
      else
      {
        a_vofsEaten.resize(0);

        retval =
          getRefinedCoarseVolume(a_vofsEaten,  a_startvof, a_geoserv, a_domain,
                                 a_dit, a_dx, a_vofsLeft, a_graph, a_printStuff);
      }
      return retval;
    }

    ///
    static inline bool
    getRegularNeighborPoint(pr_pt                            & a_regpt,
                            const ebcm_volu                  & a_smallvol,
                            const pr_hbd_int                 & a_map)
    {
      const pr_pt& startpt = a_smallvol.m_pt;
      bool found = false;
      pr_box graphregion= a_map.box();
      for(SideIterator sit; sit.ok(); ++sit)
      {
        int isign = sign(sit());
        for(int idir = 0; ((idir < DIM) && !found); ++idir)
        {
          pr_pt offset = isign*pr_pt::Basis(idir);
          auto newpt = startpt + offset;
          if(graphregion.contains(newpt))
          {
            if(a_map(newpt, 0) == -1)
            {
              a_regpt= newpt;
              found = true;
            }
          }
        }
      }
      return found;
    }


    ///put indicies of volumes in correct places in map
    static inline void
    fixVolumeEntriesInMap(shared_ptr<pr_hbd_int >                      & a_map,
                          const shared_ptr<vector<ebcm_volu> >         & a_volumes)
    {
      for(int ivol = 0; ivol < a_volumes->size(); ivol++)
      {
        const auto& curvol = (*a_volumes)[ivol];
        //these are the points that went into the volume
        const auto& volpts = curvol.m_cells;
        for(int icell = 0; icell < volpts.size(); icell++)
        {
          const pr_pt& curpt = volpts[icell];
          (*a_map)(curpt, 0) = ivol;
        }
      }
    }
    ///
    /**
       Modifies a_volumes in place by merging all of the entries in a_tooSmall
       to their neighbors found in a_volumes.   a_volumes should not change in size
       but some of its entries will get fatter.
    */
    static inline void
    mergeSmallVolumesWithNeighbors(shared_ptr<vector<ebcm_volu > >           & a_volumes,
                                   shared_ptr<pr_hbd_int>                    & a_map,
                                   shared_ptr<vector<ebcm_volu > >           & a_tooSmall,
                                   const shared_ptr<ch_geoserv>              & a_geoserv,
                                   const ch_box                              & a_domain,
                                   const ch_datind                           & a_dit,
                                   const double                              & a_dx,
                                   const shared_ptr<set<pr_vof   > >         & a_vofsLeft,
                                   const pr_ebgrap                           & a_graph,
                                   bool a_printStuff)
    {
      if(a_tooSmall->size() > 0)
      {
        //remember input volumes
        vector<ebcm_volu > inputVols = *a_volumes;
        a_volumes->resize(0);

        /**
           Remember which members of inputVols get used in merger
           so we do not double add them to a_volumes
        **/
        set<int> alreadyUsed;

        for(int ismall = 0; ismall < a_tooSmall->size(); ismall++)
        {
          const auto& smallvol = (*a_tooSmall)[ismall];
          ebcm_volu merged, bigvol;
          Point regPt;
          bool hasRegularNeighbor = getRegularNeighborPoint(regPt, smallvol,  *a_map);
          if(hasRegularNeighbor)
          {
            bigvol.defineAsRegular(regPt, a_dx);
            //set map to fake non-negative number so it will no longer get tagged as regular
            //map gets fixed up at the end of this function.
            (*a_map)(regPt, 0) = 4586;
          }
          else
          {
            int ineigh = -1; bool found = false;
            //find a neighbor in the input volumes (which are all big)
            for(int ivol = 0; ((ivol < inputVols.size()) && (!found)); ivol++)
            {
              const auto& bigvol = inputVols[ivol];
              if(volumesConnected(smallvol, bigvol))
              {
                ineigh = ivol;
                found = true;
              }
            }
            if(!found)
            {
              pout() << "mergeSmallVolumesWithNeighbors ismall = " << ismall << endl;
              MayDay::Error("small volume has no connections");
            }
            bigvol = inputVols[ineigh];
            alreadyUsed.insert(ineigh);
          }

          mergeTwoVolumesIntoThird(merged, smallvol, bigvol, a_dx, a_printStuff);
          int curIndex = a_volumes->size();
          //fix map to match merger
          for(int icell = 0; icell < merged.m_cells.size(); icell++)
          {
            auto map_pt=  merged.m_cells[icell];
            (*a_map)(map_pt, 0) = curIndex;
          }
          a_volumes->push_back(merged);
        }

        /**
           Add into a_volumes any inputs that were not used in merging.
           The merged ones are already in a_volumes.
        **/
        for(int ivol = 0; ivol < inputVols.size(); ivol++)
        {
          if(alreadyUsed.find(ivol) == alreadyUsed.end())
          {
            a_volumes->push_back(inputVols[ivol]);
          }
        }
      }

      fixVolumeEntriesInMap(a_map, a_volumes);
    }

    ///
    /**
       Translate all this mixed stuff into sparse  representation that we actually store
    */
    static inline  int
    makeTheBlobs(shared_ptr<vector<ebcm_volu > >     & a_volumes,
                 shared_ptr<pr_hbd_int >                   & a_map,
                 shared_ptr<set<pr_vof   > >               & a_inputVoFs, 
                 const shared_ptr<ch_geoserv>              & a_geoserv,
                 const ch_box                              & a_domain,
                 const ch_datind                           & a_dit,
                 const double                              & a_dx,
                 const pr_ebgrap                           & a_graph,
                 bool                                        a_mergeSmallCells,
                 bool a_printStuff)
    {
      bool  printStuff = false;
      a_volumes = shared_ptr<vector<ebcm_volu > >(new vector<ebcm_volu > ());
      a_map    =  shared_ptr<pr_hbd_int>(new pr_hbd_int());
      //define the map and set it all to covered.   the input set contains only regular and irregular
      a_map->define(a_graph.getRegion());
      a_map->setVal(-2);
      shared_ptr<set<pr_vof   > >  vofsLeft(new set<pr_vof>() );
      *vofsLeft = *a_inputVoFs;
      shared_ptr< vector<ebcm_volu > > tooSmalls(new vector<ebcm_volu >());
      ///two pass merger algorithm.   merge on coarsened cell then try to find neighbors
      //for the orphans that did not have enough volume inthe 
      double kappamin = 1.0;
      for(auto vofit = a_inputVoFs->begin(); vofit != a_inputVoFs->end(); ++vofit)
      {
        auto vof = *vofit;
        if(a_graph.isRegular(vof.m_pt))
        {
          //only reset it if it is set to covered.   It can also be reset to the vector value if merged
          //and we do not want to undo the merger work of nearby volumes.
          int mapval = (*a_map)(vof.m_pt, 0);
          if(mapval == -2)
          {
            //again, can be reset in the irreg section if nearby volumes merge
            (*a_map)(vof.m_pt, 0) = -1;
          }
        } //end if (this point regular          
        else if(a_graph.isIrregular(vof.m_pt))
        {
          if(vofsLeft->find(vof) != vofsLeft->end())
          {
            //begin debug
            printStuff = ((vof.m_pt[0]==25) && (vof.m_pt[1]== 27));
            //end debug
            vector<pr_vof> mergevec;
            ebcm_volu volume =
              createEBCM_Volume(mergevec, vof, a_geoserv, a_domain, a_dit, a_dx,
                                vofsLeft, a_graph, a_mergeSmallCells, printStuff);
            
      
            for(int ivec = 0; ivec < mergevec.size(); ivec++)
            {
              vofsLeft->erase(mergevec[ivec]);
            }
            kappamin = std::min(kappamin, volume.m_kappa);
            if( (volume.m_kappa < 1.) && (a_mergeSmallCells))
            {
              tooSmalls->push_back(volume);
            }
            else
            {
              a_volumes->push_back(volume);
            }
          } //end if this is in vofsleft
        } // end if this cell irregular
      } // end oop over a_ebvolumes

      fixVolumeEntriesInMap(a_map, a_volumes);

      if(a_printStuff)
      {
        pout() << "makeTheBlobs: number of small cells = " << tooSmalls->size() <<", min kappa = " << kappamin << endl;
      }
      
      if(tooSmalls->size() > 0)
      {
        printStuff= false;
        mergeSmallVolumesWithNeighbors(a_volumes, a_map,  tooSmalls,
                                       a_geoserv, a_domain, a_dit, a_dx,
                                       vofsLeft, a_graph, printStuff);
      }
      return 0;
    }//end  function makeTheBlobs

    //strong construction would be better but way less convenient
    EBCM_Graph()
    {
      m_defined = false;
    }

    inline void 
    define(const shared_ptr<ch_geoserv>          & a_geoserv,
           const ch_box                          & a_domain,
           const ch_datind                       & a_dit,
           const double                          & a_dx,
           int                                     a_ghost,
           bool                                    a_mergeSmallCells,
           bool                                    a_printStuff) 
    {
      m_defined = true;
      m_dx         = a_dx;
      m_domain     = a_domain;
      auto grids = a_geoserv->getDBL(a_domain);
      
      m_grid = grids[a_dit];
      m_grid.grow(a_ghost);
      m_grid &= a_domain;
      
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
        
        //translate all this mixed stuff into sparse  representation that we actually store
        makeTheBlobs(m_volumes, m_map, allEBInd, a_geoserv, a_domain, a_dit, a_dx, ebgrap,
                     a_mergeSmallCells, a_printStuff);
      }
    }

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

    ///
    bool hasIrregular() const
    {
      bool     regularOrCovered = m_allRegular  || m_allCovered;
      return (!regularOrCovered);
    }

    ///serialization routines
    

    ///Fills byte stream from with all the data from this object.
    inline void
    energize(char*                a_buf,
             size_t               a_expected_size,
             bool                 a_printStuff = false) const
    {
      size_t charsiz = 0;
      char* charbuf = (char *) a_buf;
      
      if(!m_defined)
      {
        ch_mayday::Error("EBCM_Graph::energize called when object is undefined");
      }

      //output a magic tag to show if all regular, all covered or has irregular
      int label = 4586;
      if(m_allRegular)
      {
        label = -1;
      }
      else if(m_allCovered)
      {
        label = -2;
      }
      // tag for all reg, all covered, has irreg
      {
        int* intbuf = (int *)(charbuf);
        *intbuf = label;
        charbuf += sizeof(int);
        charsiz += sizeof(int);
      }
      //m_dx
      {
        double* doubbuf = (double*)(charbuf);
        *doubbuf = m_dx;
        charbuf += sizeof(double);
        charsiz += sizeof(double);
      }
      
      // m_grid
      {
        m_grid.linearOut(charbuf);
        charbuf += m_grid.charsize();
        charsiz += m_grid.charsize();
        
        // m_domain
        m_domain.linearOut(charbuf);
        charbuf += m_domain.charsize();
        charsiz += m_domain.charsize();
      }
      
      /**
         Irregular volume information.  The stuff below are  only defined if 
         hasIrregular is true (which means allRegular and allCovered are false).
         or m_allCovered are true.  
         shared_ptr< vector<ebcm_volu > > m_volumes;
         shared_ptr< HostBoxData<int, 1> >      m_map;
      **/
      if( hasIrregular() )
      {
        size_t volsize = vector_charsize(*m_volumes);
        {
          size_t* sizbuf = (size_t*)charbuf;
          *sizbuf = volsize;
          charbuf += sizeof(size_t);
          charsiz += sizeof(size_t);
        }
        
        vector_energize(charbuf, *m_volumes, volsize);
        charbuf += volsize;
        charsiz += volsize;
        
        ch_box fabbox  = m_map->box();
        fabbox.energize(charbuf);
        charbuf += fabbox.charsize();
        charsiz += fabbox.charsize();
        
        size_t fabsize = m_map->charsize(fabbox, 0, 1, a_printStuff);;
        {
          size_t* sizbuf = (size_t*)charbuf;
          *sizbuf = fabsize; 
          charbuf += sizeof(size_t);
          charsiz += sizeof(size_t);
        }
        
        m_map->linearOut( charbuf, fabsize, fabbox, 0, 1, a_printStuff);;
        charbuf += fabsize;
        charsiz += fabsize;
      }
        
      if(a_expected_size != charsiz)
      {
        ch_mayday::Error("EBCMGraph::energize: size mismatch");
      }
    }

    ///Returns the size of the byte stream
    inline
    std::size_t
    charsize(const ch_box&   a_bx,
             unsigned int a_startcomp,
             unsigned int a_numcomps,
             bool a_printStuff = false) const
    {
      using ProtoSer::linearCharBufSize;
      size_t charsiz = 0;
      
      if(!m_defined)
      {
        ch_mayday::Error("EBCM_Graph::energize called when object is undefined");
      }

      //output a magic tag to show if all regular, all covered or has irregular
      {
        charsiz += sizeof(int);
      }
      //m_dx
      {
        charsiz += sizeof(double);
      }
      
      {
        // m_grid
        charsiz += m_grid.charsize();
        
        // m_domain
        charsiz += m_domain.charsize();
      }
      
      /**
         Irregular volume information.  The stuff below are  only defined if 
         hasIrregular is true (which means allRegular and allCovered are false).
         or m_allCovered are true.  
         shared_ptr< vector<ebcm_volu > > m_volumes;
         shared_ptr< HostBoxData<int, 1> >      m_map;
      **/
      if( hasIrregular() )
      {
        //byte length of volume vector
        charsiz += sizeof(size_t);
        
        size_t volsize = vector_charsize(*m_volumes);
        charsiz += volsize;
        
        ch_box fabbox  = m_map->box();
        charsiz += fabbox.charsize();

        //space  for fabsize
        charsiz += sizeof(size_t);
        
        size_t fabsize = m_map->charsize(fabbox, 0, 1, a_printStuff);
        //space for data
        charsiz += fabsize;
      }
      return charsiz;
    }

    ///blame it all on MPI
    int nComp() const
    {
      return 1;
    }
    
    ///
    /**
       Defines this object from a byte stream.  This has to be
       symmetric with the energize function.  If I could figure out
       how to force this, I would.
    */
    EBCM_Graph(void * a_buf,
               size_t a_expectedSize,
               bool a_printStuff = false)
    {

      m_defined = true;

      
      size_t charsiz = 0;
      char* charbuf = (char *) a_buf;
      
      //input a magic tag to show if all regular, all covered or has irregular
      {
        int* intbuf = (int *)(charbuf);
        int label = *intbuf;
        charbuf += sizeof(int);
        charsiz += sizeof(int);

        if(label == -1)
        {
          m_allRegular = true;
          m_allCovered = false;
        }
        else if(label == -2)
        {
          m_allRegular = false;
          m_allCovered = true;
        }
        else if(label == 4586)
        {
          //this one has cut cells
          m_allRegular = false;
          m_allCovered = false;
        }
        else
        {
          ch_mayday::Error("EBCMGraph bytstream constructor: bogus regirregcov flag");
        }
      }
      //m_dx
      {
        double* doubbuf = (double*)(charbuf);
        m_dx = *doubbuf;
        charbuf += sizeof(double);
        charsiz += sizeof(double);
      }
      
      // m_grid
      m_grid.linearIn(charbuf);
      charbuf += m_grid.charsize();
      charsiz += m_grid.charsize();
        
      // m_domain
      m_domain.linearIn(charbuf);
      charbuf += m_domain.charsize();
      charsiz += m_domain.charsize();

      
      /**
         Irregular volume information.  The stuff below are  only defined if 
         hasIrregular is true (which means allRegular and allCovered are false).
         or m_allCovered are true.  
         shared_ptr< vector<ebcm_volu > > m_volumes;
         shared_ptr< HostBoxData<int, 1> >      m_map;

         The function isIrregular ought to work because we have set m_allRegular and m_allCovered.
         Do not move this section of code to before that is true.
      **/
      if( hasIrregular() )
      {
        size_t volsize;
        {
          size_t* sizbuf = (size_t*)charbuf;
          volsize = *sizbuf;
          charbuf += sizeof(size_t);
          charsiz += sizeof(size_t);
        }
        
        vector_linearIn(charbuf, *m_volumes, volsize);
        charbuf += volsize;
        charsiz += volsize;
        
        ch_box fabbox  = m_map->box();
        fabbox.linearIn(charbuf);
        charbuf += fabbox.charsize();
        charsiz += fabbox.charsize();
        
        m_map = shared_ptr< HostBoxData<int, 1> >
          (new HostBoxData<int, 1>(fabbox));

        size_t fabsize;
        {
          size_t* sizbuf = (size_t*)charbuf;
          fabsize = *sizbuf;
          charbuf += sizeof(size_t);
          charsiz += sizeof(size_t);
        }
        m_map->linearIn( charbuf, fabsize, fabbox, 0, 1, a_printStuff);
        charbuf += fabsize;
        charsiz += fabsize;
      }
    
      if(a_expectedSize != charsiz)
      {
        ch_mayday::Error("EBCMGraph::EBCMGraph(char*): size mismatch");
      }
    }

    ///
    void
    copy(const EBCM_Graph              & a_src,
         const ch_box                  & a_srcBox,
         unsigned int                 a_srcComp,
         const ch_box                  & a_destBox,
         unsigned int                 a_destComp,
         unsigned int                 a_numcomp)
    {
      copy(a_src, a_srcBox);
    }
    
    ///this is for MPI communication
    inline void
    linearIn(char*                  a_buf,
             size_t                 a_expectedSize,
             const ch_box    &      a_bx,
             unsigned int           a_startcomp,
             unsigned int           a_numcomps,
             bool a_printStuff = false)
    {
      if(a_printStuff)
      {
        Chombo4::pout() << "EBCM_Graph::linearIn: entering " << endl;
      }

      EBCM_Graph srcGraph(a_buf, a_expectedSize, a_printStuff);
      copy(srcGraph, a_bx);

      if(a_printStuff)
      {
        Chombo4::pout() << "EBCM_Graph::linearIn: leaving " << endl;
      }
    }
    
    ///this is for MPI communication
    inline void
    linearOut(char*               a_buf,
              size_t              a_expectedSize,
              const ch_box    &   a_bx,
              unsigned int        a_startcomp,
              unsigned int        a_numcomps,
              bool a_printStuff = false) const
    {
      energize(a_buf, a_expectedSize, a_printStuff);
    }

      
    /// for exchange.   checks over box and adds volumes as needed from source.
    inline void
    copy(const EBCM_Graph              & a_src,
         const ch_box                  & a_region)
    {
      if((!m_defined) || (!a_src.m_defined))
      {
        ch_mayday::Error("EBCM::Graph::copy called using undefined objects");
      }
      if((m_allRegular && a_src.m_allRegular) || (m_allCovered && a_src.m_allCovered))
      {
        return;
      }
      else if((m_allRegular && a_src.m_allCovered) || (m_allCovered && a_src.m_allRegular))
      {
        Chombo4::pout() << "EBCM::Graph::copy: how is there an all regular graph trying to talk to an all covered one?";
        Chombo4::pout() << "Sanity demands there should be intervening cut cells." << endl;
        ch_mayday::Error();
      }
      else if(m_allRegular || m_allCovered)
      {
        //unsparsify recasts this object
        //so we can add irregular cells
        unsparsify();
      }
      pr_box region = ProtoCh::getProtoBox(a_region);
      //begin debug
//      pr_pt debpt; debpt[0] = 17; debpt[1]=1;
//      bool gettinWarma = region.contains(debpt);
//      if(gettinWarma)
//      {
//        Chombo4::pout() << "getting warmer" << endl;
//      }
      //end debug

      std::set<pr_pt> alreadyCopied; //have to do this because of merging
      for(auto bit = region.begin(); bit != region.end(); ++bit)
      {
        auto pt = *bit;
        Chombo4::pout() << "pt = " << pt << endl;
//      //begin debug
//        if(pt== debpt)
//        {
//          Chombo4::pout() << "pretty close" << endl;
//        }
//      //end debug  
        bool notFound = (alreadyCopied.find(pt) == alreadyCopied.end());
        if(notFound) //only go in if we have not already been there
        {
          alreadyCopied.insert(pt);
          if(a_src.isRegular(pt))
          {
            (*m_map)(pt, 0) = -1;

          }
          else if(a_src.isCovered(pt))
          {
            (*m_map)(pt, 0) = -2;
          }
          else
          {
            ebcm_volu  inputVol = a_src.getVolumeCoveringPoint(pt);

            int myindex = m_volumes->size();
            (*m_map)(pt, 0) =  myindex;
            m_volumes->push_back(inputVol);
            for(int ivec = 0; ivec < inputVol.m_cells.size(); ivec++)
            {
              auto inputpt = inputVol.m_cells[ivec];
              (*m_map)(inputpt, 0) = myindex;
              //begin debug
//              if(inputpt == debpt)
//              {
//                Chombo4::pout() << "red hot?"  << endl;
//              }
              ///end debug
              alreadyCopied.insert(inputpt);

            }
          }
        }
      }
        
    }

    ///
    /**
       The function unsparsify recasts this object as a fully mature 
       graph so  we can add irregular cells when filling in ghost cell
       informaiton.  This happens in the copy function just as it does
       in standard EB.
    **/
    inline void
    unsparsify()
    {
      if(!m_defined)
      {
        ch_mayday::Error("EBCM::Graph::unsparsify: called on undefined object");
      }
      if((!m_allRegular) && (!m_allCovered))
      {
        return;  //already not sparse
      }
      else
      {
        int mapval = 4586;
        if(m_allRegular)
        {
          mapval = -1;
        }
        else if(m_allCovered)
        {
          mapval = -2;
        }
        else
        {
          //I do not see how I can get here but I have been wrong before.
          ch_mayday::Error("EBCM::Graph::unsparsify: logic error");
        }
        m_allRegular = false;
        m_allCovered = false;
        m_volumes = shared_ptr< vector<ebcm_volu > >(new vector<ebcm_volu >());
        m_map     = shared_ptr< pr_hbd_int        > (new pr_hbd_int(m_grid));
        m_map->setVal(mapval);
      }
    }
    
    ///Class Data
    /**
       I am not bothering with fancy accessor functions.   Member data is left in the clear.  
       Use it wisely.
       Remember the map and the vector are not defined if allRegular or allCovered are true.   
       This should save a ton of memory where there is not much EB.
       Assignment ought to work but it will be shallow 
       (std::shared_ptr semantics).
    **/

    ///Whether define function has been called on this object. If m_defined == false, everything else is not defined.
    bool m_defined;

    ///
    bool m_allRegular;
    bool m_allCovered;
    //need  dx because we compute regular moment info on the fly.
    double                 m_dx;
    ch_box                 m_domain;
    ch_box                 m_grid;

    ///Irregular volume information.   
    /**
       We do not store unmerged regular cell volumes.
       The stuff below is not defined if(m_allRegular or m_allCovered are true
       (if(m_allRegular || m_allCovered)).
    **/

    ///This is the meta data and moment information.
    shared_ptr< vector<ebcm_volu > > m_volumes;
    //-1 regular, -2 covered, >= 0 is the index into m_volumes.  merged cells both hold the same index.
    shared_ptr< HostBoxData<int, 1> >      m_map;

    ///For stencil stuff, you should really use SubVolumeVector.   This is meant to be just for diagonostics.
    ebcm_volu getVolumeCoveringPoint(const pr_pt& a_pt) const
    {
      ebcm_volu retval;
      if(m_allCovered)
      {
        ch_mayday::Error("EBCMGraph:: covered points do not have volumes");
      }
      else if(m_allRegular)
      {
        retval.defineAsRegular(a_pt, m_dx);
      }
      else
      {
        int mapval = (*m_map)(a_pt, 0);
        if(mapval == -2)
        {
          ch_mayday::Error("EBCMGraph:: covered points do not have volumes");
        }
        else if(mapval == -1)
        {
          retval.defineAsRegular(a_pt, m_dx);
        }
        else
        {
          retval = (*m_volumes)[mapval];
        }
      }
      return retval;
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
    typedef Chombo4::NullConstructorDataFactory<   EBCM_Graph<order> >      cmgraph_factory_t;    
    typedef CH4_Data_Choreography::DistributedData<EBCM_Graph<order> >      cmgraph_distrib_t;
    
    MetaDataLevel(const shared_ptr<ch_geoserv>& a_geoserv,
                  const ch_box                & a_domain,
                  const double                & a_dx,
                  int                           a_ghost,
                  bool                          a_mergeSmallCells,
                  bool                          a_printStuff)
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
        (*m_graphs)[dit[ibox]].define(a_geoserv, a_domain, dit[ibox], a_dx, a_ghost,
                                      a_mergeSmallCells, a_printStuff);
      }
      m_graphs->exchange();
    }

    void outputGraphMapAsHDF5(const string & a_filename) const
    {
      LevelData<FArrayBox>  ldfabs(m_grids, 1, IntVect::Zero);
      Chombo4::DataIterator dit = m_grids.dataIterator();

      for(int ibox = 0;  ibox < dit.size(); ibox++)
      {
        auto& fab = ldfabs[dit[ibox]];
        const auto& cmgraph = (*m_graphs)[dit[ibox]];
        if(cmgraph.m_allRegular)
        {
          fab.setVal(-1);
        }
        else if(cmgraph.m_allCovered)
        {
          fab.setVal(-2);
        }
        else
        {
          ch_box chbox = m_grids[dit[ibox]];
          pr_box prbox = ProtoCh::getProtoBox(chbox);
          const auto &graphfab = *(cmgraph.m_map);
          for(auto bit = prbox.begin(); bit != prbox.end(); ++bit)
          {
            auto pt = *bit;
            auto iv = ProtoCh::getIntVect(pt);
            fab(iv, 0) = graphfab(pt, 0);
          }
        } //end has irregular 
          
      } //end loop over boxes

      writeLevelname(&ldfabs, a_filename.c_str());
    }
    shared_ptr<cmgraph_distrib_t>  m_graphs;
    ch_dbl                         m_grids;
    ch_box                         m_domain;
    double                         m_dx;
  }; //end class EBCMMetaDataLevel


  ///
  /**
     Subvolume is a a vector<ebcm_volumes> that corresponds with all the volumes in region defined by a box.
     This can be used to find the volumes within the neighborhood of a particular cell.
  */
  template <int order>
  class SubVolumeVector: public std::vector<EBCM_Volu<order> >
  {

  public:    
    virtual ~SubVolumeVector()
    {
    }
    
    SubVolumeVector(const EBCM_Graph<order> & a_graph, const ch_box& a_region, bool a_printStuff)
    {
      vector<pr_pt>   regularPoints(0);
      vector<int>     irregIndicies(0);
      Proto::Box prregion = ProtoCh::getProtoBox(a_region);
      for (auto it = prregion.begin(); it != prregion.end(); ++it)
      {
        if(a_graph.isRegular(*it))
        {
          regularPoints.push_back(*it);
        }
      }
      if(a_graph.hasIrregular())
      {
        for(int ivec = 0; ivec < a_graph.m_volumes->size(); ivec++)
        {
          const auto& volu = (*(a_graph.m_volumes))[ivec];
          if(a_region.contains(volu.m_pt))
          {
            irregIndicies.push_back(ivec);
          }
        }
      }
      int totalSize = irregIndicies.size() + regularPoints.size();
      this->resize(totalSize);
      for(int ientry = 0; ientry < totalSize; ientry++)
      {
        auto& curVolu =  (*this)[ientry];

        if(ientry < regularPoints.size())
        {
          curVolu.defineAsRegular(regularPoints[ientry], a_graph.m_dx);
        }
        else 
        {
          unsigned int irregInd = ientry - regularPoints.size();
          unsigned int graphind = irregIndicies[irregInd];
          curVolu = (*a_graph.m_volumes)[graphind];
        }
        if(a_printStuff)
        {
          static bool printedFirstOne = false;
          if(!printedFirstOne) //turn this off if you want way too much info
          {
            Chombo4::pout() << "SubVolVec: ientry = " << ientry << endl;
            curVolu.print();
            printedFirstOne = true;
          }
            
        }
      }
    }

    
  private:
    SubVolumeVector();

  }; //end class subvolumevector

}//end namespace EBCM


#endif
