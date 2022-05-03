#ifndef __hoeb_neighborhood__
#define __hoeb_neighborhood__
#include "Proto.H"
#include "Chombo_EBChombo.H"
#include "Chombo_ParmParse.H"
#include "Chombo_EBLevelBoxData.H"
#include "Chombo_GeometryService.H"
#include "Chombo_EBDictionary.H"
#include "Chombo_ParmParse.H"
#include "Hoeb_LAPACKMatrix.H"
#include "implem/Proto_Factorial.H"
#define HOEB_MAX_ORDER 4

namespace hoeb
{
  
  template <CENTERING srccent>
  class Neighborhood
  {
  public:
    //stencil information (so templated on where src not dst)
    vector<EBIndex<srccent> >       m_volumes;
    vector<EBIndex<XFACE> >      m_xfacesDiri;
    vector<EBIndex<YFACE> >      m_yfacesDiri;
    vector<EBIndex<ZFACE> >      m_zfacesDiri;
    vector<EBIndex<XFACE> >      m_xfacesNeum;
    vector<EBIndex<YFACE> >      m_yfacesNeum;
    vector<EBIndex<ZFACE> >      m_zfacesNeum;
    vector<EBIndex<BOUNDARY> > m_cutfacesDiri;
    vector<EBIndex<BOUNDARY> > m_cutfacesNeum;

    //weights used in WLS
    vector<Real>      m_volumesWgt;
    vector<Real>   m_xfacesDiriWgt;
    vector<Real>   m_yfacesDiriWgt;
    vector<Real>   m_zfacesDiriWgt;
    vector<Real>   m_xfacesNeumWgt;
    vector<Real>   m_yfacesNeumWgt;
    vector<Real>   m_zfacesNeumWgt;
    vector<Real> m_cutfacesDiriWgt;
    vector<Real> m_cutfacesNeumWgt;
    
    //number of degrees of freedom.  need to know this to make
    //sure our matricies have sufficient entries
    unsigned int size() const
    {
      unsigned int retval = 0;
      retval +=      m_volumes.size();
      retval +=   m_xfacesDiri.size();
      retval +=   m_yfacesDiri.size();
      retval +=   m_zfacesDiri.size();
      retval += m_cutfacesDiri.size();
      retval +=   m_xfacesNeum.size();
      retval +=   m_yfacesNeum.size();
      retval +=   m_zfacesNeum.size();
      retval += m_cutfacesNeum.size();
      return retval;
    }

    template<CENTERING cent>
    inline void printVec(string prefix, const vector<Real>& a_wgt, const vector<EBIndex<cent> >& a_index)  const
    {
      PR_assert(a_wgt.size() == a_index.size());
      using Chombo4::pout;
      pout()    << setprecision(8)
                << setiosflags(ios::showpoint)
                << setiosflags(ios::scientific);

      pout() << prefix << ":" << endl;
      for(int ivec = 0; ivec < a_wgt.size(); ivec++)
      {
        pout() << "(" << a_index[ivec].m_pt << "," << a_wgt[ivec] << ")" << endl;
      }
    }
  
    void poutAll()  const
    {
      printVec(string("     m_volumes"),      m_volumesWgt,      m_volumes);
      printVec(string("  m_xfacesDiri"),   m_xfacesDiriWgt,   m_xfacesDiri);
      printVec(string("  m_yfacesDiri"),   m_yfacesDiriWgt,   m_yfacesDiri);
      printVec(string("  m_zfacesDiri"),   m_zfacesDiriWgt,   m_zfacesDiri);
      printVec(string("  m_xfacesNeum"),   m_xfacesNeumWgt,   m_xfacesNeum);
      printVec(string("  m_yfacesNeum"),   m_yfacesNeumWgt,   m_yfacesNeum);
      printVec(string("  m_zfacesNeum"),   m_zfacesNeumWgt,   m_zfacesNeum);
      printVec(string("m_cutfacesDiri"), m_cutfacesDiriWgt, m_cutfacesDiri);
      printVec(string("m_cutfacesNeum"), m_cutfacesNeumWgt, m_cutfacesNeum);
    }
  };
  ///will add boundary elements if their corresponding bc is Dirichlet or Neumann"
  //doing this with sets to avoid double-counting faces
  template <CENTERING cent>
  inline void addNeighborhoodBoundaryFaces(std::set<EBIndex<cent> >          & a_neumfaces,
                                           std::set<EBIndex<cent> >          & a_dirifaces,
                                           const std::vector<EBIndex<cent> > & a_allfaces,
                                           const EBGraph                     & a_graph,
                                           const std::string                 & a_dombcname)
  {
    std::set<EBIndex<cent> >* setptr= NULL;
    if(a_dombcname == string("no_bcs"))
    {
      return;
    }
    else if(a_dombcname == string("Dirichlet"))
    {
      setptr = &a_dirifaces;
    }
    else if(a_dombcname == string("Neumann"))
    {
      setptr = &a_neumfaces;
    }

    for(int ivec = 0; ivec < a_allfaces.size(); ivec++)
    {
      if(a_allfaces[ivec].m_isBoundary)
      {
        setptr->insert(a_allfaces[ivec]);
      }
    }
  }
                                    
  template <CENTERING dstcent, CENTERING srccent>
  inline Neighborhood<srccent>
  getNeighbors(const EBIndex<dstcent>& a_start,
               const EBGraph         & a_graph,
               const string            a_dombcname[2*DIM],
               const string          & a_ebbcname,
               const int             & a_facedir);

  template<CENTERING cent>
  inline RealVect
  getNeighborLocation(EBIndex<cent> a_index)
  {
    Chombo4::MayDay::Error("not implemented");
    return RealVect::Zero();
  }

  template< >
  inline RealVect
  getNeighborLocation(EBIndex<CELL> a_index)
  {
    RealVect retval;
    for(int idir = 0; idir < DIM; idir++)
    {
      retval[idir] = Real(a_index.m_pt[idir]) + 0.5;
    }
    return retval;
  }
  template< >
  inline RealVect
  getNeighborLocation(EBIndex<BOUNDARY> a_index)
  {
    RealVect retval;
    for(int idir = 0; idir < DIM; idir++)
    {
      retval[idir] = Real(a_index.m_pt[idir]) + 0.5;
    }
    return retval;
  }
  template< >
  inline RealVect
  getNeighborLocation(EBIndex<XFACE> a_index)
  {
    RealVect retval;
    for(int idir = 0; idir < DIM; idir++)
    {
      retval[idir] = Real(a_index.m_pt[idir]) + 0.5;
    }
    retval[0] -= 0.5;
    return retval;
  }
  template< >
  inline RealVect
  getNeighborLocation(EBIndex<YFACE> a_index)
  {
    RealVect retval;
    for(int idir = 0; idir < DIM; idir++)
    {
      retval[idir] = Real(a_index.m_pt[idir]) + 0.5;
    }
    retval[1] -= 0.5;
    return retval;
  }
#if DIM==3 
  template< >
  inline RealVect
  getNeighborLocation(EBIndex<ZFACE> a_index)
  {
    RealVect retval;
    for(int idir = 0; idir < DIM; idir++)
    {
      retval[idir] = Real(a_index.m_pt[idir]) + 0.5;
    }
    retval[2] -= 0.5;
    return retval;
  }
#endif  
  ///Cartesian distance with unit dx, bounded by 0.5
  template <CENTERING startcent, CENTERING veccent>
  inline Real 
  getDistance( const EBIndex<veccent  >  & a_end,
               const EBIndex<startcent>  & a_start)
  {
    //doing Cartesian distance
    Real sqval = 0;
    RealVect startloc = getNeighborLocation<startcent>(a_start);
    RealVect   endloc = getNeighborLocation<veccent>(a_end);

    for(int idir =0; idir < DIM; idir++)
    {
      Real dist = endloc[idir] -startloc[idir];
      sqval += dist*dist;
    }
    Real retval = sqrt(sqval);
    retval = std::max(retval, 0.5);
    return retval;
  }
  //done as dx=1
  template <CENTERING veccent, CENTERING startcent>
  inline void
  setWgtSection( vector<Real >                   & a_wgt,
                 const vector<EBIndex<veccent> > & a_indices, 
                 const        EBIndex<startcent> & a_start)
  {
    int ipow = 5;
    ParmParse pp;
    pp.query("neighborhood_weight_power", ipow);
    a_wgt.resize(a_indices.size());
    for(int ivec = 0 ; ivec < a_indices.size(); ivec++)
    {
      Real distance = getDistance<startcent, veccent>(a_indices[ivec], a_start);
      Real invwgt = Proto::POW(distance, ipow);
      a_wgt[ivec] = 1.0/invwgt;
    }
  }
  //done as dx=1
  template <CENTERING srccent, CENTERING dstcent>
  inline void
  setNeighborWeights( Neighborhood<srccent>    & a_block,
                      const EBIndex<dstcent>   & a_start)
  {
    setWgtSection<srccent ,dstcent>(a_block.     m_volumesWgt, a_block.     m_volumes,  a_start);
    setWgtSection<XFACE   ,dstcent>(a_block.  m_xfacesDiriWgt, a_block.  m_xfacesDiri,  a_start);
    setWgtSection<YFACE   ,dstcent>(a_block.  m_yfacesDiriWgt, a_block.  m_yfacesDiri,  a_start);
    setWgtSection<ZFACE   ,dstcent>(a_block.  m_zfacesDiriWgt, a_block.  m_zfacesDiri,  a_start);
    setWgtSection<XFACE   ,dstcent>(a_block.  m_xfacesNeumWgt, a_block.  m_xfacesNeum,  a_start);
    setWgtSection<YFACE   ,dstcent>(a_block.  m_yfacesNeumWgt, a_block.  m_yfacesNeum,  a_start);
    setWgtSection<ZFACE   ,dstcent>(a_block.  m_zfacesNeumWgt, a_block.  m_zfacesNeum,  a_start);
    setWgtSection<BOUNDARY,dstcent>(a_block.m_cutfacesDiriWgt, a_block.m_cutfacesDiri,  a_start);
    setWgtSection<BOUNDARY,dstcent>(a_block.m_cutfacesNeumWgt, a_block.m_cutfacesNeum,  a_start);
  }
  template <CENTERING dstcent>
  inline Neighborhood<CELL>
  getCellNeighbors(const EBIndex<dstcent>& a_start,
                   const EBGraph         & a_graph,
                   const string            a_dombcname[2*DIM],
                   const string          & a_ebbcname,
                   const int             & a_facedir)
  {
    //need at least as many neighbors as we have coefficients.
    static const int minsize = IndexedMoments<DIM, HOEB_MAX_ORDER>::size();
    int radius = 3;
    int maxrad = 4;
    ParmParse pp;
    pp.get("default_neighbor_radius", radius);
    pp.get("maximum_neighbor_radius", maxrad);
    PR_assert(radius <= maxrad);
    bool hadEnough = false;
    Neighborhood <CELL> retval;
    while(!hadEnough)
    {
      Point lo = a_start.m_pt;
      int facedir = max(a_facedir, 0);
      Point hi = lo + Point::Basis(facedir, 1);
      Bx ptbx(lo, hi);
      Bx neibx = ptbx.grow(radius);

      neibx &= a_graph.getDomain();
      auto  & neighvofs = retval.m_volumes;
      std::set<EBIndex<XFACE> > xfacesetDiri;
      std::set<EBIndex<XFACE> > xfacesetNeum;
      std::set<EBIndex<YFACE> > yfacesetDiri;
      std::set<EBIndex<YFACE> > yfacesetNeum;
#if DIM==3
      std::set<EBIndex<ZFACE> > zfacesetDiri;
      std::set<EBIndex<ZFACE> > zfacesetNeum;
#endif        
        
      neighvofs.clear();
      for(auto bit = neibx.begin(); bit != neibx.end(); ++bit)
      {
        auto vofs = a_graph.getVoFs(*bit);
        //first do the volumes
        neighvofs.insert(neighvofs.end(), vofs.begin(), vofs.end());
          
        //now we do the EBBC
        for(int ivof = 0; ivof < vofs.size(); ivof++)
        {
          const auto& vof = vofs[ivof];
          if(a_graph.isIrregular(vof.m_pt))
          {
            EBIndex<BOUNDARY> cutface = vof.getCutFace();
            if(a_ebbcname == string("no_bcs"))
            {
              //do not add in this case
            }
            else if(a_ebbcname == string("Dirichlet"))
            {
              retval.m_cutfacesDiri.push_back(cutface);
            }
            else if(a_ebbcname == string("Neumann"))
            {
              retval.m_cutfacesNeum.push_back(cutface);
            }
            else
            {
              PROTO_ASSERT(false, "bogus ebbcname in get neighbors");
            }
          }

          //now for the domain boundary faces
          for(SideIterator sit; sit.ok(); ++sit)
          {
            auto xfaces = a_graph.getXFaces(vof, sit());
            auto yfaces = a_graph.getYFaces(vof, sit());
            
            auto xindex = ebp_index(0, sit());
            auto yindex = ebp_index(1, sit());
            addNeighborhoodBoundaryFaces<XFACE>(xfacesetNeum, xfacesetDiri,  xfaces,
                                                a_graph, a_dombcname[xindex]);
            addNeighborhoodBoundaryFaces<YFACE>(yfacesetNeum, yfacesetDiri,  yfaces,
                                                a_graph, a_dombcname[yindex]);
#if DIM==3                                                                    
            auto zfaces = a_graph.getZFaces(vof, sit());
            auto zindex = ebp_index(2, sit());
            addNeighborhoodBoundaryFaces<ZFACE>(zfacesetNeum, zfacesetDiri,  zfaces,
                                                a_graph, a_dombcname[zindex]);
#endif              
          }
        } //end loop over vofs in the cell
          
      } //end loop over the box

      //insert the sets into the vectors
      retval.m_xfacesDiri.assign(xfacesetDiri.begin(), xfacesetDiri.end());
      retval.m_xfacesNeum.assign(xfacesetNeum.begin(), xfacesetNeum.end());
      retval.m_yfacesDiri.assign(yfacesetDiri.begin(), yfacesetDiri.end());
      retval.m_yfacesNeum.assign(yfacesetNeum.begin(), yfacesetNeum.end());
#if DIM==3
      retval.m_zfacesDiri.assign(zfacesetDiri.begin(), zfacesetDiri.end());
      retval.m_zfacesNeum.assign(zfacesetNeum.begin(), zfacesetNeum.end());
#endif        
      
      //if we have too few neighbors, expand the radius
      hadEnough =((radius >= maxrad) || (retval.size() >= minsize));
      if(!hadEnough)
      {
        radius++;
      }
    }
    setNeighborWeights(retval, a_start);
    return retval;
  }

  inline int dummyFunc()
  {
    Neighborhood<CELL> bla;
    bla.poutAll();
    return 0;
  }
}

#endif
