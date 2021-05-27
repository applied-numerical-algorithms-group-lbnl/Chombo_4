#ifndef __hoeb_neighborhood__
#define __hoeb_neighborhood__
#include "Proto.H"
#include "Chombo_EBChombo.H"
#include "Chombo_ParmParse.H"
#include "Chombo_EBLevelBoxData.H"
#include "Chombo_GeometryService.H"
#include "Chombo_EBDictionary.H"
#include "Hoeb_LAPACKMatrix.H"
#define HOEB_MAX_ORDER 4

namespace hoeb
{
  struct Neighborhood
  {
  public:
    vector<EBIndex<CELL> >       m_volumes;
    vector<EBIndex<XFACE> >      m_xfacesDiri;
    vector<EBIndex<YFACE> >      m_yfacesDiri;
    vector<EBIndex<ZFACE> >      m_zfacesDiri;
    vector<EBIndex<XFACE> >      m_xfacesNeum;
    vector<EBIndex<YFACE> >      m_yfacesNeum;
    vector<EBIndex<ZFACE> >      m_zfacesNeum;
    vector<EBIndex<BOUNDARY> > m_cutfacesDiri;
    vector<EBIndex<BOUNDARY> > m_cutfacesNeum;
    //starting location    
    RealVect a_xbar;
    //number of degrees of freedom.  need to know this to make
    //sure our matricies have sufficient entries
    unsigned int size()
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
                                    
  template <CENTERING cent>
  inline Neighborhood
  getNeighbors(const EBIndex<cent>& a_start,
               const EBGraph      & a_graph,
               const string         a_dombcname[2*DIM],
               const string       & a_ebbcname,
               const int          & a_facedir)
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
    Neighborhood retval;
    while(!hadEnough)
    {
      Point lo = a_start.m_pt;
      Point hi = lo + Point::Basis(a_facedir, 1);
      Bx ptbx(lo, hi);
      Bx neibx = ptbx.grow(radius);

      neibx &= a_graph.getDomain();
      auto neighvofs = retval.m_volumes;
      std::set<EBIndex<XFACE> > xfacesetDiri;
      std::set<EBIndex<XFACE> > xfacesetNeum;
      std::set<EBIndex<YFACE> > yfacesetDiri;
      std::set<EBIndex<YFACE> > yfacesetNeum;
#if DIM==3
      std::set<EBIndex<ZFACE> > zfacesetDiri;
      std::set<EBIndex<ZFACE> > zfacesetneum;
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
      retval.m_yfacesDiri.assign(yfacesetNeum.begin(), yfacesetNeum.end());
      retval.m_yfacesNeum.assign(yfacesetDiri.begin(), yfacesetDiri.end());
#if DIM==3
      retval.m_zfacesDiri.assign(zfacesetDiri.begin(), zfacesetDiri.end());
      retval.m_zfacesNeum.assign(zfacesetneum.begin(), zfacesetneum.end());
#endif        
      
      //if we have too few neighbors, expand the radius
      hadEnough =((radius >= maxrad) || (retval.size() >= minsize));
      if(!hadEnough)
      {
        radius++;
      }
    }
    return retval;
  }

}
#include "Hoeb_UtilitiesImplem.H"
#endif
