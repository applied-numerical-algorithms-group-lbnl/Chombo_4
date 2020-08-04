#ifndef __PROTO_IFdata_H__
#define __PROTO_IFdata_H__

#include <map>
using std::map;
using std::vector;

#include "Proto_IndexTM.H"
#include "Proto_IFSlicer.H"
#include "Proto_CoordinateSystem.H"

using std::endl;
using std::cerr;
using std::shared_ptr;
namespace Proto
{
  template <int dim> class IFData
  {
  public:
    // data types to record vertex signs and edge intersections
    typedef IndexTM<int,dim>                      EdgeIndex;
    typedef map<EdgeIndex,double > EdgeIntersections;

    typedef IndexTM<int,dim>                      Vertex;
    typedef map<Vertex,int >        CornerSigns;

    typedef IndexTM<int,dim>                      IvDim;
    typedef IndexTM<double,dim>                     RvDim;

    typedef map<IvDim,RvDim >        NormalDerivatives;

    // member data
    CornerSigns           m_cornerSigns;
    EdgeIntersections     m_intersections;
    shared_ptr<IFSlicer<dim> >        m_function;

    CoordinateSystem<dim> m_globalCoord;
    CoordinateSystem<dim> m_cellCenterCoord;
    CoordinateSystem<dim> m_parentCoord;
    CoordinateSystem<dim> m_localCoord;

    int                   m_maxOrder;
    NormalDerivatives     m_normalDerivatives;
    bool                  m_badNormal;

    bool                  m_allVerticesIn;
    bool                  m_allVerticesOut;
    bool                  m_allVerticesOn;

    // empty constructor
    inline
    IFData();

    // makes an initial IFData from a BaseIF
    inline
    IFData(const BaseIF & a_function,
           const RvDim  & a_dx,
           const RvDim  & a_origin,
           const int    & a_maxOrder);

    // makes an lower dimensional IFData from a IFSlicer; used for refinement (or it used to be)
    inline
    IFData(IFSlicer<dim> * a_function,
           const RvDim   & a_dx,
           const RvDim   & a_origin,
           const int     & a_maxOrder);



    // makes a reduced IFData
    inline
    IFData(const IFData<dim+1> & a_hInfo,
           const int           & a_maxOrder,
           const int           & a_idir,
           const int           & a_hilo);

    // destructor
    inline
    ~IFData();

    inline void setNormalDerivatives();
    inline void makeCornerSigns();
    inline void findIntersectionPts();
    inline void defineLocalCoords();

    // used by the constructor
    inline bool isConnected(int          & a_edgeDir,
                            const Vertex & a_vertex1,
                            const Vertex & a_vertex2);

    inline void  makeEdgeKey(const int    & a_edgeDir,
                             const Vertex & a_vertex1,
                             const Vertex & a_vertex2);

    //this function calls the midPtRootFinder
    inline double rootFinder(const EdgeIndex & a_thisEdge);

    inline double BrentRootFinder(const RvDim & a_loPt,
                                const RvDim & a_hiPt,
                                const int   & a_edgeDir) const;


    inline void checkIntersection(bool       & a_hiOn,
                                  bool       & a_loOn,
                                  const double & a_pt) const;

    inline void remakeCornerSigns();

    // this function sets the changing direction given the normal direction
    inline
    void setChangingDirection();

    inline
    bool searchNormalDir(const int i) const;
  };

  template<> class IFData<1>
  {
  public:
    // true is in the fluid
    typedef double            EdgeIntersection;

    typedef int             Vertex;
    typedef map<Vertex,int> CornerSigns;

    typedef int             IvDim;
    typedef double            RvDim;

    // member data
    CornerSigns         m_cornerSigns;
    EdgeIntersection    m_intersection;

    CoordinateSystem<1> m_globalCoord;
    CoordinateSystem<1> m_cellCenterCoord;
    CoordinateSystem<1> m_parentCoord;

    bool                m_allVerticesIn;
    bool                m_allVerticesOut;
    bool                m_allVerticesOn;
    bool                m_badNormal;

    inline
    IFData();

    inline
    IFData(const IFData<2> & a_hIFData,
           const int       & a_maxOrder,
           const int       & a_idir,
           const int       & a_hilo);

    inline
    ~IFData();

  };
}

#include "Proto_IFDataImplem.H"

#endif
