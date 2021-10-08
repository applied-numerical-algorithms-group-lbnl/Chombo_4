
#ifndef _proto_COORDINATESYSTEM_H___
#define _proto_COORDINATESYSTEM_H___

#include <map>
using std::map;
using std::vector;


#include "Proto_IndexTM.H"
#include "Proto_IFSlicer.H"

using std::endl;
using std::cerr;

namespace Proto
{

  template <int dim> class CoordinateSystem
  {
  public:
    // data types for vector of int and double in dimension "dim"
    typedef IndexTM<int,dim>  IvDim;
    typedef IndexTM<double,dim> RvDim;

    // member data
    RvDim m_origin;
    RvDim m_dx;

    // empty constructor
    CoordinateSystem();

    // makes an CoordinateSystem from a origin and dx
    CoordinateSystem(const RvDim  & a_origin,
                     const RvDim  & a_dx);

    // makes a reduced CoordinateSystem
    CoordinateSystem(const CoordinateSystem<dim+1> & a_hISystem,
                     const int                     & a_idir);

    // destructor
    ~CoordinateSystem();

    // convert "a_point" in coordinate system "a_system" to a point in the
    // current coordinate system.
    inline RvDim convert(const RvDim                 & a_point,
                         const CoordinateSystem<dim> & a_system) const;

    // convert "a_coord" in direction "a_dir" in coordinate system "a_system"
    // to a value in "a_dir" in the current coordinate system.
    inline double convertDir(const double                  & a_coord,
                             const CoordinateSystem<dim> & a_system,
                             const int                   & a_dir) const;

  };
}
#include "Proto_CoordinateSystemImplem.H"

#endif
