
#ifndef _proto_COORDINATESYSTEMIMPLEM_H_
#define _proto_COORDINATESYSTEMIMPLEM_H_

namespace Proto
{
// empty constructor
  template <int dim> CoordinateSystem<dim>::CoordinateSystem()
  {
  }


// Constructor from an origin and a dx
  template <int dim> CoordinateSystem<dim>::CoordinateSystem( const RvDim  & a_origin,
                                                              const RvDim  & a_dx)
    :m_origin(a_origin),
     m_dx(a_dx)
  {
  }

  template <int dim> CoordinateSystem<dim>::CoordinateSystem(const CoordinateSystem<dim+1> & a_hiCoord,
                                                             const int                     & a_idir)
  {
    PR_TIME("CoordinateSystem::ConstructorReduceInfo");
    for (int idir = 0; idir < dim; ++idir)
    {
      if (idir < a_idir)
      {
        m_origin[idir] = a_hiCoord.m_origin[idir];
        m_dx[idir]     = a_hiCoord.m_dx[idir];
      }
      else
      {
        m_origin[idir] = a_hiCoord.m_origin[idir + 1];
        m_dx[idir]     = a_hiCoord.m_dx    [idir + 1];
      }
    }
  }

// Destructor
  template <int dim> CoordinateSystem<dim>::~CoordinateSystem()
  {
  }

  template <int dim> IndexTM<double,dim> CoordinateSystem<dim>::convert(const RvDim                 & a_point,
                                                                      const CoordinateSystem<dim> & a_system) const
  {
    RvDim retval;

    for (int idir = 0; idir < dim; ++idir)
    {
      // If   a_point = a_system.m_origin
      // then retval  =          m_origin
      //
      // If   a_point = a_system.m_dx[idir]^-1 + a_system.m_origin
      // then retval   =         m_dx[idir]^-1 +          m_origin

      retval[idir]  = a_point[idir];

      retval[idir] -= a_system.m_origin[idir];
      retval[idir] /= a_system.m_dx    [idir];
      retval[idir] *=          m_dx    [idir];
      retval[idir] +=          m_origin[idir];
    }

    return retval;
  }

  template <int dim> double CoordinateSystem<dim>::convertDir(const double                  & a_coord,
                                                            const CoordinateSystem<dim> & a_system,
                                                            const int                   & a_dir) const
  {
    double retval;

    retval  = a_coord;

    retval -= a_system.m_origin[a_dir];
    retval /= a_system.m_dx    [a_dir];
    retval *=          m_dx    [a_dir];
    retval +=          m_origin[a_dir];

    return retval;
  }



}
#endif
