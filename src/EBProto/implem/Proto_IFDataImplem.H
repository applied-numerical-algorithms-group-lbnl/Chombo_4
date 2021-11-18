#ifndef _Proto_IFDATAIMPLEM_H__
#define _Proto_IFDATAIMPLEM_H__

#include <iostream>
#include <iomanip>

#include "Proto_MultiIndex.H"
#include "Proto_NormalDerivativeNew.H"

using std::cerr;
using std::cout;
using std::endl;
namespace Proto
{
// empty constructor
  template <int dim> IFData<dim>::IFData()
  {
  }




// Constructor from the implicit function at the highest dimension
  template <int dim> IFData<dim>::IFData(const BaseIF & a_function,
                                         const RvDim  & a_dx,
                                         const RvDim  & a_cellCenter,
                                         const int    & a_maxOrder)
    :m_globalCoord    (a_cellCenter,a_dx),
     m_cellCenterCoord(IndexTM<double,dim>::Zero,a_dx),
     m_parentCoord    (IndexTM<double,dim>::Zero,a_dx),
     m_maxOrder       (a_maxOrder)
  {
    PR_TIME("ifdata constructor 2");
    // Copy the implicit function; m_function is a pointer
    // IFData always owns the data pointed to by m_function
    m_function = shared_ptr<IFSlicer<dim> >(new IFSlicer<dim>(a_function));

    // makecornerSigns sets the bools m_allVerticesIn and m_allVerticesOut,
    // and m_allVerticesOn
    makeCornerSigns();

    // find an intersection point, if it exists, of the interface with each edge
    findIntersectionPts();

    // will ultimately assign m_localCoord.m_origin to the the average of
    // intersection points.  If full or empty cell, m_localCoord.m_origin = m_parentCoord.m_origin
    defineLocalCoords();

    // set the normal derivatives
    setNormalDerivatives();
  }

// Constructor using IFSlicer, used for refinement
  template <int dim> IFData<dim>::IFData(IFSlicer<dim> * a_function,
                                         const RvDim   & a_dx,
                                         const RvDim   & a_cellCenter,
                                         const int     & a_maxOrder)
    :m_globalCoord    (a_cellCenter,a_dx),
     m_cellCenterCoord(IndexTM<double,dim>::Zero,a_dx),
     m_parentCoord    (IndexTM<double,dim>::Zero,a_dx),
     m_maxOrder       (a_maxOrder)
  {
    PR_TIME("ifdata constructor 3");
    // IFData always owns the data pointed to by m_function
    m_function = shared_ptr<IFSlicer<dim> >(new IFSlicer<dim>(*a_function));

    // makeCornerSigns sets the bools m_allVerticesIn and m_allVerticesOut,
    // and m_allVerticesOn
    makeCornerSigns();

    // Find an intersection point, if it exists, of the interface with each edge
    findIntersectionPts();

    //
    defineLocalCoords();

    // Set the normal derivatives
    setNormalDerivatives();
  }

  template <int dim> IFData<dim>::IFData(const IFData<dim+1> & a_hIFData,
                                         const int           & a_maxOrder,
                                         const int           & a_idir,
                                         const int           & a_hilo)
  {
    PR_TIME("ifdata constructor 4");
    typedef IndexTM<int,dim+1>HDEdgeIndex;
    typedef map<HDEdgeIndex,double > HDEdgeIntersections;
    typedef IndexTM<int,dim+1>HDVertex;
    typedef map<HDVertex,int > HDCornerSigns;

    // pm = -1 for lo and 1 for hi
    int pm = a_hilo;
    pm *= 2;
    pm -= 1;

    // fixedComp and fixed value are used by the slicer to pull points into one higher dimension
    int fixedComp   = a_idir;
    double fixedValue = a_hIFData.m_globalCoord.convertDir(pm*0.5*a_hIFData.m_cellCenterCoord.m_dx[fixedComp],
                                                         a_hIFData.m_cellCenterCoord,
                                                         fixedComp);

    // IFData always owns the data pointed to by m_function
    m_function = shared_ptr< IFSlicer<dim> >(new IFSlicer<dim>(a_hIFData.m_function,fixedComp,fixedValue));

    // reduce the global coordinate system
    m_globalCoord = CoordinateSystem<dim>(a_hIFData.m_globalCoord,fixedComp);

    // Construct cell centered coord sys
    m_cellCenterCoord = CoordinateSystem<dim>(a_hIFData.m_cellCenterCoord,fixedComp);

    // reduce the parent coordinate system
    m_parentCoord = CoordinateSystem<dim>(a_hIFData.m_localCoord,fixedComp);

    m_allVerticesIn  = true;
    m_allVerticesOut = true;
    m_allVerticesOn  = true;

    for (typename HDCornerSigns::const_iterator it = a_hIFData.m_cornerSigns.begin();
         it != a_hIFData.m_cornerSigns.end(); ++it)
    {
      const HDVertex& hVertex = it->first;
      if (hVertex[fixedComp] == a_hilo)
      {
        Vertex vertex;
        for (int j = 0; j < dim; ++j)
        {
          if (j < fixedComp)
          {
            vertex[j] = hVertex[j];
          }
          else
          {
            vertex[j] = hVertex[j+1];
          }
        }

        m_cornerSigns[vertex] = it->second;
        if (m_cornerSigns[vertex] == IN)
        {
          m_allVerticesOut = false;
          m_allVerticesOn  = false;
        }
        else if (m_cornerSigns[vertex] == OUT)
        {
          m_allVerticesIn = false;
          m_allVerticesOn = false;
        }
      }
    }

    for (typename HDEdgeIntersections::const_iterator it = a_hIFData.m_intersections.begin();
         it != a_hIFData.m_intersections.end(); ++it)
    {
      const HDEdgeIndex& hEdgeIndex = it->first;
      EdgeIndex edgeIndex;
      int hEdgeDir = hEdgeIndex[0];

      // only interested in edgeDir != fixedComp(ie. a_idir)
      // Among these, need to find edges in the higher dim that are a_hilo in the fixedComp component.
      // Example: fixedComp = 1; a_hilo = 0. Observe that an edge in the x (N.B 0 < 1) direction
      // might look like (edgeDir= 0:yhilo,zhilo ). In other words, y information is in index 1.
      // On the other hand, an edge in the z direction(N.B 2 > 1)
      // might look like (edgeDir= 0:xhilo,yhilo ). In other words, y information is in index 2.
      if (hEdgeDir < fixedComp)
      {
        if (hEdgeIndex[fixedComp] == a_hilo)
        {
          edgeIndex[0] = hEdgeDir;

          for (int j = 1; j < dim; ++j)
          {
            if (j < fixedComp)
            {
              edgeIndex[j] = hEdgeIndex[j];
            }
            else
            {
              edgeIndex[j] = hEdgeIndex[j+1];
            }
          }
          m_intersections[edgeIndex] = it->second;
        }
      }
      else if (hEdgeDir > fixedComp)
      {
        if (hEdgeIndex[fixedComp+1] == a_hilo)
        {
          edgeIndex[0] = hEdgeDir - 1;

          for (int j = 1; j < dim; ++j)
          {
            if (j < fixedComp+1)
            {
              edgeIndex[j] = hEdgeIndex[j];
            }
            else
            {
              edgeIndex[j] = hEdgeIndex[j+1];
            }
          }
          m_intersections[edgeIndex] = it->second;
        }
      }
    }

    // assigns m_localCoord.m_origin to the the average of intersection points
    // or to m_parentCoord.m_origin if there are no interesection points
    defineLocalCoords();

    m_maxOrder = a_maxOrder;

    if (!m_allVerticesOut && !m_allVerticesIn)
    {
      setNormalDerivatives();
    }
    else
    {
      m_badNormal = false;
    }
  }

// Destructor
  template <int dim> IFData<dim>::~IFData()
  {
//    if (m_function != NULL)
//    {
//      delete m_function;
//    }
  }

  template <int dim> void IFData<dim>::setNormalDerivatives()
  {
    PR_TIME("ifdata::setNormalDerivatives");
    NormalDerivativeNew<dim> normalDerivative;
    typename NormalDerivativeNew<dim>::NormalDerivativeMap ndMap = 
      normalDerivative.calculateAll(m_maxOrder,
                                    m_globalCoord.convert(RvDim::Zero,m_localCoord),
                                    &(*m_function));


    for (int order = 0; order <= m_maxOrder; order++)
    {
      std::vector<IvDim> derivatives;

      generateMultiIndices(derivatives,order);

      for (int i = 0; i < derivatives.size(); ++i)
      {
        const IvDim & derivative = derivatives[i];

        for (int idir = 0; idir < dim; idir++)
        {
          m_normalDerivatives[derivative][idir] = ndMap[derivative][idir];
        }
      }
    }

    //took out a bad Normal Check  --- normalDerivativeNew does not have a getMagnitudeOfGradient
    m_badNormal = false;
  }

  template <int dim> void IFData<dim>::makeCornerSigns()
  {
    PR_TIME("ifdata::makeCornerSigns");
    m_allVerticesIn  = true;
    m_allVerticesOut = true;
    m_allVerticesOn  = true;

    Vertex vertex;

    // there are 2^dim vertices
    int numVertices = 1 << dim;

    for (int i = 0; i < numVertices; ++i)
    {
      int ii = i;

      //label the ith vertex by i represented in base 2
      for (int j = 0; j < dim; ++j)
      {
        //convert to base 2 by successivly "anding" with 1 and shifting
        vertex[j] = (ii & 1);
        ii = ii >> 1;
      }

      // represent the vertex as an RvDim in cell centered coordinates
      RvDim corner;
      for (int idir = 0; idir < dim; ++idir)
      {
        corner[idir] = vertex[idir] - 0.5;
        corner[idir] *= m_cellCenterCoord.m_dx[idir];
      }

      // compute coordinates of a_vertex in global coordinates
      RvDim cornerCoord = m_globalCoord.convert(corner,m_cellCenterCoord);

      // evaluate (zeroth derivative of) the function at the corner
      double val = m_function->value(IndexTM<int,dim>::Zero,cornerCoord);

      // true = negative = in the computational domain = in the fluid
      if (val < -MACHINEPRECISION)
      {
        m_cornerSigns[vertex] = IN;
        m_allVerticesOut      = false;
        m_allVerticesOn       = false;
      }
      else if (val > MACHINEPRECISION)
      {
        m_cornerSigns[vertex] = OUT;
        m_allVerticesIn       = false;
        m_allVerticesOn       = false;
      }
      else
      {
        m_cornerSigns[vertex] = ON;
      }
    }
  }

  template <int dim> void IFData<dim>::findIntersectionPts()
  {
    PR_TIME("ifdata::findIntersectionPoints");
    // double iterate through corner signs map,writing to a edge map if corner1 is connected to corner2
    for (typename CornerSigns::const_iterator it0 = m_cornerSigns.begin();
         it0 != m_cornerSigns.end(); ++it0)
    {
      for (typename CornerSigns::const_iterator it1 = m_cornerSigns.begin();
           it1 != m_cornerSigns.end(); ++it1)
      {
        int edgeDir = LARGEINTVAL;
        if (isConnected(edgeDir,it0->first,it1->first))
        {
          // make edge key:m_intersections[thisEdge] = pt;
          makeEdgeKey(edgeDir,it0->first,it1->first);
        }
      }
    }
  }
  template <int dim> void IFData<dim>::defineLocalCoords()
  {
    PR_TIME("ifdata::defineLocalCoords");
    // default for full cells
    m_localCoord.m_origin = m_cellCenterCoord.m_origin;
    /**
       the stuff below makes convergence tests do weird things
       but makes for better answers
       dtg-- this is broken and needs to go away --dtg
    */
#if DIM==2
    if(1)
    {
      if (m_intersections.size() == 0)
      {
        m_localCoord.m_origin = m_cellCenterCoord.m_origin;
      }
      else
      {
        m_localCoord.m_origin = RvDim::Zero;

        for (typename EdgeIntersections::const_iterator it = m_intersections.begin();it != m_intersections.end(); ++it)
        {
          const EdgeIndex& edgeIndex = it->first;
          const double&      intercept = it->second;

          int varyDir = edgeIndex[0];

          RvDim intersectPt = RvDim::Zero;
          intersectPt[varyDir] = intercept;

          for (int i = 1; i < dim; i++)
          {
            int curDir;
            int loHi;

            if (i <= varyDir)
            {
              curDir = i-1;
            }
            else
            {
              curDir = i;
            }

            loHi = edgeIndex[i];
            intersectPt[curDir]  = loHi - 0.5;
            intersectPt[curDir] *= m_globalCoord.m_dx[curDir];
          }

          m_localCoord.m_origin -= intersectPt;
        }

        m_localCoord.m_origin /= m_intersections.size();
      }
    }
#endif
    /**/
    m_localCoord.m_dx =  m_globalCoord.m_dx;
  }

  template <int dim> bool IFData<dim>::isConnected(int          & a_edgeDir,
                                                   const Vertex & a_vertex1,
                                                   const Vertex & a_vertex2)
  {
    PR_TIME("ifdata::isConnected");
    bool connected = true;
    int numDif = 0;

    // connected = true if and only if  (a_vertex1 and a_vertex2 differ in
    // exactly one coordinate.)
    for (int idir = 0; idir < dim; ++idir)
    {
      if (a_vertex1[idir] != a_vertex2[idir])
      {
        // vertices differ in idir direction
        a_edgeDir = idir;
        numDif += 1;
      }
    }

    if (numDif == 1)
    {
      connected = true;
    }
    else
    {
      connected = false;
      a_edgeDir = LARGEINTVAL;
    }

    return connected;
  }

  template <int dim> void IFData<dim>::makeEdgeKey(const int    & a_edgeDir,
                                                   const Vertex & a_vertex0,
                                                   const Vertex & a_vertex1)
  {
    PR_TIME("ifdata::makeedgekey");
    EdgeIndex thisEdge;

    thisEdge[0] = a_edgeDir;

    for (int idir = 0; idir < dim; ++idir)
    {
      if (idir < a_edgeDir)
      {
        thisEdge[idir + 1] = a_vertex0[idir];
      }
      else if (idir > a_edgeDir)
      {
        thisEdge[idir] = a_vertex0[idir];
      }
    }

    int lo = LARGEINTVAL;
    int hi = LARGEINTVAL;

    if (m_cornerSigns.find(a_vertex0)!= m_cornerSigns.end())
    {
      lo = m_cornerSigns.find(a_vertex0)->second;
    }
    else
    {
      cout << "Vertex not well defined in makeEdgeKey" << endl;
      abort();
    }

    if (m_cornerSigns.find(a_vertex1)!= m_cornerSigns.end())
    {
      hi = m_cornerSigns.find(a_vertex1)->second;
    }
    else
    {
      cout << "Vertex not well defined in makeEdgeKey" << endl;
      abort();
    }

    if ((lo == OUT && hi == IN) || (lo == IN && hi == OUT))
    {
      // calculate a value between -1.0 and 1.0
      double pt = rootFinder(thisEdge);

      // check whether intersection is at the end point
      bool hiOn = false;
      bool loOn = false;

      checkIntersection(hiOn,loOn,pt);

      if (hiOn || loOn)
      {
        if (loOn)
        {
          m_cornerSigns[a_vertex0] = ON;
          lo = ON;
        }

        if (hiOn)
        {
          m_cornerSigns[a_vertex1] = ON;
          hi = ON;
        }

        remakeCornerSigns();
      }
      else
      {
        // Scale the intersection between -dx/2 and dx/2
        m_intersections[thisEdge] = 0.5*m_cellCenterCoord.m_dx[thisEdge[0]]*pt;
      }
    }

    // If the edge in full and the EB intersects a vertex, record an
    // intersection point
    if (lo == IN && hi == ON)
    {
      // Intersection at -dx/2
      m_intersections[thisEdge] = -0.5*m_cellCenterCoord.m_dx[thisEdge[0]];
    }
    else if (lo == ON && hi == IN)
    {
      // Intersection at dx/2
      m_intersections[thisEdge] =  0.5*m_cellCenterCoord.m_dx[thisEdge[0]];
    }
  }

  template <int dim> double IFData<dim>::rootFinder(const EdgeIndex& a_thisEdge)
  {
    PR_TIME("ifdata::rootfinder");
    double pt = LARGEREALVAL;

    // the first component of an EdgeIndex gives the direction of the edge
    // the other components select an edge by using the same (lo-hi) logic as Vertex
    int edgeDir = a_thisEdge[0];

    Vertex loEnd;
    loEnd[edgeDir] = 0;

    Vertex hiEnd;
    hiEnd[edgeDir] = 1;

    for (int idir = 0; idir < dim; ++idir)
    {
      if (idir < edgeDir)
      {
        loEnd[idir] = a_thisEdge[idir + 1];
        hiEnd[idir] = a_thisEdge[idir + 1];
      }
      else if (idir > edgeDir)
      {
        loEnd[idir] = a_thisEdge[idir];
        hiEnd[idir] = a_thisEdge[idir];
      }
    }

    // computes coordinates of a_vertex relative to the cell-center = (0,0...0)
    // represent the vertex as an RvDim in cell centered coordinates
    RvDim loCorner;
    RvDim hiCorner;
    for (int idir = 0; idir < dim; ++idir)
    {
      loCorner[idir] = (loEnd[idir] - 0.5) * m_cellCenterCoord.m_dx[idir];
      hiCorner[idir] = (hiEnd[idir] - 0.5) * m_cellCenterCoord.m_dx[idir];
    }

    // compute coordinates of a_vertex in global coordinates
    RvDim loPt = m_globalCoord.convert(loCorner,m_cellCenterCoord);
    RvDim hiPt = m_globalCoord.convert(hiCorner,m_cellCenterCoord);

    // range of possible roots
    // double smallestRoot = -0.5;
    // double biggestRoot  =  0.5;

    // returns a number between -1.0 and 1.0
    pt = BrentRootFinder(loPt,hiPt,edgeDir); // ,smallestRoot,biggestRoot);

    // return this result
    return pt;
  }

  template <int dim> double IFData<dim>::BrentRootFinder(const RvDim & a_loPt,
                                                       const RvDim & a_hiPt,
                                                       const int   & a_edgeDir) const
  //                                                     const double  & a_smallestRoot,
  //                                                     const double  & a_biggestRoot) const
  {
    PR_TIME("ifdata::brentrootfinder");
    // Max allowed iterations and floating point precision
    const unsigned int MAXITER = 100;
    const double         EPS   = 3.0e-15;

    unsigned int i;
    double aPt;
    double bPt;
    double c, fa, fb, fc;
    double d, e;
    double tol1, xm;
    double p, q, r, s;

    aPt = -1.0; // a_smallestRoot;
    bPt =  1.0; // a_biggestRoot;

    RvDim physCoordAPt = a_loPt;
    RvDim physCoordBPt = a_hiPt;

    fa = -m_function->value(IndexTM<int,dim>::Zero,physCoordAPt);
    fb = -m_function->value(IndexTM<int,dim>::Zero,physCoordBPt);

    //  Init these to be safe
    c = d = e = 0.0;

    if (fb*fa > 0)
    {
      std::cout << "fa " << fa << " fb " << fb << std::endl;
      PR_error("IFData::BrentRootFinder. Root must be bracketed, but instead the supplied end points have the same sign.");
    }

    fc = fb;

    for (i = 0; i < MAXITER; i++)
    {
      if (fb*fc > 0)
      {
        //  Rename a, b, c and adjust bounding interval d
        c = aPt;
        fc  = fa;
        d = bPt - aPt;
        e = d;
      }

      if (std::abs(fc) < std::abs(fb))
      {
        aPt = bPt;
        bPt = c;
        c = aPt;
        fa  = fb;
        fb  = fc;
        fc  = fa;
      }

      //  Convergence check
      tol1  = 2.0 * EPS * std::abs(bPt) + 0.5 * TOLERANCE;
      xm    = 0.5 * (c - bPt);

      if (std::abs(xm) <= tol1 || fb == 0.0)
      {
        break;
      }

      if (std::abs(e) >= tol1 && std::abs(fa) > std::abs(fb))
      {
        //  Attempt inverse quadratic interpolation
        s = fb / fa;
        if (aPt == c)
        {
          p = 2.0 * xm * s;
          q = 1.0 - s;
        }
        else
        {
          q = fa / fc;
          r = fb / fc;
          p = s * (2.0 * xm * q * (q-r) - (bPt-aPt) * (r-1.0));
          q = (q-1.0) * (r-1.0) * (s-1.0);
        }

        //  Check whether in bounds
        if (p > 0) q = -q;

        p = std::abs(p);

        if (2.0 * p < std::min(((double)3.0)*xm*q-std::abs(tol1*q), std::abs(e*q)))
        {
          //  Accept interpolation
          e = d;
          d = p / q;
        }
        else
        {
          //  Interpolation failed, use bisection
          d = xm;
          e = d;
        }
      }
      else
      {
        //  Bounds decreasing too slowly, use bisection
        d = xm;
        e = d;
      }

      //  Move last best guess to a
      aPt = bPt;
      fa  = fb;

      //  Evaluate new trial root
      if (std::abs(d) > tol1)
      {
        bPt = bPt + d;
      }
      else
      {
        if (xm < 0) bPt = bPt - tol1;
        else        bPt = bPt + tol1;
      }

      physCoordBPt[a_edgeDir] = ((1.0 - bPt)/2.0) * a_loPt[a_edgeDir] + ((1.0 + bPt)/2.0) * a_hiPt[a_edgeDir];
      fb = -m_function->value(IndexTM<int,dim>::Zero,physCoordBPt);
    }

    if (i >= MAXITER)
    {
      cerr << "BrentRootFinder: exceeding maximum iterations: "
           << MAXITER
           << "\n";
    }

    return bPt;
  }

  template <int dim> void IFData<dim>::checkIntersection(bool       & a_hiOn,
                                                         bool       & a_loOn,
                                                         const double & a_pt) const
  {
    a_hiOn = false;
    a_loOn = false;

    if (a_pt >= 1.0 - MACHINEPRECISION)
    {
      a_hiOn = true;
    }
    else if (a_pt <= -1.0 + MACHINEPRECISION)
    {
      a_loOn = true;
    }
  }

  template <int dim> void IFData<dim>::remakeCornerSigns()
  {
    m_allVerticesIn  = true;
    m_allVerticesOut = true;
    m_allVerticesOn  = true;

    for (typename CornerSigns::const_iterator it = m_cornerSigns.begin();
         it != m_cornerSigns.end(); ++it)
    {
      if (it->second == IN)
      {
        m_allVerticesOut = false;
      }

      if (it->second == OUT)
      {
        m_allVerticesIn = false;
      }

      if (it->second != ON)
      {
        m_allVerticesOn = false;
      }
    }
  }


// empty constructor (dim == 1)
  IFData<1>::IFData()
  {
  }


// Constructor from the implicit function
  IFData<1>::IFData(const IFData<2> & a_2DIFData,
                    const int       & a_maxOrder,
                    const int       & a_idir,
                    const int       & a_hilo)
    :m_globalCoord(a_2DIFData.m_globalCoord,a_idir),
     m_cellCenterCoord(a_2DIFData.m_cellCenterCoord,a_idir),
     m_parentCoord(a_2DIFData.m_localCoord,a_idir)
  {
    // we want the edge on the a_hilo side of the square with normal in the
    // a_idir direction
    IFData<2>localInfo = a_2DIFData;

    // This 2D edgeIndex locates the 1D edge in the edgeIntersection map
    IndexTM<int,2>twoDEdge;
    twoDEdge[0] = (a_idir + 1)%2;
    twoDEdge[1] = a_hilo;

    m_intersection = LARGEREALVAL;
    if (localInfo.m_intersections.find(twoDEdge) != localInfo.m_intersections.end())
    {
      m_intersection = localInfo.m_intersections[twoDEdge];
    }

    // This 2D vertex locates the hi and lo ends of the 1D segment in the
    // cornerSigns map
    IndexTM<int,2>loPt2D;
    loPt2D[(a_idir + 1)%2] = 0;
    loPt2D[a_idir] = a_hilo;

    IndexTM<int,2>hiPt2D;
    hiPt2D[(a_idir+ 1)%2] = 1;
    hiPt2D[a_idir] = a_hilo;

    if (localInfo.m_cornerSigns.find(loPt2D) != localInfo.m_cornerSigns.end())
    {
      m_cornerSigns[0] = localInfo.m_cornerSigns[loPt2D];
    }
    else
    {
      PR_error("Lo endpoint not in Map");
    }

    if (localInfo.m_cornerSigns.find(hiPt2D) != localInfo.m_cornerSigns.end())
    {
      m_cornerSigns[1] = localInfo.m_cornerSigns[hiPt2D];
    }
    else
    {
      PR_error("Hi endpoint not in Map");
    }

    // set bools
    m_allVerticesIn  = true;
    m_allVerticesOut = true;
    m_allVerticesOn  = true;

    if (m_cornerSigns[0] != ON || m_cornerSigns[1] != ON)
    {
      m_allVerticesOn  = false;
    }

    if (m_cornerSigns[0] == IN || m_cornerSigns[1] == IN)
    {
      m_allVerticesOut = false;
    }

    if (m_cornerSigns[0] == OUT || m_cornerSigns[1] == OUT)
    {
      m_allVerticesIn = false;
    }

    //there is no normal in one dimension. However, if m_badNormal = true at a lower dimension, then the higher dimension refines.
    m_badNormal = false;
  }

// Destructor (dim == 1)
  IFData<1>::~IFData()
  {
  }

}

#endif
