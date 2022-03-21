
#ifndef ___Proto_SimpleImplicitFunctions__H__
#define ___Proto_SimpleImplicitFunctions__H__

#include "Proto_RealVect.H"
#include "implem/Proto_IndexTM.H"
#include "Proto_BaseIF.H"
namespace Proto
{
  ///truly simple geometries can have some simplified bits about what is all covered or regular
  /**
     For simple geometries, this provides an easy way to encapsulate the writing of a BaseIF
   */
  class SimpleGeometry: public BaseIF
  {
  public:


    IndexTM<double,DIM> getIndexTM(const RealVect& a_rv)
    {
      IndexTM<double, DIM> retval;
      for(int idir = 0; idir < DIM; idir++)
      {
        retval[idir] = a_rv[idir];
      }
      return retval;
    }

    vector<RealVect> getCornerLocations(const Box    & a_box,
                                        const RealVect& a_origin,
                                        const double & a_dx)
    {
      vector<RealVect> cornerlocs;
      Point lopt = a_box.low();
      Point hipt = a_box.high() + Point::Ones(1);

#if DIM==2
      int numcorners = 4;
      cornerlocs.resize(numcorners);

      Point lohi = lopt;
      lohi[1] = hipt[1];   
      Point hilo = hipt;
      hilo[1] = lopt[1];   

      for(int idir = 0; idir < DIM; idir++)
      {
        cornerlocs[0][idir] = a_dx*lopt[idir] + a_origin[idir];
        cornerlocs[1][idir] = a_dx*hipt[idir] + a_origin[idir];
        cornerlocs[2][idir] = a_dx*lohi[idir] + a_origin[idir];
        cornerlocs[3][idir] = a_dx*hilo[idir] + a_origin[idir];
      }
#else
      int numcorners = 8;

      cornerlocs.resize(numcorners);
      Point lolohi = lopt;
      lolohi[2] = hipt[2];   
      Point hihilo = hipt;
      hihilo[2] = lopt[2];   

      Point lohihi = hipt;
      lohihi[0] = lopt[0];   
      Point hilolo = lopt;



      Point lohilo = lopt;
      lolohi[1] = hipt[1];   
      Point hilohi = hipt;
      hihilo[1] = lopt[1];   

      for(int idir = 0; idir < DIM; idir++)
      {
        cornerlocs[0][idir] =   a_dx*lopt[idir];
        cornerlocs[1][idir] =   a_dx*hipt[idir];
        cornerlocs[2][idir] = a_dx*lolohi[idir];
        cornerlocs[3][idir] = a_dx*hihilo[idir];
        cornerlocs[4][idir] = a_dx*hilolo[idir];
        cornerlocs[5][idir] = a_dx*lohihi[idir];
        cornerlocs[6][idir] = a_dx*hilohi[idir];
        cornerlocs[7][idir] = a_dx*lohilo[idir];
      }
#endif


      return cornerlocs;
    }
  };

///  Cylinder.
/**
   (x-xo)^2 + (y-y0)^2 + - r^2 = f(x,y,z)
*/
  class SimpleCylinderIF: public SimpleGeometry
  {
  public:
    ///
    SimpleCylinderIF(const RealVect& a_X0,
                     const double  & a_R,
                     const bool    & a_inside)
    {
      m_X0=       a_X0;    
      m_R=        a_R;     
      m_inside=   a_inside;
    }

    /// Destructor
    virtual ~SimpleCylinderIF()
    {
    }

    double
    value(const IndexTM<double,DIM>& a_point) const
    {
  
      RealVect X;
      //note the 2.   z does not figure in here
      for(int idir = 0; idir < 2; idir++)
      {
        X[idir] = a_point[idir] - m_X0[idir];
      }

      double retval = 0;
      for(int idir = 0; idir < 2; idir++)
      {
        retval += (X[idir]*X[idir]);
      }

      retval -= m_R*m_R;
      if(!m_inside)
        retval = -retval;
      return retval;
    }

    virtual double value(const IndexTM<int,DIM> & a_partialDerivative,
                         const IndexTM<double,DIM>& a_point) const
    {
      int order = a_partialDerivative.sum();
      double retval = LARGEREALVAL;

      //note the 2.   z does not figure in here
      RealVect X;
      for(int idir = 0; idir < 2; idir++)
      {
        X[idir] = a_point[idir] - m_X0[idir];
      }

      if (order == 0)
      {
        retval = value(a_point);
        return retval; //have to return here to avoid double neg
      }
      else if (order == 1)
      {
        bool found = false;
        for(int idir = 0; idir < 2; idir++)
        {
          if(a_partialDerivative[idir] == 1)
          {
            retval = 2*X[idir];
            found = true;
          }
        }
        if(!found) //z deriv
        {
          retval = 0;
        }
      }
      else if (order == 2)
      {
      
        bool found = false;
        //note the 2.   z does not figure in here
        for(int idir = 0; idir < 2; idir++)
        {
          if(a_partialDerivative[idir] == 2)
          {
            found = true;
            retval = 2;
          }
        }
        if(!found) //mixed deriv or z deriv
        {
          retval = 0;
        }
      }
      else
      {
        retval = 0.0;
      }

      if(!m_inside)
        retval = -retval;
      return retval;
    }


    virtual BaseIF* newImplicitFunction() const
    {
      SimpleCylinderIF* SimpleCylinderPtr = new SimpleCylinderIF(m_X0, m_R, m_inside);
      return static_cast<BaseIF*>(SimpleCylinderPtr);
    }

  protected:
    RealVect m_X0;
    double   m_R;
    bool     m_inside;

  private:
    SimpleCylinderIF();

  };
///  Ellipsoid.
/**
   ((x-xo)/A)^2 + ((y-y0)/B)^2 + ((z-z0)/C)^2 - r^2 = f(x,y,z)
*/
  class SimpleEllipsoidIF: public SimpleGeometry
  {
  public:
    ///
    SimpleEllipsoidIF(const RealVect& a_ABE,
                      const RealVect& a_X0,
                      const double  & a_R,
                      const bool    & a_inside)
    {
      m_ABE=      a_ABE;   
      m_X0=       a_X0;    
      m_R=        a_R;     
      m_inside=   a_inside;
      
      for(int idir = 0; idir < DIM; idir++)
      {
        PR_assert(a_ABE[idir] > 1.e-10);
        m_A2[idir] = a_ABE[idir]*a_ABE[idir];
      }
    }

    /// Destructor
    virtual ~SimpleEllipsoidIF()
    {
    }

    double
    value(const IndexTM<double,DIM>& a_point) const
    {
  
      RealVect X;
      for(int idir = 0; idir < DIM; idir++)
      {
        X[idir] = a_point[idir] - m_X0[idir];
      }

      double retval = 0;
      for(int idir = 0; idir < DIM; idir++)
      {
        retval += (X[idir]*X[idir])/m_A2[idir];
      }

      retval -= m_R*m_R;
      if(!m_inside)
        retval = -retval;
      return retval;
    }

    virtual double value(const IndexTM<int,DIM> & a_partialDerivative,
                         const IndexTM<double,DIM>& a_point) const
    {
      int order = a_partialDerivative.sum();
      double retval = LARGEREALVAL;

      RealVect X;
      for(int idir = 0; idir < DIM; idir++)
      {
        X[idir] = a_point[idir] - m_X0[idir];
      }

      if (order == 0)
      {
        retval = value(a_point);
        return retval; //have to return here to avoid double neg
      }
      else if (order == 1)
      {
        for(int idir = 0; idir < DIM; idir++)
        {
          if(a_partialDerivative[idir] == 1)
          {
            retval = 2*X[idir]/m_A2[idir];
          }
        }
      }
      else if (order == 2)
      {
      
        bool found = false;
        for(int idir = 0; idir < DIM; idir++)
        {
          if(a_partialDerivative[idir] == 2)
          {
            found = true;
            retval = 2/m_A2[idir];
          }
        }
        if(!found) //mixed deriv
        {
          retval = 0;
        }
      }
      else
      {
        retval = 0.0;
      }

      if(!m_inside)
        retval = -retval;
      return retval;
    }


    virtual BaseIF* newImplicitFunction() const
    {
      SimpleEllipsoidIF* SimpleEllipsoidPtr = new SimpleEllipsoidIF(m_ABE,m_X0, m_R, m_inside);
      return static_cast<BaseIF*>(SimpleEllipsoidPtr);
    }

  protected:
    RealVect m_ABE;
    RealVect m_X0;
    double   m_R;
    bool     m_inside;
    RealVect m_A2;
  private:
    SimpleEllipsoidIF();

  };


  


/// Sphere
/**
   (x-xo)^2 + (y-y0)^2 + (z-z0)^2 - r^2 = f(x,y,z)

*/
  class SimpleSphereIF: public SimpleEllipsoidIF
  {
  public:

    Box getBoundingBox(RealVect  a_origin,
                       double    a_dx,
                       double    a_radius,
                       RealVect  a_cent)
    {
      Point pcenter;
      double rradius =round(a_radius/a_dx);
      int    pradius =   int(rradius);
      for(int idir = 0; idir < DIM; idir++)
      {
        pcenter[idir] = int((a_origin[idir] + a_cent[idir])/a_dx);
      }
      Box retval = (Box(pcenter, pcenter)).grow(pradius);
      return  retval;
    }
    

    ///return true if the box is entirely in the fluid
    virtual bool entireBoxRegular(const Box     & a_box,
                                  const RealVect& a_origin,
                                  const double  & a_dx)
    {
      bool retval = false;
      if(m_inside)
      {
        //if m_inside, see if the entire region is inside the sphere (shrunk by a bit)
        //could be as big as (3)^-1/2 
        Box insideBox = getBoundingBox(a_origin, a_dx, 0.4*m_R, m_X0);
        retval = insideBox.contains(a_box);
      }
      else
      {
        //otherwise, see if the entire region is outside the sphere (grown by a bit)
        //could be a small as 1.0
        Box outsideBox = getBoundingBox(a_origin, a_dx, 1.1*m_R, m_X0);
        outsideBox &= a_box;
        retval = outsideBox.empty();
      }
      return retval;
    }

    ///return true if the entire box is covered by the ellipsoid
    virtual bool entireBoxCovered(const Box& a_box,
                                  const RealVect& a_origin,
                                  const double& a_dx)
    {
      bool retval = false;
      if(!m_inside)
      {
        //if not inside, see if the entire region is inside the sphere (shrunk by a bit)
        //could be as big as (3)^-1/2 
        Box insideBox = getBoundingBox(a_origin, a_dx, 0.4*m_R, m_X0);
        retval = insideBox.contains(a_box);
      }
      else
      {
        //otherwise, see if the entire region is outside the sphere (grown by a bit)
        //could be a small as 1.0
        Box outsideBox = getBoundingBox(a_origin, a_dx, 1.1*m_R, m_X0);
        outsideBox &= a_box;
        retval = outsideBox.empty();
      }
      return retval;
    }
    ///
    SimpleSphereIF(const RealVect& a_X0,
                   const double  & a_R,
                   const bool    & a_inside):
      SimpleEllipsoidIF(RealVect::Unit(), a_X0, a_R, a_inside)
    {
    }

    /// Destructor
    virtual ~SimpleSphereIF()
    {
    }


  private:
    SimpleSphereIF();

  };

  ///
  /**
     IF where all cells are regular
   */
  class AllRegularIF: public BaseIF
  {
  public:

    /// Default constructor
    AllRegularIF()
    {
    }

    /// Default destructor
    virtual ~AllRegularIF()
    {
    }

    ///
    /**
       Return the value of the function at a_point.  When delineating a domain,
       the level set value=0 represents the boundary and value<0 is inside the
       fluid.
    */
    virtual double value(const IndexTM<int   ,DIM> & a_partialDerivative,
                         const IndexTM<double,DIM> & a_point) const
    {
      double retval;
      if(a_partialDerivative.sum() == 0)
      {
        retval = -1;
      }
      else
      {
        retval = 0;
      }

      return retval;
    }

    
    //.overwrite if you have a fast way of knowing if a whole region is regular
    virtual bool entireBoxRegular(const Box     & a_box,
                                  const RealVect& a_origin,
                                  const double  & a_dx)
    {
      return true;
    }

    ///overwrite if you have a fast way of knowing if a whole region is regular
    virtual bool entireBoxCovered(const Box      & a_box,
                                  const RealVect & a_origin,
                                  const double   & a_dx)
    {
      return false;
    }

    /// factory method
    /**
       Return a newly allocated derived class.  The responsibility
       for deleting the memory is left to the calling function.
    */
    virtual BaseIF* newImplicitFunction() const
    {
      return new AllRegularIF();
    }
      

  };


  ///
  /**
     IF where all cells are regular
   */
  class PlaneIF: public SimpleGeometry
  {
  public:


  public:
    /// Default constructor
    /// Assumes normal facing outward from the interior
    PlaneIF(RealVect a_startPt, RealVect a_normal)
    {
      m_startPt = a_startPt;
      m_normal  = a_normal;
    }

    /// Default destructor
    virtual ~PlaneIF()
    {
    }

    ///
    /**
       Return the value of the function at a_point.  When delineating a domain,
       the level set value=0 represents the boundary and value<0 is inside the
       fluid.
    */
    virtual double value(const IndexTM<int   ,DIM> & a_partialDerivative,
                         const IndexTM<double,DIM> & a_point) const
    {
      double retval = 0;
      if(a_partialDerivative.sum() == 0)
      {
        for(int idir = 0; idir < DIM; idir++)
        {
          retval += -m_normal[idir]*(m_startPt[idir] - a_point[idir]);
        }
      }
      else if(a_partialDerivative.sum() == 1)
      {
        for(int idir = 0; idir < DIM; idir++)
        {
          retval += m_normal[idir]*(a_partialDerivative[idir]);
        }
      }
      else
      {
        retval = 0;
      }


      return retval;
    }

    
    ///return true if the box is entirely in the fluid
    virtual bool entireBoxRegular(const Box& a_box,
                                  const RealVect& a_origin,
                                  const double& a_dx)
    {
      vector<RealVect> cornerlocs = getCornerLocations(a_box, a_origin, a_dx);
      bool allNeg = true;

      IndexTM<int, DIM> deriv;
      deriv.setAll(0);
      for(int icorn = 0; icorn < cornerlocs.size(); icorn++)
      {
        IndexTM<double, DIM> itmloc = getIndexTM(cornerlocs[icorn]);
        double funcvalue = value(deriv, itmloc);
        if(funcvalue > 0) 
        {
          allNeg = false;
        }
      }
      return allNeg;
    }

    ///return true if the entire box is covered by the plane
    virtual bool entireBoxCovered(const Box& a_box,
                                  const RealVect& a_origin,
                                  const double& a_dx)
    {
      vector<RealVect> cornerlocs = getCornerLocations(a_box, a_origin, a_dx);

      bool allPos = true;

      IndexTM<int, DIM> deriv;
      deriv.setAll(0);
      for(int icorn = 0; icorn < cornerlocs.size(); icorn++)
      {
        IndexTM<double, DIM> itmloc = getIndexTM(cornerlocs[icorn]);
        double funcvalue = value(deriv, itmloc);
        if(funcvalue < 0) 
        {
          allPos = false;
        }
      }

      return allPos;
    }
    /// factory method
    /**
       Return a newly allocated derived class.  The responsibility
       for deleting the memory is left to the calling function.
    */
    virtual BaseIF* newImplicitFunction() const
    {
      return static_cast<BaseIF*>(new PlaneIF(m_startPt, m_normal));
    }
    
  private:
    PlaneIF();
    RealVect m_startPt;
    RealVect m_normal;

  };


///  Sin bubbles
/**
   A*sin(pi*a*(x-x0))*sin(pi*b*(y-y0))*sin(pi*c*(z-z0)) - R = f(x,y,z)
*/
  class SimpleSinProductIF: public SimpleGeometry
  {
  public:
    ///
    SimpleSinProductIF(const RealVect a_center,
                       const Real     a_R,
                       const RealVect a_freq,
                       const bool     a_inside)
    {
      m_R=        a_R;
      m_freq = a_freq;
      m_center = a_center;
      m_inside = a_inside;
    }

    /// Destructor
    virtual ~SimpleSinProductIF()
    {
    }

    double
    value(const IndexTM<double,DIM>& a_point) const
      {
        Real retval = 1.;
        for(int idir = 0; idir < DIM; idir++)
          {
            Real xdiff = a_point[idir] - m_center[idir];
            Real factor= sin(xdiff*m_PI*m_freq[idir]);
            retval *= factor;
          }
        retval -= m_R;
        return retval;
      }

    virtual double value(const IndexTM<int,DIM> & a_partialDerivative,
                         const IndexTM<double,DIM>& a_point) const
      {
        Real retval = 1.;
        for(int idir = 0; idir < DIM; idir++)
          {
            Real prefix = 1.0;
            for(int ipre = 0; ipre < a_partialDerivative[idir]; ipre++)
              {
                prefix *= m_PI*m_freq[idir];
              }
            Real arg = m_PI*m_freq[idir]*(a_point[idir]-m_center[idir]);
            if(a_partialDerivative[idir]%2 == 0)
              {
                retval *= pow(-1.,a_partialDerivative[idir]/2)*prefix*sin(arg);
              }
            else
              {
                retval *= pow(-1.,(a_partialDerivative[idir]-1)/2)*prefix*cos(arg);
              }
          }
        if (a_partialDerivative.sum() == 0)
          retval -= m_R;
        if (!m_inside)
          retval *= -1;
        return retval;
      }


    virtual BaseIF* newImplicitFunction() const
    {
      SimpleSinProductIF* SimpleSinProductPtr = new SimpleSinProductIF(m_center, m_R, m_freq, m_inside);
      return static_cast<BaseIF*>(SimpleSinProductPtr);
    }

  protected:
    Real      m_R;
    RealVect  m_freq;
    RealVect  m_center;
    Real      m_PI=3.141592653589793238462643383279502884197;
    bool      m_inside;

  private:
    SimpleSinProductIF();

  };

  ///  Cos surfaces
/**
   SUM_i( PROD_j(A_i*sin(pi*b_ij*(x_j-x0_ij))) = f(x,y,z)
*/
  class SimpleCosSurfIF: public SimpleGeometry
  {
  public:
    ///
    SimpleCosSurfIF(const std::vector<RealVect> a_center,
                    const std::vector<Real>     a_amp,
                    const std::vector<RealVect> a_freq,
                    const Real                  a_R,
                    const bool                  a_inside)
      {
        m_center = a_center;
        m_amp = a_amp;
        m_freq = a_freq;
        m_R = a_R;
        m_inside = a_inside;
        CH_assert(a_center.size() == a_amp.size());
        CH_assert(a_center.size() == a_freq.size());
      }

    /// Destructor
    virtual ~SimpleCosSurfIF()
      {
      }

    double
    value(const IndexTM<double,DIM>& a_point) const
      {
        Real retval = 0.;
        for(int i=0; i!=m_center.size(); i++)
          {
            Real factor = m_amp[i];
            for(int j=0; j!=SpaceDim; j++)
              {
                Real xdiff = a_point[j] - m_center[i][j];
                factor *= cos(xdiff*m_PI*m_freq[i][j]);
              }
            retval += factor;
          }
        retval -= m_R;
        return retval;
      }

    virtual double value(const IndexTM<int,DIM> & a_partialDerivative,
                         const IndexTM<double,DIM>& a_point) const
      {
        Real retval = 0.;
        for(int i=0; i!=m_center.size(); i++)
          {
            Real factor = m_amp[i];
            for(int j=0; j!=SpaceDim; j++)
              {
                Real prefix = 1.;
                for(int ipre = 0; ipre < a_partialDerivative[j]; ipre++)
                  {
                    prefix *= m_PI*m_freq[i][j];
                  }
                Real xdiff = a_point[j] - m_center[i][j];
                Real arg = m_PI*m_freq[i][j]*xdiff;
                if(a_partialDerivative[j]%2 == 0)
                  {
                    factor *= pow(-1.,a_partialDerivative[j]/2)*prefix*cos(arg);
                  }
                else
                  {
                    factor *= pow(-1.,(a_partialDerivative[j]-1)/2)*prefix*sin(arg);
                  }
              }
            retval += factor;
          }
        if (a_partialDerivative.sum() == 0)
          retval -= m_R;
        if (!m_inside)
          retval *= -1;
        return retval;
      }


    virtual BaseIF* newImplicitFunction() const
      {
        SimpleCosSurfIF* SimpleCosSurfPtr = new SimpleCosSurfIF(m_center, m_amp, m_freq, m_R, m_inside);
        return static_cast<BaseIF*>(SimpleCosSurfPtr);
      }

  protected:
    Real              m_R;
    vector<Real>      m_amp;
    vector<RealVect>  m_freq;
    vector<RealVect>  m_center;
    Real      m_PI=3.141592653589793238462643383279502884197;
    bool      m_inside;

  private:
    SimpleCosSurfIF();

  };
}
#endif
