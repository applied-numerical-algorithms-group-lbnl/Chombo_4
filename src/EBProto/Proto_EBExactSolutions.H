
#include <iostream>
using std::cerr;

#ifndef __EBEExactSolutions__
#define __EBEExactSolutions__

template <unsigned int order>
class BaseExactSolution
{
public:
  typedef       IndexTM<int, DIM>            IvDIM;
  typedef HostIrregData<CELL,IndMomDIM , 1>  VoluData;

  
  virtual ~BaseExactSolution()
  {
  }
  
  BaseExactSolution()
  {
  }

  ///
  Real operator()(const EBGraph & a_graph,
                  const Real    & a_dx,
                  const VoluData& a_voldat,
                  const EBIndex<CELL>   & a_vof)
  {
    Real retval = 0;
    RealVect cornerloc;
    Poitn cell = a_vof.m_pt;
    for(int idir = 0; idir < DIM; idir++)
    {
      cornerloc[idir] = a_dx*cell[idir];
    }
    IndexedMoments<DIM, order> moments;
    if(a_graph.isRegular(a_cell))
    {
      moments.setToRegular(a_dx);
    }
    else if(a_graph.isCovered(a_cell))
    {
      return 0;
    }
    else
    {
      Point pt = ProtoCh::getPoint(a_cell);
      moments = a_voldata(a_vof, 0);
    }
      
    Real volu = moments[IvDIM::Zero]
    for(MomentIterator<DIM, order> momit; momit.ok(); ++momit)
    {
      IvDIM  power = momit();
      Real   deriv = getDerivative(power, cornerloc);
      Real   denom = (Real)(deriv.factorial());
      Real   momen = moments[power];
      Real   incre = deriv*momen/denom;
      retval += incre;
    }
    Real minvol = 1.0e-10;
    if(volu > minvol)
    {
      retval /= volu;
    }
    return retval;   
  }

  virtual Real getDerivative(IvDIM    a_power,
                             RealVect a_loc) = 0;
};

///class for always returning zero (for bcs as opposed to being for testing convergence rates)
template <unsigned int order>
class ZeroSolution: public BaseSolution<order>
{
public:
  
  typedef       IndexTM<int, DIM>            IvDIM;
  typedef HostIrregData<CELL,IndMomDIM , 1>  VoluData;
  
  virtual Real getDerivative(IvDIM    a_power,
                             RealVect a_loc)
  {
    return 0.;
  }

  virtual ~ZeroSolution()
  {
  }
  
  ZeroSolution()
  {
  }

};

///phi = sin(pi(r^2 - R0^2)), where r is the distance from a center (c, c, c) and R0 is a radius
template <unsigned int order>
class SineSphere: public BaseSolution<order>
{
public:
  
  typedef       IndexTM<int, DIM>            IvDIM;
  typedef HostIrregData<CELL,IndMomDIM , 1>  VoluData;
  
  virtual Real getDerivative(IvDIM    a_power,
                             RealVect a_loc)
  {
    Real x = a_loc[0];
    Real y = a_loc[1];
    Real z = 0;
    Real pi = 4.*atan(1.0);
    Real r02 = m_R0*m_R0;
#if DIM==3
    z = a_loc[0];
#endif
    Real r2 = x*x + y*y + z*z;

    int sum = a_power.sum();
    Real retval =0;
    if(sum == 0)
    {
      //(%i2) display2d:false
      //(%o2) false
      //(%i3) r:x^2+y^2+z^2
      //(%o3) z^2+y^2+x^2
      //(%i4) phi:sin(pi*(r-R0))
      //(%o4) sin(pi*(z^2+y^2+x^2-R0))
      retval = sin(pi*(r2 - r02));
    }
    else if(sum == 1)
    {
      if(a_power[0]== 1)
      {
        //(%i5) der_x:diff(phi,x)
        //(%o5) 2*pi*x*cos(pi*(z^2+y^2+x^2-R0))
        retval = 2*pi*x*cos(pi*(r2 - r02));
      }
      else if (a_power[0]== 1)
      {
        //(%i6) der_y:diff(phi,y)
        //(%o6) 2*pi*y*cos(pi*(z^2+y^2+x^2-R0))
        retval = 2*pi*y*cos(pi*(r2 - r02));
      }
      else if(a_power[2] == 1)
      {
        //(%i7) der_z:diff(phi,z)
        //(%o7) 2*pi*z*cos(pi*(z^2+y^2+x^2-R0))
        retval = 2*pi*z*cos(pi*(r2 - r02));
      }
    }
    else if(sum == 2)
    {
      if((a_power[0]== 1) && (a_power[1]== 1))
      {
        //(%i8) der_xy:diff(der_x,y)
        //(%o8) -4*pi^2*x*y*sin(pi*(z^2+y^2+x^2-R0))
        retval = -4.*pi*pi*x*y*sin(pi*(r2 -r02));
      }
#if DIM==3      
      else if((a_power[0]== 1) && (a_power[2]== 1))
      {
        //(%i9) der_xz:diff(der_x,z)
        //(%o9) -4*pi^2*x*z*sin(pi*(z^2+y^2+x^2-R0))
        retval = -4.*pi*pi*x*z*sin(pi*(r2 -r02));
        
      }
      else if((a_power[1]== 1) && (a_power[2]== 1))
      {
        //(%i10) der_yz:diff(der_y,z)
        //(%o10) -4*pi^2*y*z*sin(pi*(z^2+y^2+x^2-R0))
        retval = -4.*pi*pi*y*z*sin(pi*(r2 -r02));
      }
      else if(a_power[2]== 2)
      {
        //(%i13) der_zz:diff(der_z,z)
        //(%o13) 2*pi*cos(pi*(z^2+y^2+x^2-R0))-4*pi^2*z^2*sin(pi*(z^2+y^2+x^2-R0))
        retval = 2.*pi*cos(pi*(r2 -r02))      -4.*pi*pi*z*z*sin(pi*(r2 - r02));

      }
#endif
      else if(a_power[0]== 2)
      {
        //(%i11) der_xx:diff(der_x,x)
        //(%o11)  2*pi*cos(pi*(z^2+y^2+x^2-R0))-4*pi^2*x^2*sin(pi*(z^2+y^2+x^2-R0))
        retval = 2.*pi*cos(pi*(r2 -r02))      -4.*pi*pi*x*x*sin(pi*(r2 - r02));
      }
      else if(a_power[1]== 2)
      {
        //(%i12) der_yy:diff(der_y,y)
        //(%o12) 2*pi*cos(pi*(z^2+y^2+x^2-R0))-4*pi^2*y^2*sin(pi*(z^2+y^2+x^2-R0))
        retval = 2.*pi*cos(pi*(r2 -r02))      -4.*pi*pi*y*y*sin(pi*(r2 - r02));

      }
      else
      {
        MayDay::Error("sinesphere ran out of cases 2");
      }
    }
    else if(sum == 3)
    {
      if(a_power[0] == 3)
      {
        //(%i14) der_xxx:diff(der_xx,x)
        //(%o14) (-12*pi^2*x*sin(pi*(z^2+y^2+x^2-R0)))      -8*pi^3*x^3*cos(pi*(z^2+y^2+x^2-R0))
        retval = (-12.*pi*pi*x*sin(pi*(r2 -r02))) -   8.*pi*pi*pi*x*x*x*cos(pi*(r2 - r02));
      }
      else if(a_power[1] == 3)
      {
        //(%i15) der_yyy:diff(der_yy,y)
        //(%o15) (-12*pi^2*y*sin(pi*(z^2+y^2+x^2-R0)))-8*pi^3*y^3*cos(pi*(z^2+y^2+x^2-R0))
        retval = (-12.*pi*pi*y*sin(pi*(r2 -r02))) -    8.*pi*pi*pi*y*y*y*cos(pi*(r2 - r02));
      }
      else if( (a_power[0] == 2) && (a_power[1]==1))
      {
        //(%i19) der_yxx:diff(der_xy,x)
        //(%o19) (-4*pi^2*y*sin(pi*(z^2+y^2+x^2-R0)))-8*pi^3*x^2*y*cos(pi*(z^2+y^2+x^2-R0))
        retval = (-4.*pi*pi*y*sin(pi*(r2 -r02))) -   8.*pi*pi*pi*x*x*y*cos(pi*(r2 - r02));
      }
      else if( (a_power[0] == 1) && (a_power[1]==2))
      {
        //(%i18) der_xyy:diff(der_xy,y)
        //(%o18) (-4*pi^2*x*sin(pi*(z^2+y^2+x^2-R0)))-8*pi^3*x*y^2*cos(pi*(z^2+y^2+x^2-R0))
        retval = (-4.*pi*pi*x*sin(pi*(r2 -r02))) -   8.*pi*pi*pi*x*y*y*cos(pi*(r2 - r02));

      }
#if DIM==3
      else if(a_power[2] == 3)
      {
        //(%i17) der_zzz:diff(der_zz,z)
        //(%o17) (-12*pi^2*z*sin(pi*(z^2+y^2+x^2-R0)))-8*pi^3*z^3*cos(pi*(z^2+y^2+x^2-R0))
        retval = (-12.*pi*pi*z*sin(pi*(r2 -r02))) -    8.*pi*pi*pi*z*z*z*cos(pi*(r2 - r02));
      }
      else if( (a_power[2] == 2) && (a_power[1]==1))
      {
        //(%i20) der_yzz:diff(der_yz,z)
        //(%o20) (-4*pi^2*y*sin(pi*(z^2+y^2+x^2-R0)))-8*pi^3*y*z^2*cos(pi*(z^2+y^2+x^2-R0))
        retval = (-4.*pi*pi*y*sin(pi*(r2 -r02))) -   8.*pi*pi*pi*y*z*z*cos(pi*(r2 - r02));
      }
      else if( (a_power[2] == 2) && (a_power[0]==1))
      {
        //(%i21) der_xzz:diff(der_xz,z)
        //(%o21) (-4*pi^2*x*sin(pi*(z^2+y^2+x^2-R0)))-8*pi^3*x*z^2*cos(pi*(z^2+y^2+x^2-R0))
        retval = (-4.*pi*pi*x*sin(pi*(r2 -r02))) -    8.*pi*pi*pi*x*z*z*cos(pi*(r2 - r02));
      }
      else if( (a_power[0] == 2) && (a_power[2]==1))
      {
        //  (%i23) der_xxz:diff(der_xx, z);
        //  (%o23) (-4*pi^2*z*sin(pi*(z^2+y^2+x^2-R0)))-8*pi^3*x^2*z*cos(pi*(z^2+y^2+x^2-R0))
        retval = (-4.*pi*pi*z*sin(pi*(r2 -r02))) -    8.*pi*pi*pi*x*x*z*cos(pi*(r2 - r02));
      }
      else if( (a_power[1] == 2) && (a_power[2]==1))
      {
        //(%i24) der_yyz:diff(der_yy, z);
        //(%o24) (-4*pi^2*z*sin(pi*(z^2+y^2+x^2-R0)))-8*pi^3*y^2*z*cos(pi*(z^2+y^2+x^2-R0))
        retval = (-4.*pi*pi*z*sin(pi*(r2 -r02))) -    8.*pi*pi*pi*y*y*z*cos(pi*(r2 - r02));
      }
#endif      
      else
      {
        MayDay::Warning("sinesphere ran out of cases 3");
        retval = 0;
      }
    }
    else
    {
      MayDay::Warning("sinesphere ran out of cases");
      retval = 0;
    }
        
    return retval;
  }

  virtual ~SineSphere()
  {
  }

  SineSphere(Real a_radius,
             Real a_center)
  {
    m_radius = a_radius;
    m_center = a_center;
  }

private:
  Real  m_radius;
  Real  m_center;
  
  SineSphere()
  {
  }

};


#endif
