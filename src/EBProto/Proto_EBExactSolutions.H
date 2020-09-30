
#include <iostream>
using std::cerr;

#ifndef __EBEExactSolutions__
#define __EBEExactSolutions__

template <unsigned int order>
class BaseExactSolution
{
public:
  typedef IndexedMoments<DIM  , order> IndMomDIM;
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
    RealVect centloc;
    Point cell = a_vof.m_pt;
    for(int idir = 0; idir < DIM; idir++)
    {
      centloc[idir] = a_dx*(Real(cell[idir]) + 0.5);
    }
    IndexedMoments<DIM, order> moments;
    if(a_graph.isRegular(cell))
    {
      moments.setToRegular(a_dx);
    }
    else if(a_graph.isCovered(cell))
    {
      return 0;
    }
    else
    {
      moments = a_voldat(a_vof, 0);
    }
      
    Real volu = moments[IvDIM::Zero];
    for(MomentIterator<DIM, order> momit; momit.ok(); ++momit)
    {
      IvDIM  power = momit();
      Real   deriv = getDerivative(power, centloc);
      int   idenom = (power.factorial());
      Real   denom = (Real)(idenom);
      Real   momen = moments[power];
      Real   incre = deriv*momen/denom;
      retval += incre;
    }
    Real minvol = 1.0e-16;
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
class ZeroSolution: public BaseExactSolution<order>
{
public:
  
  typedef       IndexTM<int, DIM>            IvDIM;
  typedef IndexedMoments<DIM  , order> IndMomDIM;
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
class SineSphere: public BaseExactSolution<order>
{
public:
  
  typedef       IndexTM<int, DIM>            IvDIM;
  typedef IndexedMoments<DIM  , order> IndMomDIM;
  typedef HostIrregData<CELL,IndMomDIM , 1>  VoluData;
  
  virtual Real getDerivative(IvDIM    a_power,
                             RealVect a_loc)
  {
    Real x = a_loc[0];
    Real y = a_loc[1];
    Real x2 = x*x;
    Real x3 = x*x*x;
    Real x4 = x*x*x*x;
    Real y2 = y*y;
    Real y3 = y*y*y;
    Real y4 = y*y*y*y;

    Real pi = 4.*atan(1.0);
    Real pi2 = pi*pi;
    Real pi3 = pi*pi*pi;
    Real pi4 = pi*pi*pi*pi;
    Real r02 = m_radius*m_radius;
    Real z2 = 0;
#if DIM==3
    Real z = a_loc[2];
    z2 = z*z;
    Real z3 = z*z*z;
    Real z4 = z*z*z*z;
#endif
    Real r2 = x2 + y2 + z2;

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
      else if (a_power[1]== 1)
      {
        //(%i6) der_y:diff(phi,y)
        //(%o6) 2*pi*y*cos(pi*(z^2+y^2+x^2-R0))
        retval = 2*pi*y*cos(pi*(r2 - r02));
      }
#if DIM==3      
      else if(a_power[2] == 1)
      {
        //(%i7) der_z:diff(phi,z)
        //(%o7) 2*pi*z*cos(pi*(z^2+y^2+x^2-R0))
        retval = 2*pi*z*cos(pi*(r2 - r02));
      }
#endif      
      else
      {
        MayDay::Error("sinesphere ran out of cases 2");
      }
    }
    else if(sum == 2)
    {
      if((a_power[0]== 1) && (a_power[1]== 1))
      {
        //(%i8) der_xy:diff(der_x,y)
        //(%o8) -4*pi^2*x*y*sin(pi*(z^2+y^2+x^2-R0))
        retval = -4.*pi2*x*y*sin(pi*(r2 -r02));
      }
#if DIM==3      
      else if((a_power[0]== 1) && (a_power[2]== 1))
      {
        //(%i9) der_xz:diff(der_x,z)
        //(%o9) -4*pi^2*x*z*sin(pi*(z^2+y^2+x^2-R0))
        retval = -4.*pi2*x*z*sin(pi*(r2 -r02));
        
      }
      else if((a_power[1]== 1) && (a_power[2]== 1))
      {
        //(%i10) der_yz:diff(der_y,z)
        //(%o10) -4*pi^2*y*z*sin(pi*(z^2+y^2+x^2-R0))
        retval = -4.*pi2*y*z*sin(pi*(r2 -r02));
      }
      else if(a_power[2]== 2)
      {
        //(%i13) der_zz:diff(der_z,z)
        //(%o13) 2*pi*cos(pi*(z^2+y^2+x^2-R0))-4*pi^2*z^2*sin(pi*(z^2+y^2+x^2-R0))
        retval = 2.*pi*cos(pi*(r2 -r02))      -4.*pi2*z2*sin(pi*(r2 - r02));

      }
#endif
      else if(a_power[0]== 2)
      {
        //(%i11) der_xx:diff(der_x,x)
        //(%o11)  2*pi*cos(pi*(z^2+y^2+x^2-R0))-4*pi^2*x^2*sin(pi*(z^2+y^2+x^2-R0))
        retval = 2.*pi*cos(pi*(r2 -r02))      -4.*pi2*x2*sin(pi*(r2 - r02));
      }
      else if(a_power[1]== 2)
      {
        //(%i12) der_yy:diff(der_y,y)
        //(%o12) 2*pi*cos(pi*(z^2+y^2+x^2-R0))-4*pi^2*y^2*sin(pi*(z^2+y^2+x^2-R0))
        retval = 2.*pi*cos(pi*(r2 -r02))      -4.*pi2*y2*sin(pi*(r2 - r02));

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
        retval = (-12.*pi2*x*sin(pi*(r2 -r02))) -   8.*pi3*x3*cos(pi*(r2 - r02));
      }
      else if(a_power[1] == 3)
      {
        //(%i15) der_yyy:diff(der_yy,y)
        //(%o15) (-12*pi^2*y*sin(pi*(z^2+y^2+x^2-R0)))-8*pi^3*y^3*cos(pi*(z^2+y^2+x^2-R0))
        retval = (-12.*pi2*y*sin(pi*(r2 -r02))) -    8.*pi3*y3*cos(pi*(r2 - r02));
      }
      else if( (a_power[0] == 2) && (a_power[1]==1))
      {
        //(%i19) der_yxx:diff(der_xy,x)
        //(%o19) (-4*pi^2*y*sin(pi*(z^2+y^2+x^2-R0)))-8*pi^3*x^2*y*cos(pi*(z^2+y^2+x^2-R0))
        retval = (-4.*pi2*y*sin(pi*(r2 -r02))) -   8.*pi3*x2*y*cos(pi*(r2 - r02));
      }
      else if( (a_power[0] == 1) && (a_power[1]==2))
      {
        //(%i18) der_xyy:diff(der_xy,y)
        //(%o18) (-4*pi^2*x*sin(pi*(z^2+y^2+x^2-R0)))-8*pi^3*x*y^2*cos(pi*(z^2+y^2+x^2-R0))
        retval = (-4.*pi2*x*sin(pi*(r2 -r02))) -   8.*pi3*x*y2*cos(pi*(r2 - r02));

      }
#if DIM==3
      else if(a_power[2] == 3)
      {
        //(%i17) der_zzz:diff(der_zz,z)
        //(%o17) (-12*pi^2*z*sin(pi*(z^2+y^2+x^2-R0)))-8*pi^3*z^3*cos(pi*(z^2+y^2+x^2-R0))
        retval = (-12.*pi2*z*sin(pi*(r2 -r02))) -    8.*pi3*z3*cos(pi*(r2 - r02));
      }
      else if( (a_power[2] == 2) && (a_power[1]==1))
      {
        //(%i20) der_yzz:diff(der_yz,z)
        //(%o20) (-4*pi^2*y*sin(pi*(z^2+y^2+x^2-R0)))-8*pi^3*y*z^2*cos(pi*(z^2+y^2+x^2-R0))
        retval = (-4.*pi2*y*sin(pi*(r2 -r02))) -   8.*pi3*y*z2*cos(pi*(r2 - r02));
      }
      else if( (a_power[2] == 2) && (a_power[0]==1))
      {
        //(%i21) der_xzz:diff(der_xz,z)
        //(%o21) (-4*pi^2*x*sin(pi*(z^2+y^2+x^2-R0)))-8*pi^3*x*z^2*cos(pi*(z^2+y^2+x^2-R0))
        retval = (-4.*pi2*x*sin(pi*(r2 -r02))) -    8.*pi3*x*z2*cos(pi*(r2 - r02));
      }
      else if( (a_power[0] == 2) && (a_power[2]==1))
      {
        //  (%i23) der_xxz:diff(der_xx, z);
        //  (%o23) (-4*pi^2*z*sin(pi*(z^2+y^2+x^2-R0)))-8*pi^3*x^2*z*cos(pi*(z^2+y^2+x^2-R0))
        retval = (-4.*pi2*z*sin(pi*(r2 -r02))) -    8.*pi3*x2*z*cos(pi*(r2 - r02));
      }
      else if( (a_power[1] == 2) && (a_power[2]==1))
      {
        //(%i24) der_yyz:diff(der_yy, z);
        //(%o24) (-4*pi^2*z*sin(pi*(z^2+y^2+x^2-R0)))-8*pi^3*y^2*z*cos(pi*(z^2+y^2+x^2-R0))
        retval = (-4.*pi2*z*sin(pi*(r2 -r02))) -    8.*pi3*y2*z*cos(pi*(r2 - r02));
      }
#endif      
      else
      {
        MayDay::Error("sinesphere ran out of cases 3");
        retval = 0;
      }
    }
    else if(sum==4)
    {
      if(a_power[0] == 4)
      {
        //(%i22) der_xxxx:diff(der_xxx,x)
        //(%o22) 16*pi^4*x^4*sin(pi*(z^2+y^2+x^2-R0))   -12*pi^2*sin(pi*(z^2+y^2+x^2-R0))-48*pi^3*x^2*cos(pi*(z^2+y^2+x^2-R0))
        retval = 16.*pi4*x4*sin(pi*(r2-r02)) -12.*pi2*sin(pi*(r2 -r02))     -48.*pi3*x2*cos(pi*(r2 - r02));
      }
      else if(a_power[1] == 4)
      {
        //(%i23) der_yyyy:diff(der_yyy,y)
        //(%o23) 16*pi^4*y^4*sin(pi*(z^2+y^2+x^2-R0))   -12*pi^2*sin(pi*(z^2+y^2+x^2-R0))-48*pi^3*y^2*cos(pi*(z^2+y^2+x^2-R0))
        retval = 16.*pi4*y4*sin(pi*(r2-r02)) -12.*pi2*sin(pi*(r2 -r02))     -48.*pi3*y2*cos(pi*(r2 - r02));
      }
      else if( (a_power[0] == 3) && (a_power[1] == 1))
      {
        //(%i25) der_xxxy:diff(der_xxx,y)
        //(%o25) 16*pi^4*x^3*y*sin(pi*(z^2+y^2+x^2-R0))-24*pi^3*x*y*cos(pi*(z^2+y^2+x^2-R0))
        retval = 16.*pi4*x3*y*sin(pi*(r2-r02))          -24.*pi3*x*y*cos(pi*(r2 -r02))     ;
      }
      else if( (a_power[0] == 1) && (a_power[1] == 3))
      {
        //(%i27) der_yyyx:diff(der_yyy,x)
        //(%o27) 16*pi^4*x*y^3*sin(pi*(z^2+y^2+x^2-R0))-24*pi^3*x*y*cos(pi*(z^2+y^2+x^2-R0))
        retval = 16.*pi4*y3*x*sin(pi*(r2-r02))          -24.*pi3*x*y*cos(pi*(r2 -r02))     ;
      }
      else if( (a_power[0] == 2) && (a_power[1] == 2))
      {
        //(%i33) der_xxyy:diff(der_xyy,x)
        //(%o33) 16*pi^4*x^2*y^2*sin(pi*(z^2+y^2+x^2-R0))-4*pi^2*sin(pi*(z^2+y^2+x^2-R0))-8*pi^3*y^2*cos(pi*(z^2+y^2+x^2-R0))-8*pi^3*x^2*cos(pi*(z^2+y^2+x^2-R0))
        retval = 16.*pi4*x2*y2*sin(pi*(r2-r02))           -4.*pi2*sin(pi*(r2 -r02))       -8.*pi3*y2*cos(pi*(r2 - r02))       -8.*pi3*x2*cos(pi*(r2 - r02));
      }
#if DIM==3
      else if(a_power[2] == 4)
      {
        //(%i24) der_zzzz:diff(der_zzz,z)
        //(%o24) 16*pi^4*z^4*sin(pi*(z^2+y^2+x^2-R0))   -12*pi^2*sin(pi*(z^2+y^2+x^2-R0))-48*pi^3*z^2*cos(pi*(z^2+y^2+x^2-R0))
        retval = 16.*pi4*z4*sin(pi(r2-r02)) -12.*pi2*sin(pi*(r2 -r02))     -48.*pi3*z2*cos(pi*(r2 - r02));
      }
      else if((a_power[2] == 3) && (a_power[0] == 1))
      {
        //(%i29) der_zzzx:diff(der_zzz,x)
        //(%o29) 16*pi^4*x*z^3*sin(pi*(z^2+y^2+x^2-R0))-24*pi^3*x*z*cos(pi*(z^2+y^2+x^2-R0))
        retval = 16.*pi4*x*z3*sin(pi*(r2-r02))  -24.*pi3*x*z*cos(pi*(r2 -r02));
      }
      else if((a_power[2] == 3) && (a_power[1] == 1))
      {
        //(%i30) der_zzzy:diff(der_zzz,y)
        //(%o30) 16*pi^4*y*z^3*sin(pi*(z^2+y^2+x^2-R0))-24*pi^3*y*z*cos(pi*(z^2+y^2+x^2-R0))
        retval = 16.*pi4*y*z3*sin(pi*(r2-r02))  -24.*pi3*x*z*cos(pi*(r2 -r02));

      }
      else if((a_power[2] == 2) && (a_power[0] == 2))
      {
        //(%i32) der_zzxx:diff(der_xzz,x)
        //(%o32) 16*pi^4*x^2*z^2*sin(pi*(z^2+y^2+x^2-R0))-4*pi^2*sin(pi*(z^2+y^2+x^2-R0))-8*pi^3*z^2*cos(pi*(z^2+y^2+x^2-R0))-8*pi^3*x^2*cos(pi*(z^2+y^2+x^2-R0))
        retval = 16.*pi4*x2*z2*sin(pi*(r2-r02))           -4.*pi2*sin(pi*(r2 -r02))       -8.*pi3*z2*cos(pi*(r2 - r02))       -8.*pi3*x2*cos(pi*(r2 - r02));
      }
      else if((a_power[2] == 2) && (a_power[1] == 2))
      {
        //(%i31) der_zzyy:diff(der_yzz,y)
        //(%o31) 16*pi^4*y^2*z^2*sin(pi*(z^2+y^2+x^2-R0))-4*pi^2*sin(pi*(z^2+y^2+x^2-R0))-8*pi^3*z^2*cos(pi*(z^2+y^2+x^2-R0))-8*pi^3*y^2*cos(pi*(z^2+y^2+x^2-R0))
        retval = 16.*pi4*y2*z2*sin(pi*(r2-r02))           -4.*pi2*sin(pi*(r2 -r02))       -8.*pi3*z2*cos(pi*(r2 - r02))       -8.*pi3*y2*cos(pi*(r2 - r02));
      }
      else if((a_power[2] == 1) && (a_power[0] == 3))
      {
        //(%i26) der_xxxz:diff(der_xxx,z)
        //(%o26) 16*pi^4*x^3*z*sin(pi*(z^2+y^2+x^2-R0))-24*pi^3*x*z*cos(pi*(z^2+y^2+x^2-R0))
        retval = 16.*pi4*x3*z*sin(pi*(r2-r02))          -24.*pi3*x*z*sin(pi*(r2 -r02))      ;
        
      }
      else if((a_power[2] == 1) && (a_power[1] == 3))
      {
        //(%i28) der_yyyz:diff(der_yyy,z)
        //(%o28) 16*pi^4*y^3*z*sin(pi*(z^2+y^2+x^2-R0))-24*pi^3*y*z*cos(pi*(z^2+y^2+x^2-R0))
        retval = 16.*pi4*y3*z*sin(pi*(r2-r02))          -24.*pi3*y*z*sin(pi*(r2 -r02))      ;
      }
#endif      
      else
      {
        MayDay::Error("sinesphere ran out of cases 4");
        retval = 0;
      }
    }
    else
    {
      MayDay::Error("sinesphere ran out of power cases");
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


///phi = cp x^p  
template <unsigned int order>
class MonomialEF: public BaseExactSolution<order>
{
public:
  
  typedef       IndexTM<int, DIM>            IvDIM;
  typedef IndexedMoments<DIM  , order>       IndMomDIM;
  typedef HostIrregData<CELL,IndMomDIM , 1>  VoluData;
  
  virtual Real getDerivative(IvDIM    a_power,
                             RealVect a_loc)
  {
  {
    Point power = m_entry.second;
    Real  C     = m_entry.first;
    
    Real x = a_loc[0];
    Real y = a_loc[1];
    int px = power[0];
    int py = power[1];
    Real x2 = x*x;
    Real x3 = x*x*x;
    Real x4 = x*x*x*x;
    Real y2 = y*y;
    Real y3 = y*y*y;
    Real y4 = y*y*y*y;

    Real pi = 4.*atan(1.0);
    Real pi2 = pi*pi;
    Real pi3 = pi*pi*pi;
    Real pi4 = pi*pi*pi*pi;
    Real r02 = m_radius*m_radius;
    Real z2 = 0;
#if DIM==3
    Real z = a_loc[2];
    z2 = z*z;
    Real z3 = z*z*z;
    Real z4 = z*z*z*z;
#endif
    Real r2 = x2 + y2 + z2;

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
      else if (a_power[1]== 1)
      {
        //(%i6) der_y:diff(phi,y)
        //(%o6) 2*pi*y*cos(pi*(z^2+y^2+x^2-R0))
        retval = 2*pi*y*cos(pi*(r2 - r02));
      }
#if DIM==3      
      else if(a_power[2] == 1)
      {
        //(%i7) der_z:diff(phi,z)
        //(%o7) 2*pi*z*cos(pi*(z^2+y^2+x^2-R0))
        retval = 2*pi*z*cos(pi*(r2 - r02));
      }
#endif      
      else
      {
        MayDay::Error("sinesphere ran out of cases 2");
      }
    }
    else if(sum == 2)
    {
      if((a_power[0]== 1) && (a_power[1]== 1))
      {
        //(%i8) der_xy:diff(der_x,y)
        //(%o8) -4*pi^2*x*y*sin(pi*(z^2+y^2+x^2-R0))
        retval = -4.*pi2*x*y*sin(pi*(r2 -r02));
      }
#if DIM==3      
      else if((a_power[0]== 1) && (a_power[2]== 1))
      {
        //(%i9) der_xz:diff(der_x,z)
        //(%o9) -4*pi^2*x*z*sin(pi*(z^2+y^2+x^2-R0))
        retval = -4.*pi2*x*z*sin(pi*(r2 -r02));
        
      }
      else if((a_power[1]== 1) && (a_power[2]== 1))
      {
        //(%i10) der_yz:diff(der_y,z)
        //(%o10) -4*pi^2*y*z*sin(pi*(z^2+y^2+x^2-R0))
        retval = -4.*pi2*y*z*sin(pi*(r2 -r02));
      }
      else if(a_power[2]== 2)
      {
        //(%i13) der_zz:diff(der_z,z)
        //(%o13) 2*pi*cos(pi*(z^2+y^2+x^2-R0))-4*pi^2*z^2*sin(pi*(z^2+y^2+x^2-R0))
        retval = 2.*pi*cos(pi*(r2 -r02))      -4.*pi2*z2*sin(pi*(r2 - r02));

      }
#endif
      else if(a_power[0]== 2)
      {
        //(%i11) der_xx:diff(der_x,x)
        //(%o11)  2*pi*cos(pi*(z^2+y^2+x^2-R0))-4*pi^2*x^2*sin(pi*(z^2+y^2+x^2-R0))
        retval = 2.*pi*cos(pi*(r2 -r02))      -4.*pi2*x2*sin(pi*(r2 - r02));
      }
      else if(a_power[1]== 2)
      {
        //(%i12) der_yy:diff(der_y,y)
        //(%o12) 2*pi*cos(pi*(z^2+y^2+x^2-R0))-4*pi^2*y^2*sin(pi*(z^2+y^2+x^2-R0))
        retval = 2.*pi*cos(pi*(r2 -r02))      -4.*pi2*y2*sin(pi*(r2 - r02));

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
        retval = (-12.*pi2*x*sin(pi*(r2 -r02))) -   8.*pi3*x3*cos(pi*(r2 - r02));
      }
      else if(a_power[1] == 3)
      {
        //(%i15) der_yyy:diff(der_yy,y)
        //(%o15) (-12*pi^2*y*sin(pi*(z^2+y^2+x^2-R0)))-8*pi^3*y^3*cos(pi*(z^2+y^2+x^2-R0))
        retval = (-12.*pi2*y*sin(pi*(r2 -r02))) -    8.*pi3*y3*cos(pi*(r2 - r02));
      }
      else if( (a_power[0] == 2) && (a_power[1]==1))
      {
        //(%i19) der_yxx:diff(der_xy,x)
        //(%o19) (-4*pi^2*y*sin(pi*(z^2+y^2+x^2-R0)))-8*pi^3*x^2*y*cos(pi*(z^2+y^2+x^2-R0))
        retval = (-4.*pi2*y*sin(pi*(r2 -r02))) -   8.*pi3*x2*y*cos(pi*(r2 - r02));
      }
      else if( (a_power[0] == 1) && (a_power[1]==2))
      {
        //(%i18) der_xyy:diff(der_xy,y)
        //(%o18) (-4*pi^2*x*sin(pi*(z^2+y^2+x^2-R0)))-8*pi^3*x*y^2*cos(pi*(z^2+y^2+x^2-R0))
        retval = (-4.*pi2*x*sin(pi*(r2 -r02))) -   8.*pi3*x*y2*cos(pi*(r2 - r02));

      }
#if DIM==3
      else if(a_power[2] == 3)
      {
        //(%i17) der_zzz:diff(der_zz,z)
        //(%o17) (-12*pi^2*z*sin(pi*(z^2+y^2+x^2-R0)))-8*pi^3*z^3*cos(pi*(z^2+y^2+x^2-R0))
        retval = (-12.*pi2*z*sin(pi*(r2 -r02))) -    8.*pi3*z3*cos(pi*(r2 - r02));
      }
      else if( (a_power[2] == 2) && (a_power[1]==1))
      {
        //(%i20) der_yzz:diff(der_yz,z)
        //(%o20) (-4*pi^2*y*sin(pi*(z^2+y^2+x^2-R0)))-8*pi^3*y*z^2*cos(pi*(z^2+y^2+x^2-R0))
        retval = (-4.*pi2*y*sin(pi*(r2 -r02))) -   8.*pi3*y*z2*cos(pi*(r2 - r02));
      }
      else if( (a_power[2] == 2) && (a_power[0]==1))
      {
        //(%i21) der_xzz:diff(der_xz,z)
        //(%o21) (-4*pi^2*x*sin(pi*(z^2+y^2+x^2-R0)))-8*pi^3*x*z^2*cos(pi*(z^2+y^2+x^2-R0))
        retval = (-4.*pi2*x*sin(pi*(r2 -r02))) -    8.*pi3*x*z2*cos(pi*(r2 - r02));
      }
      else if( (a_power[0] == 2) && (a_power[2]==1))
      {
        //  (%i23) der_xxz:diff(der_xx, z);
        //  (%o23) (-4*pi^2*z*sin(pi*(z^2+y^2+x^2-R0)))-8*pi^3*x^2*z*cos(pi*(z^2+y^2+x^2-R0))
        retval = (-4.*pi2*z*sin(pi*(r2 -r02))) -    8.*pi3*x2*z*cos(pi*(r2 - r02));
      }
      else if( (a_power[1] == 2) && (a_power[2]==1))
      {
        //(%i24) der_yyz:diff(der_yy, z);
        //(%o24) (-4*pi^2*z*sin(pi*(z^2+y^2+x^2-R0)))-8*pi^3*y^2*z*cos(pi*(z^2+y^2+x^2-R0))
        retval = (-4.*pi2*z*sin(pi*(r2 -r02))) -    8.*pi3*y2*z*cos(pi*(r2 - r02));
      }
#endif      
      else
      {
        MayDay::Error("sinesphere ran out of cases 3");
        retval = 0;
      }
    }
    else if(sum==4)
    {
      if(a_power[0] == 4)
      {
        //(%i22) der_xxxx:diff(der_xxx,x)
        //(%o22) 16*pi^4*x^4*sin(pi*(z^2+y^2+x^2-R0))   -12*pi^2*sin(pi*(z^2+y^2+x^2-R0))-48*pi^3*x^2*cos(pi*(z^2+y^2+x^2-R0))
        retval = 16.*pi4*x4*sin(pi*(r2-r02)) -12.*pi2*sin(pi*(r2 -r02))     -48.*pi3*x2*cos(pi*(r2 - r02));
      }
      else if(a_power[1] == 4)
      {
        //(%i23) der_yyyy:diff(der_yyy,y)
        //(%o23) 16*pi^4*y^4*sin(pi*(z^2+y^2+x^2-R0))   -12*pi^2*sin(pi*(z^2+y^2+x^2-R0))-48*pi^3*y^2*cos(pi*(z^2+y^2+x^2-R0))
        retval = 16.*pi4*y4*sin(pi*(r2-r02)) -12.*pi2*sin(pi*(r2 -r02))     -48.*pi3*y2*cos(pi*(r2 - r02));
      }
      else if( (a_power[0] == 3) && (a_power[1] == 1))
      {
        //(%i25) der_xxxy:diff(der_xxx,y)
        //(%o25) 16*pi^4*x^3*y*sin(pi*(z^2+y^2+x^2-R0))-24*pi^3*x*y*cos(pi*(z^2+y^2+x^2-R0))
        retval = 16.*pi4*x3*y*sin(pi*(r2-r02))          -24.*pi3*x*y*cos(pi*(r2 -r02))     ;
      }
      else if( (a_power[0] == 1) && (a_power[1] == 3))
      {
        //(%i27) der_yyyx:diff(der_yyy,x)
        //(%o27) 16*pi^4*x*y^3*sin(pi*(z^2+y^2+x^2-R0))-24*pi^3*x*y*cos(pi*(z^2+y^2+x^2-R0))
        retval = 16.*pi4*y3*x*sin(pi*(r2-r02))          -24.*pi3*x*y*cos(pi*(r2 -r02))     ;
      }
      else if( (a_power[0] == 2) && (a_power[1] == 2))
      {
        //(%i33) der_xxyy:diff(der_xyy,x)
        //(%o33) 16*pi^4*x^2*y^2*sin(pi*(z^2+y^2+x^2-R0))-4*pi^2*sin(pi*(z^2+y^2+x^2-R0))-8*pi^3*y^2*cos(pi*(z^2+y^2+x^2-R0))-8*pi^3*x^2*cos(pi*(z^2+y^2+x^2-R0))
        retval = 16.*pi4*x2*y2*sin(pi*(r2-r02))           -4.*pi2*sin(pi*(r2 -r02))       -8.*pi3*y2*cos(pi*(r2 - r02))       -8.*pi3*x2*cos(pi*(r2 - r02));
      }
#if DIM==3
      else if(a_power[2] == 4)
      {
        //(%i24) der_zzzz:diff(der_zzz,z)
        //(%o24) 16*pi^4*z^4*sin(pi*(z^2+y^2+x^2-R0))   -12*pi^2*sin(pi*(z^2+y^2+x^2-R0))-48*pi^3*z^2*cos(pi*(z^2+y^2+x^2-R0))
        retval = 16.*pi4*z4*sin(pi(r2-r02)) -12.*pi2*sin(pi*(r2 -r02))     -48.*pi3*z2*cos(pi*(r2 - r02));
      }
      else if((a_power[2] == 3) && (a_power[0] == 1))
      {
        //(%i29) der_zzzx:diff(der_zzz,x)
        //(%o29) 16*pi^4*x*z^3*sin(pi*(z^2+y^2+x^2-R0))-24*pi^3*x*z*cos(pi*(z^2+y^2+x^2-R0))
        retval = 16.*pi4*x*z3*sin(pi*(r2-r02))  -24.*pi3*x*z*cos(pi*(r2 -r02));
      }
      else if((a_power[2] == 3) && (a_power[1] == 1))
      {
        //(%i30) der_zzzy:diff(der_zzz,y)
        //(%o30) 16*pi^4*y*z^3*sin(pi*(z^2+y^2+x^2-R0))-24*pi^3*y*z*cos(pi*(z^2+y^2+x^2-R0))
        retval = 16.*pi4*y*z3*sin(pi*(r2-r02))  -24.*pi3*x*z*cos(pi*(r2 -r02));

      }
      else if((a_power[2] == 2) && (a_power[0] == 2))
      {
        //(%i32) der_zzxx:diff(der_xzz,x)
        //(%o32) 16*pi^4*x^2*z^2*sin(pi*(z^2+y^2+x^2-R0))-4*pi^2*sin(pi*(z^2+y^2+x^2-R0))-8*pi^3*z^2*cos(pi*(z^2+y^2+x^2-R0))-8*pi^3*x^2*cos(pi*(z^2+y^2+x^2-R0))
        retval = 16.*pi4*x2*z2*sin(pi*(r2-r02))           -4.*pi2*sin(pi*(r2 -r02))       -8.*pi3*z2*cos(pi*(r2 - r02))       -8.*pi3*x2*cos(pi*(r2 - r02));
      }
      else if((a_power[2] == 2) && (a_power[1] == 2))
      {
        //(%i31) der_zzyy:diff(der_yzz,y)
        //(%o31) 16*pi^4*y^2*z^2*sin(pi*(z^2+y^2+x^2-R0))-4*pi^2*sin(pi*(z^2+y^2+x^2-R0))-8*pi^3*z^2*cos(pi*(z^2+y^2+x^2-R0))-8*pi^3*y^2*cos(pi*(z^2+y^2+x^2-R0))
        retval = 16.*pi4*y2*z2*sin(pi*(r2-r02))           -4.*pi2*sin(pi*(r2 -r02))       -8.*pi3*z2*cos(pi*(r2 - r02))       -8.*pi3*y2*cos(pi*(r2 - r02));
      }
      else if((a_power[2] == 1) && (a_power[0] == 3))
      {
        //(%i26) der_xxxz:diff(der_xxx,z)
        //(%o26) 16*pi^4*x^3*z*sin(pi*(z^2+y^2+x^2-R0))-24*pi^3*x*z*cos(pi*(z^2+y^2+x^2-R0))
        retval = 16.*pi4*x3*z*sin(pi*(r2-r02))          -24.*pi3*x*z*sin(pi*(r2 -r02))      ;
        
      }
      else if((a_power[2] == 1) && (a_power[1] == 3))
      {
        //(%i28) der_yyyz:diff(der_yyy,z)
        //(%o28) 16*pi^4*y^3*z*sin(pi*(z^2+y^2+x^2-R0))-24*pi^3*y*z*cos(pi*(z^2+y^2+x^2-R0))
        retval = 16.*pi4*y3*z*sin(pi*(r2-r02))          -24.*pi3*y*z*sin(pi*(r2 -r02))      ;
      }
#endif      
      else
      {
        MayDay::Error("sinesphere ran out of cases 4");
        retval = 0;
      }
    }
    else
    {
      MayDay::Error("sinesphere ran out of power cases");
      retval = 0;
    }
        
    return retval;
  }
  virtual ~MonomialEF()
  {
  }

  MonomialEF(std::pair<Point, Real> > a_entry)
  {
    m_entry = a_entry;
  }

private:
  std::pair<Point, Real>  m_entry

  
  MonomialEF()
  {
  }
};

///phi = sum cp x^p  
template <unsigned int order>
class PolynomialEF: public BaseExactSolution<order>
{
public:
  
  typedef       IndexTM<int, DIM>            IvDIM;
  typedef IndexedMoments<DIM  , order>       IndMomDIM;
  typedef HostIrregData<CELL,IndMomDIM , 1>  VoluData;
  
  virtual Real getDerivative(IvDIM    a_power,
                             RealVect a_loc)
  {
    Real retval = 0;
    for(int ient = 0; ient < m_entries.size(); ient++)
    {
      Monomial mono(m_entries[ient]);
      Real monoval = mono.derivative(a_power, a_loc);
      retval += monoval;
    }
    return retval;

  virtual ~PolynomialEF()
  {
  }

  PolynomialEF(std::vector<std::pair<Point, Real> > a_entries)
  {
    m_entries = a_entries;
  }

private:
  std::vector<std::pair<Point, Real> > m_entries;

  
  PolynomialEF()
  {
  }
};


#endif
