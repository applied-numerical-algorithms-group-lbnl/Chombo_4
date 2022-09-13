#ifndef __SetupFunctionss__
#define __SetupFunctionss__
#define PI 3.141592653589793
#include "Chombo_NamespaceHeader.H"
using Proto::Var;
//=================================================
PROTO_KERNEL_START 
void InitializeVFaceF(int          a_p[DIM],
                      Var<Real, 1> a_vel,
                      Real         a_geomCen,
                      Real         a_geomRad,
                      Real         a_blobCen,
                      Real         a_blobRad,
                      Real         a_maxVelMag,
                      Real         a_maxVelRad,
                      Real         a_dx,
                      unsigned int a_dir)
{
  Real rlocsq = 0;
  Real xrel[DIM];
  for(int idir = 0; idir < DIM; idir++)
  {
    xrel[idir] = (a_p[idir] + 0.5)*a_dx - a_geomCen;
    rlocsq += xrel[idir]*xrel[idir];
  }

  Real velmag = sin(PI*rlocsq);
  Real vel[DIM];
#if DIM==2
  vel[0] =  velmag*xrel[1];
  vel[1] = -velmag*xrel[0];
#else
  vel[0] =  velmag*( xrel[1] + xrel[2]);
  vel[1] =  velmag*(-xrel[0] - xrel[2]);
  vel[2] = -velmag*( xrel[0] - xrel[1]);
#endif
  a_vel(0) = vel[a_dir];
}
PROTO_KERNEL_END(InitializeVFaceF, InitializeVFace)

#include "Chombo_NamespaceFooter.H"

#endif
