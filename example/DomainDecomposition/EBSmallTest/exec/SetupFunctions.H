#ifndef __SetupFunctionss__
#define __SetupFunctionss__
#define PI 3.141592653589793
#include "DebugFunctions.H"
#include "Chombo_NamespaceHeader.H"
using Proto::Var;


//=================================================
PROTO_KERNEL_START 
void InitCellF(int                 a_p[DIM],
               Var<Real, 1>       a_val)
{
  a_val(0) =  (Real(a_p[0]) + .5)*PI;
  a_val(0) += (Real(a_p[1]) - .5)*PI;
  if(a_p[1] %2 == 1) a_val(0) *= -1;
}
PROTO_KERNEL_END(InitCellF, InitCell)


//=================================================
PROTO_KERNEL_START 
void InitFaceF(int                 a_p[DIM],
               Var<Real, 1>       a_val,
               int idir)
{

  a_val(0) = Real(a_p[idir]);
  if(a_p[1] %2 == 1) a_val(0) *= -1;
}
PROTO_KERNEL_END(InitFaceF, InitFace)
#include "Chombo_NamespaceFooter.H"

#endif
