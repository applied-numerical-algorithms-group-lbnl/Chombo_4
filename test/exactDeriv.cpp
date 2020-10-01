#include <cstdio>
#include <cstring>
#include <cassert>
#include <cmath>

#include <vector>
#include <memory>

#include <iostream>
#include <fstream>
#include <sstream>

#include "Proto.H"
#include "EBProto.H"
#include "Chombo_EBChombo.H"
#include "Chombo_GeometryService.H"
#include "Chombo_EBLevelBoxData.H"
#include "Chombo_LAPACKMatrix.H"
#include "Proto_EBExactSolutions.H"
#include "Chombo_EBDictionary.H"

/***/
template< unsigned int order>
void
runderivs(int a_nx)
{

  std::pair<Point, Real> entry;
  entry.first  = Point::Ones(order);
  entry.second = order;
  RealVect loca = RealVect::Unit();
  loca *= 7.;
  MonomialEF<order> funk(entry);
  IndexTM<int, DIM> curDeriv= IndexTM<int, DIM>::Zero;
  pout() << "monomial derivs" << endl;
  for(int xderiv = 0; xderiv <= order; xderiv++)
  {
    curDeriv[0] = xderiv;
    for(int yderiv = 0; yderiv <= order; yderiv++)
    {
      curDeriv[1] = yderiv;
      int zderiv = 0;
#if DIM==3
      for(zderiv = 0; zderiv <= order; zderiv++)
      {
        curDeriv[2] = zderiv;
#endif
        int sum = curDeriv.sum();
        if(sum <= order)
        {
          Real derivval = funk.getDerivative(curDeriv, loca);
          pout() << "(xd yd zd) derivval: ("
                 << xderiv << ", "<< yderiv << ", " << zderiv << ")= "
                 << derivval << std::endl;
        }
#if DIM==3
      }
#endif    
    }
  }

  pout() << "sinesphere derivs" << endl;
  Real radius = 0.45;
  Real center = 0.5;
  SineSphereEF<order> soulChicken(radius, center);
  for(int xderiv = 0; xderiv < 5; xderiv++)
  {
    curDeriv[0] = xderiv;
    for(int yderiv = 0; yderiv < 5; yderiv++)
    {
      curDeriv[1] = yderiv;
      int zderiv = 0;
#if DIM==3
      for(zderiv = 0; zderiv < 5; zderiv++)
      {
        curDeriv[2] = zderiv;
#endif
        int sum = curDeriv.sum();
        if(sum < 5)
        {
          Real derivval = soulChicken.getDerivative(curDeriv, loca);
          pout() << "(xd yd zd) derivval: ("
                 << xderiv << ", "<< yderiv << ", " << zderiv << ")= "
                 << derivval << std::endl;
        }
#if DIM==3
      }
#endif    
    }
  }
    
}
/***/
int main(int argc, char* argv[])
{
  int nx = 64;
  
  runderivs<4>(nx);
  printf("runDerivs test PASSED with  \n");

  return 0;

}
