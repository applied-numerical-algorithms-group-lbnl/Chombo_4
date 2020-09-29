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
#include "Chombo_LAPACKMatrix.H"



/***/
int main(int argc, char* argv[])
{
  int nrow = 2;
  int ncol = 2;


  typedef Chombo4::LAPACKMatrix Matrix;
  Matrix mat(nrow, ncol);
  Real a = 1;  Real b = 2; Real c = -1; Real d = 3;
  
  mat(0, 0) = a;
  mat(0, 1) = b;
  mat(1, 0) = c;
  mat(1, 1) = d;
   
  Matrix inv = mat;
  int retval = inv.invertUsingLeastSquares();
  if(retval != 0)
  {
    printf("matrix inversion returned error code so it FAILED\n");
    return retval;
  }

  Matrix product;
  multiply(product, mat, inv);
  
  Real tol = 1.0e-6;
  for(int irow = 0; irow < nrow; irow++)
  {
    for(int icol = 0; icol < ncol; icol++)
    {
      Real rightans = 0;
      if(irow == icol) rightans = 1;
      Real diff = std::abs(product(irow, icol) - rightans);
      if(diff > tol)
      {
        printf("inversion produced difference too big = %f so it FAILED\n", diff);
        return -1;
      }
    }
  }
  printf("matrix test PASSED \n");
  return 0;

}
