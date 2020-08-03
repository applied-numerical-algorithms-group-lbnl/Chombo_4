#include "ShapeArray.H"
#include "crunchflow.h"

using namespace shape;

double aqueoustobulkconvert(const int jx, 
			    const int jy, 
			    const int jz,
			    const int nx,
			    const int ny, 
			    const int nz, 
			    double *por_3d, 
			    double *satliq_3d, 
			    double *ro_3d)
{
  ShapeArray<double, 3> por(por_3d, nz, ny, nx);
  ShapeArray<double, 3> satliq(satliq_3d, nz, ny, nx);
  ShapeArray<double, 3> ro(ro_3d, nz, ny, nx);

  double massFraction = 1.0;
  return(massFraction * por[jz][jy][jx] * satliq[jz][jy][jx] * ro[jz][jy][jx]);
}
