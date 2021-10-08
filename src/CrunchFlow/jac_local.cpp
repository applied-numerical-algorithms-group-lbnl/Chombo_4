/******************GIMRT98     ************************
      
!     Code converted using TO_F90 by Alan Miller
!     Date: 2000-07-27  Time: 09:59:43
      
!**************(C) COPYRIGHT 1995,1998,1999 ******************
!*******************C.I. Steefel      *******************
!     All Rights Reserved

!     GIMRT98 IS PROVIDED "AS IS" AND WITHOUT ANY WARRANTY EXPRESS OR IMPLIED.
!     THE USER ASSUMES ALL RISKS OF USING GIMRT98. THERE IS NO CLAIM OF THE
!     MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE.

!     YOU MAY MODIFY THE SOURCE CODE FOR YOUR OWN USE, BUT YOU MAY NOT
!     DISTRIBUTE EITHER THE ORIGINAL OR THE MODIFIED CODE TO ANY OTHER
!     WORKSTATIONS
!**********************************************************************/

#include <string.h>
#include "ShapeArray.H"
#include "crunchflow.h"

using namespace shape;

void jac_local(const int ncomp, 
	       const int nspec, 
	       const int neqn,
	       const int jx, 
	       const int jy, 
	       const int jz,
	       double *fjac_loc_2d, // [neqn][neqn],
	       const int nx, 
	       const int ny, 
	       const int nz,
	       double *muaq_2d, // [ncomp][nspec], 
	       /* const */ double *sp10_4d) // [nz][ny][nx][ncomp + nspec], 
{
  ShapeArray<double, 2> fjac_loc(fjac_loc_2d, neqn, neqn);
  ShapeArray<double, 2> muaq(muaq_2d, ncomp, nspec);
  ShapeArray<double, 4> sp10(sp10_4d, nz, ny, nx, ncomp + nspec);

  double spec_conc,mutemp;

  for(size_t i = 0; i < neqn; ++i)
    for(size_t j = 0; j < neqn; ++j)
      fjac_loc[i][j] = 0.0;

  for (int ksp = 1; ksp <= nspec; ++ksp) {
    spec_conc = sp10[jz-1][jy-1][jx-1][ksp + ncomp - 1]; // TODO: verify spec_conc
    for (int i = 1; i <= ncomp; ++i) {
      if (muaq[i-1][ksp-1] != 0.0) 
	{
	  mutemp = muaq[i-1][ksp-1];
	  for (int i2 = 1; i2 <= i-1; ++i2) 
	    {
	      fjac_loc[i-1][i2-1] = fjac_loc[i-1][i2-1] + mutemp * muaq[i2-1][ksp-1] * spec_conc;
	      fjac_loc[i2-1][i-1] = fjac_loc[i-1][i2-1];
	    }
	  fjac_loc[i-1][i-1] = fjac_loc[i-1][i-1] + mutemp * mutemp * spec_conc;
	}
    }
  }
  for (int i = 1; i <= ncomp; ++i) {
    fjac_loc[i-1][i-1] = fjac_loc[i-1][i-1] +	sp10[jz-1][jy-1][jx-1][i-1];
  }
  return;
}
