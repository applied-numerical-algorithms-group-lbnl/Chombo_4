/**
* Offsets depending on dimensions of arrays.
**/

/* 1D */
#define array1D(name,index1) ( *((name) + ((index1) - 1)) )
/* Matrix:
 * ivoume(nkin)
 * sppTMP(ncomp + nspec)
 * nreactmin(nkin)           INT
 * surf(nkin)
 * jac_sat(ncomp)
 * eqhom(nspec)              REAL
 * fxmax(neqn)               REAL
 */

/* 2D */
#define array2D(name,index1,index2,dim1) ( *((name) + ((index1) - 1) + ((index2) - 1)*(dim1)) )
/* Matrices:
 * AffinityDepend1(np,nkin)
 * si(np,nkin)
 * snorm(np,nkin)
 * ndepend(np,nkin)
 * jac_pre(ncomp,np)
 * pre_rmin(np,nkin)
 * muaq(nspec,ncomp)          REAL
 * rate0(np,nkin)
 * silog(np,nkin)             REAL
 * as1(nspec+nptot,5)         REAl
 */

/* 3D */
#define array3D(name,index1,index2,index3,dim1,dim2) ( *((name) + ((index1) - 1) + ((index2) - 1)*(dim1) + ((index3) - 1)*(dim1)*(dim2)) )
/* Matrices:
 * idepend(nd,np,nkin)
 * itot_min(nd,np,nkin)
 * depend(nd,np,nkin)
 * s(neqn,nx,ny,nz)
 * jac_rmin(ncomp,np,nkin)
 * mumin(np,nkin,ncomp)      REAL
 * por(nx,ny,nz)
 * satliq(nx,ny,nz)
 * ro(nx,ny,nz)
 * t(nx,ny,nz)               REAL
 * satliq(nx,ny,nz)          REAL
 */

/* 4D */
#define array4D(name,index1,index2,index3,index4,dim1,dim2,dim3) ( *((name) + ((index1) - 1) + ((index2) - 1)*(dim1) + ((index3) - 1)*(dim1)*(dim2) + ((index4) - 1)*(dim1)*(dim2)*(dim3)) )
/* Matrices:
 * gam(ncomp+nspec,nx,ny,nz)
 * sp(ncomp+nspec,nx,ny,nz)    REAL
 * s(neqn,nx,ny,nz)
 * sp10(ncomp+nspec,nx,ny,nz)
 * keqaq(nspec,nx,ny,nz)       REAL
 * dppt(nkin,nx,ny,nz)         REAL
 * u_rate(nkin,nx,ny,nz)
 */

/* 5D */
#define array5D(name,index1,index2,index3,index4,index5,dim1,dim2,dim3,dim4) ( *((name) + ((index1) - 1) + ((index2) - 1)*(dim1) + ((index3) - 1)*(dim1)*(dim2) + ((index4) - 1)*(dim1)*(dim2)*(dim3) + ((index5) - 1)*(dim1)*(dim2)*(dim3)*(dim4)) )
/* Matrices:
 * fjac(neqn,neqn,nx,ny,nz)
 * keqmin(np,nkin,nx,ny,nz)   REAL
 *
 *
 *
 *
 */
//( *((name) + ((index1) - 1) + ((index2) - 1)*(dim1) + ((index3) - 1)*(dim1)*(dim2) + ((index4) - 1)*(dim1)*(dim2)*(dim3) + ((index5) - 1)*(dim1)*(dim2)*(dim3)*(dim4)) )
