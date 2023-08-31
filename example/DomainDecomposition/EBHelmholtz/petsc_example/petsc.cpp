#include <cmath>
#include <cstdio>
#include <iostream>


#include "EBProto.H"
#include "Chombo_EBLevelBoxData.H"
#include "Chombo_LevelData.H"
#include "Chombo_BaseFab.H"

#include "Chombo_ParmParse.H"
#include "Chombo_LoadBalance.H"
#include "Chombo_ProtoInterface.H"
#include "Chombo_BRMeshRefine.H"
#include "Chombo_GeometryService.H"
#include "Chombo_EBDictionary.H"
#include "Chombo_EBChombo.H"
#include "EBMultigrid.H"

#include "DebugFunctions.H"
#include <iomanip>


/*
  Include "petscksp.h" so that we can use KSP solvers.
*/
#include <petscksp.h>
///
int
simplePetscTest(int argc, char **args)
{
  //from petsc's example 2:
  static char help[] = "Solves a linear system in parallel with KSP.\n\
Input parameters include:\n\
  -view_exact_sol   : write exact solution vector to stdout\n\
  -m <mesh_x>       : number of mesh points in x-direction\n\
  -n <mesh_y>       : number of mesh points in y-direction\n\n";



  Vec         x, b, u; /* approx solution, RHS, exact solution */
  Mat         A;       /* linear system matrix */
  KSP         ksp;     /* linear solver context */
  PetscReal   norm;    /* norm of solution error */
  PetscInt    i, j, Ii, J, Istart, Iend, m = 8, n = 7, its;
  PetscBool   flg;
  PetscScalar v;

  PetscFunctionBeginUser;
  PetscInitialize(&argc, &args, (char *)0, help);
  PetscOptionsGetInt(NULL, NULL, "-m", &m, NULL);
  PetscOptionsGetInt(NULL, NULL, "-n", &n, NULL);
  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Compute the matrix and right-hand-side vector that define
     the linear system, Ax = b.
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  /*
    Create parallel matrix, specifying only its global dimensions.
    When using MatCreate(), the matrix format can be specified at
    runtime. Also, the parallel partitioning of the matrix is
    determined by PETSc at runtime.

    Performance tuning note:  For problems of substantial size,
    preallocation of matrix memory is crucial for attaining good
    performance. See the matrix chapter of the users manual for details.
  */
  MatCreate(PETSC_COMM_WORLD, &A);
  MatSetSizes(A, PETSC_DECIDE, PETSC_DECIDE, m * n, m * n);
  MatSetFromOptions(A);
  MatMPIAIJSetPreallocation(A, 5, NULL, 5, NULL);
  MatSeqAIJSetPreallocation(A, 5, NULL);
  MatSeqSBAIJSetPreallocation(A, 1, 5, NULL);
  MatMPISBAIJSetPreallocation(A, 1, 5, NULL, 5, NULL);
  MatMPISELLSetPreallocation(A, 5, NULL, 5, NULL);
  MatSeqSELLSetPreallocation(A, 5, NULL);

  /*
    Currently, all PETSc parallel matrix formats are partitioned by
    contiguous chunks of rows across the processors.  Determine which
    rows of the matrix are locally owned.
  */
  MatGetOwnershipRange(A, &Istart, &Iend);

   /*
     Set matrix elements for the 2-D, five-point stencil in parallel.
     - Each processor needs to insert only elements that it owns
     locally --but any non-local elements will be sent to the
     appropriate processor during matrix assembly.
     - Always specify global rows and columns of matrix entries.

     Note: this uses the less common natural ordering that orders first
     all the unknowns for x = h then for x = 2h etc; Hence you see J = Ii +- n
     instead of J = I +- m as you might expect. The more standard ordering
     would first do all variables for y = h, then y = 2h etc.

   */
  for (Ii = Istart; Ii < Iend; Ii++)
  {
    v = -1.0;
    i = Ii / n;
    j = Ii - i * n;
    if (i > 0)
    {
      J = Ii - n;
      MatSetValues(A, 1, &Ii, 1, &J, &v, ADD_VALUES);
    }
    if (i < m - 1)
    {
      J = Ii + n;
      MatSetValues(A, 1, &Ii, 1, &J, &v, ADD_VALUES);
    }
    if (j > 0)
    {
      J = Ii - 1;
      MatSetValues(A, 1, &Ii, 1, &J, &v, ADD_VALUES);
    }
    if (j < n - 1)
    {
      J = Ii + 1;
      MatSetValues(A, 1, &Ii, 1, &J, &v, ADD_VALUES);
    }
    v = 4.0;
    MatSetValues(A, 1, &Ii, 1, &Ii, &v, ADD_VALUES);
  }

/*
  Assemble matrix, using the 2-step process:
  MatAssemblyBegin(), MatAssemblyEnd()
  Computations can be done while messages are in transition
  by placing code between these two statements.
*/
  MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);

/* A is symmetric. Set symmetric flag to enable ICC/Cholesky preconditioner */
  MatSetOption(A, MAT_SYMMETRIC, PETSC_TRUE);

/*
  Create parallel vectors.
  - We form 1 vector from scratch and then duplicate as needed.
  - When using VecCreate(), VecSetSizes and VecSetFromOptions()
  in this example, we specify only the
  vector's global dimension; the parallel partitioning is determined
  at runtime.
  - When solving a linear system, the vectors and matrices MUST
  be partitioned accordingly.  PETSc automatically generates
  appropriately partitioned matrices and vectors when MatCreate()
  and VecCreate() are used with the same communicator.
  - The user can alternatively specify the local vector and matrix
  dimensions when more sophisticated partitioning is needed
  (replacing the PETSC_DECIDE argument in the VecSetSizes() statement
  below).
*/
  VecCreate(PETSC_COMM_WORLD, &u);
  VecSetSizes(u, PETSC_DECIDE, m * n);
  VecSetFromOptions(u);
  VecDuplicate(u, &b);
  VecDuplicate(b, &x);

/*
  Set exact solution; then compute right-hand-side vector.
  By default we use an exact solution of a vector with all
  elements of 1.0;
*/
  VecSet(u, 1.0);
  MatMult(A, u, b);

/*
  View the exact solution vector if desired
*/
flg = PETSC_FALSE;
PetscOptionsGetBool(NULL, NULL, "-view_exact_sol", &flg, NULL);
if (flg) VecView(u, PETSC_VIEWER_STDOUT_WORLD);

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   Create the linear solver and set various options
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
KSPCreate(PETSC_COMM_WORLD, &ksp);

/*
  Set operators. Here the matrix that defines the linear system
  also serves as the preconditioning matrix.
*/
KSPSetOperators(ksp, A, A);

/*
  Set linear solver defaults for this problem (optional).
  - By extracting the KSP and PC contexts from the KSP context,
  we can then directly call any KSP and PC routines to set
  various options.
  - The following two statements are optional; all of these
  parameters could alternatively be specified at runtime via
  KSPSetFromOptions().  All of these defaults can be
  overridden at runtime, as indicated below.
*/
KSPSetTolerances(ksp, 1.e-2 / ((m + 1) * (n + 1)), 1.e-50, PETSC_DEFAULT, PETSC_DEFAULT);

/*
  Set runtime options, e.g.,
  -ksp_type <type> -pc_type <type> -ksp_monitor -ksp_rtol <rtol>
  These options will override those specified above as long as
  KSPSetFromOptions() is called _after_ any other customization
  routines.
*/
KSPSetFromOptions(ksp);

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   Solve the linear system
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

KSPSolve(ksp, b, x);

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   Check the solution and clean up
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
VecAXPY(x, -1.0, u);
VecNorm(x, NORM_2, &norm);
KSPGetIterationNumber(ksp, &its);

/*
  Print convergence information.  PetscPrintf() produces a single
  print statement from all processes that share a communicator.
  An alternative is PetscFPrintf(), which prints to a file.
*/
PetscPrintf(PETSC_COMM_WORLD, "Norm of error %g iterations %" PetscInt_FMT "\n", (double)norm, its);

/*
  Free work space.  All PETSc objects should be destroyed when they
  are no longer needed.
*/
KSPDestroy(&ksp);
VecDestroy(&u);
VecDestroy(&x);
VecDestroy(&b);
MatDestroy(&A);

/*
  Always call PetscFinalize() before exiting a program.  This routine
  - finalizes the PETSc libraries as well as MPI
  - provides summary and diagnostic information if certain runtime
  options are chosen (e.g., -log_view).
*/
PetscFinalize();
return 0;



}

///
class
minimalPetsc
{
public:
  minimalPetsc(int  a_numPts,
               bool a_printStuff)
  {
    createMatrixAndVectors(a_printStuff);
    setupSolver()
    formMatrix(a_printStuff)
  }
  
};



///
int
modifiedPetscTest(int argc, char **args)
{

  int numPts = 56;


  Vec         x, b, u; /* approx solution, RHS, exact solution */
  Mat         A;       /* linear system matrix */
  KSP         ksp;     /* linear solver context */
  PetscReal   norm;    /* norm of solution error */
  PetscInt    i, j, Ii, J, Istart, Iend, m = 8, n = 7, its;
  PetscBool   flg;
  PetscScalar v;

  PetscFunctionBeginUser;
  PetscInitialize(&argc, &args, (char *)0, help);
  PetscOptionsGetInt(NULL, NULL, "-m", &m, NULL);
  PetscOptionsGetInt(NULL, NULL, "-n", &n, NULL);
  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Compute the matrix and right-hand-side vector that define
     the linear system, Ax = b.
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  /*
    Create parallel matrix, specifying only its global dimensions.
    When using MatCreate(), the matrix format can be specified at
    runtime. Also, the parallel partitioning of the matrix is
    determined by PETSc at runtime.

    Performance tuning note:  For problems of substantial size,
    preallocation of matrix memory is crucial for attaining good
    performance. See the matrix chapter of the users manual for details.
  */
  MatCreate(PETSC_COMM_WORLD, &A);
  MatSetSizes(A, PETSC_DECIDE, PETSC_DECIDE, m * n, m * n);
  MatSetFromOptions(A);
  MatMPIAIJSetPreallocation(A, 5, NULL, 5, NULL);
  MatSeqAIJSetPreallocation(A, 5, NULL);
  MatSeqSBAIJSetPreallocation(A, 1, 5, NULL);
  MatMPISBAIJSetPreallocation(A, 1, 5, NULL, 5, NULL);
  MatMPISELLSetPreallocation(A, 5, NULL, 5, NULL);
  MatSeqSELLSetPreallocation(A, 5, NULL);

  /*
    Currently, all PETSc parallel matrix formats are partitioned by
    contiguous chunks of rows across the processors.  Determine which
    rows of the matrix are locally owned.
  */
  MatGetOwnershipRange(A, &Istart, &Iend);

   /*
     Set matrix elements for the 2-D, five-point stencil in parallel.
     - Each processor needs to insert only elements that it owns
     locally --but any non-local elements will be sent to the
     appropriate processor during matrix assembly.
     - Always specify global rows and columns of matrix entries.

     Note: this uses the less common natural ordering that orders first
     all the unknowns for x = h then for x = 2h etc; Hence you see J = Ii +- n
     instead of J = I +- m as you might expect. The more standard ordering
     would first do all variables for y = h, then y = 2h etc.

   */
  for (Ii = Istart; Ii < Iend; Ii++)
  {
    v = -1.0;
    i = Ii / n;
    j = Ii - i * n;
    if (i > 0)
    {
      J = Ii - n;
      MatSetValues(A, 1, &Ii, 1, &J, &v, ADD_VALUES);
    }
    if (i < m - 1)
    {
      J = Ii + n;
      MatSetValues(A, 1, &Ii, 1, &J, &v, ADD_VALUES);
    }
    if (j > 0)
    {
      J = Ii - 1;
      MatSetValues(A, 1, &Ii, 1, &J, &v, ADD_VALUES);
    }
    if (j < n - 1)
    {
      J = Ii + 1;
      MatSetValues(A, 1, &Ii, 1, &J, &v, ADD_VALUES);
    }
    v = 4.0;
    MatSetValues(A, 1, &Ii, 1, &Ii, &v, ADD_VALUES);
  }

/*
  Assemble matrix, using the 2-step process:
  MatAssemblyBegin(), MatAssemblyEnd()
  Computations can be done while messages are in transition
  by placing code between these two statements.
*/
  MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);

/* A is symmetric. Set symmetric flag to enable ICC/Cholesky preconditioner */
  MatSetOption(A, MAT_SYMMETRIC, PETSC_TRUE);

/*
  Create parallel vectors.
  - We form 1 vector from scratch and then duplicate as needed.
  - When using VecCreate(), VecSetSizes and VecSetFromOptions()
  in this example, we specify only the
  vector's global dimension; the parallel partitioning is determined
  at runtime.
  - When solving a linear system, the vectors and matrices MUST
  be partitioned accordingly.  PETSc automatically generates
  appropriately partitioned matrices and vectors when MatCreate()
  and VecCreate() are used with the same communicator.
  - The user can alternatively specify the local vector and matrix
  dimensions when more sophisticated partitioning is needed
  (replacing the PETSC_DECIDE argument in the VecSetSizes() statement
  below).
*/
  VecCreate(PETSC_COMM_WORLD, &u);
  VecSetSizes(u, PETSC_DECIDE, m * n);
  VecSetFromOptions(u);
  VecDuplicate(u, &b);
  VecDuplicate(b, &x);

/*
  Set exact solution; then compute right-hand-side vector.
  By default we use an exact solution of a vector with all
  elements of 1.0;
*/
  VecSet(u, 1.0);
  MatMult(A, u, b);

/*
  View the exact solution vector if desired
*/
flg = PETSC_FALSE;
PetscOptionsGetBool(NULL, NULL, "-view_exact_sol", &flg, NULL);
if (flg) VecView(u, PETSC_VIEWER_STDOUT_WORLD);

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   Create the linear solver and set various options
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
KSPCreate(PETSC_COMM_WORLD, &ksp);

/*
  Set operators. Here the matrix that defines the linear system
  also serves as the preconditioning matrix.
*/
KSPSetOperators(ksp, A, A);

/*
  Set linear solver defaults for this problem (optional).
  - By extracting the KSP and PC contexts from the KSP context,
  we can then directly call any KSP and PC routines to set
  various options.
  - The following two statements are optional; all of these
  parameters could alternatively be specified at runtime via
  KSPSetFromOptions().  All of these defaults can be
  overridden at runtime, as indicated below.
*/
KSPSetTolerances(ksp, 1.e-2 / ((m + 1) * (n + 1)), 1.e-50, PETSC_DEFAULT, PETSC_DEFAULT);

/*
  Set runtime options, e.g.,
  -ksp_type <type> -pc_type <type> -ksp_monitor -ksp_rtol <rtol>
  These options will override those specified above as long as
  KSPSetFromOptions() is called _after_ any other customization
  routines.
*/
KSPSetFromOptions(ksp);

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   Solve the linear system
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

KSPSolve(ksp, b, x);

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   Check the solution and clean up
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
VecAXPY(x, -1.0, u);
VecNorm(x, NORM_2, &norm);
KSPGetIterationNumber(ksp, &its);

/*
  Print convergence information.  PetscPrintf() produces a single
  print statement from all processes that share a communicator.
  An alternative is PetscFPrintf(), which prints to a file.
*/
PetscPrintf(PETSC_COMM_WORLD, "Norm of error %g iterations %" PetscInt_FMT "\n", (double)norm, its);

/*
  Free work space.  All PETSc objects should be destroyed when they
  are no longer needed.
*/
KSPDestroy(&ksp);
VecDestroy(&u);
VecDestroy(&x);
VecDestroy(&b);
MatDestroy(&A);

/*
  Always call PetscFinalize() before exiting a program.  This routine
  - finalizes the PETSc libraries as well as MPI
  - provides summary and diagnostic information if certain runtime
  options are chosen (e.g., -log_view).
*/
PetscFinalize();
return 0;

}

int main(int a_argc, char* a_argv[])
{
//  int retval = simplePetscTest(a_argc, a_argv);

  int retval = modifiedPetscTest(a_argc, a_argv);
  return retval;
}
