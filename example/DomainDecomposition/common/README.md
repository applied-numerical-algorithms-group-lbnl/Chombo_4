# Chombo_4/example/DomainDecomposition/common:
* Herein lie the common code for Chombo_4 examples.
* The cut cell tools  here do not include cell merger.
* For more detailed documentation, all header files are designed to work with doxygen.

## Notable tools herein:
* BiCGStabSolver is used a bottom solver in multigrid.
* EBParabolicIntegrators implements several parabolic algoritms, including backward Euler, Crank Nicolson, and TGA,
* DebugFunctions are some useful things you can call from gdb.  Mostly this prints out Chombo data in useful and concise ways.
* DumpArea is another debugging tool that allows you  print small areas of the domain in gdb.
* EBAdvection is the class that does second order advection in a cut cell context.
* EBCCProjector implements an cell-centered projection operator in a cut cell context.
* EBMACProjector implements a MAC projection operator in a cut cell context.
* EBMultigrid implements geometric multigrid in cut cell context.
* EBPetscSolver is an interface to the PETSc infrastructure for cut cell calculations.


