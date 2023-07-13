# Chombo_4/example/DomainDecomposition/common:
* Herein lie the common code for Chombo_4 examples.
* For more detailed documentation, all header files are designed to work with doxygen.

## Basic organization:
* As these are older examples developed in the service of larger things, each example has its own flavor.
* The important directories have  README files  which provide some guidance.
* Anything with the hoeb_ prefix should be considered experimental.
* The cut cell examples here do not include cell merger.

## The important directories:
* Chombo_4/example/DomainDecomposition/common has the shared example code.
* Chombo_4/example/DomainDecomposition/EBAdvection advects scalar given a fixed velocity in a cut cell context.
* Chombo_4/example/DomainDecomposition/EBHelmholz solves Poisson/Helmholtz equations in a cut cell context.
* Chombo_4/example/DomainDecomposition/EBMenagerie is a pedagogical example which shows how to make cut cell grids.
* Chombo_4/example/DomainDecomposition/EBINS solves the incompressible Navier Stokes equations in a cut cell context.


