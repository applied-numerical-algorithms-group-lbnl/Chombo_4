# Chombo_4/example/DomainDecomposition:
* Herein lie the Chombo_4 examples that are ports of Chombo3 examples.
* For more detailed documentation, all header files are designed to work with doxygen.

## Basic organization:
* As these are older examples developed in the service of larger things, each example has its own flavor.
* Where necessary, the important directories have  README files  which provide some guidance.
* Anything with the hoeb_ prefix should be considered experimental.

## The important directories:
* Chombo_4/example/DomainDecomposition/common has the shared example code.
* Chombo_4/example/DomainDecomposition/EBAdvection advects scalar given a fixed velocity in a cut cell context.
* Chombo_4/example/DomainDecomposition/EBHelmholz solves Poisson/Helmholtz equations in a cut cell context.  This includes projections.
* Chombo_4/example/DomainDecomposition/EBMenagerie is a pedagogical example which shows how to make cut cell geometries.
* Chombo_4/example/DomainDecomposition/EBINS solves the incompressible Navier Stokes equations in a cut cell context.


