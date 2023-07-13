# Chombo_4/example/DomainDecomposition/EBHelmholtz:
* This is where all the second-order, cut cell elliptic and parabolic stuff  lives.
* Anything with the hoeb_ prefix should be considered experimental.
* For more detailed documentation, all header files are designed to work with doxygen.

## Important Directories:
* exec solves the Helmholtz equation in a cut cell context by setting a charge distribution and solving for the field.
* ccProjExec sets a  cell-centered velocity and projects it onto its  divergence-free subspace.
* macProjExec sets a  face-centered velocity and projects it onto its  divergence-free subspace.



