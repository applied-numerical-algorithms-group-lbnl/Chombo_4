# Chombo
## Courtesy of the Applied Numerical Algorithms Group
## Lawrence Berkeley National Laboratory

## Introduction
This version of Chombo is fortran-free and depends on the Proto infrasturcture
for performance portability.   Eigen is used for linear algebra.   MPI is used for
off-node communication.  HDF5 is used for data output.

## Build instructions:
* Put machine and compiler specifics into Chombo4/mk/Make.defs.local
* Go to the example directory.
* Use the configure.example python script to create makefiles. It is simple and well-documented.
* Go to the example in which you are most interested.
* Type make.
