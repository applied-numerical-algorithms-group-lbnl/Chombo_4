# Chombo_4
# A lightweight algorithm development framework for finite volume algorithms.
## Courtesy of the Applied Numerical Algorithms Group
## Lawrence Berkeley National Laboratory

## Introduction
# Chombo_4 is a lightweight algorithm development framework used for finite volume calculations.
This version of Chombo is fortran-free and depends on the Proto infrasturcture
for performance portability.   Eigen is used for linear algebra.   MPI is used for
off-node communication.  HDF5 is used for data output.

## Build instructions:
* Put machine and compiler specifics into Chombo4/mk/Make.defs.local
* Go to the example directory.
* Use the configure.example python script to create makefiles. It is simple and well-documented.
* Go to the example in which you are most interested.
* Type make.

##Directories
* example has many established algorithms and some in progress
* src is where the common code lives
* documents has a few documents created using this software
* mk  is where compiler and machine specific flags are set
* util has a few random tools that I kept around.

