# Chombo_4  EBCM Current Release Branch
* Free software courtesy of the Applied Numerical Algorithms Group
* Lawrence Berkeley National Laboratory
* Berkeley, California, USA.

## Embedded boundary with cell merging (EBCM) is under current development.
* This branch will be where the most stable version of this technology resides.
* dtg_dev will be the relevant development branch.
* Currently working fine as of 6-27-2025:
* Meta-data for cell merging (getting all the moments worked out and so on).
* Data holders for cell merging.
* PETSc interface.
* Operator infrastructure works okay. Currently I only have implemented Helmholtz.
* Several technical documents (included) have been written using this stuf.
* What I will work on next (in order of appearance):
* Finish the writeup for the truncation error test.
* Projection operator so I can have divergence-free fields
* Advection and incompressible Navier Stokes to follow.
* Note 11: Hyperbolics is the point of this stuff, as cell merging makes hyperbolic stability much easier.
* Note 2: Gas dynamics is on hold for now as I have no idea how to do limiting with this stuff.


## Introduction
* Chombo_4 is a lightweight algorithm development framework used for finite volume calculations.
* Chombo_4 is free software (see the BSD-style license in Chombo_4/Copyright.txt).
* Chombo_4 is also blessedly  fortran-free.    
* All calculations are in C++.
* File management and build configuration tools are in Python.

## Dependencies
* proto is used for performance portability.
* Eigen is used for linear algebra.
* MPI is used for off-node communication.
* HDF5 is used for data.

## Build instructions:
* Put machine and compiler specifics into the file Chombo_4/mk/Make.defs.local.  The process is well-documented.
* Go to the Chombo_4/example directory.
* Use the configure.example python script to create makefiles. It is simple and well-documented.
* Go to the example directory in which you are most interested.
* Type make.

## Directories
* Chombo_4/example has many established algorithms and some in progress.
* Chombo_4/src is where the common code lives.
* Chombo_4/documents has a few documents created using this software.  It is also where doxygen output lives.
* Chombo_4/mk  is where compiler and machine specific flags are set.
* Chombo_4/util has a few random tools that I keep around.

