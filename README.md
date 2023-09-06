# Chombo_4
* Free software courtesy of the Applied Numerical Algorithms Group
* Lawrence Berkeley National Laboratory
* Berkeley, California, USA.

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

