# Chombo
## Courtesy of the Applied Numerical Algorithms Group
## Lawrence Berkeley National Laboratory

## Introduction
This version of Chombo is fortran-free and depends on the Proto infrasturcture  for performance portability.
On GPU devices, Proto's data holders are used.  Chombo maintains host-based data holders.   Chombo 
handles all MPI and HDF5 interactions.

## Build instructions on Cori:
### Modules
* The default CMake version is insufficient. The minumum required version is 3.17, and the minimum acceptable version available as a module on Cori is 3.18. Configuration has been tested through version 3.21.
* If doing a CUDA build, load the `cudatoolkit` module.
* If using HDF5, load the `cray-hdf5` module

### Submodules
* There are two submodules within this repository. One is the (Proto)[https://github.com/applied-numerical-algorithms-group-lbnl/proto_import] library that has headers and data types used by Chombo. The other is the (BLT)[https://github.com/LLNL/blt] library of CMake macros used during configuration.
* If you *haven't* yet cloned this repository, include the `--recurse-submodules` flag in your `git clone` command so that the two submodules will be initialized and updated.
* If you *have* already cloned this repository, run these two commands to checkout the commits in each submodule that this project links to:
   - `git submodule init`
   - `git submodule update`
   
   After doing this, the `.gitmodules` file will show the path and url of each submodule.

### Configuring
* The simplest command assumes you are at the top directory of the repo and is `cmake -S . -B build`. The argument to `-S` is the source directory containing the top-level `CMakeLists.txt` file. The argument to `-B` is where the binaries should be built. The specified directory does not need to exist when the `cmake` command is invoked. Additionally, there are various options which can be set during this step, each preceded by the `-D` flag. Valid choices are listed in brackets, with defaults in bold. They are:
   - ENABLE_CUDA=\[ON, *OFF*\]
   - ENABLE_MPI=\[ON, *OFF*\]
   - ENABLE_HDF5=\[ON, *OFF*\]
   - Build executables from the `examples` subdirectory: ENABLE_EXAMPLES=\[*ON*, OFF\]
   - Build executables from the `tests` subdirectory: ENABLE_TESTS=\[*ON*, OFF\]
   - Floating point precision: PREC=\[SINGLE, *DOUBLE*\]
   - Dimensionality of examples: DIM=\[*2*, 3\]
   - Optimization level: CMAKE_BUILD_TYPE=\[*Debug*, Release, MinSizeRel, RelWithDebInfo\]
   - Size of allocations from the stack: STACK=<int> (default is 4294967296 bytes)
   - Root of Proto repository: PROTO_HOME=<path> (default is submodule within this Chombo repository)
   
### Building
* To build all of the CMake targets in this project in parallel, run `cmake --build <build-dir> --parallel`. An integer can be provided at the end of this command to set the level of parallelization.
   
### Testing
* After moving to the build directory, run all compiled tests by giving the command `make test` to a SLURM script or `salloc` line. This command must be invoked on the proper node (eg. one with GPUs if doing a CUDA run) for the tests to run successfully. The BLT submodule used by this project has various "smoke" tests that can be used to determine if the tests are being run on the proper architecture for the desired configuration. They are automatically included in the `test` target.
