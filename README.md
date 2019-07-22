# Chombo

## Introduction
This version of Chombo is fortran-free and depends on the Proto infrasturcture  for performance portability.
On GPU devices, Proto's data holders are used.  Chombo maintains host-based data holders.   Chombo 
handles all MPI and HDF5 interactions.

## Build instructions:
1. Create the configuration file (Chombo_4/mk/Make.defs.local)
   - cp Make.defs.local.template Make.defs.local
     - There are other examples of configuration files in Chombo_4/mk/local
   - Edit the configuration file 
     - Specify your HDF5 location.
     - Specify your c++ compiler.
     - Specify your compiler flags.
     - Specify your PROTO location.
3. Configure the examples using Chombo_4/examples/configure.example
   - For example type 

```
<unix prompt>   ./configure --opt TRUE --mpi  FALSE --dim 3 --cuda TRUE 
```

   - The configuration parameters are given by:
```
<unix prompt>  ./configure.example --help
usage: configure.example [-h] [--dim DIM] [--opt {DEBUG,TRUE}] [--mpi {TRUE,FALSE}] [--hdf5 {TRUE,FALSE}] [--prec {SINGLE,DOUBLE}] [--cuda {TRUE,FALSE}]

optional arguments:
  -h, --help              show this help message and exit
  --dim DIM               dimensionality to build executables [1,2,3]
  --opt {DEBUG,TRUE}      compiler optimization [DEBUG]
  --mpi {TRUE,FALSE}      MPI on or off [FALSE]
  --hdf5 {TRUE,FALSE}     HDF5 on or off [FALSE]
  --prec {SINGLE,DOUBLE}  precision [DOUBLE]
  --cuda {TRUE,FALSE}     CUDA on or off [FALSE]
```
4.  This will create a GNUmakefile in each example directory.
5.  Go into the example dirctory and type 'make'
