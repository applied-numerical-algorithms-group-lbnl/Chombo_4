cmake_minimum_required(VERSION 3.14.2)

if(DEFINED CHOMBO_LOADED)
    return() # Chombo has already been loaded by encompassing CMake project
endif()
# checks if Chombo is submodule of larger CMake project but hasn't been loaded yet
get_directory_property(hasParent PARENT_DIRECTORY)
if(hasParent)
    set(CHOMBO_LOADED "${Chombo_VERSION_MAJOR}.${Chombo_VERSION_MINOR}" PARENT_SCOPE)
endif()

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)  # if OFF, compilation will revert to older standard if C++11 unavailable
set(CMAKE_CXX_EXTENSIONS OFF)  # if OFF, compilation flag uses c++11, not gnu++11

# setup basics
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Debug" CACHE STRING "Choose build type, \
options are: Debug Release RelWithDebInfo" FORCE)
else()
    message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")
endif()

option(BUILD_SHARED_LIBS "Build using shared libraries" OFF)
set(CMAKE_POSITION_INDEPENEDENT_CODE TRUE)
option(ENABLE_MPI "Build using MPI" OFF)
option(ENABLE_OPENMP "Build using OpenMP" OFF)
option(ENABLE_CUDA "Build using CUDA" OFF)
option(ENABLE_HIP "Build using HIP" OFF)
option(ENABLE_HDF5 "Build using HDF5" OFF)
option(ENABLE_TESTS "Build tests" ON)
option(ENABLE_EXAMPLES "Build examples" ON)

set(DIM 2 CACHE STRING "Dimensionality of examples")
set(PREC DOUBLE CACHE STRING "Precision of examples")
if(PREC STREQUAL DOUBLE)
    add_compile_definitions(PR_USE_DOUBLE)
endif()
message(STATUS "Dimensionality: ${DIM}\n   Precision: ${PREC}")

add_library(chombo_compiler_flags INTERFACE)
add_compile_definitions(DIM=${DIM})
if(NOT DEFINED PROTO_HOME)
    set(PROTO_HOME $ENV{HOME}/proto CACHE STRING "Proto home directory")
endif()
set(CHOMBO_HOME ${CMAKE_CURRENT_SOURCE_DIRECTORY} CACHE PATH "Chombo home directory")

if(CMAKE_BUILD_TYPE STREQUAL "Debug")
    set(CMAKE_CXX_FLAGS "-w" CACHE STRING "compiler flags")
endif()

project(Chombo VERSION 4.0 DESCRIPTION "Chombo examples" LANGUAGES CXX C)

# setup packages
if(ENABLE_CUDA OR ENABLE_HIP)
    if(ENABLE_CUDA AND ENABLE_HIP)
        message(FATAL_ERROR "Must specify at most one of CUDA and HIP")
    endif()
    target_compile_definitions(chombo_compiler_flags INTERFACE CHOMBO_CUDA)
    if(CMAKE_BUILD_TYPE STREQUAL Debug)
        set(CMAKE_CUDA_FLAGS " -G")
        target_compile_definitions(chombo_compiler_flags INTERFACE THRUST_DEBUG)
    endif()
    if(ENABLE_HIP)
        if(NOT DEFINED ENV{HIP_PATH})
            set(HIP_PATH "/opt/rocm/hip" CACHE PATH "Path to find HIP")
        else()
            set(HIP_PATH $ENV{HIP_PATH} CACHE PATH "Path to find HIP")
        endif()
        if(NOT DEFINED ENV{CMAKE_PREFIX_PATH})
            set(CMAKE_PREFIX_PATH "${HIP_PATH} ${HIP_PATH}/.." CACHE PATH "Path to find Rocm")
        else()
            set(CMAKE_PREFIX_PATH $ENV{CMAKE_PREFIX_PATH} CACHE PATH "Path to find Rocm")
        endif()
        set(CMAKE_MODULE_PATH "${HIP_PATH}/cmake" ${CMAKE_MODULE_PATH})
        set(LIB_DEP hip)
        target_compile_definitions(chombo_compiler_flags INTERFACE CHOMBO_HIP)
    else()
        if(NOT DEFINED ENV{CUDA_DIR} AND NOT DEFINED ENV{CUDAPATH} AND NOT DEFINED ENV{CUDA_HOME})
            message(FATAL_ERROR "No CUDA compiler found")
        elseif(DEFINED ENV{CUDA_HOME})
            set(CUDA_TOOLKIT_ROOT_DIR $ENV{CUDA_HOME} CACHE PATH "Makes BLT look in right place")
            set(CMAKE_CUDA_COMPILER $ENV{CUDA_HOME}/bin/nvcc CACHE PATH "CUDA compiler")
        endif()
        set(CMAKE_CUDA_HOST_COMPILER ${CMAKE_CXX_COMPILER} CACHE PATH "Setting for BLT")
        set(CMAKE_CUDA_SEPARABLE_COMPILATION ON CACHE BOOL "")
        enable_language(CUDA)
        if(NOT DEFINED CMAKE_CUDA_ARCHITECTURES)
            set(CMAKE_CUDA_ARCHITECTURES 70 CACHE STRING "Compute architecture needed by NVCC")
        endif()
        if(CMAKE_CUDA_COMPILER_VERSION VERSION_LESS 11)
            message(FATAL_ERROR "CUDA version >=11 needed, found ${CMAKE_CUDA_COMPILER_VERSION}")
        elseif(CMAKE_CUDA_COMPILER_VERSION VERSION_EQUAL 11)
            if(CMAKE_COMPILER_ID STREQUAL "GNU" AND NOT CMAKE_CXX_COMPILER_VERSION VERSION_EQUAL 9)
                message(FATAL_ERROR "GCC version 9 needed, found ${CMAKE_CXX_COMPILER_VERSION}")
            endif()
        endif()
        set(CMAKE_CUDA_FLAGS "--expt-relaxed-constexpr --expt-extended-lambda" CACHE STRING "used in Reduction class" FORCE)
        set(LIB_DEP cuda)
    endif()
endif()

if(ENABLE_HDF5)
    add_compile_definitions(PR_HDF5)
    find_package(HDF5 REQUIRED COMPONENTS CXX HL)
    if(HDF5_IS_PARALLEL AND NOT ENABLE_MPI)
        message(WARNING "Parallel HDF5 module is loaded -- enabling MPI")
        set(ENABLE_MPI ON)
    endif()
endif()

if(ENABLE_MPI)
    target_compile_definitions(chombo_compiler_flags INTERFACE PR_MPI)
    set(LIB_DEP ${LIB_DEP} mpi)
endif()

if(ENABLE_OPENMP)
    set(LIB_DEP ${LIB_DEP} openmp)
endif()

if (NOT BLT_LOADED)
  if (DEFINED BLT_SOURCE_DIR)
      if (NOT EXISTS ${BLT_SOURCE_DIR}/SetupBLT.cmake)
          message(FATAL_ERROR "Given BLT_SOURCE_DIR does not contain SetupBLT.cmake")
      endif()
  else ()
      set (BLT_SOURCE_DIR ${PROJECT_SOURCE_DIR}/blt CACHE PATH "")

      if (NOT EXISTS ${BLT_SOURCE_DIR}/SetupBLT.cmake)
          message(FATAL_ERROR "\
          The BLT submodule is not present. \
          If in git repository run the following two commands:\n \
          git submodule init\n \
          git submodule update")
      endif ()
  endif ()

# Default to C++11 if not set so GTest/GMock can build
  if (NOT BLT_CXX_STD)
      set(BLT_CXX_STD "c++11" CACHE STRING "")
  endif()

  include(${BLT_SOURCE_DIR}/SetupBLT.cmake)
endif()

add_subdirectory(src)

if(ENABLE_EXAMPLES)
    add_subdirectory(examples)
endif()

if(ENABLE_TESTS)
    add_subdirectory(tests)
endif()
