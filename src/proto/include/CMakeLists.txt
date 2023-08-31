blt_add_library(NAME Headers_Base HEADERS Proto.H
                INCLUDES ${CMAKE_CURRENT_SOURCE_DIR}
                ${CMAKE_CURRENT_SOURCE_DIR}/base
                ${CMAKE_CURRENT_SOURCE_DIR}/base/implem)
blt_add_library(NAME Headers_AMR HEADERS ProtoAMR.H
                INCLUDES ${CMAKE_CURRENT_SOURCE_DIR}
                ${CMAKE_CURRENT_SOURCE_DIR}/amr
                ${CMAKE_CURRENT_SOURCE_DIR}/amr/implem
                DEPENDS_ON Headers_Base)
blt_add_library(NAME Headers_MMB HEADERS ProtoMMB.H
                INCLUDES ${CMAKE_CURRENT_SOURCE_DIR}
                ${CMAKE_CURRENT_SOURCE_DIR}/mmb
                ${CMAKE_CURRENT_SOURCE_DIR}/mmb/implem
                DEPENDS_ON Headers_Base)
if(ENABLE_CUDA)
    target_link_libraries(Headers_Base INTERFACE nvToolsExt)
    target_link_libraries(Headers_AMR INTERFACE nvToolsExt)
    target_link_libraries(Headers_MMB INTERFACE nvToolsExt)
endif()
if(ENABLE_HIP)
    target_include_directories(Headers_Base INTERFACE $ENV{ROCM_PATH}/roctracer/include/)
    target_include_directories(Headers_AMR  INTERFACE $ENV{ROCM_PATH}/roctracer/include/)
    target_include_directories(Headers_MMB  INTERFACE $ENV{ROCM_PATH}/roctracer/include/)
    target_link_directories(Headers_Base INTERFACE $ENV{ROCM_PATH}/lib)
    target_link_directories(Headers_AMR  INTERFACE $ENV{ROCM_PATH}/lib)
    target_link_directories(Headers_MMB  INTERFACE $ENV{ROCM_PATH}/lib)
    target_link_libraries(Headers_Base INTERFACE roctx64)
    target_link_libraries(Headers_AMR INTERFACE roctx64)
    target_link_libraries(Headers_MMB INTERFACE roctx64)
endif()
if(ENABLE_MPI)
    target_link_libraries(Headers_Base INTERFACE MPI::MPI_CXX)
    target_link_libraries(Headers_AMR INTERFACE MPI::MPI_CXX)
    target_link_libraries(Headers_MMB INTERFACE MPI::MPI_CXX)
    target_include_directories(Headers_Base INTERFACE $ENV{MPICH_DIR}/include/)
endif()
if(ENABLE_OPENMP)
    target_link_libraries(Headers_Base INTERFACE OpenMP::OpenMP_CXX)
    target_link_libraries(Headers_AMR  INTERFACE OpenMP::OpenMP_CXX)
    target_link_libraries(Headers_MMB  INTERFACE OpenMP::OpenMP_CXX)
endif()
if(ENABLE_HDF5)
    target_link_libraries(Headers_Base INTERFACE HDF5::HDF5)
    target_link_libraries(Headers_AMR  INTERFACE HDF5::HDF5)
    target_link_libraries(Headers_MMB  INTERFACE HDF5::HDF5)
endif()
