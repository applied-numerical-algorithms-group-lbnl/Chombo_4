blt_add_executable(NAME MHD SOURCES MHD.cpp
    INCLUDES ${CMAKE_CURRENT_SOURCE_DIR}
    DEPENDS_ON Headers common MHDOpsteps)
