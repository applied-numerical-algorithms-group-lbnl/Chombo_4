blt_add_executable(
    NAME EBHelmholtz
    SOURCES HelmholtzSolve.cpp
    INCLUDES ${CMAKE_CURRENT_SOURCE_DIR}
    DEPENDS_ON Headers common)
