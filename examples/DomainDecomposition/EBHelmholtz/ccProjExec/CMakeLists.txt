message(STATUS "LIB_DEP = ${LIB_DEP}")
blt_add_executable(NAME ccProjExec
    HEADERS SetupFunctions.H
    SOURCES ccProject.cpp
    INCLUDES ${CMAKE_CURRENT_SOURCE_DIR}
    DEPENDS_ON Headers common ${LIB_DEP})
