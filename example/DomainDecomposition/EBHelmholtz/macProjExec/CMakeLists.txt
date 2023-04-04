blt_add_executable(NAME macProjExec HEADERS SetupFunctions.H
    SOURCES macProject.cpp INCLUDES ${CMAKE_CURRENT_SOURCE_DIR}
    DEPENDS_ON Headers common ${LIB_DEP})
