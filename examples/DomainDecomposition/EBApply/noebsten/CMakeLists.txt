blt_add_executable(NAME noebsten
    HEADERS HostDebugFunctions.H
    SOURCES noebstentest.cpp
    INCLUDES ${CMAKE_CURRENT_SOURCE_DIR}
    DEPENDS_ON Headers common ${LIB_DEP})
