blt_add_executable(NAME ebforall HEADERS HostDebugFunctions.H
    SOURCES ebforallTest.cpp INCLUDES ${CMAKE_CURRENT_SOURCE_DIR}
    DEPENDS_ON Headers common ${LIB_DEP})
