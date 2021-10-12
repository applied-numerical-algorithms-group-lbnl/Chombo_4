blt_add_executable(NAME exchangeTest
    HEADERS HostDebugFunctions.H SetupFunctions.H
    SOURCES exchangeTest.cpp
    INCLUDES ${CMAKE_CURRENT_SOURCE_DIR}
    DEPENDS_ON Headers common)
