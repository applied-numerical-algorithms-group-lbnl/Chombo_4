blt_add_executable(NAME EBSmallTest HEADERS SetupFunctions.H
    SOURCES smalltest.cpp INCLUDES ${CMAKE_CURRENT_SOURCE_DIR}
    DEPENDS_ON Headers common)
