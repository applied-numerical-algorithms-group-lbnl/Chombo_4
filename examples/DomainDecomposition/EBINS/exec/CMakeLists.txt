blt_add_executable(NAME navier
    HEADERS SetupFunctions.H
    SOURCES navier.cpp
    INCLUDES ${CMAKE_CURRENT_SOURCE_DIR}
    DEPENDS_ON Headers common)
