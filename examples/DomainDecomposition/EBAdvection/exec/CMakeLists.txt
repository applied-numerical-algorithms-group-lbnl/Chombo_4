blt_add_executable(NAME EBAdvection HEADERS SetupFunctions.H
    SOURCES ebadvect.cpp INCLUDES ${CMAKE_CURRENT_SOURCE_DIR}
    DEPENDS_ON Headers common)
