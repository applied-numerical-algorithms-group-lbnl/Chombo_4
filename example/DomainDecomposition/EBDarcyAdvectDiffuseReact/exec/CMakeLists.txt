blt_add_executable(NAME darcy
    HEADERS SetupFunctions.H SOURCES darcy.cpp
    INCLUDES ${CMAKE_CURRENT_SOURCE_DIR}
    DEPENDS_ON Headers common EBDarcy)
