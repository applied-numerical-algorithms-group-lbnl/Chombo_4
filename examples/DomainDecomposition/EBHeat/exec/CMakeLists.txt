blt_add_executable(NAME EBHeat
    HEADERS SetupFunctions.H SOURCES HeatEqn.cpp
    INCLUDES ${CMAKE_CURRENT_SOURCE_DIR}
    DEPENDS_ON Headers common)
