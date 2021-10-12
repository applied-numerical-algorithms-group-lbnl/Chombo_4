blt_add_executable(NAME euler SOURCES euler.cpp
    INCLUDES ${CMAKE_CURRENT_SOURCE_DIR}
    DEPENDS_ON Headers common EulerOpRK4)
