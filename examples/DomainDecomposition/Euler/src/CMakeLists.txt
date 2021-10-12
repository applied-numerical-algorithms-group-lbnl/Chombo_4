blt_add_library(NAME EulerOpRK4
    HEADERS EulerOp.H EulerRK4.H
    SOURCES EulerOp.cpp EulerRK4.cpp
    INCLUDES ${CMAKE_CURRENT_SOURCE_DIR}
    DEPENDS_ON Headers common)
