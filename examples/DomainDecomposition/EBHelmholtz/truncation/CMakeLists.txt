blt_add_executable(NAME truncation
    SOURCES truncation.cpp
    INCLUDES ${CMAKE_CURRENT_SOURCE_DIR}
    DEPENDS_ON Headers common)
