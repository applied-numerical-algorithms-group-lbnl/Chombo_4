blt_add_executable(
    NAME hoeb_truncation
    SOURCES truncation.cpp
    INCLUDES ${CMAKE_CURRENT_SOURCE_DIR}
    DEPENDS_ON Headers common hoeb)
