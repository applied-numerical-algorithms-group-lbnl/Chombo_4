blt_add_executable(NAME applyOnly SOURCES applyHelmholtz.cpp
    INCLUDES ${CMAKE_CURRENT_SOURCE_DIR} DEPENDS_ON Headers common)
