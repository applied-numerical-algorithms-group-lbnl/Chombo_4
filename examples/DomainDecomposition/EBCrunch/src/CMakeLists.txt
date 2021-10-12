blt_add_executable(NAME EBCrunch
    SOURCES initialize.cpp
    INCLUDES ${CMAKE_CURRENT_SOURCE_DIR}
    DEPENDS_ON Headers common)
