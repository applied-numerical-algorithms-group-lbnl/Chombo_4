blt_add_library(NAME EBDarcy HEADERS EBDarcy.H SOURCES EBDarcy.cpp
    INCLUDES ${CMAKE_CURRENT_SOURCE_DIR} DEPENDS_ON Headers common ${LIB_DEP})
