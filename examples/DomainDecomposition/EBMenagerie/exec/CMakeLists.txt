blt_add_executable(NAME EBMenagerie
    SOURCES drMoreau.cpp
    INCLUDES ${CMAKE_CURRENT_SOURCE_DIR}
    DEPENDS_ON Headers common ${LIB_DEP})
