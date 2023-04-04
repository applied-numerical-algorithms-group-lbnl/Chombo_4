blt_add_executable(NAME serialization SOURCES serializationTest.cpp
    INCLUDES ${CMAKE_CURRENT_SOURCE_DIR} DEPENDS_ON Headers common ${LIB_DEP})
