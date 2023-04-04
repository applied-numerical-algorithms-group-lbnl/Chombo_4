blt_add_executable(NAME EBPETSc SOURCES petscDriver.cpp
    INCLUDES ${CMAKE_CURRENT_SOURCE_DIR} DEPENDS_ON Headers common)
