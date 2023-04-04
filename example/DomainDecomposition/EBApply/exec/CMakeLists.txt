blt_add_executable(NAME PoissonStencil SOURCES PoissonStencil.cpp
    INCLUDES ${CMAKE_CURRENT_SOURCE_DIR} DEPENDS_ON Headers common ${LIB_DEP})
