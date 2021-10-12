blt_add_library(NAME MHDOpsteps HEADERS MHDOp.H MHDsteps.H
    Proto_EulerStep.H Proto_RHScalc.H
    SOURCES MHDOp.cpp MHDsteps.cpp
    INCLUDES ${CMAKE_CURRENT_SOURCE_DIR}
    DEPENDS_ON Headers common)
