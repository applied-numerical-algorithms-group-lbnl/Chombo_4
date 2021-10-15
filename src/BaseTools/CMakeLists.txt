blt_add_library(NAME BaseTools
    HEADERS
    Chombo_CHOMBO_VERSION.H
    Chombo_CH_System.H
    Chombo_CH_Thread.H
    Chombo_CH_Timer.H
    Chombo_CH_assert.H
    Chombo_CH_config.H
    Chombo_CONSTANTS.H
    Chombo_FortranNameMacro.H
    Chombo_Interval.H
    Chombo_LoHiSide.H
    Chombo_MayDay.H
    Chombo_NamespaceHeader.H
    Chombo_NamespaceFooter.H
    Chombo_ParmParse.H
    Chombo_REAL.H
    Chombo_RefCountedPtr.H
    Chombo_SPMD.H
    Chombo_SPMDI.H
    Chombo_Tuple.H
    Chombo_Vector.H
    Chombo_parstream.H
    SOURCES
    Chombo_CH_System.cpp
    Chombo_CH_Thread.cpp
    Chombo_Interval.cpp
    Chombo_MayDay.cpp
    Chombo_ParmParse.cpp
    Chombo_SPMD.cpp
    Chombo_parstream.cpp
    INCLUDES ${CMAKE_CURRENT_SOURCE_DIR}
    DEPENDS_ON Proto chombo_compiler_flags ${LIB_DEP})
