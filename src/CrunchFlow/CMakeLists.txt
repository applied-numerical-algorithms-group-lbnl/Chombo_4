if(ENABLE_CUDA)
    list(APPEND ${LIB_DEP} cuda_runtime)
endif()
blt_add_library(NAME CrunchFlow
    HEADERS
    BaseNamespaceFooter.H
    BaseNamespaceHeader.H
    CHOMBO_VERSION.H
    CH_assert.H
    GenericArithmetic.H
    GenericArithmeticI.H
    IndexTM.H
    IndexTMI.H
    IntVect.H
    MayDay.H
    Metaprograms.H
    Misc.H
    NamespaceFooter.H
    NamespaceHeader.H
    NamespaceVar.H
    REAL.H
    SPACE.H
    SPMD.H
    SPMDI.H
    ShapeArray.H
    UsingNamespace.H
    Vector.H
    memtrack.H
    parstream.H
    SOURCES
    AffinityNumerical.cpp
    AqueousToBulkConvert.cpp
    IntVect.cpp
    MayDay.cpp
    SpeciesLocal.cpp
    assemble_local.cpp
    fit.cpp
    fitgamma.cpp
    fx_local.cpp
    gamma.cpp
    jac_local.cpp
    jacmin.cpp
    keqcalc2.cpp
    os3d_newton.cpp
    parstream.cpp
    reaction.cpp
    satcalc.cpp
    totconc.cpp
    INCLUDES ${CMAKE_CURRENT_SOURCE_DIR}
    DEPENDS_ON LAPACK::LAPACK chombo_compiler_flags ${LIB_DEP})
if(ENABLE_CUDA)
    target_link_libraries(CrunchFlow PUBLIC CUDA::cusolver)
endif()
