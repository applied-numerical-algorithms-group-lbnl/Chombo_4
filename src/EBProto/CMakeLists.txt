blt_add_library(NAME EBProto HEADERS
    EBProto.H
    Proto_BaseIF.H
    Proto_EBBoxData.H
    Proto_EBDictionary.H
    Proto_EBFluxData.H
    Proto_EBGraph.H
    Proto_EBHostData.H
    Proto_EBHostFluxData.H
    Proto_EBStencil.H
    Proto_HostBoxData.H
    Proto_RealVect.H
    Proto_SimpleImplicitFunctions.H
    INCLUDES ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/implem
    DEPENDS_ON chombo_compiler_flags Proto ${LIB_DEP})
