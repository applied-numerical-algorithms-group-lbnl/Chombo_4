blt_add_library(NAME EBTools HEADERS
    Chombo_EBChombo.H
    Chombo_EBDictionary.H
    Chombo_EBDictionaryImplem.H
    Chombo_EBEncyclopedia.H
    Chombo_EBLevelBoxData.H
    Chombo_EBLevelFluxData.H
    Chombo_GeometryService.H
    Chombo_GeometryServiceImplem.H
    Chombo_ProtoFactories.H
    Chombo_SmoothAbsoluteValue.H
    Chombo_SmoothIntersection.H
    Chombo_SmoothUnion.H
    Chombo_UnionIntersection.H
    SOURCES
    Chombo_SmoothAbsoluteValue.cpp
    Chombo_SmoothIntersection.cpp
    Chombo_SmoothUnion.cpp
    Chombo_UnionIntersection.cpp
    INCLUDES ${CMAKE_CURRENT_SOURCE_DIR}
    DEPENDS_ON BaseTools BoxTools chombo_compiler_flags)
