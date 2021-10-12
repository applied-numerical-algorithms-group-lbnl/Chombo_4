blt_add_library(NAME hoeb
    HEADERS
    Hoeb_ExactSolutions.H
    Hoeb_FortranNameMacro.H
    Hoeb_LAPACKMatrix.H
    Hoeb_Lapack.H
    Hoeb_Neighborhood.H
    Hoeb_Utilities.H
    Hoeb_UtilitiesImplem.H
    SOURCES
    Hoeb_LAPACKMatrix.cpp
    Hoeb_Utilities.cpp
    DEPENDS_ON BaseTools EBTools)
