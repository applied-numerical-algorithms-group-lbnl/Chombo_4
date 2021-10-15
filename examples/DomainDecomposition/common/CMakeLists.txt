blt_add_library(NAME common HEADERS
    BiCGStabSolver.H
    DebugFunctions.H
    EBAdvection.H
    EBAdvectionFunctions.H
    EBCCProjector.H
    EBIBC.H
    EBMACProjector.H
    EBMultigrid.H
    EBMultigridFunctions.H
    EBParabolicIntegrators.H
    EBPetscSolver.H
    EBRelaxSolver.H
    SOURCES
    EBAdvection.cpp
    EBCCProjector.cpp
    EBMACProjector.cpp
    EBMultigrid.cpp
    EBParabolicIntegrators.cpp
    INCLUDES ${CMAKE_CURRENT_SOURCE_DIR}
    DEPENDS_ON EBTools BaseTools ${LIB_DEP}) # PETSc)
