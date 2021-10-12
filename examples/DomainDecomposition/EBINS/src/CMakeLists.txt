blt_add_library(NAME EBINS
    HEADERS
    BCGVelAdvect.H
    CrunchInterface.H
    EBINS.H
    SOURCES
    BCGVelAdvect.cpp
    CrunchInterface.cpp
    EBINS.cpp
    INCLUDES ${CMAKE_CURRENT_SOURCE_DIR}
    DEPENDS_ON Headers common)
