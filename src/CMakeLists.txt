add_subdirectory(BaseTools)
add_subdirectory(EBProto)
add_subdirectory(BoxTools)
add_subdirectory(CrunchFlow)
add_subdirectory(EBTools)
blt_add_library(NAME Headers
    HEADERS
    BaseTools/Chombo_BaseNamespaceHeader.H
    BaseTools/Chombo_BaseNamespaceFooter.H
    DEPENDS_ON BaseTools BoxTools CrunchFlow EBProto EBTools)
