blt_add_executable(NAME applyHelmholtz SOURCES applyHelmholtz.cpp
    INCLUDES ${CMAKE_CURRENT_SOURCE_DIR} DEPENDS_ON Headers gtest Proto ${LIB_DEP} CUDA::nvToolsExt)
blt_add_test(NAME applyHelmholtz COMMAND applyHelmholtz)
blt_add_executable(NAME applyHelmholtzStressCase
    SOURCES applyHelmholtzStressCase.cpp INCLUDES ${CMAKE_CURRENT_SOURCE_DIR}
    DEPENDS_ON Headers gtest Proto ${LIB_DEP} CUDA::nvToolsExt)
blt_add_test(NAME applyHelmholtzStressCase COMMAND applyHelmholtzStressCase)
blt_add_executable(NAME exchangeTest SOURCES exchangeTest.cpp
    INCLUDES ${CMAKE_CURRENT_SOURCE_DIR} DEPENDS_ON Headers gtest Proto ${LIB_DEP} CUDA::nvToolsExt)
blt_add_test(NAME exchangeTest COMMAND exchangeTest)
blt_add_executable(NAME forallTest SOURCES forallTest.cpp
    INCLUDES ${CMAKE_CURRENT_SOURCE_DIR} DEPENDS_ON Headers gtest Proto ${LIB_DEP} CUDA::nvToolsExt)
blt_add_test(NAME forallTest COMMAND forallTest)
blt_add_executable(NAME func_mess SOURCES func_mess.cpp
    INCLUDES ${CMAKE_CURRENT_SOURCE_DIR} DEPENDS_ON Headers gtest Proto ${LIB_DEP} CUDA::nvToolsExt)
blt_add_test(NAME func_mess COMMAND func_mess)
blt_add_executable(NAME helmholtzExact SOURCES helmholtzExact.cpp
    INCLUDES ${CMAKE_CURRENT_SOURCE_DIR} DEPENDS_ON Headers gtest Proto ${LIB_DEP} CUDA::nvToolsExt)
blt_add_test(NAME helmholtzExact COMMAND helmholtzExact)
blt_add_executable(NAME iota SOURCES iota.cpp
    INCLUDES ${CMAKE_CURRENT_SOURCE_DIR} DEPENDS_ON Headers gtest ${LIB_DEP} CUDA::nvToolsExt)
blt_add_test(NAME iota COMMAND iota)
#blt_add_executable(NAME testLibEBForAll SOURCES testLibEBForAll.cpp
#    INCLUDES ${CMAKE_CURRENT_SOURCE_DIR} DEPENDS_ON Headers gtest ${LIB_DEP} CUDA::nvToolsExt)
#blt_add_test(NAME testLibEBForAll COMMAND testLibEBForAll)
blt_add_executable(NAME testSuite SOURCES testSuite.cpp
    INCLUDES ${CMAKE_CURRENT_SOURCE_DIR} DEPENDS_ON Headers gtest Proto ${LIB_DEP} CUDA::nvToolsExt)
blt_add_test(NAME testSuite COMMAND testSuite)
if(ENABLE_CUDA)
    blt_add_executable(NAME erfSuite SOURCES erfSuite.cpp
        INCLUDES ${CMAKE_CURRENT_SOURCE_DIR} DEPENDS_ON Headers gtest ${LIB_DEP} CUDA::nvToolsExt)
    blt_add_test(NAME erfSuite COMMAND erfSuite)
endif()
