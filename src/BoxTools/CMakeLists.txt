blt_add_library(NAME BoxTools SOURCES Chombo_AMRIO.cpp Chombo_BRMeshRefine.cpp
    Chombo_BaseFab.cpp Chombo_BitSet.cpp Chombo_BoundaryIterator.cpp
    Chombo_Box.cpp Chombo_BoxIterator.cpp Chombo_BoxLayout.cpp
    Chombo_BoxLayoutData.cpp Chombo_Copier.cpp Chombo_DataIndex.cpp
    Chombo_DataIterator.cpp Chombo_DebugOut.cpp Chombo_DenseIntVectSet.cpp
    Chombo_DisjointBoxLayout.cpp Chombo_FArrayBox.cpp Chombo_IntVect.cpp
    Chombo_IntVectSet.cpp Chombo_LayoutIterator.cpp Chombo_LoadBalance.cpp
    Chombo_MeshRefine.cpp Chombo_NeighborIterator.cpp Chombo_Pool.cpp
    Chombo_ProblemDomain.cpp Chombo_ProtoInterface.cpp
    Chombo_ReductionCopier.cpp Chombo_SpreadingCopier.cpp Chombo_Stencil.cpp
    Chombo_TimedDataIterator.cpp Chombo_TreeIntVectSet.cpp Chombo_CH_HDF5.cpp
    HEADERS
    Chombo_AMRIO.H Chombo_BFM.H Chombo_BRMeshRefine.H Chombo_BaseFab.H
    Chombo_BaseFabImplem.H Chombo_BitSet.H Chombo_BoundaryIterator.H
    Chombo_Box.H Chombo_BoxIterator.H Chombo_BoxLayout.H Chombo_BoxLayoutData.H
    Chombo_BoxLayoutDataI.H Chombo_CH_HDF5.H Chombo_CH_OpenMP.H Chombo_Copier.H
    Chombo_DataIndex.H Chombo_DataIterator.H Chombo_DebugDump.H
    Chombo_DebugOut.H Chombo_DenseIntVectSet.H Chombo_DisjointBoxLayout.H
    Chombo_FABView.H Chombo_FArrayBox.H Chombo_HDF5Portable.H Chombo_IntVect.H
    Chombo_IntVectSet.H Chombo_LayoutData.H Chombo_LayoutDataI.H
    Chombo_LayoutIterator.H Chombo_LevelBoxData.H Chombo_LevelData.H
    Chombo_LevelDataI.H Chombo_LoadBalance.H Chombo_MeshRefine.H
    Chombo_NeighborIterator.H Chombo_Pool.H Chombo_ProblemDomain.H
    Chombo_ProtoInterface.H Chombo_RealTensor.H Chombo_RealVect.H
    Chombo_ReductionCopier.H Chombo_SPACE.H Chombo_SliceSpec.H
    Chombo_SpreadingCopier.H Chombo_Stencil.H Chombo_TimedDataIterator.H
    Chombo_TreeIntVectSet.H
    INCLUDES ${CMAKE_CURRENT_SOURCE_DIR}
    DEPENDS_ON BaseTools EBProto Proto chombo_compiler_flags ${LIB_DEP})
if(ENABLE_HDF5)
    target_link_libraries(BoxTools PUBLIC HDF5::HDF5)
endif()
