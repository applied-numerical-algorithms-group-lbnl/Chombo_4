#!/bin/bash

LIST="Chombo_AMRIO.H               Chombo_Box.cpp            Chombo_DataIndex.H            Chombo_FArrayBox.cpp       Chombo_LevelData.H           Chombo_ProtoInterface.H     Chombo_TimedDataIterator.cpp
Chombo_AMRIO.cpp             Chombo_BoxIterator.H      Chombo_DataIndex.cpp          Chombo_FluxBox.H           Chombo_LevelData.cpp         Chombo_ProtoInterface.cpp   Chombo_TreeIntVectSet.H
Chombo_BFM.H                 Chombo_BoxIterator.cpp    Chombo_DataIterator.H         Chombo_FluxBox.cpp         Chombo_LevelDataI.H          Chombo_RealTensor.H         Chombo_TreeIntVectSet.cpp
Chombo_BRMeshRefine.H        Chombo_BoxLayout.H        Chombo_DataIterator.cpp       Chombo_HDF5Portable.H      Chombo_LoadBalance.H         Chombo_RealVect.H           GNUmakefile
Chombo_BRMeshRefine.cpp      Chombo_BoxLayout.cpp      Chombo_DebugDump.H            Chombo_IntVect.H           Chombo_LoadBalance.cpp       Chombo_ReductionCopier.H    README
Chombo_BaseFab.H             Chombo_BoxLayoutData.H    Chombo_DebugOut.H             Chombo_IntVect.cpp         Chombo_MeshRefine.H          Chombo_ReductionCopier.cpp  multidim
Chombo_BaseFab.cpp           Chombo_BoxLayoutData.cpp  Chombo_DebugOut.cpp           Chombo_IntVectSet.H        Chombo_MeshRefine.cpp        Chombo_SPACE.H              myScriptCudaToProto.sh
Chombo_BaseFabImplem.H       Chombo_BoxLayoutDataI.H   Chombo_DenseIntVectSet.H      Chombo_IntVectSet.cpp      Chombo_NeighborIterator.H    Chombo_SliceSpec.H          script2.sh
Chombo_BitSet.H              Chombo_CH_HDF5.H          Chombo_DenseIntVectSet.cpp    Chombo_LayoutData.H        Chombo_NeighborIterator.cpp  Chombo_SpreadingCopier.H
Chombo_BitSet.cpp            Chombo_CH_HDF5.cpp        Chombo_DisjointBoxLayout.H    Chombo_LayoutDataI.H       Chombo_Pool.H                Chombo_SpreadingCopier.cpp
Chombo_BoundaryIterator.H    Chombo_CH_OpenMP.H        Chombo_DisjointBoxLayout.cpp  Chombo_LayoutIterator.H    Chombo_Pool.cpp              Chombo_Stencil.H
Chombo_BoundaryIterator.cpp  Chombo_Copier.H           Chombo_FABView.H              Chombo_LayoutIterator.cpp  Chombo_ProblemDomain.H       Chombo_Stencil.cpp
Chombo_Box.H                 Chombo_Copier.cpp         Chombo_FArrayBox.H            Chombo_LevelBoxData.H      Chombo_ProblemDomain.cpp     Chombo_TimedDataIterator.Hi"


for var in $LIST
do 
	sh myScriptCudaToProto.sh $var
done
