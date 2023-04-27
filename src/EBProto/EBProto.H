#pragma once

#include "Proto.H"
#include "base/Proto_LoHiSide.H"

namespace Proto
{

  enum CENTERING {CELL=0,BOUNDARY,XFACE,YFACE,ZFACE};

  // A way to loop over each face type
  constexpr std::array<CENTERING, DIM>   REG_FACES{
    XFACE
#if DIM>=2
    ,YFACE
#endif
#if DIM>= 3
    ,ZFACE
#endif
  };
  constexpr std::array<CENTERING, DIM+1> ALL_FACES{BOUNDARY,
    XFACE
#if DIM>=2
    ,YFACE
#endif
#if DIM>= 3
    ,ZFACE
#endif
  };

  inline int FACE_DIR(CENTERING cent)
  {
    int retval = -1;
    if (cent == XFACE) { retval = 0;}
    else if (cent == YFACE) { retval = 1;}
    else if (cent == ZFACE) { retval = 2;}
    return retval;
  }

  /// Return the index into a SpaceDim*2 vector
  /// use case is:
  /// for (int dir = 0; dir != SpaceDim; dir++)
  /// {
  ///    for (SideIterator sit; sit.ok(); ++sit)
  ///      {
  ///        int faceIndex = ebp_index(dir, side);
  ///        val[faceIndex] = side.sign() * flux[dir];
  ///      }
  /// }
  inline int ebp_index(const int a_idir,
                       const Side::LoHiSide a_side)
  {
    PR_assert(a_idir >= 0);
    //CH_assert(a_idir < Chombo4::SpaceDim); // FIXME add once dimensionality be better managed
    return a_idir + a_side*DIM;
  }

  /// Get the side and direction from a vector index
  //  use case is:
  /// for (int faceIdx = 0; faceIdx != 2*SpaceDim; faceIdx++)
  /// {
  ///    int faceDir;
  ///    Side::LoHiSide faceSide
  ///    ebp_face(faceIndx, faceDir, faceSide)
  ///    val[faceIndex] = faceSide.sign() * flux[faceDir];
  /// }
  inline void ebp_face(const int a_index,
                       int& a_idir,
                       Side::LoHiSide& a_side)
  {
    PR_assert(a_index >= 0);
    //CH_assert(a_index < 2*Chombo4::SpaceDim); // FIXME add once dimensionality be better managed
    a_idir = a_index%DIM;
    if (a_index/DIM == 0)
      {
        a_side = Side::Lo;
      }
    else
      {
        a_side = Side::Hi;
      }
    return;
  }
}

#include "implem/Proto_PointSet.H"
#include "Proto_HostBoxData.H"
#include "implem/Proto_HostIrregData.H"
#include "implem/Proto_IrregNode.H"
#include "Proto_EBGraph.H"
//#include "Proto_LevelData.H"

#include "Proto_EBBoxData.H"
#include "Proto_EBFluxData.H"
#include "implem/Proto_AggStencil.H"
#include "Proto_EBStencil.H"
#include "implem/Proto_EBStencilArchive.H"
#include "Proto_SimpleImplicitFunctions.H"
#include "Proto_RealVect.H"
#include "base/Proto_Box.H"
#include "base/Proto_Point.H"
#include "implem/Proto_LSquares.H"
#include "Proto_BaseIF.H"
#include "implem/Proto_IrregNode.H"
#include "implem/Proto_CutCellMoments.H"
#include "Proto_HostBoxData.H"
#include "implem/Proto_HostIrregData.H"
#include "implem/Proto_PointSet.H"
//#include "Proto_LevelData.H"
#include "Proto_EBGraph.H"
#include "Proto_EBStencil.H"

#include "Proto_EBHostData.H"
#include "Proto_EBHostFluxData.H"

//#include "hoeb/Proto_EBStencilHO.H"
