#pragma once

#include "Proto.H"
#include "Proto_LoHiSide.H"

namespace Proto
{

  enum CENTERING {CELL=0,BOUNDARY,XFACE,YFACE,ZFACE};
/// Return the index into the arc vector
  inline int ebp_index(int a_idir, Side::LoHiSide a_side) 
  {
    return a_idir + a_side*DIM;
  }
/// side and direction from the arc vector index
  inline ebp_face(int a_index, int& a_idir, Side::LoHiSide& a_side)
  {
    a_idir = a_index%2;
    a_side = a_index/DIM;
  }
}

#include "implem/Proto_PointSet.H"
#include "Proto_HostBoxData.H"
#include "implem/Proto_HostIrregData.H"
#include "implem/Proto_IrregNode.H"
#include "Proto_EBGraph.H"
//#include "implem/Proto_LevelData.H"

#include "Proto_EBBoxData.H"
#include "Proto_EBFluxData.H"
#include "implem/Proto_AggStencil.H"
#include "Proto_EBStencil.H"
#include "implem/Proto_EBStencilArchive.H"
#include "Proto_SimpleImplicitFunctions.H"
#include "Proto_RealVect.H"
#include "Proto_Box.H"
#include "Proto_Point.H"
#include "implem/Proto_LSquares.H"
#include "Proto_BaseIF.H"
#include "implem/Proto_IrregNode.H"
#include "implem/Proto_CutCellMoments.H"
#include "Proto_HostBoxData.H"
#include "implem/Proto_HostIrregData.H"
#include "implem/Proto_PointSet.H"
//#include "implem/Proto_LevelData.H"
#include "Proto_EBGraph.H"
#include "Proto_EBStencil.H"

#include "Proto_EBHostData.H"
#include "Proto_EBHostFluxData.H"

#include "hoeb/Proto_EBStencilHO.H"
