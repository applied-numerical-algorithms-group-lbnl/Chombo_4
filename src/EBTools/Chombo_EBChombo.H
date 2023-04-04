#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif


#ifndef _Chombo_EBChombo_H_
#define _Chombo_EBChombo_H_


#include "Proto.H"
#include "EBProto.H"
#include "Proto_EBGraph.H"
#include "implem/Proto_Notation.H"
#include "Chombo_ProtoInterface.H"
#include "Chombo_ProtoFactories.H"
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
#include "Proto_EBBoxData.H"
#include "implem/Proto_MinimalCCCM.H"


using std::shared_ptr;
using std::vector;
using Proto::IndexedMoments;
using Proto::MomentIterator;
using Proto::IrregNode;
using Proto::CENTERING;
using Proto::CELL;
using Proto::BOUNDARY;
using Proto::XFACE;
using Proto::YFACE;
using Proto::ZFACE;
using Proto::BaseIF;
using Proto::HostIrregData;
using Proto::HostBoxData;
using Proto::EBHostData;
using Proto::EBBoxData;
using Proto::EBFluxData;
using Proto::EBHostFluxData;
using Proto::BoxData;
using Proto::EBGraph;
using Proto::IndexTM;
using Proto::PointSet;
using Proto::PointSetIterator;
using Proto::CutCellMoments;
using Proto::IrregGraph;
using Proto::EBIndex;
using Proto::IFData;
using Proto::MinimalCCCM;
using Proto::NormalDerivativeNew;
using Proto::IFSlicer;
using Proto::EBStencil;
using Proto::Var;
using Proto::AggStencil;
using Proto::LocalStencil;
using Proto::RealVect;
using Proto::SimpleEllipsoidIF;
using Proto::Stencil;
using Proto::EBStencilArchive;
using Chombo4::DisjointBoxLayout;
using Chombo4::DataIterator;
using Chombo4::DataIndex;

using ProtoCh::getProtoBox;
typedef Proto::Box Bx;
typedef Chombo4::Box          ChBx;
typedef Chombo4::DataIterator ChDit;

typedef Proto::IndexTM<int,2>      Iv2;
typedef Proto::IndexTM<int,DIM>    IvgDim;
typedef Proto::IndexTM<int,DIM-1 > IvgLess1;
typedef Proto::IndexTM<double,DIM> RvgDim;

#endif
