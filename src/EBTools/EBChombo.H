#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif


#ifndef _EBChombo_GeometryService_H_
#define _EBChombo_GeometryService_H_


#include "Proto.H"
#include "EBProto.H"
#include "Proto_EBGraph.H"
#include "implem/Proto_Notation.H"
#include "ProtoInterface.H"
#include "ProtoFactories.H"


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
typedef Proto::Box Bx;
using ProtoCh::getProtoBox;
typedef Proto::IndexTM<int,2>      Iv2;
typedef Proto::IndexTM<int,DIM>    IvgDim;
typedef Proto::IndexTM<int,DIM-1 > IvgLess1;
typedef Proto::IndexTM<double,DIM> RvgDim;

#endif
