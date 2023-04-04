#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _Chombo_PROTOINTERFACE_H_
#define _Chombo_PROTOINTERFACE_H_

#include "Proto.H"
#include "Chombo_IntVect.H"
#include "Chombo_Box.H"
#include "base/Proto_Point.H"
#include "Chombo_NamespaceHeader.H"

///these functions are meant to interface with the Proto language.
namespace ProtoCh
{
  using ::Proto::Point;
///get point from intvect
  extern Point      getPoint( const IntVect& a_iv);

/// gets proto box from chombo box
  extern ::Proto::Box  getProtoBox( const Box& a_box);

///get intvect from point
  extern IntVect getIntVect(const  Point  & a_pt);

///get chombo box from proto box
  extern Box getBox(const ::Proto::Box & a_bx);

///get chombo box from proto box
  Box getChomboBox(const ::Proto::Box & a_bx);
  
} 

#include "Chombo_NamespaceFooter.H"

#endif
