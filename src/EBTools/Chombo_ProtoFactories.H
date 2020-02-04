#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif


#ifndef _ProtoFactories_H_
#define _ProtoFactories_H_


#include "Proto.H"
#include "EBProto.H"
#include "Proto_EBGraph.H"
#include "Chombo_EBChombo.H"
#include "Chombo_LevelData.H"
#include "Chombo_ProtoInterface.H"
#include "Chombo_NamespaceHeader.H"

/// factory where constructor called is the null constructor
template <typename DataHolder >
class NullConstructorDataFactory: public DataFactory<DataHolder>
{
public:

  virtual ~NullConstructorDataFactory()
  {
  }

  NullConstructorDataFactory()
  {
  }


  /// factory function.  creates a new 'T' object
  /** creates a new 'T' object and returns a pointer to it.
   */
  virtual DataHolder* create(const Box& box, int ncomps, const DataIndex& a_datInd) const
  {
    return new DataHolder();
  }

  virtual bool callDelete() const
  {
    return true;
  }

  //return true if copy and construction are thread safe
  virtual bool threadSafe() const
  {
    return false;
  }

};

using Proto::EBGraph;
/// Factory object to data members of a BoxLayoutData container where objects need a graph in the constructor
template <typename DataHolder> 
class GraphConstructorFactory: public DataFactory< DataHolder >
{
public:

  virtual ~GraphConstructorFactory()
  {
  }

  GraphConstructorFactory(shared_ptr<LevelData<EBGraph> > a_graphs)
  {
    m_graphs = a_graphs;
  }

  /// factory function.  creates a new 'T' object
  /** creates a new 'T' object and returns a pointer to it.
   */
  virtual DataHolder * create(const Box& a_box, int ncomps, const DataIndex& a_datInd) const
  {
    bool usestack = false;
    Proto::Box bx = ProtoCh::getProtoBox(a_box);
    return new DataHolder(bx, (*m_graphs)[a_datInd], usestack);
  }

  virtual bool callDelete() const
  {
    return true;
  }

  //return true if copy and construction are thread safe
  virtual bool threadSafe() const
  {
    return false;
  }

private:
  GraphConstructorFactory();
  shared_ptr<LevelData<EBGraph> > m_graphs;
};


template <typename DataHolder,  int NCOMP> 
class AliasConstructorFactory: public DataFactory< DataHolder >
{
public:

  AliasConstructorFactory(LevelData<DataHolder>* a_source, int whichcomp)
  {
    m_source = a_source;
  }
  virtual ~AliasConstructorFactory()
  {
  }

  /// factory function.  creates a new 'T' object
  /** creates a new 'T' object and returns a pointer to it.
   */
  virtual DataHolder * create(const Box& a_box, int ncomps, const DataIndex& a_datInd) const
  {
    DataHolder& vecholder = (*m_source)[a_datInd];
    DataHolder* retval = new DataHolder();
    //aliasing define
    retval->define(vecholder, m_whichcomp);
    return retval;
  }

  virtual bool callDelete() const
  {
    return false;
  }

  //return true if copy and construction are thread safe
  virtual bool threadSafe() const
  {
    return false;
  }

private:
  AliasConstructorFactory();
  LevelData<DataHolder>* m_source;
  unsigned int              m_whichcomp;

};


#include "Chombo_NamespaceFooter.H"
#endif
