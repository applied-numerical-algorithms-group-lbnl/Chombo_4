#ifndef _Chombo_UNIONIntersection_H_
#define _Chombo_UNIONIntersection_H_

#include "Proto_RealVect.H"
#include "Proto_BaseIF.H"
#include "implem/Proto_IndexTM.H"

#include "Chombo_NamespaceHeader.H"

///
/**
    This implicit function is the minimum of two (of more) implicit
    functions with which it is constructed.  This results in an implicit
    function whose interior is the union of the interiors of the given
    implicit functions.
 */
class UnionIF: public Proto::BaseIF
{
public:

  ///
  /**
     Return the value of the function at a_point.  When delineating a domain,
     the level set value=0 represents the boundary and value<0 is inside the
     fluid.
  */
  virtual double value(const Proto::IndexTM<int   ,DIM> & a_partialDerivative,
                       const Proto::IndexTM<double,DIM> & a_point) const;

    
  /// factory method
  /**
     Return a newly allocated derived class.  The responsibility
     for deleting the memory is left to the calling function.
  */
  virtual BaseIF* newImplicitFunction() const;

  ///
  /**
     Constructor specifying any number of implicit functions as inputs.
  */
  UnionIF(const std::vector<BaseIF*>& a_impFuncs);

  /// Destructor
  virtual ~UnionIF();


  void findClosest(const Proto::IndexTM<double, DIM>& a_point, int& closestIF) const;

protected:
  int                  m_numFuncs; // number of implicit functions
  std::vector<BaseIF*> m_impFuncs; // implicit functions

private:
  UnionIF();
  UnionIF(const UnionIF& a_inputIF);
  void operator=(const UnionIF& a_inputIF);

};



///
/**
    This implicit function is the maximum of two (of more) implicit
    functions with which it is constructed.  This results in an implicit
    function whose interior is the intersection of the interiors of the given
    implicit functions.
 */
class IntersectionIF: public Proto::BaseIF
{
public:

  ///
  /**
      Constructor specifying any number of implicit functions as inputs.
   */
  IntersectionIF(const std::vector<BaseIF*>& a_impFuncs);



  /// Destructor
  virtual ~IntersectionIF();

  virtual double value(const Proto::IndexTM<int,   DIM> & a_partialDerivative,
                       const Proto::IndexTM<double,DIM> & a_point) const;

  virtual BaseIF* newImplicitFunction() const;

  void findClosest(const Proto::IndexTM<double,DIM>& a_point,int& closestIF) const;


protected:
  int             m_numFuncs; // number of implicit functions
  std::vector<BaseIF*> m_impFuncs; // implicit functions

private:
  /// Copy constructor
  IntersectionIF(const IntersectionIF& a_inputIF);
  IntersectionIF();
  void operator=(const IntersectionIF& a_inputIF);

};
///
/**
    This implicit function is the negative of the implicit function with which
    it is constructed.  This results in an implicit function whose interior is
    the complement of the interior of the given implicit function (minus the
    the surface).
 */
class ComplementIF: public Proto::BaseIF
{
public:
  ///
  /**
      Constructor specifying the implicit function and whether to complement
   */
  ComplementIF(const BaseIF& a_impFunc);

  /// Destructor
  virtual ~ComplementIF();

  virtual double value(const Proto::IndexTM<int,   DIM> & a_partialDerivative,
                       const Proto::IndexTM<double,DIM> & a_point) const;

  virtual BaseIF* newImplicitFunction() const;


protected:

  BaseIF* m_impFunc; // implicit function to complement

private:
  ComplementIF(const ComplementIF& a_inputIF);
  ComplementIF();
  void operator=(const ComplementIF& a_inputIF);
};

#include "Chombo_NamespaceFooter.H"
#endif
