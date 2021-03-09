#include "Chombo_UnionIntersection.H"

#include "Chombo_NamespaceHeader.H"
using Proto::IndexTM;
using Proto::BaseIF;
///union code
///You can't scare me, I'm sticking with the union.
UnionIF::
UnionIF(const std::vector<BaseIF *>& a_impFuncs)
{
  // Number of implicit function in union
  m_numFuncs = a_impFuncs.size();

  // Vector of implicit function pointers
  m_impFuncs.resize(m_numFuncs, NULL);

  // Make copies of the implicit functions
  for (int ifunc = 0; ifunc < m_numFuncs; ifunc++)
  {
    if (a_impFuncs[ifunc] == NULL)
    {
      m_impFuncs[ifunc] = NULL;
    }
    else
    {
      m_impFuncs[ifunc] = a_impFuncs[ifunc]->newImplicitFunction();
    }
  }
}

UnionIF::~UnionIF()
{
  // Delete all the copies
  for (int ifunc = 0; ifunc < m_numFuncs; ifunc++)
  {
    if (m_impFuncs[ifunc] != NULL)
    {
      delete m_impFuncs[ifunc];
    }
  }
}


double
UnionIF::
value(const IndexTM<int,   DIM> & a_partialDerivative,
      const IndexTM<double,DIM>& a_point) const
{
  int closestIF = -1;
  findClosest(a_point,closestIF);

  if (closestIF == -1)
  {
    if (a_partialDerivative.sum() == 0)
    {
      return 1.0;
    }
    else
    {
      return 0.0;
    }
  }
  else
  {
    //this is the weirdest infinite loop if it unions all the way down.
    //not likely enough to test for. --dtg
    return m_impFuncs[closestIF]->value(a_partialDerivative,a_point);
  }
}

void UnionIF::findClosest(const IndexTM<double,DIM> & a_point,
                          int                       & a_closestIF) const
{
  double retval = 0.0;

  if (m_numFuncs > 0)
    {
      IndexTM<int, DIM> deriv = IndexTM<int, DIM>::Zero;
      retval = m_impFuncs[0]->value(deriv,a_point);
      a_closestIF = 0;

      for (int ifunc = 1; ifunc < m_numFuncs; ifunc++)
        {
          double cur;
          cur = m_impFuncs[ifunc]->value(deriv, a_point);
          if (cur < retval)
            {
              retval = cur;
              a_closestIF = ifunc;
            }
        }
    }
}

BaseIF* UnionIF::newImplicitFunction() const
{
  UnionIF* unionPtr = new UnionIF(m_impFuncs);

  return static_cast<BaseIF*>(unionPtr);
}



IntersectionIF::
IntersectionIF(const std::vector<BaseIF *>& a_impFuncs)
{
  // Number of implicit function in intersection
  m_numFuncs = a_impFuncs.size();

  // Vector of implicit function pointers
  m_impFuncs.resize(m_numFuncs);

  // Make copies of the implicit functions

  for (int ifunc = 0; ifunc < m_numFuncs; ifunc++)
  {
    if (a_impFuncs[ifunc] == NULL)
    {
      m_impFuncs[ifunc] = NULL;
    }
    else
    {
      m_impFuncs[ifunc] = a_impFuncs[ifunc]->newImplicitFunction();
    }
  }
}


IntersectionIF::~IntersectionIF()
{
  // Delete all the copies
  for (int ifunc = 0; ifunc < m_numFuncs; ifunc++)
  {
    if (m_impFuncs[ifunc] != NULL)
    {
      delete m_impFuncs[ifunc];
    }
  }
}


double
IntersectionIF::
value(const IndexTM<int ,  DIM> & a_partialDerivative,
      const IndexTM<double,DIM> & a_point) const
{
  int closestIF = -1;
  findClosest(a_point,closestIF);

  if (closestIF == -1)
  {
    if (a_partialDerivative.sum() == 0)
    {
      return -1.0;
    }
    else
    {
      return 0.0;
    }
  }
  else
  {
    return m_impFuncs[closestIF]->value(a_partialDerivative,a_point);
  }
}

BaseIF* IntersectionIF::newImplicitFunction() const
{
  IntersectionIF* intersectionPtr = new IntersectionIF(m_impFuncs);

  return static_cast<BaseIF*>(intersectionPtr);
}

void IntersectionIF::findClosest(const IndexTM<double, DIM> & a_point,
                                 int                        & a_closestIF) const
{
  double retval = 0.0;
  
  IndexTM<int, DIM> deriv = IndexTM<int, DIM>::Zero;
  if (m_numFuncs > 0)
    {
      retval = m_impFuncs[0]->value(deriv, a_point);
      a_closestIF = 0;

      for (int ifunc = 1; ifunc < m_numFuncs; ifunc++)
        {
          double cur;
          cur = m_impFuncs[ifunc]->value(deriv, a_point);
          if (cur > retval)
            {
              retval = cur;
              a_closestIF = ifunc;
            }
        }
    }
}

#include "Chombo_NamespaceFooter.H"
