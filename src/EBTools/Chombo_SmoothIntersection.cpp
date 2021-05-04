#include "Chombo_SmoothIntersection.H"
#include "Chombo_SmoothAbsoluteValue.H"
#include <iomanip>
#include "Chombo_NamespaceHeader.H"

SmoothIntersection::
SmoothIntersection(const Vector<BaseIF *>& a_impFuncs,
                   const Real            & a_delta)
{
  if (a_impFuncs.size() == 0)
    {
      MayDay::Abort("Construction of SmoothIntersection requires at least one implicit function.");
    }
  m_delta = a_delta;
  // Number of implicit function in intersection
  m_numFuncs = a_impFuncs.size();


  // Vector of implicit function pointers
  m_impFuncs.resize(m_numFuncs, NULL);

  // Make copies of the implicit functions

  for (int ifunc = 0; ifunc < m_numFuncs; ifunc++)
    {
      if (a_impFuncs[ifunc] == NULL)
        {
          MayDay::Error("you sent in a null baseif");
        }
      else
        {
          m_impFuncs[ifunc] = a_impFuncs[ifunc]->newImplicitFunction();
        }
    }
}
////
SmoothIntersection::
~SmoothIntersection()
{
  // Delete all the copies
  for (int ifunc = 0; ifunc < m_numFuncs; ifunc++)
    {
      if (m_impFuncs[ifunc] != NULL)
        {
          delete m_impFuncs[ifunc];
          m_impFuncs[ifunc] = NULL;
        }
    }
}

////
Real           
SmoothIntersection::
smoothMax(const  IntVect & a_deriv,
          const  RealVect& a_point,
          const  int     & a_closestIF,
          const  int     & a_nextClosestIF) const
          
{

  CH_assert(a_closestIF     >= 0);
  CH_assert(a_nextClosestIF >= 0);
  CH_assert(a_closestIF      < m_numFuncs);
  CH_assert(a_nextClosestIF  < m_numFuncs);
  
  const BaseIF* ffunc = m_impFuncs[a_closestIF];
  const BaseIF* gfunc = m_impFuncs[a_nextClosestIF];
  //this gets a smooth approximation of |f(x) - g(x)|
  SmoothAbsoluteValue smoothie(ffunc, gfunc, m_delta);
  Real wval;
  int icase;
  Real fval = ffunc->value(a_point);
  Real gval = gfunc->value(a_point);
  Real fder = ffunc->value(a_deriv, a_point);
  Real gder = gfunc->value(a_deriv, a_point);
  smoothie.getWCase(icase, wval, a_point);
  Real retval = 0;
  //only deal with smoothie if it is in the range of integration.
  if(icase !=  0)
    {
      if(fval > gval)
        {
          retval = fder;
        }
      else
        {
          retval = gder;
        }
    }
  else
    {
      Real absFminusG = smoothie.smoothAbsFMinusG(a_deriv, a_point);
      retval =  0.5*(fder + gder + absFminusG);
    }
  return retval;
}
////
void 
SmoothIntersection::
findClosest(int            & a_closestIF, 
            int            & a_nextClosestIF,
            int            & a_numWithinDelta,
            const RealVect & a_point) const
{
  CH_assert(m_numFuncs > 0);

  Real valueClosest        = -1.0e30;
  Real valueNextClosest    = -1.0e30;
  a_closestIF     = -1;
  a_nextClosestIF = -1;
  bool foundClosest = false;
  for (int ifunc = 0; ifunc < m_numFuncs; ifunc++)
    {
      Real cur;
      cur = m_impFuncs[ifunc]->value(a_point);
      if (cur > valueClosest)
        {
          valueClosest = cur;
          a_closestIF  = ifunc;
          foundClosest = true;
        }
    }
  if(!foundClosest) MayDay::Error("logic error smoothie0");
  //might not be a next closest so have to be careful here
  bool foundNextClosest = false;
  for (int ifunc = 0; ifunc < m_numFuncs; ifunc++)
    {
      if(ifunc != a_closestIF)
        {
          Real cur;
          cur = m_impFuncs[ifunc]->value(a_point);
          if(cur > valueNextClosest)
            {
              valueNextClosest = cur;
              a_nextClosestIF  = ifunc;
              foundNextClosest = true;
            }
        }
    }
  if(!foundNextClosest && (m_numFuncs > 1)) MayDay::Error("logic error smoothie1");

  a_numWithinDelta = 0;
  if(Abs(valueClosest)     < m_delta) a_numWithinDelta++;
  if(Abs(valueNextClosest) < m_delta) a_numWithinDelta++;
    
}


////
Real 
SmoothIntersection::
value(const RealVect& a_point) const
{

  Real retval;
  if(m_numFuncs < 2)
    {
      retval = m_impFuncs[0]->value(a_point);
    }
  else
    {
      // Maximum of the implicit functions values
      int  closestIF;
      int  nextClosestIF;
      int  numWithinDelta;
  

      findClosest(closestIF, 
                  nextClosestIF,
                  numWithinDelta,
                  a_point);


      retval = smoothMax(IntVect::Zero, a_point,
                         closestIF, nextClosestIF);
      

      
    }
  return retval;
}
////
Real 
SmoothIntersection::
derivative(const IntVect & a_deriv,
           const RealVect& a_point) const
{

  Real retval;
  //only need to do the smoothMax thing if there are two functions close enough
  if(m_numFuncs < 2)
    {
      retval = m_impFuncs[0]->value(a_deriv, a_point);
    }
  else
    {

      // Maximum of the implicit functions values
      int  closestIF;
      int  nextClosestIF;
      int  numWithinDelta;
  

      findClosest(closestIF, 
                  nextClosestIF,
                  numWithinDelta,
                  a_point);

      retval = smoothMax(a_deriv, a_point,
                         closestIF, nextClosestIF);

    }
  return retval;
}

////
BaseIF* 
SmoothIntersection::
newImplicitFunction() const
{
  SmoothIntersection* intersectionPtr = new SmoothIntersection(m_impFuncs, m_delta);

  return static_cast<BaseIF*>(intersectionPtr);
}

#include "Chombo_NamespaceFooter.H"
