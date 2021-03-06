
#ifndef _Proto_MINIMALCCCMIMPLEM_H_
#define _Proto_MINIMALCCCMIMPLEM_H_

#include <cmath>
#include <cstdlib>
#include <fstream>
#include <string>


#include "Proto_LSProblem.H"
#include "Proto_Factorial.H"

namespace Proto
{
// Null constructor
  template <int dim> MinimalCCCM<dim>::MinimalCCCM()
  {
  }

// Copy constructor
  template <int dim> MinimalCCCM<dim>::MinimalCCCM(const MinimalCCCM<dim>& a_MinimalCCCM)
    :m_cutCellMoments(a_MinimalCCCM.m_cutCellMoments)
  {
  }

// Use this for initializing
  template <int dim> MinimalCCCM<dim>::MinimalCCCM(const IFData<dim>& a_info)
    :m_cutCellMoments(a_info)
  {
    m_cutCellMoments.m_bdCCOn = false;

    for (int hilo = 0; hilo < 2; ++hilo)
    {
      for (int idir = 0; idir < dim; ++idir)
      {
        // Identifier for which boundary cutCellMoment
        Iv2 bdId;
        bdId[BDID_DIR] = idir;
        bdId[BDID_HILO] = hilo;

        IFData<dim-1> reducedInfo(a_info,a_info.m_maxOrder+1,idir,hilo);

        CutCellMoments<dim-1>bdCutCellMoments(reducedInfo);

        m_cutCellMoments.m_bdCutCellMoments[bdId] = bdCutCellMoments;


        // Notice whether at least one lower dimensional cutCell is on the
        // interface
        if (reducedInfo.m_allVerticesOn)
        {
          m_cutCellMoments.m_bdCCOn = true;
        }
      }
    }
  }

// Destructor
  template <int dim> MinimalCCCM<dim>::~MinimalCCCM()
  {
  }


// Solve
  template <int dim> void MinimalCCCM<dim>::computeMoments(const int                & a_orderPmax,
                                                           const int                & a_degreePmax)
  {
    PR_TIME("computemoments");
    int zerothOrderOfAccuracy = 0;
    int highestDegree         = a_orderPmax + a_degreePmax;

    vector<double> RNorm(3);
    for (int i = 0; i < 3; i++)
    {
      RNorm[i] = LARGEREALVAL;
    }

    for (int i = 0; i <= highestDegree; i++)
    {
      m_cutCellMoments.m_residual.push_back(RNorm);
    }

    m_boundaryMomentsComputed = false;

    computeMomentsRecursively(zerothOrderOfAccuracy,
                              highestDegree);

    IndexTM<int,dim> refineInDir;

    // Move moments from local coord to parent coord
    IndexTM<double,dim> delta = m_cutCellMoments.m_IFData.m_parentCoord.m_origin;
    delta                  -= m_cutCellMoments.m_IFData.m_localCoord .m_origin;

    PthMoment copyMoments = m_cutCellMoments.m_moments;
    for (typename PthMoment::const_iterator it = copyMoments.begin();
         it != copyMoments.end(); ++it)
    {
      IvDim mono = it->first;
      m_cutCellMoments.m_moments[mono] = m_cutCellMoments.changeMomentCoordinates(copyMoments, mono, delta);
    }

    PthMoment copyEBMoments = m_cutCellMoments.m_EBmoments;
    for (typename PthMoment::const_iterator it = copyEBMoments.begin(); it != copyEBMoments.end(); ++it)
    {
      IvDim mono = it->first;
      m_cutCellMoments.m_EBmoments[mono] = m_cutCellMoments.changeMomentCoordinates(copyEBMoments, mono, delta);
    }

    // Move bdCutCell moments from parent coord to cell center coord. From here on bdCutCell moments will
    // not be used in any least squares problem

    for (typename BdCutCellMoments::iterator it = m_cutCellMoments.m_bdCutCellMoments.begin(); it != m_cutCellMoments.m_bdCutCellMoments.end(); ++it)
    {
      it->second.changeMomentCoordinatesToCellCenter();
    }
  }

  template <int dim> void MinimalCCCM<dim>::computeMomentsRecursively(const int                & a_orderPmax,
                                                                      const int                & a_degreePmax)
  {
    PR_TIME("computemomentsRecursively");
    //pout() << "a_orderPmax  = " << a_orderPmax << endl;
    //pout() << "a_degreePmax = " << a_degreePmax << endl;
    //pout() << "m_cutCellMoments.m_IFData.m_maxOrder " << m_cutCellMoments.m_IFData.m_maxOrder  << endl;
    PR_assert(m_cutCellMoments.m_IFData.m_maxOrder >= a_orderPmax);

    if (m_cutCellMoments.m_IFData.m_allVerticesOut)
    {
      for (int iOrder = 0; iOrder <= a_orderPmax; ++iOrder)
      {
        LSProblem<dim> lsp(iOrder + a_degreePmax, false);

        // Fill moments
        const PthMomentLoc monDegreeP = lsp.getMonomialLocMapDegreeP();
        for (typename PthMomentLoc::const_iterator it = monDegreeP.begin(); it != monDegreeP.end(); ++it)
        {
          m_cutCellMoments.m_EBmoments[it->first] = 0.0;
        }

        const PthMomentLoc monDegreePLess1 = lsp.getMonomialLocMapDegreePLess1();
        for (typename PthMomentLoc::const_iterator it = monDegreePLess1.begin(); it != monDegreePLess1.end(); ++it)
        {
          m_cutCellMoments.m_moments[it->first] = 0.0;
        }
      }
    }
    else if (m_cutCellMoments.m_IFData.m_allVerticesIn && !m_cutCellMoments.m_bdCCOn)
    {
      for (int iOrder = 0; iOrder <= a_orderPmax; ++iOrder)
      {
        LSProblem<dim> lsp(iOrder + a_degreePmax, false);

        // Fill moments of degree P and P-1
        const PthMomentLoc monDegreeP = lsp.getMonomialLocMapDegreeP();
        for (typename PthMomentLoc::const_iterator it = monDegreeP.begin(); it != monDegreeP.end(); ++it)
        {
          m_cutCellMoments.m_EBmoments[it->first] = 0.0;
        }

        const PthMomentLoc monDegreePLess1 = lsp.getMonomialLocMapDegreePLess1();
        for (typename PthMomentLoc::const_iterator it = monDegreePLess1.begin(); it != monDegreePLess1.end(); ++it)
        {
          m_cutCellMoments.m_moments[it->first] = m_cutCellMoments.fullCellQuadrature(it->first,m_cutCellMoments.m_IFData.m_parentCoord);
        }
      }
    }
    else
    {
      // Only compute the boundary moments if they haven't already been
      // computed (earlier in the recursion)
      if (!m_boundaryMomentsComputed)
      {
        for (typename BdCutCellMoments::iterator it = m_cutCellMoments.m_bdCutCellMoments.begin(); it != m_cutCellMoments.m_bdCutCellMoments.end(); ++it)
        {
          MinimalCCCM<dim-1> subProblem(it->second.m_IFData);
          subProblem.computeMoments(a_orderPmax,a_degreePmax+1);
          it->second = subProblem.m_cutCellMoments;
        }

        m_boundaryMomentsComputed = true;
      }

      // Make a LSProb
      IvDim zeroDerivative = IndexTM<int,dim>::Zero;
      LSProblem<dim> lsp(a_orderPmax,a_degreePmax, false,m_cutCellMoments.m_IFData.m_normalDerivatives[zeroDerivative]);

      vector<double> rhs = computeRhs(lsp,a_orderPmax);

      // Solve the problem and return residual
      int lsCode=lsp.invertNormalEq(rhs,m_cutCellMoments.m_residual[a_degreePmax]);
      PR_assert(lsCode == 0);


      // Fill moments
      const PthMomentLoc monDegreePLess1 = lsp.getMonomialLocMapDegreePLess1();
      for (typename PthMomentLoc::const_iterator it = monDegreePLess1.begin();
           it != monDegreePLess1.end(); ++it)
      {
        m_cutCellMoments.m_moments[it->first] = lsp.getUnknown(it->second + lsp.getNumberDegP());
      }

      const PthMomentLoc monDegreeP = lsp.getMonomialLocMapDegreeP();
      for (typename PthMomentLoc::const_iterator it = monDegreeP.begin(); it != monDegreeP.end(); ++it)
      {
        m_cutCellMoments.m_EBmoments[it->first] = lsp.getUnknown(it->second);
      }

    }

    if (a_degreePmax > 0)
    {
      computeMomentsRecursively(a_orderPmax + 1,
                                a_degreePmax - 1);
    }
  }

  template <int dim> vector<double> MinimalCCCM<dim>::computeRhs(LSProblem<dim> & a_lsp,
                                                               const int      & a_orderPmax)
  {
    PR_TIME("computeRHS");
    // Resize rhs
    int numEq = dim*a_lsp.getNumberDegP();
    vector<double> rhs(numEq);

    // For each moment iterate thru bd CutCellMoments incrementing same comp
    // of rhs
    const LocPthMoment& locMap = a_lsp.getLocMonomialMapDegreeP();
    for (typename LocPthMoment::const_iterator it = locMap.begin(); it != locMap.end(); ++it)
    {
      int jth = it->first;
      IvDim mono = it->second;

      int hiSide = 1;
      int loSide = 0;
      Iv2 bdId;

      for (int idir = 0; idir < dim; ++idir)
      {
        // Which lower dimensional monomial corresponds (mono,j)
        IndexTM<int,dim-1> mono1Less;
        for (int jdir = 0; jdir < dim; ++jdir)
        {
          if (jdir < idir)
          {
            mono1Less[jdir] = mono[jdir];
          }
          else if (jdir > idir)
          {
            mono1Less[jdir-1] = mono[jdir];
          }
        }

        bdId[0] = idir;
        bdId[1] = hiSide;

        double hiMom = m_cutCellMoments.m_bdCutCellMoments[bdId].m_moments[mono1Less];

        bdId[1] = loSide;

        double loMom = m_cutCellMoments.m_bdCutCellMoments[bdId].m_moments[mono1Less];
        int exponent = it->second[idir];

        double loSideValue;
        double hiSideValue;

        loSideValue = m_cutCellMoments.m_IFData.m_localCoord.convertDir(-0.5*m_cutCellMoments.m_IFData.m_cellCenterCoord.m_dx[idir],
                                                                        m_cutCellMoments.m_IFData.m_cellCenterCoord,
                                                                        idir);

        hiSideValue = m_cutCellMoments.m_IFData.m_localCoord.convertDir( 0.5*m_cutCellMoments.m_IFData.m_cellCenterCoord.m_dx[idir],
                                                                         m_cutCellMoments.m_IFData.m_cellCenterCoord,
                                                                         idir);

        double loFactor = POW(loSideValue,exponent);
        double hiFactor = POW(hiSideValue,exponent);

        rhs[(dim*jth) + idir] = hiMom*hiFactor - loMom*loFactor;

        // Add the Taylor series terms
        for (int order = 1; order <= a_orderPmax; order++)
        {
          vector<IvDim> taylorMonomials;

          generateMultiIndices(taylorMonomials,order);

          for (int i = 0; i < taylorMonomials.size(); i++)
          {
            const IvDim & taylorMonomial = taylorMonomials[i];

            IvDim totalMonomial = mono + taylorMonomial;

            if (m_cutCellMoments.m_EBmoments.find(totalMonomial) !=
                m_cutCellMoments.m_EBmoments.end())
            {
              double normalDerivative = m_cutCellMoments.m_IFData.m_normalDerivatives[taylorMonomial][idir];
              double fact = factorial(taylorMonomial);

              double moment = m_cutCellMoments.m_EBmoments[totalMonomial];

              rhs[(dim*jth) + idir] += normalDerivative * moment / fact;
            }
            else
            {
              PR_error("Unable to find needed monomial for Taylor series");
            }
          }
        }
      }
    }

    return rhs;
  }


// Operators
  template <int dim> void MinimalCCCM<dim>::operator=(const MinimalCCCM<dim> & a_MinimalCCCM)
  {
    // Only copy if the objects are distinct
    if (this != &a_MinimalCCCM)
    {
      m_cutCellMoments = a_MinimalCCCM.m_cutCellMoments;
    }
  }


  template <int dim> double MinimalCCCM<dim>::factorial(const IvDim & a_multiIndex) const
  {
    double fact = 1.0;

    for (int i = 0; i < dim; i++)
    {
      for (int j = 2; j <= a_multiIndex[i]; j++)
      {
        fact *= j;
      }
    }

    return fact;
  }

// Null constructor
  MinimalCCCM<1>::MinimalCCCM()
  {
  }

// Copy constructor
  MinimalCCCM<1>::MinimalCCCM(const MinimalCCCM<1>& a_thisMinimalCCCM)

    :m_cutCellMoments(a_thisMinimalCCCM.m_cutCellMoments)
  {
  }

// This constructor is used in the recursion
  MinimalCCCM<1>::MinimalCCCM(const IFData<1> & a_info)
    :m_cutCellMoments(a_info)
  {
  }

// Destructor
  MinimalCCCM<1>::~MinimalCCCM()
  {
  }

// Integrate along line segments aligned in a coordinate direction
  void MinimalCCCM<1>::computeMoments(const int              & a_orderPmax,
                                      const int              & a_degreePmax)
  {
    int lo = 0;
    int hi = 1;
    int loSign = m_cutCellMoments.m_IFData.m_cornerSigns[lo];
    int hiSign = m_cutCellMoments.m_IFData.m_cornerSigns[hi];

    // If entire edge out of the fluid, then moments = 0.0
    if (loSign <= ON && hiSign<= ON)
    {
      for (int iDegree = 0; iDegree <= a_degreePmax; ++iDegree)
      {
        // Definition of m_cutCellMoments.m_moments typedef requires that we define degree as
        // a oneTuple:
        IndexTM<int,1> degree;
        degree[0] = iDegree;
        m_cutCellMoments.m_moments[degree] = 0.0;
      }
    }
    else
    {
      // Assign loPt and hiPt in m_cutCellMoments.m_IFData.m_parentCoord system

      double loPt = LARGEREALVAL;
      double hiPt = LARGEREALVAL;

      // m_origin is an IndexTM<double,1>, which implies we need to use [0] everywhere
      // m_intersection is undefined if hiSign >= ON  && loSign >= ON
      if (loSign >= ON)
      {
        loPt = m_cutCellMoments.m_IFData.m_parentCoord.convertDir(
          -0.5*m_cutCellMoments.m_IFData.m_cellCenterCoord.m_dx[0],
          m_cutCellMoments.m_IFData.m_cellCenterCoord,
          0);
      }
      else
      {
        PR_assert(-0.5*m_cutCellMoments.m_IFData.m_cellCenterCoord.m_dx[0]
                  <= m_cutCellMoments.m_IFData.m_intersection
                  &&
                  m_cutCellMoments.m_IFData.m_intersection
                  <= 0.5*m_cutCellMoments.m_IFData.m_cellCenterCoord.m_dx[0]);

        loPt = m_cutCellMoments.m_IFData.m_parentCoord.convertDir(
          m_cutCellMoments.m_IFData.m_intersection,
          m_cutCellMoments.m_IFData.m_cellCenterCoord,
          0);
      }

      if (hiSign >= ON)
      {
        hiPt = m_cutCellMoments.m_IFData.m_parentCoord.convertDir(
          0.5*m_cutCellMoments.m_IFData.m_cellCenterCoord.m_dx[0],
          m_cutCellMoments.m_IFData.m_cellCenterCoord,
          0);
      }
      else
      {
        PR_assert(-0.5*m_cutCellMoments.m_IFData.m_cellCenterCoord.m_dx[0]
                  <= m_cutCellMoments.m_IFData.m_intersection
                  &&
                  m_cutCellMoments.m_IFData.m_intersection
                  <= 0.5*m_cutCellMoments.m_IFData.m_cellCenterCoord.m_dx[0]);

        hiPt = m_cutCellMoments.m_IFData.m_parentCoord.convertDir(
          m_cutCellMoments.m_IFData.m_intersection,
          m_cutCellMoments.m_IFData.m_cellCenterCoord,
          0);
      }

      // Integrate x^degree over the line segment[loPt,hiPt]
      computeMomentsUsingBinomial(loPt,hiPt,loSign,hiSign,a_degreePmax);
    }
  }

  void MinimalCCCM<1>::simpleComputeMoments(const double & a_loPt,
                                            const double & a_hiPt,
                                            const int  & a_degreePmax)
  {
    for (int iDegree = 0; iDegree <= a_degreePmax; ++iDegree)
    {
      //definition of m_cutCellMoments.m_moments typedef requires that we define degree thus:
      IndexTM<int,1>degree;
      degree[0] = iDegree;
      m_cutCellMoments.m_moments[degree] = POW(a_hiPt,iDegree + 1) - POW(a_loPt,iDegree +1);
      //    double dxFactor = POW(m_cutCellMoments.m_IFData.m_globalCoord.m_dx[0],iDegree + 1);
      //m_cutCellMoments.m_moments[degree] *= dxFactor;
      m_cutCellMoments.m_moments[degree] /= (iDegree + 1);
    }
  }

  void MinimalCCCM<1>::computeMomentsUsingBinomial(const double & a_loPt,
                                                   const double & a_hiPt,
                                                   const int  & a_loSign,
                                                   const int  & a_hiSign,
                                                   const int  & a_degreePmax)
  {
    for (int iDegree = 0; iDegree <= a_degreePmax; ++iDegree)
    {
      double epsilon = a_hiPt - a_loPt;

      // Definition of m_cutCellMoments.m_moments typedef requires that we define degree thus:
      IndexTM<int,1> degree;
      degree[0] = iDegree;
      m_cutCellMoments.m_moments[degree] = 0.0;

      // New method without substracting higher order terms
      for (int j = 1; j <= iDegree + 1; j++)
      {
        int bigger = j;
        int smaller = j;
        if (iDegree + 1 - j > j)
        {
          bigger = iDegree + 1 - j;
        }
        else
        {
          smaller = iDegree + 1 - j;
        }

        int numerator = 1;
        for (int i = bigger + 1; i <= iDegree + 1; ++i)
        {
          numerator *= i;
        }

        int denominator = 1;
        for (int i = 1; i <= smaller; ++i)
        {
          denominator *= i;
        }

        double factor = numerator / denominator;
        if (a_loSign >= ON)
        {
          m_cutCellMoments.m_moments[degree] += factor * POW(a_loPt,iDegree + 1 - j) * POW(epsilon,j);
        }
        else if (a_hiSign >= ON)
        {
          m_cutCellMoments.m_moments[degree] -= factor * POW(a_hiPt,iDegree + 1 - j) * POW(epsilon,j) * POW(-1.0,j);
        }
      }

      //double dxFactor = POW(m_cutCellMoments.m_IFData.m_globalCoord.m_dx[0],iDegree + 1);
      //m_cutCellMoments.m_moments[degree] *= dxFactor;
      m_cutCellMoments.m_moments[degree] /= iDegree + 1;
    }
  }


  void MinimalCCCM<1>::operator=(const MinimalCCCM<1> & a_MinimalCCCM)
  {
    // Only copy of the two objects are distinct
    if (this != &a_MinimalCCCM)
    {
      m_cutCellMoments = a_MinimalCCCM.m_cutCellMoments;
    }
  }

}

#endif
