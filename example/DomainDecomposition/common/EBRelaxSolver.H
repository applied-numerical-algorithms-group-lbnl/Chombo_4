#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

// derived from:  BVS, Feb 7, 2005

#ifndef _EB_RELAXSOLVER_H_
#define _EB_RELAXSOLVER_H_

#include "EBMultigrid.H"
#include "Chombo_NamespaceHeader.H"

///
class EBRelaxSolver
{

public:
  ///
  /**
     Define the solver.   a_op is the linear operator.
  */
  EBRelaxSolver(EBMultigridLevel*                     a_op, 
                const DisjointBoxLayout             & a_grids,
                const shared_ptr<LevelData<EBGraph> > a_graphs,
                const IntVect                       & a_ghost)
  {
    m_op = a_op;
    m_res.define(a_grids, a_ghost, a_graphs);
    m_cor.define(a_grids, a_ghost, a_graphs);
  }


  ///
  /**
     Solve L(phi) = rho
  */
  void solve(EBLevelBoxData<CELL, 1>      & a_phi, 
             const EBLevelBoxData<CELL, 1>& a_rhs, 
             int   a_maxIter, 
             Real  a_tolerance)
  {
    CH_TIME("EBRelaxSolver::solve");
    m_res.setVal(0);
    m_op->residual(m_res, a_phi, a_rhs);

    Real norm = m_res.maxNorm(0);
    m_minNorm = norm;

    pout() << "      EBRelaxSolver:: initial Residual norm = "  << norm << "\n";


    if (norm > 0.) // Do not iterate if norm == 0.
    {
      Real initialNorm = norm;
      int iter =0;
      while (iter < a_maxIter)
      {
        m_cor.setVal(0.);
        m_op->preCond(m_cor, m_res);
        a_phi += m_cor;

        m_op->residual(m_res, a_phi, a_rhs);

        Real oldNorm = norm;
        norm = m_res.maxNorm(0);
        iter++;

        pout() << "      EBRelaxSolver:: iteration = " << iter
               << " residual norm = " << norm
               << ", rate = " << oldNorm/norm << endl;

        if (norm <  a_tolerance*initialNorm) break;

        if (norm < m_minNorm)
        {
          m_minNorm = norm;
        }
      }
    }

  }
private:
  EBMultigridLevel* m_op;
  EBLevelBoxData<CELL, 1>  m_res, m_cor;
  Real m_minNorm;

  EBRelaxSolver();
  EBRelaxSolver(const EBRelaxSolver& a_in);
  void operator=(const EBRelaxSolver& a_in);


};

#include "Chombo_NamespaceFooter.H"
#endif /*_RELAXSOLVER_H_*/
