#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

// BVS, June 18, 2003

#ifndef _C4_BICGSTABSOLVER_H_
#define _C4_BICGSTABSOLVER_H_

#include "Chombo_parstream.H"
#include "Chombo_NamespaceHeader.H"

///
/**
   Elliptic solver using the BiCGStab algorithm.
*/
template <class dataholder_t, class op>
class BiCGStabSolver
{
public:
  BiCGStabSolver();
  ~BiCGStabSolver();
  
  inline
  static int
  solve(dataholder_t & a_phi, const dataholder_t& a_rhs,
        const op& a_op,
        int a_verbosity = 0,
        Real a_metric = -1.0,
        Real a_tol = 1.0e-6,
        Real a_hang = 1.0e-8,
        Real a_small = 1.0e-16,
        int  a_imax = 27,
        int numRestarts = 5)
  {
    if(a_imax == 0) return 1;
    CH_TIMERS("BiCGStabSolver::solve");
    CH_TIMER("BiCGStabSolver::solve::Initialize",timeInitialize);
    CH_TIMER("BiCGStabSolver::solve::MainLoop",timeMainLoop);
    CH_TIMER("BiCGStabSolver::solve::Cleanup",timeCleanup);

    CH_START(timeInitialize);
    int exitStatus = 1;
    dataholder_t r, r_tilde, e, p, p_tilde, s_tilde, t, v;
    a_op.create(r_tilde);
    a_op.create(p_tilde);
    a_op.create(s_tilde);
    a_op.create(r);
    a_op.create(e);
    a_op.create(p);
    a_op.create(t);
    a_op.create(v);

    int recount = 0;

    a_op.setToZero(r); // added by petermc, 26 Nov 2013, to zero out ghosts

    a_op.residual(r, a_phi, a_rhs);

    a_op.assignLocal(r_tilde, r);
    a_op.setToZero(e);
    // (DFM 2/1/07) these next two need to be set to zero to prevent
    // problems in the multilevel case
    a_op.setToZero(p_tilde);
    a_op.setToZero(s_tilde);

    int i = 0;

    // rho[0] = r_i , rho[1] = r_(i-1), etc.
    Real rho[4]=
      {
        0,0,0,0
      };
    Real norm[2];
    norm[0] = a_op.norm(r);
    Real initial_norm = norm[0];
    Real initial_rnorm = norm[0];
    norm[1] = norm[0];

    Real alpha[2] =
      {
        0,0
      };
    Real  beta[2] =
      {
        0,0
      };
    Real omega[2] =
      {
        0,0
      };

    bool init = true;
    int restarts = 0;

    if (a_verbosity >= 5)
    {
      pout() << "      BiCGStab:: initial Residual norm = "
             << initial_norm << "\n";
    }

    // if a convergence metric has been supplied, replace initial residual
    // with the supplied convergence metric...
    if (a_metric > 0)
    {
      initial_norm = a_metric;
    }

    CH_STOP(timeInitialize);

    CH_START(timeMainLoop);
    while ((i<a_imax && norm[0] > a_tol*norm[1]) && (norm[1] > 0))
    {
      i++;

      norm[1] = norm[0];
      alpha[1]= alpha[0];
      beta[1] = beta[0];
      omega[1]= omega[0];

      if (a_verbosity >= 5)
      {
        pout() << "      BiCGStab::       norm[0]  = "  << norm[0] << ", "
               <<                        "norm[1]  = "  << norm[1]
               << "\n";
        pout() << "      BiCGStab::       alpha[0] = "  << alpha[0] << ", "
               <<                        "alpha[1] = "  << alpha[1]
               << "\n";
        pout() << "      BiCGStab::       beta[0]  = "  << beta[0] << ", "
               <<                        "beta[1]  = "  << beta[1]
               << "\n";
        pout() << "      BiCGStab::       omega[0] = "  << omega[0] << ", "
               <<                        "omega[1] = "  << omega[1]
               << "\n";
      }

      rho[3] = rho[2];
      rho[2] = rho[1];
      rho[1] = a_op.dotProduct(r_tilde, r);

      if (a_verbosity >= 5)
      {
        pout() << "      BiCGStab::       rho[1] = "  << rho[1] << ", "
               <<                        "rho[2] = "  << rho[2] << ", "
               <<                        "rho[3] = "  << rho[3]
               << "\n";
      }

      if (rho[1] == 0.0)
      {
        // we are finished, we will not converge anymore
        a_op.incr(a_phi, e, 1.0);

        if (a_verbosity >= 5)
        {
          pout() << "      BiCGStab:: rho = 0, returning"
                 << " -- Residual norm = "
                 << norm[0] << "\n";
        }

        CH_STOP(timeMainLoop);

        CH_START(timeCleanup);

        exitStatus = 2;

        CH_STOP(timeCleanup);

        return exitStatus;
      }

      if (init)
      {
        a_op.assignLocal(p, r);
        init = false;
      }
      else
      {
        beta[1] = (rho[1]/rho[2])*(alpha[1]/omega[1]);
        a_op.scale(p, beta[1]);
        a_op.incr(p, v, -beta[1]*omega[1]);
        a_op.incr(p, r, 1.0);
      }

      if (a_verbosity >= 5)
      {
        pout() << "      BiCGStab::       beta[1]  = "  << beta[1]
               << "\n";
      }

      a_op.preCond(p_tilde, p);
      a_op.setToZero(v); // added by petermc, 27 Nov 2013, to zero out ghosts
      a_op.applyOp(v, p_tilde);
      Real m = a_op.dotProduct(r_tilde, v);
      alpha[0] = rho[1]/m;

      if (a_verbosity >= 5)
      {
        pout() << "      BiCGStab::       rho[1] = "  << rho[1]   << ", "
               <<                             "m = "  << m        << ", "
               <<                      "alpha[0] = "  << alpha[0]
               << "\n";
      }

      if (std::abs(m) > a_small*std::abs( rho[1]))
      {
        a_op.incr(r, v, -alpha[0]);
        norm[0] = a_op.norm(r);
        a_op.incr(e, p_tilde, alpha[0]);
      }
      else
      {
        a_op.setToZero(r);
        norm[0] = 0.0;
      }

      if (a_verbosity >= 4)
      {
        pout() << "norm[0] = " << norm[0]
               << " initial_norm = " << initial_norm
               << " initial_rnorm = " << initial_rnorm << endl;
      }

      if (norm[0] > a_tol*initial_norm && norm[0])
      {
        a_op.preCond(s_tilde, r);
        a_op.setToZero(t); // added by petermc, 27 Nov 2013, to zero out ghosts
        a_op.applyOp(t, s_tilde);
        omega[0] = a_op.dotProduct(t, r)/a_op.dotProduct(t, t);
        a_op.incr(e, s_tilde, omega[0]);
        a_op.incr(r, t,      -omega[0]);
        norm[0] = a_op.norm(r);
      }

      if (a_verbosity >= 4)
      {
        pout() << "      BiCGStab::     iteration = "  << i        << ", error norm = " << norm[0] << ", rate = " << norm[1]/norm[0] << "\n";
      }

      if (norm[0] <= a_tol*initial_norm)
      {
        // converged to tolerance
        exitStatus = 1;
        break;
      }

      if (omega[0] == 0.0 || norm[0] > (1-a_hang)*norm[1])
      {
        if (recount == 0)
        {
          recount = 1;
        }
        else
        {
          recount = 0;
          a_op.incr(a_phi, e, 1.0);
          if (restarts == numRestarts)
          {
            if (a_verbosity >= 4)
            {
              pout() << "      BiCGStab: max restarts reached" << endl;
              pout() << "                init  norm = " << initial_norm << endl;
              pout() << "                final norm = " << norm[0] << endl;
            }

            CH_STOP(timeMainLoop);

            exitStatus = 4;
            return exitStatus;
          }

          {
            CH_TIME("BiCGStabSolver::solve::Restart");

            a_op.residual(r, a_phi, a_rhs);
            norm[0] = a_op.norm(r);
            rho[1]=0.0; rho[1]=0.0; rho[2]=0.0; rho[3]=0.0;
            alpha[0]=0; beta[0]=0; omega[0]=0;
            a_op.assignLocal(r_tilde, r);
            a_op.setToZero(e);

            restarts++;
          }

          if (a_verbosity >= 4)
          {
            pout() << "      BiCGStab::   restart =  " << restarts << "\n";
          }

          init = true;
        }
      }
    }
    CH_STOP(timeMainLoop);

    CH_START(timeCleanup);

    if (a_verbosity >= 4)
    {
      pout() << "      BiCGStab:: " << i << " iterations, final Residual norm = "
             << norm[0] << "\n";
    }

    a_op.incr(a_phi, e, 1.0);


    CH_STOP(timeCleanup);
    return exitStatus;
  }
};
#include "Chombo_NamespaceFooter.H"
#endif /*_BICGSTABSOLVER_H_*/
