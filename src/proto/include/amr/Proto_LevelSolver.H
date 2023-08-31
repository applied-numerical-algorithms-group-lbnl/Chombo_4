#pragma once
#ifndef _PROTO_LEVEL_SOLVER_
#define _PROTO_LEVEL_SOLVER_

#include "Proto_LevelOp.H"

namespace Proto {

/// Level Solver
/**
    Level Solver is an interface for building a solver to invert a LevelOp. 
    In particular, we would like to solve the equation L(phi) = G where L is defined by a
    LevelOp, G is an external forcing defined on a level, and phi is the solution, also on a level.

    TODO: specifically, this is an interface for *iterative* solvers.
*/

template <template<typename, MemType> class OpType, typename T, MemType MEM = MEMTYPE_DEFAULT>
class LevelSolver
{
    public:
    typedef OpType<T, MEM> BOP;
    typedef LevelOp<OpType, T, MEM> LOP;
    typedef LevelBoxData<T, BOP::numState(), MEM, PR_CELL> LevelStateData;
    typedef LevelBoxData<T, BOP::numAux(),   MEM, PR_CELL> LevelAuxData;
    static constexpr unsigned int numState() { return BOP::numState(); }
    static constexpr unsigned int numAux()   { return BOP::numAux(); }
    
    inline virtual double
    solve(
        LevelStateData& a_state,
        LevelStateData& a_forcing,
        LevelAuxData&   a_aux,
        int             a_maxIter,
        double          a_tolerance) = 0;
    
    inline virtual double
    solve(
        LevelStateData& a_state,
        LevelStateData& a_forcing,
        int             a_maxIter,
        double          a_tolerance) = 0;

};
} // end namespace Proto

#endif // end include guard
