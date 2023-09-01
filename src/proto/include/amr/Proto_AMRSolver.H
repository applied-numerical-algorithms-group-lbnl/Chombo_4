#pragma once
#ifndef _PROTO_AMR_SOLVER_
#define _PROTO_AMR_SOLVER_

#include "Proto_AMROp.H"

namespace Proto {

/// AMR Solver
/**
    AMR Solver is an interface for building a solver to invert an AMROp. 
    In particular, we would like to solve the equation L(phi) = G where L is defined by an
    AMROp, G is an external forcing defined on an AMR hierarchy, and phi is the solution,
    also on a hierarchy.

    TODO: specifically, this is an interface for *iterative* solvers.
*/


template <template<typename, MemType> class OpType, typename T, MemType MEM = MEMTYPE_DEFAULT>
class AMRSolver
{
    public:
    typedef OpType<T, MEM> BOP;
    typedef LevelOp<OpType, T, MEM> LOP;
    typedef LevelOp<OpType, T, MEM> AMROP;
    typedef AMRData<T, BOP::numState(), MEM> AMRStateData;
    typedef AMRData<T, BOP::numAux(),   MEM> AMRAuxData;
    static constexpr unsigned int numState() { return BOP::numState(); }
    static constexpr unsigned int numAux()   { return BOP::numAux(); }
    
    inline virtual double
    solve(
        AMRStateData& a_state,
        AMRStateData& a_forcing,
        AMRAuxData&   a_aux,
        int           a_maxIter,
        double        a_tolerance) = 0;
    
    inline virtual double
    solve(
        AMRStateData& a_state,
        AMRStateData& a_forcing,
        int           a_maxIter,
        double        a_tolerance) = 0;

};
} // end namespace Proto

#endif // end include guard
