#pragma once
#ifndef _PROTO_AMR_OP_
#define _PROTO_AMR_OP_

#include "Proto_AMRData.H"

namespace Proto {

/// AMR-Scope Operator
/**
    AMROp contains the necessary tools to apply a user-defined descendent of 
    BoxOp on an AMR hierarchy.

    \tparam OpType  A BoxOp class
    \tparam T       Datatype of the data holder (e.g. int, double, etc.)
    \tparam MEM     Proto::MemType of the data holder
*/

template <template<typename, MemType> class OPType,
    typename T, MemType MEM>
class AMROp
{
    public:

    typedef OPType<T, MEM> BOP;
    typedef LevelOp<OPType, T, MEM> LOP;
    typedef AMRData<T, BOP::numState(), MEM, PR_CELL> AMRStateData;
    typedef AMRData<T, BOP::numAux(),   MEM, PR_CELL> AMRAuxData;
    typedef LevelBoxData<T, BOP::numState(), MEM, PR_CELL> LevelStateData;
    typedef LevelBoxData<T, BOP::numAux(),   MEM, PR_CELL> LevelAuxData;
    typedef BoxData<T, BOP::numState(), MEM> StateData;
    typedef BoxData<T, BOP::numAux(), MEM>   AuxData;

    static constexpr unsigned int numState() { return BOP::numState(); }
    static constexpr unsigned int numAux()   { return BOP::numAux(); }
    
    /// Ghost Region Size (State)
    /**
        Amount of data needed in ghost regions of state variables to apply the operator
    */ 
    static Point ghost()      { return BOP::ghost(); }
    
    /// Ghost Region Size (Auxiliary)
    /**
        Amount of data needed in ghost regions of auxiliary variables to apply the operator
    */ 
    static Point auxGhost()
    {
        if (numAux() > 0) { return BOP::auxGhost(); }
        return Point::Zeros();
    }

    static int order()
    {
        return BOP::order();
    }
   
    /// Default Constructor 
    inline AMROp();

    /// Non-Trivial Constructor
    /**
        \param grid  AMR Grid
        \param cdx   Coarsest level grid spacing
    */
    inline AMROp(AMRGrid& a_grid, T a_cdx);
    inline AMROp(AMRGrid& a_grid, std::array<T, DIM>& a_cdx);
    
    /// Lazy Constructor
    /**
        \param grid  AMR Grid
        \param cdx   Coarsest level grid spacing
    */
    inline void define(AMRGrid& a_grid, T a_cdx);
    inline void define(AMRGrid& a_grid, std::array<T, DIM>& a_cdx);
    
    /// Apply
    /**
        Compute L(phi, rho) on all levels of the hierarchy including reflux
        corrections and boundary condition interpolation where necessary

        \param output   Evaluated operator (output)
        \param state    State variables
        \param aux      Auxiliary variables
    */
    inline void
    operator()(
        AMRStateData&       a_output,
        const AMRStateData& a_state,
        const AMRAuxData&   a_aux,
        T                   a_scale = 1.0);
    
    /// Apply
    /**
        Compute L(phi) on all levels of the hierarchy including reflux
        corrections and boundary condition interpolation where necessary

        \param output   Evaluated operator (output)
        \param state    State variables
    */
    inline void
    operator()(
        AMRStateData& a_output,
        AMRStateData& a_state,
        T             a_scale = 1.0);

    /// Level Apply
    /**
        Compute L(phi, rho) on a specified level of the hierarchy including reflux
        corrections and boundary condition interpolation.

        \param output   Evaluated operator (output)
        \param state    State variables
        \param aux      Auxiliary variables
        \param level    Specified level of application
    */
    inline void
    levelApply(
        LevelStateData& a_output,
        AMRStateData&   a_state,
        AMRAuxData&     a_aux,
        int             a_level,
        T               a_scale = 1.0); 
    
    /// Level Apply
    /**
        Compute L(phi) on a specified level of the hierarchy including reflux
        corrections and boundary condition interpolation.

        \param output   Evaluated operator (output)
        \param state    State variables
        \param level    Specified level of application
    */
    inline void
    levelApply(
        LevelStateData& a_output,
        AMRStateData&   a_state,
        int             a_level,
        T               a_scale = 1.0); 

    /// Reflux
    /**
        Given state variable information on two levels l and l+1, compute and apply
        the appropriate reflux correction to a dataset on level l.

        \param crseOut      The coarse level dataset to which refluxing will be applied
        \param crseState    coarse level state data (level l)
        \param fineState    fine level state data (level l+1)
        \param crseLevel    The level on which refluxing is computed (level l)
        \param scale        An optional value which will scale the reflux correction
    */
    inline void
    reflux(
        LevelStateData& a_crseOut,
        LevelStateData& a_crseState,
        LevelStateData& a_fineState,
        int             a_crseLevel,
        T               a_scale) const;


    inline InterpStencil<T>& interpStencil(int a_level);

    inline void setDiagScale(T a_value, int a_gridFactor = 0);
    inline void setFluxScale(T a_value, int a_gridFactor = 0);
    inline std::array<T, DIM> dx(int a_level) const { return m_dx[a_level]; }
    const LOP& levelOp(int a_level) const { return m_levelOps[a_level]; }
    const BOP& boxOp(int a_level) const { return levelOp(a_level).op(); }

    private:

    bool m_defined;    
    AMRGrid m_grid; //this is intentionally a deep copy.
    std::vector<std::array<T, DIM>>      m_dx;
    std::vector<LOP>         m_levelOps;
    std::vector<std::shared_ptr<LevelFluxRegister<T, BOP::numState(), MEM>>> m_fluxRegisters;
    std::vector<InterpStencil<T>> m_interp;
};

#include "implem/Proto_AMROpImplem.H"

} // end namespace Proto
#endif // end include guard

