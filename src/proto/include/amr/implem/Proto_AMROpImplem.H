
template <template<typename, MemType> class OPType, typename T, MemType MEM>
AMROp<OPType, T, MEM>::AMROp()
{
    m_defined = false;
}

template <template<typename, MemType> class OPType, typename T, MemType MEM>
AMROp<OPType, T, MEM>::AMROp(AMRGrid& a_grid, T a_cdx)
{
    define(a_grid, a_cdx);
}

template <template<typename, MemType> class OPType, typename T, MemType MEM>
AMROp<OPType, T, MEM>::AMROp(AMRGrid& a_grid, std::array<T, DIM>& a_cdx)
{
    define(a_grid, a_cdx);
}

template <template<typename, MemType> class OPType, typename T, MemType MEM>
void AMROp<OPType, T, MEM>::define(AMRGrid& a_grid, T a_cdx)
{
    std::array<T, DIM> cdx;
    cdx.fill(a_cdx);
    define(a_grid, cdx);
}

template <template<typename, MemType> class OPType, typename T, MemType MEM>
void AMROp<OPType, T, MEM>::define(AMRGrid& a_grid, std::array<T, DIM>& a_cdx)
{
    int numLevels = a_grid.numLevels();
    m_dx.clear();
    m_levelOps.clear();
    m_fluxRegisters.clear();

    m_dx.resize(numLevels);
    m_levelOps.resize(numLevels);
    m_fluxRegisters.resize(numLevels - 1);
    m_interp.resize(numLevels - 1); 
    m_grid = a_grid;
    auto dx = a_cdx;
    for (int lvl = 0; lvl < numLevels; lvl++)
    {
        m_dx[lvl] = dx;
        m_levelOps[lvl].define(dx);
        if (lvl < numLevels - 1)
        {
            Point refRatio = m_grid.refRatio(lvl);
            m_fluxRegisters[lvl] =
                std::make_shared<LevelFluxRegister<T, BOP::numState(), MEM>>(
                    a_grid[lvl], a_grid[lvl+1],
                    refRatio, dx);
            m_interp[lvl] = 
                InterpStencil<T>::Build(BOP::order() + 1, refRatio);
            for (int dir = 0; dir < DIM; dir++)
            {
                dx[dir] /= refRatio[dir];
            }
        }
    }
    m_defined = true;
}

template <template<typename, MemType> class OPType, typename T, MemType MEM>
void AMROp<OPType, T, MEM>::operator()(
        AMRStateData& a_output,
        const AMRStateData& a_state,
        const AMRAuxData&   a_aux,
        T                   a_scale)
{
    PROTO_ASSERT(a_output.layout().compatible(a_state.layout()),
        "AMROp::operator() | Error: output and state have incompatible layouts.");
    PROTO_ASSERT(a_output.layout().compatible(a_aux.layout()),
        "AMROp::operator() | Error: output and aux have incompatible layouts.");
    
    a_state.averageDown();
    for (int lvl = 0; lvl < m_grid.numLevels(); lvl++)
    {
        auto& output = a_output[lvl];
        auto& state  = a_state[lvl];
        auto& aux    = a_aux[lvl];
        auto& op     = m_levelOps[lvl];
        
        // interpolate boundary conditions
        if (lvl > 0)
        {
            auto& crseState = a_state[lvl-1];
            interpBoundaries(crseState, state, m_interp);
        }
        // apply the operator
        op(output, state, aux, a_scale);
        // reflux coarse level
        // FIXME: This method computes all of the fluxes three times and should
        // definitely be replaced with something more clever.
        if (lvl > 0)
        {
            auto& crseState = a_state[lvl-1];
            auto& crseOutput = a_output[lvl-1];
            reflux(crseOutput, crseState, state, lvl-1, a_scale);
        }
    }
}

template <template<typename, MemType> class OPType, typename T, MemType MEM>
void AMROp<OPType, T, MEM>::operator()(
        AMRStateData& a_output,
        AMRStateData& a_state,
        T             a_scale)
{
    PROTO_ASSERT(a_output.grid().compatible(a_state.grid()),
        "AMROp::operator() | Error: output and state have incompatible grids.");
    
    // average down all levels recursively
    a_state.averageDown(); 
  
    for (int lvl = 0; lvl < m_grid.numLevels(); lvl++)
    {
        auto& output = a_output[lvl];
        auto& state  = a_state[lvl];
        auto& op     = m_levelOps[lvl];
        
        // if a coarser level exists, interpolate boundaries
        if (lvl > 0)
        {
            auto& crseState = a_state[lvl-1];
            interpBoundaries(crseState, state, m_interp[lvl-1]);
        }
        
        // apply the level operator
        op(output, state, a_scale);
        
        // reflux coarse level
        // FIXME: This method computes all of the fluxes three times and should
        // definitely be replaced with something more clever.
        if (lvl > 0)
        {
            auto& crseState = a_state[lvl-1];
            auto& crseOutput = a_output[lvl-1];
            reflux(crseOutput, crseState, state, lvl-1, a_scale);
        }
    }
}

template <template<typename, MemType> class OPType, typename T, MemType MEM>
void AMROp<OPType, T, MEM>::levelApply(
        LevelStateData& a_output,
        AMRStateData& a_state, 
        AMRAuxData& a_aux, 
        int a_level,
        T   a_scale)
{
    PROTO_ASSERT(a_output.layout().compatible(a_state[a_level].layout()),
        "AMROp::operator() | Error: output and state[%i] have incompatible layouts.", a_level);
    PROTO_ASSERT(a_output.layout().compatible(a_aux[a_level].layout()),
        "AMROp::operator() | Error: output and aux[%i] have incompatible layouts.", a_level);

    auto& state = a_state[a_level];
    auto& aux   = a_aux[a_level];
    auto& op = m_levelOps[a_level];
    
    // if there is a finer level, average down to this level
    if (a_level < m_grid.numLevels() - 1)
    {
        // aux data shouldn't need averaging
        Point refRatio = m_grid.refRatio(a_level);
        auto& fineState = a_state[a_level + 1];
        averageDown(state, fineState, refRatio);
    }
    
    // if there is a coarser level, average down to it and then
    // interpolate boundary conditions from it
    if (a_level > 0)
    {
        // aux data shouldn't need averaging / interpolation
        Point refRatio = m_grid.refRatio(a_level-1);
        auto& crseState = a_state[a_level-1];
        averageDown(crseState, state, refRatio);
        interpBoundaries(crseState, state, m_interp[a_level-1]);
    }
     
    op(a_output, state, aux, a_scale);
    
    if (a_level < m_grid.numLevels()-1)
    {
        auto& fineState = a_state[a_level+1];
        reflux(a_output, state, fineState, a_level, a_scale);
    }
}


template <template<typename, MemType> class OPType, typename T, MemType MEM>
void AMROp<OPType, T, MEM>::levelApply(
        LevelStateData& a_output,
        AMRStateData& a_state, 
        int a_level,
        T   a_scale)
{
    PROTO_ASSERT(a_output.layout().compatible(a_state[a_level].layout()),
        "AMROp::operator() | Error: output and state[%i] have incompatible layouts.", a_level);
     
    auto& state = a_state[a_level];
    auto& op = m_levelOps[a_level];
    
    // if there is a finer level, average down to this level
    if (a_level < m_grid.numLevels() - 1)
    {
        Point refRatio = m_grid.refRatio(a_level);
        auto& fineState = a_state[a_level + 1];
        averageDown(state, fineState, refRatio);
    }
    
    // if there is a coarser level, average down to it and then
    // interpolate boundary conditions from it
    if (a_level > 0)
    {
        Point refRatio = m_grid.refRatio(a_level-1);
        auto& crseState = a_state[a_level-1];
        averageDown(crseState, state, refRatio);
        interpBoundaries(crseState, state, m_interp[a_level-1]);
    }
    
    // apply the level operator
    op(a_output, state, a_scale);

    // if there is a finer level, use it to compute reflux correction
    if (a_level < a_state.numLevels() - 1)
    {
        auto& fineState = a_state[a_level + 1];
        reflux(a_output, state, fineState, a_level, a_scale);
    }
}

template <template<typename, MemType> class OPType, typename T, MemType MEM>
void AMROp<OPType, T, MEM>::reflux(
        LevelStateData& a_crseOut,
        LevelStateData& a_crseState,
        LevelStateData& a_fineState,
        int a_crseLevel,
        T   a_scale) const
{
    auto& fluxRegister = *(m_fluxRegisters[a_crseLevel]);
    auto& crseOp = m_levelOps[a_crseLevel];
    auto& fineOp = m_levelOps[a_crseLevel+1];
    
    PROTO_ASSERT(a_crseOut.layout().compatible(fluxRegister.crseLayout()),
            "AMROp::reflux | Error coarse output has an incompatible layout.");
    PROTO_ASSERT(a_crseState.layout().compatible(fluxRegister.crseLayout()),
            "AMROp::reflux | Error coarse input has an incompatible layout.");
    PROTO_ASSERT(a_fineState.layout().compatible(fluxRegister.fineLayout()),
            "AMROp::reflux | Error fine input has an incompatible layout.");
    
    fluxRegister.reset();
   
    for (auto iter = fluxRegister.crseLayout().begin(); iter.ok(); ++iter)
    {
        for (int dir = 0; dir < DIM; dir++)
        {
            auto& crse_i = a_crseState[*iter];
            StateData flux_id(fluxRegister.crseLayout()[*iter].grow(dir, Side::Hi, 1));
            crseOp.op().flux(flux_id, crse_i, dir);
            fluxRegister.incrementCoarse(flux_id, *iter, 1.0, dir);
        }
    }
    for (auto iter = fluxRegister.fineLayout().begin(); iter.ok(); ++iter)
    {
        for (int dir = 0; dir < DIM; dir++)
        {
            auto& fine_i = a_fineState[*iter];
            StateData flux_id(fluxRegister.fineLayout()[*iter].grow(dir, Side::Hi, 1));
            fineOp.op().flux(flux_id, fine_i, dir);
            fluxRegister.incrementFine(flux_id, *iter, 1.0, dir);
        }
    }

    double scale = a_scale * crseOp.fluxScale();
    fluxRegister.reflux(a_crseOut, scale);
}   

template <template<typename, MemType> class OPType, typename T, MemType MEM>
InterpStencil<T>& AMROp<OPType, T, MEM>::interpStencil(int a_level)
{
    PROTO_ASSERT(a_level < m_grid.numLevels()-1,
        "AMROp::interpStencil | Error: level %i is out of bounds.", a_level);
    return m_interp[a_level];
}

template <template<typename, MemType> class OPType, typename T, MemType MEM>
void AMROp<OPType, T, MEM>::setDiagScale(T a_value, int a_gridFactor)
{
    T value = a_value;
    for (int lvl = 0; lvl < m_grid.numLevels(); lvl++)
    {
        m_levelOps[lvl].setDiagScale(value);
        if (lvl < m_grid.numLevels() - 1)
        {
            Point refRatio = m_grid.refRatio(lvl);
            if (a_gridFactor != 0)
            {
                PROTO_ASSERT(refRatio.isotropic(), 
                        "AMROp::setDiagScale | Error: \
                        non-zero gridFactor input is not implemented for anisotropic \
                        refinement ratios.");
            }
            double factor = pow(refRatio[0], a_gridFactor);
            value /= factor;
        }
    }
}

template <template<typename, MemType> class OPType, typename T, MemType MEM>
void AMROp<OPType, T, MEM>::setFluxScale(T a_value, int a_gridFactor)
{
    T value = a_value;
    for (int lvl = 0; lvl < m_grid.numLevels(); lvl++)
    {
        m_levelOps[lvl].setFluxScale(value);
        if (lvl < m_grid.numLevels() - 1)
        {
            Point refRatio = m_grid.refRatio(lvl);
            if (a_gridFactor != 0)
            {
                PROTO_ASSERT(refRatio.isotropic(), 
                        "AMROp::setFluxScale | Error: \
                        non-zero gridFactor input is not implemented for anisotropic \
                        refinement ratios.");
            }
            double factor = pow(refRatio[0], a_gridFactor);
            value /= factor;
        }
    }
}
