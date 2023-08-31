
template <template<typename, MemType> class OpType, typename T, MemType MEM>
LevelOp<OpType, T, MEM>::LevelOp() {}

template <template<typename, MemType> class OpType, typename T, MemType MEM>
LevelOp<OpType, T, MEM>::LevelOp(T a_dx)
{
    define(a_dx);
}

template <template<typename, MemType> class OpType, typename T, MemType MEM>
LevelOp<OpType, T, MEM>::LevelOp(std::array<T, DIM> a_dx)
{
    define(a_dx);
}

template <template<typename, MemType> class OpType, typename T, MemType MEM>
void LevelOp<OpType, T, MEM>::define(T a_dx)
{
    std::array<T, DIM> dx;
    dx.fill(a_dx);
    define(dx);
}

template <template<typename, MemType> class OpType, typename T, MemType MEM>
void LevelOp<OpType, T, MEM>::define(std::array<T, DIM> a_dx)
{
    m_op.define(a_dx);
}

template <template<typename, MemType> class OpType, typename T, MemType MEM>
void LevelOp<OpType, T, MEM>::operator()(
        LevelStateData& a_output,
        const LevelStateData& a_state,
        const LevelAuxData&   a_aux,
        T                     a_scale) const
{
    PROTO_ASSERT(a_output.layout().compatible(a_state.layout()),
            "LevelOp::operator() | Error: Output and state have incompatible layouts.");
    PROTO_ASSERT(a_output.layout().compatible(a_aux.layout()),
            "LevelOp::operator() | Error: Output and aux have incompatible layouts.");
    
    for (auto iter : a_state)
    {
        auto& out_i = a_output[iter];
        const auto& state_i = a_state[iter];
        const auto& aux_i   = a_aux[iter];
        m_op(out_i, state_i, aux_i, a_scale);
    }
}

template <template<typename, MemType> class OpType, typename T, MemType MEM>
void LevelOp<OpType, T, MEM>::operator()(
        LevelStateData&       a_output,
        const LevelStateData& a_state,
        T                     a_scale) const
{
    PROTO_ASSERT(a_output.layout().compatible(a_state.layout()),
            "LevelOp::operator() | Error: Output and state have incompatible layouts.");
    
    for (auto iter : a_state)
    {
        auto& out_i = a_output[iter];
        const auto& state_i = a_state[iter];
        std::array<StateData, DIM> fluxes;
        for (int dir = 0; dir < DIM; dir++)
        {
            fluxes[dir].define(out_i.box().grow(dir, Side::Hi, 1));
        }
        m_op(out_i, fluxes, state_i, a_scale);
    }
}
template <template<typename, MemType> class OpType, typename T, MemType MEM>
void LevelOp<OpType, T, MEM>::setDiagScale(T a_value)
{
    m_op.setDiagScale(a_value);
}

template <template<typename, MemType> class OpType, typename T, MemType MEM>
void LevelOp<OpType, T, MEM>::setFluxScale(T a_value)
{
    m_op.setFluxScale(a_value);
}

template <template<typename, MemType> class OpType, typename T, MemType MEM>
void LevelOp<OpType, T, MEM>::setTime(T a_time)
{
    m_op.setTime(a_time);
}

template <template<typename, MemType> class OpType, typename T, MemType MEM>
std::array<T, DIM> LevelOp<OpType, T, MEM>::dx() const
{
    return m_op.dx();
}
