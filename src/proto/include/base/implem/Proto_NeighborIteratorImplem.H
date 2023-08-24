NeighborIterator::NeighborIterator(const DisjointBoxLayout& a_dbl, const Box& a_bx)
{
    PROTO_ASSERT(!a_bx.empty(),
        "NeighborIterator::Constructor | Error: Input Box cannot be empty.");
    m_layout = a_dbl;
    m_bx = a_bx;
    m_current = DataIndex<BoxPartition>(a_dbl.m_partition,0);
    m_end = DataIndex<BoxPartition>();
    m_dProbC = m_layout.domain().coarsen(m_layout.boxSize());    

    // Set up BoxIterator.
    Box bxC = a_bx.coarsen(m_layout.boxSize());
    Box bxPoints = m_dProbC & bxC;
    m_bxit = BoxIterator(bxPoints);
}

NeighborIterator&
NeighborIterator::begin()
{
    if (!m_bxit.ok()) { m_current = m_end; }
    m_bxit.begin();
    Point ptCoarseShifted = m_dProbC.image(*m_bxit);
    while ((m_layout.find(ptCoarseShifted) == *m_layout.end()))
    {
        ++m_bxit;
        if (!m_bxit.ok()){ break; }
        ptCoarseShifted =m_dProbC.image(*m_bxit);
    }
    if (m_bxit.ok())
    {
        m_shift = (m_dProbC.image(*m_bxit) - *m_bxit)*m_layout.boxSize();
        m_current = m_layout.find(ptCoarseShifted);
    }
    else
    {
        m_current = m_end;
    }
    return *this;
}

NeighborIterator&
NeighborIterator::operator++()
{
    ++m_bxit;
    if (!m_bxit.ok())
    {
        m_current = m_end;
        return *this;
    }
    Point ptCoarseShifted = m_dProbC.image(*m_bxit);
    while ((m_layout.find(ptCoarseShifted) == *m_layout.end()))
    {
        ++m_bxit;
        if (!m_bxit.ok()){break;}
        ptCoarseShifted = m_dProbC.image(*m_bxit);
    }
    if (m_bxit.ok())
    {
        m_shift = (m_dProbC.image(*m_bxit) - *m_bxit)*m_layout.boxSize();
        m_current = m_layout.find(ptCoarseShifted);
    }
    else
    {
        m_current = m_end;
    }
    return *this;
}

Point 
NeighborIterator::point() const
{
    return m_dProbC.image(*m_bxit);
}

Box
NeighborIterator::srcBox() const
{
    Box bx(*m_bxit,*m_bxit);
    return (bx.refine(m_layout.boxSize()) & m_bx).shift(m_shift);
}

Box
NeighborIterator::destBox() const
{
    Box bx(*m_bxit,*m_bxit);
    return (bx.refine(m_layout.boxSize()) & m_bx);
}

Point
NeighborIterator::shift()
{
    return m_shift;
}

