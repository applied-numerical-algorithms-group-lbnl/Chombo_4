MBDisjointBoxLayout::MBDisjointBoxLayout(
        const MBProblemDomain& a_domain,
        Point a_boxSize)
{
    std::vector<Point> boxSizes(a_domain.numBlocks(), a_boxSize);
    define(a_domain, boxSizes);
}

MBDisjointBoxLayout::MBDisjointBoxLayout(
        const MBProblemDomain& a_domain,
        const std::vector<Point>& a_boxSizes)
{
    define(a_domain, a_boxSizes);
}

MBDisjointBoxLayout::MBDisjointBoxLayout(
        const MBProblemDomain& a_domain,
        const std::vector<MBPatchID_t>& a_patches,
        const std::vector<Point>& a_boxSizes)
{
    define(a_domain, a_patches, a_boxSizes);
}

void MBDisjointBoxLayout::define(
        const MBProblemDomain& a_domain, 
        const std::vector<Point>& a_boxSizes)
{
    PROTO_ASSERT(a_boxSizes.size() == a_domain.size(),
            "MBDisjointBoxLayout::define | Error: incorrect number of boxSize inputs.");
    
    std::vector<MBPatchID_t> patches;
    for (unsigned int bi = 0; bi < a_domain.numBlocks(); bi++)
    {
        Box patchBox = a_domain.blockDomain(bi).box().coarsen(a_boxSizes[bi]);
        for (auto pi : patchBox)
        {
            patches.push_back(MBPatchID_t(pi, bi));
        }
    }
    define(a_domain, patches, a_boxSizes);
}

void MBDisjointBoxLayout::define(
        const MBProblemDomain& a_domain,
        const std::vector<MBPatchID_t>& a_patches,
        const std::vector<Point>& a_boxSizes)
{
    PROTO_ASSERT(a_boxSizes.size() == a_domain.size(),
            "MBDisjointBoxLayout::define | Error: incorrect number of boxSize inputs.");

    m_domain = a_domain;
    m_partition = std::make_shared<MBBoxPartition>(a_domain.coarsen(a_boxSizes), a_patches);
    m_layouts.clear();
    m_layouts.resize(a_domain.numBlocks());
    
    for (unsigned int bi = 0; bi < a_domain.numBlocks(); bi++)
    {
        m_layouts[bi].define(m_partition->blockPartition(bi), a_boxSizes[bi]);
    }
}

MBProblemDomain MBDisjointBoxLayout::domain() const
{
    return m_domain;
}

unsigned int MBDisjointBoxLayout::numBlocks() const
{
    return m_partition->numBlocks();
}

unsigned int MBDisjointBoxLayout::numBoxes() const 
{
    return m_partition->numBoxes();
}

unsigned int MBDisjointBoxLayout::numBoxes(unsigned int a_proc) const 
{
    return m_partition->numBoxes(a_proc);
}

unsigned int MBDisjointBoxLayout::numProcs() const
{
    return m_partition->numProcs();
}

DataIndex<BoxPartition>
MBDisjointBoxLayout::blockIndex(const index_t& a_index) const
{
    Point p = point(a_index);
    unsigned int b = block(a_index);
    return m_layouts[b].find(p); 
}

Box MBDisjointBoxLayout::operator[](const index_t& a_index)
{
    return box(a_index);
}

MBDisjointBoxLayout::iterator_t MBDisjointBoxLayout::begin()
{
    iterator_t iter(m_partition);
    return iter.begin();
}

MBDisjointBoxLayout::iterator_t MBDisjointBoxLayout::end()
{
    iterator_t iter(m_partition);
    return iter.end();
}

Point MBDisjointBoxLayout::point(const index_t& a_index) const
{
    auto data = m_partition->partition()[a_index.global()];
    auto block = data.first.second;
    auto patch = data.first.first;
    return patch;
}

Box MBDisjointBoxLayout::box(const index_t& a_index) const
{
    auto data = m_partition->partition()[a_index.global()];
    auto block = data.first.second;
    auto patch = data.first.first;
    auto boxSize = m_layouts[block].boxSize();
    return Box(patch, patch).refine(boxSize);
}

unsigned int MBDisjointBoxLayout::block(const index_t& a_index) const
{
    auto data = m_partition->partition()[a_index.global()];
    auto block = data.first.second;
    return block;
}

DisjointBoxLayout& MBDisjointBoxLayout::layout(unsigned int a_block)
{
    PROTO_ASSERT(a_block < numBlocks(),
            "MBDisjointBoxLayout::layout | Error: Block index is out of bounds \
            (numBlocks() = %u)", numBlocks());
    return m_layouts[a_block];
}
