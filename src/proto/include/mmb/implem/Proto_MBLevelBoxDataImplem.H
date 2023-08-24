
template<typename T, unsigned int C, MemType MEM, Centering CTR>
MBLevelBoxData<T, C, MEM, CTR>::MBLevelBoxData(
        MBDisjointBoxLayout& a_layout,
        Point a_ghost)
{
    std::array<Point, DIM+1> ghost;
    ghost.fill(a_ghost);
    define(a_layout, ghost);
}

template<typename T, unsigned int C, MemType MEM, Centering CTR>
MBLevelBoxData<T, C, MEM, CTR>::MBLevelBoxData(
        MBDisjointBoxLayout&     a_layout,
        std::array<Point, DIM+1> a_ghost)
{
    define(a_layout, a_ghost);
}

template<typename T, unsigned int C, MemType MEM, Centering CTR>
void MBLevelBoxData<T, C, MEM, CTR>::define(
        MBDisjointBoxLayout&      a_layout,
        std::array<Point, DIM+1> a_ghost)
{
    // a_ghost convention:
    // [interior, codim1, codim2, ...]
    m_ghost = a_ghost;
    m_layout = &a_layout;
    m_data.clear();
    m_data.resize(a_layout.numBlocks());
    for (int bi = 0; bi < a_layout.numBlocks(); bi++)
    {
        m_data[bi] = std::make_shared<LevelBoxData<T, C, MEM, CTR>>(
                a_layout.layout(bi), a_ghost[0]);
    }
    for (auto iter : a_layout)
    {
        auto block = a_layout.block(iter);
        auto patch = a_layout.point(iter);
        auto localLayout = a_layout.layout(block);
        auto patchDomain = localLayout.patchDomain().box();
        auto localIndex = localLayout.find(patch);
        auto patchBox = localLayout[localIndex];
        
        //std::cout << std::endl << "Getting boundaries for patch " << patch << " in block " << block << std::endl;

        Box K = Box::Kernel(1);
        for (int cc = 1; cc <= DIM; cc++)
        {
            //std::cout << "\tChecking codimension " << cc << " boundaries." << std::endl;
            for (auto di : K)
            {
                if (codim(di) != cc) { continue; }
                if (patchDomain.contains(patch+di)) { continue; }
                Point shrinkDirs;
                for (int d = 0; d < DIM; d++)
                {
                    shrinkDirs[d] = (di[d] != 0) ? -1 : 0;
                }
                Box k = K.grow(shrinkDirs).shift(di);
                PROTO_ASSERT(k.contains(di),
                        "MBLevelBoxData::define | Error: Data corruption");
                auto bounds = a_layout.domain().boundaries(block, di);
                Box boundaryPatches = patchDomain.adjacent(di,1);
                std::cout << "\tdir: " << di << " | k: " << k << " | boundaryPatches: " << boundaryPatches << std::endl;
                for (auto bound : bounds)
                {
                    PROTO_ASSERT(bound.srcToDst == di,
                        "MBLevelBoxData::define | Error: Data corruption");
                    auto adjBlock = bound.dstBlock;
                    Point adjDir = bound.dstToSrc;
                    auto adjLocalLayout = a_layout.layout(adjBlock);
                    Box adjPatchDomain = adjLocalLayout.patchDomain().box();
                    Box adjBoundaryPatches = adjPatchDomain.edge(adjDir);
                    auto R = a_layout.domain().graph().rotation(block, di, adjBlock);
                    //std::cout << "\t\tadjBlock: " << adjBlock << " | adjDir: " << adjDir << " | adjBoundaryPatches: " << adjBoundaryPatches << std::endl;
                    for (auto ki : k)
                    {
                        boundsKey_t key(patch, block, ki);
                        Point neighbor = patch + ki;
                        //std::cout << "\t\t\tChecking ki: " << ki << " | neighbor: " << neighbor << std::endl;
                        if (boundaryPatches.contains(neighbor))
                        {
                            Point L = boundaryPatches.low();
                            auto adjPatch = R.rotateCell(neighbor - L, boundaryPatches.shift(-L))
                                + adjBoundaryPatches.low();
                            auto adjIndex = adjLocalLayout.find(adjPatch);
                            auto kiAdj = R(ki);
                            // FIXME: this should use anisotropic ghost value for codim cc
                            //std::cout << "\t\t\t\tPatchBox: " << patchBox << " | patchBox.adjacent(ki,1): " << patchBox.adjacent(ki,1) << std::endl;
                            Box patchBoundary = patchBox.adjacent(ki, 1);
                            Box adjPatchBoundary = adjLocalLayout[adjIndex].edge(-kiAdj,1);
                            
                            //std::cout << "\t\t\t\tFound boundary | adjPatch: " << adjPatch << " | patchBoundary: " << patchBoundary << " | adjPatchBoundary: " << adjPatchBoundary << std::endl;
                            MBBoundaryData<T, C, MEM> boundData(
                                    block, adjBlock, patchBoundary, adjPatchBoundary);
                            m_bounds[key].push_back(boundData);
                        }
                    }
                }
            }
        }
    }
}

template<typename T, unsigned int C, MemType MEM, Centering CTR>
template<typename Func, typename... Srcs>
void MBLevelBoxData<T, C, MEM, CTR>::initialize(Func& a_func, Srcs... a_srcs) {}

template<typename T, unsigned int C, MemType MEM, Centering CTR>
template<typename Func, typename... Srcs>
void MBLevelBoxData<T, C, MEM, CTR>::initConvolve(Func& a_func, Srcs... a_srcs){}

template<typename T, unsigned int C, MemType MEM, Centering CTR>
BoxData<T, C, MEM>&
MBLevelBoxData<T, C, MEM, CTR>::operator[](index_t a_index)
{
    unsigned int block = m_layout->block(a_index);
    auto localIndex = m_layout->blockIndex(a_index);
    return (*m_data[block])[localIndex];
}

template<typename T, unsigned int C, MemType MEM, Centering CTR>
const BoxData<T, C, MEM>&
MBLevelBoxData<T, C, MEM, CTR>::operator[](index_t a_index) const
{
    unsigned int block = m_layout->block(a_index);
    auto localIndex = m_layout->blockIndex(a_index);
    return (*m_data[block])[localIndex];
}

template<typename T, unsigned int C, MemType MEM, Centering CTR>
LevelBoxData<T, C, MEM, CTR>&
MBLevelBoxData<T, C, MEM, CTR>::blockData(unsigned int a_block)
{
    return m_data[a_block];
}

template<typename T, unsigned int C, MemType MEM, Centering CTR>
const LevelBoxData<T, C, MEM, CTR>&
MBLevelBoxData<T, C, MEM, CTR>::blockData(unsigned int a_block) const
{
    return m_data[a_block];
}

template<typename T, unsigned int C, MemType MEM, Centering CTR>
std::vector<MBBoundaryData<T, C, MEM>>&
MBLevelBoxData<T, C, MEM, CTR>::bounds(index_t a_index, Point a_dir)
{
    auto patch = m_layout->point(a_index);
    auto block = m_layout->block(a_index);
    return m_bounds[boundsKey_t(patch, block, a_dir)];
}

template<typename T, unsigned int C, MemType MEM, Centering CTR>
const std::vector<MBBoundaryData<T, C, MEM>>&
MBLevelBoxData<T, C, MEM, CTR>::bounds(index_t a_index, Point a_dir) const
{
    auto patch = m_layout->point(a_index);
    auto block = m_layout->block(a_index);
    return m_bounds[boundsKey_t(patch, block, a_dir)];
}

template<typename T, unsigned int C, MemType MEM, Centering CTR>
typename MBLevelBoxData<T, C, MEM, CTR>::iterator_t
MBLevelBoxData<T, C, MEM, CTR>::begin()
{
    return m_layout->begin();
}

template<typename T, unsigned int C, MemType MEM, Centering CTR>
typename MBLevelBoxData<T, C, MEM, CTR>::iterator_t
MBLevelBoxData<T, C, MEM, CTR>::end()
{
    return m_layout->end();
}
