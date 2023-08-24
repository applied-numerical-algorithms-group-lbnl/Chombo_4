
MBProblemDomain::MBProblemDomain()
{
    m_domains.clear();
    m_domainsDefined.clear();
}

MBProblemDomain::MBProblemDomain(unsigned int a_numBlocks)
{
    m_domains.clear();
    m_domains.resize(a_numBlocks);
    m_graph = std::make_shared<MBGraph>(a_numBlocks);
    m_domainsDefined.clear();
    m_domainsDefined.resize(a_numBlocks, false);
}

void MBProblemDomain::defineDomain(unsigned int a_blockID, Point a_domainSize)
{
    std::array<bool, DIM> periodicity;
    periodicity.fill(false);
    Box domainBox(a_domainSize);
    m_domains[a_blockID].define(domainBox, periodicity);
    m_domainsDefined[a_blockID] = true;
    validate();
}

void MBProblemDomain::defineBoundary(
        unsigned int        a_srcBlock,
        unsigned int        a_dstBlock,
        unsigned int        a_dir,
        Side::LoHiSide      a_side,
        CoordPermutation&   a_rotation)
{
    m_graph->defineBoundary(a_srcBlock, a_dstBlock, a_dir, a_side, a_rotation);
    validate();
}

void MBProblemDomain::defineBoundary(
        unsigned int        a_srcBlock,
        unsigned int        a_dstBlock,
        Point               a_dir,
        CoordPermutation&   a_rotation)
{
    m_graph->defineBoundary(a_srcBlock, a_dstBlock, a_dir, a_rotation);
    validate();
}

const ProblemDomain& MBProblemDomain::blockDomain(unsigned int a_blockID) const
{
    return m_domains[a_blockID];
}

std::vector<MBGraphArc> MBProblemDomain::boundaries(
        unsigned int a_srcBlock,
        Point        a_dir) const
{
    return m_graph->boundaries(a_srcBlock, a_dir);
}

Point MBProblemDomain::convert(
        Point a_point,
        unsigned int a_srcBlock,
        unsigned int a_dstBlock) const
{
    PROTO_ASSERT(a_srcBlock < size(),
            "MBProblemDomain::convert | Error: srcBlock is out of bounds.");
    PROTO_ASSERT(a_dstBlock < size(),
            "MBProblemDomain::convert | Error: dstBlock is out of bounds.");
    Point n01 = m_graph->connectivity(a_srcBlock, a_dstBlock);
    PROTO_ASSERT(codim(n01) == 1,
            "MBProblemDomain::convert | Error: src and dst blocks are do not \
share a codimension 1 boundary.");
    Box srcBox = m_domains[a_srcBlock].box();
    Box dstBox = m_domains[a_dstBlock].box();
    auto R01 = m_graph->rotation(a_srcBlock, n01);

    Point d01, d10, v01, v10, s10;
    
    // convert the normal n01 to a codimension-DIM boundary in dst
    // which is contained in the codimension-1 boundary n01
    for (int dir = 0; dir < DIM; dir++)
    {
        if (n01[dir] == 0)
        {
            d01[dir] = -1;
        } else {
            d01[dir] = n01[dir];
        }
    }
    //std::cout << "MBProblemDomain::convert(" << a_point << ", " << a_srcBlock << ", " << a_dstBlock << ")" << std::endl; 
    // find the corresponding codimension-DIM boundary in src
    d10 = m_graph->reverseArc(a_srcBlock, a_dstBlock, d01);
    //std::cout << "n01: " << n01 << " | d10: " << d10 << " | d01: " << d01 << std::endl;
    
    // construct position vectors to the common corner defined
    // by d01 and d10
    for (int dir = 0; dir < DIM; dir++)
    {
        v01[dir] = d01[dir] > 0 ? srcBox.size(dir) : 0;
        v10[dir] = d10[dir] > 0 ? dstBox.size(dir) : 0;
    }

    // compute the displacement between origins in dst coordinates
    s10 = v10 - R01(v01);

    // convert the input point to dst coordinates and add to the origin
    // displacement
    //std::cout << "v10: " << v10 << " | v01: " << v01 << " | R(v01): " << R01(v01) << " | s10: " << s10 << std::endl;
    return s10 + R01(a_point);
}

Box MBProblemDomain::convert(
        Box a_box,
        unsigned int a_srcBlock,
        unsigned int a_dstBlock) const
{
    Point L = convert(a_box.low(), a_srcBlock, a_dstBlock);
    Point H = convert(a_box.high() + Point::Ones(), a_srcBlock, a_dstBlock);
    Box B(L, L);
    B = B & H;
    return Box(B.low(), B.high()-Point::Ones());
}

MBProblemDomain MBProblemDomain::coarsen(Point a_refRatio) const
{
    PROTO_ASSERT(coarsenable(a_refRatio),
            "MBProblemDomain::coarsen | Error: Problem domain is not coarsenable.");
    MBProblemDomain coarsened(size());
    for (int ii = 0; ii < size(); ii++)
    {
        coarsened.m_domains[ii] = m_domains[ii].coarsen(a_refRatio);
        coarsened.m_domainsDefined[ii] = true;
    }
    coarsened.m_graph = m_graph;
    return coarsened;
}

MBProblemDomain MBProblemDomain::coarsen(std::vector<Point> a_refRatios) const
{
    MBProblemDomain coarsened(size());
    for (int ii = 0; ii < size(); ii++)
    {
        PROTO_ASSERT(coarsenable(ii, a_refRatios[ii]),
                "MBProblemDomain::coarsen | Error: Problem domain is not coarsenable by inputs.");
        coarsened.m_domains[ii] = m_domains[ii].coarsen(a_refRatios[ii]);
        coarsened.m_domainsDefined[ii] = true;
    }
    coarsened.m_graph = m_graph;
    return coarsened;
}

MBProblemDomain MBProblemDomain::refine(Point a_refRatio) const
{
    MBProblemDomain refined(size());
    for (int ii = 0; ii < size(); ii++)
    {
        refined.m_domains[ii] = m_domains[ii].refine(a_refRatio);
        refined.m_domainsDefined[ii] = true;
    }
    refined.m_graph = m_graph;
    return refined;
}

bool MBProblemDomain::coarsenable(Point a_refRatio) const
{
    for (auto d : m_domains)
    {
        if (!d.coarsenable(a_refRatio)) { return false; }
    }
    return true;
}

bool MBProblemDomain::coarsenable(unsigned int a_block, Point a_refRatio) const
{
    return m_domains[a_block].coarsenable(a_refRatio);
}

bool MBProblemDomain::defined()
{
    for (auto def : m_domainsDefined)
    {
        if (!def) { return false; }
    }
    return true;
}

void MBProblemDomain::validate() const
{
    // TODO: validate conformity of blocks
}
