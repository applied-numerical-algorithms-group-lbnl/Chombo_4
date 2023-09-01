BoxPartition::BoxPartition(
        const ProblemDomain& a_patchDomain,
        const std::vector<Point>& a_patches)
{
    define(a_patchDomain, a_patches, 0, numProc());
}

BoxPartition::BoxPartition(
        const ProblemDomain& a_patchDomain,
        const std::vector<Point>& a_patches,
        unsigned int a_startProc,
        unsigned int a_endProc)
{
    define(a_patchDomain, a_patches, a_startProc, a_endProc);
}

void BoxPartition::define(
        const ProblemDomain& a_patchDomain,
        const std::vector<Point>& a_patches,
        unsigned int a_startProc,
        unsigned int a_endProc)
{
    for (auto patch : a_patches)
    {
        PROTO_ASSERT(a_patchDomain.box().contains(patch),
                "BoxPartion::Constructor | Error: \
                All patches must be contained in the patchDomain.");
    }
    
    m_patchDomain = a_patchDomain;
    std::vector<Point> sortedPatches = a_patches;
    Morton::sort(sortedPatches);
    loadBalance(sortedPatches, a_startProc, a_endProc);
}

void BoxPartition::assign(
        std::vector<Point>& a_patches,
        unsigned int a_globalIndex,
        int          a_proc,
        unsigned int a_num)
{
    int global = a_globalIndex;
    for (int ii = 0; ii < a_num; ii++, global++)
    {
        Point patch = a_patches[global];
        uint64_t mortonIndex = Morton::index(patch);
        m_indexMap.insert(std::pair<uint64_t, int>(mortonIndex, global));
        auto procAssign = pair<Point, unsigned int>(patch, a_proc);
        m_partition.push_back(procAssign);
    }
}

template<typename... vals>
void BoxPartition::loadAssign(
        std::vector<Point>& a_patches,
        vals... a_args)
{
    m_procMap.clear();
    m_indexMap.clear();
    m_partition.clear();
    unpack(a_patches, 0, a_args...);
}

void BoxPartition::loadAssign(
        std::vector<std::pair<int, unsigned int>>& a_assignment)
{
    std::vector<Point> patches;
    for (auto data : m_partition)
    {
        patches.push_back(data.first);
    }
    loadAssign(a_assignment, patches);
}
void BoxPartition::loadAssign(
        std::vector<std::pair<int, unsigned int>>& a_assignment,
        std::vector<Point>& a_patches)
{
    m_procMap.clear();
    m_indexMap.clear();
    m_partition.clear();
    int globalIndex = 0;
    for (auto item : a_assignment)
    {
        int proc = item.first;
        unsigned int numBoxes = item.second;
        m_procMap[proc] = std::pair<int, int>(globalIndex, globalIndex + numBoxes);
        assign(a_patches, globalIndex, proc, numBoxes);
        globalIndex += numBoxes;
    }
    if (globalIndex < a_patches.size())
    {
        MayDay<void>::Warning("BoxPartition::loadAssign | Warning: \
                Incomplete load assignment.");
    }
}

void BoxPartition::loadBalance(
        std::vector<Point>& a_patches,
        unsigned int a_startProc,
        unsigned int a_endProc)
{
    m_procMap.clear();
    m_indexMap.clear();
    m_partition.clear();
    if (a_startProc >= numProc() || a_startProc >= a_endProc) { return; } //nothing else to do
    
    int nsegs  = a_endProc - a_startProc;
    int length = a_patches.size() / nsegs;
    int remainder = a_patches.size() % nsegs;
    int globalIndex = 0;
    for (int procNum = a_startProc; procNum < a_startProc + nsegs; procNum++)
    {
        int segLength = length;
        if (procNum < remainder)
        {
            segLength = length+1;
        }
        m_procMap[procNum] = std::pair<int, int>(globalIndex, globalIndex + segLength);
        assign(a_patches, globalIndex, procNum, segLength);
        globalIndex += segLength;
    }
}


bool BoxPartition::compatible(const BoxPartition& a_rhs)
{
    if (this == &a_rhs) { return true; }
    auto lhsPartition = m_partition;
    auto rhsPartition = a_rhs.partition();
    if (lhsPartition.size() != rhsPartition.size())
    {
        return false;
    } 
    for (int ii = 0; ii < lhsPartition.size(); ii++)
    {
        auto& lhs_i = lhsPartition[ii];
        auto& rhs_i = rhsPartition[ii];
        if (lhs_i.first != rhs_i.first) 
        {
            return false;
        }
        if (lhs_i.second != rhs_i.second) 
        {
            return false; 
        }
    }
    return true;
}

unsigned int BoxPartition::numProcs() const
{
    return m_procMap.size();
}

unsigned int BoxPartition::numBoxes() const
{
    return m_partition.size();
}

unsigned int BoxPartition::numBoxes(unsigned int a_proc) const
{
    return procEndIndex(a_proc) - procStartIndex(a_proc);
}

unsigned int BoxPartition::procStartIndex(unsigned int a_proc) const
{
    auto data = m_procMap.find(a_proc);
    if (data == m_procMap.end()) { return numBoxes(); }
    return (*data).second.first;
}

unsigned int BoxPartition::procEndIndex(unsigned int a_proc) const
{
    auto data = m_procMap.find(a_proc);
    if (data == m_procMap.end()) { return numBoxes(); }
    return (*data).second.second;
}

/*
unsigned int BoxPartition::patchIndex(Point a_pt)
{
    PROTO_ASSERT(m_patchDomain.contains(a_pt),
            "BoxPartition::patchIndex | Error: \
            The input patch is not contained in the domain or any of its periodic images.");
    Point image = m_patchDomain.image(a_pt);
    return m_indexMap[Morton::index(image)];
}
*/

unsigned int BoxPartition::find(Point a_pt)
{
    PROTO_ASSERT(m_patchDomain.contains(a_pt),
            "BoxPartition::patchIndex | Error: \
            The input patch is not contained in the domain or any of its periodic images.");
    Point image = m_patchDomain.image(a_pt);
    auto mortonIndex = Morton::index(image);
    auto iter = m_indexMap.find(mortonIndex);
    if (iter == m_indexMap.end())
    {
        return numBoxes();
    }
    return m_indexMap[mortonIndex];
}

void BoxPartition::print() const
{
    pout() << "Box Partition: " << std::endl;
    pout() << "\tPatch Domain: " << m_patchDomain.box() << std::endl;
    pout() << "\tNum Boxes: " << numBoxes() << std::endl;
    pout() << "\tNum Procs: " << numProcs() << std::endl;
    pout() << "\tData By Processor: " << std::endl;
    for (auto procData : m_procMap)
    {
        unsigned int proc = procData.first;
        pout() << "\t\tProc: " << proc << " | NumBoxes: " << numBoxes(proc);
        pout() << " | Start Index: " << procStartIndex(proc);
        pout() << " | End Index: " << procEndIndex(proc) << std::endl;
    }
    pout() << "\tData Partition: " << std::endl;
    for (auto item : m_partition)
    {
        pout() << "\t\tPatch: " << item.first << " | Proc: " << item.second << std::endl; 
    }
}

template<typename... vals>
inline void BoxPartition::unpack(
        std::vector<Point>& a_patches,
        unsigned int a_globalIndex,
        int          a_proc,
        unsigned int a_num,
        vals... a_args)
{
    m_procMap[a_proc] = std::pair<int, int>(a_globalIndex, a_globalIndex + a_num);
    assign(a_patches, a_globalIndex, a_proc, a_num); 
    unpack(a_patches, a_globalIndex + a_num, a_args...);
}

template<>
inline void BoxPartition::unpack(
        std::vector<Point>& a_patches,
        unsigned int a_globalIndex,
        int          a_proc,
        unsigned int a_num)
{
    m_procMap[a_proc] = std::pair<int, int>(a_globalIndex, a_globalIndex + a_num);
    assign(a_patches, a_globalIndex, a_proc, a_num); 
}
