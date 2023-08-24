MBBoxPartition::MBBoxPartition(
        const MBProblemDomain& a_patchDomain,
        const std::vector<MBPatchID_t>& a_patches)
{
    define(a_patchDomain, a_patches);
}

void MBBoxPartition::define(
        const MBProblemDomain& a_patchDomain,
        const std::vector<MBPatchID_t>& a_patches)
{
    m_patchDomain = a_patchDomain;
    std::vector<std::vector<Point>> blockPatches;
    blockPatches.resize(a_patchDomain.numBlocks());
    for (auto pi : a_patches)
    {
        blockPatches[pi.second].push_back(pi.first);
    }
    
    m_blockPartitions.clear();
    for (unsigned int bi = 0; bi < a_patchDomain.size(); bi++)
    {
        auto partition = std::make_shared<BoxPartition>(a_patchDomain.blockDomain(bi), blockPatches[bi]);
        m_blockPartitions.push_back(partition);
    }
    loadBalance();
}

void MBBoxPartition::buildLocalMaps()
{
    uint64_t N = numBoxes();
    uint64_t n = N / Proto::numProc();
    uint64_t r = N % Proto::numProc();
    std::vector<uint64_t> unassignedBoxes(Proto::numProc(), n);
    for (int ii = 0; ii < r; ii++)
    {
        unassignedBoxes[ii]++;
    }
    unsigned int proc = 0;
    for (unsigned int bi = 0; bi < numBlocks(); bi++)
    {
        uint64_t numAssigned = 0;
        uint64_t blockSize = m_blockPartitions[bi]->numBoxes();
        std::vector<std::pair<int, unsigned int>> blockAssignment;
        while (numAssigned < blockSize)
        {
            uint64_t ni = min(unassignedBoxes[proc], blockSize - numAssigned);
            blockAssignment.push_back(std::pair<int, unsigned int>(proc,ni));
            unassignedBoxes[proc] -= ni;
            numAssigned += ni;
            if (unassignedBoxes[proc] == 0) { proc++; }
        }
        m_blockPartitions[bi]->loadAssign(blockAssignment);
    }
}

void MBBoxPartition::buildGlobalMaps()
{
    m_partition.clear();
    m_indexMap.clear();
    m_procMap.clear();
    uint64_t globalIndex = 0;
    unsigned int currentProc = 0;
    uint64_t procRangeStart = 0;
    for (unsigned int bi = 0; bi < numBlocks(); bi++)
    {
        auto patchData = m_blockPartitions[bi]->partition();
        for (auto patchData_i : patchData)
        {
            Point patch = patchData_i.first;
            unsigned int  proc  = patchData_i.second;
            if (proc != currentProc)
            {
                std::pair<uint64_t, uint64_t> procRange(procRangeStart, globalIndex);
                m_procMap[currentProc] = procRange;
                currentProc = proc;
                procRangeStart = globalIndex;
            }
            MBPatchID_t mbPatch(patch, bi);
            m_partition.push_back(std::pair<MBPatchID_t, unsigned int>(mbPatch, proc));
            m_indexMap[mbPatch] = globalIndex;
            globalIndex++;
        }
        std::pair<uint64_t, uint64_t> procRange(procRangeStart, numBoxes());
        m_procMap[currentProc] = procRange;
    }
}

void MBBoxPartition::loadBalance()
{
    buildLocalMaps();
    buildGlobalMaps();
}

bool MBBoxPartition::compatible(const MBBoxPartition& a_rhs)
{
    if (this == &a_rhs) { return true; }
    if (numBlocks() != a_rhs.numBlocks()) { return false; }
    for (unsigned int bi = 0; bi < numBlocks(); bi++)
    {
        if (!m_blockPartitions[bi]->compatible(*a_rhs.m_blockPartitions[bi]))
        {
            return false;
        }
    }
    return true;
}
unsigned int MBBoxPartition::numProcs() const
{
    return m_procMap.size();
}
unsigned int MBBoxPartition::numBoxes() const
{
    int N = 0;
    for (auto p : m_blockPartitions)
    {
        N += p->numBoxes();
    }
    return N;
}
unsigned int MBBoxPartition::numBoxes(unsigned int a_proc)
{
    return m_procMap[a_proc].second - m_procMap[a_proc].first;
}

unsigned int MBBoxPartition::numBlocks() const
{
    return m_blockPartitions.size();
}

uint64_t MBBoxPartition::procStartIndex(unsigned int a_proc) const
{
    auto data = m_procMap.find(a_proc);
    if (data == m_procMap.end()) { return numBoxes(); }
    return (*data).second.first;
}
uint64_t MBBoxPartition::procEndIndex(unsigned int a_proc) const
{
    auto data = m_procMap.find(a_proc);
    if (data == m_procMap.end()) { return numBoxes(); }
    return (*data).second.second;
}


inline std::shared_ptr<BoxPartition> MBBoxPartition::blockPartition(unsigned int a_block) const
{
    return m_blockPartitions[a_block];
}

uint64_t MBBoxPartition::find(Point a_patch, unsigned int a_block)
{
    auto blockDomain = m_patchDomain.blockDomain(a_block);
    PROTO_ASSERT(blockDomain.contains(a_patch),
            "MBBoxPartition::find | Error: \
            Patch is not contained within the problem domain.");
    Point image = blockDomain.image(a_patch);
    MBPatchID_t patchID(image, a_block);
    auto data = m_indexMap.find(patchID);
    if (data == m_indexMap.end()) { return numBoxes(); }
    return m_indexMap[patchID];
}


void MBBoxPartition::print()
{
    pout() << "MBBoxPartition" << std::endl;
    for (unsigned int bi = 0; bi < numBlocks(); bi++)
    {
        pout() << "Block " << bi << std::endl;
        m_blockPartitions[bi]->print();
    }
    pout() << "Global Partition Data" << std::endl;
    pout() << "\tNum Boxes:  " << numBoxes() << std::endl;
    pout() << "\tNum Blocks: " << numBlocks() << std::endl; 
    pout() << "\tNum Procs:  " << numProcs() << std::endl;
    pout() << "\tData By Processor: " << std::endl;
    for (auto procData : m_procMap)
    {
        unsigned int proc = procData.first;
        pout() << "\t\tProc: " << proc << " | Num Boxes: " << numBoxes(proc);
        pout() << " | Start Index: " << procStartIndex(proc);
        pout() << " | End Index: " << procEndIndex(proc) << std::endl;
    }
    pout() << "\tData Partition: " << std::endl;
    for (auto patchData : m_partition)
    {
        auto patchID = patchData.first;
        auto proc = patchData.second;
        pout() << "\t\tBlock: " << patchID.second << " | Patch: " << patchID.first;
        pout() << " | Proc: " << proc << std::endl;
    }
}

