using std::shared_ptr;
using std::array;

// Default Constructor
DisjointBoxLayout::DisjointBoxLayout()
{
}

// Simple Constructor
DisjointBoxLayout::DisjointBoxLayout(const ProblemDomain   & a_problemDomain, 
                                     const Point           & a_boxSize)
{
    define(a_problemDomain, a_boxSize);
}

// General Constructor
DisjointBoxLayout::DisjointBoxLayout(const ProblemDomain   & a_problemDomain, 
                                     const vector<Point>   & a_coarsenedPatches,
                                     const Point           & a_boxSize)
{
    define(a_problemDomain, a_coarsenedPatches, a_boxSize);
}

// Sub-Region Constructor
DisjointBoxLayout::DisjointBoxLayout(const ProblemDomain   & a_problemDomain, 
                                     const Box             & a_region,
                                     const Point           & a_boxSize)
{
    define(a_problemDomain, a_region, a_boxSize);
}

void DisjointBoxLayout::define(
        std::shared_ptr<BoxPartition> a_partition,
        const Point&  a_boxSize)
{

    m_partition = a_partition;
    m_problemDomain = a_partition->domain().refine(a_boxSize);
    m_boxSize = a_boxSize;
    //m_end = DataIndex<BoxPartition>(m_partition, a_partition->numBoxes());
}

// Copy Constructor
DisjointBoxLayout::DisjointBoxLayout(const DisjointBoxLayout& a_input)
{
    m_partition = a_input.m_partition;
    m_problemDomain = a_input.m_problemDomain;
    m_boxSize = a_input.m_boxSize;
    //m_end = a_input.m_end;
}

// Define (Simple)
void 
DisjointBoxLayout::define(const ProblemDomain   & a_problemDomain, 
                          const Point           & a_boxSize)
{
    Box bxCoarse = a_problemDomain.box().coarsen(a_boxSize);
    BoxIterator bxit(bxCoarse);
    vector<Point> allPoints;
    for (bxit.begin();bxit.ok();++bxit)
    {
        allPoints.push_back(*bxit);
    }
    define(a_problemDomain,allPoints,a_boxSize);
}

// Define (Sub-Domain)
void
DisjointBoxLayout::define(const ProblemDomain   & a_problemDomain,
                          const Box             & a_region,
                          const Point           & a_boxSize)
{
    PROTO_ASSERT(a_problemDomain.box().contains(a_region),
        "DisjointBoxLayout::define | Error: Sub-region is not a subset of the domain.");
    Box patchBox = a_region.coarsen(a_boxSize);
    std::vector<Point> patches;
    for (auto iter = patchBox.begin(); iter.ok(); ++iter)
    {
        patches.push_back(*iter);
    }
    define(a_problemDomain, patches, a_boxSize);
}

// Define (General)
void
DisjointBoxLayout::define(const ProblemDomain   & a_problemDomain,
                          const vector<Point>   & a_patches,
                          const Point           & a_boxSize)
{
    PROTO_ASSERT(a_problemDomain.box().coarsenable(a_boxSize),
        "DisjointBoxLayout::define | Error: Sub-region cannot be tiled by input box size.");    
    auto patchDomain = a_problemDomain.coarsen(a_boxSize);
    m_problemDomain = a_problemDomain;
    m_partition = std::make_shared<BoxPartition>(patchDomain, a_patches);
    m_boxSize = a_boxSize;
    //m_end = DataIndex<BoxPartition>(m_partition, a_patches.size());
}

void
DisjointBoxLayout::loadBalance(
        unsigned int a_startProc,
        unsigned int a_endProc)
{
    std::vector<Point> patches;
    for (auto p : boxes())
    {
        patches.push_back(p.first);
    }
    m_partition->loadBalance(patches, a_startProc, a_endProc);
}

template<typename... Args>
void
DisjointBoxLayout::loadAssign(Args... a_args)
{
    std::vector<Point> patches;
    for (auto p : boxes())
    {
        patches.push_back(p.first);
    }
    m_partition->loadAssign(patches, a_args...);
}

void DisjointBoxLayout::loadAssign(std::vector<std::pair<int, unsigned int>>& a_assignment)
{
    std::vector<Point> patches;
    for (auto p : boxes())
    {
        patches.push_back(p.first);
    }
    m_partition->loadAssign(a_assignment, patches);
}

// Assignment Operator
DisjointBoxLayout&
DisjointBoxLayout::operator=(const DisjointBoxLayout& a_input)
{
    m_partition = a_input.m_partition;
    m_problemDomain = a_input.m_problemDomain;
    m_boxSize = a_input.m_boxSize;
    //m_end = a_input.m_end;
    return *this;
}

// Equality Operator
bool
DisjointBoxLayout::operator==(const DisjointBoxLayout& a_input) const
{
    return ((a_input.m_partition == m_partition)&&
            (a_input.m_problemDomain == m_problemDomain)&&
            (a_input.m_boxSize == m_boxSize));
}

const std::vector<std::pair<Point, unsigned int>>
DisjointBoxLayout::boxes() const
{
    return m_partition->partition(); 
}

unsigned int
DisjointBoxLayout::offset(int a_proc) const
{
    return m_partition->procStartIndex(a_proc);
}

///  
// number of boxes in grid (over all procs)
unsigned int
DisjointBoxLayout::size() const
{
    return m_partition->numBoxes();
}  

unsigned int
DisjointBoxLayout::localSize() const
{
    return m_partition->numBoxes(Proto::procID());
}

Box
DisjointBoxLayout::operator[](const DataIndex<BoxPartition>& a_index) const
{
    PROTO_ASSERT(compatible(a_index),
        "DisjointBoxLayout::operator[] | Error: \
        Input DataIndex is not compatible with this layout.");
    return box(a_index);
}

Point
DisjointBoxLayout::point(const DataIndex<BoxPartition>& a_index) const
{
    PROTO_ASSERT(compatible(a_index),
        "DisjointBoxLayout::point | Error: Input DataIndex is not compatible with this layout.");
    return boxes()[a_index.global()].first;
}

Box
DisjointBoxLayout::box(const DataIndex<BoxPartition>& a_index) const
{
    PROTO_ASSERT(compatible(a_index),
        "DisjointBoxLayout::box | Error: Input DataIndex is not compatible with this layout.");
    Point p = point(a_index);
    //return Box(p*m_boxSize,(p + Point::Ones())*m_boxSize - Point::Ones());
    return Box(p,p).refine(m_boxSize);
}

int
DisjointBoxLayout::procID(const DataIndex<BoxPartition>& a_index) const
{
    PROTO_ASSERT(compatible(a_index),
        "DisjointBoxLayout::procID | Error: Input DataIndex is not compatible with this layout.");
    return boxes()[a_index.global()].second;
}

DataIndex<BoxPartition>
DisjointBoxLayout::index(Point& a_patchPoint) const
{
    return DataIndex<BoxPartition>(m_partition, m_partition->find(a_patchPoint));
}

DataIndex<BoxPartition>
DisjointBoxLayout::index(unsigned int a_indexInt) const
{
    PROTO_ASSERT(a_indexInt < size(),
        "DisjointBoxLayout::index | Error: \
        Global integer index %u is out of bounds.", a_indexInt);
    return DataIndex<BoxPartition>(m_partition,a_indexInt);
}

DataIndex<BoxPartition>
DisjointBoxLayout::localIndex(unsigned int a_indexInt) const
{
    PROTO_ASSERT(a_indexInt < localSize(),
        "DisjointBoxLayout::index | Error: \
        Global integer index %u is out of bounds.", a_indexInt);
    return DataIndex<BoxPartition>(m_partition, m_partition->procStartIndex(Proto::procID()) + a_indexInt);
}

DataIndex<BoxPartition>
DisjointBoxLayout::find(const Point& a_pt) const
{
    if (!patchDomain().contains(a_pt))
    {
        return *end();
    }
    auto globalIndex = m_partition->find(a_pt);
    if (globalIndex == size())
    {
        return *end();
    } else {
        return DataIndex<BoxPartition>(m_partition, globalIndex);
    }
}

bool
DisjointBoxLayout::contains(const Point& a_pt) const
{
    auto findIndex = find(a_pt);
    auto endIndex = *end();
    return find(a_pt) != *end();
}

bool
DisjointBoxLayout::onDomainBoundary(const Point& a_pt) const
{
    PROTO_ASSERT(find(a_pt) != *end(),
        "DisjointBoxLayout::onDomainBoundary | Error: \
        Point is not in this layout.");
    Box K = Box::Kernel(1).shift(a_pt);
    
    for (auto biter = K.begin(); biter != K.end(); ++biter)
    {
        if (*biter == a_pt) {continue;}
        if (!patchDomain().contains(*biter)) { return true; }
    }
    return false;
}

bool
DisjointBoxLayout::onLevelBoundary(const Point& a_pt) const
{
    PROTO_ASSERT(find(a_pt) != *end(),
        "DisjointBoxLayout::onDomainBoundary | Error: \
        Point is not in this layout.");
    Box K = Box::Kernel(1).shift(a_pt);
    
    for (auto biter = K.begin(); biter != K.end(); ++biter)
    {
        if (*biter == a_pt) {continue;}
        if (find(*biter) == *end()) { return true; }
    }
    return false;
}

DataIterator<BoxPartition>
DisjointBoxLayout::end() const
{
    DataIterator<BoxPartition> iter(m_partition);
    return iter.end();
}

DataIterator<BoxPartition>
DisjointBoxLayout::begin() const
{
    DataIterator<BoxPartition> iter(m_partition);
    return iter.begin();
}

bool 
DisjointBoxLayout::coarsenable(const Point& a_refRatio) const
{
    return (((m_problemDomain.box().coarsenable(a_refRatio))) &&
            ((m_boxSize%a_refRatio) == (Point::Zeros())));
}

bool 
DisjointBoxLayout::compatible(const DisjointBoxLayout& a_layout) const
{
    return m_partition->compatible(*a_layout.m_partition);
}

bool 
DisjointBoxLayout::compatible(const DataIndex<BoxPartition>& a_index) const
{
    return a_index.compatible(*m_partition);
}

bool 
DisjointBoxLayout::compatible(const DataIterator<BoxPartition>& a_iter) const
{
    return a_iter.compatible(*m_partition);
}

DisjointBoxLayout 
DisjointBoxLayout::coarsen(const Point& a_refRatio) const
{
    PROTO_ASSERT(coarsenable(a_refRatio),
        "DisjointBoxLayout::coarsen | Error: Invalid coarsening ratio (not coarsenable)");
    DisjointBoxLayout retval;
    retval.m_partition = m_partition;
    retval.m_problemDomain = m_problemDomain.coarsen(a_refRatio);
    retval.m_boxSize = m_boxSize/a_refRatio;
    return retval;
} 

DisjointBoxLayout 
DisjointBoxLayout::coarsen(const Point& a_refRatio, const Point& a_boxSize) const
{
    // The simple case where there are the same number of boxes 
    if (boxSize() / a_refRatio == a_boxSize)
    {
        return coarsen(a_refRatio); 
    }
    MayDay<void>::Error("DisjointBoxLayout::coarsen | Error: Function not implemented.");
    PROTO_ASSERT(coarsenable(a_refRatio),
        "DisjointBoxLayout::coarsen | Error: Invalid coarsening ratio (not coarsenable)");
    for (auto iter : m_partition->partition())
    {
        Point p = iter.first;
    }
    return DisjointBoxLayout();
}

DisjointBoxLayout 
DisjointBoxLayout::refine(const Point& a_refRatio) const
{
    DisjointBoxLayout retval;
    retval.m_partition = m_partition;
    retval.m_problemDomain = m_problemDomain.refine(a_refRatio);
    retval.m_boxSize = m_boxSize*a_refRatio;
    return retval;    
}

Box
DisjointBoxLayout::boundingBox() const
{
    if (size() == 0) {return Box(); }
    Box boundingBox;
    for (auto iter : boxes())
    {
        Point p = iter.first;
        if (boundingBox.empty())
        {
            boundingBox = Box(p, p);
        } else {
            boundingBox &= p;
        }
    }
    return boundingBox.refine(boxSize());
}

bool
DisjointBoxLayout::radialSymmetry(Point a_origin) const
{
    std::vector<Point> assymmetries; 
    for (auto iter : boxes())
    {
        Point shifted = iter.first - a_origin;
        Point mirrorImage = -(shifted + Point::Ones());
        Point searchPoint = mirrorImage + a_origin;
        if (find(searchPoint) == *end())
        {
            assymmetries.push_back(iter.first);
        }
    }
    if (assymmetries.size() == 0){ return true; }
    else {return false;}
}

bool
DisjointBoxLayout::mirrorSymmetry(Point a_origin, int a_coord) const
{
    if (a_coord < 0 || a_coord >= DIM) { return true; }
    
    std::vector<Point> assymmetries; 
    for (auto iter : boxes())
    {
        Point shifted = iter.first - a_origin;
        Point mirrorImage = shifted;
        mirrorImage[a_coord] = -(mirrorImage[a_coord] + 1);
        Point searchPoint = mirrorImage + a_origin;
        if (find(searchPoint) == *end())
        {
            assymmetries.push_back(iter.first);
        }
    }
    if (assymmetries.size() == 0){return true;}
    else {return false;}
}

void
DisjointBoxLayout::print(std::string a_name) const
{
    using Proto::pout;

    pout() << "Disjoint Box Layout: " << a_name << std::endl;
    pout() << std::setw(40) << std::setfill('-') << "-" << std::endl;
    pout() << std::setfill(' ');
    pout() << std::setw(10) << std::left << "Domain: "   << m_problemDomain.box() << std::endl;
    pout() << std::setw(10) << std::left << "Box Size: " << m_boxSize << std::endl;
    pout() << std::setw(40) << std::setfill('-') << "-" << std::endl;
    pout() << std::setfill(' ');
    if (size() == 0)
    {
        pout() << "( EMPTY LAYOUT ) " << std::endl;
        pout() << setfill('\0');
        return;
    }
    pout() << std::right;
    for (int ii = 0; ii < size(); ii++)
    {
        Point p = boxes()[ii].first;
        Box b = Box(p,p).refine(m_boxSize);
        pout() << std::setw(4) << boxes()[ii].second << " | ";
        pout() << b << " | ";
        pout() << p << std::endl;
    }
    pout() << setfill(' ') << std::endl;
}

inline std::ostream& operator<< (std::ostream& os, const DisjointBoxLayout& a_dbl)
{
    os << "DisjointBoxLayout: \n" ;
    os << a_dbl.domain() << "\n";
    os << "coarsened Domain = " << a_dbl.patchDomain() << "\n" ;
    os << "Patch Partition: \n";
    for (int k = 0; k < a_dbl.size(); k++)
    {
        os << a_dbl.boxes()[k].first 
            << " , " << a_dbl.boxes()[k].second << "\n";
    }
    return os;
}



