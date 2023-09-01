
FinitePointSet::FinitePointSet(Box a_domain, std::array<bool, DIM> a_periodic)
{
    m_domain = ProblemDomain(a_domain, a_periodic);
#ifdef PR_MPI
    m_sendBufferSize = 0;
    m_recvBufferSize = 0;
    m_sendBuffer = nullptr;
    m_recvBuffer = nullptr;
    
    m_globalBufferSize = 0;
    m_localBufferSize = 0;
    //m_bufferSizes = (int*)malloc(numProc()*sizeof(int));
    //m_bufferOffsets = (int*)malloc(numProc()*sizeof(int));
    protoMalloc(MEMTYPE_DEFAULT, m_bufferSizes, numProc()*sizeof(int)); 
    protoMalloc(MEMTYPE_DEFAULT, m_bufferOffsets, numProc()*sizeof(int)); 
    m_localBuffer = nullptr;
    m_globalBuffer = nullptr;
#endif
}

FinitePointSet::FinitePointSet(Box a_domain, bool a_periodic)
{
    m_domain = ProblemDomain(a_domain, a_periodic);
#ifdef PR_MPI
    m_sendBufferSize = 0;
    m_recvBufferSize = 0;
    m_sendBuffer = nullptr;
    m_recvBuffer = nullptr;
    
    m_globalBufferSize = 0;
    m_localBufferSize = 0;
    //m_bufferSizes = (int*)malloc(numProc()*sizeof(int));
    //m_bufferOffsets = (int*)malloc(numProc()*sizeof(int));
    protoMalloc(MEMTYPE_DEFAULT, m_bufferSizes, numProc()*sizeof(int)); 
    protoMalloc(MEMTYPE_DEFAULT, m_bufferOffsets, numProc()*sizeof(int)); 
    m_localBuffer = nullptr;
    m_globalBuffer = nullptr;
#endif
}

FinitePointSet::FinitePointSet(ProblemDomain a_domain)
{
    m_domain = a_domain;
#ifdef PR_MPI
    m_sendBufferSize = 0;
    m_recvBufferSize = 0;
    m_sendBuffer = nullptr;
    m_recvBuffer = nullptr;
    
    m_globalBufferSize = 0;
    m_localBufferSize = 0;
    //m_bufferSizes = (int*)malloc(numProc()*sizeof(int));
    //m_bufferOffsets = (int*)malloc(numProc()*sizeof(int));
    protoMalloc(MEMTYPE_DEFAULT, m_bufferSizes, numProc()*sizeof(int)); 
    protoMalloc(MEMTYPE_DEFAULT, m_bufferOffsets, numProc()*sizeof(int)); 
    m_localBuffer = nullptr;
    m_globalBuffer = nullptr;
#endif
}

FinitePointSet::~FinitePointSet()
{
#ifdef PR_MPI
    if (m_sendBuffer != nullptr)
    {
        protoFree(MEMTYPE_DEFAULT, m_sendBuffer);
    }
    if (m_recvBuffer != nullptr)
    {
        protoFree(MEMTYPE_DEFAULT, m_recvBuffer);
    }
    protoFree(MEMTYPE_DEFAULT, m_bufferSizes);
    protoFree(MEMTYPE_DEFAULT, m_bufferOffsets);
    if (m_localBufferSize > 0)
    {
        protoFree(MEMTYPE_DEFAULT, m_localBuffer);
    }
    if (m_globalBufferSize > 0)
    {
        protoFree(MEMTYPE_DEFAULT, m_globalBuffer);
    }
#endif
}

bool FinitePointSet::add(Point a_point)
{
    if (m_points.count(a_point) > 0) {return false; }
    Point p = m_domain.image(a_point); 
    if (!m_domain.box().contains(p))
    {
        return false;
    }
#ifdef PR_MPI
    m_pointBuffer.push_back(a_point);
#else
    m_points.insert(a_point);
#endif
    return true;
}

void FinitePointSet::exchange()
{
#ifdef PR_MPI
    MPI_Barrier(MPI_COMM_WORLD);
    m_points.clear();

    // GET BUFFER SIZES
    int localBufferSize = DIM*m_pointBuffer.size();
    MPI_Allgather(&localBufferSize, 1, MPI_INT, m_bufferSizes, 1, MPI_INT, MPI_COMM_WORLD);  

    // ALLOCATE LOCAL BUFFER
    if (localBufferSize > m_localBufferSize)
    {
        if (m_localBufferSize > 0)
        {
            protoFree(MEMTYPE_DEFAULT, m_sendBuffer);
        }
        m_localBufferSize = localBufferSize;
        //m_localBuffer = (int*)malloc(m_localBufferSize*sizeof(int));
        protoMalloc(MEMTYPE_DEFAULT, m_localBuffer, m_localBufferSize*sizeof(int));
    }

    // WRITE TO LOCAL BUFFER
    int writeIndex = 0;
    for (auto iter = m_pointBuffer.begin(); iter != m_pointBuffer.end(); ++iter)
    {
        char* buf = ((char*)m_localBuffer) + writeIndex;
        (*iter).linearOut(buf);
        writeIndex += DIM*sizeof(int);
    }

    // ALLOCATE GLOBAL BUFFER AND OFFSETS
    int globalBufferSize = 0;
    for (int ii = 0; ii < numProc(); ii++)
    {
        m_bufferOffsets[ii] = globalBufferSize;
        globalBufferSize += m_bufferSizes[ii];
    }
    if (globalBufferSize > m_globalBufferSize)
    {
        if (m_globalBufferSize > 0)
        {
            protoFree(MEMTYPE_DEFAULT, m_globalBuffer);
        }
        m_globalBufferSize = globalBufferSize;
        //m_globalBuffer = (int*)malloc(m_globalBufferSize*sizeof(int));
        protoMalloc(MEMTYPE_DEFAULT, m_globalBuffer, m_globalBufferSize*sizeof(int));
    }

    // EXCHANGE DATA
    MPI_Allgatherv(m_localBuffer, localBufferSize, MPI_INT, 
            m_globalBuffer, m_bufferSizes, m_bufferOffsets, MPI_INT, MPI_COMM_WORLD);

    // READ DATA FROM BUFFER
    int readIndex = 0;
    for (int ii = 0; ii < globalBufferSize / DIM; ii++)
    {
        char* buf = ((char*)m_globalBuffer) + readIndex;
        Point p;
        p.linearIn(buf);
        m_points.insert(p);
        readIndex += DIM*sizeof(int);
    }
#endif
}

std::vector<Point>
FinitePointSet::points()
{
    exchange();
    std::vector<Point> points;
    for (auto iter = m_points.begin(); iter != m_points.end(); ++iter)
    {
        points.push_back(*iter);
    }
    return points;
}

void FinitePointSet::clear()
{
    m_points.clear();
    m_pointBuffer.clear();
}

#ifdef PR_MPI
/*
int
FinitePointSet::recvBufferSize()
{
    // assume the largest message is twice the perfectly load balanced value
    return (m_domain.size() / Proto::numProc())*2*DIM;
}
*/
#endif

