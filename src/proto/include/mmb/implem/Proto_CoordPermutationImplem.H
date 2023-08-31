CoordPermutation::CoordPermutation()
{
    for (int ii = 0; ii < DIM; ii++)
    {
        m_matrix[ii].fill(0);
        m_matrix[ii][ii] = 1;
    }
}

CoordPermutation::CoordPermutation(std::array<std::array<int, DIM>, DIM> a_matrix)
{
    m_matrix = a_matrix;
}

CoordPermutation::CoordPermutation(
        std::vector<std::array<int, 3>>& a_mapping)
{
    define(a_mapping);
}

CoordPermutation::CoordPermutation(
        std::initializer_list<std::array<int, 3>> a_mapping)
{
    std::vector<std::array<int, 3>> v(a_mapping);
    define(v);
}

void CoordPermutation::define(std::vector<std::array<int, 3>>& a_mapping)
{
    for (auto map : a_mapping)
    {
        PROTO_ASSERT(map[0] >= 0 && map[0] < DIM,
                "CoordPermutation::Constructor | Invalid mapping input: [%i, %i, %i].",
                map[0], map[1], map[2]);
        PROTO_ASSERT(map[1] >= 0 && map[1] < DIM,
                "CoordPermutation::Constructor | Invalid mapping input: [%i, %i, %i].",
                map[0], map[1], map[2]);
        PROTO_ASSERT(map[2] == 1 || map[2] == -1,
                "CoordPermutation::Constructor | Invalid mapping input: [%i, %i, %i].",
                map[0], map[1], map[2]);
    }
    for (int ii = 0; ii < DIM; ii++)
    {
        m_matrix[ii].fill(0);
        m_matrix[ii][ii] = 1;
    }
    for (auto map : a_mapping)
    {
        m_matrix[map[1]].fill(0);
        m_matrix[map[1]][map[0]] = map[2];
    }
    PROTO_ASSERT(this->isValid(),
            "CoordPermutation::Constructor | Resulting object is not a valid permutation.");
}

CoordPermutation CoordPermutation::identity()
{
    return CoordPermutation();
}

CoordPermutation CoordPermutation::reverse()
{
    CoordPermutation ret;
    for (int ii = 0; ii < DIM; ii++)
    {
        ret.m_matrix[ii][ii] = -1;
    }
    return ret;
}

CoordPermutation CoordPermutation::cw()
{
#if DIM >= 2
    return CoordPermutation{{0,1,1},{1,0,-1}};
#else
    return CoordPerutation();
#endif
}

CoordPermutation CoordPermutation::rotatePlane(unsigned int a_c0, unsigned int a_c1)
{
#if DIM >= 2
    return CoordPermutation({{a_c0,a_c1,1},{a_c1,a_c0,-1}});
#else
    return CoordPerutation();
#endif
}

CoordPermutation CoordPermutation::ccw()
{
#if DIM >= 2
    return CoordPermutation{{0,1,-1},{1,0,1}};
#else
    return CoordPerutation();
#endif
}

Point CoordPermutation::operator()(Point a_pt) const
{
    Point mapped;
    for (int jj = 0; jj < DIM; jj++)
    {
        for (int ii = 0; ii < DIM; ii++)
        {
            if (m_matrix[jj][ii] != 0)
            {
                mapped[jj] = m_matrix[jj][ii]*a_pt[ii];
            }
        }
    }
    return mapped;
}

CoordPermutation CoordPermutation::inverse() const
{
    std::array<std::array<int, DIM>, DIM> invMatrix;
    for (int jj = 0; jj < DIM; jj++)
    {
        for (int ii = 0; ii < DIM; ii++)
        {
            invMatrix[ii][jj] = m_matrix[jj][ii];
        }
    }
    return CoordPermutation(invMatrix);
}

Point CoordPermutation::rotatePoint(Point a_point, Box a_box) const
{
    PROTO_ASSERT(a_box.low() == Point::Zeros(),
            "CoordPermutation::rotatePoint | Error: input Box must have an origin at (0,...,0)");
    Point a = -(a_box.sizes());
    Point b = (*this)(a);
    Point s;
    for (int d = 0; d < DIM; d++)
    {
        s[d] = max(0,b[d]);
    }
    Point r = (*this)(a_point) + s;
    //std::cout << "Rotated " << a_point << " in box " << a_box << " and got " << r << std::endl;
    return r;
}

Point CoordPermutation::rotateCell(Point a_point, Box a_box) const
{
    PROTO_ASSERT(a_box.low() == Point::Zeros(),
            "CoordPermutation::rotatePoint | Error: input Box must have an origin at (0,...,0)");
    Point a = -(a_box.sizes() - Point::Ones());
    Point b = (*this)(a);
    Point s;
    for (int d = 0; d < DIM; d++)
    {
        s[d] = max(0,b[d]);
    }
    Point r = (*this)(a_point) + s;
    //std::cout << "Rotated " << a_point << " in box " << a_box << " and got " << r << std::endl;
    return r;
}

bool CoordPermutation::operator==(const CoordPermutation& a_rhs) const
{
    for (int jj = 0; jj < DIM; jj++)
    {
        for (int ii = 0; ii < DIM; ii++)
        {
            if (m_matrix[ii][jj] != a_rhs.m_matrix[ii][jj]) { return false; }
        }
    }
    return true;
}

void CoordPermutation::print() const
{
    for (int jj = 0; jj < DIM; jj++)
    {
        for (int ii = 0; ii < DIM; ii++)
        {
            if (m_matrix[jj][ii] != 0)
            {
                pout() << ii << " --> " << jj;
                pout() << " (" << m_matrix[jj][ii] << ")" << std::endl;
            }
        }
    }
    this->printMatrix();
}

void CoordPermutation::printMatrix() const
{
    for (auto row : m_matrix)
    {
        for (auto val : row)
        {
            pout() << val << ", ";
        }
        pout() << std::endl;
    }
}

bool CoordPermutation::isValid() const
{
    for (auto row : m_matrix)
    {
        int nonZeros = 0;
        for (auto val : row)
        {
            if (val != 0)
            {
                nonZeros++;
                if (val < -1 || val > 1)
                {
                    return false;
                }
            }
        }
        if (nonZeros != 1) { return false; }
    }
    return true;
}

CoordPermutation operator*(const CoordPermutation& a_A, const CoordPermutation& a_B)
{
    auto A = a_A.matrix();
    auto B = a_B.matrix();
    std::array<std::array<int, DIM>, DIM> M;
    for (int ii = 0; ii < DIM; ii++)
    {
        for (int jj = 0; jj < DIM; jj++)
        {
            int Mij = 0;
            for (int kk = 0; kk < DIM; kk++)
            {
                Mij += A[kk][ii]*B[jj][kk];
            }
            M[jj][ii] = Mij;
        }
    }
    return CoordPermutation(M);
}

