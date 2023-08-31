CInterval::CInterval(
        unsigned int a_c0,  unsigned int a_c1,
        unsigned int a_d0,  unsigned int a_d1,
        unsigned int a_e0,  unsigned int a_e1)
{
    PROTO_ASSERT(a_c1 >= a_c0,
            "CInterval Constructor | Error: \
            CInterval(...) invalid for a_c0 = %i and a_c1 = %i. \
            lower bound must be less than or equal to high bound.",
            a_c0, a_c1);
    PROTO_ASSERT(a_d1 >= a_d0,
            "CInterval Constructor | Error: \
            CInterval(...) invalid for a_d0 = %i and a_d1 = %i. \
            lower bound must be less than or equal to high bound.",
            a_d0, a_d1);
    PROTO_ASSERT(a_e1 >= a_e0,
            "CInterval Constructor | Error: \
            CInterval(...) invalid for a_e0 = %i and a_e1 = %i. \
            lower bound must be less than or equal to high bound.",
            a_e0, a_e1);
    m_comps[0][0] = a_c0; m_comps[0][1] = a_c1;
    m_comps[1][0] = a_d0; m_comps[1][1] = a_d1;
    m_comps[2][0] = a_e0; m_comps[2][1] = a_e1;
}

CInterval::CInterval(
        std::initializer_list<std::initializer_list<unsigned int>> a_lst)
{
    int comp = 0;
    for (auto& c : a_lst)
    {
        if (c.size() == 0)
        {
            m_comps[comp][0] = 0;
            m_comps[comp][1] = 0;
            comp++;
            continue;
        }
        int bound = 0;
        for (auto& b : c)
        {
            m_comps[comp][bound] = b;
            bound++;
        }
        comp++;
    }
    PROTO_ASSERT(m_comps[0][1] >= m_comps[0][0],
            "CInterval Constructor | Error: \
            CInterval(...) invalid for a_c0 = %i and a_c1 = %i. \
            lower bound must be less than or equal to high bound.",
            m_comps[0][0], m_comps[0][1]);
    PROTO_ASSERT(m_comps[1][1] >= m_comps[1][0],
            "CInterval Constructor | Error: \
            CInterval(...) invalid for a_d0 = %i and a_d1 = %i. \
            lower bound must be less than or equal to high bound.",
            m_comps[1][0], m_comps[1][1]);
    PROTO_ASSERT(m_comps[2][1] >= m_comps[2][0],
            "CInterval Constructor | Error: \
            CInterval(...) invalid for a_e0 = %i and a_e1 = %i. \
            lower bound must be less than or equal to high bound.",
            m_comps[2][0], m_comps[2][1]);
}

bool CInterval::operator==(CInterval& a_rhs)
{
    bool equal = true;
    equal &= (m_comps[0][0] == a_rhs.m_comps[0][0]);
    equal &= (m_comps[1][0] == a_rhs.m_comps[1][0]);
    equal &= (m_comps[2][0] == a_rhs.m_comps[2][0]);
    equal &= (m_comps[0][1] == a_rhs.m_comps[0][1]);
    equal &= (m_comps[1][1] == a_rhs.m_comps[1][1]);
    equal &= (m_comps[2][1] == a_rhs.m_comps[2][1]);
    return equal;
}

unsigned int CInterval::low(unsigned int a_comp) const
{
    PROTO_ASSERT((a_comp < 3),
        "CInterval::low | Error: \
        invalid for a_comp = %i. a_comp must be in [0,3).", a_comp);
    return m_comps[a_comp][0];
}

unsigned int CInterval::high(unsigned int a_comp) const
{
    PROTO_ASSERT((a_comp < 3),
        "CInterval::high | Error: \
         Invalid for a_comp = %i. a_comp must be in [0,3).", a_comp);
    return m_comps[a_comp][1];
}

bool CInterval::contains( unsigned int a_index, unsigned int a_comp) const
{
    PROTO_ASSERT((a_comp < 3),
        "CInterval::contains | Error: \
        invalid for a_comp = %i. a_comp must be in [0,3).", a_comp);
    return ((a_index >= m_comps[a_comp][0])
            && (a_index <= m_comps[a_comp][1]));
}

unsigned int CInterval::size(unsigned int a_comp) const
{
    PROTO_ASSERT((a_comp < 3),
        "CInterval::size | Error: \
        invalid for a_comp = %i. a_comp must be in [0,3).", a_comp);
    return (m_comps[a_comp][1] - m_comps[a_comp][0] + 1);
}

void CInterval::print() const
{
    std::cout << "{ ";
    std::cout << "{" << m_comps[0][0] << ", " << m_comps[0][1] << "}";
    std::cout << ", {" << m_comps[1][0] << ", " << m_comps[1][1] << "}";
    std::cout << ", {" << m_comps[2][0] << ", " << m_comps[2][1] << "}";
    std::cout << " }";
}

std::ostream& operator<<(std::ostream& a_os, const CInterval& a_int)
{
    a_os << "{ ";
    a_os << "{" << a_int.low(0) << ", " << a_int.high(0) << "}";
    a_os << ", {" << a_int.low(1) << ", " << a_int.high(1) << "}";
    a_os << ", {" << a_int.low(2) << ", " << a_int.high(2) << "}";
    a_os << " }";
    return a_os;
}
