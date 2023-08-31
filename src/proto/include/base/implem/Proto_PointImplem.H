
CUDA_DECORATION
Point::Point()
{
    for (int j = 0; j < DIM; j++)
    {
        m_tuple[j]=0;
    }
}
CUDA_DECORATION
Point::Point(const int (&a_tuple)[DIM])
{
    for (int j = 0; j < DIM; j++)
    {
        m_tuple[j]=a_tuple[j];
    }
}

CUDA_DECORATION
Point::Point(const Point& a_pt)
{
    *this = a_pt;
}

CUDA_DECORATION
Point Point::Ones(int a_scale)
{
    int ones[DIM];
    for (int k = 0 ; k < DIM; k++)
    {
        ones[k] = a_scale;
    }
    Point p(ones);
    return p;
}

CUDA_DECORATION
Point Point::Zeros()
{
    int zeros[DIM];
    for (int k = 0 ; k < DIM; k++)
    {
        zeros[k] = 0;
    }
    return Point(zeros);
}

CUDA_DECORATION
Point Point::Basis(int idir, int a_scale)
{
    Point retval = Point::Zeros();
    retval[idir] = a_scale;
    return retval;
}

CUDA_DECORATION
Point Point::Basis(int idir, const Side::LoHiSide a_side)
{
    if (a_side == Side::Lo)
    {
        return Point::Basis(idir, -1);
    } else if (a_side == Side::Hi)
    {
        return Point::Basis(idir, 1);
    } else {
        return Point::Zeros();
    }
}

CUDA_DECORATION
const int& Point::operator[](unsigned char a_index) const
{
    PROTO_ASSERT((a_index < DIM),
            "Point::operator[](uchar a_index) invalid for a_index = %i.\
            a_index must be in [0,DIM = %i)",
            a_index, DIM);
    return m_tuple[a_index];
}

CUDA_DECORATION
int& Point::operator[](unsigned char a_index)
{
    PROTO_ASSERT(((a_index < DIM)),
            "Point::operator[](uchar a_index) invalid for a_index = %i.\
            a_index must be in [0,DIM = %i)",
            a_index, DIM);
    return m_tuple[a_index];
}



CUDA_DECORATION
Point& Point::operator=(const Point& a_rhs)
{
    for (int ii = 0; ii < DIM; ii++)
    {
        m_tuple[ii] = a_rhs[ii];
    }
    return *this;
}

CUDA_DECORATION
bool Point::operator<(const Point& a_rhs) const
{
    for (int ii = DIM-1; ii >= 0; ii--)
    {
        if (m_tuple[ii] > a_rhs[ii]){return false;}
        else if (m_tuple[ii] < a_rhs[ii]){return true;}
        //else they are equal, check next slowest index
    }
    //Points are equal; return false
    return false;
}

CUDA_DECORATION
bool Point::operator==(const Point& a_rhsPoint) const 
{
    bool ret = true;
    for (int j = 0; j < DIM; j++)
    {
        if (m_tuple[j] != a_rhsPoint[j]) ret=false;
    }
    return ret;
}

CUDA_DECORATION
bool Point::operator!=(const Point& a_rhsPoint) const 
{
    return !(*this == a_rhsPoint);
}

CUDA_DECORATION
Point Point::operator+(const Point& a_rhsPoint) const
{
    int tuple[DIM];
    for (int j = 0; j < DIM; j++)
    {
        tuple[j] = m_tuple[j]+a_rhsPoint.m_tuple[j]; //a_tuple
    }
    Point ret(tuple);
    return ret;
}

CUDA_DECORATION
Point Point::operator-(const Point& a_rhsPoint) const
{
    int tuple[DIM];
    for (int j = 0; j < DIM; j++)
    {
        tuple[j] = m_tuple[j]-a_rhsPoint.m_tuple[j]; 
    }
    Point ret(tuple);
    return ret;
}

CUDA_DECORATION
Point Point::operator*(const Point& a_rhsPoint) const
{
    int tuple[DIM];
    for (int j = 0; j < DIM; j++)
    {
        tuple[j] = m_tuple[j]*a_rhsPoint.m_tuple[j]; //a_tuple
    }
    Point ret(tuple);
    return ret;
} 
CUDA_DECORATION
Point Point::operator/(const Point& a_rhsPoint) const
{
    #ifndef NDEBUG
    for (int ii = 0; ii < DIM; ii++)
    {
        PROTO_ASSERT(a_rhsPoint[ii] != 0, 
        "Point::operator/(Point a_pt) invalid: division by zero (a_pt[%i] = 0)",ii);
    }
    // Division rounds towards -infinity, consistent with the geometric meaning coarsening.
    #endif
    int tuple[DIM];
    for (int j = 0; j < DIM; j++)
    {
        if (m_tuple[j]*a_rhsPoint[j] >= 0)
        {
            tuple[j] = m_tuple[j]/a_rhsPoint[j];
        }
        else if (-m_tuple[j]%a_rhsPoint[j] == 0)
          {
            tuple[j] = -(-m_tuple[j]/a_rhsPoint[j]);
          }
        else
          {
            tuple[j] = -(-m_tuple[j]/a_rhsPoint[j] + 1);
          }
    }
    Point ret(tuple);
    return ret;
}

CUDA_DECORATION
Point Point::operator%(const Point& a_pt) const
{
    int arr[DIM];
    for (int ii = 0; ii < DIM; ii++)
    {
      if (m_tuple[ii] < 0)
        {  
          arr[ii] = (m_tuple[ii]+(-m_tuple[ii]/a_pt[ii]+1)*a_pt[ii]) % a_pt[ii];
        }
      else
        {
          arr[ii] = m_tuple[ii]%a_pt[ii];
        }
    }
    return Point(arr);
}

Point Point::operator+(int a_nref) const
{
    return (*this)+(Point::Ones(a_nref));
}    

Point Point::operator-(int a_nref) const
{
    return (*this)-(Point::Ones(a_nref));
}    

Point Point::operator*(int a_nref) const
{
    int tuple[DIM];
    for (int j = 0; j < DIM; j++)
    {
        tuple[j] = m_tuple[j]*a_nref;
    }
    Point ret(tuple);
    return ret;
}
    
Point Point::operator/(int a_nref) const
{
    PROTO_ASSERT(a_nref != 0,
    "Point::operator/(int a_nref) invalid for a_nref = %i. Division by zero.",
    a_nref);
    return (*this)/(Point::Ones()*a_nref);
}

Point Point::operator%(int a_nref) const
{
    PROTO_ASSERT(a_nref != 0,
    "Point::operator%%(int a_nref) invalid for a_nref = %i. Division by zero.",
    a_nref);
    return (*this)%(Point::Ones()*a_nref);
}

CUDA_DECORATION
void Point::operator+=(const Point& a_rhsPoint)
{
    (*this) = (*this)+a_rhsPoint;
}

CUDA_DECORATION
void Point::operator-=(const Point& a_rhsPoint)
{
    (*this) = (*this)-a_rhsPoint;
}

CUDA_DECORATION
void Point::operator*=(const Point& a_rhsPoint)
{
    (*this) = (*this)*a_rhsPoint;
}

CUDA_DECORATION
void Point::operator/=(const Point& a_rhsPoint)
{
    (*this) = (*this)/a_rhsPoint;
}

CUDA_DECORATION
void Point::operator%=(const Point& a_pt)
{
    (*this) = (*this)%a_pt;
}

void Point::operator+=(int a_nref)
{
    (*this) = (*this)+a_nref;
}

void Point::operator-=(int a_nref)
{
    (*this) = (*this)-a_nref;
}
void Point::operator*=(int a_nref)
{
    (*this) = (*this)*a_nref;
}

void Point::operator/=(int a_nref)
{
    (*this) = (*this)/a_nref;
}

void Point::operator%=(int a_n)
{
    (*this) = (*this)%a_n;
}

int Point::dot(const Point& a_rhs) const
{
    return ((*this)*a_rhs).sum();
}

Point Point::coarsen(unsigned int a_refRatio) const
{
    PROTO_ASSERT(a_refRatio != 0, 
        "Error in Proto::Point::coarsen | refinement ratio must be non-zero.");
    Point retval = *this;
    for(int idir = 0; idir < DIM; idir++)
    {
        retval[idir] /= a_refRatio; 
    }
    return retval;
}

CUDA_DECORATION
inline int Point::sum() const
{
    int s = 0;
    for (int ii = 0; ii < DIM; ii++) {
        s += m_tuple[ii];
    }
    return s;
}

int Point::min() const
{
    int val = m_tuple[0];
    for (int dir = 1; dir < DIM; dir++)
    {
        val = std::min(val, m_tuple[dir]);
    }
    return val;
}

int Point::max() const
{
    int val = m_tuple[0];
    for (int dir = 1; dir < DIM; dir++)
    {
        val = std::max(val, m_tuple[dir]);
    }
    return val;
}

int Point::sumAbs() const
{
    int s = 0;
    for (int ii = 0; ii < DIM; ii++) {
        s += std::abs(m_tuple[ii]);
    }
    return s;
}

bool Point::isotropic() const
{
    for (int dir = 1; dir < DIM; dir++)
    {
        if (m_tuple[dir] != m_tuple[0]) {return false;}
    }
    return true;
}

size_t Point::linearSize() const
{
    size_t retval = DIM*sizeof(int);
    return retval;
}

void Point::linearIn(const char* a_buf, size_t expected_size)
{
    const int* intbuf = (int *) a_buf;
    for(int idir = 0; idir < DIM; idir++)
    {
        m_tuple[idir] = intbuf[idir];
    }
}


void Point::linearOut(char* a_buf) const
{
    int* intbuf = (int *) a_buf;
    for(int idir = 0; idir < DIM; idir++)
    {
        intbuf[idir] = m_tuple[idir];
    }
}

void Point::print() const
{
    using std::cout;
    using std::endl;
    cout << "(";
    for (int dir = 0; dir < DIM ; dir++)
    {
        cout << m_tuple[dir];
        if(dir < DIM-1) cout << ",";
    }
    cout <<   ")"  << endl;
}

std::ostream& operator <<(std::ostream& a_os, const Point& a_pt)
{
    a_os << "(" ;
    for (int dir = 0; dir < DIM; dir++) 
    {
        a_os << a_pt[dir];
        if(a_pt[dir] < 10)
        {
          a_os << " ";
        }
        if (dir < DIM-1) 
        {
            a_os << ",";
        }
    }
    a_os << ")";
    return a_os;
}

std::istream& operator >>(std::istream& a_is,  Point& a_pt)
{
  char dum;
  a_is >> dum;
  for (int dir = 0; dir < DIM; dir++) 
  {
    a_is >> a_pt[dir];
    if (dir < DIM-1) 
    {
      a_is >> dum;
    }
  }
  a_is >> dum;
  return a_is;
}
