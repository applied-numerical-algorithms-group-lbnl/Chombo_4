//=================================================================================================
// Constructors ||
//==============++ 

Box::Box()
{
    m_low = Point::Zeros();
    m_high= Point::Ones()*(-1);
    m_size = 0;
}

CUDA_DECORATION
Box::Box(const Box& a_box)
{
    m_low = a_box.m_low;
    m_high = a_box.m_high;
    m_size = a_box.m_size;
}

Box::Box(const Point& a_lowCorner,const Point& a_highCorner)
{
    m_low = a_lowCorner;
    m_high = a_highCorner;
    recomputeSize();
}

Box::Box(const Point& a_sizes)
{
    m_low = Point::Zeros();
    m_high = a_sizes - Point::Ones();
    recomputeSize();
}

void
Box::recomputeSize()
{
    m_size = 1;
    for(int idir = 0; idir < DIM; idir++)
    {
        m_size *= size(idir);
    }
}

void
Box::define(const int* a_lo, const int* a_hi)
{
    Point lo, hi;
    for(int idir = 0; idir < DIM; idir++)
    {
        lo[idir] = a_lo[idir];
        hi[idir] = a_hi[idir];
    }
    *this = Box(lo, hi);
}


//=================================================================================================
// Static Methods ||
//================++ 

Box Box::Cube(int a_size)
{
    return Box(Point::Zeros(),Point::Ones()*(a_size-1));
}

Box Box::Kernel(int a_radius)
{
    return Box::Cube(2*a_radius+1).shift(Point::Ones(-a_radius));
}

CUDA_DECORATION
inline unsigned int Box::sindex(Point a_pt, Point a_low, Point a_high)
{
    unsigned int factor = 1;
    unsigned int linIndex = a_pt[0] - a_low[0];
    for (unsigned char i = 1;i < DIM;i++)
    { 
        factor = factor*(a_high[i-1] - a_low[i-1] + 1);
        linIndex = linIndex + (a_pt[i] - a_low[i])*factor;
    }
    return linIndex;
}

/*
//TODO: Doesn't appear to be used anywhere -CLG
CUDA_DECORATION
inline unsigned int Box::soffset(const Point& a_pt, const Point& a_low, const Point& a_high)
{
    int factor = 1;
    int linIndex = a_pt[0];
    for(unsigned char i=1; i<DIM; i++)
    {
        factor = factor*(a_high[i-1] - a_low[i-1] + 1);
        linIndex = linIndex + (a_pt[i])*factor;
    }
    return linIndex;
}     
*/

//=================================================================================================
// Accessors and Queries ||
//=======================++ 

CUDA_DECORATION
inline std::size_t Box::size(unsigned char a_dim) const
{
    PROTO_ASSERT((a_dim < DIM),
        "Box::size(int a_dim) with a_dim = %i is invalid. a_dim should be in [0,DIM=%i).\n",
        a_dim, DIM);
    if (m_high[a_dim] < m_low[a_dim]){return 0;}
    else {
        return (m_high[a_dim] - m_low[a_dim] + 1);
    }
}

CUDA_DECORATION
inline std::size_t Box::size() const
{
    if (m_size <= 0)
    {
        return 0;
    }
    else 
    {
        return m_size;
    }
}

CUDA_DECORATION
inline Point Box::sizes() const
{
    return high() - low() + Point::Ones();
}

CUDA_DECORATION
bool Box::contains(const Point& a_pt) const
{
    if (this->empty()){return false;}
    for(int idir = 0; idir < DIM; idir++)
    {
        if(a_pt[idir] < m_low[idir])
        {
            return false;
        }
        if(a_pt[idir] > m_high[idir])
        {
            return false;
        }
    }
    return true;
}
  
CUDA_DECORATION
bool Box::contains(const Box& a_rhs) const
{
    return (this->contains(a_rhs.m_low) && this->contains(a_rhs.m_high));
}

bool Box::onBoundary(const Point& a_p) const
{
    if (!(this->contains(a_p))){return false;}
    for (int ii = 0; ii < DIM; ii++)
    {
        if ((a_p[ii] == m_low[ii]) || (a_p[ii] == m_high[ii]))
        {
            return true;
        }
    }
    return false;
}

CUDA_DECORATION
inline unsigned int Box::index(const Point& a_pt) const
{
    PROTO_ASSERT(contains(a_pt),
        "Box::index(Point a_pt) argument invalid. a_pt should be contained in *this Box.\n");
    return Box::sindex(a_pt, m_low, m_high);
}


CUDA_DECORATION
inline Point Box::operator[](unsigned int a_index) const
{
    PROTO_ASSERT(!empty(),
            "Box::operator[] | Error: Box is empty.");
    PROTO_ASSERT(a_index < m_size,
            "Box::operator[] | Error: argument %i invalid. a_index should be in [0,%i).",
            a_index, m_size);
    
    int tuple[DIM];
    for (unsigned char i = 0;i < DIM; i++)
    {      
        int factor = (m_high[i] - m_low[i] + 1);
        int kred = a_index % factor;
        tuple[i] = kred + m_low[i];
        a_index = (a_index - kred)/factor;
    }
    return Point(tuple);
}

CUDA_DECORATION
inline Point Box::operator()(unsigned int a_idx, unsigned int a_idy, unsigned int a_idz) const
{
    
    int tuple[DIM];
    tuple[0] = m_low[0] + a_idx;
    tuple[1] = m_low[1] + a_idy;
#if DIM == 3
    tuple[2] = m_low[2] + a_idz;
#endif 
    return Point(tuple);
}

//=================================================================================================
// Operators ||
//===========++ 

Box Box::operator&(const Box& a_rightBox) const
{
    int newLo[DIM];
    int newHi[DIM];
    for (int i = 0 ; i < DIM; i++)
    {
        newLo[i] = m_low[i];
        if (m_low[i] < a_rightBox.m_low[i])
        {
            newLo[i] = a_rightBox.m_low[i];
        }
        newHi[i] = m_high[i];
        if (m_high[i] > a_rightBox.m_high[i])
        {
            newHi[i] = a_rightBox.m_high[i];
        }
    }
    Point newLowCorner(newLo);
    Point newHighCorner(newHi); 
    for (int dir = 0; dir < DIM; dir++)
    {
        if (newHi[dir] < newLo[dir])
        {
            //if any of the dimensions are inside-out, return an empty (default) Box
            Box ret0;
            return ret0;
        }
    }
    Box ret(newLowCorner,newHighCorner);
    return ret;
}

void Box::operator&=(const Box& a_rhs)
{
    Box intersect = a_rhs & (*this);
    *this = intersect;
}

Box Box::operator&(const Point& a_pt) const
{
    if (contains(a_pt)){ return *this; }
    else if(empty()) {return Box(a_pt, a_pt); }
    else {
        Point L = m_low;
        Point H = m_high;
        for (int ii = 0; ii < DIM; ii++)
        {
            L[ii] = std::min(a_pt[ii], L[ii]);
            H[ii] = std::max(a_pt[ii], H[ii]);
        }
        return Box(L,H);
    }
}

void Box::operator&=(const Point& a_pt)
{
    if (this->contains(a_pt)){return;}
    Box intersect = (*this) & a_pt;
    *this = intersect;
}

bool Box::operator==(const Box& a_rhsBox) const
{
    return ((m_low == a_rhsBox.m_low) 
       && (m_high == a_rhsBox.m_high));
}

CUDA_DECORATION
Point Box::operator%(const Point& a_pt) const
{
    return mod(a_pt);
}

CUDA_DECORATION
inline Point Box::mod(const Point& a_pt) const
{
    PROTO_ASSERT(m_size > 0,
    "Box::mod(Point a_pt) invalid for empty Box.");
    Point p = a_pt - m_low;
    Point q = p % (m_high - m_low + Point::Ones());
    for (int ii = 0; ii < DIM; ii++)
    {
        if (q[ii] < 0)
        {
            q[ii] += this->size(ii);
        }
    }
    return q + m_low;
}

//=================================================================================================
// Transformations ||
//=================++ 

Box Box::shift(int a_direction, int a_offset) const
{
    PROTO_ASSERT(((a_direction >= 0) && (a_direction < DIM)),
        "Box::shift(int a_dir, int a_dist) with a_dir = %i is invalid. a_dir should be in [0,DIM=%i).\n",
        a_direction, DIM);
    Box returnBox = Box(*this);
    returnBox.m_low += Point::Basis(a_direction)*a_offset;
    returnBox.m_high += Point::Basis(a_direction)*a_offset;
    return returnBox;
}

Box Box::shift(const Point& a_pt) const
{
    Box returnBox = Box(*this);
    returnBox.m_low += a_pt;
    returnBox.m_high += a_pt;
    return returnBox;
}

Box Box::grow(int a_offset) const
{
    return grow(Point::Ones()*a_offset);
}

Box Box::grow(const Point& a_offset) const
{
    Point lo = m_low;
    Point hi = m_high;
    lo -= a_offset;
    hi += a_offset;
    Box returnBox(lo, hi);

    return returnBox;
}

Box Box::grow(int a_dir, int a_offset) const
{
    return grow(Point::Basis(a_dir, a_offset));
}

Box Box::grow(int a_dir, Side::LoHiSide a_side, int a_dist) const
{
    Point p = Point::Basis(a_dir, a_side);
    return extrude(p, a_dist);
}

Box Box::extrude(const Point& a_dir, int a_dist) const
{
    Point lo = m_low;
    Point hi = m_high;

    for (int ii = 0; ii < DIM; ii++)
    {
        if (a_dir[ii] > 0)
        {
            hi[ii] += a_dist;
        } else if (a_dir[ii] < 0)
        {
            lo[ii] -= a_dist;
        }
    }
    return Box(lo, hi);
}

Box Box::extrude(int a_dir, int a_dist, bool a_upper) const
{
  PROTO_ASSERT(((a_dir >= 0) && (a_dir < DIM)),
    "error in Box::extrude | a_dir = %i is invalid. a_dir should be in [0,DIM=%i).",
    a_dir, DIM);
    if (a_upper)
    {
        return extrude(Point::Basis(a_dir), a_dist);
    } else {
        return extrude(-Point::Basis(a_dir), a_dist);
    }
}

Box Box::growHi(int a_dir, int a_offset) const
{
    return grow(a_dir, Side::Hi, a_offset);
}

Box Box::growLo(int a_dir, int a_offset) const
{
    return grow(a_dir, Side::Lo, a_offset);
}

Point Box::boundary(const Side::LoHiSide a_side) const
{
    Point retval = low();
    if(a_side == Side::Hi)
    {
        retval =  high();
    }
    return retval;
}

/* TODO: Redundant, and doesn't appear to be used anywhere
Box Box::faceBox(int a_dir, const Side::LoHiSide a_side) const
{
    Box retval = adjCellSide(a_idir, 1, a_side);
    if(a_side == Side::Lo)
    {
        retval =  retval.shift(Point::Basis(a_idir, 1));
    }
    return retval;
}
*/

Box Box::coarsen(unsigned int a_nref) const
{
    PROTO_ASSERT(a_nref > 0,
    "Box::coarsen(uint a_nref) with a_nref = %i is invalid. a_nref should be positive.", a_nref);

    return coarsen(Point::Ones()*a_nref);
}

Box Box::coarsen(const Point& a_pt) const
{
    #ifndef NDEBUG
    for (int ii = 0; ii < DIM; ii++)
    {
        PROTO_ASSERT(a_pt[ii] > 0,
    "Box::coarsen(Point a_pt) with a_pt[%i]=%i is invalid. All components of a_pt must be positive",
        ii,a_pt[ii]);
    }
    #endif
    Point lo = m_low;
    Point hi = m_high;
    for (int d=0; d< DIM; d++)
      {
        if (lo[d] < 0)
          {
            if ((-lo[d])%a_pt[d] == 0)
              {
                lo[d] = -(-lo[d]/a_pt[d]);
              }
            else
              {
                lo[d] = -(-lo[d]/a_pt[d] + 1);
              }
          }
        else
          {
            lo[d] = lo[d]/a_pt[d];
          }
        if (hi[d] < 0)
          {
            if (hi[d]%a_pt[d] == 0)
              {
                hi[d] = -(-hi[d]/a_pt[d]);
              }
            else
              {
                hi[d] = -(-hi[d]/a_pt[d] + 1);
              }
          }
        else
          {
            hi[d] = hi[d]/a_pt[d];
          }    
      } 
    Box returnBox(lo, hi);
    return returnBox;
}

Box Box::taperCoarsen(const Point& a_ref) const
{
    Point H = m_high/a_ref;
    Point L;
    for (int ii = 0; ii < DIM; ii++)
    {
        L[ii] = ceil(m_low[ii]/((double)a_ref[ii]));
    }
    return Box(L,H);
}

bool Box::coarsenable(const int& a_ratio) const
{
    Point ratio = Point::Ones(a_ratio);
    return coarsenable(ratio);
}

bool Box::coarsenable(const Point& a_ratio) const
{
    Box coarBox  =  this->coarsen(a_ratio);
    Box refCoBox = coarBox.refine(a_ratio);
    return(refCoBox == (*this));
}

Box Box::refine(unsigned int a_nref) const
{
    PROTO_ASSERT(a_nref > 0,
    "Box::refine(uint a_nref) with a_nref = %i is invalid. a_nref should be positive.", a_nref);

    return refine(Point::Ones()*a_nref);
}

Box Box::refine(const Point& a_pt) const
{
    #ifndef NDEBUG
    for (int ii = 0; ii < DIM; ii++)
    {
        PROTO_ASSERT(a_pt[ii] > 0,
            "Box::refine(Point a_pt) with a_pt[%i]=%i is invalid. All components of a_pt must be positive.",
        ii, a_pt[ii]);
    }
    #endif
    Point lo = m_low;
    Point hi = m_high;
    lo *= a_pt;
    hi += Point::Ones();
    hi *= a_pt;
    hi -= Point::Ones();
    Box returnBox(lo, hi);

    return returnBox;
}

Box Box::edge(const Point& a_dir, int a_dist) const
{
    int d = a_dist - 1;
    Point H = m_high; Point L = m_low;
    for (int ii = 0; ii < DIM; ii++)
    {
        if (a_dir[ii] > 0){L[ii] = H[ii] - d;}
        if (a_dir[ii] < 0){H[ii] = L[ii] + d;}
    }
    return Box(L,H);
}

Box Box::face(int a_dir, Side::LoHiSide a_side, int a_dist) const
{
    PROTO_ASSERT(((a_dir >= 0) && (a_dir < DIM)),
        "Box::face(int a_dir) with a_dir = %i is invalid. a_dir should be in [0,DIM=%i).",
        a_dir, DIM);
    int s = size(a_dir);
    int dist = std::min(a_dist, s);
    if (a_side == Side::Lo)
    {
        return edge(Point::Basis(a_dir, -1), dist);
    } else if (a_side == Side::Hi)
    {
        return edge(Point::Basis(a_dir, +1), dist);
    } else {
        return Box();
    }
}

Box Box::flatten(const int a_dir, bool a_upper) const
{
    if (a_upper)
    {
        return face(a_dir, Side::Hi, 1);
    } else {
        return face(a_dir, Side::Lo, 1);
    }
}

Box Box::adjacent(const Point& a_dir, int a_dist) const
{
    Point L = low();
    Point H = high();
    if (a_dist > 0)
    {
        for (int ii = 0; ii < DIM; ii++)
        {
            if (a_dir[ii] > 0)
            {
                L[ii] = H[ii] + 1;
                H[ii] = H[ii] + a_dist;
            } else if (a_dir[ii] < 0)
            {
                H[ii] = L[ii] - 1;
                L[ii] = L[ii] - a_dist;
            }
        }
    } else {
        for (int ii = 0; ii < DIM; ii++)
        {
            if (a_dir[ii] > 0)
            {
                L[ii] = H[ii] + 1;
                H[ii] = H[ii] + size(ii);
            } else if (a_dir[ii] < 0)
            {
                H[ii] = L[ii] - 1;
                L[ii] = L[ii] - size(ii);
            }
        }
    }
    return Box(L,H);
}

Box Box::adjacent(int a_dir, Side::LoHiSide a_side, int a_dist) const
{
    PROTO_ASSERT(((a_dir >= 0) && (a_dir < DIM)),
        "Box::flatten(int a_dir) with a_dir = %i is invalid. a_dir should be in [0,DIM=%i).",
        a_dir, DIM);
    if (a_side == Side::Lo)
    {
        return adjacent(Point::Basis(a_dir, -1), a_dist);
    } else if (a_side == Side::Hi)
    {
        return adjacent(Point::Basis(a_dir, +1), a_dist);
    }
    else return Box();
}

Box Box::toOrigin() const
{
    return this->shift(-m_low);
}

Box Box::adjCellLo(int a_idir, int a_length) const
{
    return adjacent(Point::Basis(a_idir, -1), a_length);
}


Box Box::adjCellHi(int a_idir, int a_length) const
{
    return adjacent(Point::Basis(a_idir, +1), a_length);
}


//=================================================================================================
// Utility ||
//=========++ 

BoxIterator Box::begin() const
{
    return BoxIterator(*this,0);
}

BoxIterator Box::end() const
{
    return BoxIterator(*this, this->size());
}

BoxIterator Box::rbegin() const
{
    return BoxIterator(*this,this->size()-1);
}

BoxIterator Box::rend() const
{
    return BoxIterator(*this, -1);
}

void Box::print() const 
{
    std::cout<<*this<<std::endl;
}

size_t Box::linearSize() const
{
    size_t retval = 0;
    retval += m_low.linearSize();
    retval += m_high.linearSize();
    retval += sizeof(int);
    return retval;
}

void Box::linearIn(const char* a_buf)
{
    size_t charsiz = 0;
    const char* charbuf = a_buf;

    m_low.linearIn(charbuf);
    charbuf += m_low.linearSize();
    charsiz += m_low.linearSize();

    m_high.linearIn(charbuf);
    charbuf += m_high.linearSize();
    charsiz += m_high.linearSize();

    int* intbuf = (int* ) charbuf;
    m_size = *intbuf;
    charsiz += sizeof(int);
    charbuf += sizeof(int);
}

void Box::linearOut(char* a_buf) const
{
    size_t charsiz = 0;
    char* charbuf = a_buf;

    m_low.linearOut(charbuf);
    charbuf += m_low.linearSize();
    charsiz += m_low.linearSize();

    m_high.linearOut(charbuf);
    charbuf += m_high.linearSize();
    charsiz += m_high.linearSize();

    int* intbuf = (int* ) charbuf;
    *intbuf = m_size;
    charsiz += sizeof(int);
    charbuf += sizeof(int);
}

//=================================================================================================
// BoxIterator ||
//=============++ 

BoxIterator::BoxIterator(const Box& a_box, int a_pos)
{
    m_box = a_box;
    m_pos = a_pos;
}

BoxIterator::BoxIterator(const Box& a_box, const Point& a_pt)
{
    m_box = a_box;
    m_pos = a_box.index(a_pt);
}

bool BoxIterator::ok() const
{
    return (m_pos >= 0 && m_pos < m_box.size()); 
}

bool BoxIterator::operator==(const BoxIterator& a_iter) const
{
    return ((m_pos == a_iter.m_pos) && (m_box == a_iter.m_box));
}

BoxIterator& BoxIterator::operator++()
{
    m_pos++;
    return *this;
}

BoxIterator& BoxIterator::operator--()
{
    m_pos--;
    return *this;
}

BoxIterator BoxIterator::operator++(int)
{
    BoxIterator ptr(*this);
    m_pos++;
    return ptr;
}

BoxIterator BoxIterator::operator--(int)
{
    BoxIterator ptr(*this);
    m_pos--;
    return ptr;
}
