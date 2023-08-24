template<class T>
void InterpStencil<T>::define(Point a_ratio)
{
    for (int ii = 0; ii < DIM; ii++)
    {
        PROTO_ASSERT(a_ratio[ii] > 0,
                "InterpStencil(Point ratio) invalid. All ratios must be 1 or greater");
    }
    m_r = a_ratio;
    m_k = Box(a_ratio);
    m_s.resize(m_k.size());
    m_closed = false;
}

template<class T>
void InterpStencil<T>::close()
{
    if (m_closed == true) { return; }
    for (auto iter = m_k.begin(); iter != m_k.end(); ++iter)
    {
        Stencil<T>& S = (*this)(*iter);
        S.destRatio() = m_k.high() + Point::Ones();
        S.destShift() = *iter;
    }
    m_closed = true;
}

template<class T>
Stencil<T>& InterpStencil<T>::operator()(Point a_p)
{
    PROTO_ASSERT(!m_closed,
        "InterpStencil::operator() | Error: \
        Components of InterpStencil are read-only once closed. Use InterpStencil::get instead.");
    PROTO_ASSERT(m_k.contains(a_p),
        "InterpStencil::operator() | Error: \
        p is not in the InterStencil kernel");
    return m_s[m_k.index(a_p)];
}

template<class T>
const Stencil<T>& InterpStencil<T>::get(Point a_p) const
{
    PROTO_ASSERT(m_k.contains(a_p),
        "InterpStencil::operator() | Error: \
        p is not a member of the InterStencil kernel");
    return m_s[m_k.index(a_p)];
}

template<class T>
Box InterpStencil<T>::span() const
{
    Box span;
    for (int ii = 0; ii < m_s.size(); ii++)
    {
        span = span & m_s[ii].span().low();
        span = span & m_s[ii].span().high();
    }
    return span;
}

template<class T>
Point InterpStencil<T>::ghost() const
{
    Box spanBox = span();
    int v[DIM];
    for (int ii = 0; ii < DIM; ii++)
    {
        v[ii] = std::max(std::abs(spanBox.low()[ii]), std::abs(spanBox.high()[ii]));
    }
    return Point(v);
}

template<class T>
void InterpStencil<T>::print() const
{
    for (auto iter = begin(); iter != end(); ++iter)
    {
        (*iter).print();
    }
}

template<class T>
template<unsigned int C, MemType MEMTYPE, unsigned char D, unsigned char E>
LazyStencil<T,C,MEMTYPE,D,E> InterpStencil<T>::operator()(
        const BoxData<T,C,MEMTYPE,D,E>& a_src,
        Box                             a_box,
        T                               a_scale)
{
    if (!m_closed){this->close();}

    LazyStencil<T,C,MEMTYPE,D,E> ret;
    if (m_s.size() <= 0){return ret;}
    Box bsrc = a_src.box();
    for (int ii = 0; ii < m_s.size(); ii++)
    {
        Box b = Box(bsrc.low() - m_s[ii].span().low(), bsrc.high() - m_s[ii].span().high());
        b = b.taperCoarsen(m_s[ii].srcRatio());
        ret.m_stencil.push_back(&(m_s[ii]));
        ret.m_src.push_back(const_cast<BoxData<T,C,MEMTYPE,D,E>*>(&a_src));
        ret.m_box.push_back(a_box & b);
        ret.m_scale.push_back(a_scale);
    }
    return ret;
}

template<class T>
template<unsigned int C, MemType MEMTYPE, unsigned char D, unsigned char E>
inline LazyStencil<T,C,MEMTYPE,D,E> InterpStencil<T>::operator()(
        const BoxData<T,C,MEMTYPE,D,E>& a_src,
		T                               a_scale)
{
    if (!m_closed){this->close();}
    
    LazyStencil<T,C,MEMTYPE,D,E> ret;
    if (m_s.size() <= 0){return ret;}
    Box bsrc = a_src.box();

    for (int ii = 0; ii < m_s.size(); ii++)
    {
        //compute range for building the output BoxData, if necessary
        if (ii == 0)
        {
            ret.m_range = m_s[ii].range(bsrc); //initialize
        } else {
            ret.m_range = ret.m_range & m_s[ii].range(bsrc).low();
            ret.m_range = ret.m_range & m_s[ii].range(bsrc).high();
        }

        //compute the largest possible computation box (different from range) 
        Box b = Box(bsrc.low() - m_s[ii].span().low(), bsrc.high() - m_s[ii].span().high());
        b = b.taperCoarsen(m_s[ii].srcRatio());
        ret.m_stencil.push_back(&(m_s[ii]));
        ret.m_src.push_back(const_cast<BoxData<T,C,MEMTYPE,D,E>*>(&a_src));
        ret.m_box.push_back(b);
        ret.m_scale.push_back(a_scale);
    }
    return ret;
}
