
template<typename T, unsigned int C, MemType MEM, Centering CTR>
void AMRData<T, C, MEM, CTR>::define(AMRGrid& a_grid, Point a_ghost)
{
    m_ghost = a_ghost;
    m_grid = a_grid;
    m_data.clear();
    for (int ii = 0; ii < a_grid.numLevels(); ii++)
    {
        auto level = std::make_shared<LevelBoxData<T, C, MEM, CTR>>(a_grid[ii], a_ghost);
        m_data.push_back(level);
    }

    m_defined = true;
}
template<typename T, unsigned int C, MemType MEM, Centering CTR>
void AMRData<T, C, MEM, CTR>::regrid(AMRGrid& a_newgrid,int a_level,int a_order)
{
  if (m_defined)
    {
      int ghostInterp = 2;
        for (int ii = a_level;ii + 1 < a_newgrid.numLevels();ii++)          
        {
          Point refratio = a_newgrid.refRatio(ii);
          InterpStencilTP<double > stencil(a_order,refratio[0]);
          PR_assert(refratio == Point::Ones()*refratio[0]);
          
        //PC: Need to extend InterpStencilTP to non-isotropic case !!
          
          auto newleveldata =std::make_shared<LevelBoxData<T, C, MEM, CTR> >
            (a_newgrid[ii+1], m_ghost);
          auto dblCoarse =
            a_newgrid[ii+1].coarsen(a_newgrid.refRatio(ii));
          LevelBoxData<T, C, MEM, CTR > dataCoarse(dblCoarse,
                                                   Point::Ones(ghostInterp));
          (*this)[ii].copyTo(dataCoarse);
          for (auto dit=dataCoarse.begin();dit != dit.end();++dit)
            {
              stencil.apply((*newleveldata)[*dit],dataCoarse[*dit]);
            }
          // copy old data to new data on intersection.
          (*this)[ii+1].copyTo(*newleveldata);
          // replace the data pointer with the new data pointer.
          m_data[ii+1] = newleveldata;
        }
        m_grid = a_newgrid;
    }
  else
    {
       MayDay<void>::Error("AMRData::regrid | Error: attempting to regrid undefined AMRData.");
    }
  // Check to see if grids are consistent.
  PR_assert(m_grid.numLevels() == (*this).numLevels()); 
  for (int ii = 0; ii < m_grid.numLevels();ii++)
    {
      if (!(m_data[ii]->layout() == m_grid[ii]))
        {
          (*m_data[ii]).layout().print();
          m_grid[ii].print();
          PR_assert(m_data[ii]->layout() == m_grid[ii]); 
        }
    }
}
template<typename T, unsigned int C, MemType MEM, Centering CTR>
LevelBoxData<T, C, MEM, CTR>&
AMRData<T, C, MEM, CTR>::operator[](unsigned int a_level)
{
    PROTO_ASSERT(a_level < m_data.size(),
        "AMRData::operator[] | Error: level %u is out of bounds.", a_level);
    return *m_data[a_level];
}

template<typename T, unsigned int C, MemType MEM, Centering CTR>
const LevelBoxData<T, C, MEM, CTR>&
AMRData<T, C, MEM, CTR>::operator[](unsigned int a_level) const
{
    PROTO_ASSERT(a_level < m_data.size(),
        "AMRData::operator[] | Error: level %u is out of bounds.", a_level);
    return *m_data[a_level];
}

template<typename T, unsigned int C, MemType MEM, Centering CTR>
void
AMRData<T, C, MEM, CTR>::setToZero()
{
    for (int ii = 0; ii < m_data.size(); ii++)
    {
        m_data[ii]->setToZero();
    }
}

template<typename T, unsigned int C, MemType MEM, Centering CTR>
template<typename Func, typename... Srcs>
void
AMRData<T, C, MEM, CTR>::initialize(T a_dx0, Func& a_func, Srcs... a_srcs)
{
    auto dx = a_dx0;
    for (int ii = 0; ii < m_data.size(); ii++)
    {
        m_data[ii]->initialize(a_func, dx, a_srcs...);
        if (ii < m_data.size()-1)
        {
            dx /= m_grid.refRatio(ii)[0];
            PROTO_ASSERT(m_grid.refRatio(ii).isotropic(),
                    "AMRData::initialize | Error: Non-isotropic refinement ratio used in \
                    isotropic initialization function.");
        }
    }
}

template<typename T, unsigned int C, MemType MEM, Centering CTR>
template<typename Func, typename... Srcs>
void
AMRData<T, C, MEM, CTR>::initialize(std::array<T, DIM>& a_dx0, Func& a_func, Srcs... a_srcs)
{
    auto dx = a_dx0;
    for (int ii = 0; ii < m_data.size(); ii++)
    {
        m_data[ii]->initialize(a_func, dx, a_srcs...);
        if (ii < m_data.size()-1)
        {
            for (int dir = 0; dir < DIM; dir++)
            {
                dx[dir] /= m_grid.refRatio(ii)[dir];
            }
        }
    }
}

template<typename T, unsigned int C, MemType MEM, Centering CTR>
template<typename Func, typename... Srcs>
void
AMRData<T, C, MEM, CTR>::initConvolve(double a_dx0, Func& a_func, Srcs... a_srcs)
{
    double dx = a_dx0;
    for (int ii = 0; ii < m_data.size(); ii++)
    {
        m_data[ii]->initConvolve(a_func, dx, a_srcs...);
        if (ii < m_data.size()-1)
        {
            // FIXME: Assumes isotropic refinement
            dx /= m_grid.refRatio(ii)[0];
        }
    }
}

template<typename T, unsigned int C, MemType MEM, Centering CTR>
void
AMRData<T, C, MEM, CTR>::averageDown()
{
    if (numLevels() < 2) { return; }
    for (int lvl = numLevels() - 2; lvl >= 0; lvl--)
    {
        auto& crse = operator[](lvl);
        auto& fine = operator[](lvl+1);
        Point refRatio = m_grid.refRatio(lvl);
        auto  cfLayout = fine.layout().coarsen(refRatio);
        LevelBoxData<T, C, MEM, CTR> temp(cfLayout, Point::Zeros());
        //crse.layout().print();
        //fine.layout().print();
        //cfLayout.print();

        Proto::averageDown(crse, fine, temp, refRatio);
    }
}

template<typename T, unsigned int C, MemType MEM, Centering CTR>
void
AMRData<T, C, MEM, CTR>::increment(
    AMRData<T, C, MEM, CTR>& a_data,
    T a_scale)
{
    int maxLevels = std::min(numLevels(), a_data.numLevels());
    for (int lvl = 0; lvl < maxLevels; lvl++)
    {
        auto& lhs = (*this)[lvl];
        auto& rhs = a_data[lvl];
        lhs.increment(rhs, a_scale);  
    }
}

template<typename T, unsigned int C, MemType MEM, Centering CTR>
T
AMRData<T, C, MEM, CTR>::integrate(T a_cdx, unsigned int a_c)
{
    std::array<T, DIM> dx;
    dx.fill(a_cdx);
    return integrate(dx, a_c);
}

template<typename T, unsigned int C, MemType MEM, Centering CTR>
T
AMRData<T, C, MEM, CTR>::integrate(const std::array<T, DIM>& a_cdx, unsigned int a_c)
{
    averageDown();
    double integral = (*this)[0].integrate(a_cdx, a_c);
    return integral;
}

template<typename T, unsigned int C, MemType MEM, Centering CTR>
T
AMRData<T, C, MEM, CTR>::integrateAbs(T a_cdx, unsigned int a_c)
{
    std::array<T, DIM> dx;
    dx.fill(a_cdx);
    return integrateAbs(dx, a_c);
}

template<typename T, unsigned int C, MemType MEM, Centering CTR>
T
AMRData<T, C, MEM, CTR>::integrateAbs(const std::array<T, DIM>& a_cdx, unsigned int a_c)
{
    averageDown();
    double integral = (*this)[0].integrateAbs(a_cdx, a_c);
    return integral;
}

template<typename T, unsigned int C, MemType MEM, Centering CTR>
T
AMRData<T, C, MEM, CTR>::absMax(unsigned int a_c)
{
    averageDown();
    T maxValue = 0;
    for (int lvl = 0; lvl < maxLevels(); lvl++)
    {
        // I *think* the mean value theorem guarantees that if an extremum
        // occurs on a coarse patch covered by fine data, that extremum
        // will be *more* extreme on the finer level. Hence, we don't need to
        // worry about valid/invalid data. -CLG
        maxValue = std::max(m_data[lvl]->absMax(), maxValue);
    }
    return maxValue;
}

template<typename T, unsigned int C, MemType MEM, Centering CTR>
void
AMRData<T, C, MEM, CTR>::copyTo(AMRData<T, C, MEM, CTR>& a_rhs)
{
    int minLevels = min(numLevels(), a_rhs.numLevels());
    for (int lvl = 0; lvl < minLevels; lvl++)
    {
        (*this)[lvl].copyTo(a_rhs[lvl]);
    }
}

template<typename T, unsigned int C, MemType MEM, Centering CTR>
void
AMRData<T, C, MEM, CTR>::exchange()
{
    for (int lvl = 0; lvl < numLevels(); lvl++)
    {
        m_data[lvl]->exchange();
    }
}
