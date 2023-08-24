PROTO_KERNEL_START
void
maxTagFcnF(Var<short,1>& a_tag, Var<short,1>& a_coarse)
{
  if (a_coarse(0) > 0) a_tag(0) = 1;
}
PROTO_KERNEL_END(maxTagFcnF,maxTagFcn)

template<typename T, unsigned int C, MemType MEM, Centering CTR>
void
AMRGrid::computeTags(
        LevelTagData&                 a_tags, 
        LevelBoxData<T, C, MEM, CTR>& a_data,
        Point                         a_bufferSize,
        T                             a_threshold,
        unsigned int                  a_comp)
{
    a_tags.define(a_data.layout(), a_bufferSize);

    std::vector<Stencil<T>> D;
    for (int dir = 0; dir < DIM; dir++)
    {
        Stencil<T> Si = ((T)1)*Shift::Basis(dir, 1) - ((T)1)*Shift::Basis(dir, -1);
        D.push_back(Si);
    }

    for (auto iter = a_data.begin(); iter.ok(); ++iter)
    {
        const auto& data_i = slice(a_data[*iter], a_comp);
        auto& tags_i = a_tags[*iter];
        BoxData<T, DIM, MEM> diffs(iter.box()); 
        for (int dir = 0; dir < DIM; dir++)
        {
            BoxData<T, 1> diff = slice(diffs, dir);
            diff |= D[dir](data_i);
        }
        forallInPlace(
        [] PROTO_LAMBDA (Var<short, 1>& v_tags, Var<T, DIM>& v_diffs, T v_threshold)
        {
            T diffSq = 0;
            for (int dir = 0; dir < DIM; dir++)
            {
                diffSq += v_diffs(dir)*v_diffs(dir);
            }
            T diff = sqrt(diffSq);
            if (diff > v_threshold)
            {
                v_tags(0) = 1;
            } else {
                v_tags(0) = 0;
            }
        }, tags_i, diffs, a_threshold);
    }
    
    AMRGrid::buffer(a_tags, a_bufferSize);
}

template<typename T, unsigned int C, MemType MEM, Centering CTR>
void
AMRGrid::computeTags(
        LevelTagData&                 a_tags, 
        LevelBoxData<T, C, MEM, CTR>& a_data,
        int                           a_bufferSize,
        T                             a_threshold,
        unsigned int                  a_comp)
{
    PROTO_ASSERT(a_bufferSize > 0,
        "AMRGrid::computeTags | Error: Buffer size %i is invalid. Must be non-negative.",
        a_bufferSize);
    computeTags(a_tags, a_data, Point::Ones(a_bufferSize), a_threshold, a_comp);
}

void
AMRGrid::buffer(LevelTagData& a_tags, int a_bufferSize)
{
    AMRGrid::buffer(a_tags, Point::Ones(a_bufferSize));
}
void
AMRGrid::buffer(LevelTagData& a_tags, Point a_bufferSize)
{
  
    Box K(-a_bufferSize, a_bufferSize);
    Stencil<short> S;
    short w = 1;
    for (auto iter = K.begin(); iter.ok(); ++iter)
    {
        S += w*Shift(*iter);
    }
    a_tags.exchange();
    for (auto iter = a_tags.begin(); iter.ok(); ++iter)
    {
        auto& patch = a_tags[*iter];
        TagData tmp(a_tags.layout()[*iter]);
        tmp.setVal(0);
        tmp |= S(patch);
        forallInPlace(
                [=] PROTO_LAMBDA (Var<short>& a_data)
                {
                a_data(0) = min((int)a_data(0), 1);
                }, tmp);
        tmp.copyTo(patch);
    }
}

AMRGrid::AMRGrid(
    std::vector<DisjointBoxLayout>& a_layouts,
    std::vector<Point>& a_refRatios,
    int a_maxLevels)
{
    define(a_layouts, a_refRatios, a_maxLevels);
}

AMRGrid::AMRGrid(
    DisjointBoxLayout& a_layout,
    std::vector<Point>& a_refRatios,
    int a_maxLevels)
{
    define(a_layout, a_refRatios, a_maxLevels);
}

void
AMRGrid::define(
    std::vector<DisjointBoxLayout>& a_layouts,
    std::vector<Point>& a_refRatios,
    int a_maxLevels)
{
    PROTO_ASSERT(a_maxLevels > 0,
        "AMRGrid::define | Error: maxLevels should not be zero.");
    PROTO_ASSERT(a_layouts.size() <= a_maxLevels,
        "AMRGrid::define | Error: Input layout array has %lu levels, but maxLevels is set to %i.",
        a_layouts.size(), a_maxLevels);
    PROTO_ASSERT(a_refRatios.size() == a_maxLevels - 1,
        "AMRGrid::define | Error: the number of ref ratios should equal maxLevels - 1");
    for (int lvl = 1; lvl < a_layouts.size(); lvl++)
    {
        auto crseDomainBox = a_layouts[lvl-1].domain().box();
        auto fineDomainBox = a_layouts[lvl].domain().box();
        PROTO_ASSERT(crseDomainBox == fineDomainBox.coarsen(a_refRatios[lvl-1]),
            "AMRGrid::define | Error: Incorrect refinement ratio on level %i.", lvl-1);
        /*
        PROTO_ASSERT(a_refRatios[lvl-1].isotropic(),
            "AMRGrid::define | Error: Anisotropic refinement ratios are not yet supported.");
        */
    }
    
    m_maxLevels = a_maxLevels;
    m_layouts.clear();
    m_refRatios.clear();
    for (int ii = 0; ii < a_layouts.size(); ii++)
    {
        m_layouts.push_back(a_layouts[ii]); //copying here is fine; dbl internals do their job
    }
    for (int ii = 1; ii < a_maxLevels; ii++)
    {
        m_refRatios.push_back(a_refRatios[ii-1]);
    }
    m_defined = true;
}

void
AMRGrid::define(
    DisjointBoxLayout& a_layout,
    std::vector<Point>& a_refRatios,
    int a_maxLevels)
{
    std::vector<DisjointBoxLayout> layouts;
    layouts.push_back(a_layout);
    define(layouts, a_refRatios, a_maxLevels);
}

const DisjointBoxLayout&
AMRGrid::operator[](unsigned int a_level) const
{
    PROTO_ASSERT(a_level < m_layouts.size(),
        "AMRGrid::operator[] | Error: level index %i out of bounds. Should be in [0, %lu)",
        a_level, m_layouts.size());

    return m_layouts[a_level];
}

DisjointBoxLayout&
AMRGrid::operator[](unsigned int a_level)
{
    PROTO_ASSERT(a_level < m_layouts.size(),
        "AMRGrid::operator[] | Error: level index %i out of bounds. Should be in [0, %lu)",
        a_level, m_layouts.size());

    return m_layouts[a_level];
}

void
AMRGrid::regrid(LevelTagData& a_bufferedTags, unsigned int a_level)
{
    PROTO_ASSERT(a_level < m_layouts.size(),
        "AMRGrid::regrid | Error: Input level %u exceeds size %lu of AMRGrid.", a_level, m_layouts.size());
    //TODO: Maybe regridding on the finest level should be a null-op instead of an error.  
    PROTO_ASSERT(a_level+1 < m_maxLevels,
        "AMRGrid::regrid | Error: Attempt to regrid on the finest level.");
    Point boxSize = operator[](a_level).boxSize();
    regrid(a_bufferedTags, a_level, boxSize);

}
void
AMRGrid::regrid(LevelTagData& a_bufferedTags, unsigned int a_level, Point a_fineBoxSize)
{
    PROTO_ASSERT(a_level < m_layouts.size(),
        "AMRGrid::regrid | Error: Input level %u exceeds size %lu of AMRGrid.", a_level, m_layouts.size());
    //TODO: Maybe regridding on the finest level should be a null-op instead of an error.  
    PROTO_ASSERT(a_level+1 < m_maxLevels,
        "AMRGrid::regrid | Error: Attempt to regrid on the finest level.");
    auto& crseLayout = operator[](a_level);
    PROTO_ASSERT(crseLayout == a_bufferedTags.layout(),
        "AMRGrid::regrid | Error: Tag data layout does not match AMRGrid layout on level %u.", a_level);
   
    Point patchRefRatio = (crseLayout.boxSize() / a_fineBoxSize) * m_refRatios[a_level];
    Box crseBitDomain = crseLayout.patchDomain().box();
    Box fineBitDomain = crseBitDomain.refine(patchRefRatio);
    
    FinitePointSet taggedPatches(fineBitDomain, crseLayout.domain().periodicity());
    Point tagRefRatio = crseLayout.boxSize() / patchRefRatio;
    for ( auto iter = a_bufferedTags.begin(); iter.ok(); ++iter)
    {
        auto& patch = a_bufferedTags[*iter];
        TagData tags(a_bufferedTags.layout()[*iter].coarsen(tagRefRatio));
        tags |= Stencil<short>::Sum(tagRefRatio)(patch);
         
        for (auto biter = tags.box().begin(); biter.ok(); ++biter)
        {
            if (tags(*biter) != 0)
            {
                taggedPatches.add(*biter);
            }
        }
    }
  
    //ProblemDomain fineProblemDomain = crseLayout.domain().refine(Point::Ones(PR_AMR_REFRATIO));
    ProblemDomain fineProblemDomain = crseLayout.domain().refine(m_refRatios[a_level]);
    DisjointBoxLayout fineLayout(fineProblemDomain, taggedPatches.points(), a_fineBoxSize);
    if (a_level + 1 == m_layouts.size())
    {
        // adding a new level of refinement with the next coarser level's boxSize
        m_layouts.push_back(fineLayout);
    } else {
        // regridding an existing level of refinement
        m_layouts[a_level + 1] = fineLayout;
    }
}

Point
AMRGrid::refRatio(int a_level) const
{
    PROTO_ASSERT(a_level < m_refRatios.size(),
        "AMRGrid::refRatio | Error: level %i is out of bounds. Only *coarse* levels have refinement ratios.", 
        a_level)
    return m_refRatios[a_level];
}

bool
AMRGrid::compatible(const AMRGrid& a_grid) const
{
    int n = min(numLevels(), a_grid.numLevels());
    for (int lvl = 0; lvl < n; lvl++)
    {
        if (!m_layouts[lvl].compatible(a_grid[lvl])) { return false; }
    }
    return true;
}

void
AMRGrid::enforceNesting(unsigned int a_level, int a_nestingDistance)
{
    // input level is the "middle" level (e.g. the one being modified)
    PROTO_ASSERT(a_level+1 < m_layouts.size(),
        "AMRGrid::enforceNesting | Error: Invalid input level %u", a_level);
    PROTO_ASSERT(a_level > 0,
        "AMRGrid::enforceNesting | Error: Invalid input level %u", a_level);

    auto& layout = m_layouts[a_level];
    auto& fineLayout = m_layouts[a_level+1];
    
    ProblemDomain patchDomain = layout.domain().coarsen(layout.boxSize());

    const auto& patches = layout.boxes();
    const auto& finePatches = fineLayout.boxes();
    
    std::set<Point> points; 
    
    // patch-space refinement ratio between this level and fine level
    Point finePatchRatio = layout.boxSize() / fineLayout.boxSize() * m_refRatios[a_level];

    // Add this layout's patches in the coarse point space
    for (auto iter = patches.begin(); iter != patches.end(); ++iter)
    {
        points.insert((*iter).first);
    }
   
    // Union with the fine layout's patches in the coarse point space 
    int nestingVect[DIM];
    for (int ii = 0; ii < DIM; ii++)
    {
        int n = ceil((double)(a_nestingDistance * m_refRatios[a_level][ii]) / (double)(fineLayout.boxSize()[ii]));
        nestingVect[ii] = n;
    }
    Point growRadius(nestingVect);
    for (auto iter = finePatches.begin(); iter != finePatches.end(); ++iter)
    {
        Point finePoint = (*iter).first;
        auto K = Box(-growRadius, growRadius);
        K = K.shift(finePoint);
        for (auto biter = K.begin(); biter.ok(); ++biter)
        {
            Point pi = (*biter) / finePatchRatio;
            if (patchDomain.contains(pi))
            {
                points.insert(patchDomain.image(pi));
            }
        }
    }

    std::vector<Point> pointVect;
    for (auto iter = points.begin(); iter != points.end(); ++iter)
    {
        pointVect.push_back(*iter);
    }
    
    // Redefine the layout 
    layout.define(layout.domain(), pointVect, layout.boxSize());
}

void AMRGrid::addFinerGrids(LevelTagData& a_tags,int a_level)
{

    PR_assert(a_level +2 < m_maxLevels);
    PR_assert(a_tags.layout() == (*this)[a_level]);
    // midTags2: Tags on a_level + 1 tagging a_level+2 grids.
    LevelTagData midTags2((*this)[a_level+2].coarsen(m_refRatios[a_level+1]),
            Point::Zeros());
    //midTags2.setToZero();
    for (auto dit = (*this)[a_level+2].begin();dit != dit.end();++dit)
    {
        midTags2[*dit].setVal(1);
    }
    // midTags1:midTags2 tags on a_level + 1 dbl.
    LevelTagData midTags1((*this)[a_level+1],m_refRatios[a_level]);
    midTags1.setToZero();
    midTags2.copyTo(midTags1);
    // coarseTags1: midTags1 projected onto a_level grid.
    LevelTagData coarseTags1((*this)[a_level+1].coarsen(m_refRatios[a_level]),
            m_refRatios[a_level+1]);
    coarseTags1.setToZero();
    for (auto dit = (*this)[a_level+1].begin();dit != dit.end();++dit)
    {
        coarseTags1[*dit] |= Stencil<short>::Sum(m_refRatios[a_level])(midTags1[*dit]);
        forallInPlace(maxTagFcn,coarseTags1[*dit],coarseTags1[*dit]);
    }
    // coarseTags0: coarseTags1 defined on the a_level dbl.

    Point ptCoarsen2 =
        (*this)[a_level+2].boxSize()/m_refRatios[a_level+1]/m_refRatios[a_level];
    if (ptCoarsen2*m_refRatios[a_level+1]*m_refRatios[a_level] !=
            (*this)[a_level+2].boxSize()) ptCoarsen2 = ptCoarsen2 + Point::Ones();
    LevelTagData coarseTags0((*this)[a_level],ptCoarsen2);
    coarseTags0.setToZero();
    coarseTags1.copyTo(coarseTags0);
    coarseTags0.exchange();
    // We buffer so that enforceNesting will not throw away any of the
    // a_level+2 grids.
    buffer(coarseTags0, ptCoarsen2);
    for (auto dit = (*this)[a_level].begin();dit != dit.end();++dit)
    {
        forallInPlace(maxTagFcn,a_tags[*dit],coarseTags0[*dit]);
    }
}
void
AMRGrid::enforceNesting2(unsigned int a_level,int a_nestingDistance)
{
  PROTO_ASSERT(a_level > 0,
        "AMRGrid::enforceNesting | Error: Invalid input level %u", a_level);
  auto dbl = m_layouts[a_level];
  auto dblCoarse = m_layouts[a_level-1];
  auto oldvector = dbl.boxes();
  Point refRatio = m_refRatios[a_level-1];
  Point coarsen = refRatio * dblCoarse.boxSize()/dbl.boxSize();
  Box bx(Point::Ones(-a_nestingDistance),Point::Ones(a_nestingDistance));
  std::vector<Point> newPoints;
  for (int kold = 0; kold < oldvector.size();kold++)
    {
      bool isNested = true;
      Point pt = oldvector[kold].first;
      for (auto bxit = bx.begin();bxit.ok();++bxit)
        {
          Point ptShiftCoarse = (pt + *bxit)/coarsen;
          auto di = dblCoarse.find(ptShiftCoarse);
          if (di == *dblCoarse.end()) isNested=false;
        }
      if (isNested) newPoints.push_back(pt);
    }
  DisjointBoxLayout dblNew(dbl.domain(),newPoints,dbl.boxSize());
  m_layouts[a_level] = dblNew;
}

void AMRGrid::addFinerTags(
        LevelTagData& a_tags,
        unsigned int  a_level)
{
    PROTO_ASSERT(a_tags.layout().compatible((*this)[a_level]),
            "AMRGrid::addFinerTags | Error: \
            Input tag layout is incompatible with the layout on level %i.", a_level);
    if (a_level+2 >= numLevels()) { return;}

    Point refRatio_12 = m_refRatios[a_level+1];
    Point refRatio_01 = m_refRatios[a_level];
    DisjointBoxLayout cfLayout_12 = (*this)[a_level+2].coarsen(refRatio_12);
    DisjointBoxLayout cfLayout_01 = (*this)[a_level+1].coarsen(refRatio_01);
    DisjointBoxLayout& layout_2 = (*this)[a_level+1];
    DisjointBoxLayout& layout_1 = (*this)[a_level+1];
    DisjointBoxLayout& layout_0 = (*this)[a_level];
    auto SUM_01 = Stencil<short>::Sum(refRatio_01);

    LevelTagData cfTags_12(cfLayout_12, Point::Zeros());
    cfTags_12.setVal(1);
    LevelTagData tags_1(layout_1, Point::Zeros());
    tags_1.setToZero();
    cfTags_12.copyTo(tags_1);
    LevelTagData cfTags_01(cfLayout_01, Point::Zeros());
    for (auto iter = tags_1.begin(); iter.ok(); ++iter)
    {
        auto& crseTags = cfTags_01[*iter];
        auto& fineTags = tags_1[*iter];
        crseTags |= SUM_01(fineTags);
    }
    //cfTags_01.clampVal(0,1);

    Point refRatio_02 = refRatio_12*refRatio_01;
    Point bufferSize = layout_2.boxSize()/refRatio_02;
    if (bufferSize*refRatio_02 != layout_2.boxSize())
    {
        bufferSize = bufferSize + Point::Ones();
    }

    LevelTagData tags_0(layout_0, bufferSize);
    tags_0.setToZero();
    cfTags_01.copyTo(tags_0);
    buffer(tags_0, bufferSize);
    a_tags.increment(tags_0);
    //a_tags.clampVal(0,1);
}

void AMRGrid::print() const
{
    pout() << "AMRGrid | numLevels: " << numLevels() << " | maxLevels: " << maxLevels() << std::endl;
    for (int lvl = 0; lvl < numLevels(); lvl++)
    {
        (*this)[lvl].print();
    }
}

