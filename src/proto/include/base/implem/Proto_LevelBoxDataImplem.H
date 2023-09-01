template<typename T, unsigned int C, MemType MEM, Centering CTR>
LevelBoxData<T, C, MEM, CTR>::
LevelBoxData()
{
    m_isDefined = false;
}

template<typename T, unsigned int C, MemType MEM, Centering CTR>
LevelBoxData<T, C, MEM, CTR>::
LevelBoxData(const DisjointBoxLayout& a_layout, 
             const Point& a_ghost)
{
    define(a_layout, a_ghost);
}

template<typename T, unsigned int C, MemType MEM, Centering CTR>
template<unsigned int CC>
LevelBoxData<T, C, MEM, CTR>::
LevelBoxData(const LevelBoxData<T, CC, MEM, CTR>& a_src,
             unsigned int a_comp)
{
    PROTO_ASSERT(C+a_comp <= CC,
            "LevelBoxData alias constructor | Error: \
            Invalid number of components.");
    m_isDefined = true;
    m_ghost = a_src.ghost();
    m_layout = a_src.layout();
    for (auto iter : m_layout)
    {
        auto& patch = a_src[iter];
        Point p = patch.box().low();
        // Calling BoxData alias constructor
        m_data.push_back(
                std::make_shared<BoxData<T, C, MEM>>(
                    patch.data(p, a_comp), patch.box())); 
    }
    
    if (m_ghost != Point::Zeros())
    {
        m_exchangeCopier.define(LevelCopierOp<T, C, MEM, MEM, CTR>(*this, *this));
    }
}


template<typename T, unsigned int C, MemType MEM, Centering CTR>
void LevelBoxData<T, C, MEM, CTR>::define(
                 const DisjointBoxLayout& a_layout, 
                 const Point& a_ghost)
{
    m_isDefined = true;
    m_ghost = a_ghost;
    m_layout = a_layout;
    for (auto dit = begin(); dit.ok(); ++dit)
    {
        auto patch = std::make_shared<BoxData<T,C,MEM>>(a_layout[*dit].grow(m_ghost));
        m_data.push_back(patch);
    }

    if (m_ghost != Point::Zeros())
    {
        m_exchangeCopier.define(LevelCopierOp<T, C, MEM, MEM, CTR>(*this, *this));
    }
}

template<typename T, unsigned int C, MemType MEM, Centering CTR>
template<typename Func, typename... Srcs>
void LevelBoxData<T, C, MEM, CTR>::initialize(Func& a_func, Srcs... a_srcs)
{
    for (auto iter : m_layout)
    {
        auto& patch = (*this)[iter];
        forallInPlace_p(a_func, patch, a_srcs...);
    }
}

template<typename T, unsigned int C, MemType MEM, Centering CTR>
template<typename Func, typename... Srcs>
void LevelBoxData<T, C, MEM, CTR>::initConvolve(Func& a_func, Srcs... a_srcs)
{
    for (auto iter : m_layout)
    {
        auto& patch = (*this)[iter];
        BoxData<T, C, MEM> temp(patch.box().grow(1));
        forallInPlace_p(a_func, temp, a_srcs...);
        Operator::convolve(patch, temp);
    }
}

template<typename T, unsigned int C, MemType MEM, Centering CTR>
unsigned int
LevelBoxData<T, C, MEM, CTR>::patchSize() const
{
    int size = 1;
    Point boxSize = layout().boxSize() + 2*ghost();
    for (int ii = 0; ii < DIM; ii++)
    {
        size *= boxSize[ii];
    }
    return C*size;
}

template<typename T, unsigned int C, MemType MEM, Centering CTR>
unsigned int
LevelBoxData<T, C, MEM, CTR>::offset(int a_proc) const
{
    Point boxSizeVect = layout().boxSize() + 2*ghost();
    int boxSize = 1;
    for (int ii = 0; ii < DIM; ii++)
    {
        boxSize *= boxSizeVect[ii];
    }
    return layout().offset(a_proc)*boxSize*C;
}

template<typename T, unsigned int C, MemType MEM, Centering CTR>
unsigned int 
LevelBoxData<T, C, MEM, CTR>::size() const
{
    return m_layout.size();
}

template<typename T, unsigned int C, MemType MEM, Centering CTR>
void   
LevelBoxData<T, C, MEM, CTR>::setToZero(int a_comp)
{
    setVal(0, a_comp);
}

template<typename T, unsigned int C, MemType MEM, Centering CTR>
void   
LevelBoxData<T, C, MEM, CTR>::setVal(T a_value, int a_comp)
{
    PROTO_ASSERT((a_comp >= -1) && (a_comp < DIM),
            "LevelBoxData::setVal | Error: %i is not a valid component specification.", a_comp);
    for (auto iter : m_layout)
    {
        auto& patch = (*this)[iter];
        if (a_comp == -1)
        {
            patch.setVal(a_value);
        } else {
            auto patch_i = slice(patch, a_comp);
            patch_i.setVal(a_value);
        }
    }
}

template<typename T, unsigned int C, MemType MEM, Centering CTR>
void   
LevelBoxData<T, C, MEM, CTR>::setBoundary(T a_value, int a_comp)
{
    PROTO_ASSERT((a_comp >= -1) && (a_comp < DIM),
            "LevelBoxData::setBoundary | Error: \
            %i is not a valid component specification.", a_comp);
    if (ghost() == Point::Zeros()){ return; } //nothing to do
    for (auto iter : m_layout)
    {
        if (!layout().onLevelBoundary(m_layout.point(iter))) { continue; } //no boundary
        auto& patch = (*this)[iter];
        if (a_comp == -1)
        {
            BoxData<T, C, MEM> tmp(layout()[iter]);
            patch.copyTo(tmp);
            patch.setVal(a_value);
            tmp.copyTo(patch);
        } else {
            auto patch_i = slice(patch, a_comp);
            BoxData<T, 1, MEM> tmp(layout()[iter]);
            patch_i.copyTo(tmp);
            patch_i.setVal(a_value);
            tmp.copyTo(patch_i);
        }
    }
    //TODO: this could be coded more elegantly / efficiently 
    exchange();
}

/*
template<typename T, unsigned int C, MemType MEM, Centering CTR>
void   
LevelBoxData<T, C, MEM, CTR>::clampVal(T a_min, T a_max, int a_comp)
{
    PROTO_ASSERT((a_comp >= -1) && (a_comp < DIM),
            "LevelBoxData::clampVal | Error: %i is not a valid component specification.", a_comp);
    for (auto iter : m_layout)
    {
        auto& patch = (*this)[iter];
        if (a_comp == -1)
        {
            patch.clampVal(a_min, a_max);
        } else {
            auto patch_i = slice(patch, a_comp);
            patch_i.clampVal(a_min, a_max);
        }
    }
}
*/

template<typename T, unsigned int C, MemType MEM, Centering CTR>
void 
LevelBoxData<T, C, MEM, CTR>::exchange()
{
    PR_TIME("LevelBoxData::exchange");
#ifdef PR_MPI
    {
        PR_TIME("MPI_Barrier exchange");
        barrier();
    }
#endif
    m_exchangeCopier.execute();
}

// TODO: The following definition of copyTo is a stopgap until we teach Copier
// and CopierOp how to move data between HOST and DEVICE
template<typename T, unsigned int C, MemType SRC_MEM, Centering CTR>
template<MemType DST_MEM>
void 
LevelBoxData<T, C, SRC_MEM, CTR>::copyToSimple(LevelBoxData<T, C, DST_MEM, CTR>& a_dest) const
{

    PROTO_ASSERT(a_dest.layout().compatible(this->layout()),
        "LevelBoxData::copyToSimple | Error: inputs are not compatible.");
    for (auto iter : m_layout)
    {
        (*this)[iter].copyTo(a_dest[iter]);
    }
    a_dest.exchange();
}

template<typename T, unsigned int C, MemType SRC_MEM, Centering CTR>
template<MemType DST_MEM>
void 
LevelBoxData<T, C, SRC_MEM, CTR>::copyTo(LevelBoxData<T, C, DST_MEM, CTR>& a_dest) const
{
    if (a_dest.layout().compatible(this->layout()))
    {
        this->copyToSimple(a_dest);
        return;
    }
    if (SRC_MEM == DST_MEM)
    {
        auto mutableSrc = const_cast<LevelBoxData<T, C, SRC_MEM, CTR>*> (this);
        LevelCopierOp<T, C, SRC_MEM, DST_MEM, CTR> op(*mutableSrc, a_dest);
        LevelCopier<T, C, SRC_MEM, DST_MEM, CTR> copier;
        copier.define(op);
        copier.execute();
    } else {
        LevelBoxData<T, C, DST_MEM, CTR> tmpSrc(this->layout(), this->ghost());
        this->copyToSimple(tmpSrc);
#ifdef PROTO_CUDA
        //cudaDeviceSynchronize();
	protoThreadSynchronize();
#endif
        tmpSrc.copyTo(a_dest);
    }
}

template<typename T, unsigned int C, MemType MEM, Centering CTR>
void 
LevelBoxData<T, C, MEM, CTR>::coarsenTo(
    LevelBoxData<T, C, MEM, CTR>&   a_dest,
    Point                           a_refRatio) const
{
    if (CTR == PR_CELL)
    {
        auto AVG = Stencil<T>::AvgDown(a_refRatio);
        DisjointBoxLayout cfLayout = m_layout.coarsen(a_refRatio);
        LevelBoxData<T, C, MEM, CTR> cfLevel(cfLayout, Point::Zeros());
        for (auto iter : m_layout)
        {
            auto& fine = (*this)[iter]; 
            auto& crse = cfLevel[iter];
            
            crse |= AVG(fine);
        }
        cfLevel.copyTo(a_dest);
    } else {
        MayDay<void>::Warning("LevelBoxData::coarsenTo | Warning: Non-Cell centered implementation not yet implemented");
    }
}

template<typename T, unsigned int C, MemType MEM, Centering CTR>
void 
LevelBoxData<T, C, MEM, CTR>::linearIn(
    void* a_buffer,
    const Box& a_box,
    const DataIndex<BoxPartition>& a_index)
{
    auto& patch = this->operator[](a_index);
    CInterval comps(0, C-1);
    patch.linearIn(a_buffer, a_box, comps);
}

template<typename T, unsigned int C, MemType MEM, Centering CTR>
void 
LevelBoxData<T, C, MEM, CTR>::linearIn(
    void* a_buffer,
    const DataIndex<BoxPartition>& a_index)
{
    Box b = this->operator[](a_index).box();
    linearIn(a_buffer, b, a_index);
}

template<typename T, unsigned int C, MemType MEM, Centering CTR>
void 
LevelBoxData<T, C, MEM, CTR>::linearIn( void* a_buffer)
{
    char* readPoint = (char*)a_buffer;
    for (auto iter : m_layout)
    {
        linearIn(readPoint, iter);
        readPoint += linearSize(iter); 
    }
}

template<typename T, unsigned int C, MemType MEM, Centering CTR>
void 
LevelBoxData<T, C, MEM, CTR>::linearOut(
    void* a_buffer,
    const Box& a_box,
    const DataIndex<BoxPartition>& a_index) const
{
    auto& patch = this->operator[](a_index);
    CInterval comps(0, C-1);
    patch.linearOut(a_buffer, a_box, comps);
}

template<typename T, unsigned int C, MemType MEM, Centering CTR>
void 
LevelBoxData<T, C, MEM, CTR>::linearOut(
    void* a_buffer,
    const DataIndex<BoxPartition>& a_index) const
{
    Box b = this->operator[](a_index).box();
    linearOut(a_buffer, b, a_index);
}

template<typename T, unsigned int C, MemType MEM, Centering CTR>
void 
LevelBoxData<T, C, MEM, CTR>::linearOut( void* a_buffer) const
{
    char* writePoint = (char*)a_buffer;
    for (auto iter : m_layout)
    {
        linearOut(writePoint, iter);
        writePoint += linearSize(iter); 
    }
}

template<typename T, unsigned int C, MemType MEM, Centering CTR>
size_t
LevelBoxData<T, C, MEM, CTR>::linearSize(
    const Box& a_box,
    const DataIndex<BoxPartition>& a_index) const
{
    return C*a_box.size()*sizeof(T);
}

template<typename T, unsigned int C, MemType MEM, Centering CTR>
size_t
LevelBoxData<T, C, MEM, CTR>::linearSize( const DataIndex<BoxPartition>& a_index) const
{
    Box b = this->operator[](a_index).box();
    return this->linearSize(b, a_index);
}

template<typename T, unsigned int C, MemType MEM, Centering CTR>
size_t
LevelBoxData<T, C, MEM, CTR>::linearSize() const
{
    size_t totalSize = 0;
    for (auto iter : m_layout)
    {
        totalSize += (*this)[iter].linearSize();
    }
    return totalSize;
}

//// Interaction with iterators.
template<typename T, unsigned int C, MemType MEM, Centering CTR>
DataIterator<BoxPartition> LevelBoxData<T, C, MEM, CTR>::begin() const
{
    return m_layout.begin();
} 

template<typename T, unsigned int C, MemType MEM, Centering CTR>
DataIterator<BoxPartition> LevelBoxData<T, C, MEM, CTR>::end() const
{
    return m_layout.end();
} 

template<typename T, unsigned int C, MemType MEM, Centering CTR>
BoxData<T,C,MEM>& 
LevelBoxData<T, C, MEM, CTR>::operator[](const DataIndex<BoxPartition>& a_di) 
{
    //return (*(m_data[m_layout.myIntIndex(a_di)]));
    return (*(m_data[a_di]));
}

template<typename T, unsigned int C, MemType MEM, Centering CTR>
const BoxData<T,C,MEM>& 
LevelBoxData<T, C, MEM, CTR>::operator[](const DataIndex<BoxPartition>& a_di) const
{
    //return (*(m_data[m_layout.myIntIndex(a_di)]));
    return (*(m_data[a_di]));
}


template<typename T, unsigned int C, MemType MEM, Centering CTR>
void
LevelBoxData<T, C, MEM, CTR>::increment(LevelBoxData<T, C, MEM, CTR>& a_data, T a_scale)
{
    PROTO_ASSERT(m_layout.compatible(a_data.layout()),
        "LevelBoxData::increment | Error: Incompatible layouts.");
    PR_TIME("LevelBoxData::increment");
    Stencil<T> INCR = a_scale*Shift::Zeros();
    for (auto iter : m_layout)
    {
        auto& lhs_i = (*this)[iter];
        auto& rhs_i = a_data[iter];
        lhs_i += INCR(rhs_i); // using a Stencil here fuses the add and multiply in a single kernel launch -clg
    }
}

template<typename T, unsigned int C, MemType MEM, Centering CTR>
void
LevelBoxData<T, C, MEM, CTR>::operator*=(T a_scale)
{
    PR_TIME("LevelBoxData::operator*=");
    for (auto iter : m_layout)
    {
        auto& lhs_i = (*this)[iter];
        lhs_i *= a_scale;
    }
}

template<typename T, unsigned int C, MemType MEM, Centering CTR>
template <Proto::Operation OP>
double
LevelBoxData<T, C, MEM, CTR>::reduce(unsigned int a_comp) const
{
    PR_TIME("LevelBoxData::reduce");
    Reduction<T, OP, MEM> rxn;
    for (auto iter : m_layout)
    {
        auto& data = (*this)[iter];
        if (ghost() == Point::Zeros())
        {
            // avoid the copyTo call if there are no ghost cells
            data.reduce(rxn, a_comp);
        } else {
            BoxData<T, 1, MEM> temp(m_layout[iter]);
            auto comp = slice(data, a_comp);
            comp.copyTo(temp);
            temp.reduce(rxn);
        }
    }
    return rxn.fetch();
}


template<typename T, unsigned int C, MemType MEM, Centering CTR>
double
LevelBoxData<T, C, MEM, CTR>::absMax(unsigned int a_comp) const
{
    return reduce<Abs>(a_comp);
}

template<typename T, unsigned int C, MemType MEM, Centering CTR>
double
LevelBoxData<T, C, MEM, CTR>::sum(unsigned int a_comp) const
{
    return reduce<Sum>(a_comp);
}

template<typename T, unsigned int C, MemType MEM, Centering CTR>
double
LevelBoxData<T, C, MEM, CTR>::integrate(double a_dx, unsigned int a_comp) const
{
    return reduce<Sum>(a_comp) * pow(a_dx, DIM);
}

template<typename T, unsigned int C, MemType MEM, Centering CTR>
double
LevelBoxData<T, C, MEM, CTR>::integrate(std::array<double, DIM> a_dx, unsigned int a_comp) const
{
    double factor = 1.0;
    for (int dir = 0; dir < DIM; dir++) { factor *= a_dx[dir]; }
    return reduce<Sum>(a_comp) * factor;
}

template<typename T, unsigned int C, MemType MEM, Centering CTR>
double
LevelBoxData<T, C, MEM, CTR>::integrateAbs(double a_dx, unsigned int a_comp) const
{
    return reduce<SumAbs>(a_comp) * pow(a_dx, DIM);
}

template<typename T, unsigned int C, MemType MEM, Centering CTR>
double
LevelBoxData<T, C, MEM, CTR>::integrateAbs(std::array<double, DIM> a_dx, unsigned int a_comp) const
{
    double factor = 1.0;
    for (int dir = 0; dir < DIM; dir++) { factor *= a_dx[dir]; }
    return reduce<SumAbs>(a_comp) * factor;
}

template<typename T, unsigned int C, MemType MEM, Centering CTR>
double
LevelBoxData<T, C, MEM, CTR>::max(unsigned int a_comp) const
{
    return reduce<Max>(a_comp);
}

template<typename T, unsigned int C, MemType MEM, Centering CTR>
double
LevelBoxData<T, C, MEM, CTR>::min(unsigned int a_comp) const
{
    return reduce<Min>(a_comp);
}


// =======================================================================
// LEVEL COPIER OP

template<typename T, unsigned int C, MemType SRC_MEM, MemType DST_MEM, Centering CTR>
LevelCopierOp<T, C, SRC_MEM, DST_MEM, CTR>::LevelCopierOp(
        LevelBoxData<T, C, SRC_MEM, CTR>& a_src,
        LevelBoxData<T, C, DST_MEM, CTR>& a_dst)
{
    PROTO_ASSERT(SRC_MEM == DST_MEM,
        "LevelCopierOp::constructor | Error: \
        LevelCopierOp with different MemType src and dst is not yet implemented.");
    m_src = &a_src;
    m_dst = &a_dst;
}

template<typename T, unsigned int C, MemType SRC_MEM, MemType DST_MEM, Centering CTR>
int
LevelCopierOp<T, C, SRC_MEM, DST_MEM, CTR>::linearSize(const Box& a_box,
        const DataIndex<BoxPartition>& a_index) const
{
    //TODO: Does this do the right thing?
    return C*a_box.size()*sizeof(T);
}

template<typename T, unsigned int C, MemType SRC_MEM, MemType DST_MEM, Centering CTR>
void
LevelCopierOp<T, C, SRC_MEM, DST_MEM, CTR>::linearOut(void* a_buf, 
        const Box& a_bx,
        const DataIndex<BoxPartition>& a_index) const
{
    MemType mem = pointerMemType(a_buf);
    PROTO_ASSERT(mem == SRC_MEM,
        "LevelCopierOp::linearOut | Error: Memtype mismatch. Detected %i but should be %i",
        mem, SRC_MEM);
    auto& patch = m_src->operator[](a_index);
    CInterval cint(0,C-1);
    patch.linearOut(a_buf, a_bx, cint);
    // BEGIN TEST CODE
    /*
    pout() << "LevelCopierOp::linearOut " << std::endl;
    int memsize = a_bx.size()*C*sizeof(T);
    T* tmp = (T*)proto_malloc<HOST>(memsize);
    proto_memcpy<SRC_MEM, HOST>(a_buf, tmp, memsize);
    for (int ii = 0; ii < memsize/sizeof(T); ii++)
    {
        pout() << tmp[ii] << ", ";
    }
    pout() << std::endl;
    */
    // END TEST CODE
}

template<typename T, unsigned int C, MemType SRC_MEM, MemType DST_MEM, Centering CTR>
void
LevelCopierOp<T, C, SRC_MEM, DST_MEM, CTR>::linearIn(void* a_buf, 
        const Box& a_bx,
        const DataIndex<BoxPartition>& a_index) const
{
    MemType mem = pointerMemType(a_buf);
    PROTO_ASSERT(mem == DST_MEM,
        "LevelCopierOp::linearIn | Error: Memtype mismatch. Detected %i but should be %i",
        mem, DST_MEM);
    auto& patch = m_dst->operator[](a_index);
    CInterval cint(0,C-1);
    patch.linearIn(a_buf, a_bx, cint);
    // BEGIN TEST CODE
    /*
    pout() << "LevelCopierOp::linearIn " << std::endl;
    int memsize = a_bx.size()*C*sizeof(T);
    T* tmp = (T*)proto_malloc<HOST>(memsize);
    proto_memcpy<DST_MEM, HOST>(a_buf, tmp, memsize);
    for (int ii = 0; ii < memsize/sizeof(T); ii++)
    {
        pout() << tmp[ii] << ", ";
    }
    pout() << std::endl;
    */
    // END TEST CODE
}

template<typename T, unsigned int C, MemType SRC_MEM, MemType DST_MEM, Centering CTR>
void
LevelCopierOp<T, C, SRC_MEM, DST_MEM, CTR>::localCopy(
        const Box& a_domain, const DataIndex<BoxPartition>& a_domainIndex,
        const Box& a_range,  const DataIndex<BoxPartition>& a_rangeIndex) const
{
    int LinearizationTest = 0;
    if (LinearizationTest == 0)
    {
        PR_TIMERS("LOP_localCopy_1");
        auto& src = m_src->operator[](a_domainIndex);
        auto& dst = m_dst->operator[](a_rangeIndex);
        Point shift = a_range.low() - a_domain.low();
        src.copyTo(dst, a_domain, shift);
    }
    else
    {
        PR_TIMERS("LOP_localCopy_2");
        int sizeSrc = linearSize(a_domain, a_domainIndex);
        int sizeDst = linearSize(a_range, a_rangeIndex);

        if (sizeSrc != sizeDst)
        {
            MayDay<void>::Abort("LinearizationTest failure: dest and source have different sizes");
        }
        
        void* srcBuffer = proto_malloc<SRC_MEM>(sizeSrc);
        void* dstBuffer = proto_malloc<DST_MEM>(sizeDst);
        {
            PR_TIMERS("linearIn / linearOut inside LDOP::op");
            linearOut(srcBuffer, a_domain, a_domainIndex);
            proto_memcpy<SRC_MEM, DST_MEM>(srcBuffer, dstBuffer, sizeSrc);
            if (SRC_MEM != DST_MEM)
            {
#ifdef PROTO_CUDA
                //cudaDeviceSynchronize();
	        protoThreadSynchronize();
#endif
            }
            linearIn (dstBuffer, a_range,  a_rangeIndex);
        }
    }
}

// =======================================================================
// LEVEL COPIERS
template<typename T, unsigned int C, MemType SRC_MEM, MemType DST_MEM, Centering CTR>
void
LevelCopier<T, C, SRC_MEM, DST_MEM, CTR>::buildMotionPlans(LevelCopierOp<T, C, SRC_MEM, DST_MEM, CTR>& a_op)
{
    PROTO_ASSERT(SRC_MEM == DST_MEM, 
        "LevelCopier::buildMotionPlans | Error: \
        Copy operations between HOST and DEVICE are not yet implemented.");
    // Construct Motion items for copyTo and insert into MotionPlans.
    // For each motion item, the source Box is specified relative to the location
    // of the source patch, and the dest Box relative to the location of the dest
    // patch. This is easy for the fromMotionPlan / localMotionPlan. In constructing
    // the toMotionPlan, though, must explicitly shift back and forth between local
    // coordinate systems defined by the source and destination.
    int destProcID = procID();

    DisjointBoxLayout srclayout = a_op.m_src->layout();
    DisjointBoxLayout destlayout = a_op.m_dst->layout();
    ProblemDomain pd = srclayout.domain();
    Point destGhost = a_op.m_dst->ghost();
    PR_assert(pd==destlayout.domain());
    // source pd and dest pd must be the same
    // For each Box on this processor, compute all the pairs of motionItems that move data
    // to this Box.
    for (auto diter : destlayout) 
    {
        Box bx = destlayout[diter].grow(destGhost);

        // Compute all the patches (as Points in a bitmap) that might move data to bx.

        Box dProbC = bx.coarsen(srclayout.boxSize());
        // Iterate over the possible source patches to contruct MotionItems.

        for (auto pt : dProbC)
        {
            // Find the possible source patch, or its periodic image.
            auto di = srclayout.find(pt);
            if (di != *srclayout.end())
            {
                Box destBx = Box(pt, pt).refine(srclayout.boxSize()) & bx;
                // Compute srcBx to be the shifted image of destBx. If destBx is in
                // the ProblemDomain box, srcBx == destBx.
                //Box srcBx = Box(pd.image(destBx.low()),pd.image(destBx.high()));

                if (!destBx.empty())
                {
                    // Create MotionItem.
                    // Compute srcBx to be the shifted image of destBx. If destBx 
                    // is in the ProblemDomain box, srcBx == destBx.
                    Box srcBx = Box(pd.image(destBx.low()),
                                    pd.image(destBx.high()));
                    MotionItem<BoxPartition> toMotionItem(di, diter, srcBx, destBx);
                    int srcProcID = srclayout.procID(di);
                    if (srcProcID == procID())
                    {
                        PROTO_ASSERT(srcProcID == destlayout.procID(diter),
                            "LevelCopier::buildMotionPlan | Error: Data corruption.");
                        this->m_localMotionPlan.push_back(toMotionItem);
                    }
                    else
                    {
                        PROTO_ASSERT(srcProcID != destlayout.procID(diter),
                            "LevelCopier::buildMotionPlan | Error: Data corruption.");
                        toMotionItem.procID = srcProcID;
                        this->m_toMotionPlan.push_back(toMotionItem);
                    }
                }
            }
        }
    }

    // Iterate over Boxes in the source to construct fromMotionPlan  
    for (auto siter : srclayout)
    {
        // We grow the source Box by dstGhost to capture all the patches that the
        // source will copy to. dProbC contains all the points in the bitmap for
        // the destination.
        Box srcBox = srclayout[siter];
        Box dProbC = srcBox.grow(destGhost).coarsen(destlayout.boxSize());
        for (auto pt : dProbC)
        {
            auto di = destlayout.find(pt);
            if (di != *destlayout.end())
            {
                Box bxDestShifted = Box(pt, pt).refine(destlayout.boxSize()).grow(destGhost);
                Box bxDest = destlayout[di].grow(destGhost);
                Point shiftSrcDest = bxDest.low() - bxDestShifted.low();
                bxDest = bxDest & (srcBox.shift(shiftSrcDest));
                Box bxSrc = bxDestShifted & srcBox;
                if (!bxSrc.empty())
                {
                    int destProcID = destlayout.procID(di);           
                    if (destProcID != procID())
                    {
                        MotionItem<BoxPartition> fromMotionItem(siter, di, bxSrc, bxDest);
                        fromMotionItem.procID =destProcID;
                        this->m_fromMotionPlan.push_back(fromMotionItem);
                    }
                }
            }
        }
    } // end src loop
    this->sort();
    /*
    pout() << std::endl << "=======================================================================" << std::endl;
    pout() << "Printing LevelCopier" << std::endl;
    pout() << "Source Layout" << std::endl;
    srclayout.print();
    pout() << "Dest Layout" << std::endl;
    destlayout.print();
    this->print();
    */
}

template<typename T, unsigned int C, MemType MEM, Centering CTR>
void
LevelExchangeCopier<T, C, MEM, CTR>::buildMotionPlans(LevelCopierOp<T, C, MEM, MEM, CTR>& a_op)
{
    auto layout = a_op.m_src->layout();
    auto ghost  = a_op.m_src->ghost();

    int myProcID = Proto::procID();
    for (auto dit : layout)
    {
        Box bx0 = layout[dit];
        Box bx = layout.domain() & bx0.grow(ghost);
        NeighborIterator nit(layout, bx);
        for (nit.begin(); nit.ok(); ++nit)
        {
            int nbrProcID = nit.procID();
            Point nbrShift = nit.shift();
            // For each pair of interacting boxes, we construct a
            // Motion to, from, and local.
            // Compute the boxes corresponding to filling from *dit to *nit
            // The source box is shifted, the dest box is not.
            
            Box destBx = nit.destBox();
            if (!destBx.empty() && !bx0.contains(destBx))
            {
                Box srcBx = nit.srcBox();
                MotionItem<BoxPartition> toMotionItem(*nit, dit, srcBx, destBx);
                if (nbrProcID == myProcID)
                {
                    // local motion
                    this->m_localMotionPlan.push_back(toMotionItem);
                } else {
                    toMotionItem.procID = nbrProcID;
                    this->m_toMotionPlan.push_back(toMotionItem);
                }
            }
            
            // Construct the motion for filling the ghost region for *nit from *dit

            if (nbrProcID != myProcID)
            {
                Box nbrbx = Box(nit.point(), nit.point()).refine(layout.boxSize()).grow(ghost);
                Box nbrbxShift = nbrbx.shift(-nbrShift);
                Box srcFromBx = nbrbxShift & bx0;
                Box destFromBx = srcFromBx.shift(nbrShift);
                MotionItem<BoxPartition> fromMotionItem(dit, *nit, srcFromBx, destFromBx);
                fromMotionItem.procID = nbrProcID;
                this->m_fromMotionPlan.push_back(fromMotionItem);
            }
        } // end neighbor loop
    } // end patch loop
    /*
    this->sort();
    pout() << std::endl << "=======================================================================" << std::endl;
    pout() << "Printing ExchangeCopier" << std::endl;
    pout() << "Layout" << std::endl;
    layout.print();
    pout() << "Copier" << std::endl;
    this->print();
    */
}

template<typename T, unsigned int C, MemType MEM, Centering CTR>
void interpBoundaries(
        LevelBoxData<T, C, MEM, CTR>& a_crse,
        LevelBoxData<T, C, MEM, CTR>& a_fine,
        InterpStencil<T>&       a_interp)
{
    //FIXME: This implementation doesn't address refined regions directly
    //       adjacent to domain boundaries. -CLG
    Point refRatio = a_interp.ratio();
    auto  cfLayout = a_fine.layout().coarsen(refRatio);
    Point cfGhost = a_fine.ghost() / refRatio + Point::Ones();
    Point interpGhost = a_interp.ghost();
    LevelBoxData<T, C, MEM, CTR> temp(cfLayout, cfGhost + interpGhost);
    interpBoundaries(a_crse, a_fine, temp, a_interp);
}

template<typename T, unsigned int C, MemType MEM, Centering CTR>
void interpBoundaries(
        LevelBoxData<T, C, MEM, CTR>& a_crse,
        LevelBoxData<T, C, MEM, CTR>& a_fine,
        LevelBoxData<T, C, MEM, CTR>& a_crseFine,
        InterpStencil<T>&       a_interp)
{
    a_crse.exchange();
    a_crse.copyTo(a_crseFine);
    const auto& fineLayout = a_fine.layout();
    for (auto iter : fineLayout)
    {
        if (!fineLayout.onLevelBoundary(fineLayout.point(iter))){ continue; }
        auto& fine_i = a_fine[iter];
        auto& crse_i = a_crseFine[iter];
        BoxData<T, C, MEM> fine_0(fineLayout[iter]);
        fine_i.copyTo(fine_0);  // save a copy of non-interior data
        fine_i |= a_interp(crse_i);
        fine_0.copyTo(fine_i);
        // This method seems like a shortcut / hack, but each interpolation call
        // results in refRatio^DIM stencil applies, each of which is a kernel
        // launch for GPU targets. I suspect that two BoxData::copyTo calls
        // are faster, but this should be tested. -CLG
    }
    a_fine.exchange();
}

template<typename T, unsigned int C, MemType MEM, Centering CTR>
void averageDown(
        LevelBoxData<T, C, MEM, CTR>& a_crse,
        LevelBoxData<T, C, MEM, CTR>& a_fine,
        Point                         a_refRatio)
{
    auto  cfLayout = a_fine.layout().coarsen(a_refRatio);
    LevelBoxData<T, C, MEM, CTR> temp(cfLayout, Point::Zeros());
    averageDown(a_crse, a_fine, temp, a_refRatio);
}

template<typename T, unsigned int C, MemType MEM, Centering CTR>
void averageDown(
        LevelBoxData<T, C, MEM, CTR>& a_crse,
        LevelBoxData<T, C, MEM, CTR>& a_fine,
        int                           a_refRatio)
{
    averageDown(a_crse, a_fine, Point::Ones(a_refRatio));
}

template<typename T, unsigned int C, MemType MEM, Centering CTR>
void averageDown(
        LevelBoxData<T, C, MEM, CTR>& a_crse,
        LevelBoxData<T, C, MEM, CTR>& a_fine,
        LevelBoxData<T, C, MEM, CTR>& a_crseFine,
        Point                         a_refRatio)
{
    a_fine.exchange();
    auto AVG = Stencil<T>::AvgDown(a_refRatio);
    for (auto iter : a_fine.layout())
    {
        auto& fine_i = a_fine[iter];
        auto& crse_i = a_crseFine[iter];
        crse_i |= AVG(fine_i);
    }
    a_crseFine.copyTo(a_crse);
}

template<typename T, unsigned int C, MemType MEM, Centering CTR>
void averageDown(
        LevelBoxData<T, C, MEM, CTR>& a_crse,
        LevelBoxData<T, C, MEM, CTR>& a_fine,
        LevelBoxData<T, C, MEM, CTR>& a_crseFine,
        int                           a_refRatio)
{
    averageDown(a_crse, a_fine, a_crseFine, Point::Ones(a_refRatio));
}

template<typename T, unsigned int C, MemType MEM, Centering CTR>
template<typename _T, unsigned int _C, MemType _MEM, Centering _CTR>
bool LevelBoxData<T, C, MEM, CTR>::compatible(const LevelBoxData<_T, _C, _MEM, _CTR>& a_data) const
{
    return m_layout.compatible(a_data.layout());
}

template<typename T, unsigned int C, MemType MEM, Centering CTR>
bool LevelBoxData<T, C, MEM, CTR>::compatible(const DisjointBoxLayout& a_layout) const
{
    return m_layout.compatible(a_layout);
}

template<typename T, unsigned int C, MemType MEM, Centering CTR>
bool LevelBoxData<T, C, MEM, CTR>::compatible(const DataIterator<BoxPartition>& a_iter) const
{
    return m_layout.compatible(a_iter);
}

template<typename T, unsigned int C, MemType MEM, Centering CTR>
bool LevelBoxData<T, C, MEM, CTR>::compatible(const DataIndex<BoxPartition>& a_index) const
{
    return m_layout.compatible(a_index);
}
