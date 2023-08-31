template<typename T, unsigned int C, MemType MEM>
FluxRegisterCopierOp<T,C,MEM>::FluxRegisterCopierOp(LevelFluxRegister<T, C, MEM>& a_register)
{
    m_register = &a_register;
}

template<typename T, unsigned int C, MemType MEM>
int FluxRegisterCopierOp<T,C,MEM>::linearSize(
                                   const Box& a_bx,
                                   const DataIndex<BoxPartition>& a_index) const
{
    return C*a_bx.size()*sizeof(T);
}
template<typename T, unsigned int C, MemType MEM>
void FluxRegisterCopierOp<T,C,MEM>::linearOut(
                                 void* a_buf,
                                 const Box& a_bx,
                                 const DataIndex<BoxPartition>& a_index) const
{
    // linearOut copies from source data holders into buffers.
    BoxData<T,C,MEM>& data = m_register->sourceData(a_bx, a_index);
    CInterval cint(0,C-1);
    data.linearOut(a_buf, a_bx, cint);
}

template<typename T, unsigned int C, MemType MEM>
void FluxRegisterCopierOp<T,C,MEM>::linearIn(
                                 void* a_buf,
                                 const Box& a_bx,
                                 const DataIndex<BoxPartition>& a_index) const
{
    // linearIn copies from buffers into destination data holders.
    BoxData<T,C,MEM>& data = m_register->destData(a_bx,a_index);
    CInterval cint(0,C-1);
    data.linearIn(a_buf, a_bx,cint);
}

template<typename T, unsigned int C, MemType MEM>
void FluxRegisterCopierOp<T,C,MEM>::localCopy(
        const Box& a_domain,
        const DataIndex<BoxPartition>& a_domainIndex,
        const Box& a_range,
        const DataIndex<BoxPartition>& a_rangeIndex) const
{
#define noLinearizationTest 1
  //PROTO_ASSERT(a_domain == a_range, 
  //     "FluxRegisterCopierOp::localCopy | Error: domain and range should be the same.");
#if noLinearizationTest

    PR_TIMERS("LOP_localCopy_1");
    BoxData<T,C,MEM>& src  = m_register->sourceData(a_domain,a_domainIndex);
    BoxData<T,C,MEM>& dest = m_register->destData(a_range,a_rangeIndex);
    Point shift = a_range.low() - a_domain.low();
    src.copyTo(dest, a_domain, shift);
#else
    // code for testing linearIn / linearOut while running serial.
    PR_TIMERS("LOP_localCopy_2");
    int srcSize = linearSize(a_domain, a_domainIndex);
    int dstSize = linearSize(a_range,  a_rangeIndex);

    if (srcSize != dstSize)
    {
        MayDay<void>::Abort("LinearizationTest failure: dest and source have different sizes");
    }

    vector<char> buffer(sizeSource);
    {
        PR_TIMERS("linearIn / linearOut inside LDOP::op");
        void* b = (void*)&(buffer[0]);
        linearOut(b, a_domain, a_domainIndex);
        linearIn (b, a_range,  a_rangeIndex);
    }
#endif
}

template<typename T, unsigned int C, MemType MEM>
void FluxRegisterCopier<T,C,MEM>::buildMotionPlans (FluxRegisterCopierOp<T, C, MEM >& a_op)
{
    auto & lfr = *(a_op.m_register);
    auto & refRatio = lfr.m_refRatio;
    auto & fineLayout = lfr.m_fineLayout;
    auto & crseLayout = lfr.m_crseLayout;
    auto & crseRegisters = lfr.m_crseRegisters;
    auto & crseIncrement = lfr.m_crseIncrement;
    auto & fineRegisters = lfr.m_fineRegisters;
    auto & fineIndices = lfr.m_fineIndices;

    //std::cout << "Building Destination Motion Items" << std::endl;
    for (auto citer = crseLayout.begin(); citer.ok(); ++citer)
    {
        auto& localRegisters = crseRegisters[*citer];
        for (int k = 0; k < localRegisters.size(); k++)
        {
            // Find the index of the fineRegister source.
            // The index is associated with the fine patch contributing the averaged flux
            auto& reg     = localRegisters[k];
            Box   destBox = reg.m_data->box();
            auto  side    = reg.m_side;
            int   dir     = reg.m_dir;
            Point finePatchPoint = destBox.low()*refRatio/fineLayout.boxSize() + Point::Basis(dir,side);
            //std::cout << "\tCoarse Box: " << destBox << " | Dir: " << dir << " | Side: " << side << std::endl;
            //std::cout << "\tAdjacent Fine Patch: " << finePatchPoint << std::endl;
            auto fineIndex = fineLayout.find(finePatchPoint);
            
            PROTO_ASSERT(fineIndex != *fineLayout.end(),
                "FluxRegisterCopier::buildMotionPlans | Error: Could not find fine register index.");
            // Compute fine register box (nontrivial b/c of possibility of periodic BCs. 
            Box cfBox  = fineLayout.box(fineIndex).coarsen(refRatio);
            Box srcBox = cfBox.adjacent(dir, flip(side), 1);
            // Boxes for the fine registers are defined relative to the fine patch
            // coordinates, even if the register Boxes are periodic images outside the
            // ProblemDomain Box - PC.
            MotionItem<BoxPartition> toMotionItem(fineIndex, *citer, srcBox, destBox);
            auto fineID = fineLayout.procID(fineIndex);
            if (fineID == procID())
            { 
                this->m_localMotionPlan.push_back(toMotionItem);
            }
            else
            {
                toMotionItem.procID = fineID;
                this->m_toMotionPlan.push_back(toMotionItem);
            }
        } // for each local coarse register
    } // end for coarse patch

    // Iterate over fine grids that I own to obtain entries in m_fromMotionPlan.           
    int destProcID = procID();
    // std::cout << "Building Source Motion Items" << std::endl;
    for (auto fiter = fineLayout.begin(); fiter.ok(); ++fiter)
    {
      //Box cfBox = fiter.box().coarsen(refRatio);
        for (int dir = 0; dir < DIM; dir++)
        {
            SideIterator siter;
            for (siter.begin(); siter.ok(); ++siter)
            {
                int k = fineIndices[*fiter][dir][*siter];
                if (k >= 0) //otherwise there is no fine register for this (dir, side)
                {
                  Box srcBox = fineRegisters[*fiter][k].m_data->box();
                  Point crsePatchPoint = srcBox.low()/crseLayout.boxSize();
                    // std::cout << "\tCoarse Box: " << srcBox << " | Dir: " << dir << " | Side: " << *siter << std::endl;
                    
                    auto crseIndex  = crseLayout.find(crsePatchPoint);
                    int  crseProcID = crseLayout.procID(crseIndex);
                    //Automatically defines dstBox to be the correct periodic image. NO IT DOESN'T.
                    Box dstBox = crseLayout.domain() & srcBox;
                    dstBox = crseLayout.domain().image(dstBox); 
                    if (crseProcID != procID())
                    {
                        MotionItem<BoxPartition> fromMotionItem(*fiter, crseIndex, srcBox, dstBox);
                        fromMotionItem.procID = crseProcID;
                        this->m_fromMotionPlan.push_back(fromMotionItem);
                    }
                }
            }
        }
    }     
    this->sort();
    /*
    pout() << std::endl << "=======================================================================" << std::endl;
    pout() << "Printing FluxRegisterCopier" << std::endl;
    pout() << "Coarse Layout" << std::endl;
    crseLayout.print();
    pout() << "Fine Layout" << std::endl;
    fineLayout.print();
    this->print();
    */
}

template<typename T, unsigned int C, MemType MEM>
LevelFluxRegister<T,C,MEM>::LevelFluxRegister(
        const DisjointBoxLayout& a_fineLayout,
        const DisjointBoxLayout& a_crseLayout,
        const Point& a_refRatio,
        std::array<T,DIM> a_dxCrse)
{
    define(a_fineLayout,a_crseLayout,a_refRatio,a_dxCrse);
}

template<typename T, unsigned int C, MemType MEM>
void LevelFluxRegister<T,C,MEM>::define (
        const DisjointBoxLayout& a_crseLayout,
        const DisjointBoxLayout& a_fineLayout,
        const Point& a_refRatio,
        std::array<T, DIM> a_dxCrse)
{
    // Fine registers are defined relative to the fine patch coordinates, and coarse 
    // registers are defined relative to the coarse patch coordinates. What happens in 
    // the presence of periodic boundaries is mediated in LevelFluxRegister<T,C,MEM>Operator.
    
    m_refRatio = a_refRatio;
    m_dxCrse = a_dxCrse;
    m_fineLayout = a_fineLayout;
    m_crseLayout = a_crseLayout;

    // Check to see whether Coarse patch is exactly covered by fine patches.
    Point patchRefRatio = m_crseLayout.boxSize()*m_refRatio / m_fineLayout.boxSize();
    PROTO_ASSERT(patchRefRatio*m_fineLayout.boxSize() == m_crseLayout.boxSize()*m_refRatio,
        "LevelFluxRegister::define | Error: coarse region is not covered exactly by fine region.");

    // Check to see whether coarsened fine patch is more than one coarse cell wide in all directions that are being refined. 
    for (int dir = 0; dir < DIM; dir++)
    {
        if (m_refRatio[dir] > 1)
        {
            PROTO_ASSERT(m_fineLayout.boxSize()[dir]/m_refRatio[dir] > 1,
                "LevelFluxRegister::define | Error: fine patch cell size is smaller than the refinement ratio.");
        }
    }
    m_crseRegisters.clear();
    m_crseIncrement.clear();
    m_crseIndices.clear();
    m_fineRegisters.clear();
    m_fineIndices.clear();

    int numCrsePatches = m_crseLayout.localSize();
    //std::cout << "FluxRegister::Define | Building sets of registers for " << numCrsePatches << " coarse patches." << std::endl;
    m_crseRegisters.resize(numCrsePatches);
    m_crseIncrement.resize(numCrsePatches);
    m_crseIndices.resize(numCrsePatches);
    
    // Set up coarse data holders.  
    // std::cout << "Building Coarse Data Holders" << std::endl;
    for (auto citer = m_crseLayout.begin(); citer.ok(); ++citer)
    {
        Box crseBox = m_crseLayout[*citer];
        m_crseRegisters[*citer].clear();
        m_crseIncrement[*citer].clear();
        m_crseIndices[*citer].clear();
        
        //setup: compute box corresponding to possible fine patches covered by coarse patch.
        Box finePatchPoints = crseBox.refine(m_refRatio).coarsen(m_fineLayout.boxSize());
        // Iterate over the fine patches contained in the coarse box, and see whether patches 
        // adjacent to them in each direction can contribute to a flux register in the coarse 
        // patch. 
        
        //std::cout << "\tCoarse Patch: " << crseBox << " | Fine Patch Points: " << finePatchPoints << std::endl;
        int k = 0;
        for (auto biter = finePatchPoints.begin(); biter.ok(); ++biter)
        {
            Point pt = *biter;
            auto fineIndex = m_fineLayout.find(pt);
            // only check fine patches NOT in part of the fine layout
            if (fineIndex != *m_fineLayout.end()) { continue; } 
            Box cfBox = Box(pt, pt).refine(m_fineLayout.boxSize()).coarsen(m_refRatio);
            //std::cout << "\t\tFound External Patch: " << pt << " | Coarsened Box: " << cfBox << std::endl;
            for (int dir = 0; dir < DIM; dir++)
            {
                SideIterator siter;
                for (siter.begin(); siter.ok(); ++siter)
                {
                    Point ptAdj = *biter + Point::Basis(dir,*siter);
                    //std::cout << "\t\t\tChecking Point: " << ptAdj << std::endl;
                    // find accounts for periodic boundaries
                    auto adjIndex = m_fineLayout.find(ptAdj);
                    if (adjIndex != *m_fineLayout.end())
                    {
                        // fineIndex is the index of a fine patch NOT in the fine layout
                        // adjIndex is the index of a fine patch IN the fine layout
                        // hence, the interface between the two is a coarse-fine boundary
                        Box destBox = cfBox.face(dir, *siter);
                        // Intersect with the Problem domain to shift the Box to its proper 
                        // periodic image in the coarse patch.
                        destBox = m_crseLayout.domain() & destBox; 
                        
                        //std::cout << "\t\t\tFound Boundary Patch: " << ptAdj << " | Shifted Coarsened Box: " << destBox << std::endl;

                        // Create register and insert it into the coarse register data structure.
                        // the register is not uniquely determined by it's (dir, side) 
                        // the key is generated using destBox which incodes (dir, side, tile position)
                        int thisKey = key(destBox, *citer);
                        shared_ptr<BoxData<T,C,MEM> > temp1(new BoxData<T,C,MEM>(destBox));
                        shared_ptr<BoxData<T,C,MEM> > temp2(new BoxData<T,C,MEM>(destBox));
                        Register<T,C,MEM> reg(temp1, dir, *siter);
                        m_crseRegisters[*citer].push_back(reg);
                        m_crseIncrement[*citer].push_back(temp2);
                        m_crseIndices  [*citer][thisKey] = k;
                        k++;
                    }
                } // end for side
            } // end for DIM
        } // end for fine in coarse
    }  // end for coarse

    // Set up fine data holders.
    //std::cout << "LevelFluxReguster::define | Building Fine Registers" << std::endl;
    int numFinePatches = m_fineLayout.localSize();
    //std::cout << "FluxRegister::Define | Building sets of registers for " << numFinePatches << " coarse patches." << std::endl;
    m_fineRegisters.resize(numFinePatches);
    m_fineIndices.resize(numFinePatches);
    for (auto fiter = m_fineLayout.begin(); fiter.ok(); ++fiter)
    {
        Point finePatchPoint = m_fineLayout.point(*fiter);
        Box fineBox = m_fineLayout[*fiter];
        Box cfBox = fineBox.coarsen(m_refRatio);
        //std::cout << "\tFine Patch Point: " << finePatchPoint << " | Coarse-Fine Box: " << cfBox << std::endl;
        int k = 0;
        for (int d = 0; d < DIM; d++)
        {
            SideIterator siter;
            for (siter.begin(); siter.ok(); ++siter)
            {
                Point srcPt = finePatchPoint + Point::Basis(d, *siter);
                auto adjIndex = m_fineLayout.find(srcPt);
                //std::cout << "\t\tChecking Adjacent Point: " << srcPt << std::endl;
                // look for adjacent patches NOT in the fine layout
                if (adjIndex == *m_fineLayout.end())
                {
                    // Fine registers are AVERAGED fine fluxes
                    // hence, domains are coarse
                    Box srcBox = cfBox.adjacent(d, *siter, 1);
                    // Deal with possible periodic images
                    srcBox = m_crseLayout.domain() & srcBox;
                    shared_ptr<BoxData<T,C,MEM> > temp(new BoxData<T,C,MEM>(srcBox));
                    Register<T,C,MEM> reg(temp, d, *siter);
                    m_fineRegisters[*fiter].push_back(reg);
                    m_fineIndices[*fiter][d][(int)(*siter)] = k;
                    //std::cout << "\t\t\tFound CF Boundary. Register Box: " << srcBox << " | Index: " << k << std::endl;
                    k++;
                } else {
                    m_fineIndices[*fiter][d][(int)(*siter)] = -1;
                }
            } // end for side
        } // end for DIM
    } // end for fine

    // We now have enough information to build the motion plans for the copier.
    FluxRegisterCopierOp<T,C,MEM> op(*this);
    m_copier.define(op);
    // m_copier.buildMotionPlans(op);

    // Create stencils for averaging down fluxes.
    for (int d = 0; d < DIM; d++)
    {   
        SideIterator siter;
        for (siter.begin(); siter.ok(); ++siter)
        {
            // Side::Lo is not a typo -CLG
            auto avg = Stencil<T>::AvgDownFace(d, Side::Lo, m_refRatio);
            if (*siter == Side::Hi)
            {
                avg = avg*(-1);
            } else if (*siter == Side::Lo)
            {
                // Side::Lo is not a typo here either -CLG
                avg.destShift() = Point::Basis(d, Side::Lo);
            }
            m_avgFlux[d][*siter] = avg;
        }
    }

    reset();
    //std::cout << "Built Flux Register: " << std::endl;
    //std::cout << "Coarse Registers: " << std::endl;
    for (auto iter = m_crseLayout.begin(); iter.ok(); ++iter)
    {
      //std::cout << "\tCoarse Box: " << m_crseLayout[*iter] << std::endl;
        int numRegisters = m_crseRegisters[*iter].size();
        //std::cout << "\tNumber of Registers: " << numRegisters << std::endl;
        for (int ii = 0; ii < numRegisters; ii++)
        {
            auto& reg = m_crseRegisters[*iter][ii];
            //std::cout << "\t\tBox: " << reg.m_data->box() << " | Dir: " << reg.m_dir << " | Side: " << reg.m_side << std::endl;
        }
    }
    // std::cout << "Fine Registers: " << std::endl;
    for (auto iter = m_fineLayout.begin(); iter.ok(); ++iter)
    {
      //std::cout << "\tFine Box: " << m_fineLayout[*iter] << std::endl;
        int numRegisters = m_fineRegisters[*iter].size();
        //std::cout << "\tNumber of Registers: " << numRegisters << std::endl;
        for (int ii = 0; ii < m_fineRegisters[*iter].size(); ii++)
        {
            auto& reg = m_fineRegisters[*iter][ii];
            //std::cout << "\t\tBox: " << reg.m_data->box() << " | Dir: " << reg.m_dir << " | Side: " << reg.m_side << std::endl;
        }
    }
}

template<typename T, unsigned int C, MemType MEM>
void LevelFluxRegister<T,C,MEM>::incrementCoarse (
        const BoxData<T,C,MEM>& a_flux,
        const DataIndex<BoxPartition>& a_crseIndex,
        const T& a_weight,
        unsigned int a_dir)
{
    //std::cout << "FluxRegister::incrementCoarse | patch index: " << (int)a_crseIndex << " | box: " << a_flux.box() << std::endl;
    auto& crseRegisters = m_crseRegisters[a_crseIndex];
    T weight = a_weight / m_dxCrse[a_dir];
    Stencil<T> scaleHigh = weight*Shift::Basis(a_dir, +1);
    Stencil<T> scaleLow = -weight*Shift::Zeros();
    for (int k = 0; k < crseRegisters.size(); k++)
    {   
        auto& reg = crseRegisters[k];
        if (reg.m_dir == a_dir)
        {
            auto& data = reg.m_data; 
            auto  side = reg.m_side;
            if (side == Side::Lo)
            {
                (*data) += scaleLow(a_flux);
            }
            else if (side == Side::Hi)
            {
                (*data) += scaleHigh(a_flux);
            }
        }
    }
}
template<typename T, unsigned int C, MemType MEM>
void LevelFluxRegister<T,C,MEM>::incrementFine (
        const BoxData<T, C , MEM>& a_flux,
        const DataIndex<BoxPartition>& a_fineIndex,
        const T& a_weight,
        unsigned int a_dir)
{
    auto & fineRegister = m_fineRegisters[a_fineIndex];
    T weight = a_weight / m_dxCrse[a_dir];
    for (int k = 0; k < fineRegister.size(); k++)
    {
        auto& reg = fineRegister[k];
        if (reg.m_dir == a_dir)
        {
            auto& data = reg.m_data;
            auto  side = reg.m_side;
            (*data) += m_avgFlux[a_dir][side](a_flux, weight);
        }
    }
}

template<typename T, unsigned int C, MemType MEM>
void LevelFluxRegister<T,C,MEM>::reflux(
        LevelBoxData<T,C,MEM>& a_coarseData,
        const T& a_weight)
{
    PROTO_ASSERT(m_crseLayout.compatible(a_coarseData.layout()),
            "LevelFluxRegister::reflux | Error: \
            Input data has incompatible layout.");
    
    m_copier.execute();

    // fine register data has been copied to m_crseIncrement.
    // We compute the difference between the coarse and fine registers, and increment a_coarseData.
    Stencil<T> increment = (-1.0)*Shift::Zeros();
    for (auto iter = a_coarseData.begin(); iter.ok(); ++iter)
    {
        int numRegisters = m_crseIncrement[*iter].size();
        auto& dest = m_crseIncrement[*iter];
        auto& creg = m_crseRegisters[*iter];
        // subtract off the coarse fluxes
        for (int k = 0; k < numRegisters; k++)
        {
            *(dest[k]) += increment(*(creg[k].m_data));
            *(dest[k]) *= a_weight;
        }
        // add increment to input data
        for (int k = 0; k < numRegisters; k++)
        {
            a_coarseData[*iter] += *(dest[k]);
        }
    }          
}

template<typename T, unsigned int C, MemType MEM>
void LevelFluxRegister<T,C,MEM>::reset(T a_val)
{
    // set values in the registers to zero.
    for (auto citer = m_crseLayout.begin(); citer.ok(); ++citer)
    {
        for (int k = 0; k < m_crseIndices[*citer].size(); k++)
        {
            m_crseRegisters[*citer][k].m_data->setVal(a_val);
            m_crseIncrement[*citer][k]->setToZero();
        }
    }
    for (auto fiter = m_fineLayout.begin(); fiter.ok(); ++fiter)
    {
        for (int k = 0; k < m_fineRegisters[*fiter].size();k++)
        {
            m_fineRegisters[*fiter][k].m_data->setToZero();
        }
    }
}

// a_bx is the coarse box updated by refluxing
template<typename T, unsigned int C, MemType MEM>
int LevelFluxRegister<T,C,MEM>::key(const Box& a_bx, const DataIndex<BoxPartition>& a_di)
{
    // Find which direction is normal to the register. We use the fact that the register
    // is of size 1 only in the normal direction.
   
    Box crseBox = m_crseLayout[a_di];
    Point origin = crseBox.low();
    Box tileBox = a_bx.shift(-origin);
    Point sizes = crseBox.high() - crseBox.low() + Point::Ones();

    int k = tileBox.low()[0];
    int p = 1;
    for (int ii = 1; ii < DIM; ii++)
    {
        p *= sizes[ii-1];
        k += p*tileBox.low()[ii];
    }
    p *= sizes[DIM-1];
    k += p*tileBox.high()[0];
    for (int ii = 1; ii < DIM; ii++)
    {
        p *= sizes[ii-1];
        k += p*tileBox.high()[ii];
    }
    
    return k;

    /* 
    int dir = -1;
    int side = -1;
    Box crseBox = m_crseLayout[a_di];
    for (int d = 0; d < DIM; d++)
    {
        SideIterator siter;
        for (siter.begin(); siter.ok(); ++siter)
        {
            Box faceBox = crseBox.face(d, *siter, 1);
            if (faceBox.contains(a_bx))
            {
                Point tileSize = a_bx.high() - a_bx.low() + Point::Ones() - Point::Basis(d);
                Point shiftPoint = faceBox.low();
                Box tile = a_bx.shift(-shiftPoint).coarsen(tileSize);
                Box tileDomain = faceBox.shift(-shiftPoint).coarsen(tileSize);
                PROTO_ASSERT(tileDomain.contains(tile), "LevelFluxRegister::key | Error: Could not compute tile index");
                PROTO_ASSERT(tileDomain.contains(tile), "LevelFluxRegister::key | Error: Invalid tile box");
                int numTiles = tileDomain.size();
                int tileIndex = tileDomain.index(tile.low());

                return tileIndex + numTiles*d + numTiles*DIM*(*siter);
            }
        }
    }
    MayDay<void>::Abort("LevelFluxRegister::key | Error: could not generate key.");
    */
}

template<typename T, unsigned int C, MemType MEM>
BoxData<T,C,MEM>& LevelFluxRegister<T,C,MEM>::sourceData(const Box& a_bx, const DataIndex<BoxPartition>& a_fineIndex)
{
    Box fineBox = m_fineLayout.box(a_fineIndex);
    Box cfBox = fineBox.coarsen(m_refRatio);
    SideIterator siter;
    for (siter.begin(); siter.ok(); ++siter)
    {
        for (int dir = 0; dir < DIM; dir++)
        {
            if (a_bx == cfBox.adjacent(dir, *siter, 1))
            {
                int k = m_fineIndices[a_fineIndex][dir][(int)(*siter)];
                return *(m_fineRegisters[a_fineIndex][k].m_data);
            }
        }
    }
    std::cout << "Invalid Box: " << a_bx << " or invalid index: " << (int)a_fineIndex << std::endl;
    PROTO_ASSERT(false,
        "LevelFluxRegister::sourceData | Error: Could not locate data.");
}

template<typename T, unsigned int C, MemType MEM>
BoxData<T,C,MEM>& LevelFluxRegister<T,C,MEM>::destData(const Box& a_bx, const DataIndex<BoxPartition>& a_crseIndex)
{
    int crseKey       = key(a_bx, a_crseIndex);
    int mapIndex      = m_crseIndices[a_crseIndex][crseKey];
    return *(m_crseIncrement[a_crseIndex][mapIndex]);
}

template<typename T, unsigned int C, MemType MEM>
void LevelFluxRegister<T,C,MEM>::print() const
{
    using Proto::pout;
    pout() << "Level Flux Register: " << std::endl;
    pout() << "Coarse Registers: " << std::endl;
    for (int ii = 0; ii < m_crseRegisters.size(); ii++)
    {
        for (int jj = 0; jj < m_crseRegisters[ii].size(); jj++)
        {
            auto item = m_crseRegisters[ii][jj];
            pout() << "\t box: " << item.m_data->box();
            pout() << " | dir: " << item.m_dir;
            pout() << " | side: " << item.m_side << std::endl;
        }
    }
    pout() << "Fine Registers: " << std::endl;
    for (int ii = 0; ii < m_fineRegisters.size(); ii++)
    {
        for (int jj = 0; jj < m_fineRegisters[ii].size(); jj++)
        {
            auto item = m_fineRegisters[ii][jj];
            pout() << "\t box: " << item.m_data->box();
            pout() << " | dir: " << item.m_dir;
            pout() << " | side: " << item.m_side << std::endl;
        }
    }
}
