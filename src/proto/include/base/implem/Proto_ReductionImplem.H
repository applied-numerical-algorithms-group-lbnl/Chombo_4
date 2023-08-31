//// CONSTRUCTOR
template<typename T, Operation OP, MemType MEM>
Reduction<T,OP,MEM>::Reduction(bool dynamic) {
    m_hostTotal = (T*)proto_malloc<HOST>(sizeof(T));
    m_hostTemp = (T*)proto_malloc<HOST>(sizeof(T));
    *m_hostTotal = init();
    *m_hostTemp = init();
    if (MEM == DEVICE)
    {
#ifdef PROTO_CUDA
        // TODO: get number of streams in runtime
        //m_dynamic = dynamic;
        m_dynamic = dynamic;
        protoDeviceProp prop;
        protoGetDeviceProperties(&prop, 0); // assumes system has identical GPUs
        m_numThreads = prop.maxThreadsPerBlock; // threads in a block
        m_warpSize = prop.warpSize; // threads in a warp
        m_deviTotal   = (T*)proto_malloc<MEM>(sizeof(T));
        initKernel<T,OP><<<1,1>>>(m_deviTotal);
        if (!m_dynamic) {
            // filling device with blocks for this kernel. TODO: divide by streams.
            m_numBlocks = (prop.maxThreadsPerMultiProcessor/m_numThreads)*prop.multiProcessorCount;
            m_deviTemp = (T*)proto_malloc<DEVICE>(m_numBlocks*sizeof(T)); //TODO: mult by streams
        }
#else
#endif
    } else {
        // Nothing else to do for HOST
    }
}

/// DESTRUCTOR
template<typename T, Operation OP, MemType MEM>
Reduction<T,OP,MEM>::~Reduction() {
    proto_free<HOST>(m_hostTotal);
    proto_free<HOST>(m_hostTemp);
    if (MEM == DEVICE)
    {
#ifdef PROTO_CUDA
        proto_free<MEM>(m_deviTotal);
        proto_free<DEVICE>(m_deviTemp);
#endif
    }
}

// RESET
template<typename T, Operation OP, MemType MEM>
void Reduction<T,OP,MEM>::reset() {
    if (MEM == DEVICE)
    {
#ifdef PROTO_CUDA
        initKernel<T,OP><<<1,1>>>(m_deviTotal);
#else
#endif
    } else {
        *m_hostTotal = init();
    }
}

/// INITIALIZE VALUE
template<typename T, Operation OP, MemType MEM>
CUDA_DECORATION
T Reduction<T,OP,MEM>::init()
{
    switch (OP) {
        case Max:
            return std::numeric_limits<T>::min();
        case Min:
            return std::numeric_limits<T>::max();
        default: // Abs, Sum, SumAbs
            return static_cast<T>(0);
    }
}

/// UPDATE VALUE
template<typename T, Operation OP, MemType MEM>
CUDA_DECORATION
void Reduction<T,OP,MEM>::update(T& a_v1, const T a_v2)
{
    switch(OP) {
        case Max:
            a_v1 = max<T>(a_v1,a_v2);
            break;
        case Min:
            a_v1 = min<T>(a_v1,a_v2);
            break;
        case Abs:
            a_v1 = (a_v2 > 0 ? max<T>(a_v1,a_v2) : max<T>(a_v1,-a_v2));
            break;
        case Sum:
            a_v1 += a_v2;
            break;
        case SumAbs:
            a_v1 += (a_v2 > 0 ? a_v2 : -a_v2);
            break;
    }
}

/// FETCH
// TODO: on-DEVICE fetch?
template<typename T, Operation OP, MemType MEM>
T Reduction<T,OP,MEM>::fetch() {
    if (MEM == DEVICE)
    {
#ifdef PROTO_CUDA
        //cudaDeviceSynchronize(); // make sure kernels are finished
	protoThreadSynchronize();
        proto_memcpy<MEM, HOST>(m_deviTotal, m_hostTemp, sizeof(T));
       // cudaDeviceSynchronize(); // make sure copy to HOST is done (unnecessary?)
	protoThreadSynchronize();
        update(*m_hostTotal, *m_hostTemp);
        // TODO: synchronize multiple streams
#endif
    } else {
        // nothing to do for HOST
    }
#ifdef PR_MPI
    T global;
    auto datatype = mpiDatatype<T>();
    barrier();
    switch (OP)
    {
        case Abs:
        case Max:
            MPI_Allreduce(m_hostTotal, &global, 1, datatype, MPI_MAX, MPI_COMM_WORLD);
            break;
        case Min:
            MPI_Allreduce(m_hostTotal, &global, 1, datatype, MPI_MIN, MPI_COMM_WORLD);
            break;
        case Sum:
        case SumAbs:
            MPI_Allreduce(m_hostTotal, &global, 1, datatype, MPI_SUM, MPI_COMM_WORLD);
            break;
        default:
            MayDay<void>::Abort("Reduction::fetch | Error: Unknown reduction operator.");
            break;
    }
    barrier();
    return global;
#else
    return *m_hostTotal;
#endif
}

/// REDUCE
template<typename T, Operation OP, MemType MEM>
void Reduction<T,OP,MEM>::reduce(const T *a_data, const size_t a_size) {
    if (MEM == DEVICE)
    {
#ifdef PROTO_CUDA
        if (m_dynamic)
        {
            m_numThreads = min(a_size,(size_t)m_numThreads); // threads in a block
            m_numBlocks = (a_size+m_numThreads-1)/m_numThreads; // blocks in a kernel
            m_deviTemp = (T*)proto_malloc<DEVICE>(m_numBlocks*sizeof(T));
        }
        T* temp_h = (T*)proto_malloc<HOST>(m_numBlocks*sizeof(T));
        int shmem = (m_numThreads+m_warpSize-1)/m_warpSize; // warps in a block
        protoLaunchKernelMemAsyncGPU((kernel<T,OP>),m_numBlocks,shmem*m_warpSize,shmem*sizeof(T),
                (protoStream_t) 0, a_size, a_data, m_deviTemp, m_deviTotal);
        if (m_numBlocks > 1)
        {
            // each block in first kernel left partial reduction
            shmem = (m_numBlocks+m_warpSize-1)/m_warpSize;
            // update m_deviTotal
            protoLaunchKernelMemAsyncGPU((kernel<T,OP>),1,shmem*m_warpSize,shmem*sizeof(T),
                    (protoStream_t) 0, m_numBlocks, m_deviTemp, nullptr, m_deviTotal);
        }
#else
        MayDay<void>::Abort("Reduction::reduce | Error: Unrecognized acceleration flag.");
#endif
    } else if (MEM == HOST)
    {
        T val = *m_hostTotal;
        for(int it=0 ; it < a_size ; it++)
        {
            Reduction<T,OP>::update(val,a_data[it]);
        }
        *m_hostTotal = val;
    }
}
