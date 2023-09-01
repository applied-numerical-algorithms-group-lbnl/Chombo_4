
#ifdef PROTO_CUDA
#ifdef PROTO_HIP
inline void hipErrorCheck(hipError_t a_err, std::string a_message)
{
    if (a_err != hipSuccess)
    {
        pout() << a_message;
        pout() << hipGetErrorString(a_err) << std::endl;
    }
}
#else
inline void cudaErrorCheck(cudaError_t a_err, std::string a_message)
{
    if (a_err != cudaSuccess)
    {
        pout() << a_message;
        pout() << cudaGetErrorString(a_err) << std::endl;
    }
}
#endif
#endif

MemType pointerMemType(const void* a_ptr)
{
#ifdef PROTO_CUDA
#ifdef PROTO_HIP
    hipPointerAttribute_t attributes;
    hipPointerGetAttributes(&attributes, a_ptr);
    if (attributes.devicePointer != NULL) { return DEVICE; }
    else { return HOST; }
#else
    cudaPointerAttributes attributes;
    cudaPointerGetAttributes(&attributes, a_ptr);
    if (attributes.devicePointer != NULL) { return DEVICE; }
    else { return HOST; }
#endif
#endif
    return HOST;
}

//TODO: Add an asynchronous option
template<MemType SRC_MEM, MemType DST_MEM>
void proto_memcpy(const void* a_src, void* a_dst, size_t a_size)
{
    auto srcMem = pointerMemType(a_src);
    auto dstMem = pointerMemType(a_dst);
    PROTO_ASSERT(srcMem == SRC_MEM,
        std::string("protoMemcpy | Error: Source pointer is of type " + parseMemType(srcMem) 
        + "but should be " + parseMemType(SRC_MEM) + "\n").c_str());
    PROTO_ASSERT(dstMem == DST_MEM,
        std::string("protoMemcpy | Error: Destination pointer is of type " + parseMemType(dstMem) 
        + "but should be " + parseMemType(DST_MEM) + "\n").c_str());
    if (SRC_MEM == HOST && DST_MEM == HOST)
    {
        std::memcpy(a_dst, a_src, a_size);
    } else {
#ifdef PROTO_CUDA
#ifdef PROTO_HIP
        hipMemcpyKind copyDir;
        if (     SRC_MEM == HOST   && DST_MEM == DEVICE) { copyDir = hipMemcpyHostToDevice; }
        else if (SRC_MEM == DEVICE && DST_MEM == HOST)   { copyDir = hipMemcpyDeviceToHost; }
        else if (SRC_MEM == DEVICE && DST_MEM == DEVICE) { copyDir = hipMemcpyDeviceToDevice; }
        hipErrorCheck(hipMemcpy(a_dst, a_src, a_size, copyDir), "proto_memcpy | Error: ");
#else
        cudaMemcpyKind copyDir;
        if (     SRC_MEM == HOST   && DST_MEM == DEVICE) { copyDir = cudaMemcpyHostToDevice; }
        else if (SRC_MEM == DEVICE && DST_MEM == HOST)   { copyDir = cudaMemcpyDeviceToHost; }
        else if (SRC_MEM == DEVICE && DST_MEM == DEVICE) { copyDir = cudaMemcpyDeviceToDevice; }
        cudaErrorCheck(cudaMemcpy(a_dst, a_src, a_size, copyDir), "proto_memcpy | Error: ");
#endif    
#endif    
    }
}

template<MemType MEM>
void* proto_malloc(size_t a_size)
{
    void* buffer;
    if (MEM == HOST)
    {
        buffer = std::malloc(a_size);
    } else if (MEM == DEVICE)
    {
#ifdef PROTO_CUDA
#ifdef PROTO_HIP
        hipErrorCheck(hipMalloc(&buffer, a_size), "proto_malloc | Error: ");
#else
        cudaErrorCheck(cudaMalloc(&buffer, a_size), "proto_malloc | Error: ");
#endif
#endif
    }
    return buffer;
}

template<MemType MEM>
void proto_free(void* a_buffer)
{
    if (a_buffer == nullptr) 
    {
        return;
    }
    auto mem = pointerMemType(a_buffer);
    PROTO_ASSERT(mem == MEM,
        "proto_free | Error: Pointer is of type %s but should be %s.",
        parseMemType(mem).c_str(), parseMemType(MEM).c_str());
    if (MEM == HOST)
    {
        std::free(a_buffer);
    } else if (MEM == DEVICE)
    {
#ifdef PROTO_CUDA
#ifdef PROTO_HIP
        hipErrorCheck(hipFree(a_buffer), "proto_free | Error: ");
#else
        cudaErrorCheck(cudaFree(a_buffer), "proto_free | Error: ");
#endif
#endif
    }
}
