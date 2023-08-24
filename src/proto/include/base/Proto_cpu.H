#pragma once
#include <type_traits>

#define cpuMalloc(ptr,nbBytes)     ptr = (decltype(ptr)) malloc(nbBytes)

#define protoMallocCPU(PTR,NBYTES)    storeMemInfo(Proto::HOST,NBYTES); countMallocHost(cpuMalloc(PTR, NBYTES))
#define protoFreeCPU(PTR)             free(PTR)

#ifndef PROTO_CUDA
#define protoHostAlloc(a,b)        storeMemInfo(Proto::HOST,b); countMallocHost(cpuMalloc(a, b))
#define protoFreeHost(PTR)         free(PTR)
#else
#ifdef PROTO_HIP
#define protoHostAlloc(a,b)        storeMemInfo(Proto::HOST,b); countMallocHost(hipHostMalloc(&a, b))
#define protoFreeHost(PTR)         hipHostFree(PTR)
#else
#define protoHostAlloc(a,b)        storeMemInfo(Proto::HOST,b); countMallocHost(cudaMallocHost(&a, b))
#define protoFreeHost(PTR)         cudaFreeHost(PTR)
#endif
#endif
// trick

#ifndef PROTO_CUDA
//static bool isGpu=false;
#define protoMemcpyDeviceToDevice false
#define protoMemcpyHostToDevice   false
#define protoMemcpyDeviceToHost   false
#endif
// end trick


// COPY


#define protoMemcpyCPU(to,from,size,copyType)             std::memcpy(to,from,size)
#define protoMemcpyAsyncCPU(to,from,size,copyType,stream) std::memcpy(to,from,size)
#define protoMemcpyFromSymbolCPU(a,b,c,d,e)               std::memcpy(a,b,c)
#define protoMemcpyToSymbolCPU(a,b,c,d,e)                 std::memcpy(a,b,c)

#define protoDeviceSynchronizeCPU()
