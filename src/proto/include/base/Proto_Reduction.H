#ifndef __PROTO_REDUCTION_H__
#define __PROTO_REDUCTION_H__

//#include "Proto_gpu.H"
#include "Proto_MemType.H"
#include "Proto_cuda.H"
#include "Proto_macros.H"
#include "Proto_SPMD.H"
#include <limits>
#include <typeinfo>

namespace Proto {

enum Operation { Max, Min, Abs, Sum, SumAbs };
enum Atomic { Warp, Block, None };

constexpr int line = 128; // bytes in a cache line

#ifdef PROTO_CUDA // namespace collision in host builds
template<typename T>
CUDA_DECORATION
static constexpr T max(T a, T b) { return (a > b ? a : b); }
template<typename T>
CUDA_DECORATION
static constexpr T min(T a, T b) { return (a < b ? a : b); }
#endif

template<typename T, Operation OP, MemType MEM=MEMTYPE_DEFAULT>
class Reduction 
{
public:
    Reduction<T,OP,MEM>() : Reduction<T,OP,MEM>(false) {}

    /// Constructor
    /**
        Creates a Reduction operator with or without dynamic block allocation.
        If <code>a_dynamic==true</code>, the minimum number of GPU blocks
        needed to compute a reduction will be determined on each call to
        <code>Reduction::reduce</code>. Otherwise, the maximum bumber of blocks
        available on the device will be used instead.

        \param a_dynamic    Determine number of GPU blocks dynamically?
    */
    Reduction<T,OP,MEM>(bool a_dynamic);

    /// Destructor
    ~Reduction<T,OP,MEM>(); 
    
    /// Initialize Value
    /**
        Initialize a value based on OP.
        Minima Ops: <code>a_value</code> -> max value of T
        Maxima Ops: <code>a_value</code> -> min value of T
        Sum Ops:    <code>a_value</code> -> 0
    */
    CUDA_DECORATION
    static T init();

    /// Update Value
    /**
        Computed an updated value by executing OP on two values.

        \param a_v1 A value compared with OP
        \param a_v2 A value compared with OP
    */
    CUDA_DECORATION
    static void update(T& a_v1, const T a_v2);
    
    /// Get Reduction
    /**
        Returns the computed reduction, waiting for kernel completion on GPUs.
        If MPI is enabled, the result is communicated to all processes.
    */
    T fetch();

    /// Compute Reduction
    /**
        Calculates a reduction on the buffer <code>a_data</code> of size <code>a_size</code>.
        Here <code>a_size</code> is the number of elements of <code>a_data</code>.
        Subsequent calls to <code>reduce</code> will update the reduction until
        <code>Reduction::reset()</code> is called. 

        \param a_data   A buffer of data.
        \param a_size   The number of elements in <code>a_data</code>.
    */
    void reduce(const T *a_data, const size_t a_size); // configures and calls the kernel

    /// Reset Reduction
    /**
        Reinitializes the reduction.
    */
    void reset(); // initializes pointer based on operation

private:
    T *m_hostTotal;
    T* m_hostTemp;
#ifdef PROTO_CUDA
    bool m_dynamic;
    T *m_deviTemp, *m_deviTotal;
    int m_numThreads, m_numBlocks, m_warpSize;
#endif
};

template<typename T, Operation OP>
CUDA_KERNEL
void initKernel(T* ptr) 
{
    *ptr = Reduction<T,OP>::init();
}

#ifdef PROTO_CUDA

//=================================================================================================
// KERNEL AND DEVICE CODE

template<typename T, Operation OP>
__device__ // reduction within a warp, only lane 0 gets right result
void warpOp(T& val, size_t idx, size_t size) {
    unsigned mask = 0xffffffff; // FULL_MASK
    for (unsigned int delta = warpSize/2; delta > 0; delta /= 2)
    {
#ifdef PROTO_HIP
        // __shfl_down depreciated since CUDA 9.0 but __shfl_down_sync isn't available with hipcc (17/07/2020)
        Reduction<T,OP>::update(val, (T)__shfl_down(val,delta));
#else
        Reduction<T,OP>::update(val, (T)__shfl_down_sync(mask,val,delta));
#endif
    }
}

template<typename T, Operation OP>
__device__ // reduce within a block by storing partial warp reductions in shmem
void blockOp(T& val, int idx, int size) {
    PR_assert(blockDim.x <= warpSize*warpSize); // otherwise block won't be reduced by end of function
    extern __shared__ __align__(sizeof(T)) unsigned char shdata[];
    T *shmem = reinterpret_cast<T*>(shdata);
    int lane = threadIdx.x % warpSize;
    int wid = threadIdx.x / warpSize;
    int warps = (blockDim.x+warpSize-1)/warpSize;

    warpOp<T,OP>(val, idx, size);
    if (warps == 1) return;

    if (!lane) shmem[wid] = val; // first lane of each warp fills memory
    __syncthreads();
    if (!wid) {
        val = (threadIdx.x < warps ? shmem[lane] : Reduction<T,OP>::init());
    // only first lane of first warp ends up with real value
        warpOp<T,OP>(val, threadIdx.x, warps);
    }
}

template<typename T, Operation OP>
__global__ // every block's sum is written to out[]
void kernel(size_t size, const T* in, T* out, T* val)
{ // if called twice by reduce, out=nullptr
    PR_assert(gridDim.x <= line/sizeof(T)); // each block writes to unique out[] index
    int idx = blockIdx.x*blockDim.x + threadIdx.x;
    T ret = Reduction<T,OP>::init();
    for (size_t i = idx; i < size; i += blockDim.x*gridDim.x)
    {
        Reduction<T,OP>::update(ret,in[i]);
    }
    blockOp<T,OP>(ret, idx, size);
    if (!threadIdx.x)
    {
        if (gridDim.x > 1)
        {
            out[blockIdx.x] = ret; // block result is reused in 2nd call
        }
        else // only one block in this kernel -- always true for 2nd call
        {
            Reduction<T,OP>::update(*val, ret);
        }
    }
}
#endif

#include "implem/Proto_ReductionImplem.H"
} // end namespace Proto

#endif  // __PROTO_REDUCTION_H__
