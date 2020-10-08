#include "Proto.H"
#include "EBProto.H"
#include "Proto_GeometryService.H"

using namespace Proto;



template<typename T>
inline T&
expGetIrregData(T& a_s)
{
  return a_s;
}

template <CENTERING cent, typename  data_t, unsigned int ncomp>
inline IrregData<cent, data_t, ncomp>&
expGetIrregData(EBBoxData<cent, data_t, ncomp>& a_s)
{
  return a_s.getIrregData();
}


template<CENTERING cent,  typename data_t, unsigned int ncomp>
struct
uglyStruct
{
  data_t*        m_startPtr;
  unsigned int   m_varsize;
  unsigned int   m_offset;
  Point          m_index;
};

template <CENTERING cent, typename T>
inline T
getUglyStruct(const vector<EBIndex<cent> >& a_indices,
              T& a_T)
{
  return a_T;
}

template <CENTERING cent, typename  data_t, unsigned int ncomp>
inline vector< uglyStruct<cent, data_t, ncomp> >
getUglyStruct(const vector<EBIndex<cent> >& a_indices,
              IrregData<cent, data_t, ncomp>& a_s )
{
  data_t*      debPtr  = a_s.data();
//  printf("irreg data ptr = %p\n", debPtr);
  vector< uglyStruct<cent, data_t, ncomp> > retval;
  for(int ivec = 0; ivec < a_indices.size(); ivec++)
  {
    uglyStruct<cent, data_t, ncomp>  vecval;
    vecval.m_startPtr = a_s.data();
    vecval.m_varsize  = a_s.vecsize();
    vecval.m_offset   = a_s.index(a_indices[ivec], 0);
    vecval.m_index    = a_indices[ivec].m_pt;
    retval.push_back(vecval);
  }

//  printf("host return data ptr = %p\n", retval.data());
  return retval;
}

template <CENTERING cent, typename T>
inline T
cleanUpPtrr(T& a_T)
{
  return a_T;
}



template<typename T>
inline T&
expGetBoxData(T& a_s)
{
  return a_s;
}

template <CENTERING cent, typename  data_t, unsigned int ncomp>
inline BoxData<data_t, ncomp>&
expGetBoxData(EBBoxData<cent, data_t, ncomp>& a_s)
{
  return a_s.getRegData();
}

#ifdef PROTO_CUDA
///cuda specific functions
///
template<typename T>
__device__ __host__
inline T
expCudaGetVar(unsigned int ivec,  T a_s)
{
  return a_s;
}

///
template<CENTERING cent, typename data_t, unsigned int ncomp>
__device__ 
inline Var<data_t, ncomp>
expCudaGetVar(unsigned int a_ivec,
           uglyStruct<cent, data_t, ncomp>* a_dst)
{
  Var<data_t, ncomp> retval;

  const uglyStruct<cent, data_t, ncomp>*  rawptr = a_dst;
  const uglyStruct<cent, data_t, ncomp>&  ugly   = rawptr[a_ivec];
//  printf("expCudaGetVar rawptr   = %p \n", rawptr);
//  printf("expCudaGetVar startptr = %p \n", ugly.m_startPtr);
  for(int icomp = 0; icomp < ncomp; icomp++)
  {
    retval.m_ptrs[icomp] = ugly.m_startPtr + ugly.m_offset + (ugly.m_varsize*icomp);
  }
  return retval;
}




///going into this srcs are thrust_device_pointer<uglystruct> and other stuff
template<CENTERING cent, typename data_t,unsigned int ncomp,  typename Func, typename... Srcs>
__global__
void
vec_indexer(unsigned int a_begin, unsigned int a_end,Func a_body, 
            uglyStruct<cent, data_t, ncomp>*  a_dst, Srcs... a_srcs)
{
  int idx = threadIdx.x + blockIdx.x*blockDim.x;
  if (idx >= a_begin && idx < a_end)
  {
    a_body(expCudaGetVar(idx, a_dst), expCudaGetVar(idx, a_srcs)...);
  }
}


template<CENTERING cent, typename data_t,unsigned int ncomp,  typename Func, typename... Srcs>
__global__
void
vec_indexer_i(unsigned int a_begin, unsigned int a_end,Func a_body, 
              uglyStruct<cent, data_t, ncomp>*  a_dst, Srcs... a_srcs)
{
  int idx = threadIdx.x + blockIdx.x*blockDim.x;
  if (idx >= a_begin && idx < a_end)
  {
    a_body(a_dst[idx].m_index.m_tuple, expCudaGetVar(idx, a_dst), expCudaGetVar(idx, a_srcs)...);
  }
}

///going into this srcs are uglystruct* and other stuff
template <typename T>
inline int
cleanUpPtrs(T a_T)
{
  return 0;
}
///going into this srcs are uglystruct* and other stuff
template <CENTERING cent, typename data_t,unsigned int ncomp>
inline int
cleanUpPtrs(uglyStruct<cent, data_t, ncomp> * a_ptr)
{
  protoFree(a_ptr);
  return 0;
}

///
template<typename... Srcs>
inline void
expEmptyFunc(Srcs... a_srcs)
{
}

template<typename... Srcs>
__global__
inline void
expEmptyFuncGPU(Srcs... a_srcs)
{
}

///going into this srcs are uglystruct* and other stuff
template<CENTERING cent, typename data_t,unsigned int ncomp,  typename Func, typename... Srcs>
void
cudaExpVectorFunc(const Func& a_F, unsigned int a_Nvec, 
               uglyStruct<cent, data_t, ncomp> * a_dst,Srcs... a_srcs)
{
//  printf("cudavecf: dst  = %p\n", a_dst);
  //printf("cudavecf: src  = %p\n", a_firstsrc);
  protoStream_t curstream = DisjointBoxLayout::getCurrentStream();
  unsigned int N = a_Nvec;
  unsigned int stride = a_Nvec;
  unsigned int blocks = 1;
  size_t smem = 0;
  protoLaunchKernelMemAsync( vec_indexer, blocks, stride, smem, curstream,
    0, N, a_F, a_dst, a_srcs...);

  //there is a cudaMalloc that happens above so we have to delete
  expEmptyFunc(cleanUpPtrs(a_dst ), (cleanUpPtrs(a_srcs))...); 
}


///going into this srcs are uglystruct* and other stuff
template<CENTERING cent, typename data_t,unsigned int ncomp,  typename Func, typename... Srcs>
void
cudaExpVectorFunc_i(const Func& a_F, unsigned int a_Nvec, 
                 uglyStruct<cent, data_t, ncomp> * a_dst,Srcs... a_srcs)
{
//  printf("cudavecf: dst  = %p\n", a_dst);
  //printf("cudavecf: src  = %p\n", a_firstsrc);
  protoStream_t curstream = DisjointBoxLayout::getCurrentStream();
  const int N = a_Nvec;
  unsigned int stride = a_Nvec;
  unsigned int blocks = 1;
  size_t smem = 0;
  protoLaunchKernelMemAsync(vec_indexer_i, blocks, stride, smem, curstream,
    0, N, a_F, a_dst, a_srcs...);

  //there is a cudaMalloc that happens above so we have to delete
  expEmptyFunc(cleanUpPtrs(a_dst ), (cleanUpPtrs(a_srcs))...); 
}


template <CENTERING cent, typename T>
inline T
cudaGetUglyStruct(const vector<EBIndex<cent> >& a_indices,
                  T& a_T)
{
  return a_T;
}
///
template <CENTERING cent, typename  data_t, unsigned int ncomp>
inline  uglyStruct<cent, data_t, ncomp>*
cudaGetUglyStruct(const vector<EBIndex<cent> >& a_indices,
                  IrregData<cent, data_t, ncomp>& a_s )
{
  vector< uglyStruct<cent, data_t, ncomp> > hostvec = getUglyStruct(a_indices, a_s);

  size_t memsize = hostvec.size()*sizeof(uglyStruct<cent, data_t, ncomp>);
  uglyStruct<cent, data_t, ncomp>* retval;
  protoMalloc(&retval, memsize);
  protoMemcpy(retval, hostvec.data(), memsize, protoMemcpyHostToDevice);

  //this copies from the host to the device
//  printf("cgus: device host vector ptr = %p\n", hostvec.data());
//  printf("cgus: device return data ptr = %p\n", retval);
  return retval;
}
///going into this srcs are IrregDatas and other stuff
template<CENTERING cent, typename  data_t, unsigned int ncomp, typename Func, typename... Srcs>
inline void
cudaExpebforAllIrreg(const Func& a_F, const Box& a_box,
                  IrregData<cent, data_t, ncomp>& a_dst,
                  Srcs&...  a_srcs)
{
  //indicies into irreg vector that correspond to input box
  const vector<EBIndex<cent> >& dstvofs = a_dst.getIndices(a_box);
  unsigned int vecsize = a_dst.vecsize();
  if(vecsize > 0)
  {
//    printf("cudaexpebforall: dst  = %p\n", a_dst.data());
    cudaExpVectorFunc(a_F, vecsize, cudaGetUglyStruct(dstvofs, a_dst), 
                   (cudaGetUglyStruct(dstvofs, a_srcs))...);

   }
}


///going into this srcs are IrregDatas and other stuff
template<CENTERING cent, typename  data_t, unsigned int ncomp, typename Func, typename... Srcs>
inline void
cudaExpebforAllIrreg_i(const Func& a_F, const Box& a_box,
                    IrregData<cent, data_t, ncomp>& a_dst,
                    Srcs&...  a_srcs)
{
  //indicies into irreg vector that correspond to input box
  const vector<EBIndex<cent> >& dstvofs = a_dst.getIndices(a_box);
  unsigned int vecsize = a_dst.vecsize();
  if(vecsize > 0)
  {
//    printf("cudaexpebforall: dst  = %p\n", a_dst.data());
    cudaExpVectorFunc_i(a_F, vecsize, cudaGetUglyStruct(dstvofs, a_dst), 
                     (cudaGetUglyStruct(dstvofs, a_srcs))...);

   }
}

template<typename Func, typename... Srcs>
inline void
cudaExpebforall(const Func & a_F,  Box a_box, Srcs&... a_srcs)
{
//call regular forall
  forallInPlaceBase(a_F, a_box, (expGetBoxData(a_srcs))...);
  
//do the same thing for the irregular data
  cudaExpebforAllIrreg(a_F, a_box, expGetIrregData(a_srcs)...);
}


template<typename Func, typename... Srcs>
inline void
cudaExpebforall_i(const Func & a_F,  Box a_box, Srcs&... a_srcs)
{
//call regular forall
  forallInPlace_i(a_F, a_box, (expGetBoxData(a_srcs))...);
  
//do the same thing for the irregular data
  cudaExpebforAllIrreg_i(a_F, a_box, expGetIrregData(a_srcs)...);
}

#else
///cpu-only specific functions

///
template<typename T>
inline T
expGetVar(unsigned int ivec,  T a_s)
{
  return a_s;
}

///
template<CENTERING cent, typename data_t, unsigned int ncomp>
inline Var<data_t, ncomp>
expGetVar(unsigned int a_ivec,
       vector< uglyStruct<cent, data_t, ncomp> > a_dst)
{
  Var<data_t, ncomp> retval;
  const uglyStruct<cent, data_t, ncomp> ugly = a_dst[a_ivec];
  for(int icomp = 0; icomp < ncomp; icomp++)
  {
    retval.m_ptrs[icomp] = ugly.m_startPtr + ugly.m_offset + (ugly.m_varsize*icomp);
  }
  return retval;
}

///going into this srcs are vector<uglystruct> and other stuff
template<CENTERING cent, typename data_t,unsigned int ncomp,  typename Func, typename... Srcs>
void
hostExpVectorFunc(const Func& a_F, vector< uglyStruct<cent, data_t, ncomp> > a_dst, Srcs... a_srcs)
{
  for(unsigned int ivec = 0; ivec < a_dst.size(); ivec++)
  {
    a_F(expGetVar(ivec, a_dst), (expGetVar(ivec, a_srcs))...);
  }
       
}


///going into this srcs are vector<uglystruct> and other stuff
template<CENTERING cent, typename data_t,unsigned int ncomp,  typename Func, typename... Srcs>
void
hostExpVectorFunc_i(const Func& a_F, vector< uglyStruct<cent, data_t, ncomp> > a_dst, Srcs... a_srcs)
{
  for(unsigned int ivec = 0; ivec < a_dst.size(); ivec++)
  {
    Point pt = a_dst[ivec].m_index;
    a_F(pt.m_tuple, expGetVar(ivec, a_dst), (expGetVar(ivec, a_srcs))...);
  }
       
}

///going into this srcs are IrregDatas and other stuff
template<CENTERING cent, typename  data_t, unsigned int ncomp, typename Func, typename... Srcs>
inline void
hostExpebforAllIrreg(const Func& a_F, const Box& a_box,
                  IrregData<cent, data_t, ncomp>& a_dst,
                  Srcs&...  a_srcs)
{
//indicies into irreg vector that correspond to input box
  vector<EBIndex<cent> > dstvofs = a_dst.getIndices(a_box);
  hostExpVectorFunc(a_F, getUglyStruct(dstvofs, a_dst), (getUglyStruct(dstvofs, a_srcs))...);
}


///going into this srcs are IrregDatas and other stuff
template<CENTERING cent, typename  data_t, unsigned int ncomp, typename Func, typename... Srcs>
inline void
hostExpebforAllIrreg_i(const Func& a_F, const Box& a_box,
                    IrregData<cent, data_t, ncomp>& a_dst,
                    Srcs&...  a_srcs)
{
//indicies into irreg vector that correspond to input box
  vector<EBIndex<cent> > dstvofs = a_dst.getIndices(a_box);
  hostExpVectorFunc_i(a_F, getUglyStruct(dstvofs, a_dst), (getUglyStruct(dstvofs, a_srcs))...);
}


///going into this srcs are EBBoxDatas and other stuff
template<typename Func, typename... Srcs>
inline void
hostExpebforall(const Func & a_F,  Box a_box, Srcs&... a_srcs)
{
//call regular forall
  forallInPlaceBase(a_F, a_box, (expGetBoxData(a_srcs))...);
  
//do the same thing for the irregular data
  hostExpebforAllIrreg(a_F, a_box, expGetIrregData(a_srcs)...);
}

///going into this srcs are EBBoxDatas and other stuff
template<typename Func, typename... Srcs>
inline void
hostExpebforall_i(const Func & a_F,  Box a_box, Srcs&... a_srcs)
{
//call regular forall
  forallInPlace_i(a_F, a_box, (expGetBoxData(a_srcs))...);
  
//do the same thing for the irregular data
  hostExpebforAllIrreg_i(a_F, a_box, expGetIrregData(a_srcs)...);
}

#endif

///version that does not send the point to the function
template<typename Func, typename... Srcs>
inline void ExpebforallInPlace(unsigned long long int a_num_flops_point,
                            const char*            a_timername,
                            const Func & a_F,  Box a_box, Srcs&... a_srcs)
{
  PR_TIME(a_timername);

//  printf("in expebforall function pointer = %p\n", &a_F);
  unsigned long long int boxfloops = a_num_flops_point*a_box.size();

#ifdef PROTO_CUDA
  cudaExpebforall(a_F, a_box, a_srcs...);
  protoDeviceSynchronize();
#else
  hostExpebforall(a_F, a_box, a_srcs...);
#endif
  PR_FLOPS(boxfloops);
}


/////version that sends the point to the function
template<typename Func, typename... Srcs>
inline void ExpebforallInPlace_i(unsigned long long int a_num_flops_point,
                              const char*            a_timername,
                              const Func & a_F,  Box a_box, Srcs&... a_srcs)
{
  PR_TIME(a_timername);

  unsigned long long int boxfloops = a_num_flops_point*a_box.size();

#ifdef PROTO_CUDA
  cudaExpebforall_i(a_F, a_box, a_srcs...);
#else
  hostExpebforall_i(a_F, a_box, a_srcs...);
#endif


  PR_FLOPS(boxfloops);
}


///  after this are specific to the test
typedef Var<double,DIM> V;


PROTO_KERNEL_START 
void UsetUF(V a_U, double  a_val)
{
//  printf("in set U\n");
//  printf("setu: uptr[0] = %p, uptr[1] = %p\n",a_U.m_ptrs[0],a_U.m_ptrs[1]);
  for(int idir = 0; idir < DIM; idir++)
  {
    a_U(idir) = a_val;
    if(a_U(idir) != a_val)
    {
      printf("p1: values do not match \n");
      printf("setu: val = %f, uval = %f\n",a_val, a_U(idir));
    }
  }
}
PROTO_KERNEL_END(UsetUF, UsetU)


PROTO_KERNEL_START 
void setUptF(int  a_p[DIM], V a_U, double  a_val)
{
  for(int idir = 0; idir < DIM; idir++)
  {
    a_U(idir) = a_val;
    if(a_U(idir) != a_val)
    {
      printf("upt: values do not match \n");
      printf("upt: val = %f, uval = %f\n",a_val, a_U(idir));
    }
  }
}
PROTO_KERNEL_END(setUptF, setUpt)


PROTO_KERNEL_START 
void VsetVF(V a_V, double  a_val, int a_intvar)
{
//  printf("setv: vptr[0] = %p, vptr[1] = %p\n",a_V.m_ptrs[0],a_V.m_ptrs[1]);
//  printf("in set V\n");
 for(int idir = 0; idir < DIM; idir++)
 {
   a_V(idir) = a_val;
   if(a_V(idir) != a_val)
   {
     printf("setv: values do not match \n");
     //     printf("setv: val = %f, vval = %f\n",a_val, a_V(idir));
   }
 }
}
PROTO_KERNEL_END(VsetVF, VsetV)

PROTO_KERNEL_START 
void setVptF(int  a_p[DIM], V a_V, double  a_val, int a_vvar)
{
  for(int idir = 0; idir < DIM; idir++)
  {
   a_V(idir) = a_val;
   if(a_V(idir) != a_val)
   {
     printf("vpt: values do not match \n");
//     printf("setv: val = %f, vval = %f\n",a_val, a_V(idir));
   }
  }
}
PROTO_KERNEL_END(setVptF, setVpt)


PROTO_KERNEL_START 
void WsetWtoUplusVF(V a_W,
                    V a_U,
                    V a_V,
                    double  a_val)
{
  for(int idir = 0; idir < DIM; idir++)
  {
    a_W(idir) = a_U(idir) + a_V(idir);
    if(a_W(idir) != a_val)
    {
      printf("w: values do not match\n");
    }
  }

}
PROTO_KERNEL_END(WsetWtoUplusVF, WsetWtoUplusV)

PROTO_KERNEL_START 
void setWtoUplusVptF(int a_p[DIM],
                     V a_W,
                     V a_U,
                     V a_V,
                     double  a_val)
{
//  printf("setw: uptr[0] = %p, uptr[1] = %p\n" ,a_U.m_ptrs[0],a_U.m_ptrs[1]);
//  printf("setw: vptr[0] = %p, vptr[1] = %p\n" ,a_V.m_ptrs[0],a_V.m_ptrs[1]);
//  printf("setw: wptr[0] = %p, wptr[1] = %p\n" ,a_W.m_ptrs[0],a_W.m_ptrs[1]);
//  printf("in set W \n");
  for(int idir = 0; idir < DIM; idir++)
  {
    a_W(idir) = a_U(idir) + a_V(idir);
    if(a_W(idir) != a_val)
    {
      printf("wpt: values do not match\n");
    }
  }
}
PROTO_KERNEL_END(setWtoUplusVptF, setWtoUplusVpt)

int main(int argc, char* argv[])
{

  int nx = 16;
  int maxGrid = nx; //this test is for a single grid
  Box domain(Point::Zeros(), Point::Ones(nx-1));
  double dx = 1.0/domain.size(0);
  std::array<bool, DIM> periodic;
  for(int idir = 0; idir < DIM; idir++) periodic[idir]=true;
  DisjointBoxLayout grids(domain, maxGrid, periodic);

  Point dataGhost = Point::Zeros();
  Point geomGhost = Point::Zeros();

  RealVect ABC = RealVect::Unit();
  RealVect X0 = RealVect::Unit();
//  X0 *= 100.;
  X0 *= 0.5;
  RealVect origin= RealVect::Zero();
  double R = 0.25;
  shared_ptr<BaseIF>              impfunc(new SimpleEllipsoidIF(ABC, X0, R, false));
  shared_ptr<GeometryService<2> > geoserv(new GeometryService<2>(impfunc, origin, dx, domain, grids, geomGhost, 0));
  shared_ptr<LevelData<EBGraph> > graphs = geoserv->getGraphs(domain);

  for(int ibox = 0; ibox < grids.size(); ibox++)
  {
    Box grid = grids[ibox];
    EBBoxData<CELL, double, DIM> U(grid,(*graphs)[ibox]);
    EBBoxData<CELL, double, DIM> V(grid,(*graphs)[ibox]);
    EBBoxData<CELL, double, DIM> W(grid,(*graphs)[ibox]);
    unsigned long long int numFlopsPt = 0;
    double uval = 1;
    double vval = 2;
    //printf("going into setU\n");
    ExpebforallInPlace(numFlopsPt, "setU", UsetU, grid, U, uval);
    //printf("going into setV\n");
    int vvar = -1;
    ExpebforallInPlace(numFlopsPt, "setV", VsetV, grid, V, vval, vvar);  //tweaking signature to clarify compilers job
    double wval = 3;
    //printf("going into setWtoUPlusV\n");
    ExpebforallInPlace(numFlopsPt, "setWtoUPlusV", WsetWtoUplusV, grid, W, U, V, wval);

    uval = 2;
    vval = 5;
    wval = 7;
    printf("going into setUpt\n");
    ExpebforallInPlace_i(numFlopsPt, "setU", setUpt, grid, U, uval);
    printf("going into setVpt\n");
    ExpebforallInPlace_i(numFlopsPt, "setV", setVpt, grid, V, vval, vvar);
    printf("going into setWpt\n");
    ExpebforallInPlace_i(numFlopsPt, "setWtoUPlusV", setWtoUplusVpt, grid, W, U, V, wval);

  }

}








