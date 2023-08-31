#pragma once


template<typename T, BoxDataOp op>
struct opKernel{

#ifdef PROTO_CUDA
  CUDA_DECORATION
  void launchKernel(T* x, T* y, const unsigned int size) 
  {
    int id = blockIdx.x*blockDim.x+threadIdx.x;
    if (id < size)
    switch(op)
    {
      case BoxDataOp::Add: 
        x[id] += y[id];
        break;
      case BoxDataOp::Subtract: 
        x[id] -= y[id];
        break;
      case BoxDataOp::Divide: 
        x[id] /= y[id];
        break;
      case BoxDataOp::Multiply: 
        x[id] *= y[id];
        break;
    }
  }
#else
  void launchKernel(T* x, T* y, const unsigned int size) 
  {
    switch(op)
    {
      case BoxDataOp::Add: 
      for(int id = 0 ; id < size ; id++)
        x[id] += y[id];
        break;
      case BoxDataOp::Subtract: 
      for(int id = 0 ; id < size ; id++)
        x[id] -= y[id];
        break;
      case BoxDataOp::Divide: 
      for(int id = 0 ; id < size ; id++)
        x[id] /= y[id];
        break;
      case BoxDataOp::Multiply: 
      for(int id = 0 ; id < size ; id++)
        x[id] *= y[id];
        break;
    }
  }
#endif
};

