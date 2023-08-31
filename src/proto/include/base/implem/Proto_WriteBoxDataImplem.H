#include <cstdio>
#include <iostream>


//using namespace std;

inline const char* WriteBoxData(BoxData<double>& a_array,
                    double dx)
{
    static int fileCount = 0;
    static char nameBuffer[20];

    if (!a_array.defined())
    {
        return nameBuffer;
    }

    sprintf(nameBuffer, "data%06d",fileCount);
    WriteBoxData(nameBuffer, a_array, dx);

    fileCount++;

    return nameBuffer;
}

inline void WriteBoxData(const char*           a_filename,
             BoxData<double>&  a_array,
             double dx)
{
    if (!a_array.defined()){return;}
    if(DIM==1)
    {
            std::string file_str(a_filename);
            file_str+=".curve";
            std::ofstream fs;
            fs.open(file_str);
            fs << "#var0" << std::endl;
            for(int i=0; i<a_array.box().size(); i++)
                    fs<< i*dx << " " << *a_array[i] << std::endl;
            fs.close();
    }
    else
    {
            const char* const vars[] = {a_filename};
            double origin[DIM];
            for (int ii = 0; ii < DIM; ii++)
            {
                    origin[ii] = 0.0;
            }

            WriteBoxData(a_filename, a_array, vars, origin, dx);
    }
}

template<class T, unsigned int C, MemType MEMTYPE, unsigned char D, unsigned char E>
inline void WriteBoxData(const char*       a_filename,
                         BoxData<T,C,MEMTYPE,D,E>&      a_array,
        const char* const      a_varnames[C*D*E],
        const double*          a_corner,
        const double&          a_h)
{
//  static StackAlloc allocator;
    int dim[3] = {1,1,1};

    int vardims[C*D*E];
    int centering[C*D*E];

    double* x[3] {nullptr,nullptr,nullptr};
    T* hostPtr[C*D*E];
    T* wholeHostPtr = nullptr;
#ifdef PROTO_CUDA
    int nsize = a_array.size();
    wholeHostPtr = new T[nsize];
    T* wholeDevPtr  = a_array.data();
    //get memory to host
    protoMemcpyGPU(wholeHostPtr, wholeDevPtr, nsize*sizeof(T), protoMemcpyDeviceToHost);
    protoError err = protoGetLastError();
    if (err != protoSuccess)
    {
        fprintf(stderr, "protoGetLastError() failed at %s:%i : %s\n",
                __FILE__, __LINE__, protoGetErrorString(err));
    }

#else
    wholeHostPtr = a_array.data();
#endif
    int ilin = 0;
    for(unsigned int ee = 0; ee < E; ++ee)
        for(unsigned int dd = 0; dd < D; ++dd)
            for(unsigned int cc = 0; cc < C; ++cc)
            {
                unsigned int iindex = a_array.index(a_array.box().low(), cc, dd, ee);
                hostPtr  [ilin] = &(wholeHostPtr[iindex]);
                vardims  [ilin] = 1;
                centering[ilin] = 0;
                ilin++;
            }


    Point lo, hi;
    const Box& box = a_array.box();
    lo = box.low();
    hi = box.high();
    for (int i = 0; i < DIM; i++)
    {
        dim[i] = hi[i] - lo[i] + 2;
        //x[i] = (double*)allocator.alloc((dim[i]+1)*sizeof(double));
        x[i] = new double[dim[i]+1];
        for (int k = 0; k < dim[i]+1; k++)
        {
            x[i][k] = (double)(a_corner[i] + a_h * k);
        }
    }

    if (DIM == 2)
    {
      //x[2] = (double*)allocator.alloc(2*sizeof(double));
      x[2] = new double[2];
        x[2][0] = 0.0;
        x[2][1] = 1.0;
    }


    write_rectilinear_mesh( a_filename,
            dim,
            x[0],x[1],x[2],
            C*D*E,
            vardims,
            centering,
            a_varnames,
            hostPtr);

    // allocator.free(x[2]); // it's stack, have to do it backwards :-)
    //allocator.free(x[1]);
    //allocator.free(x[0]);
    delete[] x[2];
    delete[] x[1];
    delete[] x[0];
#ifdef PROTO_CUDA
    delete[] wholeHostPtr;
#endif

}
