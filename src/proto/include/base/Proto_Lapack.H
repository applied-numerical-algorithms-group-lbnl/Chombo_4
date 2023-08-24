#pragma once
#ifndef _PROTO_LAPACK_
#define _PROTO_LAPACK_

extern "C" {
    void dgemv_(char*, int*, int*, double*, double*,
            int*, double*, int*, double*, double*, int*);
    void dgemm_(char*, char*, int*, int*, int*, double*,
            double*, int*, double*, int*, double*, double*, int*);
    void dgetri_(int*, double*, int*, int*, double*, int*, int*);
    void dgetrf_(int*, int*, double*, int*, int*, int*);
}

#endif // end include guard
