#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif
#ifndef _Chombo_BFM_H_
#define _Chombo_BFM_H_

#include "Chombo_Box.H"
#include "Chombo_NamespaceHeader.H"





template<unsigned char D, typename T1, typename T2, typename F>
struct BFMImpl
{
  static void f(T1* __restrict t1, const T2* __restrict t2, int* lo1, int* lo2, int* stride1, int* stride2, int* counting,F f)
  {
    t1+=lo1[D]*stride1[D-1];
    t2+=lo2[D]*stride2[D-1];
    for(int i=0; i<counting[D]; ++i, t1+=stride1[D-1], t2+=stride2[D-1]) 
      BFMImpl<D-1,T1,T2,F>::f(t1, t2, lo1, lo2, stride1, stride2, counting, f);
  }
};
 
template<typename T1, typename T2, typename F>
struct BFMImpl<0,T1,T2,F>
{
  static void f(T1* __restrict__ t1, const T2* __restrict__ t2, int* lo1, int* lo2, int* stride1, int* stride2, int* counting,F f)
  {
    int trips=counting[0];
    int l1=lo1[0];
    int l2=lo2[0];
    for(int i=0; i!=trips; ++i) f(t1+i+l1,t2[i+l2]);
  }
};

template<unsigned char D, typename T1, typename F>
struct BFMImpl1
{
  static void f(T1* __restrict t1, int* lo1, int* stride1, int* counting,F f)
  {
    t1+=lo1[D]*stride1[D-1];
    for(int i=0; i<counting[D]; ++i, t1+=stride1[D-1]) 
      BFMImpl1<D-1,T1,F>::f(t1, lo1, stride1, counting, f);
  }
};
 
template<typename T1, typename F>
struct BFMImpl1<0,T1,F>
{
  static void f(T1* __restrict__ t1, int* lo1, int* stride1, int* counting,F f)
  {
    int trips=counting[0];
    int l1=lo1[0];
    for(int i=0; i!=trips; ++i) f(t1+i+l1);
  }
};

template<unsigned char R> inline void rcount(int* c, const Box& r)
{
  rcount<R-1>(c,r);
  c[R] = r.bigEnd()[R]-r.smallEnd()[R]+1;
}
template<> inline void rcount<0>(int* c, const Box& r)
{
  c[0]=r.bigEnd()[0]-r.smallEnd()[0]+1;
}
//helper template functions for Grid
template<unsigned char R> inline void gridStride(int* stride, int* lo, const Box& b, const Box& r)
{
  gridStride<R-1>(stride, lo, b, r);
  stride[R] = (b.bigEnd()[R]-b.smallEnd()[R]+1)*stride[R-1];
  lo[R] = r.smallEnd()[R]-b.smallEnd()[R];
}

template<> void inline gridStride<0>(int* stride, int* lo, const Box& b, const Box& r)
{
  stride[0] = b.bigEnd()[0]-b.smallEnd()[0]+1;
  lo[0] = r.smallEnd()[0]-b.smallEnd()[0];
}

/// This is the user template meta program function, everything above here is to support this function
/**

Example usage of BFM template meta-program.  the performCopy operation that
I made for CFArrayBox in Chombo 3

void CFArrayBox::performCopy(const BaseFab<Complex>& a_src,
                             const Box&           a_srcbox,
                             int                  a_srccomp,
                             const Box&           a_destbox,
                             int                  a_destcomp,
                             int                  a_numcomp)
{
  //BaseFab<Complex>::performCopy(a_src, a_srcbox, a_srccomp, a_destbox, a_destcomp, a_numcomp);
  Box r = a_srcbox&a_destbox;
  BFM<CH_SPACEDIM>(dataPtr(),a_src.dataPtr(), box(), a_src.box(), a_destcomp, a_srccomp,a_numcomp,r,
                   [](Complex* __restrict__ d, Complex s){(*d)=s;});
}
*/
template<unsigned char D, typename T1, typename T2, typename F>
void BFM(T1* __restrict t1, const T2* __restrict t2, const Box& b1, const Box& b2, int c1, 
         int c2, int ncomp, const Box& r, F f)
{
  CH_assert(b1.contains(r));
  CH_assert(b2.contains(r));
  int lo1[D+1], lo2[D+1]; lo1[D]=c1; lo2[D]=c2;
  int stride1[D], stride2[D];
  int counting[D+1]; counting[D]=ncomp;
  rcount<D-1>(counting,r);
  gridStride<D-1>(stride1, lo1, b1, r);
  gridStride<D-1>(stride2, lo2, b2, r);

  BFMImpl<D,T1,T2,F>::f(t1,t2,lo1,lo2,stride1,stride2,counting,f);
}

template<unsigned char D, typename T1, typename F>
void BFM1(T1* __restrict t1, const Box& b1, int c1, 
         int ncomp, const Box& r, F f)
{
  CH_assert(b1.contains(r));
  int lo1[D+1]; lo1[D]=c1;
  int stride1[D];
  int counting[D+1]; counting[D]=ncomp;
  rcount<D-1>(counting,r);
  gridStride<D-1>(stride1, lo1, b1, r);

  BFMImpl1<D,T1,F>::f(t1,lo1,stride1,counting,f);
}




#include "Chombo_NamespaceFooter.H"
#endif
