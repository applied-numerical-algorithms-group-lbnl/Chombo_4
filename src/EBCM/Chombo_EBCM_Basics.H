#ifndef __ebcm__basics__
#define __ebcm__basics__
#include "Proto.H"
#include "Chombo_EBChombo.H"


namespace ebcm_basics
{

  ///
  IndexTM<int, DIM>
  getVolIndexFromFace(const  IndexTM<int, DIM-1>& a_facind, int a_facedir)
  {
    IndexTM<int, DIM> retval;
    int ifac = 0;
    for(int idir = 0; idir < DIM; idir++)
    {
      if(idir != a_facedir)
      {
        retval[idir] = a_facind[ifac];
        ifac++;
      }
      else
      {
        retval[idir] = 0;
      }
    }
    return retval;
  }
  IndexTM<int, DIM-1>
  getFaceIndexFromVolIndex(const  IndexTM<int, DIM>& a_volind,
                           int                       a_facedir)
  {
    IndexTM<int, DIM-1> retval;
    int ifac = 0;
    for(int idir = 0; idir < DIM; idir++)
    {
      if(idir != a_facedir)
      {
        retval[ifac] = a_volind[idir];
        ifac++;
      }
    }
    return retval;
  }

  template <int ebcm_order>
  IndexedMoments<DIM, ebcm_order>
  getVolMomFromFacMom(const IndexedMoments<DIM-1, ebcm_order>& a_facmom,
                      int                                      a_facedir)
  {
    IndexedMoments<DIM, ebcm_order> retval;
    for(MomentIterator< DIM, ebcm_order> momit; momit.ok(); ++momit)
    {
      IndexTM<int, DIM> volpow = momit();
      Real momval;
      if(volpow[a_facedir] > 0) //no derivs in this direction
      {
        momval = 0;
      }
      else
      {
        IndexTM<int, DIM-1> facpow = getFaceIndexFromVolIndex(volpow, a_facedir);
        momval = a_facmom[facpow];
      }
      retval[volpow] = momval;
    }
    return retval;
  }

}

#endif
