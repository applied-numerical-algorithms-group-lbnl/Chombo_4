#pragma once

namespace Proto
{
    CUDA_DECORATION
    inline bool contains(const Point& pt, Point& a_low, Point& a_high)
    {
      for(int idir = 0; idir < DIM; idir++)
      {
        if(pt[idir] < a_low[idir])
        {
            return false;
        }
        if(pt[idir] > a_high[idir])
        {
            return false;
        }
      }
      return true;
    }


    template<CENTERING cent>
    CUDA_DECORATION
    inline bool contains(EBIndex<cent>& a_in, Point& a_low, Point& a_high)
    {
      auto pt = a_in.m_pt;
      return contains(pt,a_low,a_high);
    }
}


