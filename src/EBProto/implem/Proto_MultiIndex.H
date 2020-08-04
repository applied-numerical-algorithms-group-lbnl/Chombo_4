
#ifndef _Proto_MULTIINDEX_H_
#define _Proto_MULTIINDEX_H_

#include <map>
using std::map;
using std::vector;


#include "Proto_IndexTM.H"

namespace Proto
{

// Generate and return a vector of all multi-indices whose magnitude (the sum
// of their entries) equals "a_magnitude"
  template <int dim> inline void generateMultiIndices(vector<IndexTM<int,dim> > & a_indices,
                                                      const int                 & a_magnitude)
  {
    if (a_magnitude >= 0)
    {
      IndexTM<int,dim> monomial;

      for (int idir = 0; idir < dim; ++ idir)
      {
        monomial[idir] = 0;
      }

      while (true)
      {
        for (int j = 1; j < dim-1; ++j)
        {
          int t = a_magnitude;
          for (int k = j+1; k < dim; ++k)
          {
            t -= monomial[k];
          }

          if (monomial[j] > t)
          {
            monomial[j+1] += 1;
            monomial[j] = 0;
          }
          else
          {
            break;
          }
        }

        if (monomial[dim-1] > a_magnitude)
        {
          break;
        }

        monomial[0] = a_magnitude;

        for (int j = 1; j < dim; ++j)
        {
          monomial[0] -= monomial[j];
        }

        a_indices.push_back(monomial);

        monomial[1] += 1;
      }
    }
  }

}


#endif
