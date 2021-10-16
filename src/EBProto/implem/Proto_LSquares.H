
#ifndef _Proto_LSQUARES_H_
#define _Proto_LSQUARES_H_

using std::vector;

namespace Proto
{
///
/**
   This class solves least squares problems
*/
  class LSquares
  {
  public:
    inline void LeastSquares(double** A,
                             vector<double>&x,
                             const vector<double>&rhs);

    /// gaussian elimination with partial pivoting
    inline int gaussElim(double**A,
                         vector<double>& rhs);

    inline void swapRows(double** A,
                         const int& rowi,
                         const int& rowj,
                         const int&numberOfCols);

    inline void swapRows(vector<double>& rhs,
                         const int& currRow,
                         const int& pivot);

    inline int findPivot(double** A,
                         const int& currCol,
                         const int& currRow,
                         const int& numRows,
                         int& pivot);

    inline void addRows(vector<double>& rhs,
                        const int& rowi,
                        const double& alpha,
                        const int& rowj);

    inline void addRows(double** A,
                        const int& rowi,
                        const double& alpha,
                        const int& rowj,
                        const int& numberOfCols);

    inline void timesBeta(vector<double>&rhs,
                          const int& currRow,
                          const double& Beta);

    inline void timesBeta(double** A,
                          const int& rowi,
                          const double& Beta,
                          const int& numberOfcols);

    inline void transpose(double** a_A,
                          double ** a_Atrans,
                          const int& a_numRows,
                          const int& a_numCols);

    inline void matMul(double** a_A,
                       double** a_B,
                       double** a_C,
                       const int& a_numRowsA,
                       const int& a_numColsA,
                       const int& a_numColsB);

    inline void backSolve(double** a_A,
                          const vector<double>& a_rhs,
                          const int& a_numArows,
                          vector<double>& a_x);

    inline void AtimesX(double** A,
                        const vector<double>&x,
                        const int& numRowsA,
                        vector<double>& Ax);

    inline void allocArray(const int& rows,
                           const int& cols,
                           double**& A);

    inline void freeArray(const int& rows,
                          const int& cols,
                          double**& A);

  };

}

#include "Proto_LSquaresImplem.H"

#endif
