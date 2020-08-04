
namespace Proto
{
  void LSquares::LeastSquares(double** A,
                              vector<double>&x,
                              const vector<double>&rhs)

  {

    int numColsA = x.size();
    int numRowsA = rhs.size();

    double** Atrans;
    allocArray(numColsA,numRowsA,Atrans);

    double** LS; // LS = least squares matrix; Ax = rhs is an overdetermined system
    allocArray(numColsA,numColsA,LS);

    transpose(A,Atrans,numRowsA,numColsA);

    matMul(Atrans,A,LS,numColsA,numRowsA,numColsA);

    // pattern is (A,B,A*B,numrowsA,numcolsA,numcolsB)
    vector<double> ATrhs(numColsA);
    AtimesX(Atrans,rhs,numColsA,ATrhs); // numCols = number of rows of Atranspose*A

    gaussElim(LS,ATrhs);
    // char* bug = "LS";
    // output(numColsA,numColsA,LS,bug);

    backSolve(LS,ATrhs,numColsA,x); // x is the answer

    freeArray(numColsA,numRowsA,Atrans);
    freeArray(numColsA,numColsA,LS);
  }

  void LSquares::AtimesX(double** A,
                         const vector<double>& x,
                         const int& numRowsA,
                         vector<double>& Ax)
  {
    for (int i = 0; i < numRowsA; i++)
    {
      double* scanA = A[i];

      double sum = 0.0;
      int xSize = x.size();

      for (int j = 0; j < xSize; j++)
      {
        sum += *(scanA++) * x[j];
      }

      Ax[i] = sum;
    }
  }

// this is written for square A
  int  LSquares::gaussElim(double** A,
                           vector<double>& rhs)
  {
    //   char* name = "A";
    int currRow = 0;
    int numRows = rhs.size();
    int numCols = rhs.size();

    for (int currCol = 0; currCol < numCols; ++currCol)
    {
      int pivot;

      findPivot(A,currCol, currRow,numRows,pivot);
      // output(numRows,numCols,A,name);
      //
      if (std::abs(A[pivot][currCol]) > 1.0e-15)
      {
        swapRows(A,currRow,pivot,numCols);
        swapRows(rhs,currRow,pivot);

        double Beta = A[currRow][currCol];
        Beta = 1.0/Beta;

        timesBeta(A,currRow,Beta, numCols);
        timesBeta(rhs,currRow,Beta);

        // output(numRows,numCols,A,name);
        for (int rows = currRow; rows < numRows-1; ++rows)
        {
          double alpha = -A[rows+1][currCol];

          addRows(A,rows+1,alpha,currRow,numCols);
          addRows(rhs,rows+1,alpha,currRow);
          // output(numRows,numCols,A,name);
        }
        currRow += 1;
      }
      else
      {
        // MayDay::Warning("small pivot in gaussElim");
      }
    }

    return 0;
  }

  void LSquares::swapRows(vector<double>& rhs,
                          const int& currRow,
                          const int& pivot)
  {
    double    temp = rhs[currRow];
    rhs[currRow] = rhs[pivot];
    rhs[pivot]   = temp;
  }

  void LSquares::swapRows(double** A,
                          const int& rowi,
                          const int& rowj,
                          const int& numCols)
  {
    double *scani = A[rowi];
    double *scanj = A[rowj];

    for (int count = 0; count < numCols; ++count)
    {
      double temp;

      temp     = (*scani);
      (*scani) = (*scanj);
      (*scanj) = temp;

      scani++;
      scanj++;
    }
  }

  int LSquares::findPivot(double** A,
                          const int& currCol,
                          const int& currRow,
                          const int& numRows,
                          int& pivot)
  {
    double max = 0;
    pivot = currRow;

    for (int count = currRow; count < numRows; count++)
    {
      if (std::abs(A[count][currCol]) > max)
      {
        max = std::abs(A[count][currCol]);
        pivot = count;
      }
    }

    return 0;
  }

  void LSquares::addRows(vector<double>& rhs,
                         const int& rowi,
                         const double& alpha,
                         const int& rowj)
  {
    rhs[rowi] += alpha * rhs[rowj];
  }

  void LSquares::addRows(double** A,
                         const int& rowi,
                         const double& alpha,
                         const int& rowj,
                         const int& numCols) // rowi += alpha * rowj
  {
    double *scani = A[rowi];
    double *scanj = A[rowj];

    for (int count = 0; count < numCols; ++count)
    {
      (*scani) += alpha * (*scanj);
      scani++;
      scanj++;
    }
  }

  void LSquares::timesBeta(vector<double>&rhs,
                           const int& currRow,
                           const double& Beta)
  {
    rhs[currRow] *= Beta;
  }

  void  LSquares::timesBeta(double** A,
                            const int& rowi,
                            const double& Beta,
                            const int& numCols)
  {
    double *scanA = A[rowi];

    for (int count = 0; count < numCols; ++count)
    {
      (*scanA) *= Beta;
      scanA++;
    }
  }

  void LSquares::transpose(double** a_A,
                           double** a_Atrans,
                           const int& a_numRowsA,
                           const int& a_numColsA)
  {
    for (int irow = 0; irow < a_numColsA; ++irow)
    {
      double *scanAtrans = a_Atrans[irow];

      for (int icol = 0; icol < a_numRowsA; ++icol)
      {
        (*scanAtrans) = a_A[icol][irow];
        scanAtrans++;
      }
    }
  }

  void LSquares::matMul(double** a_A,
                        double** a_B,
                        double** a_C,
                        const int& a_numRowsA,
                        const int& a_numColsA,
                        const int& a_numColsB)
  {
    for (int i = 0; i < a_numRowsA; ++i)
    {
      double* scanC = a_C[i];

      for (int j = 0; j < a_numColsB; ++j)
      {
        double* scanA = a_A[i];
        double sum = 0.0;

        for (int k = 0; k < a_numColsA; ++k)
        {
          sum += (*scanA) * a_B[k][j];
          scanA++;
        }

        (*scanC) = sum;
        scanC++;
      }
    }
  }

  void LSquares::backSolve(double** a_A,
                           const vector <double>& a_rhs,
                           const int& a_numArows,
                           vector<double>& a_x)
  {
    int N = a_numArows;
    for (int n = 0; n < N; ++n)
    {
      a_x[n] = 0.0;
    }

    for (int n = N-1; n >= 0; --n)
    {
      for (int m = 1; m < N-n; ++m)
      {
        a_x[n] -= a_x[n+m] * a_A[n][n+m];
      }
      if (std::abs(a_A[n][n]) > 0.0)
      {
        a_x[n] += a_rhs[n] / a_A[n][n]; // this only works for square A
      }
    }
  }

  void LSquares::allocArray(const int& rows,
                            const int& cols,
                            double**& A)
  {
    A = new double* [rows];

    for (int i = 0; i < rows;i++)
    {
      A[i] = new double [cols];
      double* scanA = A[i];

      for (int j = 0; j < cols; j++)
      {
        *(scanA++) = 0.0;
      }
    }
  }

  void LSquares::freeArray(const int& rows,
                           const int& cols,
                           double**& A)
  {
    for (int i = 0; i < rows; i++)
    {
      delete[] A[i];
    }

    delete[] A;
  }

}
