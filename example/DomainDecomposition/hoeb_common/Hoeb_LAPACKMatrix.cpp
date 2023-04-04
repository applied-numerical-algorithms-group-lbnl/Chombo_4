#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "Hoeb_LAPACKMatrix.H"
#include <iomanip>
#include "Chombo_parstream.H"
#include "Hoeb_Lapack.H"

#include <cstddef> 
#include "Chombo_CH_Timer.H"
#include <cmath>
#include <istream>
#include <cblas.h>
namespace hoeb
{

  bool LAPACKMatrix::s_checkConditionNumber = false;
  bool LAPACKMatrix::s_verbose = false;
  bool LAPACKMatrix::s_outputStenData = false;


  ostream&
  operator<< (ostream&            a_os,
              const LAPACKMatrix& a_mat)
  {
    a_os << a_mat.m_nrow << " ";
    a_os << a_mat.m_ncol << endl;
    for(int irow =0; irow < a_mat.m_nrow; irow++)
    {
      for(int icol = 0; icol < a_mat.m_ncol; icol++)
      {

        Real val = a_mat(irow, icol);
        a_os      << setprecision(10)
                  << setiosflags(ios::showpoint)
                  << setiosflags(ios::scientific);
        if(val >= 0.) a_os << " ";
        a_os << val;
        a_os << " ";
      }
      a_os << endl;
    }
    if (a_os.fail())
      Chombo4::MayDay::Error("operator<<(ostream&,LAPACKMatrix&) failed");
    return a_os;
  }


  std::istream&
  operator>> (std::istream      & a_is,
              LAPACKMatrix      & a_mat)
  {
    int nrow, ncol;
    a_is >> nrow >> ncol;
    a_mat.define(nrow, ncol);
    for(int irow =0; irow < a_mat.m_nrow; irow++)
    {
      for(int icol = 0; icol < a_mat.m_ncol; icol++)
      {
        Real val;
        a_is >> val;
        a_mat(irow, icol)= val;
      }
    }
    return a_is;
  }
  LAPACKMatrix::
  LAPACKMatrix(const LAPACKMatrix& a_input):
    m_nrow(a_input.m_nrow),
    m_ncol(a_input.m_ncol),
    m_data(new Real[a_input.m_nrow*a_input.m_ncol])
  {
    *this=a_input;
  }

///
  Real 
  LAPACKMatrix::
  normLTwo() const
  {
    CH_TIME("LAPACKMatrix::normLTwo");
    Real retval = 0;

    for(int irow = 0; irow < m_nrow; irow++)
    {
      for(int icol = 0; icol < m_ncol; icol++)
      {
        Real val = (*this)(irow, icol);
        retval += val*val;
      }
    }  
    retval = sqrt(retval);
    return retval;
  }
///
  void
  LAPACKMatrix::
  setVal(const Real& a_val)
  {
    CH_TIME("LAPACKMatrix::setval");
    for(int irow = 0; irow < m_nrow; irow++)
    {
      for(int icol = 0; icol < m_ncol; icol++)
      {
        (*this)(irow, icol) = a_val;
      }
    }
  }

///
  LAPACKMatrix::
  ~LAPACKMatrix()
  {
    if(!m_alias && m_data) delete[] m_data;
  }
///
  int 
  LAPACKMatrix::offset(int a_irow, int a_icol) const
  {
    CH_assert((a_irow >= 0) && (a_irow < m_nrow));
    CH_assert((a_icol >= 0) && (a_icol < m_ncol));

    int retval = a_irow + a_icol*m_nrow;
    return retval;
  }

  const Real& 
  LAPACKMatrix::
  operator() (int a_irow, int a_icol) const
  {
    int ioffset = offset(a_irow, a_icol);
    return m_data[ioffset];
  }
///
  Real& 
  LAPACKMatrix::
  operator() (int a_irow, int a_icol)
  {
    int ioffset = offset(a_irow, a_icol);
    return m_data[ioffset];
  }
///
  Real* 
  LAPACKMatrix::
  dataPtr()
  {
    return m_data;
  }
///
  const Real* 
  LAPACKMatrix::
  dataPtr() const
  {
    return m_data;
  }
///
  Real
  LAPACKMatrix::
  maxNorm() const
  {
    CH_TIME("LAPACKMatrix::maxNorm");
    Real maxval = 0;
    for(int irow =0; irow < m_nrow; irow++)
    {
      for(int icol = 0; icol < m_ncol; icol++)
      {
        maxval = std::max(maxval, (*this)(irow, icol));
      }
    }
    return maxval;
  }
///
  void 
  LAPACKMatrix::
  poutAll() const
  {
    for(int irow =0; irow < m_nrow; irow++)
    {
      Chombo4::pout() << "(";
      if(irow < 10) Chombo4::pout() << " ";
      Chombo4::pout() << irow  <<  "): ";
      for(int icol = 0; icol < m_ncol; icol++)
      {
        //          Chombo4::pout() << "(" << irow  << "," << icol << ")";
        Real val = (*this)(irow, icol);
        Chombo4::pout()    << setprecision(4)
                  << setiosflags(ios::showpoint)
                  << setiosflags(ios::scientific);
        if(val >= 0.) Chombo4::pout() << " ";
        Chombo4::pout() << val;
        Chombo4::pout() << " ";
      }
      Chombo4::pout() << endl;
    }
  }
  void 
  LAPACKMatrix::
  poutMaxMins() const
  {
    Real maxDiag = -1.0e10;
    Real minDiag =  1.0e10;
    Real maxOffd = -1.0e10;
    Real minOffd =  1.0e10;
    for(int irow =0; irow < m_nrow; irow++)
    {
      for(int icol = 0; icol < m_ncol; icol++)
      {
        Real val = (*this)(irow, icol);
        if(irow == icol)
        {
          if(val > maxDiag)
          {
            maxDiag = val;
          }
          if(val < maxDiag)
          {
            minDiag = val;
          }
        }
        else
        {
          if(val > maxOffd)
          {
            maxOffd = val;
          }
          if(val < minOffd)
          {
            minOffd = val;
          }
        }
      }
    }
    Chombo4::pout() << "max value ON  diagonal = " << maxDiag;
    Chombo4::pout() << "min value ON  diagonal = " << minDiag;
    Chombo4::pout() << "max value OFF diagonal = " << maxOffd;
    Chombo4::pout() << "min value OFF diagonal = " << minOffd;
  }
///
  void
  LAPACKMatrix::
  poutMatlab() const
  {
    for(int irow =0; irow < m_nrow; irow++)
    {
      for(int icol = 0; icol < m_ncol; icol++)
      {
        Real val = (*this)(irow, icol);
        Chombo4::pout()    << setprecision(16)
                  << setiosflags(ios::showpoint)
                  << setiosflags(ios::scientific);
        Chombo4::pout() << val;
        Chombo4::pout() << " ";
      }
      Chombo4::pout() << endl;
    }
  }
///
  void 
  LAPACKMatrix::
  poutDiag() const
  {
    for(int irow =0; irow < m_nrow; irow++)
    {
      if(irow < m_ncol)
      {
        Chombo4::pout() << "(" << irow  <<  "): ";
        Real val = (*this)(irow, irow);
        Chombo4::pout()    << setprecision(4)
                  << setiosflags(ios::showpoint)
                  << setiosflags(ios::scientific);
        Chombo4::pout() << val;
        Chombo4::pout() << endl;
      }
    }
  }
///
  void
  LAPACKMatrix::
  poutDiagMatlab() const
  {
    for(int irow =0; irow < m_nrow; irow++)
    {
      if(irow < m_ncol)
      {
        Real val = (*this)(irow, irow);
        Chombo4::pout()    << setprecision(16)
                  << setiosflags(ios::showpoint)
                  << setiosflags(ios::scientific);
        Chombo4::pout() << val;
        Chombo4::pout() << endl;
      }
    }
  }
///
  LAPACKMatrix& 
  LAPACKMatrix::operator=(const LAPACKMatrix& a_matrix)
  {
    CH_TIME("LAPACKMatrix::operator=");
    if(&a_matrix != this)
    {
      if(m_nrow*m_ncol < a_matrix.m_nrow*a_matrix.m_ncol)
      {
        // existing storage not big enough. need to delete and make a new one
        if(!m_alias && m_data) delete[] m_data;
        m_data = new Real[ a_matrix.m_nrow*a_matrix.m_ncol ];
        m_alias = false;
      }
      m_nrow=a_matrix.m_nrow;
      m_ncol=a_matrix.m_ncol;
      for(int irow =0; irow < m_nrow; irow++)
      {
        for(int icol = 0; icol < m_ncol; icol++)
        {
          (*this)(irow, icol)  = a_matrix(irow, icol);
        }
      }
    }
    return *this;
  }

///
  LAPACKMatrix& 
  LAPACKMatrix::operator+=(const LAPACKMatrix& a_matrix)
  {
    CH_TIME("LAPACKMatrix::operator+=");
    for(int irow =0; irow < m_nrow; irow++)
    {
      for(int icol = 0; icol < m_ncol; icol++)
      {
        (*this)(irow, icol)  += a_matrix(irow, icol);
      }
    }
    return *this;
  }

///
  LAPACKMatrix& 
  LAPACKMatrix::operator-=(const LAPACKMatrix& a_matrix)
  {
    CH_TIME("LAPACKMatrix::operator-=");
    for(int irow =0; irow < m_nrow; irow++)
    {
      for(int icol = 0; icol < m_ncol; icol++)
      {
        (*this)(irow, icol)  -= a_matrix(irow, icol);
      }
    }
    return *this;
  }

///
  LAPACKMatrix& 
  LAPACKMatrix::operator*=(const Real& a_scalingFactor)
  {
    CH_TIME("LAPACKMatrix::operator*=");
    for(int irow =0; irow < m_nrow; irow++)
    {
      for(int icol = 0; icol < m_ncol; icol++)
      {
        (*this)(irow, icol)  *= a_scalingFactor;
      }
    }
    return *this;
  }
///
  void
  LAPACKMatrix::setToIdentity()
  {
    for(int irow =0; irow < m_nrow; irow++)
    {
      for(int icol = 0; icol < m_ncol; icol++)
      {
        Real value = 0;
        if(irow == icol) 
        {
          value = 1;
        }
        (*this)(irow, icol)  = value;
      }
    }
  }
///
  void
  LAPACKMatrix::transpose()
  {
    CH_TIME("LAPACKMatrix::transpose");
    LAPACKMatrix copy(*this);
    m_nrow = copy.m_ncol;
    m_ncol = copy.m_nrow;
    for(int irow =0; irow < m_nrow; irow++)
    {
      for(int icol = 0; icol < m_ncol; icol++)
      {
        (*this)(irow, icol)  = copy(icol, irow);
      }
    }
  }
////
  void 
  multiply(LAPACKMatrix& a_product, 
           const LAPACKMatrix& a_left,
           const LAPACKMatrix& a_right)
  {
    CH_TIME("LAPACKMatrix::multiply");
    if(a_left.m_ncol != a_right.m_nrow)
    {
      Chombo4::MayDay::Error("LAPACKMatrix multiply needs colums of left to match the rows on right");
    }

    int nr1 =  a_left.m_nrow;
    int nc2 = a_right.m_ncol;
    int nc1 =  a_left.m_ncol; //== a_right.m_nrow
    // clear existing m_data and allocate new (or use a_productData if non-NULL)
    a_product = LAPACKMatrix(nr1, nc2);

    // lapack_gemm sets C = alpha*A*B + beta*C
    char TRANS = 'N'; // no transposes
    Real ALPHA = 1.;
    Real BETA = 0.;
    // lapack_gemm(&TRANS, &TRANS, &nr1, &nc2, &nc1, &ALPHA,
    //             a_left.dataPtr(), &nr1,
    //             a_right.dataPtr(), &nc1, &BETA,
    //             a_product.dataPtr(), &nr1);
    HOEB_LAPACK(GEMM,gemm)(&TRANS, &TRANS, &nr1, &nc2, &nc1, &ALPHA,
                           (Real*)a_left.dataPtr(), &nr1,
                           (Real*)a_right.dataPtr(), &nc1, &BETA,
                           a_product.dataPtr(), &nr1);
    // for(int irow = 0; irow < nr1; irow++)
    //   {
    //     for(int icol = 0; icol < nc2; icol++)
    //       {
    //         Real value = 0;
    //         for(int jvec = 0; jvec < nc1; jvec++)
    //           {
    //             value += a_left(irow, jvec)*a_right(jvec, icol);
    //           }
    //         a_product(irow, icol) = value;
    //       }
    //   }
  }
/// returns the info value 
  int
  LAPACKMatrix::
  invert()
  {
    CH_TIME("LAPACKMatrix::invert");
    int INFO = 0;
    //only makes sense if square
    if(m_nrow != m_ncol)
    {
      Chombo4::MayDay::Error("trying to invert a non-square matrix");
    }

    Real* A = m_data;
    int   N = m_ncol; //==m_nrow

    int *IPIV = new int[N+1];
    int LWORK = N*N;
    Real *WORK = new Real[LWORK];


    HOEB_LAPACK(GETRF,getrf)(&N,&N,A,&N,IPIV,&INFO);
    HOEB_LAPACK(GETRI,getri)(&N,A,&N,IPIV,WORK,&LWORK,&INFO);
    if(INFO != 0)
    {
      Chombo4::MayDay::Warning(" info flag from lapack");
      Chombo4::pout() << "getri matrix may be singular---info = " << INFO << endl;
    }
    else if(s_checkConditionNumber)
    {
      checkConditionNumber();
    }
  
    delete[] IPIV;
    delete[] WORK;

    return INFO;

  }
///
// LAPACKMatrix::
// LAPACKMatrix( int a_nrow, int a_ncol, const Real* const a_data)
// {
//   CH_TIME("Matrix::define1");
//   setDefaultValues();
//   define(a_nrow, a_ncol); 
//   for(int irow = 0; irow < m_nrow; irow++)
//     {
//       for(int icol = 0; icol < m_ncol; icol++)
//         {
//           int ioff = offset(irow, icol);
//           (*this)(irow, icol) = *(a_data + ioff);
//         }
//     }  
// }

///
  void
  LAPACKMatrix::
  checkConditionNumber() const
  {
    Real inverse = getInverseOfConditionNumber(*this);
    Real small = 1.0e-6;
    Real reallysmall = 1.0e-15;
    if(inverse < reallysmall)
    {
      Chombo4::pout() << "matrix is poorly conditioned: 1/condition_number = " << inverse << endl;
      if(s_verbose)
      {
        this->poutAll();
      }
    }
    else if(inverse < small)
    {
      Chombo4::pout() << "matrix is marginally conditioned: 1/condition_number = " << inverse << endl;
    }
    else
    {
      Chombo4::pout() << "matrix might be OK: 1/condition_number = " << inverse << endl;
    }
  }
///
  void
  LAPACKMatrix::
  checkUpperTriangularConditionNumber() const
  {
    Real inverse = getInverseOfUpperTriangularConditionNumber(*this);
    Real small = 1.0e-6;
    Real reallysmall = 1.0e-15;
    if(inverse < reallysmall)
    {
      Chombo4::pout() << "matrix is poorly conditioned: 1/condition_number = " << inverse << endl;
      if(s_verbose)
      {
        this->poutAll();
      }
    }
    else if(inverse < small)
    {
      Chombo4::pout() << "matrix is marginally conditioned: 1/condition_number = " << inverse << endl;
    }
    else
    {
      // Chombo4::pout() << "matrix might be OK: 1/condition_number = " << inverse << endl;
    }
  }
/**
 *  Get Ainverse using  least squares with svd
 */
  int
  LAPACKMatrix::invertUsingSVD(int a_maxiter, Real a_tol)
  {
    CH_TIME("LAPACKMatrix::invertUsingSVD");
    int M = m_nrow;
    int N = m_ncol;

    LAPACKMatrix rhs(M, N);
    rhs.setToIdentity();
    LAPACKMatrix Acopy = *this;

    //  Acopy.transpose();
    int retval = solveLSTSVD(Acopy, rhs, a_maxiter, a_tol);
    //int retval = solveLSTSVDOnce(Acopy, rhs);
    *this = rhs;
    return retval;
  }
/***/
  int 
  LAPACKMatrix::
  pseudoInvertUsingSVD(int a_maxiter, Real a_tol)
  {
    int retval = invertUsingSVD(a_maxiter, a_tol);
  
    return retval; 
  }
/***/
  int
  LAPACKMatrix::invertUsingLeastSquares()
  {
    CH_TIME("LAPACKMatrix::pseudoInvertUsingLeastSquares");
    int M = m_nrow;
    int N = m_ncol;

    LAPACKMatrix rhs(M, N);
    rhs.setToIdentity();
    LAPACKMatrix Acopy = *this;
    int retval = 0;
    if(M== N)
    {
      retval = solveLeastSquares(Acopy, rhs);
    }
    else
    {
      CH_assert(M>N);
      for(int irow = N; irow < M; irow++)
      {
        rhs(irow,N-1) = 1;
      }
      retval = solveReducedRankLS(Acopy, rhs);
    }
    *this = rhs;
    return retval;
  }

/***/
  int
  LAPACKMatrix::
  pseudoInvertUsingQR()
  {
    CH_TIME("LAPACKMatrix::pseudoInvertUsingQR");
    int M = m_nrow; // number of equations
    int N = m_ncol; // number of variables
    CH_assert(M >= N);
    LAPACKMatrix TAU(N, 1);
    int LWORK = 2*M*N; // allowing block size of 2*M, which is plenty
    LAPACKMatrix WORK(LWORK, 1);
    int INFO;
    HOEB_LAPACK(GEQRF,geqrf)(&M, &N, dataPtr(), &M,
                             TAU.dataPtr(), WORK.dataPtr(), &LWORK, &INFO);
    if (INFO != 0)
    {
      Chombo4::MayDay::Warning(" info flag from lapack");
      Chombo4::pout() << "geqrf call has bad value at argument " << (-INFO) << endl;
    }

    if(LAPACKMatrix::s_verbose)
    {
      checkUpperTriangularConditionNumber();
    }

    LAPACKMatrix Xt(*this);

    // Generate the M by N matrix Q, and store it in Xt.
    HOEB_LAPACK(ORGQR,orgqr)(&M, &N, &N, Xt.dataPtr(),
                             &M, TAU.dataPtr(), WORK.dataPtr(), &LWORK, &INFO);
    if (INFO != 0)
    {
      Chombo4::MayDay::Warning(" info flag from lapack");
      Chombo4::pout() << "orgqr call has bad value at argument " << (-INFO) << endl;
    }

    // Solve for X^T in X^T * R^T = Q.
    char SIDE = 'R';
    char UPLO = 'U';
    char TRANSA = 'T';
    char DIAG = 'N';
    Real ALPHA = 1.;
    HOEB_LAPACK(TRSM,trsm)(&SIDE, &UPLO, &TRANSA, &DIAG, &M, &N, &ALPHA,
                           dataPtr(), &M, Xt.dataPtr(), &M);
  
    *this = Xt;
    transpose();

    return INFO;

  }
  void
  LAPACKMatrix::truncate(int a_nrow,int a_ncol)
  {
    CH_assert(a_nrow <= m_nrow);
    CH_assert(a_ncol <= m_ncol);
    CH_assert(a_nrow > 0);
    CH_assert(a_ncol > 0);
    m_nrow = a_nrow;
    m_ncol = a_ncol;
  }

///
  Real getInverseOfConditionNumber(const LAPACKMatrix& A)
  {
    CH_TIME("LAPACKMatrix::getInverseOfConditionNumber");
    // 1 or 0 gets L-1  norm.  I is infinity norm
    char NORM = 'I';
    int M = A.m_nrow;
    int N = A.m_ncol;
    int LDA = M;


    LAPACKMatrix tempA;
    if(M==N)
    {
      tempA = A;
    }
    else
    {
      //send in AAT to make sure matrix is square
      LAPACKMatrix MatM = A;
      LAPACKMatrix MatMT= A;
      MatMT.transpose();
      multiply(tempA, MatM, MatMT);
      M = tempA.m_nrow;
      N = tempA.m_ncol;
    }

  
    Real rcond;
    LAPACKMatrix WORK(2*M*N, 1);
    WORK.setVal(0.);
    int INFO;

    int *IPIV = new int[N+1];
    HOEB_LAPACK(GETRF,getrf)(&N,&N,tempA.m_data,&N,IPIV,&INFO);
  
    Real ANORM = 0;
    for(int icol = 0; icol < N; icol++)
    {
      Real colsum = 0;
      for(int irow = 0; irow < M; irow++)
      {
        colsum += std::abs(tempA(irow, icol));
      }
      ANORM = std::max(ANORM, colsum);
    }
      
    const int minMN = 20;
    int IWORK[25*minMN];
    HOEB_LAPACK(GECON,gecon)(&NORM, &N, tempA.m_data, &LDA, &ANORM, &rcond, WORK.dataPtr(), IWORK, &INFO);

    delete[] IPIV;

    return rcond;
  }

///
  Real getInverseOfUpperTriangularConditionNumber(const LAPACKMatrix& A)
  {
    CH_TIME("LAPACKMatrix::getInverseOfUpperTriangularConditionNumber");
    // 1 or 0 gets L-1  norm.  I is infinity norm
    char NORM = '1';
    char UPLO = 'U';
    char DIAG = 'N';
    int M = A.m_nrow;
    int N = A.m_ncol;

    Real rcond;
    LAPACKMatrix WORK(3*N, 1);
    int INFO;

    int *IPIV = new int[N];
    HOEB_LAPACK(TRCON,trcon)(&NORM, &UPLO, &DIAG, &N, A.dataPtr(), &M,
                             &rcond, WORK.dataPtr(), IPIV, &INFO);
    if (INFO != 0)
    {
      Chombo4::MayDay::Warning(" info flag from lapack");
      Chombo4::pout() <<  "trcon call has bad value at argument " << (-INFO) << endl;
    }

    delete[] IPIV;

    return rcond;
  }

///
/**
 *  Solves A*X = B using general least squares, for each column of B
 Answer goes 
 back into B I think
*/
  int solveLeastSquares(LAPACKMatrix& A, LAPACKMatrix& B)
  {
    CH_TIME("LAPACKMatrix::solveLeastSquares");
    // TODO - check that the sizes of A, B and C are compatible
    int M = A.m_nrow;
    int N = A.m_ncol;
    int NRHS = B.m_ncol;
    int LDA = M;
    int LDB = std::max(M,N);
    CH_assert(B.m_nrow == M);

    int LWORK[2] = {1,1};
    LWORK[0] = 2*M*N;

    LAPACKMatrix WORK(2*M*N, 1);
    WORK.setVal(0.);
  
    char TRANS = 'N';
    int INFO;

    HOEB_LAPACK(GELS,gels)(&TRANS, &M, &N, &NRHS, A.dataPtr(), &LDA, 
                           B.dataPtr(), &LDB, WORK.dataPtr(), LWORK, &INFO);

    if(INFO != 0)
    {
      Chombo4::MayDay::Warning(" info flag from lapack");
      Chombo4::pout() << "gels matrix may be singular---info = " << INFO << endl;
    } 
    else if(LAPACKMatrix::s_checkConditionNumber)
    {
      A.checkConditionNumber();
    }

    return INFO;

  }


/**
 *  Solves A'*X = B using least squares, for vector b
 Answer goes 
 back into B I think
*/
  int solveLeastSquaresTranspose(LAPACKMatrix& A, LAPACKMatrix& B)
  {
    CH_TIME("LAPACKMatrix::solveLeastSquaresTranspose");
    // TODO - check that the sizes of A, B and C are compatible
    int M = A.m_nrow;
    int N = A.m_ncol;
    int NRHS = B.m_ncol;
    int LDA = A.m_nrow;
    int LDB = std::max(M, N);
    int LWORK[2] = {1,1};
    LWORK[0] = 2*M*N;
    LAPACKMatrix WORK(2*M, N);
    WORK.setVal(0.);

    char TRANS = 'T';
    int INFO;

    HOEB_LAPACK(GELS,gels)(&TRANS, &M, &N, &NRHS, A.dataPtr(), &LDA, 
                           B.dataPtr(), &LDB, WORK.dataPtr(), LWORK, &INFO);

    if(INFO != 0)
    {
      Chombo4::MayDay::Warning(" info flag from lapack");
      Chombo4::pout() << "gels matrix may be singular---info = " << INFO << endl;
    }
    else if(LAPACKMatrix::s_checkConditionNumber)
    {
      A.checkConditionNumber();
    }
  
    return INFO;

  }


/**
 *  Solves A^T X = B using least squares with SVD, for vector b
 */
  int solveLSTSVD(LAPACKMatrix& A, LAPACKMatrix& B, int a_maxiter, Real a_tol)
  {
    CH_TIME("LAPACKMatrix::solveLSTSVD");
    int INFO = 0;
    Real residinit = B.maxNorm();
    Real tol = a_tol*residinit;
    Real residnorm = residinit;

    int iter = 0;
    LAPACKMatrix resid = B;
    int nrows = A.m_nrow;
    int ncols = A.m_ncol;
    //looks backward because this gets multiplied by A
    LAPACKMatrix incr(ncols, nrows);
    LAPACKMatrix    X(ncols, nrows);
    incr.setVal(0);
    X.setVal(0);
    while((iter== 0)|| ((residnorm >= tol) && (iter < a_maxiter) && (residnorm > a_tol)))
    {
      //this puts the answer into resid
      INFO = solveLSTSVDOnce(A, resid);
      incr = resid;
      incr.transpose();
      X -= incr;
 
      //recalculate residual = ax - B
      resid.setVal(0.);
      multiply(resid, A, X);
      resid.truncate(B.dims().first, B.dims().second);
      resid -= B;
      residnorm = resid.maxNorm();
      
      incr.setVal(0.);
      iter++;
    }
    if((iter >= a_maxiter) && LAPACKMatrix::s_verbose)
    {
      Chombo4::pout() << "matrix::solve warning" << endl;
      Chombo4::pout() << "matrix: maximum number of iterations reached " << endl;
      Chombo4::pout() << "residnorm = " << residnorm << ", iter = " << iter << ", maxresid = " << tol <<  ", residinit = " << residinit <<  endl;
    }
    //because this is where the answer goes
    B = X;
    return INFO;
  }

/**
 *  Solves A^T X = B using least squares with SVD, for vector b
 */
  int solveLSTSVDOnce(LAPACKMatrix& A, LAPACKMatrix& B)
  {
    CH_TIME("LAPACKMatrix::solveLSTSVDOnce");
    // TODO - check that the sizes of A, B and C are compatible
    int M = A.m_nrow;
    int N = A.m_ncol;
    int NRHS = B.m_ncol;
    int LDA = N;
    int LDB = std::max(M,N);

    // Scratch and output allocated via high-water mark
    const int maxMxN = 4000;
    const int minMN = 20;
    int LWORK[2] = {1,1};
    LWORK[0] = maxMxN;
    LAPACKMatrix WORK(1, maxMxN);
    WORK.setVal(0);
    Real S[minMN];
    int IWORK[25*minMN];

    Real RCOND = -1;
    int INFO;
    int RANK;

    //transpose of A
    LAPACKMatrix At = A;

    At.transpose();
    // - save a copy of the RHS to check the residual
    LAPACKMatrix Bcopy = B;

    HOEB_LAPACK(GELSD,gelsd)(&N, &M, &NRHS, At.dataPtr(), &LDA, 
                             B.dataPtr(), &LDB, S, &RCOND, &RANK,
                             WORK.dataPtr(), LWORK, IWORK, &INFO);

    if(INFO != 0)
    {
      Chombo4::MayDay::Warning(" info flag from lapack");
      Chombo4::pout() << "gelsd matrix may be singular---info = " << INFO << endl;
      return INFO;
    }
    else if(LAPACKMatrix::s_checkConditionNumber)
    {
      A.checkConditionNumber();
    }



    // - double-check the residual, toss a warning if it's non-zero
    WORK.setVal(0);
    char TRANS = 'T';
    Real ALPHA = 1;
    Real BETA = -1;
    int INCX = 1;
    int INCY = 1;
    HOEB_LAPACK(GEMV,gemv)(&TRANS, &M, &N, &ALPHA, A.dataPtr(), &M, 
                           B.dataPtr(), &INCX, &BETA, Bcopy.dataPtr(), &INCY);
  
    //  CH_assert(Bcopy.m_ncol == 1);
    Real resid = 0.0;
    for (int i=0; i < N; ++i)
    {
      Real val = Bcopy(i,0);
      resid += val*val;
    }
    resid = sqrt(resid);
    Real tol = 1e-3;
    if (resid > tol)
    {
      Chombo4::pout() << "solveLSTSVD residual = " << resid << endl;
      Chombo4::pout() << "WARNING: solveLSTSVD residual above tolerance " << tol << endl;
      for (int i=0; i < N; ++i)
      {
        Chombo4::pout() << "r(" << i << ") = " << Bcopy(i,0) << endl;
      }
    }

    return INFO;
  }

///
/**
   Solves A*X = B using least squares with SVD and iterative refinement
*/
  int solveLSTSVD(LAPACKMatrix      & X, 
                  const LAPACKMatrix& A, 
                  const LAPACKMatrix& B, 
                  int           a_maxiter, 
                  Real          a_tol)
  {
    CH_TIME("LAPACKMatrix::solveLSTSVD2");
    int INFO = 0;
    Real residinit = B.maxNorm();
    Real tol = a_tol*residinit;
    /**  if(tol < 1.0e-16)
         {
         tol = 1.0e-16;
         }**/

    Real residnorm = residinit;

    int iter = 0;
    LAPACKMatrix resid = B;
    int ncols = A.m_ncol;
    int nrhs  = B.m_ncol;
    X = LAPACKMatrix(ncols, nrhs);
    // incr is the solution to the residual equation
    LAPACKMatrix incr(ncols, nrhs); 
    X.setVal(0);
    incr.setVal(0);
    while((residnorm >= tol) && (iter < a_maxiter))
    {
      INFO = solveLSTSVDOnce(incr, A, resid);
 
      X += incr;

      //recalculate residual = B - AX
      LAPACKMatrix AX;
      multiply(AX, A, X);
      resid = B;
      resid -= AX;
      residnorm = resid.maxNorm();
      
      iter++;
    }
    if((iter >= a_maxiter) && LAPACKMatrix::s_verbose)
    {
      Chombo4::pout() << "matrix::solve warning" << endl;
      Chombo4::pout() << "matrix: maximum number of iterations (" << a_maxiter << ") reached " << endl;
      Chombo4::pout() << "residnorm = " << residnorm << ", iter = " << iter << ", maxresid = " << tol <<  ", residinit = " << residinit <<  endl;
    }

    return INFO;
  }

/**
 *  Solves A*X = B using least squares with SVD, for X
 */
  int solveLSTSVDOnce(LAPACKMatrix      & X, 
                      const LAPACKMatrix& A, 
                      const LAPACKMatrix& B)
  {
    CH_TIME("LAPACKMatrix::solveLSTSVDOnce2");
    int M = A.m_nrow;
    int N = A.m_ncol;
    int NRHS = B.m_ncol;
    int LDA = M;
    int LDB = std::max(M,N);

    int maxMxN = 40000;

    int minMN = std::min(M,N);
    //int LWORK = maxMxN;
    int LWORK[2] = {1,1};
    LWORK[0] = maxMxN;
    LAPACKMatrix WORK(1, maxMxN);
    WORK.setVal(0);
    Real* S = new Real[minMN];
    int* IWORK = new int[25*minMN];
  

    Real RCOND = -1;
    int INFO;
    int RANK;

    X = LAPACKMatrix(N,NRHS);
    // HOEB_LAPACK(GELSD,gelsd) destroys A, so send in a copy instead
    LAPACKMatrix Acopy = A;
    if(M >= N)
    {
      LAPACKMatrix temp = B;
      HOEB_LAPACK(GELSD,gelsd)(&M, &N, &NRHS, Acopy.dataPtr(), &LDA, 
                               temp.dataPtr(), &LDB, S, &RCOND, &RANK,
                               WORK.dataPtr(), LWORK, IWORK, &INFO);

      for(int irow = 0; irow < N; irow++)
      {
        for(int icol = 0; icol < NRHS; icol++)
        {
          X(irow, icol) = temp(irow, icol);
        }
      }
    }
    else
    {
      LAPACKMatrix temp(N,NRHS);
      temp.setVal(0.0);
      for(int irow = 0; irow < M; irow++)
      {
        for(int icol = 0; icol < NRHS; icol++)
        {
          temp(irow, icol) = B(irow, icol);
        }
      }
      HOEB_LAPACK(GELSD,gelsd)(&M, &N, &NRHS, Acopy.dataPtr(), &LDA, 
                               temp.dataPtr(), &LDB, S, &RCOND, &RANK,
                               WORK.dataPtr(), LWORK, IWORK, &INFO);

      X = temp;
    }


    if(INFO != 0)
    {
      Chombo4::MayDay::Warning(" info flag from lapack");
      Chombo4::pout() <<  "gelsd matrix may be singular---info = " << INFO << endl;
      return INFO;
    }
    else if(LAPACKMatrix::s_checkConditionNumber)
    {
      A.checkConditionNumber();
    }


    // TO DO: put in check on residuals (see LSTSVDOnce with two arguments)
    // - double-check the residual, toss a warning if it's non-zero
    WORK.setVal(0);
    // HOEB_LAPACK(GEMV,gemv) sets Y = alpha*A*X + beta*Y
    char TRANS = 'N'; // compute A*X and not A'*X
    Real ALPHA = 1;
    Real BETA = -1;
    int INCX = 1;
    int INCY = 1;
    LAPACKMatrix Bcopy = B;
    LAPACKMatrix Acast = A; // copy of A that is not const
    HOEB_LAPACK(GEMV,gemv)(&TRANS, &M, &N, &ALPHA, Acast.dataPtr(), &M, 
                           X.dataPtr(), &INCX, &BETA, Bcopy.dataPtr(), &INCY);
  
    //  CH_assert(Bcopy.m_ncol == 1);
    /**Real resid = 0.0;
       for (int i=0; i < N; ++i)
       {
       Real val = Bcopy(i,0);
       resid += val*val;
       }
       resid = sqrt(resid);
       Real tol = 1e-3;
       if (resid > tol)
       {
       Chombo4::pout() << "solveLSTSVD residual = " << resid << endl;
       Chombo4::pout() << "WARNING: solveLSTSVD residual above tolerance " << tol << endl;
       for (int i=0; i < N; ++i)
       {
       Chombo4::pout() << "r(" << i << ") = " << Bcopy(i,0) << endl;
       }
       }**/

    delete[] S;
    delete[] IWORK;
    return INFO;
  }

/**
 *  Solves equality constrained least squares problem
 *    Find x, s.t. min norm(A x - c) with B x = d
 */
  int solveEqualityConstrainedLS(LAPACKMatrix& A, LAPACKMatrix& c, LAPACKMatrix& B, LAPACKMatrix& d, LAPACKMatrix& x)
  {
    // Check that the sizes of A, c, B, d and x, for compatibility
    int M = A.m_nrow; // # rows of A, c
    int N = A.m_ncol; // # cols of A, B, rows of x
    int P = B.m_nrow; // # rows of B, d
    CH_assert(B.m_ncol == N);
    CH_assert(c.m_nrow == M);
    CH_assert(c.m_ncol == 1);
    CH_assert(d.m_nrow == P);
    CH_assert(d.m_ncol == 1);
    CH_assert(x.m_nrow == N);
    CH_assert(x.m_ncol == 1);

    // Make sure the problem is well-posed
    CH_assert(N >= P);
    CH_assert(M+P >= N);

    // Scratch and other LAPACK parameters
    int LDA = M;
    int LDB = P;
    const int maxMxN = 4000;
    int lwork = maxMxN;
    LAPACKMatrix work(lwork,1);
    work.setVal(0);
    int info;
    HOEB_LAPACK(GGLSE,gglse)(&M, &N, &P, A.dataPtr(), &LDA, B.dataPtr(), &LDB, 
                             c.dataPtr(), d.dataPtr(), x.dataPtr(), work.dataPtr(), &lwork, &info);
    int INFO = info;
    if(INFO != 0)
    {
      Chombo4::MayDay::Warning(" info flag from lapack");
      Chombo4::pout() << "gglse matrix may be singular---info = " << INFO << endl;
      return INFO;
    }
    else if(LAPACKMatrix::s_checkConditionNumber)
    {
      A.checkConditionNumber();
    }

    // Calculate the residual
    Real resid = 0.0;
    for (int i=N-P; i < M; ++i)
    {
      Real val = c(i,0);
      resid += val*val;
    }
    resid = sqrt(resid);
    Real tol = 1e-3;
    if (resid > tol)
    {
      Chombo4::pout() << "solveequalityconstrainedLS residual = " << resid << endl;
      Chombo4::pout() << "WARNING: solveLSTSVD residual above tolerance " << tol << endl;
    }
    return INFO;
  }
///
/**
 *  Solves A'*X = B using reduced rank least squares, for vector b
 */
  int solveReducedRankLS(LAPACKMatrix& A, LAPACKMatrix& b)
  {
    // - check that the sizes of A, b are compatible
    int M = A.m_nrow;
    int N = A.m_ncol;
    CH_assert(b.m_nrow == M);
    //  CH_assert(b.m_ncol == 1);
    CH_assert(M >= N); // A is over-determined
    // Chombo4::pout() << "A is " << M << " x " << N << endl;

    // - create a transpose
    LAPACKMatrix At = A;
    At.transpose();

    // - Calculate QR factorization of A' with DGEQRF
    int LWORK[2] = {1,1};
    LWORK[0] = 2*M*N;
    LAPACKMatrix WORK(LWORK[0], 1);
    WORK.setVal(0);

    LAPACKMatrix TAU(N, 1);
    TAU.setVal(0);

    int INFO;
    HOEB_LAPACK(GEQRF,geqrf)(&N, &M, At.dataPtr(), &N,
                             TAU.dataPtr(), WORK.dataPtr(), LWORK, &INFO);

    if(INFO != 0)
    {
      Chombo4::MayDay::Warning(" info flag from lapack");
      Chombo4::pout() << "geqrf matrix may be singular---info = " << INFO << endl;
    }
    else if(LAPACKMatrix::s_checkConditionNumber)
    {
      At.checkConditionNumber();
    }
    CH_assert(INFO == 0);

    char SIDE = 'L';
    char TRANS = 'T';
    int NRHS =  b.m_ncol;
    HOEB_LAPACK(ORMQR,ormqr)(&SIDE, &TRANS, &N, &NRHS, &N, 
                             At.dataPtr(), &N, TAU.dataPtr(), b.dataPtr(), &N,
                             WORK.dataPtr(), LWORK, &INFO);

    if(INFO != 0)
    {
      Chombo4::MayDay::Warning(" info flag from lapack");
      Chombo4::pout() << "ormqr matrix may be singular---info = " << INFO << endl;
      return INFO;
    }

  
    // - Solve R x = (Q' * b) using DTRTRS
    char UPLO = 'U';
    TRANS = 'N';
    char DIAG = 'N';
    HOEB_LAPACK(TRTRS,trtrs)(&UPLO, &TRANS, &DIAG, &N, &NRHS, 
                             At.dataPtr(), &N, b.dataPtr(), &N, &INFO);
    if(INFO != 0)
    {
      Chombo4::MayDay::Warning(" info flag from lapack");
      Chombo4::pout() << "trtrs matrix may be singular---info = " << INFO << endl;
    }
    
    return INFO;
  }
}
