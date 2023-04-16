#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _EigenMATRIX_H_
#define _EigenMATRIX_H_

#include "Eigen/Core"
#include "Eigen/SVD"
#include "Eigen/LU"
#include <utility>
#include <string>
#include <memory>
namespace ch_eigen
{

    typedef Chombo4::MayDay  ch_mayday;
    typedef Eigen::MatrixXd  eigen_mat;
///
/**
   Very simple wrapper around Eigen::matrix class designed to 
   be a transition from between the LAPACKMatrix class and
   simply using the excellent Eigen infrastructure for our
   linear algebra needs.
*/
  class Matrix
  {
  public:
    
    ///This is the actual matrix.  The class has no other data.
    std::shared_ptr<eigen_mat> m_mat_p;

    Matrix()
    {
    }

    Matrix(int a_nrow, int a_ncol)
    {
      m_mat_p = shared_ptr<eigen_mat>(new eigen_mat(a_nrow, a_ncol));
    }


    ///  deep copy constructor
    Matrix(const Matrix& a_input)
    {
      int nrow = a_input.m_mat_p->rows();
      int ncol = a_input.m_mat_p->cols();
      m_mat_p = shared_ptr<eigen_mat>(new eigen_mat(nrow, ncol));
      
      (*m_mat_p) =  (*a_input.m_mat_p);
    }


    ///returns the L-2 norm (what Eigen calls the Euclidean norm).
    inline double norm() const
    {
      return m_mat_p->norm();
    }
  
    /// Eigen probably has way faster ways to do this 
    inline void setToIdentity()
    {
      
      for(int irow = 0; irow < m_mat_p->rows(); irow++)
      {
        for(int icol = 0; icol < m_mat_p->cols(); icol++)
        {
          double val = 0;
          if(irow == icol)
          {
            val = 1;
          }
          (*m_mat_p)(irow, icol) = val;
        }
      }
    }
  

    /// Eigen probably has way faster ways to do this 
    inline void  setVal(const double& a_val)
    {
      for(int irow = 0; irow < m_mat_p->rows(); irow++)
      {
        for(int icol = 0; icol < m_mat_p->cols(); icol++)
        {
          (*m_mat_p)(irow, icol) = a_val;
        }
      }
    }

    ///
    const inline double& operator() (int irow, int icol) const
    {
      return (*m_mat_p)(irow, icol);
    }
  
    ///
    inline double& operator() (int irow, int icol)
    {
      return (*m_mat_p)(irow, icol);
    }
      

    ///return pair <nrow, ncol>
    inline std::pair<int, int>  dims() const
    {
    
      std::pair<int, int>  retval;
      retval.first = m_mat_p->rows();
      retval.second= m_mat_p->cols();
      return retval;
    }

    /// deep assign. 
    inline Matrix& operator=(const Matrix& a_input)
    {
      int nrow = a_input.m_mat_p->rows();
      int ncol = a_input.m_mat_p->cols();
      m_mat_p = shared_ptr<eigen_mat>(new eigen_mat(nrow, ncol));
      
      (*m_mat_p) = (*a_input.m_mat_p);
      return *this;
    }
    
    ///increment
    inline Matrix& operator+=(const Matrix& a_matrix)
    {
      (*m_mat_p) += (*a_matrix.m_mat_p);
      return *this;
    }

    ///decrement
    inline Matrix& operator-=(const Matrix& a_matrix)
    {
      (*m_mat_p) -= (*a_matrix.m_mat_p);
      return *this;
    }

    ///scale
    inline Matrix& operator*=(const double& a_scalingFactor)
    {
      (*m_mat_p) *= a_scalingFactor;
      return *this;
    }

    ///get 1/condition number.  Goes to zero if poorly conditioned.
    inline double invConditionNumber() const
    {
      double bigeig, liteig;
      getLargestAndSmallestEigenvalues(bigeig, liteig);

      return liteig/bigeig;
    }

    ///
    inline void
    getLargestAndSmallestEigenvalues(double & a_biggestEigenvalue,
                                     double & a_smallestEigenvalue) const
    {
      Eigen::JacobiSVD<Eigen::MatrixXd> svd(*m_mat_p);
      int ifirs = 0;
      int ilast = svd.singularValues().size()-1;

      double alpha = svd.singularValues()(ifirs); 
      double omega = svd.singularValues()(ilast);
      a_biggestEigenvalue  = alpha;
      a_smallestEigenvalue = omega;
    }

    ///Eigen probably has way faster ways to do this 
    void writeToFile(std::string a_filename) const
    {
      ofstream outfile;
      outfile    << setprecision(15)
                 << setiosflags(ios::showpoint)
                 << setiosflags(ios::scientific);
      
      outfile.open(a_filename.c_str(), ios::out);
      int nrow = m_mat_p->rows();
      int ncol = m_mat_p->cols();
      outfile << "matrix size = " << nrow << " x " << ncol << endl;
      for(int irow = 0; irow < nrow; irow++)
      {
        for(int icol = 0; icol < ncol; icol++)
        {
          outfile << (*this)(irow, icol) << " " ;
        }
        outfile << endl;
      }
      outfile.close();
    }
    
    /// debugging aid
    void poutMaxMins() const
    {
      Chombo4::pout() << setprecision(15)
                      << setiosflags(ios::showpoint)
                      << setiosflags(ios::scientific);
      
      Chombo4::pout() << "Matrix: maxCoef = " << m_mat_p->maxCoeff()
                      <<       ", minCoef = " << m_mat_p->minCoeff()
                      << endl;
    }
    
    ///debugging aid
    void poutAll() const
    {
      Chombo4::pout() << setprecision(3)
                      << setiosflags(ios::showpoint)
                      << setiosflags(ios::scientific);
      int nrow = m_mat_p->rows();
      int ncol = m_mat_p->cols();
      Chombo4::pout() << "matrix size = " << nrow << " x " << ncol << endl;
      for(int irow = 0; irow < nrow; irow++)
      {
        for(int icol = 0; icol < ncol; icol++)
        {
          double val = (*this)(irow, icol);
          if(val >= 0)
          {
            Chombo4::pout() << " ";
          }
          Chombo4::pout() << val << " " ;
        }
        Chombo4::pout() << endl;
      }
    }

    ///debugging aid
    void poutRow(int a_row) const
    {
      Chombo4::pout() << setprecision(3)
                      << setiosflags(ios::showpoint)
                      << setiosflags(ios::scientific);
      int nrow = m_mat_p->rows();
      int ncol = m_mat_p->cols();
      Chombo4::pout() << "matrix size = " << nrow << " x " << ncol << endl;
      int irow = a_row;
      for(int icol = 0; icol < ncol; icol++)
      {
        double val = (*this)(irow, icol);
        if(val >= 0)
        {
          Chombo4::pout() << " ";
        }
        Chombo4::pout() << val << " " ;
      }
      Chombo4::pout() << endl;

    }
    
  
    ///debugging aid
    void poutDiag() const
    {
      Chombo4::pout() << setprecision(15)
                      << setiosflags(ios::showpoint)
                      << setiosflags(ios::scientific);
      int nrow = m_mat_p->rows();
      int ncol = m_mat_p->cols();
      if(nrow != ncol)
      {
        ch_mayday::Error("ch_eigen::poutDiag written for square matricies");
      }
      
      for(int irow = 0; irow < nrow; irow++)
      {
        Chombo4::pout() << "M("<< irow << ","  << irow << ") = "
                        << (*this)(irow, irow) << endl;
      }
    }

 
    ///
    /**
       This function inverts this matix in place.  
       I only want this to work when the matrix is square
       as I prefer to deal with Moore-Penrose manually.
    */
    inline void invert()
    {
      if(m_mat_p->rows() != m_mat_p->cols())
      {
        ch_mayday::Error("ch_eigen::Matrix::invert called for non-square matrix.");
      }
      (*m_mat_p) = m_mat_p->inverse();
    }

    //for tiny cells, set moments to their limits (1 0 0 0...)
    inline void setSmallCellRow(const int& a_irow)
    {
      Matrix& Mvol = *this;
      Mvol(a_irow, 0) = 1.0;
      for(int icol = 1; icol < m_mat_p->cols(); icol++)
      {
        Mvol(a_irow, icol) = 0.0;
      }
    }


    ///
    inline void transpose()
    {
      (*m_mat_p) = (m_mat_p->transpose());
    }

  };//end class ch_eigen::matrix

///
/**
   This function sets a_prod = a_left* a_right.
   This fails if left matrix column count is not the same as the right matrix row count
   (i.e. if a_left.m_mat_p->cols() != a_right.m_mat_p->rows()).
   Eigen may do something smarter but I find a lot of bugs this way.
*/
  void multiply(Matrix      & a_prod, 
                const Matrix& a_left,
                const Matrix& a_righ)
  {
    int leftrow = a_left.m_mat_p->rows();
    int leftcol = a_left.m_mat_p->cols();
    int righrow = a_righ.m_mat_p->rows();
    int righcol = a_righ.m_mat_p->cols();
    if(leftcol != righrow)
    {
      Chombo4::MayDay::Error("ch_eigen::multiply needs left n_columns to match right n_rows");
    }
    a_prod.m_mat_p = shared_ptr<eigen_mat>(new eigen_mat(leftrow, righcol));
    
    (*a_prod.m_mat_p) = (*a_left.m_mat_p) * (*a_righ.m_mat_p);
  }



} //end ch_eigen namespace
#endif
