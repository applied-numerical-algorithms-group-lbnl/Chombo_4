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

#include "Chombo_REAL.H"
//#include "Proto_LoHiSide.H"
#include "Eigen/Core"
#include <utility>
#include <string>
namespace ch_eigen
{

///
/**
   Very simple wrapper around Eigen::matrix class designed to 
   be a transition from using the LAPACKMatrix class toward
   simply using the excellent Eigen infrastrucutre for our
   linear algebra needs.
*/
  class Matrix
  {
  public:
    typedef Chombo4::MayDay  ch_mayday;
    typedef Eigen::MatrixXd  eigen_mat;
    
    ///This is the actual matrix.  The class has no other data
    shared_ptr<eigen_mat> m_mat_p;
    

    Matrix(int a_nrow, int a_ncol)
    {
      m_mat_p = shared_ptr<eigen_mat>(new eigen_mat(a_nrow, a_ncol));
    }


    ///  deep copy constructor
    Matrix(const Matrix& a_input)
    {
      (*m_mat_p) =  (*a_input.m_mat_p)
    }



    inline Real normLTwo() const
    {
      return 0;
    }
  
    inline Real normMax() const
    {
      return 0;
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


    const Real& operator() (int irow, int icol) const
    {
      return (*m_mat_p)(irow, icol);
    }
  
    ///
    double& operator() (int irow, int icol)
    {
      return (*m_mat_p)(irow, icol);
    }
      

    ///return pair <nrow, ncol>
    std::pair<int, int>  dims() const
    {
    
      std::pair<int, int>  retval;
      retval.first = m_mat_p->rows();
      retval.second= m_mat_p->cols();
      return retval;
    }

    /// deep assign. 
    Matrix& operator=(const Matrix& a_matrix)
    {
      (*m_mat_p) = (*a_matrix.m_mat_p);
      return *this;
    }
    
    ///
    Matrix& operator+=(const Matrix& a_matrix)
    {
      (*m_mat_p) += (*a_matrix.m_mat_p);
      return *this;
    }

    ///
    Matrix& operator-=(const Matrix& a_matrix)
    {
      (*m_mat_p) -= (*a_matrix.m_mat_p);
      return *this;
    }

    ///
    Matrix& operator*=(const double& a_scalingFactor)
    {
      (m_mat_p) *= a_scalingFactor;
      return *this;
    }

    double invConditionNumber() const
    {
      double bigeig, liteig;
      getBiggestAndSmallestEigenvalues(bigeig, liteig);
      double inverseConditionNumber = liteig/bigeig;
      return inverseCondtionNumber;
    }

    void
    getBiggestAndSmallestEigenvalues(double & a_biggestEigenvalue,
                                     double & a_smallestEigenvalue)
    {
      JacobiSVD<MatrixXd> svd(*m_mat_p);
      int ifirs = 0;
      int ilast = svd.singularValues().size()-1;
      
      double alpha = svd.singularValues()(ifirst); 
      double omega = svd.singularValues()(ilast ); 
      double cond = svd.singularValues()();
    }

    void writeToFile(std::string a_filename) const;
    ///
    void poutMaxMins() const;
  
    ///
    void poutAll() const;
    
  
    ///
    void poutDiag() const;

 
    ///inverts this matix.  I only want this to work when the matrix is square.
    int invert()
    {
      if(m_mat_p->rows() != m_mat_p->cols())
      {
        ch_mayday::Error("ch_eigen::Matrix::invert called for non-square matrix.");
      }
      (*m_mat_p) = m_mat_p->inverse();
    }

    //for tiny cells, set moments to their limits (1 0 0 0...)
    void setSmallCellRow(const int& a_irow)
    {
      Matrix& Mvol = *this;
      Mvol(a_irow, 0) = 1.0;
      for(int icol = 1; icol < m_mat_p->cols(); icol++)
      {
        Mvol(a_irow, icol) = 0.0;
      }
    }


    ///
    void transpose()
    {
      m_mat_p->transpose();
    }

  };//end class ch_eigen::matrix

///
/**
   sets product = a_left* a_right
   fails if a_left.m_col != a_right.m_rows
*/
  void multiply(Matrix& a_product, 
                const Matrix& a_left,
                const Matrix& a_righ)
  {
    int leftcolumns = a_left.m_mat_p->cols();
    int rightrows   = a_righ.m_mat_p->rows();
    if(leftcolumns != rightrows)
    {
      ch_mayday::Error("ch_eigen::multiply needs left n_columns to match right n_rows")
    }
    (*a_product.m_mat_p) = (*a_left.m_mat_p) * (*a_righ.m_mat_p);
  }



} //end ch_eigen namespace
#endif
