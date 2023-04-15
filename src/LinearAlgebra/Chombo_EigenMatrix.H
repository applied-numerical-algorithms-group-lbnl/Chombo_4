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
   Very simple wrapper around Eigen::matrix class.
*/
  class Matrix
  {
  public:
    typedef Chombo4::MayDay  ch_mayday;
    typedef Eigen::MatrixXd  eigen_mat;
    
    ///this is the actual matrix
    shared_ptr<eigen_mat> m_mat_p;
    

    Matrix(int a_nrow, int a_ncol)
    {
      m_mat_p = shared_ptr<eigen_mat>(new eigen_mat(a_nrow, a_ncol));
    }


    ///  deep copy constructor
    Matrix(const Matrix& a_input)
    {
      ch_mayday::Error("not implemented");
    }



    inline Real normLTwo() const
    {
      return 0;
    }
  
    inline Real normMax() const
    {
      return 0;
    }

    inline void setToIdentity()
    {
    }
  

    ///
    inline void  setVal(const double& a_val);


    const Real& operator() (int irow, int icol) const;

      
  
    ///
    Real & operator() (int irow, int icol);

    //return row, col
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
      //
      return *this;
    }


    ///
    void poutAll() const;
  
    void writeToFile(std::string a_filename) const;
    ///
    void poutMaxMins() const;
  
    ///
    void poutDiag() const;

    ///
    void poutDiagMatlab() const;

    ///
    void poutMatlab() const;
 
    ///inverts this matix 
    /**
       fails if matrix is not square
       if return value != 0, probably a singular matrix
    */
    int invert();

    //for tiny cells, set moments to their limits (1 0 0 0...)
    void setSmallCellRow(const int& irow)
    {
      Matrix& Mvol = *this;
      Mvol(irow, 0) = 1.0;
      int ncol = m_mat_p->cols();
      for(int icol = 1; icol < ncol; icol++)
      {
        Mvol(irow, icol) = 0.0;
      }
    }


    ///
    void transpose();



  };//end class eigenmatrix

///
/**
   sets product = a_left* a_right
   fails if a_left.m_col != a_right.m_rows
*/
  void multiply(Matrix& a_product, 
                const Matrix& a_left,
                const Matrix& a_right);



} //end ch_eigen namespace
#endif
