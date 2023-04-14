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
#include <utility>
#include <string>
namespace ch_eigen
{

///
/**
   Getting sick of writing the basics here over and over.
   Silly class but it will cut down on the typing.
*/
  class EigenMatrix
  {
    typedef Chombo4::MayDay  ch_mayday; 
  public:
    bool m_defined;

    EigenMatrix()
    {
      m_defined = false;
    }
  
    EigenMatrix(int a_nrow, int a_ncol)
    {
      define(a_irow, a_icol);
    }

    ///
    inline void define(int nrow, int ncol)
    {
      m_defined = true;
      ch_mayday::Error("not implemented");
    }


    ///  deep copy constructor
    EigenMatrix(const EigenMatrix& a_input)
    {
      ch_mayday::Error("not implemented");
    }



    inline Real normLTwo() const
    {
    }
  
    inline Real normMax() const
    {
    }

    inline void setToIdentity()
    {
    }
  

    ///
    ~EigenMatrix()
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
      retval.first = m_nrow;
      retval.second= m_ncol;
      return retval;
    }

    /// deep assign. 
    EigenMatrix& operator=(const EigenMatrix& a_matrix)
    {
      return *this;
    }
    
    ///
    EigenMatrix& operator+=(const EigenMatrix& a_matrix);
    {
      return *this;
    }

    ///
    EigenMatrix& operator-=(const EigenMatrix& a_matrix);
    {
      return *this;
    }

    ///
    EigenMatrix& operator*=(const double& a_scalingFactor);
    {
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
      EigenMatrix& Mvol = *this;
      Mvol(irow, 0) = 1.0;
      for(int icol = 1; icol < m_ncol; icol++)
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
  void multiply(EigenMatrix& a_product, 
                const EigenMatrix& a_left,
                const EigenMatrix& a_right);



} //end ch_eigen namespace
#endif
