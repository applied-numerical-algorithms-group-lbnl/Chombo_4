#include <cmath>
#include <cstdio>
#include <iostream>

#include "EBProto.H"
#include "Chombo_EBLevelBoxData.H"
#include "Chombo_LevelData.H"
#include "Chombo_BaseFab.H"

#include "Chombo_ParmParse.H"
#include "Chombo_LoadBalance.H"
#include "Chombo_ProtoInterface.H"
#include "Chombo_BRMeshRefine.H"
#include "Chombo_GeometryService.H"
#include "Chombo_EBDictionary.H"
#include "Chombo_EBChombo.H"
#include "EBMultigrid.H"
#include "Proto_DebugHooks.H"
#include "DebugFunctions.H"
#include "Hoeb_ExactSolutions.H"
//this one defines HOEB_MAX_ORDER
#include "Hoeb_Utilities.H"
#include "Hoeb_LAPACKMatrix.H"
#include <iomanip>


using std::endl;
using hoeb::LAPACKMatrix;
using Chombo4::pout;
/**/
struct dotprod_t
{
  Real          m_prodval;
  unsigned int  m_rownum0;
  unsigned int  m_rownum1;
  Real          m_weight0;
  Real          m_weight1;
  bool operator<(const dotprod_t& a_input) const
  {
    return m_prodval < a_input.m_prodval;
  }
};

std::ostream& operator<< (std::ostream&       a_os,
                          const dotprod_t&    a_iv)
{
  static bool printedhead = false;
  if(!printedhead)
  {
    a_os << "prodval    \t col0\t col1\t weight0\t weight1\t" << endl;
    printedhead = true;
  }
  a_os      << setprecision(4)
            << setiosflags(ios::showpoint)
            << setiosflags(ios::scientific);
  a_os << a_iv.m_prodval <<  "\t"
       << a_iv.m_rownum0 <<  "\t"
       << a_iv.m_rownum1 <<  "\t"
       << a_iv.m_weight0 <<  "\t"
       << a_iv.m_weight1 <<  "\t" ;
  return a_os;

}
  
std::set<dotprod_t, std::less<dotprod_t> > 
getRowDotProducts(const LAPACKMatrix& a_Amat, const LAPACKMatrix& a_weight)
{
  unsigned int nrow = a_Amat.dims().first;
  std::set<dotprod_t, std::less<dotprod_t> > allprod;
  
  for(unsigned int irow = 0; irow < nrow; irow++)
  {
    LAPACKMatrix iRowVec = a_Amat.rowVector(irow);
    for(unsigned int jrow = 0; jrow < nrow; jrow++)
    {
      if(irow < jrow)
      {
        LAPACKMatrix jRowVec = a_Amat.rowVector(jrow);
        Real dotprodval = iRowVec.dotProduct(jRowVec);
        dotprod_t setval;
        setval.m_prodval = dotprodval;
        setval.m_rownum0 = irow;
        setval.m_rownum1 = jrow;
        setval.m_weight0 = a_weight(irow, irow);
        setval.m_weight1 = a_weight(jrow, jrow);
        allprod.insert(setval);
      }
    }
  }
  
  return allprod;
}
///
void
removeLinearlyDependentEquations(LAPACKMatrix                              & a_Areduced,
                                 LAPACKMatrix                              & a_Wreduced,
                                 const LAPACKMatrix                        & a_origA,
                                 const LAPACKMatrix                        & a_weight,
                                 Real a_eps)
{
//figure out which which rows are most linearly independent
  using Chombo4::pout;
  using std::endl;
  
  std::set<dotprod_t, std::less<dotprod_t> > allprod = getRowDotProducts(a_origA, a_weight);
  vector<dotprod_t> vecprod;
  for(auto it = allprod.begin(); it != allprod.end(); ++it)
  {
    vecprod.push_back(*it);
  }

  pout() << "all dot products: " << endl;
  std::sort(vecprod.begin(), vecprod.end());
//  for(unsigned int ivec = 0; ivec < vecprod.size(); ivec++)
//  {
//    pout() << vecprod[ivec] << endl;
//  }

  
  //now let's try taking out the smallest
  //min number rows == number of columns
  unsigned int numColumns    = a_origA.dims().second;
  unsigned int totalRows     = a_origA.dims().first;
  unsigned int minNumberRows = numColumns;
  int maxRemove = totalRows - minNumberRows -1; //-1 so weights still have meaning
  
  unsigned int vecsize = vecprod.size();
  //because sorted
  Real maxval = vecprod[vecsize-1].m_prodval;
  Real minval = vecprod[0].m_prodval;
  //dot products should always be positive
  PR_assert(maxval >= 0);
  PR_assert(minval >= 0);
  Real magval= maxval - minval;
  Real toosmall = a_eps*magval;


  std::set<unsigned int> rowsToRemove;
  for(unsigned int ivec = 0; ivec < vecprod.size(); ivec++)
  {
    auto& thisprod = vecprod[ivec];
    if((thisprod.m_prodval < toosmall) && (rowsToRemove.size() < maxRemove))
    {
      if(thisprod.m_weight0 < thisprod.m_weight1)
      {
        rowsToRemove.insert(thisprod.m_rownum0);
      }
      else
      {
        rowsToRemove.insert(thisprod.m_rownum1);
      }
    }      
  }
  unsigned int newRowTotal = totalRows - rowsToRemove.size();
  a_Areduced.define(newRowTotal, numColumns);
  a_Wreduced.define(newRowTotal, newRowTotal);
  a_Wreduced.setToIdentity();
  unsigned int currentRow = 0;
  for(int irow = 0; irow < totalRows; irow++)
  {
    const bool skipThisOne = rowsToRemove.find(irow) != rowsToRemove.end();
    if(!skipThisOne)
    {
      for(int icol = 0; icol < numColumns; icol++)
      {
        a_Areduced(currentRow, icol) = a_origA(irow, icol);
      }
      a_Wreduced(currentRow, currentRow) = a_weight(currentRow, currentRow);
      currentRow++;
    }
  }
}
///
void
checkWeightedMoorePenroseInverse(const LAPACKMatrix  & a_M,
                                 const LAPACKMatrix  & a_W,
                                 const string        & a_prefix)
{
  using Chombo4::pout;
  using std::endl;
  pout() << a_prefix << "Moment matrix M condition number:" << endl;
  a_M.checkConditionNumber();
  LAPACKMatrix WM;
  multiply(WM,a_W, a_M);
  pout() << a_prefix << "Weighted Moment matrix WM condition number:" << endl;
  WM.checkConditionNumber();
  
  
  LAPACKMatrix WMT = WM;
  WMT.transpose();
  LAPACKMatrix WMTWM;
  multiply(WMTWM, WMT, WM);
  pout() << a_prefix << "Squared Weighted Moment matrix WMTWM condition number:" << endl;
  WMTWM.checkConditionNumber();

  LAPACKMatrix WMTWMInv = WMTWM;
  WMTWMInv.invert();
  LAPACKMatrix checkMat;
  multiply(checkMat, WMTWM, WMTWMInv);
  pout() << a_prefix << "checkmat WMTWM WMTWMInv (should be the identity matrix):" << endl;
  checkMat.poutAll();
}
///
unsigned int lapackTest()
{
  using Chombo4::pout;
  using std::endl;
  LAPACKMatrix Amat, weight;
  Amat.defineFromFile(  string("Mmat.0.matrix"));
  weight.defineFromFile(string("weight0.matrix"));
  
  LAPACKMatrix Areduced, Wreduced;
  removeLinearlyDependentEquations(Areduced, Wreduced, Amat, weight, 1.e-4);

  checkWeightedMoorePenroseInverse(Amat, weight, string("full system "));

  checkWeightedMoorePenroseInverse(Areduced, Wreduced, string("reduced system "));

//check dot products 
  return 0;
}
/****/

int main(int a_argc, char* a_argv[])
{
  int retval = lapackTest();
  return retval;
}
