#ifndef _EBIBC_H_
#define _EBIBC_H_
#include <cmath>
#include <memory>
#include <string>
#include "Proto.H"
#include "EBProto.H"
#include "Chombo_EBChombo.H"

#include "Chombo_NamespaceHeader.H"
using std::string;
///
class EBIBC
{
public:

  /// 
  EBIBC(string  a_veloIC,
        string  a_scalIC,
        string  a_loSideDomainBC[DIM],
        string  a_hiSideDomainBC[DIM],
        string  a_ebbc = string("NoSlipWall"))
  {
    for(int idir = 0;  idir < DIM; idir++)
    {
      unsigned int loIndex = ebp_index(idir, Side::Lo);
      unsigned int hiIndex = ebp_index(idir, Side::Hi);
      m_domainBC[loIndex] = a_loSideDomainBC[idir];
      m_domainBC[hiIndex] = a_loSideDomainBC[idir];
    }

    m_ebbc   =   a_ebbc  ;
    m_veloIC =   a_veloIC;
    m_scalIC =   a_scalIC;
  }

  string getDomainIndex(unsigned int idir, Side::LoHiSide a_sd)
  {
    unsigned int index = ebp_index(idir, a_sd);
    return m_domainBC[index];
  }
  
  string         m_domainBC[2*DIM];
  string         m_ebbc;
  string         m_veloIC;
  string         m_scalIC;

private:

  EBIBC();

};
#include "Chombo_NamespaceFooter.H"

#endif
