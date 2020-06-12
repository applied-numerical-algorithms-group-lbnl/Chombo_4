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

  EBIBC()
  {
    m_isDefined = false;
  }
  /// 
  EBIBC(string  a_veloIC,
        string  a_scalIC,
        string  a_loSideDomainBC[DIM],
        string  a_hiSideDomainBC[DIM],
        string  a_ebbc = string("NoSlipWall"))
  {
    define( a_veloIC,                     
            a_scalIC,                     
            a_loSideDomainBC,
            a_hiSideDomainBC,
            a_ebbc);
  }
  ///
  void
  define(string  a_veloIC,
         string  a_scalIC,
         string  a_loSideDomainBC[DIM],
         string  a_hiSideDomainBC[DIM],
         string  a_ebbc = string("NoSlipWall"))
  {
    m_isDefined = true;
    for(int idir = 0;  idir < DIM; idir++)
    {
      unsigned int loIndex = ebp_index(idir, Side::Lo);
      unsigned int hiIndex = ebp_index(idir, Side::Hi);
      m_domainBC[loIndex] = a_loSideDomainBC[idir];
      m_domainBC[hiIndex] = a_hiSideDomainBC[idir];
    }

    m_ebbc   =   a_ebbc  ;
    m_veloIC =   a_veloIC;
    m_scalIC =   a_scalIC;
  }

  string getDomainIndex(unsigned int idir, Side::LoHiSide a_sd) const
  {
    unsigned int index = ebp_index(idir, a_sd);
    return m_domainBC[index];
  }

  void helmholtzStencilStrings(string a_bcnames[2*DIM]) const
  {
    for(int ivec = 0; ivec < 2*DIM; ivec++)
    {
      const string& faceBC = m_domainBC[ivec];
      if(faceBC == "no_slip_wall")
      {
        a_bcnames[ivec] = StencilNames::Dirichlet;
      }
      else if(faceBC == "slip_wall")
      {
        a_bcnames[ivec] = StencilNames::Neumann;
      }
      else if(faceBC == "inflow")
      {
        a_bcnames[ivec] = StencilNames::Neumann;
      }
      else if(faceBC == "outflow")
      {
        a_bcnames[ivec] = StencilNames::Neumann;
      }
      else
      {
        MayDay::Error("EBIBC: unknown bc string");
      }
    }
  }


  void projectionStencilStrings(string a_bcnames[2*DIM]) const
  {
    for(int ivec = 0; ivec < 2*DIM; ivec++)
    {
      const string& faceBC = m_domainBC[ivec];
      if(faceBC == "no_slip_wall")
      {
        a_bcnames[ivec] = StencilNames::Neumann;
      }
      else if(faceBC == "slip_wall")
      {
        a_bcnames[ivec] = StencilNames::Neumann;
      }
      else if(faceBC == "inflow")
      {
        a_bcnames[ivec] = StencilNames::Neumann;
      }
      else if(faceBC == "outflow")
      {
        a_bcnames[ivec] = StencilNames::Dirichlet;
      }
      else
      {
        MayDay::Error("EBIBC: unknown bc string");
      }
    }
  }
  string         m_domainBC[2*DIM];
  string         m_ebbc;
  string         m_veloIC;
  string         m_scalIC;

private:
  bool m_isDefined;


};
#include "Chombo_NamespaceFooter.H"

#endif
