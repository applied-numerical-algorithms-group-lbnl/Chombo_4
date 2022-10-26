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
        string  a_ebbc = string("no_slip_wall"))
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
         string  a_ebbc = string("no_slip_wall"))
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

  string& operator()(unsigned int idir, Side::LoHiSide a_sd)
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
        MayDay::Error((("EBIBC: unknown bc string " + a_bcnames[ivec]).c_str()));
      }
    }
  }


  void scalarDiffusionStencilStrings(string a_bcnames[2*DIM]) const
  {
    for(int ivec = 0; ivec < 2*DIM; ivec++)
    {
      a_bcnames[ivec] = StencilNames::Neumann;
    }
  }
  

  void projectionStencilStrings(string a_bcnames[2*DIM],
                                string& a_ebbcname) const
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
        MayDay::Error((("EBIBC: unknown bc string " + a_bcnames[ivec]).c_str()));
      }
    }
    a_ebbcname = StencilNames::Neumann;
    //a_ebbcname = StencilNames::Dirichlet;
  }

  void insStencilStrings(vector<string> a_bcnames[2*DIM],
                         string& a_ebbcname,
                         bool a_allwallslip=false) const
    {
      for(int ivec = 0; ivec < 2*DIM; ivec++)
        {
          for(int comp = 0; comp!=DIM; comp++)
            {
              a_bcnames[ivec].resize(DIM);
              int idir;
              Side::LoHiSide iside;
              ebp_face(ivec, idir, iside);
              const string& faceBC = m_domainBC[ivec];
              if(faceBC == "no_slip_wall" && !a_allwallslip)
                {
                  a_bcnames[ivec][comp] = StencilNames::Dirichlet;
                }
              else if(faceBC == "slip_wall" ||
                      (faceBC == "no_slip_wall" && a_allwallslip))
                {
                  if(idir == comp) // normal constrained
                    {
                      a_bcnames[ivec][comp] = StencilNames::Dirichlet;
                    }
                  else // tangent is free
                    {
                      a_bcnames[ivec][comp] = StencilNames::NoBC;
                    }
                }
              else if(faceBC == "inflow")
                {
                  a_bcnames[ivec][comp] = StencilNames::Dirichlet;
                }
              else if(faceBC == "outflow")
                {
                  if (a_allwallslip)
                    {
                      a_bcnames[ivec][comp] = StencilNames::NoBC;
                    }
                  else
                    {
                      a_bcnames[ivec][comp] = StencilNames::Neumann;
                    }
                }
              else
                {
                  MayDay::Error((("EBIBC: unknown bc string " + a_bcnames[ivec][comp]).c_str()));
                }
              //a_bcnames[ivec][comp] = StencilNames::NoBC;
            }
        }

      if(m_ebbc == "no_slip_wall" && !a_allwallslip)
        {
          a_ebbcname = StencilNames::Dirichlet;
        }
      else if(m_ebbc == "slip_wall" ||
              (m_ebbc == "no_slip_wall" && a_allwallslip))
        {
          a_ebbcname = StencilNames::NoBC;
        }
      else
        {
          MayDay::Error("EBIBC: unknown bc string");
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
