#pragma once

#ifndef _PROTO_EB_DOMAIN_BC_
#define _PROTO_EB_DOMAIN_BC_

#include "EBProto.H"
#include "Proto_EBExactSolutions.H"

namespace Proto
{
  // A simple register that stores things on each domain face
  template<typename T, unsigned int Comps = 1>
  class EBDomainBCRegister
  {
  public:

    inline T
    getRegComp(int a_dir,
               Side::LoHiSide a_side,
               int a_comp=0) const
      {
        CH_assert((a_dir >= 0) && (a_dir < DIM));
        CH_assert(((int)a_side == 0) || (a_side == 1));
        CH_assert((a_comp >= 0) && (a_comp < Comps));
        return m_regbcData[a_dir][(int)a_side][a_comp];
      }

    inline T
    getEBComp(int a_comp=0) const
      {
        CH_assert((a_comp >= 0) && (a_comp < Comps));
        return m_ebbcData[a_comp];
      }

    inline void
    setRegComp(const T& a_value,
               int a_dir,
               Side::LoHiSide a_side,
               int a_comp=0)
      {
        CH_assert((a_comp >= 0) && (a_comp < Comps));
        CH_assert((a_dir >= 0) && (a_dir < DIM));
        CH_assert(((int)a_side == 0) || ((int)a_side == 1));
        m_regbcData[a_dir][(int)a_side][a_comp] = a_value;
      }

    inline void
    setEBComp(const T& a_value,
              int a_comp=0)
      {
        CH_assert((a_comp >= 0) && (a_comp < Comps));
        m_ebbcData[a_comp] = a_value;
      }

    inline void
    setAllRegComp(const T& a_value,
                  int a_comp)
      {
        for(int dir = 0; dir!=DIM; dir++)
          {
            for(SideIterator sit; sit.ok(); ++sit)
              {
                setRegComp(a_value, dir, sit(), a_comp);
              }
          }
      }

    inline void
    setAllComp(const T& a_value,
               int a_comp)
      {
        setAllRegComp(a_value, a_comp);
        setEBComp(a_value, a_comp);
      }

    inline void
    setAll(const T& a_value)
      {
        for(int comp=0; comp!=Comps; comp++)
          {
            setAllComp(a_value, comp);
          }
      }

    // Sometimes we need to work with a single component object.
    // Hopefully typename T is cheap to copy
    inline EBDomainBCRegister<T, 1>
    makeSingleComponent(int a_comp)
      {
        EBDomainBCRegister<T, 1> newSingleComp;
        // copy all the regular values
        for(int dir = 0; dir!=DIM; dir++)
          {
            for(SideIterator sit; sit.ok(); ++sit)
              {
                newSingleComp.setRegComp(this-getRegComp(dir, sit(), a_comp),
                                         dir, sit(), 0);
              }
          }
        // copy the EB
        newSingleComp.setEBComp(this->getEBComp(a_comp));

        return newSingleComp;
      }

    // THIS ONLY EXISTS FOR BACKWARD COMPATIBILITY REASONS
    // DO NOT USE UNLESS NECECCARY
    void
    fillFlatArray(string a_bcnames[2*DIM],
                  int a_comp=0) const
      {
        for(int ivec = 0; ivec != 2*DIM; ivec++)
          {
            int idir;
            Side::LoHiSide iside;
            ebp_face(ivec, idir, iside);
            a_bcnames[ivec] = getRegComp(idir, iside, a_comp);
          }
      }

    void
    fillFlatArray(vector<string> a_bcnames[2*DIM]) const
      {
        for(int ivec = 0; ivec != 2*DIM; ivec++)
          {
            for(int comp = 0; comp!=DIM; comp++)
              {
                a_bcnames[ivec].resize(DIM);
                int idir;
                Side::LoHiSide iside;
                ebp_face(ivec, idir, iside);
                a_bcnames[ivec][comp] = getRegComp(idir, iside, comp);
              }
          }
      }

  private:
    std::array<std::array<std::array<T, Comps>, 2>, DIM> m_regbcData;
    std::array<T, Comps> m_ebbcData;
  };

  // aliases for common use cases
  template<unsigned int Comps> using EBDomainBCtypes = EBDomainBCRegister<std::string, Comps>;
  template<unsigned int Comps, unsigned int order> using EBDomainBCfuncs = EBDomainBCRegister<std::shared_ptr<BaseExactSolution<order> >, Comps>;

  // couple a BC type (Dirichlet, Neuman, ...) with the values on that boundary
  template<unsigned int Comps, unsigned int order>
  class EBDomainBC
  {
  public:

    inline
    EBDomainBCtypes<Comps>
    getBCtypes() const
      {
        return m_bcTypes;
      }

    inline
    EBDomainBCfuncs<Comps, order>
    getBCfuncs() const
      {
        return m_bcFuncs;
      }

    inline
    EBDomainBCtypes<Comps>&
    updateBCtypes()
      {
        return m_bcTypes;
      }

    inline
    EBDomainBCfuncs<Comps, order>&
    updateBCfuncs()
      {
        return m_bcFuncs;
      }

  private:
    EBDomainBCtypes<Comps> m_bcTypes;
    EBDomainBCfuncs<Comps, order> m_bcFuncs;
  };


  // // Utilities to build EBDomainBC from old array storage formats
  // void
  // fillDomBCFromArrays(
  //   EBDomainBCtypes & a_domBC,
  //   int a_comp,
  //   const string & a_ebbcName,
  //   const string & a_dombcName[2*DIM])
  // {
  //   a_domBC.setEBComp(a_ebbcName, a_comp);
  //   for (int face=0; face!=2*DIM; face++)
  //     {
  //       int idir; Side::LoHiSide iside;
  //       ebp_face(face, idir, iside);
  //       a_domBC.setRegComp(a_domainBCName[face][idir], idir, iside, comp);
  //     }
  // }

  // void
  // fillArraysFromDomBC(
  //   cost EBDomainBCtypes & a_domBC,
  //   int a_comp,
  //   string & a_ebbcName,
  //   string & a_dombcName[2*DIM])
  // {
  //   a_domBC.setEBComp(a_ebbcName, a_comp);
  //   for (int face=0; face!=2*DIM; face++)
  //     {
  //       int idir; Side::LoHiSide iside;
  //       ebp_face(face, idir, iside);
  //       a_domBC.setRegComp(a_domainBCName[face][idir], idir, iside, comp);
  //     }
  // }

} // end namespace Proto
#endif // end include guard
