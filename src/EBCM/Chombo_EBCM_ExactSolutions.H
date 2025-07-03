
#include <iostream>
using std::cerr;

#ifndef __EBCM_ExactSolutions__
#define __EBCM_ExactSolutions__
#include "Chombo_EBCM_Basics.H"
// autodiff include
#include <autodiff/forward/real.hpp>
#include <autodiff/forward/real/eigen.hpp>
namespace EBCM_Exact
{
  template <unsigned int order>
  class BaseExactSolution
  {
  public:
    typedef    EBCM::neig_dom_face<                       order > ebcm_dom_face;
    typedef    EBCM::neig_cut_face<                       order > ebcm_cut_face;
    typedef    EBCM::neighborhood<                        order > ebcm_neighborhood;
    typedef    EBCM::EBCM_Volu<                           order > ebcm_volu;
    typedef    EBCM::EBCM_Face<              XFACE      , order > ebcm_xface;
    typedef    EBCM::EBCM_Face<              YFACE      , order > ebcm_yface;
    typedef    EBCM::EBCM_Face<              ZFACE      , order > ebcm_zface;
    typedef    EBCM::IrregFace<                           order > ebcm_irreg_face;
    typedef   Proto::IndexedMoments<         DIM-1,       order > pr_ind_mom_dmo;
    typedef   Proto::IndexedMoments<         DIM  ,       order > pr_ind_mom_dim;
    typedef   Proto::MomentIterator<         DIM  ,       order > pr_momit_dim;
    typedef   Proto::MomentIterator<         DIM-1,       order > pr_momit_dmo;
    typedef   Proto::IndexTM<       int   ,  DIM                > pr_ind_int_dim;
    typedef   Proto::IndexTM<       int   ,  DIM-1              > pr_ind_int_dmo;
    typedef   Proto::IndexTM<       double,  DIM                > pr_ind_dou_dim;
    typedef   Proto::IndexTM<       double,  DIM-1              > pr_ind_dou_dmo;
    typedef   Proto::RealVect                                     pr_rv;
    typedef   Proto::Point                                        pr_pt;
    typedef Chombo4::MayDay                                       ch_mayday;
    virtual ~BaseExactSolution()
    {
    }
  
    BaseExactSolution()
    {
    }

    ///get full face analytical in grad phi dot n dA.
    template <CENTERING  cent>
    double get_analytical_gphi_dot_n(const EBCM::EBCM_Face<cent, order>  & a_ebcm_face,
                                     double                                a_dx) const
    {
      int direction = a_ebcm_face.getFaceDir();
      double retval = 4586.;
      if(direction == 0)
      {
        retval = getAnalytical_IntGPhiX_DA(a_ebcm_face.m_hi.m_pt, a_dx);
      }
      else if(direction == 1)
      {
        retval = getAnalytical_IntGPhiY_DA(a_ebcm_face.m_hi.m_pt, a_dx);
      }
      else if(direction == 2)
      {
        retval = getAnalytical_IntGPhiZ_DA(a_ebcm_face.m_hi.m_pt, a_dx);
      }
      else
      {
        ch_mayday::Error("get_analytical_gphi_dot_n:bogus face direction");
      }
      double area = a_ebcm_face.m_facmom[pr_ind_int_dmo::Zero];
      retval /= area;
      return retval;
    }
    ///get full face analytical in grad phi dot n dA.
    template <CENTERING  cent>
    double get_exact_gphi_dot_n(const EBCM::EBCM_Face<cent, order>  & a_ebcm_face,
                                double                                a_dx) const
    {
      double retval = getIntGradPhiDotNDA(a_ebcm_face, a_dx);

      double area = a_ebcm_face.m_facmom[pr_ind_int_dmo::Zero];
      retval /= area;
      return retval;
    }

    ///
    virtual 
    double getAnalytical_IntGPhiX_DA(const pr_pt & a_ebcm_face_hi,
                                     double        a_dx) const
    {
      static bool s_printed_once = false;
      if(!s_printed_once)
      {
        Chombo4::pout() << "found default analytical flux X" << endl;
        s_printed_once = true;
      }
      return 4586;
    }

    ///
    virtual 
    double getAnalytical_IntGPhiY_DA(const   pr_pt & a_ebcm_face_hi,
                                     double          a_dx) const
    {
      static bool s_printed_once = false;
      if(!s_printed_once)
      {
        Chombo4::pout() << "found default analytical flux Y" << endl;
        s_printed_once = true;
      }
      return 4586;
    }

    ///
    virtual 
    double getAnalytical_IntGPhiZ_DA(const  pr_pt  & a_ebcm_face_hi,
                                     double          a_dx) const
    {
      static bool s_printed_once = false;
      if(!s_printed_once)
      {
        Chombo4::pout() << "found default analytical flux Z" << endl;
        s_printed_once = true;
      }
      return 4586;
    }
    ///
    virtual double
    getAnalyticalLaplacian(const ebcm_volu&  a_volu)
    {
      static bool printed_once = false;
      if(!printed_once)
      {
        Chombo4::pout() << "no analytical laplacian available" << endl;
        printed_once = true;
      }
      return 4586;
    }
    ///
    virtual double
    getAnalyticalPhi(const ebcm_volu&  a_volu)
    {
      static bool printed_once = false;
      if(!printed_once)
      {
        Chombo4::pout() << "no analytical phi available" << endl;
        printed_once = true;
      }
      return 4586;
    }
    /// These functions are necessary for convergence tests
    
    ///
    /**
       Returns taylor series of laplacian_phi integrated over the volume
    **/
    double getCellAveLaplacian(const ebcm_volu&  a_volu) const
    {
      double retval = 4586;
      const pr_ind_mom_dim & volmom = a_volu.m_volmom;
      retval = 0;
      ///this gets int(laplacian) dV
      pr_rv centloc;
      auto pt = a_volu.m_pt;
      auto dx = a_volu.m_dx;
      for(int idir = 0; idir < DIM; idir++)
      {
        centloc[idir] = dx*(double(pt[idir]) + 0.5);
      }
      for(pr_momit_dim momit; momit.ok(); ++momit)
      {
        pr_ind_int_dim mom_power = momit();
        double rfactorial  = (double)(mom_power.factorial());
        for(int idir = 0; idir < DIM; idir++)
        {
          if(mom_power[idir] > 1)
          {
            pr_ind_int_dim lap_power = mom_power;
            lap_power[idir] -= 2;
            double deriv       = getDerivative(mom_power, centloc);
            double cp          = deriv/rfactorial;
            double factor      = cp*(mom_power[idir])*(mom_power[idir]-1);

            double lap_moment  = volmom[lap_power];
            double increment   = factor*lap_moment;
            
            retval += increment;
          } //end if this power is high enough
        } //end loop over directions
      } //end loop over moments
      //need to divide out volume
      double momzero    = volmom[pr_ind_int_dim::Zero];
      retval /= momzero; //Is cell merging grand or what?
      return retval;
    } //end function getCellAveLaplacian
    
    ///
    /**
       Returns |phi|
    **/
    double getCellAveSolution(const ebcm_volu&  a_volu) const
    {
      double retval = 4586;
      const pr_ind_mom_dim & volmom = a_volu.m_volmom;
      ///this gets int(phi) dV
      retval = getIntPhiDV(a_volu);

      ///need to div
      double momzero    = volmom[pr_ind_int_dim::Zero];
      retval /= momzero; //Is cell merging grand or what?
      return retval;
    }

    ///
    /**
       Returns int( phi dV) 
    **/
    double getIntPhiDV(const ebcm_volu&  a_volu) const
    {
      double retval = 0; //needs to be zero
      const pr_ind_mom_dim & volmom = a_volu.m_volmom;
      auto centloc                 = a_volu.m_centroid;
      for(pr_momit_dim momit; momit.ok(); ++momit)
      {
        pr_ind_int_dim power = momit();
        int   idenom = (power.factorial());
        double denom = (Real)(idenom);
        
        double deriv = getDerivative(power, centloc);
        double momen = volmom[power];
        double incre = deriv*momen/denom;
        retval += incre;
      }
      return retval;
      
    }

    /// These functions are necessary for inhomogeneous bcs

    ///
    double getIntPhiDACommonBit(const pr_ind_mom_dmo & a_fac_mom,
                                const pr_ind_dou_dim & a_cent_loc,
                                const int            & a_face_dir,
                                bool a_killing_me_smalls) const
    {
      double retval = 0;  //incremental
      for(pr_momit_dmo momit_dmo; momit_dmo.ok(); ++momit_dmo)
      {
        /// this is p
        pr_ind_int_dmo  taylor_power_dmo = momit_dmo();
        pr_ind_int_dim  taylor_power_dim;
        int dmo_index = 0;
        for(int idir = 0; idir < DIM; idir++)
        {
          if(idir == a_face_dir)
          {
            taylor_power_dim[idir] = 0;
          }
          else
          {
            taylor_power_dim[idir] = taylor_power_dmo[dmo_index];
            dmo_index++;
          }
        }
        int   idenom = (taylor_power_dim.factorial());
        double denom = (Real)(idenom);
        // this is phi^(p)
        double deriv = getDerivative(taylor_power_dim, a_cent_loc);
        double momen =     a_fac_mom[taylor_power_dmo];
        ///increment = phi^(p) m^p/p!
        double incre = deriv*momen/denom;
        retval += incre;
      }
      if(a_killing_me_smalls)
      {
        //Chombo4::pout() << "getIntPhiDACB: retval = " << retval << endl;
      }
      return retval;
      
    }

    
    ///  
    double getIntPhiDA(const ebcm_dom_face     & a_dom_face,
                       const ebcm_neighborhood & a_neigh) const
    {

      const pr_ind_mom_dmo&   fac_mom = a_dom_face.m_facmom;
      const pr_ind_dou_dim&  cent_loc = a_dom_face.m_location;
      double retval = getIntPhiDACommonBit(fac_mom, cent_loc, a_dom_face.m_direction, false);
      return retval;
    }
    ///  
    double getIntPhiDADomainFace(const ebcm_dom_face     & a_dom_face,
                                 const ebcm_neighborhood & a_neigh,
                                 bool a_killing_me_smalls) const
    {

      const pr_ind_mom_dmo&   fac_mom = a_dom_face.m_facmom;
      const pr_ind_dou_dim&  cent_loc = a_dom_face.m_location;
      double retval = getIntPhiDACommonBit(fac_mom, cent_loc, a_dom_face.m_direction, a_killing_me_smalls);
      return retval;
    }

    /// 
    pr_ind_int_dmo  getDMOPower(pr_ind_int_dim taylor_power_dim,int a_face_dir) const
    {
      pr_ind_int_dmo taylor_power_dmo;
      int i_dmo = 0;
      for(int i_dim =0; i_dim < DIM; i_dim++)
      {
        if(i_dim != a_face_dir)
        {
          taylor_power_dmo[i_dmo] = taylor_power_dim[i_dim];
          i_dmo++;
        }
      }
      return taylor_power_dmo;
    }
    
    ///
    double getIntGradPhiDotNDACommonBit(const pr_ind_mom_dmo  & a_fac_mom,
                                        const pr_ind_dou_dim & a_cent_loc,
                                        const int            & a_face_dir) const
    {
      double retval = 0;  //incremental
      for(pr_momit_dim momit_dim; momit_dim.ok(); ++momit_dim)
      {
        ///taylor_power =  p
        pr_ind_int_dim  taylor_power_dim = momit_dim();
        
        int pd       =           taylor_power_dim[a_face_dir];
        if( pd == 1 ) // > 1 has x-x0 = 0
        {
          int   idenom = (taylor_power_dim.factorial());
          double denom = (Real)(idenom);
          ///deriv_power = p-e^d
          pr_ind_int_dim moment_power_dim = taylor_power_dim;
          moment_power_dim[a_face_dir] -= 1;
          pr_ind_int_dmo moment_power_dmo = getDMOPower(moment_power_dim, a_face_dir);
          
          double momen = a_fac_mom[moment_power_dmo];
          double deriv = getDerivative(taylor_power_dim, a_cent_loc);
          ///increment = pd  phi^(p-ed) m^p/p!
          double incre = pd*deriv*momen/denom;
          retval += incre;

        } // end if pd > 0
      }// end loop over moments
      return retval;
    } //end function getGradPhiDotNDACommonBit
    ///
    template <CENTERING  cent>
    double getIntGradPhiDotNDA(const EBCM::EBCM_Face<cent, order>  & a_ebcm_face,
                               double                                a_dx) const
    {
      pr_ind_mom_dmo   fac_mom  = a_ebcm_face.m_facmom;
      pr_ind_dou_dim  cent_loc = a_ebcm_face.location(a_dx);
      int grad_dir = a_ebcm_face.getFaceDir();
      double retval = getIntGradPhiDotNDACommonBit(fac_mom, cent_loc, grad_dir);
      return retval;
    }

    ///
    double getIntGradPhiDotNDA(const ebcm_dom_face     & a_dom_face,
                               const ebcm_neighborhood & a_neigh) const
    {
      const pr_ind_mom_dmo&   fac_mom = a_dom_face.m_facmom;
      const pr_ind_dou_dim&  cent_loc = a_dom_face.m_location;

      double retval = getIntGradPhiDotNDACommonBit(fac_mom, cent_loc, a_dom_face.m_direction);
      return retval;

    }
    ///
    double getIntGradPhiDotNDADomainFace(const ebcm_dom_face     & a_dom_face,
                                         const ebcm_neighborhood & a_neigh) const
    {
      const pr_ind_mom_dmo&   fac_mom = a_dom_face.m_facmom;
      const pr_ind_dou_dim&  cent_loc = a_dom_face.m_location;

      double retval = getIntGradPhiDotNDACommonBit(fac_mom, cent_loc, a_dom_face.m_direction);
      return retval;

    }
    ///
    template <CENTERING  cent>
    double getIntPhiDA(const EBCM::EBCM_Face<cent, order>  & a_ebcm_face,
                       const ebcm_neighborhood             & a_neigh,
                       double                                a_dx) const
    {
      pr_ind_mom_dmo   fac_mom  = a_ebcm_face.m_facmom;
      pr_ind_dou_dim  cent_loc  = a_ebcm_face.location(a_dx);
      int  face_dir = a_ebcm_face.getFaceDir();
      double retval = getIntPhiDACommonBit(fac_mom, cent_loc, face_dir);
      return retval;
    }
    ///
    double getIrregIntPhiDACommonBit(const pr_ind_mom_dim & a_ebf_mom,
                                     const pr_ind_dou_dim & a_cent_loc) const
    {
      double retval = 0;  //incremental
      for(pr_momit_dim momit_dim; momit_dim.ok(); ++momit_dim)
      {
        ///taylor_power =  p
        pr_ind_int_dim  taylor_power_dim = momit_dim();
        
        double momen = a_ebf_mom[taylor_power_dim];
        pr_ind_int_dim deriv_power_dim = taylor_power_dim;
        double deriv = getDerivative(deriv_power_dim, a_cent_loc);
        int   idenom = (taylor_power_dim.factorial());
        double denom = (Real)(idenom);
          
        ///increment =  phi^(p) m^p/p!
        double incre =  deriv*momen/denom;
        retval += incre;
      }// end loop over moments
      return retval;
    } //end function getGradPhiDotNDACommonBit
    ///
    double getIntPhiDA(const EBCM::IrregFace<order>  & a_irreg_face,
                       const ebcm_neighborhood       & a_neigh,
                       double                          a_dx) const
    {
      pr_ind_mom_dim ebf_mom  = a_irreg_face.m_ebfmom;
      pr_ind_dou_dim cent_loc = a_irreg_face.location(a_dx);
      double retval = getIrregIntPhiDACommonBit(ebf_mom, cent_loc);
      return retval;
    }

    ///
    double getIntPhiDA(const ebcm_cut_face     & a_cut_face,
                       const ebcm_neighborhood & a_neigh) const
    {
      pr_ind_mom_dim    ebf_mom   =    a_cut_face.m_ebfmom;
      pr_ind_dou_dim    cent_loc  =    a_cut_face.m_location;
      double retval = getIrregIntPhiDACommonBit(ebf_mom, cent_loc);
      return retval;
    }
    ///
    double getIntPhiDACutFace(const ebcm_cut_face     & a_cut_face,
                              const ebcm_neighborhood & a_neigh) const
    {
      pr_ind_mom_dim    ebf_mom   =    a_cut_face.m_ebfmom;
      pr_ind_dou_dim    cent_loc  =    a_cut_face.m_location;
      double retval = getIrregIntPhiDACommonBit(ebf_mom, cent_loc);
      return retval;
    }
    
    ///
    double getIrregIntGradPhiDotNDACommonBit(const pr_ind_mom_dim & a_ebf_norm_mom_x,
                                             const pr_ind_mom_dim & a_ebf_norm_mom_y,
                                             const pr_ind_mom_dim & a_ebf_norm_mom_z,
                                             const pr_ind_dou_dim & a_cent_loc) const
    {
      double retval = 0;  //incremental
      for(pr_momit_dim momit_dim; momit_dim.ok(); ++momit_dim)
      {
        ///taylor_power =  p
        pr_ind_int_dim  taylor_power_dim = momit_dim();
        for(int grad_dir = 0; grad_dir < DIM; grad_dir++)
        {
          int pd = taylor_power_dim[grad_dir];
          if(pd > 0)
          {
            
            ///remember
            pr_ind_mom_dim ebf_norm_mom;
            if(grad_dir == 0)
            {
              ebf_norm_mom = a_ebf_norm_mom_y;
            }
            else if(grad_dir == 1)
            {
              ebf_norm_mom = a_ebf_norm_mom_y;
            }
            else
            {
              ebf_norm_mom = a_ebf_norm_mom_z;
            }

            // int((x-xbar)^P n_xd dA)
            double norm_mom = ebf_norm_mom[taylor_power_dim];
            //p - e^d
            pr_ind_int_dim deriv_power_dim = taylor_power_dim;
            deriv_power_dim[grad_dir] -= 1;
            double deriv = getDerivative(deriv_power_dim, a_cent_loc);
            int   idenom = (taylor_power_dim.factorial());
            double denom = (Real)(idenom);
          
            ///increment =  phi^(p) m^p/p!
            double incre =  pd*deriv*norm_mom/denom;
            retval += incre;
          } //end if pd>0 (if this power has a non-zero gradient in this direction)
        } //end loop over gradient directions
      } //end loop over moments

      return retval;
    }// end function getIrregIntGradPhiDotNDACommonBit

    ///
    double getIntGradPhiDotNDA(const EBCM::IrregFace<order>  & a_irreg_face,
                               double                          a_dx) const
    {
      pr_ind_dou_dim  cent_loc = a_irreg_face.location(a_dx);
      double retval = getIrregIntGradPhiDotNDACommonBit(a_irreg_face.m_norm_mom_x,
                                                        a_irreg_face.m_norm_mom_y,
                                                        a_irreg_face.m_norm_mom_z,
                                                        cent_loc);
      return retval;
    }
    
    
    ///
    double getIntGradPhiDotNDA(const ebcm_cut_face     & a_cut_face,
                               const ebcm_neighborhood & a_neigh,
                               double                    a_dx) const
    {
      double retval = getIrregIntGradPhiDotNDACommonBit(a_cut_face.m_norm_mom_x,
                                                        a_cut_face.m_norm_mom_y,
                                                        a_cut_face.m_norm_mom_z,
                                                        a_cut_face.m_location);

      return retval;
    }

    ///
    virtual double
    getIntGradPhiDotNDACutFace(const ebcm_cut_face     & a_cut_face,
                               const ebcm_neighborhood & a_neigh,
                               double                    a_dx) const
    {
      double retval = getIrregIntGradPhiDotNDACommonBit(a_cut_face.m_norm_mom_x,
                                                        a_cut_face.m_norm_mom_y,
                                                        a_cut_face.m_norm_mom_z,
                                                        a_cut_face.m_location);

      return retval;
    }


    ///Where all the real work gets done (in the derived class)
    virtual Real getDerivative(pr_ind_int_dim  a_power,
                               pr_rv           a_loc) const
    {
      Chombo4::MayDay::Error("EBCM::BaseExact:: not implemented");
    }


  }; //end class baseExactSolution  
  ///class for always returning zero (for bcs as opposed to being for testing convergence rates)
  template <unsigned int order>
  class ZeroSolution: public BaseExactSolution<order>
  {
  public:
  
    typedef    EBCM::EBCM_Volu<                           order > ebcm_volu;
    typedef       Proto::IndexTM<int, DIM>            pr_ind_int_dim;
    typedef       Proto::RealVect                     pr_rv;
    virtual double
    getAnalyticalLaplacian(const ebcm_volu&  a_volu)
    {
      return 0;
    }
    virtual Real getDerivative(pr_ind_int_dim  a_power,
                               pr_rv           a_loc) const
    {
      return 0.;
    }

    virtual ~ZeroSolution()
    {
    }
  
    ZeroSolution()
    {
    }

  };




  ///phi = cp x^p  
  template <unsigned int order>
  class MonomialEF: public BaseExactSolution<order>
  {
  public:
  
    typedef       IndexTM<int, DIM>            IvDIM;
    typedef IndexedMoments<DIM  , order>       IndMomDIM;
    typedef HostIrregData<CELL,IndMomDIM , 1>  VoluData;
  
    virtual Real getDerivative(IndexTM<int, DIM>    a_deriv,
                               RealVect a_loc) const
    {
      Real  C     = m_entry.second;  //coefficient
      Point expon = m_entry.first; //exponent
    
      //first deal with the zeros.  If the deriv in a particular direction is higher
      //than the original exponent, the derivative is zero
      for(int idir = 0; idir < DIM; idir++)
      {
        PROTO_ASSERT((a_deriv[idir] >= 0), "invalid derivative");
        PROTO_ASSERT((  expon[idir] >= 0), "invalid exponenent");
        if(a_deriv[idir] > expon[idir])
        {
          return 0;
        }
      }

      //Past this point, we can assume by construction that the exponent is bigger than the derivative
      // (exp >= deriv),  in all directions.   We can also assume that both >= 0.
      //int px = expon[0];
      //int py = expon[1];
      //using std::pow;
#if DIM==3
      //int pz = expon[2];
#endif
      //int sum = a_deriv.sum();
      //all the return vals are constants multiplied by this stuff
      Real ppowbit = 1;
      for(int idir = 0; idir < DIM; idir++)
      {
        ppowbit *= (pow(a_loc[idir], expon[idir]-a_deriv[idir]));
      }
      //this is to denote the cascading derivative factors ...(expon-3)*(expon-2)(expon)
      Real retval_simp = C*ppowbit;
      for(int idir = 0; idir < DIM; idir++)
      {
        int curDeriv = a_deriv[idir];
        int curExpon =   expon[idir];
        Real factor = 1;
        while(curDeriv > 0)
        {
          factor *= Real(curExpon);
          curDeriv--;
          curExpon--;
        }
        retval_simp *= factor;
      }

      return retval_simp;
    }
    virtual ~MonomialEF()
    {
    }

    MonomialEF(std::pair<Point, Real>  a_entry)
    {
      m_entry = a_entry;
    }

  private:
    std::pair<Point, Real>  m_entry;

  
    MonomialEF()
    {
    }
  };
  ///phi = amp* sin( pi*freq x)
  template <unsigned int order>
  class OneDimSineEF: public BaseExactSolution<order>
  {
  public:
    typedef    EBCM::EBCM_Volu<                           order > ebcm_volu;
    int    m_direction;
    double m_amplitude;
    int    m_frequency;
    OneDimSineEF(int    a_direction,
                 double a_amplitude,
                 int    a_frequency):
      m_direction(a_direction),
      m_amplitude(a_amplitude),
      m_frequency(a_frequency)
    {
    }

    virtual Real getDerivative(IndexTM<int, DIM>  a_deriv,
                               RealVect           a_loc) const
    {
      //any  derivatives in cross directions  = 0
      for(int idir = 0; idir < DIM; idir++)
      {
        if((idir != m_direction) && (a_deriv[idir] != 0))
        {
          return 0;
        }
      }
      double x_loc     =   a_loc[m_direction];

      int    deriv_mag = a_deriv[m_direction];
      bool   sin_now = true;
      bool   cos_now = false;
      static const double s_pi = 4.*atan(1.0);
      double trig_factor =   1;
      for(int i_der = 1; i_der <= deriv_mag; i_der++)
      {
        if(sin_now)
        {
          trig_factor *= s_pi*m_frequency;
        }
        if(cos_now)
        {
          trig_factor *= -s_pi*m_frequency;
        }
        //swap stuff for next round
        sin_now = !sin_now;
        cos_now = !cos_now;
      }
      double base_trig = 4586;
      if(sin_now)
         base_trig = sin(s_pi*x_loc*m_frequency);
      else if(cos_now)
         base_trig = cos(s_pi*x_loc*m_frequency);
      double retval = m_amplitude*trig_factor*base_trig;

      return retval;
    }

    virtual double
    getAnalyticalLaplacian(const ebcm_volu&  a_volu)
    {
      static const double pi = 4.*atan(1.0);
      double dx = a_volu.m_dx;
      double xhi = dx*(double(a_volu.m_pt[0]) +   1);
      double xlo = dx*(double(a_volu.m_pt[0])      );
      double fhi = pi*cos(pi*xhi);
      double flo = pi*cos(pi*xlo);
        
      double retval = (fhi-flo)/dx;

      return retval;
    }
    virtual ~OneDimSineEF()
    {
    }
    private:
    OneDimSineEF();
  };

  ///phi = amp* sin( pi*freq_x x)) * sin(pi*freq_y y) * sin(pi * freq_z * z)
  template <unsigned int order>
  class SineProductEF: public BaseExactSolution<order>
  {
  public:
    typedef    Proto::RealVect         pr_rv;
    typedef    Proto::Point            pr_pt;
    typedef    EBCM::EBCM_Volu<order > ebcm_volu;

    double m_amplitude;
    vector<int> m_frequency;
    
    ///
    virtual void
    getFullCellBounds(pr_rv& a_lo, pr_rv& a_hi, const ebcm_volu&  a_volu)
    {
      for(int idir = 0; idir < DIM; idir++)
      {
        double cent_ind = a_volu.m_pt[idir];
        a_lo[idir] = (a_volu.m_dx)*(cent_ind);
        a_hi[idir] = (a_volu.m_dx)*(cent_ind + 1);
      }
    }

    /// return full cell value of (1/V)(int(phi) dV)
    virtual double
    getAnalyticalPhi(const ebcm_volu&  a_volu)
    {
      static const double  s_pi = 4.*atan(1.);
      static bool printed_once = false;
      if(!printed_once)
      {
        Chombo4::pout() << "found sineproduct analytical phi" << endl;
        printed_once = true;
      }
      pr_rv lo_bound, hi_bound;
      getFullCellBounds(lo_bound, hi_bound, a_volu);
      double integral      = m_amplitude; // result is a product
      double full_cell_vol = 1;
      
      for( int idir = 0; idir < DIM; idir++)
      {
        full_cell_vol *= a_volu.m_dx;
        double freq_pi = s_pi*double(m_frequency[idir]);
        integral *= (-cos(freq_pi*hi_bound[idir]) 
                     +cos(freq_pi*lo_bound[idir]) )/freq_pi;
      }

      double average = integral/full_cell_vol;
      return average;
    }
    
    /// return full cell value of (1/V)(int(lap phi) dV)
    virtual double
    getAnalyticalLaplacian(const ebcm_volu&  a_volu)
    {
      static bool s_printed_once = false;
      if(!s_printed_once)
      {
        Chombo4::pout() << "found sineproduct analytical laplacian" << endl;
        s_printed_once = true;
      }
      double phi_val = getAnalyticalPhi(a_volu);
      static const double  s_pi = 4.*atan(1.);
      double lap_val = 0; //additive
      for(int idir = 0; idir < DIM; idir++)
      {
        double freq_pi = s_pi*double(m_frequency[idir]);
        lap_val -= freq_pi*freq_pi*phi_val;
      }
      lap_val *= m_amplitude;
      return lap_val;
    }
    
    ///
    virtual 
    double getAnalytical_IntGPhiX_DA(const pr_pt  & a_ebcm_face_hi,
                                     double         a_dx) const
    {
      static bool s_printed_once = false;
      if(!s_printed_once)
      {
        Chombo4::pout() << "found sineproduct analytical flux X" << endl;
      }
      double retval = 4586.;
#if DIM==3
      if(!s_printed_once)
      {
        Chombo4::pout() << "3d version of ineproduct analytical flux X goes here" << endl;
        s_printed_once = true;
      }

#else
      
      double x_c = a_dx*(double(a_ebcm_face_hi[0]) + 0.5);
      double y_l = a_dx*(double(a_ebcm_face_hi[1])      );
      double y_h = a_dx*(double(a_ebcm_face_hi[1]) + 1.0);
      static const double  s_pi = 4.*atan(1.);
      pr_rv freq_pi;
      for(int idir = 0; idir < DIM; idir++)
      {
        freq_pi[idir] = s_pi*double(m_frequency[idir]);
      }
      if(m_frequency[1] == 0)
      {
        retval = -a_dx*freq_pi[0]*cos(x_c*freq_pi[0]);
      }
      else
      {
        double factor = -freq_pi[0]/freq_pi[1];
        retval = factor*cos(x_c*freq_pi[0])*(cos(y_h*freq_pi[1]) - cos(y_l*freq_pi[1]));
      }
#endif
      retval *= m_amplitude;
      s_printed_once = true;
      return retval;
    }

    ///
    virtual 
    double getAnalytical_IntGPhiY_DA(const pr_pt & a_ebcm_face_hi,
                                     double        a_dx) const
    {
      static bool s_printed_once = false;
      if(!s_printed_once)
      {
        Chombo4::pout() << "found sineproduct analytical flux y" << endl;
      }
      double retval = 4586.;
#if DIM==3
      if(!s_printed_once)
      {
        Chombo4::pout() << "3d version of ineproduct analytical flux Y goes here" << endl;
        s_printed_once = true;
      }
#else
      
      double y_c = a_dx*(double(a_ebcm_face_hi[1]) + 0.5);
      double x_l = a_dx*(double(a_ebcm_face_hi[0])      );
      double x_h = a_dx*(double(a_ebcm_face_hi[0]) + 1.0);
      static const double  s_pi = 4.*atan(1.);
      pr_rv freq_pi;
      for(int idir = 0; idir < DIM; idir++)
      {
        freq_pi[idir] = s_pi*double(m_frequency[idir]);
      }
      if(m_frequency[0] == 0)
      {
        retval = -a_dx*freq_pi[1]*cos(y_c*freq_pi[1]);
      }
      else
      {
        double factor = -freq_pi[1]/freq_pi[0];
        retval = factor*cos(y_c*freq_pi[1])*(cos(x_h*freq_pi[0]) - cos(x_l*freq_pi[0]));
      }
#endif
      retval *= m_amplitude;
      s_printed_once = true;
      return retval;
    }

    ///
    virtual 
    double getAnalytical_IntGPhiZ_DA(const pr_pt & a_ebcm_face_hi,
                                     double        a_dx) const
    {
      static bool s_printed_once = false;
      if(!s_printed_once)
      {
        Chombo4::pout() << "found the function for sineproduct analytical flux Z (but it is not written yet)" << endl;
        s_printed_once = true;
      }
      return 4586;
    }
    ///
    SineProductEF(double a_amplitude,
                  vector<int> a_frequency)
    {
      m_amplitude =  a_amplitude;
      m_frequency =  a_frequency;
    }

    virtual Real getDerivative(IndexTM<int, DIM>  a_deriv,
                               RealVect           a_loc) const
    {
      //compute the derivatives of each direction separately
      // and multiply ?
      double directional_factors[DIM];

      for(int deriv_dir = 0; deriv_dir < DIM; deriv_dir++)
      {
        //set amplitude to 1 to avoid double-counting.
        OneDimSineEF<order> one_dim_sine(deriv_dir, 1.0, m_frequency[deriv_dir]);
        IndexTM<int, DIM>   one_dim_deriv;
        one_dim_deriv.setToBasisVector(deriv_dir);
        one_dim_deriv *= a_deriv[deriv_dir];
        directional_factors[deriv_dir] = one_dim_sine.getDerivative(one_dim_deriv, a_loc);
      }

      double retval = m_amplitude;
      for(int idir = 0; idir < DIM; idir++)
      {
        retval *= directional_factors[idir];
      }
      return retval;
    }
    virtual ~SineProductEF()
    {
    }
    private:
    SineProductEF();
  };
  ///phi = sum cp x^p  
  template <unsigned int order>
  class PolynomialEF: public BaseExactSolution<order>
  {
  public:
  
    typedef       IndexTM<int, DIM>            IvDIM;
    typedef IndexedMoments<DIM  , order>       IndMomDIM;
    typedef HostIrregData<CELL,IndMomDIM , 1>  VoluData;
  
    virtual Real getDerivative(IndexTM<int, DIM>  a_deriv,
                               RealVect a_loc) const
    {
      Real retval = 0;
      for(int ient = 0; ient < m_entries.size(); ient++)
      {
        MonomialEF<order> mono(m_entries[ient]);
        Real monoval = mono.getDerivative(a_deriv, a_loc);
        retval += monoval;
      }
      return retval;
    }

    virtual ~PolynomialEF()
    {
    }

    PolynomialEF(std::vector<std::pair<Point, Real> > a_entries)
    {
      m_entries = a_entries;
    }

  private:
    std::vector<std::pair<Point, Real> > m_entries;

  
    PolynomialEF()
    {
    }
  };
}

#endif
