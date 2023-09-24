#ifndef __Chombo_EBCM_Elliptic_Framework_
#define __Chombo_EBCM_Elliptic_Framework_

///Nothing in this file is usable without PETSc.

#ifdef CH_USE_PETSC 
#include "petsc.h"
#include "petscmat.h"
#include "petscvec.h"
#include "petscksp.h"
#include "petscviewer.h"
#include "petscmat.h"

///  Chombo_EBCM_Graph.H is where MetaDataLevel lives.  It made sense at the time.
#include "Chombo_EBCM_Graph.H"
#include "Chombo_EBCM_PETSc_Framework.H"
#include "Chombo_EBCM_ExactSolutions.H"
#include "Chombo_EBLevelBoxData.H"

/// 
/**
   This is for elliptic PDEs.   This contains HelmholtzOp a constant coefficient Helmholtz operator that conforms to the 
   PETScFramework::BasePETScOp.    
**/
namespace EBCM
{

  
  template <int ebcm_order>
  class Elliptic_Framework
  {
  public:

    typedef ch_eigen::Matrix                                                eigen_mat;
    
    typedef Chombo4::Box                                                    ch_box;
    typedef Chombo4::DisjointBoxLayout                                      ch_dbl;
    typedef Chombo4::BoxIterator                                            ch_bit;
    typedef Chombo4::IntVect                                                ch_iv;
    typedef Chombo4::DataIterator                                           ch_dit;
    typedef Chombo4::DataIndex                                              ch_dat_ind;
    typedef Chombo4::ParmParse                                              ch_parmparse;
                                                                            
    typedef Proto::Point                                                    pr_pt;
    typedef Proto::RealVect                                                 pr_realvect;
    typedef Proto::IndexTM<double, DIM>                                     pr_itm_r_dim;
    typedef Proto::IndexTM<int ,   DIM>                                     pr_itm_i_dim;
    typedef Proto::IndexTM<double, DIM-1>                                   pr_itm_r_dmo;
    typedef Proto::IndexTM<int ,   DIM-1>                                   pr_itm_i_dmo;
    typedef Proto::IndexedMoments<DIM  ,      ebcm_order >                  pr_mom_dim;
    typedef Proto::IndexedMoments<DIM-1,      ebcm_order >                  pr_mom_dmo;
                                                                            
    typedef EBCM::Algorithm_Framework<        ebcm_order >                  algo_framework;
    typedef EBCM::PETSc_Framework<            ebcm_order >                  petsc_framework;
    
    typedef EBCM::neighborhood   <            ebcm_order >                  ebcm_neighborhood;
    typedef EBCM::MetaDataLevel<              ebcm_order >                  ebcm_meta;
    typedef EBCM::EBCM_Graph<                 ebcm_order >                  ebcm_graph;
    typedef EBCM::HostLevelData<int,    1,    ebcm_order >                  ebcm_int_sca_data;
    typedef EBCM::HostLevelData<double, 1,    ebcm_order >                  ebcm_dou_sca_data;
    typedef EBCM::HostLevelData<int,    DIM,  ebcm_order >                  ebcm_int_vec_data;
    typedef EBCM::HostLevelData<double, DIM,  ebcm_order >                  ebcm_dou_vec_data;
    typedef EBCM::SubVolumeVector<            ebcm_order >                  ebcm_subvol_vec;
    typedef EBCM::EBCM_Volu<                  ebcm_order >                  ebcm_volu;
    typedef typename EBCM::PETSc_Framework<   ebcm_order >::local_stencil_t ebcm_local_stencil;
    typedef typename EBCM::PETSc_Framework<   ebcm_order >::Base_PETSc_Op   ebcm_base_op;
    typedef typename EBCM::Analytic_Framework<ebcm_order >                  ebcm_analytic;
    typedef typename ebcm_analytic::BaseExactSolution                       ebcm_base_exact;

    

    ///
    template <CENTERING cent>
    static shared_ptr<eigen_mat>
    getHelmholtzQMatrix( const EBCM_Face<cent, ebcm_order>  & a_face,
                         shared_ptr<ebcm_neighborhood>        a_locality,
                         const int                          & a_facedir,
                         const double                       & a_dx,
                         bool a_print)
    {
      pr_mom_dmo facmom =  a_face.m_facmom;
      pr_mom_dim volmom = ebcm_basics::getVolMomFromFacMom(facmom, a_facedir);
      double momzero    = facmom[pr_itm_i_dmo::Zero];

      pr_itm_r_dim xbar = a_locality->getXBar(a_face, a_dx);
      volmom.shift(xbar);

      unsigned int Pcols = pr_mom_dim::size();
      shared_ptr<eigen_mat> retval(new eigen_mat(1, Pcols));


      for(MomentIterator<DIM, ebcm_order> momit; momit.ok(); ++momit)
      {
        auto momind = momit();
        double matval = 0;
        if(momind[a_facedir] > 0)
        {
          int pd = momind[a_facedir];
          pr_itm_i_dim derind = momind - pr_itm_i_dim::BASISV(a_facedir);
          Real derval = volmom[derind];
          if(std::abs(momzero) > 1.0e-10)
          {
            matval = pd*derval/momzero;
          }
        }
        unsigned int currentCol = pr_mom_dim::indexOf(momind);
        (*retval)(0,currentCol) = matval;
      }
    
      return retval;
    }

    ///
    template<CENTERING cent>
    static ebcm_local_stencil
    getHelmholtzFDotEhatStencil(const EBCM_Face<cent, ebcm_order>  & a_face,
                                const int                          & a_stenrad,
                                const int                          & a_weightpower,
                                const EBCM::EBCM_Graph<ebcm_order> & a_graph,
                                const string                       & a_dom_bc,
                                const string                       & a_eb_bc,
                                const int                          & a_facedir,
                                const double                       & a_dx,
                                bool a_print = false)
    {
      /**
         From the document (Lawrence Berkeley National Laboratory Technical Report LBNL-2001537):  
         Q : Q matrix (Unfortunately, I did not give it a fancy mathy name in the document.
         It is the one that converts  Taylor coefficients to fluxes).   
         S : stencil matrix.  
         P : phi vector.  
         W : weighting matrix.
         M : moment matrix.
         A: system matrix =  ((WM)^T(WM))
         (S P) = <\nabla \psi> = Q C= (Q ((WM)^T(WM))^{-1} (WM)^T  P) Q A^{-1} WM^T P.
         S = Q C = (Q ((WM)^T(WM))^{-1} (WM)^T  ) = Q A^{-1} (WM)^T.
      **/
      shared_ptr< EBCM::neighborhood< ebcm_order > > locality;
      shared_ptr< eigen_mat > Amat_p =
        EBCM::Algorithm_Framework<ebcm_order>::
        getSystemMatrix(locality, a_face, a_graph,
                        a_stenrad, a_weightpower,
                        a_dom_bc, a_eb_bc, a_print);
      shared_ptr<eigen_mat> Wmat_p = algo_framework::getWeightMatrix(locality, a_print);
      shared_ptr<eigen_mat> Mmat_p = algo_framework::getMomentMatrix(locality, a_dom_bc, a_eb_bc,a_print);

      eigen_mat WMmat;
      multiply(WMmat, *Wmat_p, *Mmat_p);
      eigen_mat WMTmat = WMmat;
      WMTmat.transpose();

      eigen_mat Ainvmat = *Amat_p;
      Ainvmat.invert();
      
      shared_ptr<eigen_mat> Qmat_p = getHelmholtzQMatrix(a_face, locality, a_facedir, a_dx, a_print);
      eigen_mat AinvWMTmat;
      multiply( AinvWMTmat, Ainvmat, WMTmat);
      eigen_mat Smat;
      multiply( Smat, *Qmat_p, AinvWMTmat);
      ebcm_local_stencil retval = petsc_framework::getStencilFromMat(Smat, locality);
      return retval;
    }

    static double
    getEBQMatrixIncrement(pr_itm_i_dim  a_momind,
                          int           a_direction,
                          double        a_momzerofactor,
                          pr_mom_dim    a_norm_mom)
    {
      double increment = 0;
      if(a_momind[a_direction] > 0)
      {
        int pd = a_momind[a_direction];
        pr_itm_i_dim derind = a_momind - pr_itm_i_dim::BASISV(a_direction);
        double ndda  = a_norm_mom[derind];
        increment = a_momzerofactor*pd*ndda;
      }
      return increment;
    }
    ///
    static shared_ptr<eigen_mat>
    getHelmholtzEBQMatrix( const IrregFace< ebcm_order >      & a_face,
                           shared_ptr<ebcm_neighborhood>        a_locality,
                           const double                       & a_dx,
                           bool a_print)
    {
      
      pr_mom_dim ebf_mom    =  a_face.m_ebfmom;
      pr_mom_dim norm_mom_x =  a_face.m_norm_mom_x;
      pr_mom_dim norm_mom_y =  a_face.m_norm_mom_y;
      pr_mom_dim norm_mom_z =  a_face.m_norm_mom_y;

      double momzero    = ebf_mom[pr_itm_i_dim::Zero];
      double momzerofactor = 0;
      if(std::abs(momzero) > 1.0e-10) 
      {
        momzerofactor = 1./momzero;
      }

      pr_itm_r_dim xbar = a_locality->getXBarEB(a_face, a_dx);
      //only needed for zeroth moment (which is not altered by shifts)
      //ebf_mom.shift(xbar); 
      norm_mom_x.shift(xbar);
      norm_mom_y.shift(xbar);
      norm_mom_z.shift(xbar);
      unsigned int Pcols = pr_mom_dim::size();
      shared_ptr<eigen_mat> retval(new eigen_mat(1, Pcols));

      for(MomentIterator<DIM, ebcm_order> momit; momit.ok(); ++momit)
      {
        pr_itm_i_dim momind = momit();
        double xincr = getEBQMatrixIncrement(momind, 0, momzerofactor, norm_mom_x);
        double yincr = getEBQMatrixIncrement(momind, 1, momzerofactor, norm_mom_y);
#if DIM== 3      
        double zincr = getEBQMatrixIncrement(momind, 2, momzerofactor, norm_mom_z);
#else       
        double zincr = 0;
#endif      
        double matval =  xincr + yincr + zincr;
        unsigned int currentCol = pr_mom_dim::indexOf(momind);
        (*retval)(0,currentCol) = matval;
      }
    
      return retval;

    }
    ///Returns a local stencil for \int (\grad \phi) \cdot (\ehat^facedir)
    /**
       Get stencil for integrated gradient in a  coordinate direction. 
       This is the same as the average gradient * face area.
       \int (\grad \phi) \cdot (\ehat^d) dA = A_f < \frac{\partial \phi}{\partial x_d} >
    */
    static ebcm_local_stencil
    getIrregHelmholtzFDotNStencil(const IrregFace<ebcm_order> & a_face,
                                  const int                          & a_stenrad,
                                  const int                          & a_weightpower,
                                  const EBCM::EBCM_Graph<ebcm_order> & a_graph,
                                  const string                       & a_dom_bc,
                                  const string                       & a_eb_bc,
                                  const double                       & a_dx,
                                  bool a_print = false)
    {

      shared_ptr< EBCM::neighborhood< ebcm_order > > locality;
      shared_ptr< eigen_mat > Amat_p =
        EBCM::Algorithm_Framework<ebcm_order>::
        getEBSystemMatrix(locality, a_face, a_graph,
                          a_stenrad, a_weightpower,
                          a_dom_bc, a_eb_bc, a_print);

      shared_ptr<eigen_mat> Wmat_p = algo_framework::getWeightMatrix(locality, a_print);
      shared_ptr<eigen_mat> Mmat_p = algo_framework::getMomentMatrix(locality, a_dom_bc, a_eb_bc, a_print);

      eigen_mat WMmat;
      multiply(WMmat, *Wmat_p, *Mmat_p);
      eigen_mat WMTmat = WMmat;
      WMTmat.transpose();

      eigen_mat Ainvmat = *Amat_p;
      Ainvmat.invert();
      
      shared_ptr<eigen_mat> Qmat_p = getHelmholtzEBQMatrix(a_face, locality,  a_dx, a_print);
      eigen_mat AinvWMTmat;
      multiply( AinvWMTmat, Ainvmat, WMTmat);
      eigen_mat Smat;
      multiply( Smat, *Qmat_p, AinvWMTmat);
      ebcm_local_stencil retval = petsc_framework::getStencilFromMat(Smat, locality);
      return retval;
    }
    
    ///
    /**
       This class is  a Helmholtz operator.
       For a volume-averaged field phi, 
       this creates a volume-averaged operator of the form
       L(phi) = alpha phi + beta * div(G(phi))
       I am jettisoning the kappa-weighting trick for this operator at least.
       Ergo, this should only be used with cell merging.  
       Without it, we have other operators.
       <Div(Grad(phi))> = (1/vol) (int (Div Grad(phi) ) dV)  = (1./vol)(int (Grad(phi) \cdot \nhat dA)).
       This means that the RHS for a projection will have to be Div(v) not kappa(Div(v)).
    **/
    class HelmholtzOp: public ebcm_base_op
    {
    public:

      typedef Proto::Side              pr_sd;
      typedef Proto::Box               pr_box;
      typedef Proto::IndexTM<int, DIM> pr_int_dim;

      ///member data is in the clear.
      double                                 m_alpha;
      double                                 m_beta;
      int                                    m_weight_power;
      int                                    m_sten_rad;
      string m_dom_bc;
      string m_eb_bc;
      /// 
      /**
         alpha        = identity coeff
         beta         = coeff of laplcian
         ghost        = neighbor ghost size (grow the face box by this much)
         weight_power = goes to neighborhood to weight vofs.
         meta         = ebcm's very data-intensive meta data
      **/
      HelmholtzOp(const double                                 & a_alpha,
                  const double                                 & a_beta,
                  const int                                    & a_sten_rad,
                  const int                                    & a_weight_power,
                  const shared_ptr< ebcm_meta  >               & a_meta_ptr,
                  string a_dom_bc, string a_eb_bc, 
                  bool                                           a_print = false)
      :ebcm_base_op(a_meta_ptr, a_sten_rad,  a_print)
      {
        m_alpha        = a_alpha;       
        m_beta         = a_beta;
        m_sten_rad     = a_sten_rad;
        m_weight_power = a_weight_power;
        m_dom_bc       = a_dom_bc;
        m_eb_bc        = a_eb_bc;
      }

      virtual ~HelmholtzOp()
      {
      }

      ///this is the function that the driver cod calls
      virtual shared_ptr<ebcm_local_stencil> getLocalStencil(const pr_pt& a_pt, const ebcm_graph& a_graph, bool a_print = false)
      {
        shared_ptr<ebcm_local_stencil> retval =
          privateStencil(a_pt,          
                         m_sten_rad,     
                         m_weight_power, 
                         a_graph,       
                         m_dom_bc,      
                         m_eb_bc,       
                         a_print);        
        return retval;
      }

    private:
      ///Weak construction makes everything harder in the end (and not in a good way).
      HelmholtzOp(); 
      
      ///
      /**
         This is where the canal locks are dug in the basement.
         I am jettisoning the kappa-weighting trick for this operator at least.
         Div(Grad(phi)) = (1/vol) (int (div F) dV)  = (1./vol)(int (F dot n dA)).
         This means that the RHS for a projection will have to be div v not kappa div v.
      **/
      virtual shared_ptr<ebcm_local_stencil> privateStencil(const pr_pt                         & a_pt, 
                                                            const int                          & a_stenrad,
                                                            const int                          & a_weightpower,
                                                            const EBCM::EBCM_Graph<ebcm_order> & a_graph,
                                                            const string                       & a_dom_bc,
                                                            const string                       & a_eb_bc,
                                                            bool                                 a_print        = false)
      {
        typedef Proto::Box pr_bx;
        
        shared_ptr<ebcm_local_stencil> retval(new ebcm_local_stencil());
        pr_bx pointbx(a_pt, a_pt);
        pr_bx grbx = pointbx.grow(pr_pt::Ones(this->m_ghost));
        
        auto volume = a_graph.getVolumeCoveringPoint(a_pt);
        auto volmag = volume.m_volmom[pr_int_dim::Zero];
        auto dx = volume.m_dx;
        
        {
          for(int iface = 0; iface < volume.m_xfaces.size(); iface++)
          {
            const auto& face = volume.m_xfaces[iface];
            double fluxsign = faceSign(volume, face);
            ebcm_local_stencil flux = getHelmholtzFDotEhatStencil(face,
                                                                  a_stenrad,     
                                                                  a_weightpower, 
                                                                  a_graph,       
                                                                  a_dom_bc,      
                                                                  a_eb_bc,       
                                                                  dx, 0, a_print);        
            flux *= fluxsign;
            (*retval) += flux;
          }
        }
        {
          for(int iface = 0; iface < volume.m_yfaces.size(); iface++)
          {
            const auto& face = volume.m_yfaces[iface];
            double fluxsign = faceSign(volume, face);
            ebcm_local_stencil flux = getHelmholtzFDotEhatStencil(face,
                                                                  a_stenrad,     
                                                                  a_weightpower, 
                                                                  a_graph,       
                                                                  a_dom_bc,      
                                                                  a_eb_bc,       
                                                                  dx, 1, a_print);        

                                                               
            flux *= fluxsign;
            (*retval) += flux;
          }
        }
        {
          for(int iface = 0; iface < volume.m_zfaces.size(); iface++)
          {
            const auto& face = volume.m_zfaces[iface];
            double fluxsign = faceSign(volume, face);
            ebcm_local_stencil flux = getHelmholtzFDotEhatStencil(face,
                                                                  a_stenrad,     
                                                                  a_weightpower, 
                                                                  a_graph,       
                                                                  a_dom_bc,      
                                                                  a_eb_bc,       
                                                                  dx, 2, a_print);        
                                                               
            flux *= fluxsign;
            (*retval) += flux;
          }
        }
        {
          for(int iface = 0; iface < volume.m_ifaces.size(); iface++)
          {
            const auto& face = volume.m_ifaces[iface];
            ebcm_local_stencil flux = getIrregHelmholtzFDotNStencil(face,
                                                                    a_stenrad,     
                                                                    a_weightpower, 
                                                                    a_graph,       
                                                                    a_dom_bc,      
                                                                    a_eb_bc,
                                                                    dx, a_print);
                                                                 
            (*retval) += flux;
          }
        }

        ///I think this should be dividing by dx^d
        (*retval) *= (1.0/volmag);
        return retval;
      }

    }; //end class HelmholtzOp
    /**
       just parse stuff.  This decouples geometry generation with parmparse, allowing me to make varying refinements.
    **/ 
    static void parseFineGridParams(int            & a_nx,
                                    double         & a_dx,
                                    ch_probdom     & a_domain,
                                    int            & a_maxGrid,
                                    bool           & a_mergeSmallCells,
                                    double         & a_alpha,
                                    double         & a_beta,
                                    int            & a_sten_rad,
                                    int            & a_weight_power,
                                    string         & a_dom_bc,
                                    string         & a_eb_bc, 
                                    bool           a_print = false)
    {
      ParmParse pp("parseFineGridParams");

      pp.get("nx"             , a_nx);
      pp.get("maxGrid"        , a_maxGrid);
      pp.get("mergeSmallCells", a_mergeSmallCells);
      pp.get("alpha"          , a_alpha);
      pp.get("beta"           , a_beta );
      pp.get("stencil_rad"    , a_sten_rad);
      pp.get("weight_power"   , a_weight_power);
      pp.get("dom_bc"         , a_dom_bc);
      pp.get("eb_bc"          , a_eb_bc);

      ch_iv domLo = ch_iv::Zero;
      ch_iv domHi  = (a_nx - 1)*ch_iv::Unit;
      a_domain = ch_probdom(domLo, domHi);
      a_dx = 1./(double(a_nx));
      
    }
    /**
       Make geometry at specified refinment.
       This decouples geometry generation with parmparse, allowing me to make varying refinements.
    **/ 
    static
    shared_ptr< ebcm_meta  >
    makeGeometryAtSpecificNX( const int            & a_nx,
                              const double         & a_dx,
                              const ch_probdom     & a_domain,
                              const int            & a_maxGrid,
                              const bool           & a_mergeSmallCells,
                              const bool           & a_print)
    {
      typedef Chombo4::GeometryService<ebcm_order> ch_geoserv;
      
      vector<ch_dbl> vecgrids;
      Chombo4::pout() << "making grids" << endl;
      ch_geoserv::generateGrids(vecgrids, a_domain.domainBox(), a_maxGrid);

      int geomGhost = 6;
      shared_ptr<BaseIF>    impfunc = algo_framework::getImplicitFunction();
      Chombo4::pout() << "defining geometry in EB land" << endl;

      pr_rv origin = pr_rv::Zero();
      shared_ptr< ch_geoserv > geoserv
        (new ch_geoserv(impfunc, origin, a_dx, a_domain.domainBox(), vecgrids, geomGhost));

      int ghost = 0;
      shared_ptr< ebcm_meta  >
        metaDataPtr(new ebcm_meta(geoserv, a_domain.domainBox(), a_dx, ghost,
                                  a_mergeSmallCells, a_print));
      return metaDataPtr;
    }

    ///
    static void
    setExact(shared_ptr<ebcm_dou_sca_data>     a_phi_exact_ptr,
             shared_ptr<ebcm_dou_sca_data>     a_lph_exact_ptr,
             shared_ptr< ebcm_meta       >     a_meta_ptr,
             shared_ptr<ebcm_base_exact  >     a_exact_ptr,
             bool a_print = false)
    {
      const auto& grids  = a_meta_ptr->m_grids;
      const auto& graphs = *(a_meta_ptr->m_graphs);
      const auto&    dx  = a_meta_ptr->m_dx;
      auto           dit = grids.dataIterator();
      for(int ibox = 0; ibox < dit.size(); ibox++)
      {
        auto      & phifab = (*a_phi_exact_ptr)[dit[ibox]];
        auto      & lphfab = (*a_lph_exact_ptr)[dit[ibox]];
        const auto&  valid =  grids[dit[ibox]];
        const auto&  graph = graphs[dit[ibox]];
        ebcm_subvol_vec volvec(graph, valid, a_print);

        for(int ivec = 0; ivec < volvec.size(); ivec++)
        {
          const auto& volume = volvec[ivec];
          for(int icell = 0; icell < volume.m_cells.size(); icell++)
          {
            const auto& cell = volume.m_cells[icell];
            const pr_mom_dim& volmom = volume.m_volmom;
            pr_mom_dim    moments = volmom;
            pr_realvect   location = volume.m_centroid;

            double phi_value = a_exact_ptr->getPhi(      moments, location);
            double lph_value = a_exact_ptr->getLaplacian(moments, location);
            
            phifab(cell, 0) = phi_value;
            lphfab(cell, 0) = lph_value;
          }
        }
      }
    }
    ///
    /**
       Define exact solution data from input file.  I tried to this
       make this self-documenting.  I came up with this so I could use
       maxima to have nice smooth solutions that meet the boundary
       conditions.  Inhomogeneous boundary conditions are messy.
       Since EBCM_INS (Remember Alice? This is a song about Alice)
       will not need inhomogeneous boundary conditions
       I would like to avoid dealing with that complexity just for
       testing purposes.  
    **/
    static shared_ptr<ebcm_base_exact>
    getExactSolution()
    {
      string whichphi;
      //print out the diagnostics only once
      static bool printedonce = false;
      ch_parmparse  pp("getExactSolution");
      pp.get("which_phi", whichphi);
      shared_ptr<ebcm_base_exact>  retval;
      if(whichphi == string("SineSphere"))
      {
        typedef typename ebcm_analytic::SineSphereEF analytic_sine;

        bool   inside = false;
        double radius = 4586;
        double center = 4586;
        pp.get("radius", radius);
        pp.get("inside", inside);
        pp.get("center", center);

        if(!printedonce)
          Chombo4::pout() << "getExactSolution: exact solution = sinesphere with radius = " << radius << ", center = " << center << endl;
        ebcm_base_exact*  raw_ptr = static_cast<ebcm_base_exact*>(new analytic_sine(radius, center));
        retval =  shared_ptr<ebcm_base_exact>(raw_ptr);
      }
      else if(whichphi == string("Polynomial"))
      {
        typedef typename ebcm_analytic::PolynomialEF analytic_poly;
        vector<pair<pr_pt, double> > entries;
        int num_terms;
        pp.get("polynomial_num_terms", num_terms);
        entries.resize(num_terms);
        if(!printedonce)
          Chombo4::pout() << "getExactSolution: polynomial chosen" << endl;
        for(int iterm = 0; iterm < num_terms; iterm++)
        {
          string coef_str   = string("polynomial_coefficient_") + to_string(iterm);
          string powers_str = string("polynomial_powers_")      + to_string(iterm);
          vector<int> entry_powers;
          Real        entry_coef;
          pp.getarr(powers_str.c_str(), entry_powers, 0, DIM);
          pp.get(     coef_str.c_str(), entry_coef);
          pair<Point, Real> entry;
          entry.second = entry_coef;
          for(int idir = 0; idir < DIM; idir++)
          {
            entry.first[idir] = entry_powers[idir];
          }
          entries.push_back(entry);
          if(!printedonce)
            Chombo4::pout() << "getExactSolution: iterm = " << iterm << ", power = " << entry.first << ", weight = " << entry.second << endl;
        }
        ebcm_base_exact *  raw_ptr = static_cast< ebcm_base_exact * >(new analytic_poly(entries));
        retval = shared_ptr<ebcm_base_exact>(raw_ptr);
      }
      else if(whichphi == string("OneDimSineProduct"))
      {
        typedef Proto::RealVect                            pr_realvect;
        typedef typename  ebcm_analytic::OneDimSineProduct analytic_poly;
        vector<double> vec_freq(DIM, 4586.);
        vector<double> vec_cent(DIM, 4586.);
        double amplitude;
        pp.getarr("sine_product_freq", vec_freq, 0, DIM);
        pp.getarr("sine_product_cent", vec_cent, 0, DIM);
        pp.get(   "sine_product_amp", amplitude);
        pr_realvect  freq;
        pr_realvect  cent;
        for(int idir = 0; idir < DIM; idir++)
        {
          freq[idir] = vec_freq[idir] ;
          cent[idir] = vec_cent[idir] ;
        }
        
        if(!printedonce)
        {
          Chombo4::pout() << "getExactSolution: OneDimSineProduct  chosen."                  << endl;
          Chombo4::pout() << "getExactSolution: OneDimSineProduct frequency = " << freq      << endl;
          Chombo4::pout() << "getExactSolution: OneDimSineProduct center    = " << cent      << endl;
          Chombo4::pout() << "getExactSolution: OneDimSineProduct amplitude = " << amplitude << endl;
        }
        ebcm_base_exact *  raw_ptr = static_cast< ebcm_base_exact * >(new analytic_poly(amplitude, freq, cent));
        retval = shared_ptr<ebcm_base_exact>(raw_ptr);
      }
      else if(whichphi == string("OneDimCosineProduct"))
      {
        typedef Proto::RealVect                              pr_realvect;
        typedef typename  ebcm_analytic::OneDimCosineProduct analytic_poly;
        vector<double> vec_freq(DIM, 4586.);
        vector<double> vec_cent(DIM, 4586.);
        double amplitude;
        pp.getarr("cosine_product_freq", vec_freq, 0, DIM);
        pp.getarr("cosine_product_cent", vec_cent, 0, DIM);
        pp.get(   "cosine_product_amp", amplitude);
        pr_realvect freq;
        pr_realvect cent;
        for(int idir = 0; idir < DIM; idir++)
        {
          freq[idir] = vec_freq[idir] ;
          cent[idir] = vec_cent[idir] ;
        }
        
        if(!printedonce)
        {
          Chombo4::pout() << "getExactSolution: OneDimCosineProduct chosen."                  << endl;
          Chombo4::pout() << "getExactSolution: OneDimCosineProduct frequency = " << freq      << endl;
          Chombo4::pout() << "getExactSolution: OneDimCosineProduct center    = " << cent      << endl;
          Chombo4::pout() << "getExactSolution: OneDimCosineProduct amplitude = " << amplitude << endl;
        }
        ebcm_base_exact *  raw_ptr = static_cast< ebcm_base_exact * >(new analytic_poly(amplitude, freq, cent));
        retval = shared_ptr<ebcm_base_exact>(raw_ptr);
      }
      else
      {
        Chombo4::MayDay::Error("getExactSolution: Doh! (string for whichphi unrecognized)");
      }

      printedonce = true;
      return retval;
    }
    ///
    static double
    getMaxAbsTruncation(double                     a_alpha,
                        double                     a_beta,
                        int                        a_sten_rad,
                        int                        a_weight_power,
                        shared_ptr< ebcm_meta  >   a_meta_ptr,
                        string                     a_dom_bc,
                        string                     a_eb_bc,
                        bool                       a_printToFile,
                        string                     a_file_prefix,
                        bool                       a_print)

    {
      int num_ghost = 4;
      shared_ptr<ebcm_dou_sca_data> phi_exac_ptr(new ebcm_dou_sca_data(a_meta_ptr, num_ghost));
      shared_ptr<ebcm_dou_sca_data> lph_exac_ptr(new ebcm_dou_sca_data(a_meta_ptr, num_ghost));
      shared_ptr<ebcm_dou_sca_data> lph_calc_ptr(new ebcm_dou_sca_data(a_meta_ptr, num_ghost));
      shared_ptr<ebcm_dou_sca_data> err_calc_ptr(new ebcm_dou_sca_data(a_meta_ptr, num_ghost));

      shared_ptr<ebcm_base_exact> exactly_ptr = getExactSolution();
      
      setExact(phi_exac_ptr, lph_exac_ptr, a_meta_ptr, exactly_ptr, a_print);

      HelmholtzOp discreteOp(a_alpha, a_beta, a_sten_rad,
                             a_weight_power, a_meta_ptr,
                             a_dom_bc, a_eb_bc, a_print);

      discreteOp.applyOp(lph_calc_ptr, phi_exac_ptr);
      //settting  error to lph_calc - lph_exac
      double acoef =  1;
      double bcoef = -1;
      err_calc_ptr->setVal(0.);
      err_calc_ptr->setToAxBy(acoef, *lph_calc_ptr, bcoef, *lph_exac_ptr);

      if(a_printToFile)
      {
        string phi_exac_file = a_file_prefix +  string("phi_exac.hdf5");
        string lph_exac_file = a_file_prefix +  string("lph_exac.hdf5");
        string lph_calc_file = a_file_prefix +  string("lph_calc.hdf5");
        string err_calc_file = a_file_prefix +  string("err_calc.hdf5");
      
        phi_exac_ptr->writeToHDF5(phi_exac_file);
        lph_exac_ptr->writeToHDF5(lph_exac_file);
        lph_calc_ptr->writeToHDF5(lph_calc_file);
        err_calc_ptr->writeToHDF5(err_calc_file);
      }

      double retval = err_calc_ptr->maxNorm();
      return retval;
    }

    
    ///
    /**
       One of the tests called by  example/EBCM/hoeb_truncation.
    **/ 
    static  void run_hoeb_truncation_tests()
    {


      static const int num_ghost = 4;
      
      int    nx_fine      = 4586;
      int    maxGrid      = 4586;    
      int    sten_rad     = 4586; 
      int    weight_power = 4586; 
      double dx_fine      = 4586.;
      double alpha        = 4586.;
      double beta         = 4586.;
      string dom_bc(       "4586");
      string  eb_bc(       "4586");
      bool mergeSmallCells;
      ch_probdom domain_fine;
      bool print = false;

      parseFineGridParams(nx_fine, dx_fine, domain_fine,
                          maxGrid, mergeSmallCells,
                          alpha, beta, sten_rad, weight_power,
                          dom_bc, eb_bc, print);
      
      int    nx_coar = nx_fine/2;
      double dx_coar = dx_fine*2.;
      ch_probdom domain_coar = domain_fine;
      domain_coar.coarsen(2);

      shared_ptr< ebcm_meta  > meta_fine = makeGeometryAtSpecificNX(nx_fine, dx_fine, domain_fine, maxGrid, mergeSmallCells, print);
      shared_ptr< ebcm_meta  > meta_coar = makeGeometryAtSpecificNX(nx_coar, dx_coar, domain_coar, maxGrid, mergeSmallCells, print);

      double max_err_fine = getMaxAbsTruncation(alpha, beta, sten_rad, weight_power, meta_fine, dom_bc, eb_bc, true, string("fine_"), print); 
      double max_err_coar = getMaxAbsTruncation(alpha, beta, sten_rad, weight_power, meta_coar, dom_bc, eb_bc, true, string("coar_"), print);

      Chombo4::pout() << "hoeb_truncation: Lmax(coarse error) =  " << max_err_coar << endl;
      Chombo4::pout() << "hoeb_truncation: Lmax(fine   error) =  " << max_err_fine << endl;
      double tol = 1.0e-12;
      double err = std::abs(max_err_fine); 
      if(err > tol)
      {
        double order = (log(max_err_coar/max_err_fine))/log(2.0);
        Chombo4::pout() << "hoeb_truncation: order of accuracy = " <<  order  << endl;
      }
      else
      {
        Chombo4::pout() << "hoeb_truncation: Fine error is too small for safe order computation." << endl;
        
      }

    }
   ///
   /**      
       get the worst flux in the given direction
    **/
    template <CENTERING cent>
    static double
    getFaceError(const vector<EBCM_Face<cent, ebcm_order> >  & a_faces,
                 shared_ptr<ebcm_dou_sca_data>  &              a_phi_exac_ptr,
                 ch_dat_ind                                    a_dit,
                 double                                        a_alpha,
                 double                                        a_beta,
                 int                                           a_sten_rad,
                 int                                           a_weight_power,
                 shared_ptr< ebcm_meta  >                      a_meta_ptr,
                 string                                        a_dom_bc,
                 string                                        a_eb_bc,
                 int                                           a_face_dir,
                 bool                                          a_print)
    {
      double retval = -4586;
      ch_mayday::Error("Doh! not implemented");
      

      return retval;
    }
    ///
    /**
       The points returned are the cells where the bad face happened.
       These will be different on every MPI process.   The max errors returned
       (xface, yface, zface and overall) are  properly computed with MPI scatter and gather.
    **/
    static double
    getMaxAbsFluxError(pr_pt                     & a_worst_xface_pt,
                       pr_pt                     & a_worst_yface_pt,
                       pr_pt                     & a_worst_zface_pt,
                       double                    & a_error_xface, 
                       double                    & a_error_yface, 
                       double                    & a_error_zface, 
                       double                     a_alpha,
                       double                     a_beta,
                       int                        a_sten_rad,
                       int                        a_weight_power,
                       shared_ptr< ebcm_meta  >   a_meta_ptr,
                       string                     a_dom_bc,
                       string                     a_eb_bc,                    
                       bool                       a_print) 

    {

      int num_ghost = 4;
      
      /// For this test, do not need lph_exac is necessary but this reuses code.
      shared_ptr<ebcm_dou_sca_data>  phi_exac_ptr(new ebcm_dou_sca_data(a_meta_ptr, num_ghost));
      shared_ptr<ebcm_dou_sca_data>  lph_exac_ptr(new ebcm_dou_sca_data(a_meta_ptr, num_ghost));
      shared_ptr<ebcm_dou_sca_data> xface_err_ptr(new ebcm_dou_sca_data(a_meta_ptr, num_ghost));
      shared_ptr<ebcm_dou_sca_data> yface_err_ptr(new ebcm_dou_sca_data(a_meta_ptr, num_ghost));
      shared_ptr<ebcm_dou_sca_data> zface_err_ptr(new ebcm_dou_sca_data(a_meta_ptr, num_ghost));
      shared_ptr<ebcm_base_exact> exactly_ptr = getExactSolution();
      setExact(phi_exac_ptr, lph_exac_ptr, a_meta_ptr, exactly_ptr, a_print);

      const auto& grids  = a_meta_ptr->m_grids;
      const auto& graphs = *(a_meta_ptr->m_graphs);
      const auto&    dx  = a_meta_ptr->m_dx;
      auto           dit = grids.dataIterator();
      for(int ibox = 0; ibox < dit.size(); ibox++)
      {
        const auto&  valid =  grids[dit[ibox]];
        const auto&  graph = graphs[dit[ibox]];
        ebcm_subvol_vec volvec(graph, valid, a_print);
        
        for(int ivec = 0; ivec < volvec.size(); ivec++)
        {
          const auto& volume = volvec[ivec];
          double xface_err = -4586;
          double yface_err = -4586;
          double zface_err = -4586;
          xface_err = getFaceError<XFACE>
            (volume.m_xfaces, phi_exac_ptr, dit[ibox],
             a_alpha,a_beta, a_sten_rad, a_weight_power,
             a_meta_ptr, a_dom_bc, a_eb_bc, 0, a_print);
          yface_err = getFaceError<YFACE>
            (volume.m_yfaces,phi_exac_ptr, dit[ibox],
             a_alpha,a_beta, a_sten_rad, a_weight_power,
             a_meta_ptr, a_dom_bc, a_eb_bc, 1, a_print);
                                                 
#if DIM==3
          zface_err = getFaceError<ZFACE>
            (volume.m_zfaces, phi_exac_ptr, dit[ibox],
             a_alpha,a_beta, a_sten_rad, a_weight_power,
             a_meta_ptr, a_dom_bc, a_eb_bc, 2, a_print);
                                                 
#endif
          for(int icell = 0; icell < volume.m_cells.size(); icell++)
          {
            auto pt = volume.m_cells[icell];
            (*xface_err_ptr)[dit[ibox]](pt, 0) = xface_err;
            (*yface_err_ptr)[dit[ibox]](pt, 0) = yface_err;
            (*yface_err_ptr)[dit[ibox]](pt, 0) = zface_err;
          }
        }
      }
      double max_xface_err =  xface_err_ptr->maxNorm();
      double max_yface_err =  yface_err_ptr->maxNorm();

      double max_err_all_dirs = std::max(max_xface_err, max_yface_err);
#if DIM==3
      double max_zface_err =  zface_err_ptr->maxNorm();
      max_err_all_dirs = std::max(max_zface_err, max_err_all_dirs);
#endif       
      return max_err_all_dirs;
    }
    ///
    /**
       Another one of the tests called by  example/EBCM/hoeb_truncation.
    **/ 
    static  void run_open_flux_tests()
    {

      static const int num_ghost = 4;
      
      int    nx_fine      = 4586;
      int    maxGrid      = 4586;    
      int    sten_rad     = 4586; 
      int    weight_power = 4586; 
      double dx_fine      = 4586.;
      double alpha        = 4586.;
      double beta         = 4586.;
      string dom_bc(       "4586");
      string  eb_bc(       "4586");
      bool mergeSmallCells;
      ch_probdom domain_fine;
      bool print = false;

      parseFineGridParams(nx_fine, dx_fine, domain_fine,
                          maxGrid, mergeSmallCells,
                          alpha, beta, sten_rad, weight_power,
                          dom_bc, eb_bc, print);
      
      int    nx_coar = nx_fine/2;
      double dx_coar = dx_fine*2.;
      ch_probdom domain_coar = domain_fine;
      domain_coar.coarsen(2);

      shared_ptr< ebcm_meta  > meta_fine = makeGeometryAtSpecificNX(nx_fine, dx_fine, domain_fine, maxGrid, mergeSmallCells, print);
      shared_ptr< ebcm_meta  > meta_coar = makeGeometryAtSpecificNX(nx_coar, dx_coar, domain_coar, maxGrid, mergeSmallCells, print);

      double max_err_fine = 4586.; 
      double max_err_coar = 4586.;
      /**
         The points returned are the cells where the bad face happened.
         These will be different on every MPI process.   The max errors returned
         (xface, yface, zface and overall) are  properly computed with MPI scatter and gather.
       **/
      pr_pt          worst_xface_fine, worst_xface_coar;
      pr_pt          worst_yface_fine, worst_yface_coar;
      pr_pt          worst_zface_fine, worst_zface_coar;
      double         error_xface_fine, error_xface_coar;
      double         error_yface_fine, error_yface_coar;
      double         error_zface_fine, error_zface_coar;
      
      max_err_fine = getMaxAbsFluxError(worst_xface_fine, worst_yface_fine, worst_zface_fine,
                                        error_xface_fine, error_yface_fine, error_zface_fine,
                                        alpha, beta, sten_rad, weight_power,
                                        meta_fine, dom_bc, eb_bc, print);

      max_err_coar = getMaxAbsFluxError(worst_xface_coar, worst_yface_coar, worst_zface_coar,
                                        error_xface_coar, error_yface_coar, error_zface_coar,
                                        alpha, beta, sten_rad, weight_power,
                                        meta_coar, dom_bc, eb_bc, print);
      
      
      Chombo4::pout() << "flux_truncation: Lmax(X flux coarse error) =  " << error_xface_coar << endl;
      Chombo4::pout() << "flux_truncation: Lmax(X flux fine   error) =  " << error_xface_fine << endl;

      Chombo4::pout() << "flux_truncation: Lmax(Y flux coarse error) =  " << error_yface_coar << endl;
      Chombo4::pout() << "flux_truncation: Lmax(Y flux fine   error) =  " << error_yface_fine << endl;
      
#if DIM==3
      Chombo4::pout() << "flux_truncation: Lmax(Z flux coarse error) =  " << error_zface_coar << endl;
      Chombo4::pout() << "flux_truncation: Lmax(Z flux fine   error) =  " << error_zface_fine << endl;
#endif

      double tol = 1.0e-12;
      double err = std::abs(max_err_fine); 
      if(err > tol)
      {
        
        Chombo4::pout() << "max_err_fine =  =  " << max_err_fine << endl;
        Chombo4::pout() << "max_err_coar =  =  " << max_err_coar << endl;
        double order = (log(max_err_coar/max_err_fine))/log(2.0);
        Chombo4::pout() << "flux_truncation: overall order of accuracy = " <<  order  << endl;
      }
      else
      {
        Chombo4::pout() << "flux_truncation: Fine error is too small for any meaningful order computation." << endl;
        
      }

    }

    
  }; //end elliptic_framwork class
}
#endif //end ifdef ch_use_petsc(none of this stuff makes sense without petsc).

#endif 
