#ifndef _EBCM_ALGO_DEV_FRAMEWORK__
#define _EBCM_ALGO_DEV_FRAMEWORK__
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
#include "Chombo_EBCM_Graph.H"
#include "Chombo_EBCM_HostLevelData.H"
#include "Chombo_EigenMatrix.H"
#include "Chombo_EBCM_Basics.H"
#include "Chombo_EBCM_Neighborhood.H"
#include <iomanip>

namespace EBCM
{
  typedef Chombo4::Box                                   ch_box;
  typedef Chombo4::DisjointBoxLayout                     ch_dbl;
  typedef Chombo4::ProblemDomain                         ch_probdom;
  typedef Chombo4::IntVect                               ch_iv;
  typedef Chombo4::MayDay                                ch_mayday;
  typedef Proto::RealVect                                pr_rv;
  typedef Proto::IndexTM<double, DIM>                    pr_itm_r_dim;
  typedef Proto::IndexTM<int , DIM>                      pr_itm_i_dim;
  typedef Proto::IndexTM<double, DIM-1>                  pr_itm_r_dmo;
  typedef Proto::IndexTM<int , DIM-1>                    pr_itm_i_dmo;
  typedef Proto::BaseIF                                  pr_baseif;
  typedef Proto::Point                                   pr_pt;
  typedef ch_eigen::Matrix                               eigen_mat;
///
/**
   Purely functional framework for doing tests in EBCM.
   The framework has no data but lots of common  typedefs.
   All functions are static.
*/
  template <int ebcm_order>
  class Algorithm_Framework
  {
  public:
    typedef    EBCM::HostLevelData<double, 1,    ebcm_order >   ebcm_leveldata;
    typedef    EBCM::MetaDataLevel<              ebcm_order >   ebcm_meta;
    typedef    EBCM::SubVolumeVector  <          ebcm_order >   ebcm_subvol_vec;
    typedef    EBCM::EBCM_Volu<                  ebcm_order >   ebcm_volu;
    typedef    EBCM::EBCM_Graph<                 ebcm_order >   ebcm_graph;
    typedef    EBCM::neighborhood<               ebcm_order >   ebcm_neighborhood;
    typedef Chombo4::GeometryService<            ebcm_order >   ch_geoserv;
    typedef   Proto::IndexedMoments<DIM  ,       ebcm_order >   pr_mom_dim;
    typedef   Proto::IndexedMoments<DIM-1,       ebcm_order >   pr_mom_dmo;

    class condition_t
    {
    public:
      double m_inv_condition;
      double m_max_eigenvalue;
      double m_min_eigenvalue;
      pr_pt  m_pt;
      inline void poutLatexTable() const
      {

        using Chombo4::pout;

        ParmParse ppmain("main");
        ParmParse ppsten("stencil");
        int stnrad, wgtpow, taymax;
        ppmain.get("polynomial_order", taymax);
        ppsten.get("weight_power"    , wgtpow);
        ppsten.get("radius"          , stnrad);
      
        pout() << setw(12)
               << setprecision(6)
               << setiosflags(ios::showpoint)
               << setiosflags(ios::scientific);
        pout() << "\\begin{table}" << endl;
        pout() << "\\begin{center}" << endl;
        pout() << "\\begin{tabular}{|ccc|ccc|} \\hline" << endl;

        pout() << " $D$ & $P^T$  & $P^W$ & $R_s$ & \\lambda_{max} & \\lambda_{min}   & $I$ \\\\"<< endl;
        pout() << "\\hline " << endl;
        int dim = DIM;
        pout() <<  "-GREP_HOOK-" << dim << " & " << taymax << " & " << wgtpow << " & " << m_max_eigenvalue << " & " << m_min_eigenvalue << " & "  << m_inv_condition;
        pout() << " \\\\ " << endl;

        pout() << "\\hline " << endl;
        pout() << "\\end{tabular}" << endl;
        pout() << "\\end{center}" << endl;
        pout() << "\\end{table}" << endl;
      }
    };
    ///
    static inline void
    getConditionNumberData(shared_ptr<ebcm_leveldata>      & a_data,
                           shared_ptr<condition_t>         & a_worst,
                           const shared_ptr< ebcm_meta  >  & a_ebcm,
                           const string                    & a_dom_bc,
                           const string                    & a_eb_bc)
    {
    
      CH_TIME("EBCM_Framework::getConditionNumberData");
      a_data  = shared_ptr<ebcm_leveldata>(new ebcm_leveldata(a_ebcm, 0));
      a_worst = shared_ptr<condition_t>(new condition_t());

      double minInvCond = 1.0e30;
    
      const auto& ldgraph = *(a_ebcm->m_graphs);
      const auto&   grids =   a_ebcm->m_grids;
      ch_dit dit = a_ebcm->m_grids.dataIterator();
      auto& mpidata = *(a_data->m_data);
      for(int ibox = 0; ibox < dit.size(); ibox++)
      {
        ///set to zero everywhere so covered cells get something
        auto&     datafab = mpidata[dit[ibox]];
        const auto& graph = ldgraph[dit[ibox]];
        ch_box       grid =   grids[dit[ibox]];
        ebcm_subvol_vec allVols(graph, grid, false);
        for(int ivec = 0; ivec < allVols.size(); ivec++)
        {
          const auto& volu = allVols[ivec];
          auto allpts= volu.m_cells;
          double invCondNum, maxEig, minEig;

          ///info about where the worst one was
          int worst_face_dir, worst_face_index;
          getWorstConditionNumberInfo(invCondNum, maxEig, minEig,
                                      worst_face_dir, worst_face_index,
                                      volu, graph, a_dom_bc, a_eb_bc, false);
          if(invCondNum < minInvCond)
          {
            minInvCond = invCondNum;
            a_worst->m_inv_condition  = minInvCond;
            a_worst->m_min_eigenvalue = minEig;
            a_worst->m_max_eigenvalue = maxEig;
            a_worst->m_pt = volu.m_pt;
          }
          //because of merger, a volume can span multiple cells.
          for(int ipt = 0; ipt < allpts.size(); ipt++)
          {
            datafab(allpts[ipt], 0) = invCondNum;
          }
        }
      }
    } //end function getConditionNumberData
    
    ///gets the worst condition number info for any face in the vof
    static inline void
    getWorstConditionNumberInfo(double & a_invCondNumber,
                                double & a_maxEigenvalue,
                                double & a_minEigenvalue,
                                int    & a_face_dir,
                                int    & a_vec_index,
                                const ebcm_volu     &  a_volu,
                                const ebcm_graph    &  a_graph,
                                const string        &  a_dom_bc,
                                const string        &  a_eb_bc,
                                bool a_print)
    {
      CH_TIME("EBCM::Algorithm_Framework::getWorstConditionNumberInfo");


      //this stuff is probably more general than just this function
      //but the perfect API is eluding me at the moment so we will just
      //give the parameter a generic-ish base.
      ParmParse pp("stencil");
      int stenrad, weightpower;
      pp.get("radius", stenrad);
      pp.get("weight_power", weightpower);
      if(a_volu.m_xfaces.size() > 0)
      {
        double worstCond, maxEigen, minEigen;
        int vecIndex;
          
        getWorstConditionNumberInfo<XFACE>(worstCond, maxEigen, minEigen, vecIndex,
                                           a_volu.m_xfaces, a_graph, stenrad, weightpower,
                                           a_dom_bc, a_eb_bc, a_print);
        a_face_dir       = 0;
        a_invCondNumber = worstCond;
        a_maxEigenvalue = maxEigen;
        a_minEigenvalue = minEigen;
        a_vec_index     = vecIndex;
      }
      if(a_volu.m_yfaces.size() > 0)
      {
        double worstCond, maxEigen, minEigen;
        int vecIndex;
        getWorstConditionNumberInfo<YFACE>(worstCond, maxEigen, minEigen, vecIndex,
                                           a_volu.m_yfaces, a_graph, stenrad, weightpower,
                                           a_dom_bc, a_eb_bc, false);
        if(worstCond < a_invCondNumber)
        {
          a_face_dir       = 1;
          a_invCondNumber = worstCond;
          a_maxEigenvalue = maxEigen;
          a_minEigenvalue = minEigen;
          a_vec_index     = vecIndex;

        }
      }

      
      if(a_volu.m_zfaces.size() > 0)
      {
        double worstCond, maxEigen, minEigen;
        int    vecIndex;
        getWorstConditionNumberInfo<ZFACE>(worstCond, maxEigen, minEigen, vecIndex,
                                           a_volu.m_zfaces, a_graph, stenrad, weightpower, 
                                           a_dom_bc, a_eb_bc, false);
        if(worstCond < a_invCondNumber)
        {
          a_face_dir      = 2;
          a_invCondNumber = worstCond;
          a_maxEigenvalue = maxEigen;
          a_minEigenvalue = minEigen;
          a_vec_index     = vecIndex;
        }

      }

      bool ebFluxKnownAProri = (a_eb_bc == string("neumann"));
      if(  ebFluxKnownAProri)
      {
        static bool printedMessage = false;
        if(!printedMessage)
        {
          Chombo4::pout() << "getWorstConditionNumberInfo: Neumann EB flux needs no system matrix" << endl;
          printedMessage = true;
        }
      }
      else if(a_volu.m_ifaces.size() > 0)
      {
        double worstCond, maxEigen, minEigen;
        int    vecIndex;
        getWorstEBConditionNumberInfo(worstCond, maxEigen, minEigen, vecIndex,
                                      a_volu.m_ifaces, a_graph, stenrad, weightpower, 
                                      a_dom_bc, a_eb_bc, false);
        if(worstCond < a_invCondNumber)
        {
          a_face_dir      = -1;
          a_invCondNumber = worstCond;
          a_maxEigenvalue = maxEigen;
          a_minEigenvalue = minEigen;
          a_vec_index     = vecIndex;
        }

      }

    } //end function vof function



    ///
    /**
       The neighborhood was shoehorned into the interface here for reuse.
     */
    template <CENTERING cent>
    static inline shared_ptr<eigen_mat>
    getSystemMatrix(shared_ptr<ebcm_neighborhood>         &  a_locality,
                    const EBCM_Face<cent, ebcm_order>  &  a_face,
                    const ebcm_graph                   &  a_graph,
                    int a_nghost, int a_weight_power,
                    const string& a_dom_bc, const string& a_eb_bc,
                    bool a_print)
    {
      CH_TIME("EBCM_Framework::getSystemMatrix");
      
      a_locality = shared_ptr<ebcm_neighborhood>
        (new ebcm_neighborhood(a_graph,
                            a_face.m_lo.m_pt,
                            a_face.m_hi.m_pt,
                            a_nghost, a_weight_power, a_print));

      
      shared_ptr<eigen_mat> Wmat_p = getWeightMatrix(a_locality, a_print);
      shared_ptr<eigen_mat> Mmat_p = getMomentMatrix(a_locality, a_dom_bc, a_eb_bc, a_print);

      eigen_mat WMmat;
      multiply(WMmat, *Wmat_p, *Mmat_p);
      eigen_mat WMTmat = WMmat;
      WMTmat.transpose();

      shared_ptr<eigen_mat> Amat_p(new eigen_mat());
      multiply(*Amat_p, WMTmat, WMmat);
    
      if(a_print)
      {
        int icent = cent;
        Chombo4::pout() << "getSystemMatrix: icent = " << icent << endl;
        Chombo4::pout() << "getSystemMatrix: Weight matrix:" << endl;
        Wmat_p->poutDiag();
        Chombo4::pout() << "getSystemMatrix: Moment matrix:" << endl;
        Mmat_p->poutAll();
      
        Chombo4::pout() << "getSystemMatrix: WM   matrix: :" << endl;
        WMmat.poutAll();
      
        Chombo4::pout() << "getSystemMatrix: WM^T matrix: :" << endl;
        WMTmat.poutAll();

        Chombo4::pout() << "System matrix:" << endl;
        Amat_p->poutAll();
      }

      return Amat_p;

    }

    ///
    /**
       The locality gets defined within is passed out for reuse.
     **/
    static inline shared_ptr<eigen_mat>
    getEBSystemMatrix( shared_ptr<ebcm_neighborhood>           & a_locality,
                       const IrregFace< ebcm_order>         &  a_face,
                       const ebcm_graph                     &  a_graph,
                       int a_nghost, int a_weight_power,
                       const string& a_dom_bc, const string & a_eb_bc,
                       bool a_print)
    {
      CH_TIME("EBCM_Framework::getEBSystemMatrix");

      a_locality = shared_ptr<ebcm_neighborhood>
        (new ebcm_neighborhood(a_graph,
                            a_face.m_pt,
                            a_face.m_pt,
                            a_nghost, a_weight_power, a_print));

      shared_ptr<eigen_mat> Wmat_p = getWeightMatrix(a_locality, a_print);
      shared_ptr<eigen_mat> Mmat_p = getMomentMatrix(a_locality, a_dom_bc, a_eb_bc, a_print);

      eigen_mat WMmat;
      multiply(WMmat, *Wmat_p, *Mmat_p);
      eigen_mat WMTmat = WMmat;
      WMTmat.transpose();

      shared_ptr<eigen_mat> Amat_p(new eigen_mat());
      multiply(*Amat_p, WMTmat, WMmat);
    
      if(a_print)
      {
        Chombo4::pout() << "getEBSystemMatrix:Weight matrix:" << endl;
        Wmat_p->poutDiag();
        Chombo4::pout() << "getEBSystemMatrix:Moment matrix:" << endl;
        Mmat_p->poutAll();
      
        Chombo4::pout() << "getEBSystemMatrix:WM   matrix: :" << endl;
        WMmat.poutAll();
      
        Chombo4::pout() << "getEBSystemMatrix:WM^T matrix: :" << endl;
        WMTmat.poutAll();

        Chombo4::pout() << "getEBSystemMatrix:System matrix:" << endl;
        Amat_p->poutAll();
      }

      return Amat_p;

    }
    ///gets the worst condition number info for any face in the vector
    template <CENTERING cent>
    static inline void
    getWorstConditionNumberInfo(double & a_invCondNumber,
                                double & a_maxEigenvalue,
                                double & a_minEigenvalue,
                                int    & a_vec_index,
                                const vector<EBCM_Face<cent, ebcm_order> >  &  a_faces,
                                const ebcm_graph                            &  a_graph,
                                int a_stenrad, int a_weightpower,
                                const string& a_dom_bc, const string& a_eb_bc, bool a_print)
    {
      a_invCondNumber = 4586;
      a_maxEigenvalue = 4586;
      a_minEigenvalue = 4586;
      a_vec_index     = 4586;
      for(int ivec = 0; ivec < a_faces.size(); ivec++)
      {
        double localCond, localMax, localMin;
        shared_ptr<ebcm_neighborhood> locality; //putting this into the interface for reuse
        shared_ptr<eigen_mat> Amat = getSystemMatrix(locality, a_faces[ivec], a_graph, a_stenrad, a_weightpower,
                                                     a_dom_bc, a_eb_bc, a_print);
        Amat->getLargestAndSmallestEigenvalues(localMax, localMin);
        localCond = localMin/localMax;
        if(localCond < a_invCondNumber)
        {
          a_invCondNumber = localCond;
          a_maxEigenvalue = localMax;
          a_minEigenvalue = localMin;
          a_vec_index     = ivec;
        }
      }
      if(a_print)
      {
        Chombo4::pout() << "getWorstCondtionNumberInfo:"<< endl;
        Chombo4::pout() << "a_invCondNumber = " << a_invCondNumber << endl;
        Chombo4::pout() << "a_maxEigenvalue = " << a_maxEigenvalue << endl;
        Chombo4::pout() << "a_minEigenvalue = " << a_minEigenvalue << endl;
      }
    }
    ///same thing but for the cut face.  once again, code reuse eludes me.
    static inline void
    getWorstEBConditionNumberInfo(double & a_invCondNumber,
                                  double & a_maxEigenvalue,
                                  double & a_minEigenvalue,
                                  int    & a_vec_index,
                                  const vector<IrregFace< ebcm_order > >   &  a_faces,
                                  const ebcm_graph                         &  a_graph,
                                  int a_stenrad, int a_weightpower,
                                  const string& a_dom_bc, const string& a_eb_bc, bool a_print)
    {
      a_invCondNumber = 4586;
      a_maxEigenvalue = 4586;
      a_minEigenvalue = 4586;
      a_vec_index     = 4586;
      for(int ivec = 0; ivec < a_faces.size(); ivec++)
      {
        double localCond, localMax, localMin;
        shared_ptr<ebcm_neighborhood>  locality; 
        shared_ptr<eigen_mat> Amat =
          getEBSystemMatrix(locality, a_faces[ivec],
                            a_graph, a_stenrad, a_weightpower,
                            a_dom_bc, a_eb_bc, a_print);
        
        Amat->getLargestAndSmallestEigenvalues(localMax, localMin);
        localCond = localMin/localMax;
        if(localCond < a_invCondNumber)
        {
          a_invCondNumber = localCond;
          a_maxEigenvalue = localMax;
          a_minEigenvalue = localMin;
          a_vec_index     = ivec;
        }
      }
      if(a_print)
      {
        Chombo4::pout() << "getWorstCondtionNumberInfo:"<< endl;
        Chombo4::pout() << "a_invCondNumber = " << a_invCondNumber << endl;
        Chombo4::pout() << "a_maxEigenvalue = " << a_maxEigenvalue << endl;
        Chombo4::pout() << "a_minEigenvalue = " << a_minEigenvalue << endl;
      }
    }
    ///
    /**
       Print out max and min volume fraction in a meta.
    **/
    static void
    checkKappa(double                  & a_maxKappa,
               double                  & a_minKappa,
               shared_ptr< ebcm_meta  >  a_meta,
               const ch_dbl            & a_grids)
    {

      CH_TIME("EBCM_Framework::checkKappa");
      double maxKappa = -1.0e10;
      double minKappa =  1.0e10;
      using Chombo4::pout;
      ch_dit dit = a_grids.dataIterator();
      for(int ibox = 0; ibox < dit.size(); ibox++)
      {
        const auto& valid = a_grids[dit[ibox]];
        const auto& graph = (*(a_meta->m_graphs))[dit[ibox]];
        double maxKappaBox = -1.0e10;
        double minKappaBox =  1.0e10;
        ebcm_subvol_vec volumes(graph, valid, false);
        for(int ivec = 0; ivec < volumes.size(); ivec++)
        {
          const auto& volu = volumes[ivec];
          maxKappaBox = std::max(volu.m_kappa, maxKappaBox);
          minKappaBox = std::min(volu.m_kappa, minKappaBox);
        }
        Chombo4::pout() << "maximum (non covered) volume fraction for grids[" << ibox << "] = " << maxKappaBox << endl;
        Chombo4::pout() << "minimum (non covered) volume fraction for grids[" << ibox << "] = " << minKappaBox << endl;

        maxKappa = std::max(maxKappaBox, maxKappa);
        minKappa = std::min(minKappaBox, minKappa);
      }
      Chombo4::pout() << "maximum volume fraction overall = " << maxKappa << endl;
      Chombo4::pout() << "minimum volume fraction overall = " << minKappa << endl;
    }
    ///test function for example/EBCM/hoeb_petsc_test
    static  void run_hoeb_petsc_tests(shared_ptr< ebcm_meta  >&  a_ebcm)
    {
      CH_TIME("EBCM_Framework::runTests");

      shared_ptr< ebcm_meta  > ebcm;
      ch_dbl                   grids;
      double                   dx;
      ch_probdom               domain;
      int coarseningFactor = 1; //do not coarsen
      makeMergedGeometry(ebcm, grids, dx, domain, coarseningFactor, true, true);

      double maxKapp, minKapp;
      //this tries out iteration and checks kappa
      checkKappa(maxKapp, minKapp, ebcm, grids);

      //makes data of inv condition number info and plots it out
      shared_ptr<ebcm_leveldata> invCondition;
      
      shared_ptr<condition_t>    worst;
      string dom_bc("all_neumann");
      string  eb_bc("all_neumann");
      ParmParse ppmain("main");  //this function is a templated standin for main.
      ppmain.get("dom_bc", dom_bc);
      ppmain.get( "eb_bc",  eb_bc);
      getConditionNumberData(invCondition, worst, ebcm, dom_bc, eb_bc);

      worst->poutLatexTable();
      Chombo4::pout() << "worst condition number at point " << worst->m_pt << endl;
    
      string plot_filename;
      ppmain.get("inv_cond_filename",   plot_filename);
      invCondition->writeToHDF5(plot_filename);

      ///getting meta data out for reuse
      a_ebcm = ebcm;      
    }
    Algorithm_Framework()
    {}    
  

  
  static shared_ptr<pr_baseif> getImplicitFunction()
  {
    CH_TIME("EBCM_Framework::getImplicitFunction");
    shared_ptr<pr_baseif>  retval;
    ParmParse pp("getImplicitFunction");
    string which_geom;
    pp.get("which_geom", which_geom);
    using Chombo4::pout;
    if(which_geom == string("sphere"))
    {
      RealVect center = 0.5*RealVect::Unit();
      double radius = 0.1;
      bool inside = false;
      std::vector<double> centvec;
      pp.get("radius", radius);
      pp.get("inside", inside);
      pp.getarr("center", centvec, 0, DIM);
      for(int idir = 0; idir < DIM; idir++)
      {
        center[idir] = centvec[idir];
      }
      Proto::SimpleSphereIF* sphereptr = new Proto::SimpleSphereIF(center, radius, inside);
      Chombo4::pout() << "sphere implicit function with radius = " << radius << ", center = " << center << ", and inside = " << inside << endl;
      retval = shared_ptr<BaseIF>(static_cast<BaseIF*>(sphereptr));
    }
    else if(which_geom == string("all_regular"))
    {
      Chombo4::pout() << "all regular geometry" << endl;
      retval = shared_ptr<BaseIF>(new Proto::AllRegularIF());
    }
    else if(which_geom == string("plane"))
    {
      using Proto::PlaneIF;
      Chombo4::pout() << "plane implicit function" << endl;
      RealVect normal, startPt;
      vector<double> v_norm, v_start;
      pp.getarr("geom_normal", v_norm, 0, DIM);
      pp.getarr("geom_start_pt", v_start, 0, DIM);
      for(int idir = 0; idir < DIM; idir++)
      {
        normal[ idir] = v_norm[ idir];
        startPt[idir] = v_start[idir];
        Chombo4::pout() << "normal ["<< idir << "] = " << normal [idir]  << endl;
        Chombo4::pout() << "startPt["<< idir << "] = " << startPt[idir]  << endl;
      }
      retval = shared_ptr<BaseIF>(new PlaneIF(startPt, normal));
    }
    else
    {
      Chombo4::pout() << "bogus geometry" << endl;
      Chombo4::MayDay::Error("Doh!");
    }
    return retval;
  }
  

/****/
    static inline void
    makeMergedGeometry(   shared_ptr< ebcm_meta  >  & a_ebcm,
                          ch_dbl                    & a_grids,
                          double                    & a_dx,
                          ch_probdom                & a_domain,
                          int  a_coarseningFactor,
                          bool a_print)
    {
      CH_TIME("EBCM_Framework::makeMergedGeometry");
      using Chombo4::pout;

      int nx               = 32;
      int maxGrid          = 32;
      bool mergeSmallCells = true;
    
      ParmParse pp("makeMergedGeometry");

      pp.get("nx"             , nx);
      pp.get("maxGrid"        , maxGrid);
      pp.get("mergeSmallCells", mergeSmallCells);

      nx /= a_coarseningFactor;

      Chombo4::pout() << "makeMergedGeometry: nx"        << " = " <<  nx         << endl;
      Chombo4::pout() << "makeMergedGeometry: max_grid"  << " = " <<  maxGrid    << endl;
      if(mergeSmallCells)
      {
        Chombo4::pout() << "makeMergedGeometry: Cell merging is turned ON."  << endl;
      }
      else
      {
        Chombo4::pout() << "makeMergedGeometry: Cell merging is turned OFF."  << endl;
      }
  

      ch_iv domLo = ch_iv::Zero;
      ch_iv domHi  = (nx - 1)*ch_iv::Unit;

      ch_probdom domain(domLo, domHi);

      vector<ch_dbl> vecgrids;
      Chombo4::pout() << "makeMergedGeometry: making grids" << endl;
      GeometryService<2>::generateGrids(vecgrids, domain.domainBox(), maxGrid);

      int geomGhost = 6;
      shared_ptr<BaseIF>    impfunc = getImplicitFunction();
      Chombo4::pout() << "makeMergedGeometry: defining geometry in EB land" << endl;
      double dx = 1.0/(double(nx));
      pr_rv origin = pr_rv::Zero();
      shared_ptr< ch_geoserv > geoserv
        (new ch_geoserv(impfunc, origin, dx, domain.domainBox(), vecgrids, geomGhost));

      int ghost = 0;
      shared_ptr< ebcm_meta  >
        metaDataPtr(new ebcm_meta(geoserv, domain.domainBox(), dx, ghost,
                                  mergeSmallCells, a_print));

      a_ebcm   = metaDataPtr;
      a_grids  = vecgrids[0];
      a_dx     = dx;
      a_domain = domain;
    }


    ///
    static inline
    shared_ptr<eigen_mat>
    getWeightMatrix(shared_ptr<ebcm_neighborhood> a_locality, bool a_print)
    {
      int nrows_volumes   = a_locality->m_weig_volumes.size();
      int nrows_dom_faces = a_locality->m_weig_dom_faces.size();
      int nrows_cut_faces = a_locality->m_weig_cut_faces.size();
      int nrows = nrows_volumes + nrows_dom_faces + nrows_cut_faces;
      
      shared_ptr<eigen_mat> weight_p(new eigen_mat(nrows, nrows));
      weight_p->setVal(0.);//this has to be zero because we are only filling in the diagonal.
      int irow_global = 0;
      for(int irow = 0; irow < nrows_volumes; irow++)
      {
        (*weight_p)(irow_global, irow_global) = a_locality->m_weig_volumes[irow];
        irow_global++;
      }
      for(int irow = 0; irow < nrows_dom_faces; irow++)
      {
        (*weight_p)(irow_global, irow_global) = a_locality->m_weig_dom_faces[irow];
        irow_global++;
      }
      for(int irow = 0; irow < nrows_cut_faces; irow++)
      {
        (*weight_p)(irow_global, irow_global) = a_locality->m_weig_cut_faces[irow];
        irow_global++;
      }
      return weight_p;
    }


    ///
    static inline void
    shiftMomentAndFillRow(eigen_mat         & a_mat,
                          const pr_mom_dim  & a_moment,
                          const pr_rv       & a_distance,
                          const int         & a_currentRow,
                          bool a_print)
    {
      CH_TIME("Algorithm_Framework::shiftMomentAndFillRow");
      if(a_print)
      {
        Chombo4::pout() << "shiftMomentAndFillRow: a_currentRow = " << a_currentRow << endl;
        Chombo4::pout() << "shiftMomentAndFillRow: a_distance   = " <<   a_distance << endl;
        Chombo4::pout() << "shiftMomentAndFillRow: input mat row: " << endl;
        a_mat.poutRow(a_currentRow);
        Chombo4::pout() << "shiftMomentAndFillRow:   a_moment =" << endl;
        a_moment.print();
      }

      pr_itm_r_dim itm_dist;
      for(int idir = 0; idir < DIM; idir++)
      {
        itm_dist[idir] = a_distance[idir];
      }

      pr_mom_dim shiftedMom = a_moment;
      shiftedMom.shift(itm_dist);
      if(a_print)
      {
        Chombo4::pout() << "shiftMomentAndFillRow: shiftedMom =" << endl;
        shiftedMom.print();
      }

      double volume = shiftedMom[pr_itm_i_dim::Zero];
      bool divide = (volume > 1.0e-16);
      if(!divide)
      {
        if(a_print)
        {
          Chombo4::pout() << "shiftMomentAndFillRow: smallCellRow with volume = " << volume << endl;
        }
        a_mat.setSmallCellRow(a_currentRow);
      }
      else
      {
        for(MomentIterator<DIM, ebcm_order> momit; momit.ok(); ++momit)
        {
          auto momind = momit();
          unsigned int currentCol = pr_mom_dim::indexOf(momind);
          a_mat(a_currentRow, currentCol) = shiftedMom[momind]/volume;
        }
      }
      if(a_print)
      {
        Chombo4::pout() << "shiftMomentAndFillRow: amat leaving = "  << endl;
        a_mat.poutRow(a_currentRow);
      }
    }

    ///
    static inline void
    fillNeumannDomBCRow(eigen_mat         & a_mat,
                        const pr_mom_dim  & a_moment,
                        const pr_rv       & a_distance,
                        const int         & a_face_dir,
                        const int         & a_currentRow,
                        bool a_print)
    {
      CH_TIME("EBCM::Algorithm_Framework::fillNeumannDomBCRow");
      if(a_print)
      {
        Chombo4::pout() << "fillNeumannDomBCRow: a_currentRow = " << a_currentRow << endl;
        Chombo4::pout() << "fillNeumannDomBCRow: a_distance   = " <<   a_distance << endl;
        Chombo4::pout() << "fillNeumannDomBCRow: input mat row: " << endl;
        a_mat.poutRow(a_currentRow);
        Chombo4::pout() << "fillNeumannDomBCRow:  a_moment =" << endl;
        a_moment.print();
      }

      pr_itm_r_dim itm_dist;
      for(int idir = 0; idir < DIM; idir++)
      {
        itm_dist[idir] = a_distance[idir];
      }

      pr_mom_dim shiftedMom = a_moment;
      shiftedMom.shift(itm_dist);
      if(a_print)
      {
        Chombo4::pout() << "fillNeumannDomBCRow:: shiftedMom =" << endl;
        shiftedMom.print();
      }

      //because this is homogeneous, do not need to divde by area 
      //(so we do not have to check for it being too small)
      for(MomentIterator<DIM, ebcm_order> momit; momit.ok(); ++momit)
      {
        auto momind = momit();
        auto pd = momind[a_face_dir];
        double matval = 0;
        unsigned int currentCol = pr_mom_dim::indexOf(momind);
        if(pd > 0)
        {
          auto derivind  = momind;
          derivind[a_face_dir] -= 1;  //the moment  indexed by p - e^d for neumann bcs
          matval = pd*shiftedMom[derivind];
        }
        a_mat(a_currentRow, currentCol) = matval;
      }

      if(a_print)
      {
        Chombo4::pout() << "fillNeumannDomBCRow: amat leaving = "  << endl;
        a_mat.poutRow(a_currentRow);
      }
    } //end function fillNeuamnnDomBCRow


    ///
    static inline void
    fillDirichletDomBCRow(eigen_mat         & a_mat,
                          const pr_mom_dim  & a_moment,
                          const pr_rv       & a_distance,
                          const int         & a_face_dir,
                          const int         & a_currentRow,
                          bool a_print)
    {
      CH_TIME("EBCM::Algorithm_Framework::fillDirichletDomBCRow");
      if(a_print)
      {
        Chombo4::pout() << "fillDirichletDomBCRow: a_currentRow = " << a_currentRow << endl;
        Chombo4::pout() << "fillDirichletDomBCRow: a_distance   = " <<   a_distance << endl;
        Chombo4::pout() << "fillDirichletDomBCRow: input mat row: " << endl;
        a_mat.poutRow(a_currentRow);
        Chombo4::pout() << "fillDirichletDomBCRow:  a_moment =" << endl;
        a_moment.print();
      }

      pr_itm_r_dim itm_dist;
      for(int idir = 0; idir < DIM; idir++)
      {
        itm_dist[idir] = a_distance[idir];
      }

      pr_mom_dim shiftedMom = a_moment;
      shiftedMom.shift(itm_dist);
      if(a_print)
      {
        Chombo4::pout() << "fillNeumannDomBCRow:: shiftedMom =" << endl;
        shiftedMom.print();
      }
      for(MomentIterator<DIM, ebcm_order> momit; momit.ok(); ++momit)
      {
        auto momind = momit();
        double matval = shiftedMom[momind];
        unsigned int currentCol = pr_mom_dim::indexOf(momind);
        a_mat(a_currentRow, currentCol) = matval;
      }

      if(a_print)
      {
        Chombo4::pout() << "fillDirichletDomBCRow: amat leaving = "  << endl;
        a_mat.poutRow(a_currentRow);
      }
    } //end function fillNeuamnnDomBCRow
    
    ///
    static inline void
    fillNeumannEBBCRow(eigen_mat                      &  a_mat,
                       const IrregFace<ebcm_order>    &  a_cut_face,
                       const pr_rv                    &  a_distance,
                       const int                      &  a_currentRow,
                       bool a_print)
    {
      CH_TIME("EBCM::Algorithm_Framework::fillNeumannEBBCRow");
      if(a_print)
      {
        Chombo4::pout() << "fillNeumannEBBCRow: a_currentRow = " << a_currentRow << endl;
        Chombo4::pout() << "fillNeumannEBBCRow: a_distance   = " <<   a_distance << endl;
        Chombo4::pout() << "fillNeumannEBBCRow: input mat row: " << endl;
        a_mat.poutRow(a_currentRow);
      }

      pr_itm_r_dim itm_dist;
      for(int idir = 0; idir < DIM; idir++)
      {
        itm_dist[idir] = a_distance[idir];
      }

      pr_mom_dim shifted_normal_moments[DIM];
      for(int idir = 0;  idir < DIM; idir++)
      {
        shifted_normal_moments[idir] = a_cut_face.getNormalMoment(idir);
        shifted_normal_moments[idir].shift(itm_dist);
      }

      //because this is homogeneous, do not need to divde by area 
      //(so we do not have to check for it being too small)
      for(MomentIterator<DIM, ebcm_order> momit; momit.ok(); ++momit)
      {
        auto momind = momit();
        unsigned int currentCol = pr_mom_dim::indexOf(momind);

        double momsum = 0;
        for(int idir = 0; idir < DIM; idir ++)
        {
          auto pd        = momind[idir];
          auto derivind  = momind;
          if(pd > 0)
          {
            derivind[idir] -= 1;  //the moment  indexed by p - e^d for neumann bcs
            momsum += pd*shifted_normal_moments[idir][derivind];
          }
        }
        a_mat(a_currentRow, currentCol) = momsum;
      }

      if(a_print)
      {
        Chombo4::pout() << "fillNeumannEBBCRow: amat leaving = "  << endl;
        a_mat.poutRow(a_currentRow);
      }
    } //end function fillNeuamnnEBBCRow

    static inline void
    fillDirichletEBBCRow(eigen_mat                      &  a_mat,
                       const IrregFace<ebcm_order>    &  a_cut_face,
                       const pr_rv                    &  a_distance,
                       const int                      &  a_currentRow,
                       bool a_print)
    {
      CH_TIME("EBCM::Algorithm_Framework::fillDirichletEBCRow");
      if(a_print)
      {
        Chombo4::pout() << "fillDirichletEBBCRow: a_currentRow = " << a_currentRow << endl;
        Chombo4::pout() << "fillDirichletEBBCRow: a_distance   = " <<   a_distance << endl;
        Chombo4::pout() << "fillDirichletEBBCRow: input mat row: " << endl;
        a_mat.poutRow(a_currentRow);
      }

      pr_itm_r_dim itm_dist;
      for(int idir = 0; idir < DIM; idir++)
      {
        itm_dist[idir] = a_distance[idir];
      }

      pr_mom_dim shifted_moments = a_cut_face.m_ebfmom;
      shifted_moments.shift(itm_dist);

      for(MomentIterator<DIM, ebcm_order> momit; momit.ok(); ++momit)
      {
        auto momind = momit();
        unsigned int currentCol = pr_mom_dim::indexOf(momind);

        double momval =  shifted_moments[momind];

        a_mat(a_currentRow, currentCol) = momval;
      }

      if(a_print)
      {
        Chombo4::pout() << "fillDirichletEBBCRow: amat leaving = "  << endl;
        a_mat.poutRow(a_currentRow);
      }
    } //end function fillDirichletEBBCRow
    
    ///
    static shared_ptr<eigen_mat>
    getMomentMatrix(shared_ptr<ebcm_neighborhood> a_locality,
                    const string& a_dom_bc, const string& a_eb_bc,
                    bool a_print)
    {
      CH_TIME("Algorithm_Framework::getMomentMatrix");
      int nrows_volumes   = a_locality->m_weig_volumes.size();
      int nrows_dom_faces = a_locality->m_weig_dom_faces.size();
      int nrows_cut_faces = a_locality->m_weig_cut_faces.size();
      int nrows           = nrows_volumes + nrows_dom_faces + nrows_cut_faces;
      int n_equations     = nrows;
      int n_unknowns      = pr_mom_dim::size();
      if(a_print)
      {
        Chombo4::pout() << "getMomentMatrix: nrows_volumes   = "  << nrows_volumes   << endl;
        Chombo4::pout() << "getMomentMatrix: nrows_dom_faces = "  << nrows_dom_faces << endl;
        Chombo4::pout() << "getMomentMatrix: nrows_cut_faces = "  << nrows_cut_faces << endl;
        Chombo4::pout() << "getMomentMatrix: nrows           = "  << nrows           << endl;
        Chombo4::pout() << "getMomentMatrix: n_equations     = "  << n_equations     << endl;
        Chombo4::pout() << "getMomentMatrix: n_unknowns      = "  << n_unknowns      << endl;
      }
      if(n_unknowns > n_equations)
      {
        Chombo4::pout() << "getMomentMatrix::isufficient equations to solve system" << endl;
        ch_mayday::Error("Doh!");
      }
      int n_rows = n_equations;
      int n_cols = n_unknowns;
      shared_ptr<eigen_mat> Mmat_p(new eigen_mat(n_rows, n_cols));
      Mmat_p->setVal(4586.);
      
      int irow_global = 0;
      for(int irow = 0; irow < nrows_volumes; irow++)
      {
        const auto& volmom   =   (*(a_locality->m_volumes ))[irow].m_volmom;
        pr_rv & xbar =      ( (a_locality->m_xbar_volumes ))[irow];
        shiftMomentAndFillRow(*Mmat_p, volmom, xbar, irow_global, false);
        irow_global++;
      }
      /**/
      if(a_dom_bc == string("neumann"))
      {
        const auto& allfac = (*a_locality->m_dom_faces);
        const auto& allvol = (*a_locality->m_volumes);
        for(int irow = 0; irow < nrows_dom_faces; irow++)
        {
          const auto& curfac = allfac[irow]; //indexing info for domain boundary face

          pr_mom_dmo facmom;
          if(     curfac.m_direction == 0)
          {
            facmom = allvol[curfac.m_subvol_index].m_xfaces[curfac.m_volume_index].m_facmom;
          }
          else if(curfac.m_direction == 1)
          {
            facmom = allvol[curfac.m_subvol_index].m_yfaces[curfac.m_volume_index].m_facmom;
          }
          else if(curfac.m_direction == 2)
          {
            facmom = allvol[curfac.m_subvol_index].m_zfaces[curfac.m_volume_index].m_facmom;
          }
          else
          {
            Chombo4::pout() << "getMomentMatrix: invalid curfac.m_direction" << endl;
            ch_mayday::Error("Doh!");
          }
          
          ///This is the face mom with all derivatives in a_neigh_dir set to zero.
          pr_mom_dim volmomfromfac = ebcm_basics::getVolMomFromFacMom(facmom, curfac.m_direction);
          auto xbar = a_locality->m_xbar_dom_faces[irow];   //need to know xbar to do x-xbar^p
          //auto dist = a_locality->m_dist_dom_faces[irow]; //only used for weight
          //auto weig = a_locality->m_weig_dom_faces[irow]; //only used for weight

          fillNeumannDomBCRow(*Mmat_p, volmomfromfac, xbar, curfac.m_direction, irow_global, a_print);

          irow_global++;
        }
      }
      else if(a_dom_bc == string("dirichlet"))
      {
        const auto& allfac = (*a_locality->m_dom_faces);
        const auto& allvol = (*a_locality->m_volumes);
        for(int irow = 0; irow < nrows_dom_faces; irow++)
        {
          const auto& curfac = allfac[irow]; //indexing info for domain boundary face

          pr_mom_dmo facmom;
          if(     curfac.m_direction == 0)
          {
            facmom = allvol[curfac.m_subvol_index].m_xfaces[curfac.m_volume_index].m_facmom;
          }
          else if(curfac.m_direction == 1)
          {
            facmom = allvol[curfac.m_subvol_index].m_yfaces[curfac.m_volume_index].m_facmom;
          }
          else if(curfac.m_direction == 2)
          {
            facmom = allvol[curfac.m_subvol_index].m_zfaces[curfac.m_volume_index].m_facmom;
          }
          else
          {
            Chombo4::pout() << "getMomentMatrix: invalid curfac.m_direction" << endl;
            ch_mayday::Error("Doh!");
          }
          
          ///This is the face mom with all derivatives in a_neigh_dir set to zero.
          pr_mom_dim volmomfromfac = ebcm_basics::getVolMomFromFacMom(facmom, curfac.m_direction);
          auto xbar = a_locality->m_xbar_dom_faces[irow];   //need to know xbar to do x-xbar^p
          //auto dist = a_locality->m_dist_dom_faces[irow]; //only used for weight
          //auto weig = a_locality->m_weig_dom_faces[irow]; //only used for weight

          fillDirichletDomBCRow(*Mmat_p, volmomfromfac, xbar, curfac.m_direction, irow_global, a_print);

          irow_global++;
        }
        
      }   
      else
      {
        Chombo4::pout() << "getMomentMatrix: invalid dom_bc" << endl;
        ch_mayday::Error("Doh!");
      }
      if(a_eb_bc == string("neumann"))
      {
        for(int irow = 0; irow < nrows_cut_faces; irow++)
        {
          //indexing information to get to the cut face volumes
          const auto &  face_loc = (*(a_locality->m_cut_faces))[irow];  
          const auto &  cut_face = (*(a_locality->m_volumes))[face_loc.m_subvol_index].m_ifaces[face_loc.m_iface_index];
          auto xbar = a_locality->m_xbar_cut_faces[irow];   //need to know xbar to do x-xbar^p
          
          fillNeumannEBBCRow(*Mmat_p, cut_face, xbar, irow_global, a_print);
          irow_global++;
        }
      }
      else if(a_eb_bc == string("dirichlet"))
      {
        for(int irow = 0; irow < nrows_cut_faces; irow++)
        {
          //indexing information to get to the cut face volumes
          const auto &  face_loc = (*(a_locality->m_cut_faces))[irow];  
          const auto &  cut_face = (*(a_locality->m_volumes))[face_loc.m_subvol_index].m_ifaces[face_loc.m_iface_index];
          auto xbar = a_locality->m_xbar_cut_faces[irow];   //need to know xbar to do x-xbar^p
          
          fillDirichletEBBCRow(*Mmat_p, cut_face, xbar, irow_global, a_print);
          irow_global++;
        }
      }
      else
      {
        Chombo4::pout() << "getMomentEquaations: unknown EBBC" << endl;
        ch_mayday::Error("Doh!");
      }
      return Mmat_p;
    }//end function getMomentMatrix


    

  }; //end class Algorithm_Framework 

} //end namespace ebcm

#endif
