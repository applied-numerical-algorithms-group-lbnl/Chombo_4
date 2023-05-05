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
#include "Chombo_EBCM_Algorithm_Framework.H"
#include "Chombo_EigenMatrix.H"
#include "EBMultigrid.H"

#include "DebugFunctions.H"
#include "Hoeb_ExactSolutions.H"

#include "Chombo_LAPACKMatrix.H"
#include <iomanip>


typedef Chombo4::Box                                   ch_box;
typedef Chombo4::DisjointBoxLayout                     ch_dbl;
typedef Chombo4::ProblemDomain                         ch_probdom;
typedef Chombo4::DataIterator                          ch_dit;
typedef Chombo4::IntVect                               ch_iv;
typedef Chombo4::MayDay                                ch_mayday;
//typedef Chombo4::LAPACKMatrix                          ch_mat;
typedef Proto::RealVect                                pr_rv;
typedef Proto::IndexTM<Real, DIM>                      pr_itm_r_dim;
typedef Proto::IndexTM<int , DIM>                      pr_itm_i_dim;
typedef Proto::IndexTM<Real, DIM-1>                    pr_itm_r_dmo;
typedef Proto::IndexTM<int , DIM-1>                    pr_itm_i_dmo;
typedef Proto::BaseIF                                  pr_baseif;
typedef ch_eigen::Matrix                               eigen_mat;

/****/
///
/**
   Print out max and min volume fraction in a meta.
**/
namespace EBCM
{
  typedef Chombo4::Box                                   ch_box;
  typedef Chombo4::DisjointBoxLayout                     ch_dbl;
  typedef Chombo4::ProblemDomain                         ch_probdom;
  typedef Chombo4::IntVect                               ch_iv;
  typedef Chombo4::MayDay                                ch_mayday;
  typedef Proto::RealVect                                pr_rv;
  typedef Proto::IndexTM<Real, DIM>                      pr_itm_r_dim;
  typedef Proto::IndexTM<int , DIM>                      pr_itm_i_dim;
  typedef Proto::IndexTM<Real, DIM-1>                    pr_itm_r_dmo;
  typedef Proto::IndexTM<int , DIM-1>                    pr_itm_i_dmo;
  typedef Proto::BaseIF                                  pr_baseif;
  typedef ch_eigen::Matrix                               eigen_mat;
  
  template<int ebcm_order>
  class test_framework
  {
  public:
    typedef    EBCM::MetaDataLevel<         ebcm_order >   ebcm_meta;
    typedef    EBCM::SubVolumeVector  <     ebcm_order >   ebcm_subvol_vec;
    typedef    EBCM::EBCM_Volu<             ebcm_order >   ebcm_volu;
    typedef    EBCM::EBCM_Graph<            ebcm_order >   ebcm_graph;
    typedef    EBCM::HostLevelData<Real, 1, ebcm_order >   ebcm_leveldata;

    typedef Chombo4::GeometryService<       ebcm_order >   ch_geoserv;
    typedef   Proto::IndexedMoments<DIM  ,  ebcm_order >   pr_mom_dim;
    typedef   Proto::IndexedMoments<DIM-1,  ebcm_order >   pr_mom_dmo;
    typedef   Proto::Point pr_pt;
  
    static shared_ptr<pr_baseif> getImplicitFunction()
    {
      CH_TIME("Algorithm_Framework::getImplicitFunction");
      shared_ptr<pr_baseif>  retval;
      ParmParse pp("getImplicitFunction");
      string which_geom;
      pp.get("which_geom", which_geom);
      using Chombo4::pout;
      if(which_geom == string("sphere"))
      {
        RealVect center = 0.5*RealVect::Unit();
        Real radius = 0.1;
        bool inside = false;
        std::vector<Real> centvec;
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
        Chombo4::MayDay::Error("bogus geometry");
      }
      return retval;
    }
  
    static inline void
    makeMergedGeometry(   shared_ptr< ebcm_meta  >  & a_ebcm,
                          ch_dbl                    & a_grids,
                          double                    & a_dx,
                          ch_probdom                & a_domain,
                          int  a_ghost,
                          bool a_outputGraph,
                          bool a_printStuff)
    {
      CH_TIME("Algorithm_Framework::makeMergedGeometry");
      using Chombo4::pout;

      int nx               = 32;
      int maxGrid          = 32;
      bool mergeSmallCells = true;
    
      ParmParse pp("makeMergedGeometry");

      pp.get("nx"             , nx);
      pp.get("maxGrid"        , maxGrid);
      pp.get("mergeSmallCells", mergeSmallCells);


      Chombo4::pout() << "nx"        << " = " <<  nx         << endl;
      Chombo4::pout() << "max_grid"  << " = " <<  maxGrid    << endl;
      if(mergeSmallCells)
      {
        Chombo4::pout() << "Cell merging is turned ON."  << endl;
      }
      else
      {
        Chombo4::pout() << "Cell merging is turned OFF."  << endl;
      }
  

      ch_iv domLo = ch_iv::Zero;
      ch_iv domHi  = (nx - 1)*ch_iv::Unit;

      ch_probdom domain(domLo, domHi);

      vector<ch_dbl> vecgrids;
      Chombo4::pout() << "making grids" << endl;
      GeometryService<2>::generateGrids(vecgrids, domain.domainBox(), maxGrid);

      int geomGhost = 6;
      shared_ptr<BaseIF>    impfunc = getImplicitFunction();
      Chombo4::pout() << "defining geometry in EB land" << endl;
      Real dx = 1.0/(Real(nx));
      pr_rv origin = pr_rv::Zero();
      shared_ptr< ch_geoserv > geoserv
        (new ch_geoserv(impfunc, origin, dx, domain.domainBox(), vecgrids, geomGhost));

  
      shared_ptr< ebcm_meta  >
        metaDataPtr(new ebcm_meta(geoserv, domain.domainBox(), dx, a_ghost,
                                  mergeSmallCells, a_printStuff));

      string graph_filename;
      pp.get("graph_filename"         , graph_filename);
      metaDataPtr->outputGraphMapAsHDF5(graph_filename);
    
      //smuggle stuff out to make this a little bit useful.
      a_ebcm   = metaDataPtr;
      a_grids  = vecgrids[0];
      a_dx     = dx;
      a_domain = domain;
    }

    static void
    checkKappa(double                  & a_maxKappa,
               double                  & a_minKappa,
               shared_ptr< ebcm_meta  >  a_meta,
               const ch_dbl            & a_grids)
    {

      CH_TIME("Algorithm_Framework::checkKappa");
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


    static  void runTests()
    {
      CH_TIME("Algorithm_Framework::runTests");

      //the important stuff
      //see if we can make a geometry
      shared_ptr< ebcm_meta  > ebcm;
      ch_dbl                   grids;
      double                   dx;
      ch_probdom               domain;
      int nghost = 0;
      makeMergedGeometry(ebcm, grids, dx, domain, nghost, true, true);

      Real maxKapp, minKapp;
      //this tries out iteration and checks kappa
      checkKappa(maxKapp, minKapp, ebcm, grids);
    }
  };//end class test_framework
}//end namespace EBCM
int main(int a_argc, char* a_argv[])
{
#ifdef CH_MPI
  MPI_Init(&a_argc, &a_argv);
  Chombo4::pout() << "MPI INIT called" << std::endl;
#endif

  using Chombo4::pout;
  CH_TIMER_SETFILE("ebcm_framework_time_table");
  {
    
    if (a_argc < 2)
    {
      cerr<< " usage " << a_argv[0] << " <input_file_name> " << endl;
      return -2;
    }
    char* in_file = a_argv[1];
    //We need to declare the first  parmparse this way
    //(the reason for this is lost in time).
    ParmParse  ppdecl(a_argc-2,a_argv+2,NULL,in_file);

    int order;
    ParmParse pp("main");
    pp.get("polynomial_order", order);
    Chombo4::pout() <<  "Running test for polynomial order = " << order << endl;
    if(order == 6)
    {
      EBCM::test_framework<6>::runTests();
    }
    else if(order == 5)
    {
      EBCM::test_framework<5>::runTests();
    }
    if(order == 4)
    {
      EBCM::test_framework<4>::runTests();
    }
    else if(order == 3)
    {
      EBCM::test_framework<3>::runTests();
    }
    else if(order == 2)
    {
      EBCM::test_framework<2>::runTests();
    }
    else if(order == 1)
    {
      EBCM::test_framework<1>::runTests();
    }
    else
    {
      Chombo4::pout() << "main: unknown order = " << order << endl;
      return -1;
    }
  }

  Chombo4::pout() << "printing time table " << endl;
  CH_TIMER_REPORT();
#ifdef CH_MPI
  Chombo4::pout() << "about to call MPI Finalize" << std::endl;
  MPI_Finalize();
#endif
  return 0;
}
