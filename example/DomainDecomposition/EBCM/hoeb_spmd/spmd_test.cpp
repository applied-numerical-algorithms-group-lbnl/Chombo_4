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
typedef Proto::doubleVect                                pr_rv;
typedef Proto::IndexTM<double, DIM>                      pr_itm_r_dim;
typedef Proto::IndexTM<int , DIM>                      pr_itm_i_dim;
typedef Proto::IndexTM<double, DIM-1>                    pr_itm_r_dmo;
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
  typedef Proto::IndexTM<double, DIM>                    pr_itm_r_dim;
  typedef Proto::IndexTM<int ,   DIM>                    pr_itm_i_dim;
  typedef Proto::IndexTM<double, DIM-1>                  pr_itm_r_dmo;
  typedef Proto::IndexTM<int ,   DIM-1>                  pr_itm_i_dmo;
  typedef Proto::BaseIF                                  pr_baseif;
  typedef ch_eigen::Matrix                               eigen_mat;
  
  template<int ebcm_order>
  class test_framework
  {
  public:
    typedef    EBCM::MetaDataLevel<           ebcm_order >   ebcm_meta;
    typedef    EBCM::SubVolumeVector  <       ebcm_order >   ebcm_subvol_vec;
    typedef    EBCM::EBCM_Volu<               ebcm_order >   ebcm_volu;
    typedef    EBCM::EBCM_Graph<              ebcm_order >   ebcm_graph;
    typedef    EBCM::HostLevelData<double, 1, ebcm_order >   ebcm_leveldata;

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
        pr_rv center = 0.5*pr_rv::Unit();
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
        pr_rv normal, startPt;
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

    //if (a_forceSingleGridOnProc){ EVERY PROC WILL HOLD THE WHOLE DOMAIN.  THIS WILL NOT SCALE WELL.}
    static inline void
    getGridInfo(ch_dbl          & a_grid,
                ch_probdom      & a_domain,
                double          & a_dx,
                int a_nx, int a_maxGrid,
                bool a_forceSingleGridOnEveryProc,
                string a_prefix, bool a_printStuff)
    {
      ch_iv domLo = ch_iv::Zero;
      ch_iv domHi  = (nx - 1)*ch_iv::Unit;

      a_domain= ch_probdom(domLo, domHi);
      a_dx = 1./a_nx;
      if(a_forceSingleGridOnEveryProc)
      {
        if(a_printStuff)
        {
          Chombo4::pout() << "For comparison meta, put everything on every proc" << endl;
        }
        int procID = CH4_SPMD::procID();
        vector<ch_box> boxes(1, a_domain);
        vector<ch_box> procs(1, procID);
        a_grid = ch_dbl(boxes, procs);
      }
      else
      {
        vector<ch_dbl> vecgrids;
        if(a_printStuff)
        {
          Chombo4::pout() << "making grids using GeometryService for distributed data" << endl;
        }
        GeometryService<2>::generateGrids(vecgrids, domain.domainBox(), maxGrid);
        a_grid = vecgrids[0];
      }
      if(a_printStuff)
      {
        Chombo4::pout() << "grid = " << a_grid << endl;
      }
    }
    static inline void
    makeMergedGeometry(   shared_ptr< ebcm_meta  >        & a_ebcm,
                          const ch_dbl                    & a_grids,
                          const double                    & a_dx,
                          const ch_probdom                & a_domain,
                          int  a_ghost,  int a_nx, 
                          bool a_mergeSmalls, string a_prefix, bool a_printStuff)
    {
      CH_TIME("Algorithm_Framework::makeMergedGeometry");
      using Chombo4::pout;

    
      int nx               = a_nx;
      int maxGrid          = a_maxGrid;
      bool mergeSmallCells = a_mergeSmalls;

      if(a_printStuff)
      {
        Chombo4::pout() << "making ebcm_meta for " << a_prefix << " case." << endl;
        Chombo4::pout() << "nx"        << " = " <<  nx         << endl;
        if(mergeSmallCells)
        {
          Chombo4::pout() << "Cell merging is turned ON."  << endl;
        }
        else
        {
          Chombo4::pout() << "Cell merging is turned OFF."  << endl;
        }
      }

      vector<ch_dbl> vecgrids(1, a_grids);
      int geomGhost = 6;
      shared_ptr<BaseIF>    impfunc = getImplicitFunction();
      Chombo4::pout() << "defining geometry in EB land" << endl;
      double dx = 1.0/(double(nx));
      pr_rv origin = pr_rv::Zero();
      shared_ptr< ch_geoserv > geoserv
        (new ch_geoserv(impfunc, origin, dx, domain.domainBox(), vecgrids, geomGhost));

  
      shared_ptr< ebcm_meta  >
        metaDataPtr(new ebcm_meta(geoserv, domain.domainBox(), dx, a_ghost,
                                  mergeSmallCells, a_printStuff));

      a_ebcm   = metaDataPtr;
    }

    static void
    checkKappa(double                  & a_maxKappa,
               double                  & a_minKappa,
               shared_ptr< ebcm_meta  >  a_meta,
               const ch_dbl            & a_grids,
               string a_prefix)
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
        Chombo4::pout() << a_prefix << " maximum (non covered) volume fraction for grids[" << ibox << "] = " << maxKappaBox << endl;
        Chombo4::pout() << a_prefix << " minimum (non covered) volume fraction for grids[" << ibox << "] = " << minKappaBox << endl;

        maxKappa = std::max(maxKappaBox, maxKappa);
        minKappa = std::min(minKappaBox, minKappa);
      }
      Chombo4::pout() << a_prefix << " maximum volume fraction overall = " << maxKappa << endl;
      Chombo4::pout() << a_prefix << " minimum volume fraction overall = " << minKappa << endl;
    }

    ///see that meta data matches a single-grid distribution
    static void compareMeta(
      const shared_ptr< ebcm_meta  > & a_meta_single_grid,
      const shared_ptr< ebcm_meta  > & a_meta_distributed,
      const ch_dbl                   & a_grid_single_grid,
      const ch_dbl                   & a_grid_distributed,
      const double                   & a_dx,
      const ch_probdom               & a_domain,
      const int                      & a_nghost )
    {
    }
    ///fill data with known values
    static void fillData(
      ebcm_leveldata                 & a_data
      const shared_ptr< ebcm_meta  > & a_meta,
      const ch_dbl                   & a_grid,
      const double                   & a_dx,
      const ch_probdom               & a_domain,
      const int                      & a_nghost )
    {
    }
    /// compare data (ghost and otherwise) 
    static void compareData(
      const ebcm_leveldata           & a_data_single_grid,
      const ebcm_leveldata           & a_data_distributed,
      const shared_ptr< ebcm_meta  > & a_meta_single_grid,
      const shared_ptr< ebcm_meta  > & a_meta_distributed,
      const ch_dbl                   & a_grid_single_grid,
      const ch_dbl                   & a_grid_distributed,
      const double                   & a_dx,
      const ch_probdom               & a_domain,
      const int                      & a_nghost )
    {
    }
    ///compare distribted data vs. single grid data
    static  void spmd_tests()
    {
      CH_TIME("Algorithm_Framework::runTests");

      //the important stuff
      shared_ptr< ebcm_meta  > meta_single_grid;
      shared_ptr< ebcm_meta  > meta_distributed;
      ch_dbl                   grid_single_grid;
      ch_dbl                   grid_distributed;
      double                   dx;
      ch_probdom               domain;
      int nghost = 4;
      ParmParse pp("makeMergedGeometry");


      pp.get("nx"             , nx);
      pp.get("maxGrid"        , maxGrid);
      pp.get("mergeSmallCells", mergeSmallCells);
      int maxMaxGrid = 4096;

      getGridInfo(grid_distributed, domain, dx, nx,  maxGrid   , false, string("distributed"), true);
      getGridInfo(grid_single_grid, domain, dx, nx,  maxMaxGrid, false, string("single_grid"), true);

      makeMergedGeometry( meta_distributed,
                          grid_distributed,
                          dx, domain, nghost,
                          nx, mergeSmallCells,
                          string("distributed"), true);

      makeMergedGeometry( meta_single_grid,
                          grid_single_grid,
                          dx, domain, nghost,
                          nx,  mergeSmallCells,
                          string("single_grid"), true);


      double maxKapp, minKapp;
      //this tries out iteration and checks kappa
      checkKappa(maxKapp, minKapp, meta_distributed, grid_distributed, string("distributed"));
      checkKappa(maxKapp, minKapp, meta_single_grid, grid_single_grid, string("single-grid"));

      compareMeta( meta_single_grid, meta_distributed,
                   grid_single_grid, grid_distributed,
                   dx,        domain,        nghost );

      ebcm_leveldata data_distributed(meta_distributed, nghost);
      ebcm_leveldata data_single_grid(meta_single_grid, nghost);
      
      fillData(data_distributed, meta_distributed, grid_distributed, dx, domain, nghost);
      fillData(data_single_grid, meta_single_grid, grid_single_grid, dx, domain, nghost);

      data_distributed.exchange();
      data_single_grid.exchange();
        
      compareData( data_single_grid, data_distributed,
                   meta_single_grid, meta_distributed,
                   grid_single_grid, grid_distributed,
                   dx,        domain,        nghost );

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
  CH_TIMER_SETFILE("spmd.time.table");
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
