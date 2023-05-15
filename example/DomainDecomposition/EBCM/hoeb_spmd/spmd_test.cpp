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
typedef Chombo4::BoxIterator                           ch_boxit;
typedef Chombo4::DisjointBoxLayout                     ch_dbl;
typedef Chombo4::DataIterator                          ch_dit;
typedef Chombo4::IntVect                               ch_iv;
typedef Chombo4::MayDay                                ch_mayday;
typedef Proto::RealVect                                pr_rv;
typedef Proto::IndexTM<double, DIM>                    pr_itm_r_dim;
typedef Proto::IndexTM<int , DIM>                      pr_itm_i_dim;
typedef Proto::IndexTM<double, DIM-1>                  pr_itm_r_dmo;
typedef Proto::IndexTM<int , DIM-1>                    pr_itm_i_dmo;
typedef Proto::BaseIF                                  pr_baseif;
typedef ch_eigen::Matrix                               eigen_mat;

///
/**
   This is the testing framework for EBCM SPMD.
**/
namespace EBCM
{
  typedef Chombo4::Box                                   ch_box;
  typedef Chombo4::DisjointBoxLayout                     ch_dbl;
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
    typedef   Proto::MomentIterator<DIM  ,  ebcm_order >   pr_momit_dim;
    typedef   Proto::MomentIterator<DIM-1,  ebcm_order >   pr_momit_dmo;
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
                ch_box          & a_domain,
                double          & a_dx,
                int a_nx, int a_maxGrid,
                bool a_forceSingleGridOnEveryProc,
                string a_prefix, bool a_printStuff)
    {
      ch_iv domLo = ch_iv::Zero;
      ch_iv domHi  = (a_nx - 1)*ch_iv::Unit;

      ch_box dombox;
      a_domain= ch_box(domLo, domHi);
      a_dx = 1./a_nx;
      if(a_forceSingleGridOnEveryProc)
      {
        if(a_printStuff)
        {
          Chombo4::pout() << "For comparison meta, put everything on every proc" << endl;
        }
        int procID = CH4_SPMD::procID();
        vector<ch_box> boxes(1, a_domain);
        vector<int   > procs(1, procID);
        a_grid = ch_dbl(boxes, procs);
      }
      else
      {
        vector<ch_dbl> vecgrids;
        if(a_printStuff)
        {
          Chombo4::pout() << "making grids using GeometryService for distributed data" << endl;
        }
        GeometryService<2>::generateGrids(vecgrids, a_domain, a_maxGrid);
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
                          const ch_box                    & a_domain,
                          int  a_ghost,  int a_nx, 
                          bool a_mergeSmalls, string a_prefix, bool a_printStuff)
    {
      CH_TIME("Algorithm_Framework::makeMergedGeometry");
      using Chombo4::pout;

    
      if(a_printStuff)
      {
        Chombo4::pout() << "making ebcm_meta for " << a_prefix << " case." << endl;
        Chombo4::pout() << "nx"        << " = " <<  a_nx         << endl;
        if(a_mergeSmalls)
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
      if(a_printStuff)
      {
        Chombo4::pout() << "Defining geometry in EB land." << endl;
      }
      double dx = 1.0/(double(a_nx));
      pr_rv origin = pr_rv::Zero();
      shared_ptr< ch_geoserv > geoserv
        (new ch_geoserv(impfunc, origin, a_dx, a_domain, vecgrids, geomGhost));

  
      if(a_printStuff)
      {
        Chombo4::pout() << "Defining geometry in our new EBCM beachhead." << endl;
      }
      shared_ptr< ebcm_meta  >
        metaDataPtr(new ebcm_meta(geoserv, a_domain, a_dx, a_ghost,
                                  a_mergeSmalls, a_printStuff));

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
    ///
    static int
    compareIrregFaces(const vector<IrregFace<ebcm_order> > & a_vec_sing,
                      const vector<IrregFace<ebcm_order> > & a_vec_dist,
                      bool a_cellMerged)
    {
      if(!a_cellMerged)
      {
        int sing_size  = a_vec_sing.size();
        int dist_size  = a_vec_dist.size();
        if(sing_size != dist_size)
        {
          Chombo4::pout() << "EBCM::compareIrregFaces: vector size mismatch" << endl;
          return -1;
        }
        for(int ivec = 0; ivec < sing_size; ivec++)
        {
          const auto& face_sing = a_vec_sing[ivec];
          const auto& face_dist = a_vec_dist[ivec];

          if(face_sing.m_pt != face_dist.m_pt)
          {
            Chombo4::pout() << "EBCM::compareIrregFaces: face mismatch at ivec = " << ivec << endl;
            return -2;
          }
          
          //see if the moments are in the ballpark.
          auto mom_diff = face_sing.m_ebfmom;
          mom_diff     -= face_dist.m_ebfmom;
          double tol = 1.0e-9;
          double biggestDiffest = mom_diff.maxAbs();  //with apologies to Mel Brooks
          if(biggestDiffest > tol)
          {
            Chombo4::pout() << "EBCM::compareIrregFaces: m_facmom mismatch" << endl;
            return -6;
          }
        }
              
      }
      else
      {
        static bool printedOnce = false;
        if(printedOnce)
        {
          Chombo4::pout() << "***EBCM::compareFaces: vector test for merged case not complete***" << endl;
          printedOnce = true;
        }
      }
      return 0;
    }
    ///
    template<CENTERING cent>
    static int
    compareFaces(const vector<EBCM_Face<cent, ebcm_order> > & a_vec_sing,
                 const vector<EBCM_Face<cent, ebcm_order> > & a_vec_dist,
                 bool a_cellMerged)
    {
      if(!a_cellMerged)
      {
        int sing_size  = a_vec_sing.size();
        int dist_size  = a_vec_dist.size();
        if(sing_size != dist_size)
        {
          Chombo4::pout() << "EBCM::compareFaces: vector size mismatch" << endl;
          return -1;
        }
        for(int ivec = 0; ivec < sing_size; ivec++)
        {
          const auto& face_sing = a_vec_sing[ivec];
          const auto& face_dist = a_vec_dist[ivec];
          //this just checks low vof and high vof (not any of the moment information)
          if(face_sing.m_lo != face_dist.m_lo)
          {
            Chombo4::pout() << "EBCM::compareFaces: lo face mismatch at ivec = " << ivec << endl;
            return -7;
          }
          if(face_sing.m_hi != face_dist.m_hi)
          {
            Chombo4::pout() << "EBCM::compareFaces: hi face mismatch at ivec = " << ivec << endl;
            return -2;
          }
          
          //see if the moments are in the ballpark.
          auto mom_diff = face_sing.m_facmom;
          mom_diff     -= face_dist.m_facmom;
          double tol = 1.0e-9;
          double biggestDiffest = mom_diff.maxAbs();  //with apologies to Mel Brooks
          if(biggestDiffest > tol)
          {
            Chombo4::pout() << "EBCM::compareFaces: m_facmom mismatch" << endl;
            return -6;
          }
          
        }
              
      }
      else
      {
        static bool printedOnce = false;
        if(printedOnce)
        {
          Chombo4::pout() << "***EBCM::compareFaces: vector test for merged case not complete***" << endl;
          printedOnce = true;
        }
      }
      return 0;
    }
    ///
    static int
    compareVolumes(const ebcm_volu& a_volu_sing,
                   const ebcm_volu& a_volu_dist,
                   bool a_cellMerged)
    {
      double tol = 1.0e-9;
      double dxdiff = std::abs(a_volu_sing.m_dx - a_volu_dist.m_dx);
      if(a_volu_sing.m_regular != a_volu_dist.m_regular)
      {
        Chombo4::pout() << "EBCM::compareVolumes: m_regular mismatch" << endl;
        return -1;
      }
      if(dxdiff > tol)
      {
        Chombo4::pout() << "EBCM::compareVolumes: m_dx mismatch" << endl;
        return -2;
      }
      int csize_sing = a_volu_sing.m_cells.size();
      int csize_dist = a_volu_dist.m_cells.size();
      //this checks to see if this volume comes from a merger
      if((csize_sing==1)  && (csize_dist==1))//everything should match in this case.
      {
        if(a_volu_sing.m_pt != a_volu_dist.m_pt)
        {
          Chombo4::pout() << "EBCM::compareVolumes: m_pt mismatch" << endl;
          return -3;
        }
        for(int idir = 0; idir < DIM; idir++)
        {
          double centdiff = a_volu_sing.m_centroid[idir] - a_volu_dist.m_centroid[idir];
          if(std::abs(centdiff) > tol)
          {
            Chombo4::pout() << "EBCM::compareVolumes: m_centroid mismatch" << endl;
            return -4;
          }
        }
        if(std::abs(a_volu_sing.m_kappa - a_volu_dist.m_kappa) > tol)
        {
            Chombo4::pout() << "EBCM::compareVolumes: m_kappa mismatch" << endl;
            return -5;
        }

        auto mom_diff = a_volu_sing.m_volmom;
        mom_diff     -= a_volu_dist.m_volmom;
        double tol = 1.0e-9;
        double biggestDiffest = mom_diff.maxAbs();  //with apologies to Mel Brooks
        if(biggestDiffest > tol)
        {
            Chombo4::pout() << "EBCM::compareVolumes: m_volmom mismatch" << endl;
            return -6;
        }
      }
      else
      {
        static bool printed_once = false;
        if(!printed_once)
        {
          Chombo4::pout()  << "***EBCM::compareVolumes: Merged case does not yes test communcation of pt, centroid or kappa***" << endl;
          printed_once = true;
        }
      }
      //if these fail, compareFaces will bark
      int xeek = compareFaces<XFACE>(a_volu_sing.m_xfaces, a_volu_dist.m_xfaces, a_cellMerged); if(!xeek) return xeek;
      int yeek = compareFaces<YFACE>(a_volu_sing.m_yfaces, a_volu_dist.m_yfaces, a_cellMerged); if(!yeek) return yeek;
      int zeek = compareFaces<ZFACE>(a_volu_sing.m_zfaces, a_volu_dist.m_zfaces, a_cellMerged); if(!zeek) return zeek;
      int ieek = compareIrregFaces(  a_volu_sing.m_ifaces, a_volu_dist.m_ifaces, a_cellMerged); if(!ieek) return ieek;
      return 0;
    }
    
    ///see that meta data matches a single-grid distribution
    static int compareMeta(
      const shared_ptr< ebcm_meta  > & a_meta_single_grid,
      const shared_ptr< ebcm_meta  > & a_meta_distributed,
      const ch_dbl                   & a_grid_single_grid,
      const ch_dbl                   & a_grid_distributed,
      const double                   & a_dx,
      const ch_box                   & a_domain,
      const int                      & a_nghost,
      bool a_cellMerged, bool a_printStuff)
    {
      int eekflag = 0;
      ch_datind ind_sing = a_grid_single_grid.dataIterator()[0];
      const auto& ld_graph_sing = *(a_meta_single_grid->m_graphs);
      const auto& ld_graph_dist = *(a_meta_distributed->m_graphs);

      ch_dit    dit_dist = a_grid_distributed.dataIterator();
      for(int ibox = 0; ibox < dit_dist.size(); ibox++)
      {
        ch_datind  ind_dist    =    dit_dist[ibox];
        ch_box   valid_dist    = a_grid_distributed[ind_dist];
        const auto& graph_sing =      ld_graph_sing[ind_sing];
        const auto& graph_dist =      ld_graph_dist[ind_dist];
        ch_box grown = grow(valid_dist, a_nghost);
        grown &= a_domain;
        pr_box prgrown = ProtoCh::getProtoBox(grown);
        for(auto bit = prgrown.begin(); bit != prgrown.end(); ++bit)
        {
          auto pt = *bit;
          if(!graph_sing.isCovered(pt))
          {
            auto volu_sing = graph_sing.getVolumeCoveringPoint(pt);
            auto volu_dist = graph_dist.getVolumeCoveringPoint(pt);
            int eek = compareVolumes(volu_sing, volu_dist, a_cellMerged);
            if(eek  != 0)
            {
              Chombo4::pout() << "EBCM::test_framework::compareMeta problem at point " << pt << endl;
              return eek;
            }
          }
        }
        
      }
      return eekflag;
    }
    ///fill data with known values
    static void fillData(
      ebcm_leveldata                 & a_data,
      const shared_ptr< ebcm_meta  > & a_meta,
      const ch_dbl                   & a_grid,
      const double                   & a_dx,
      const ch_box                   & a_domain,
      const int                      & a_nghost,
      bool a_printStuff)
    {
    }
    /// compare data (ghost and otherwise) 
    static int compareData(
      const ebcm_leveldata           & a_data_single_grid,
      const ebcm_leveldata           & a_data_distributed,
      const shared_ptr< ebcm_meta  > & a_meta_single_grid,
      const shared_ptr< ebcm_meta  > & a_meta_distributed,
      const ch_dbl                   & a_grid_single_grid,
      const ch_dbl                   & a_grid_distributed,
      const double                   & a_dx,
      const ch_box                   & a_domain,
      const int                      & a_nghost,
      bool a_cellMerged, bool a_printStuff)
    {
      int eekflag = 0;
      return eekflag;
    }
    ///compare distribted data vs. single grid data
    static  int runTests()
    {
      CH_TIME("Algorithm_Framework::runTests");

      //define both single-box and multi box (distributed) objects 
      shared_ptr< ebcm_meta  > meta_single_grid;
      shared_ptr< ebcm_meta  > meta_distributed;
      ch_dbl                   grid_single_grid;
      ch_dbl                   grid_distributed;
      double                   dx;
      ch_box                   domain;
      int nghost = 4;
      ParmParse pp("makeMergedGeometry");

      int nx, maxGrid;
      bool mergeSmallCells;
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

      Chombo4::pout() << "Outputting map graphs as hdf5 data." << endl;
      meta_single_grid->outputGraphMapAsHDF5(string("single_grid.map.hdf5"));
      meta_distributed->outputGraphMapAsHDF5(string("distributed.map.hdf5"));

      //compare both single-box and multi box (distributed) objects 
      double max_distributed, min_distributed;
      double max_single_grid, min_single_grid;
      //this tries out iteration and checks kappa
      checkKappa(max_distributed, min_distributed, meta_distributed, grid_distributed, string("distributed"));
      checkKappa(max_single_grid, min_single_grid, meta_single_grid, grid_single_grid, string("single-grid"));

      if(mergeSmallCells)
      {
        Chombo4::pout() << "These might not match!   Merging can cause this." << endl;
      }
      else
      {
        double diffmax = std::abs(max_distributed - max_single_grid);
        double diffmin = std::abs(min_distributed - min_single_grid);
        double tol = 1.0e-8;
        if( (diffmax > tol) || (diffmin > tol) )
        {
          Chombo4::pout() << "Even without merging, we have a serial/parallel mismatch in volume fraction comparison." << endl;
          return -1;
        }
      }

      int eekflag = compareMeta( meta_single_grid, meta_distributed,
                                 grid_single_grid, grid_distributed,
                                 dx, domain, nghost, mergeSmallCells, true );
      
      if(eekflag != 0)
      {
        Chombo4::pout() << "test_framework::compareMeta returned " << eekflag << endl;
        return eekflag;
      }

      ebcm_leveldata data_distributed(meta_distributed, nghost);
      ebcm_leveldata data_single_grid(meta_single_grid, nghost);
      
      fillData(data_distributed, meta_distributed, grid_distributed, dx, domain, nghost, true);
      fillData(data_single_grid, meta_single_grid, grid_single_grid, dx, domain, nghost, true);

      data_distributed.exchange(true);
      data_single_grid.exchange(true);
        
      eekflag = compareData( data_single_grid, data_distributed,
                             meta_single_grid, meta_distributed,
                             grid_single_grid, grid_distributed,
                             dx, domain, nghost, mergeSmallCells, true );

      if(eekflag != 0)
      {
        Chombo4::pout() << "test_framework::compareData returned " << eekflag << endl;
        return eekflag;
      }
      return 0;
    } //end function runTests
  };//end class test_framework
}//end namespace EBCM
int main(int a_argc, char* a_argv[])
{
    int eekflag = 0;
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
      eekflag = EBCM::test_framework<6>::runTests();
    }
    else if(order == 5)
    {
      eekflag = EBCM::test_framework<5>::runTests();
    }
    if(order == 4)
    {
      eekflag = EBCM::test_framework<4>::runTests();
    }
    else if(order == 3)
    {
      eekflag = EBCM::test_framework<3>::runTests();
    }
    else if(order == 2)
    {
      eekflag = EBCM::test_framework<2>::runTests();
    }
    else if(order == 1)
    {
      eekflag = EBCM::test_framework<1>::runTests();
    }
    else
    {
      Chombo4::pout() << "main: unknown order = " << order << endl;
      return -1;
    }
    Chombo4::pout() << "EBCM::test_framework::runTest  returned " << eekflag << endl;
    if(eekflag == 0)
    {
      Chombo4::pout() << "Woo hoo! \n ---H.J.S." << endl;
    }
    else if(eekflag == 1)
    {
      Chombo4::pout() << "Doh!     \n ---H.J.S." << endl;
    }
  }

  Chombo4::pout() << "printing time table " << endl;
  CH_TIMER_REPORT();
#ifdef CH_MPI
  Chombo4::pout() << "about to call MPI Finalize" << std::endl;
  MPI_Finalize();
#endif
  return eekflag;
}
