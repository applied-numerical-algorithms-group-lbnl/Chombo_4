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

///
/**
   Purely functional framework for doing tests in EBCM.
   The framework has no data but lots of common  typedefs.
   All functions are static.
 */
template <int ebcm_order>
class EBCM_Framework
{
public:
  EBCM_Framework()
  {}    
  

  typedef    EBCM::MetaDataLevel<         ebcm_order >   ebcm_meta;
  typedef    EBCM::SubVolumeVector  <     ebcm_order >   ebcm_subvol_vec;
  typedef    EBCM::EBCM_Volu<             ebcm_order >   ebcm_volu;
  typedef    EBCM::EBCM_Graph<            ebcm_order >   ebcm_graph;
  typedef    EBCM::HostLevelData<Real, 1, ebcm_order >   ebcm_leveldata;

  typedef Chombo4::GeometryService<       ebcm_order >   ch_geoserv;
  typedef Proto::IndexedMoments<DIM  ,    ebcm_order >   pr_mom_dim;
  typedef Proto::IndexedMoments<DIM-1,    ebcm_order >   pr_mom_dmo;

  
  static shared_ptr<pr_baseif> getImplicitFunction()
  {
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
/****/
  static inline void
  makeMergedGeometry(   shared_ptr< ebcm_meta  >  & a_ebcm,
                        ch_dbl                    & a_grids,
                        double                    & a_dx,
                        ch_probdom                & a_domain,
                        bool a_outputGraph,
                        bool a_printStuff)
  {
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
      metaDataPtr(new ebcm_meta(geoserv, domain.domainBox(), dx,
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


  ///
  struct  neighborhood_t
  {
  
    ch_iv                            m_startiv;
    int                              m_nghost;
    pr_rv                            m_startloc;
    shared_ptr<ebcm_subvol_vec>      m_volumes;
    vector<Real>                     m_distance;
    vector<Real>                     m_eqweight;
  
    neighborhood_t(const ebcm_graph    & a_graph,
                   const ch_iv         & a_start,
                   int a_nghost, int a_weightPower, bool a_printStuff)
    {
      m_startiv = a_start;
      m_nghost  = a_nghost;
      m_startloc.setToCCLocation(a_start, a_graph.m_dx);

      ch_box region(a_start, a_start);
      region.grow(a_nghost);
      region &= a_graph.m_domain;

      m_volumes = shared_ptr<ebcm_subvol_vec>(new ebcm_subvol_vec(a_graph, region, a_printStuff));
      //begin debug
      Chombo4::pout() << "first volmom in the neighborhood" << endl;
      (*m_volumes)[0].m_volmom.print();
      //end debug
      int n_equations = m_volumes->size();
      int n_unknowns  = pr_mom_dim::size();
      if(n_unknowns > n_equations)
      {
        Chombo4::pout() << a_nghost <<
          " ghost cells not getting enough volumes with start = "
               << a_start << endl;
        ch_mayday::Error();
      }
      m_distance.resize(m_volumes->size());
      m_eqweight.resize(m_volumes->size());
    
      for(int ivec = 0; ivec < m_volumes->size(); ivec++)
      {
        Real xbar = distanceMetric((*m_volumes)[ivec]);
      
        //putting this check because someone is going to hack the distance function      
        if(xbar < 1.0e-16) 
        {
          if(xbar < 0)
            Chombo4::pout() << "neighborhood_t::distanceMetric returned a negtive number." << endl;
          else
            Chombo4::pout() << "neighborhood_t::distanceMetric returned too small a positive number." << endl;
          ch_mayday::Error();
        }
        m_distance[ivec] = xbar;
        m_eqweight[ivec] = 1;
        for(int iweight = 0; iweight < a_weightPower; iweight++)
        {
          m_eqweight[ivec] /= xbar;
        }
      }
    }


    //Cartesian makes more sense in this context than anything else (since which cell stuff lives in is undefined)
    Real distanceMetric(const ebcm_volu & a_volu)  const
    {
      pr_rv  vectDist = a_volu.m_centroid - this->m_startloc;
      Real distance = vectDist.vectorLength();
      return distance;
    }

  };
  ///
  static inline
  shared_ptr<eigen_mat>
  getWeightMatrix(shared_ptr<neighborhood_t> a_locality, bool a_printStuff)
  {
    int nrows = a_locality->m_eqweight.size();
    shared_ptr<eigen_mat> weight_p(new eigen_mat(nrows, nrows));
    weight_p->setVal(0.);
    for(int irow = 0; irow < nrows; irow++)
    {
      (*weight_p)(irow, irow) = a_locality->m_eqweight[irow];
    }
    return weight_p;
  }


  ///
  static inline void
  shiftMomentAndFillRow(eigen_mat         & a_mat,
                        const pr_mom_dim  & a_moment,
                        const pr_rv       & a_distance,
                        const int         & a_currentRow,
                        bool a_printStuff)
  {
    if(a_printStuff)
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
    if(a_printStuff)
    {
      Chombo4::pout() << "shiftMomentAndFillRow: shiftedMom =" << endl;
      shiftedMom.print();
    }

    Real volume = shiftedMom[pr_itm_i_dim::Zero];
    bool divide = (volume > 1.0e-16);
    if(!divide)
    {
      if(a_printStuff)
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
    if(a_printStuff)
    {
      Chombo4::pout() << "shiftMomentAndFillRow: amat leaving = "  << endl;
      a_mat.poutRow(a_currentRow);
    }
  }

  ///
  static shared_ptr<eigen_mat>
  getMomentMatrix(shared_ptr<neighborhood_t> a_locality,
                  bool a_printStuff)
  {
    //number of equations
    int n_equations = a_locality->m_volumes->size();
    int n_unknowns  = pr_mom_dim::size();
    if(n_unknowns > n_equations)
    {
      ch_mayday::Error("getMomentMatrix::bring me more equations");
    }
    int n_rows = n_equations;
    int n_cols = n_unknowns;
    
    shared_ptr<eigen_mat> Mmat_p(new eigen_mat(n_rows, n_cols));
    for(int irow = 0; irow < n_rows; irow++)
    {
      const auto& volmom   =   (*(a_locality->m_volumes ))[irow].m_volmom;
      if(a_printStuff)
      {
        volmom.print();
      }
      const auto& distance =   ( (a_locality->m_distance))[irow];
      shiftMomentAndFillRow(*Mmat_p, volmom, distance, irow, a_printStuff);
    }
    return Mmat_p;
  }

  ///
  static inline shared_ptr<eigen_mat>
  getAMatrix(const ebcm_volu     &  a_volu,
             const ebcm_graph    &  a_graph,
             int a_nghost, int a_weightPower,
             bool a_printStuff)
  {
  
    shared_ptr<neighborhood_t>
      locality(new neighborhood_t(a_graph, a_volu.m_pt, a_nghost, a_weightPower, a_printStuff));
  
    shared_ptr<eigen_mat> Wmat_p = getWeightMatrix(locality, a_printStuff);
    shared_ptr<eigen_mat> Mmat_p = getMomentMatrix(locality, a_printStuff);

    eigen_mat WMmat;
    multiply(WMmat, *Wmat_p, *Mmat_p);
    eigen_mat WMTmat = WMmat;
    WMTmat.transpose();
  
    shared_ptr<eigen_mat> Amat_p(new eigen_mat());
    multiply(*Amat_p, WMTmat, WMmat);
    
    if(a_printStuff)
    {
      Chombo4::pout() << "Weight matrix: " << endl;
      Wmat_p->poutAll();
      Chombo4::pout() << "Moment matrix: " << endl;
      Mmat_p->poutAll();
      
      Chombo4::pout() << "WM   matrix: " << endl;
      WMmat.poutAll();
      
      Chombo4::pout() << "WM^T matrix: " << endl;
      WMTmat.poutAll();

      Chombo4::pout() << "A matrix: " << endl;
      Amat_p->poutAll();
    }

    return Amat_p;
  }
  ///
  static inline Real
  getInvCondNumber(const ebcm_volu     &  a_volu,
                   const ebcm_graph    &  a_graph)
  {
    //this stuff is probably more general than just this function
    //but the perfect API is eluding me at the moment so we will just
    //give the parameter a generic-ish base.
    ParmParse pp("stencil");
    int stenrad, weightpower;
    pp.get("radius", stenrad);
    pp.get("weight_power", weightpower);
    //just print everything for the first one
    static bool printedAlready = false;
    shared_ptr<eigen_mat> Amat = getAMatrix(a_volu, a_graph, stenrad, weightpower, !printedAlready);
    printedAlready = true;
  
    Real retval = Amat->invConditionNumber();

    return retval;
  }
  ///
  static inline shared_ptr<ebcm_leveldata>
  plotInverseConditionNumberData(const shared_ptr< ebcm_meta  >  & a_ebcm,
                                 const string                    & a_filename)
  {
    shared_ptr<ebcm_leveldata>
      data(new ebcm_leveldata(a_ebcm, Chombo4::IntVect::Zero));
  
    const auto& ldgraph = *(a_ebcm->m_graphs);
    const auto&   grids =   a_ebcm->m_grids;
    ch_dit dit = a_ebcm->m_grids.dataIterator();
    auto& mpidata = *(data->m_data);
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
        Real invCondNum = getInvCondNumber(volu, graph);
        //because of merger, a volume can span multiple cells.
        for(int ipt = 0; ipt < allpts.size(); ipt++)
        {
          datafab(allpts[ipt], 0) = invCondNum;
        }
      }
    }
    data->writeToHDF5(a_filename);
    return data;
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


  static void runEigenTests()
  {

    int n =4;
    eigen_mat A(n, n);
    Real val = 1;

    for(int irow = 0; irow < n; irow++)
    {
      for(int icol = 0; icol < n; icol++)
      {
        bool odd = (irow%2 == 1);
        if(odd && (irow== icol))
          A(irow, icol) = -val;
        else
          A(irow, icol) =  val;
        val += 1;
      }
    }

    eigen_mat Ainv = A;
    Ainv.invert();


    eigen_mat AAinv;
    multiply(AAinv, A, Ainv);
    using Chombo4::pout;
    pout() << "eigen inverse test" << endl;
    pout() << "A = " << endl;
    A.poutAll();
    pout() << endl;

    pout() << "Ainv = " << endl;
    Ainv.poutAll();
    pout() << endl;

    pout() << "A * Ainv = " << endl;
    AAinv.poutAll();
    pout() << endl;

  }
  static  void runTests()
  {
    //linear algebra machinery check
    runEigenTests();
    //the important stuff
    //see if we can make a geometry
    shared_ptr< ebcm_meta  > ebcm;
    ch_dbl                   grids;
    double                   dx;
    ch_probdom               domain;
    makeMergedGeometry(ebcm, grids, dx, domain, true, true);

    Real maxKapp, minKapp;
    //this tries out iteration and checks kappa
    checkKappa(maxKapp, minKapp, ebcm, grids);

    //makes data of inv condition number info and plots it out
    string plot_filename;
    ParmParse pp("main");  //kinda pendantic to use anything else
    pp.get("inv_cond_filename", plot_filename);
    shared_ptr<ebcm_leveldata> data = 
      plotInverseConditionNumberData(ebcm, plot_filename);
  }
};
/****/

int main(int a_argc, char* a_argv[])
{
#ifdef CH_MPI
  MPI_Init(&a_argc, &a_argv);
  Chombo4::pout() << "MPI INIT called" << std::endl;
#endif

  using Chombo4::pout;
  CH_TIMER_SETFILE("ebcellmerge.time.table");
  {
    
    if (a_argc < 2)
    {
      cerr<< " usage " << a_argv[0] << " <input_file_name> " << endl;
      return -2;
    }
    char* in_file = a_argv[1];
    //We need to declare the first  parmparse this way because that is the way PP works.  
    ParmParse  ppdecl(a_argc-2,a_argv+2,NULL,in_file);

    int order;
    ParmParse pp("main");
    pp.get("polynomial_order", order);
    Chombo4::pout() <<  "Running test for polynomial order = " << order << endl;
    if(order == 4)
    {
      EBCM_Framework<4>::runTests();
    }
    else if(order == 3)
    {
      EBCM_Framework<3>::runTests();
    }
    else if(order == 2)
    {
      EBCM_Framework<2>::runTests();
    }
    else if(order == 1)
    {
      EBCM_Framework<1>::runTests();
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
