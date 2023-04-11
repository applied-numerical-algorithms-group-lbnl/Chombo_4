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
#include "EBMultigrid.H"

#include "DebugFunctions.H"
#include "Hoeb_ExactSolutions.H"
//this one defines HOEB_MAX_ORDER
#include "Hoeb_Utilities.H"
#include "Hoeb_LAPACKMatrix.H"
#include <iomanip>

typedef Chombo4::GeometryService<      HOEB_MAX_ORDER >   ch_geoserv;
typedef    EBCM::MetaDataLevel<        HOEB_MAX_ORDER >   ebcm_meta;
typedef    EBCM::SubVolumeVector  <    HOEB_MAX_ORDER >   ebcm_subvol_vec;
typedef    EBCM::EBCM_Volu<            HOEB_MAX_ORDER >   ebcm_volume;
typedef    EBCM::EBCM_Graph<           HOEB_MAX_ORDER >   ebcm_graph;
typedef    EBCM::HostLevelData<Real,1, HOEB_MAX_ORDER>    ebcm_leveldata;


typedef Chombo4::Box               ch_box;
typedef Chombo4::DisjointBoxLayout ch_dbl;
typedef Chombo4::ProblemDomain     ch_probdom;
typedef Chombo4::DataIterator      ch_dit;
typedef Chombo4::IntVect           ch_iv;
typedef Chombo4::MayDay            ch_mayday;
/****/
void
makeMergedGeometry(   shared_ptr< ebcm_meta  >  & a_ebcm,
                      ch_dbl                    & a_grids,
                      double                    & a_dx,
                      ch_probdom                & a_domain)
{
  using Chombo4::pout;

  int nx               = 32;
  int maxGrid          = 32;
  bool mergeSmallCells = true;
    
  ParmParse pp;

  pp.get("nx"             , nx);
  pp.get("maxGrid"        , maxGrid);
  pp.get("mergeSmallCells", mergeSmallCells);


  pout() << "nx"        << " = " <<  nx         << endl;
  pout() << "max_grid"  << " = " <<  maxGrid    << endl;
  if(mergeSmallCells)
  {
    pout() << "Cell merging is turned ON."  << endl;
  }
  else
  {
    pout() << "Cell merging is turned OFF."  << endl;
  }
  

  IntVect domLo = IntVect::Zero;
  IntVect domHi  = (nx - 1)*IntVect::Unit;

  Chombo4::ProblemDomain domain(domLo, domHi);

  vector<ch_dbl> vecgrids;
  pout() << "making grids" << endl;
  GeometryService<2>::generateGrids(vecgrids, domain.domainBox(), maxGrid);

  int geomGhost = 6;
  shared_ptr<BaseIF>    impfunc = hoeb::getImplicitFunction();
  pout() << "defining geometry in EB land" << endl;
  Real dx = 1.0/(Real(nx));
  RealVect origin = RealVect::Zero();
  shared_ptr< ch_geoserv > geoserv
    (new ch_geoserv(impfunc, origin, dx, domain.domainBox(), vecgrids, geomGhost));

  
  bool printStuff = true;
  shared_ptr< ebcm_meta  >
    metaDataPtr(new ebcm_meta(geoserv, domain.domainBox(), dx,
                              mergeSmallCells, printStuff));

  string filename("graph_pic.hdf5");
  metaDataPtr->outputGraphMapAsHDF5(filename);
  //smuggle stuff out to make this a little bit useful.
  a_ebcm   = metaDataPtr;
  a_grids  = vecgrids[0];
  a_dx     = dx;
  a_domain = domain;

  
}


///
class  neighborhood_t
{
public:
  
  ch_iv                            m_startiv;
  int                              m_nghost;
  RealVect                         m_startloc;
  shared_ptr<ebcm_subvol_vec>      m_volumes;
  vector<Real>                     m_distance;
  
  neighborhood_t(const ebcm_graph    & a_graph,
                 const ch_iv         & a_start,
                 int a_nghost)
  {
    m_startiv = a_start;
    m_nghost  = a_nghost;
    m_startloc.setToCCLocation(a_start, a_graph.m_dx);

    ch_box region(a_start, a_start);
    region.grow(a_nghost);
    region &= a_graph.m_domain;

    m_volumes = shared_ptr<ebcm_subvol_vec>(new ebcm_subvol_vec(a_graph, region));
    m_distance.resize(m_volumes->size());
    for(int ivec = 0; ivec < m_volumes->size(); ivec++)
    {
      m_distance[ivec] = distanceMetric((*m_volumes)[ivec]);
    }
  }

  //Cartesian makes more sense in this context than anything else (since which cell stuff lives in is undefined)
  Real distanceMetric(const ebcm_volume& a_volu)  const
  {
    RealVect  vectDist = a_volu.m_centroid - this->m_startloc;
    Real distance = vectDist.vectorLength();
    return distance;
  }


private:

  neighborhood_t();
};

///
Real
getInvCondNumber(const ebcm_volume   &  a_volu,
                 const ebcm_graph    &  a_graph)
{
  Real retval;
  //HERE
  ch_mayday::Error("not implemented");
  return retval;
}
///
/**
   Returns 1/invCondNumber
**/
shared_ptr<ebcm_leveldata>
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
    ebcm_subvol_vec allVols(graph, grid);
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

void
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
    ebcm_subvol_vec volumes(graph, valid);
    for(int ivec = 0; ivec < volumes.size(); ivec++)
    {
      const auto& volu = volumes[ivec];
      maxKappaBox = std::max(volu.m_kappa, maxKappaBox);
      minKappaBox = std::min(volu.m_kappa, minKappaBox);
    }
    pout() << "maximum (non covered) volume fraction for grids[" << ibox << "] = " << maxKappaBox << endl;
    pout() << "minimum (non covered) volume fraction for grids[" << ibox << "] = " << minKappaBox << endl;

    maxKappa = std::max(maxKappaBox, maxKappa);
    minKappa = std::min(minKappaBox, minKappa);
  }
  pout() << "maximum volume fraction overall = " << maxKappa << endl;
  pout() << "minimum volume fraction overall = " << minKappa << endl;
}
/****/

int main(int a_argc, char* a_argv[])
{
  //the stuff before the important stuff
  using Chombo4::pout;
#ifdef CH_MPI
  MPI_Init(&a_argc, &a_argv);
  pout() << "MPI INIT called" << std::endl;
#endif

  using Chombo4::pout;
  //needs to be called after MPI_Init
  CH_TIMER_SETFILE("ebcellmerge.time.table");
  {
    
    if (a_argc < 2)
    {
      cerr<< " usage " << a_argv[0] << " <input_file_name> " << endl;
      exit(0);
    }
    char* in_file = a_argv[1];
    ParmParse  pp(a_argc-2,a_argv+2,NULL,in_file);

    //the important stuff
    //see if we can make a geometry
    shared_ptr< ebcm_meta  > ebcm;
    ch_dbl                   grids;
    double                   dx;
    ch_probdom               domain;
    makeMergedGeometry(ebcm, grids, dx, domain);

    Real maxKapp, minKapp;
    //this tries out iteration and checks kappa
    checkKappa(maxKapp, minKapp, ebcm, grids);

    //makes data of inv condition number info and plots it out
    string condition_plot_filename("condition.hdf5");
    pp.query("condition_plot_filename", condition_plot_filename);
    shared_ptr<ebcm_leveldata> data = 
      plotInverseConditionNumberData(ebcm, condition_plot_filename);
  }

  //the stuff after the important stuff
  pout() << "printing time table " << endl;
  CH_TIMER_REPORT();
#ifdef CH_MPI
  pout() << "about to call MPI Finalize" << std::endl;
  MPI_Finalize();
#endif
  return 0;
}
