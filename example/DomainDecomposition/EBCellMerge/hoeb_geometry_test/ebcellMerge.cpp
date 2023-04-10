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
typedef    EBCM::SubVolumeIterator<    HOEB_MAX_ORDER >   ebcm_subvol_it;
typedef Chombo4::DisjointBoxLayout ch_dbl;
typedef Chombo4::ProblemDomain ch_probdom;
typedef Chombo4::DataIterator  ch_dit;
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

/**
**/
void
makeEBCMData(   shared_ptr< ebcm_meta  >  & a_ebcm)
{

  EBCM::HostLevelData<double ,1, HOEB_MAX_ORDER>(a_ebcm, Chombo4::IntVect::Zero);
}

/**
**/
void 

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
    ebcm_subvol_it iterator(graph, valid);
    //Yes the syntax here is weird and old.  it was easier and I am just one guy.
    for(iterator.begin(); iterator.ok(); ++iterator)
    {
      auto volu = *iterator;
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

    makeEBCMData(ebcm);
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
