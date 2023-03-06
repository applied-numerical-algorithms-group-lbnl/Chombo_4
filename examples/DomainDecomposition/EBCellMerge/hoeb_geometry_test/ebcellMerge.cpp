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
#include "Chombo_EBCMGraph.H"
#include "EBMultigrid.H"

#include "DebugFunctions.H"
#include "Hoeb_ExactSolutions.H"
//this one defines HOEB_MAX_ORDER
#include "Hoeb_Utilities.H"
#include "Hoeb_LAPACKMatrix.H"
#include <iomanip>


/****/
unsigned int
makeMergedGeometry()
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

  vector<Chombo4::DisjointBoxLayout> vecgrids;
  pout() << "making grids" << endl;
  GeometryService<2>::generateGrids(vecgrids, domain.domainBox(), maxGrid);

  int geomGhost = 6;
  shared_ptr<BaseIF>    impfunc = hoeb::getImplicitFunction();
  pout() << "defining geometry in EB land" << endl;
  Real dx = 1.0/(Real(nx));
  RealVect origin = RealVect::Zero();
  typedef Chombo4::GeometryService<  HOEB_MAX_ORDER >   ch_geoserv;
  typedef    EBCM::MetaDataLevel<    HOEB_MAX_ORDER >   ebcm_meta;
  shared_ptr< ch_geoserv > geoserv
    (new ch_geoserv(impfunc, origin, dx, domain.domainBox(), vecgrids, geomGhost));

  

  shared_ptr< ebcm_meta  >
    metaDataPtrPtr(new ebcm_meta(geoserv, domain.domainBox(), dx, mergeSmallCells));

  return 0;
}
  
/****/

int main(int a_argc, char* a_argv[])
{
  //the stuff before the important stuff
#ifdef CH_USE_PETSC  
  //because of some kind of solipsistic madness, PetscInitialize calls MPI_INIT
  PetscInt ierr = PetscInitialize(&a_argc, &a_argv, "./.petscrc",PETSC_NULL); CHKERRQ(ierr);
#else
  using Chombo4::pout;
#ifdef CH_MPI
  MPI_Init(&a_argc, &a_argv);
  pout() << "MPI INIT called" << std::endl;
#endif
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
    int retval = makeMergedGeometry();
    pout() << "make mergedGeometry returned "  << retval  << endl;

  }

  //the stuff after the important stuff
  pout() << "printing time table " << endl;
  CH_TIMER_REPORT();
#ifdef CH_MPI
#ifdef CH_USE_PETSC
  pout() << "about to call petsc Finalize" << std::endl;
  PetscFinalize();
#else  
  pout() << "about to call MPI Finalize" << std::endl;
  MPI_Finalize();
#endif
#endif
  return 0;
}
