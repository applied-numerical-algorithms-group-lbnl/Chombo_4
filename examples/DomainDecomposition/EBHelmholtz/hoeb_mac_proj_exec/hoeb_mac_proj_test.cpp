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
#include "Chombo_EBEncyclopedia.H"
#include "Chombo_EBDictionary.H"
#include "Chombo_EBChombo.H"
#include "Chombo_EBLevelFluxData.H"
#include "Hoeb_MAC_Projector.H"
#include "EBMACProjector.H"
#include "SetupFunctions.H"

#include <iomanip>

#define MAX_ORDER 4

using std::cout;
using std::endl;
using std::shared_ptr;
using Proto::Var;
using Proto::SimpleEllipsoidIF;

inline 
EBIBC getIBCs()
{
  string veloIC;
  string scalIC("blob");
  string loDomBC[DIM];
  string hiDomBC[DIM];
  ParmParse pp;
  pp.get("interior_velo", veloIC);

  using std::to_string;
  for(unsigned int idir = 0; idir < DIM; idir++)
  {
    string lostr = "domain_bc_lo_" + to_string(idir);
    string histr = "domain_bc_hi_" + to_string(idir);
    pp.get(lostr.c_str(), loDomBC[idir]);
    pp.get(histr.c_str(), hiDomBC[idir]);
  }
  string ebbc("NoSlipWall");
  EBIBC retval(veloIC, scalIC, loDomBC, hiDomBC, ebbc);
  return retval;
}


//=================================================
void initializeVelocity(EBLevelFluxData<1>               &  a_velo,
                        const Chombo4::DisjointBoxLayout &  a_grids,
                        const Real                       &  a_dx)
{
  string interior_velo;
  ParmParse pp;
  pp.get("interior_velo", interior_velo);
  
  if(interior_velo == string("zero"))
  {
    a_velo.setVal(0.);
    Real inflow_value;
    pp.get("velocity_inflow_value", inflow_value);
    Chombo4::DataIterator dit = a_grids.dataIterator();
    for(int ibox = 0; ibox < dit.size(); ibox++)
    {
      int idir = 0; Side::LoHiSide sit = Side::Lo;  //only place I am setting a velocity
      Bx valbx = ProtoCh::getProtoBox(a_grids[dit[ibox]]);
      //I have to think of better names for functions.
      EBMACProjector::setFaceStuff(idir, sit, a_velo[dit[ibox]], valbx, inflow_value);
    }
    
    
  }
  else
  {
    Chombo4::MayDay::Error("initializeVelocity:unknown interior velocity");
  }
      
}


void defineGeometry(shared_ptr<GeometryService<MAX_ORDER> >&  a_geoserv,
                    const std::vector<Chombo4::DisjointBoxLayout>& a_grids,
                    const Chombo4::Box     & a_finestDomain,
                    const Real             & a_dx)
{
  Chombo4::pout() << "defining geometry" << endl;

  int geomGhost = 6;
  RealVect origin = RealVect::Zero();

  Chombo4::pout() << "creating implicit function" << endl;
  shared_ptr<BaseIF>  impfunc = hoeb::getImplicitFunction();

  Chombo4::pout() << "creating geometry service" << endl;
  Chombo4::Box domain = a_finestDomain;
  a_geoserv  = shared_ptr<GeometryService<MAX_ORDER> >
    (new GeometryService<MAX_ORDER>(impfunc, origin, a_dx, domain, a_grids, geomGhost));
}

int
runProjection(int a_argc, char* a_argv[])
{
  int nx      = 32;
  ParmParse pp;
  int maxGrid = 32;
  pp.get("maxGrid", maxGrid); 

  pp.get("nx"           , nx);

  Chombo4::pout() << "nx       = " << nx     << endl;
  Chombo4::pout() << "maxGrid  = " << maxGrid  << endl;

  Real dx = 1.0/nx;
  std::vector<Chombo4::DisjointBoxLayout> vecgrids;

  Chombo4::pout() << "making grids" << endl;
  IntVect domLo = IntVect::Zero;
  IntVect domHi  = (nx - 1)*IntVect::Unit;

  Chombo4::ProblemDomain domain(domLo, domHi);
  GeometryService<MAX_ORDER>::generateGrids(vecgrids, domain.domainBox(), maxGrid);

  Chombo4::pout() << "defining geometry" << endl;
  shared_ptr<GeometryService<MAX_ORDER> >  geoserv;

  Chombo4::Box domainb = domain.domainBox();
  defineGeometry(geoserv, vecgrids, domainb, dx);

  IntVect dataGhostIV =   6*IntVect::Unit;
  Point   dataGhostPt = ProtoCh::getPoint(dataGhostIV); 

  Chombo4::pout() << "making dictionary" << endl;
  vector<Chombo4::Box>    vecdomain(vecgrids.size(), domain.domainBox());
  vector<Real>   vecdx    (vecgrids.size(), dx);
  for(int ilev = 1; ilev < vecgrids.size(); ilev++)
  {
    vecdomain[ilev] = coarsen(vecdomain[ilev-1], 2);
    vecdx    [ilev] =           2*vecdx[ilev-1];
  }

  shared_ptr<EBEncyclopedia<MAX_ORDER, Real> > 
    brit(new EBEncyclopedia<MAX_ORDER, Real>(geoserv, vecgrids, vecdomain, vecdx, dataGhostPt));


  Chombo4::pout() << "inititializing data"   << endl;
  
  auto graphs = geoserv->getGraphs(domain.domainBox());
  Chombo4::DisjointBoxLayout& grids = vecgrids[0];
  shared_ptr< EBLevelFluxData<1>  > velo( new EBLevelFluxData<1>(grids, dataGhostIV, graphs));

  initializeVelocity(*velo, grids, dx);

  EBIBC bc = getIBCs();
  bool printStuff = true;
  pp.get("print_stuff", printStuff);
  hoeb::Hoeb_MAC_Projector<MAX_ORDER> proj(brit, geoserv, grids, domain.domainBox(), dx, bc, dataGhostIV, printStuff);

//begin debug
  shared_ptr< EBLevelBoxData<CELL, 1>  > testPhi( new EBLevelBoxData<CELL, 1>(grids, dataGhostIV, graphs));
  shared_ptr< EBLevelBoxData<CELL, 1>  > testLph( new EBLevelBoxData<CELL, 1>(grids, dataGhostIV, graphs));
  testPhi->setVal(1.);
  testLph->setVal(0.);
  proj.m_solver->applyOp(*testLph, *testPhi);
  testPhi->writeToFileHDF5(string("testPhi1.hdf5"), 0.0);
  testLph->writeToFileHDF5(string("testLph1.hdf5"), 0.0);
  
  proj.m_solver->applyOpPointwiseHomogeneous(*testLph, *testPhi, printStuff);
  testPhi->writeToFileHDF5(string("testPhi2.hdf5"), 0.0);
  testLph->writeToFileHDF5(string("testLph2.hdf5"), 0.0);
  Chombo4::pout() << "returning before projeciton" << endl;
  return 0;
//end debug
  
  proj.project(velo, printStuff);

  auto divergence = proj.kappaDivergence(velo, printStuff);
  
  Real divnorm = divergence->maxNorm(0);
  
  Chombo4::pout() << "max norm of post projection divergence(vel) = " << divnorm << endl;
   
  return 0;
}
//void stupidPetscTest()
//{
//  Mat A;
//#ifdef CH_MPI
//    MPI_Comm wcomm = Chombo_MPI::comm;
//#else
//    MPI_Comm wcomm = PETSC_COMM_SELF;
//#endif
////    int nc = 1;
//    PetscInt ierr;
//    ierr = MatCreate(wcomm,&A);CHKERRQ(ierr);
// 
//}

int main(int a_argc, char* a_argv[])
{
  //because of some kind of solipsistic madness, PetscInitialize calls MPI_INIT
   PetscInt ierr = PetscInitialize(&a_argc, &a_argv, "./.petscrc",PETSC_NULL); CHKERRQ(ierr);

//   stupidPetscTest();
//   return 0;
  //needs to be called after MPI_Init
  CH_TIMER_SETFILE("ebadvect.time.table");
  {
    if (a_argc < 2)
    {
      cerr<< " usage " << a_argv[0] << " <input_file_name> " << endl;
      exit(0);
    }
    char* in_file = a_argv[1];
    ParmParse  pp(a_argc-2,a_argv+2,NULL,in_file);
    runProjection(a_argc, a_argv);
  }

  Chombo4::pout() << "printing time table " << endl;
  CH_TIMER_REPORT();

  Chombo4::pout() << "about to call petsc Finalize" << std::endl;
  PetscFinalize();

  return 0;
}


