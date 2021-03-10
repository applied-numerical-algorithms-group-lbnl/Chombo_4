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
#include "Chombo_GeometryService.H"
#include "Chombo_EBDictionary.H"
#include "Chombo_EBChombo.H"
#include "Proto_SimpleImplicitFunctions.H"
#include "Chombo_UnionIntersection.H"
#include <iomanip>

#include "Chombo_NamespaceHeader.H"

#define MAX_ORDER 2

using std::cout;
using std::endl;
using std::shared_ptr;


void
dumpPPS(const PointSet* a_ivs)
{
  for(PointSetIterator ivsit(*a_ivs);  ivsit.ok(); ++ivsit)
  {
    std::cout << ivsit() << " " ;
  }
  std::cout << std::endl;
 }

void  
fillKappa(EBLevelBoxData<CELL, 1>&  a_kappa,
          const shared_ptr<GeometryService<2> >   & a_geoserv,
          const shared_ptr<LevelData<EBGraph> >   & a_graphs,
          const DisjointBoxLayout a_grids,
          const Box& a_domain)
{
  CH_TIME("fillkappa");
  DataIterator dit = a_grids.dataIterator();
  int ideb = 0;
  for(int ibox = 0; ibox < dit.size(); ++ibox)
  {
    auto& kappdat = a_kappa[dit[ibox]];
    auto grid =a_grids[dit[ibox]];
    Bx  grbx = kappdat.inputBox();
    const EBGraph  & graph = (*a_graphs)[dit[ibox]];
    EBHostData<CELL, Real, 1> hostdat(grbx, graph);
    //fill kappa on the host then copy to the device
    a_geoserv->fillKappa(hostdat, grid, dit[ibox], a_domain);
    // now copy to the device
    EBLevelBoxData<CELL, 1>::copyToDevice(kappdat, hostdat);
  }
}
int
makeTwoSpheres(int a_argc, char* a_argv[])
{
  Real coveredval = -1;
  int nx      = 32;
  int maxGrid = 32;
  ParmParse pp("two_spheres");
    
  pp.get("nx"     , nx);
  pp.get("maxGrid", maxGrid);
  
  Real rad_one = 0.1;
  Real rad_two = 0.1;
  vector<Real> cen_one(DIM, 0.0);
  vector<Real> cen_two(DIM, 1.0);

  pp.get(   "rad_one"     , rad_one);
  pp.get(   "rad_two"     , rad_two);
  pp.getarr("cen_one"     , cen_one, 0, DIM);
  pp.getarr("cen_two"     , cen_two, 0, DIM);

  pout() << "two sphere case= "      << endl;
  pout() << "nx      = " << nx       << endl;
  pout() << "maxGrid = " << maxGrid  << endl;
  pout() << "rad_one = " << rad_one  << endl;
  pout() << "rad_two = " << rad_two  << endl;
  RealVect ABC_one, X0_one;
  RealVect ABC_two, X0_two;
  for(int idir = 0; idir < DIM; idir++)
  {
    pout() << "cen_one[ " << idir << "] = " << cen_one[idir]  << endl;
    pout() << "cen_two[ " << idir << "] = " << cen_two[idir]  << endl;

    ABC_one[idir] = 1.0; //makes it a sphere
    ABC_two[idir] = 1.0; //makes it a sphere

    X0_one[idir] = cen_one[idir];
    X0_two[idir] = cen_two[idir];
  }


  IntVect domLo = IntVect::Zero;
  IntVect domHi  = (nx - 1)*IntVect::Unit;

// EB and periodic do not mix
  ProblemDomain domain(domLo, domHi);

  Vector<Box> boxes;
  unsigned int blockfactor = 8;
  domainSplit(domain, boxes, maxGrid, blockfactor);
  
  Vector<int> procs;
  pout() << "making grids" << endl;
  LoadBalance(procs, boxes);
  DisjointBoxLayout grids(boxes, procs, domain);
  grids.printBalance();

  IntVect dataGhostIV =   IntVect::Zero;
  int geomGhost = 4;
  RealVect origin = RealVect::Zero();
  Real dx = 1.0/nx;
  SimpleEllipsoidIF* sphere_one = new Proto::SimpleEllipsoidIF(ABC_one, X0_one, rad_one, false);
  SimpleEllipsoidIF* sphere_two = new Proto::SimpleEllipsoidIF(ABC_two, X0_two, rad_two, false);
  std::vector<BaseIF*> both_spheres(2);
  both_spheres[0] = static_cast<BaseIF*>(sphere_one);
  both_spheres[1] = static_cast<BaseIF*>(sphere_two);
  shared_ptr<BaseIF>  impfunc(new UnionIF(both_spheres));
//  Bx domainpr = getProtoBox(domain.domainBox());
  pout() << "defining geometry" << endl;
  shared_ptr<GeometryService<MAX_ORDER> >  geoserv(new GeometryService<MAX_ORDER>(impfunc, origin, dx, domain.domainBox(), grids, geomGhost));
  shared_ptr<LevelData<EBGraph> > graphs = geoserv->getGraphs(domain.domainBox());

  pout() << "making data" << endl;
  EBLevelBoxData<CELL,   1>  kappa(grids, dataGhostIV, graphs);
  fillKappa(kappa, geoserv, graphs, grids, domain.domainBox());
  Real coveredVal = -1;
  kappa.writeToFileHDF5(string("two_spheres_kappa.hdf5"), coveredVal);
  
  pout() << "exiting twoSpheres" << endl;
  return 0;
}

#include "Chombo_NamespaceFooter.H"

int main(int a_argc, char* a_argv[])
{
#ifdef CH_MPI
  MPI_Init(&a_argc, &a_argv);
  pout() << "MPI INIT called" << std::endl;
#endif

//   Dr. Moreau : What is the law?
//   Sayer of the Law : Not to eat meat, that is the law. Are we not men?
//   Beasts (in unison) : Are we not men?
//   
//   Dr. Moreau : What is the law?
//   Sayer of the Law : Not to go on all fours, that is the law. Are we not men?
//   Beasts (in unison) : Are we not men?
//   
//   Dr. Moreau : What is the law?
//   Sayer of the Law : Not to spill blood, that is the law. Are we not men?
//   Beasts (in unison) : Are we not men?


  //needs to be called after MPI_Init
  CH_TIMER_SETFILE("dr.moreau.time.table");
  {
    if (a_argc < 2)
    {
      cerr<< " usage " << a_argv[0] << " <input_file_name> " << endl;
      exit(0);
    }
    char* in_file = a_argv[1];
    ParmParse  pp(a_argc-2,a_argv+2,NULL,in_file);

    pout() << "making two spheres" << endl;
    makeTwoSpheres(a_argc, a_argv);
//    pout() << "making the dreaded mollified box" << endl;
//    Chombo4::makeMollifiedBox(a_argc, a_argv);
//    pout() << "making packed spheres" << endl;
//    Chombo4::make(a_argc, a_argv);
  }

  pout() << "printing time table " << endl;
  CH_TIMER_REPORT();
#ifdef CH_MPI
  pout() << "about to call MPI Finalize" << std::endl;
  MPI_Finalize();
#endif
  return 0;
}
