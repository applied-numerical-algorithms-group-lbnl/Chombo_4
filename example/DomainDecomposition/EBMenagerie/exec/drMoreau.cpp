#include <cmath>
#include <cstdio>
#include <iostream>


#include "EBProto.H"
#include "Chombo_EBLevelBoxData.H"

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
#include "Chombo_SmoothIntersection.H"
#include "Chombo_SmoothUnion.H"
#include <iomanip>

using Proto::BASISREALV; using Proto::SimpleSphereIF; 
using Proto::PlaneIF; using Proto::SimpleCylinderIF;
typedef CH4_Data_Choreography::DistributedData<EBGraph>  graph_distrib_t;
#include "Chombo_NamespaceHeader.H"

#define MAX_ORDER 2
///we print out data = the volume fraction field
void  
fillKappa(EBLevelBoxData<CELL, 1>&  a_kappa,
          const shared_ptr<GeometryService<2> >   & a_geoserv,
          const shared_ptr<graph_distrib_t>   & a_graphs,
          const DisjointBoxLayout a_grids,
          const Box& a_domain)
{
  CH_TIME("fillkappa");
  DataIterator dit = a_grids.dataIterator();
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

/**************/
int
makeTwoSpheres(int a_argc, char* a_argv[])
{

  int nx      = 32;
  int maxGrid = 32;
  ParmParse pp("two_spheres");
    
  pp.get("nx"     , nx);
  pp.get("maxGrid", maxGrid);
  
  Real rad_one = 0.1;
  Real rad_two = 0.1;
  vector<Real> cen_one(DIM, 0.0);
  vector<Real> cen_two(DIM, 1.0);

  Chombo4::pout() << "two sphere case: "      << endl;
  pp.get(   "rad_one"     , rad_one);
  pp.get(   "rad_two"     , rad_two);
  pp.getarr("cen_one"     , cen_one, 0, DIM);
  pp.getarr("cen_two"     , cen_two, 0, DIM);

  Chombo4::pout() << "nx      = " << nx       << endl;
  Chombo4::pout() << "maxGrid = " << maxGrid  << endl;
  Chombo4::pout() << "rad_one = " << rad_one  << endl;
  Chombo4::pout() << "rad_two = " << rad_two  << endl;
  RealVect ABC_one, X0_one;
  RealVect ABC_two, X0_two;
  for(int idir = 0; idir < DIM; idir++)
  {
    Chombo4::pout() << "cen_one[ " << idir << "] = " << cen_one[idir]  << endl;
    Chombo4::pout() << "cen_two[ " << idir << "] = " << cen_two[idir]  << endl;

    ABC_one[idir] = 1.0; //makes it a sphere
    ABC_two[idir] = 1.0; //makes it a sphere

    X0_one[idir] = cen_one[idir];
    X0_two[idir] = cen_two[idir];
  }


  IntVect domLo = IntVect::Zero;
  IntVect domHi  = (nx - 1)*IntVect::Unit;

// EB and periodic do not mix
  ProblemDomain domain(domLo, domHi);

  std::vector<Box> boxes;
  unsigned int blockfactor = 8;
  domainSplit(domain.domainBox(), boxes, maxGrid, blockfactor);
  
  std::vector<int> procs;
  Chombo4::pout() << "making grids" << endl;
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
  shared_ptr<BaseIF>     intersect(new IntersectionIF(both_spheres));

  Chombo4::pout() << "defining geometry" << endl;
  shared_ptr<GeometryService<MAX_ORDER> >  geoserv(new GeometryService<MAX_ORDER>(intersect, origin, dx, domain.domainBox(), grids, geomGhost));
  auto graphs = geoserv->getGraphs(domain.domainBox());

  Chombo4::pout() << "making data" << endl;
  EBLevelBoxData<CELL,   1>  kappa(grids, dataGhostIV, graphs);
  fillKappa(kappa, geoserv, graphs, grids, domain.domainBox());
  Real coveredVal = -1;
  kappa.writeToFileHDF5(string("two_spheres_kappa.hdf5"), coveredVal);

  delete sphere_one;
  delete sphere_two;
  Chombo4::pout() << "exiting twoSpheres" << endl;
  return 0;
}

/**************/
int
makeIntersectingSpheres(int a_argc, char* a_argv[])
{

  int nx      = 32;
  int maxGrid = 32;
  ParmParse pp("intersecting_spheres");
    
  pp.get("nx"     , nx);
  pp.get("maxGrid", maxGrid);

  //put one sphere in each quadrant or octant
#if DIM==2
  const unsigned int num_spheres = 5;
#else
  const unsigned int num_spheres = 9;
#endif  
  Real rad = 0.1;
  Real offsetmag = 1.1*rad; 
  pp.get("radius", rad);
  pp.get("offsetmag", offsetmag);
  Chombo4::pout() << "nx      = " << nx       << endl;
  Real dx = 1.0/nx;
  
  Chombo4::pout() << "dx      = " << dx       << endl;
  Chombo4::pout() << "maxGrid = " << maxGrid  << endl;
  Chombo4::pout() << "radius = " << rad << endl;
  Chombo4::pout() << "offset = " << offsetmag << endl;
  vector<RealVect> centers(num_spheres);
  RealVect centercenter = RealVect::Unit();
  centercenter *= 0.5;
  RealVect offsetx = offsetmag*BASISREALV(0);
  RealVect offsety = offsetmag*BASISREALV(1);
#if DIM==2
  centers[0] = centercenter - offsetx - offsety;
  centers[1] = centercenter + offsetx - offsety;
  centers[2] = centercenter - offsetx + offsety;
  centers[3] = centercenter + offsetx + offsety;
  centers[4] = centercenter;
#else
  RealVect offsetz = offsetmag*BASISREALV(2);
  centers[0] = centercenter - offsetx - offsety - offsetz;
  centers[1] = centercenter + offsetx - offsety - offsetz;
  centers[2] = centercenter - offsetx + offsety - offsetz;
  centers[3] = centercenter + offsetx + offsety - offsetz;
  centers[4] = centercenter - offsetx - offsety + offsetz;
  centers[5] = centercenter + offsetx - offsety + offsetz;
  centers[6] = centercenter - offsetx + offsety + offsetz;
  centers[7] = centercenter + offsetx + offsety + offsetz;
  centers[8] = centercenter;
#endif  

  Chombo4::pout() << "intersecting sphere case with rad = : " << rad   << endl;
  for(unsigned int icen = 0; icen < num_spheres; icen++)
  {
    Chombo4::pout() << "center[" << icen << "] = " << centers[icen] << endl;
  }

  IntVect domLo = IntVect::Zero;
  IntVect domHi  = (nx - 1)*IntVect::Unit;

  // EB and periodic do not mix
  ProblemDomain domain(domLo, domHi);

  std::vector<Box> boxes;
  unsigned int blockfactor = 8;
  domainSplit(domain.domainBox(), boxes, maxGrid, blockfactor);
  
  std::vector<int> procs;
  Chombo4::pout() << "making grids" << endl;
  LoadBalance(procs, boxes);
  DisjointBoxLayout grids(boxes, procs, domain);
  grids.printBalance();

  IntVect dataGhostIV =   IntVect::Zero;
  int geomGhost = 4;
  RealVect origin = RealVect::Zero();

  std::vector<BaseIF*> spheres(num_spheres);
  bool inside = false;
  for(unsigned int icen = 0; icen < num_spheres; icen++)
  {
    SimpleSphereIF* sphereptr = new SimpleSphereIF(centers[icen], rad, inside);
    spheres[icen] = static_cast<BaseIF*>(sphereptr);
  }
  shared_ptr<BaseIF>     intersect(new IntersectionIF(spheres));
  Chombo4::pout() << "defining geometry" << endl;
  shared_ptr<GeometryService<MAX_ORDER> >
    geoserv(new GeometryService<MAX_ORDER>(intersect, origin, dx, domain.
                                           domainBox(), grids, geomGhost));
  
  auto graphs = geoserv->getGraphs(domain.domainBox());

  Chombo4::pout() << "making data" << endl;
  EBLevelBoxData<CELL,   1>  kappa(grids, dataGhostIV, graphs);
  fillKappa(kappa, geoserv, graphs, grids, domain.domainBox());
  Real coveredVal = -1;
  kappa.writeToFileHDF5(string("intersecting_spheres_kappa.hdf5"), coveredVal);

  for(int icen = 0; icen < num_spheres; icen++)
  {
    delete spheres[icen];
    spheres[icen] = NULL;
  }
  Chombo4::pout() << "exiting intersectingSpheres" << endl;
  return 0;
}
/**************/
int
makeSmoothedIntersectingSpheres(int a_argc, char* a_argv[])
{
  int nx      = 32;
  int maxGrid = 32;
  ParmParse pp("smoothed_intersecting_spheres");
    
  pp.get("nx"     , nx);
  pp.get("maxGrid", maxGrid);

  //put one sphere in each quadrant or octant
#if DIM==2
  const unsigned int num_spheres = 5;
#else
  const unsigned int num_spheres = 9;
#endif  
  Real rad = 0.1;
  Real offsetmag = 1.1*rad; 
  Real dx = 1.0/nx;
  Real deltamag = 0.1;
  Real delta = deltamag*dx;
  pp.get("radius", rad);
  pp.get("offsetmag", offsetmag);
  pp.get("deltamag" , deltamag);
  Chombo4::pout() << "nx      = " << nx       << endl;
  Chombo4::pout() << "dx      = " << dx       << endl;
  Chombo4::pout() << "maxGrid = " << maxGrid  << endl;
  Chombo4::pout() << "radius  = " << rad << endl;
  Chombo4::pout() << "offset  = " << offsetmag << endl;
  Chombo4::pout() << "delta   = " << deltamag << "*dx = " << delta << endl;
  
  vector<RealVect> centers(num_spheres);
  RealVect centercenter = RealVect::Unit();
  centercenter *= 0.5;
  RealVect offsetx = offsetmag*BASISREALV(0);
  RealVect offsety = offsetmag*BASISREALV(1);
  
#if DIM==2
  centers[0] = centercenter - offsetx - offsety;
  centers[1] = centercenter + offsetx - offsety;
  centers[2] = centercenter - offsetx + offsety;
  centers[3] = centercenter + offsetx + offsety;
  centers[4] = centercenter;
#else
  RealVect offsetz = offsetmag*BASISREALV(2);
  centers[0] = centercenter - offsetx - offsety - offsetz;
  centers[1] = centercenter + offsetx - offsety - offsetz;
  centers[2] = centercenter - offsetx + offsety - offsetz;
  centers[3] = centercenter + offsetx + offsety - offsetz;
  centers[4] = centercenter - offsetx - offsety + offsetz;
  centers[5] = centercenter + offsetx - offsety + offsetz;
  centers[6] = centercenter - offsetx + offsety + offsetz;
  centers[7] = centercenter + offsetx + offsety + offsetz;
  centers[8] = centercenter ;
#endif  

  Chombo4::pout() << "smoothed intersecting sphere case with rad = : " << rad     << endl;
  for(unsigned int icen = 0; icen < num_spheres; icen++)
  {
    Chombo4::pout() << "center[" << icen << "] = " << centers[icen] << endl;
  }

  IntVect domLo = IntVect::Zero;
  IntVect domHi  = (nx - 1)*IntVect::Unit;

  // EB and periodic do not mix
  ProblemDomain domain(domLo, domHi);

  std::vector<Box> boxes;
  unsigned int blockfactor = 8;
  domainSplit(domain.domainBox(), boxes, maxGrid, blockfactor);
  
  std::vector<int> procs;
  Chombo4::pout() << "making grids" << endl;
  LoadBalance(procs, boxes);
  DisjointBoxLayout grids(boxes, procs, domain);
  grids.printBalance();

  IntVect dataGhostIV =   IntVect::Zero;
  int geomGhost = 4;
  RealVect origin = RealVect::Zero();

  std::vector<BaseIF*> spheres(num_spheres);
  bool inside = false;
  for(unsigned int icen = 0; icen < num_spheres; icen++)
  {
    SimpleSphereIF* sphereptr = new SimpleSphereIF(centers[icen], rad, inside);
    spheres[icen] = static_cast<BaseIF*>(sphereptr);
  }
  shared_ptr<BaseIF>     intersect(new SmoothIntersection(spheres, dx));
  Chombo4::pout() << "defining geometry" << endl;
  shared_ptr<GeometryService<MAX_ORDER> >
    geoserv(new GeometryService<MAX_ORDER>(intersect, origin, dx, domain.
                                           domainBox(), grids, geomGhost));
  
  auto graphs = geoserv->getGraphs(domain.domainBox());

  Chombo4::pout() << "making data" << endl;
  EBLevelBoxData<CELL,   1>  kappa(grids, dataGhostIV, graphs);
  fillKappa(kappa, geoserv, graphs, grids, domain.domainBox());
  Real coveredVal = -1;
  kappa.writeToFileHDF5(string("smoothed_intersecting_spheres_kappa.hdf5"), coveredVal);
  
  for(int icen = 0; icen < num_spheres; icen++)
  {
    delete spheres[icen];
    spheres[icen] = NULL;
  }

  Chombo4::pout() << "exiting smoothed intersectingSpheres" << endl;
  return 0;
}
/**************/
int
makeMollifiedBox(int a_argc, char* a_argv[])
{
  int nx      = 32;
  int maxGrid = 32;
  ParmParse pp("mollified_box");
    
  pp.get("nx"     , nx);
  pp.get("maxGrid", maxGrid);
  
  Real width = 0.1;
  vector<Real> center(DIM, 0.0);

  Chombo4::pout() << "mollified_box:"      << endl;
  pp.getarr("center", center, 0, DIM);
  Chombo4::pout() << "nx      = " << nx       << endl;
  Chombo4::pout() << "maxGrid = " << maxGrid  << endl;
  Chombo4::pout() << "width   = " << nx       << endl;
  for(int idir = 0; idir < DIM; idir++)
  {
    Chombo4::pout() << "center[ " << idir << "] = " << center[idir]  << endl;
  }


  IntVect domLo = IntVect::Zero;
  IntVect domHi  = (nx - 1)*IntVect::Unit;

// EB and periodic do not mix
  ProblemDomain domain(domLo, domHi);

  std::vector<Box> boxes;
  unsigned int blockfactor = 8;
  domainSplit(domain.domainBox(), boxes, maxGrid, blockfactor);
  
  std::vector<int> procs;
  Chombo4::pout() << "making grids" << endl;
  LoadBalance(procs, boxes);
  DisjointBoxLayout grids(boxes, procs, domain);
  grids.printBalance();

  IntVect dataGhostIV =   IntVect::Zero;
  int geomGhost = 4;
  RealVect origin = RealVect::Zero();
  Real dx = 1.0/nx;
  std::vector<BaseIF*> planes(2*DIM);
  RealVect centerrv;
  for(int idir = 0; idir < DIM; idir++)
  {
    centerrv[idir]= center[idir];
  }
  for(int idir = 0; idir < DIM; idir++)
  {
    RealVect normal = RealVect::Zero();
    RealVect point  = centerrv;

    normal[idir] = -1;
    point[ idir] -= width/2.0;
    planes[2*idir  ] = new PlaneIF(point, normal);
    normal[idir] = 1;
    point[ idir] += width;
    planes[2*idir+1] = new PlaneIF(point, normal);
  }

  std::shared_ptr<BaseIF> box_planes(new UnionIF(planes));
  Chombo4::pout() << "defining geometry" << endl;
  shared_ptr<GeometryService<MAX_ORDER> >
    geoserv(new GeometryService<MAX_ORDER>(box_planes, origin, dx, domain.domainBox(), grids, geomGhost));
  
  auto graphs = geoserv->getGraphs(domain.domainBox());

  Chombo4::pout() << "making data" << endl;
  EBLevelBoxData<CELL,   1>  kappa(grids, dataGhostIV, graphs);
  fillKappa(kappa, geoserv, graphs, grids, domain.domainBox());
  Real coveredVal = -1;
  kappa.writeToFileHDF5(string("mollified_box.hdf5"), coveredVal);

  //clean up memory
  for(int iplane = 0; iplane < 2*DIM; iplane++)
  {
    delete planes[iplane];
  }
  Chombo4::pout() << "exiting mollified box" << endl;
  return 0;
}

/**************/
int
makeCylinderedSphere(int a_argc, char* a_argv[])
{
  int nx      = 32;
  int maxGrid = 32;
  ParmParse pp("cylindered_sphere");
    
  pp.get("nx"     , nx);
  pp.get("maxGrid", maxGrid);
  
  Real sphereRadius = 0.1;
  Real cylinderRadius = 0.4;

  Chombo4::pout() << "cylindered sphere (sphere/cyl center is domain center):"      << endl;
  pp.get("sphere_radius"  , sphereRadius);
  pp.get("cylinder_radius", cylinderRadius);
  Chombo4::pout() << "nx              = " << nx       << endl;
  Chombo4::pout() << "maxGrid         = " << maxGrid  << endl;
  Chombo4::pout() << "sphere radius   = " << sphereRadius       << endl;
  Chombo4::pout() << "cylinder radius = " << sphereRadius       << endl;
  


  IntVect domLo = IntVect::Zero;
  IntVect domHi  = (nx - 1)*IntVect::Unit;

// EB and periodic do not mix
  ProblemDomain domain(domLo, domHi);

  std::vector<Box> boxes;
  unsigned int blockfactor = 8;
  domainSplit(domain.domainBox(), boxes, maxGrid, blockfactor);
  
  std::vector<int> procs;
  Chombo4::pout() << "making grids" << endl;
  LoadBalance(procs, boxes);
  DisjointBoxLayout grids(boxes, procs, domain);
  grids.printBalance();

  IntVect dataGhostIV =   IntVect::Zero;
  int geomGhost = 4;
  RealVect origin = RealVect::Zero();
  Real dx = 1.0/nx;
  std::vector<BaseIF*> planes(2*DIM);
  RealVect center = RealVect::Unit();
  center *= 0.5;
  
  SimpleCylinderIF cylinder(center, cylinderRadius, true);
  SimpleSphereIF     sphere(center,   sphereRadius, false);
  vector<BaseIF*> bothIF(2, NULL);
  bothIF[0] = static_cast<BaseIF*>(&cylinder);
  bothIF[1] = static_cast<BaseIF*>(&sphere);

  std::shared_ptr<BaseIF> magilla(new IntersectionIF(bothIF));
  Chombo4::pout() << "defining geometry" << endl;
  shared_ptr<GeometryService<MAX_ORDER> >
    geoserv(new GeometryService<MAX_ORDER>(magilla, origin, dx, domain.domainBox(), grids, geomGhost));
  
  auto graphs = geoserv->getGraphs(domain.domainBox());

  Chombo4::pout() << "making data" << endl;
  EBLevelBoxData<CELL,   1>  kappa(grids, dataGhostIV, graphs);
  fillKappa(kappa, geoserv, graphs, grids, domain.domainBox());
  Real coveredVal = -1;
  kappa.writeToFileHDF5(string("cylindered_sphere.hdf5"), coveredVal);

  //clean up memory
  for(int iplane = 0; iplane < 2*DIM; iplane++)
  {
    delete planes[iplane];
  }
  Chombo4::pout() << "exiting mollified box" << endl;
  return 0;
}

#include "Chombo_NamespaceFooter.H"

int main(int a_argc, char* a_argv[])
{
#ifdef CH_MPI
  MPI_Init(&a_argc, &a_argv);
  Chombo4::pout() << "MPI INIT called" << std::endl;
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

    Chombo4::pout() << "making two spheres" << endl;
    makeTwoSpheres(a_argc, a_argv);
    
    Chombo4::pout() << "making the dreaded mollified box" << endl;
    Chombo4::makeMollifiedBox(a_argc, a_argv);

    Chombo4::pout() << "making intersecting spheres (unsmoothed)" << endl;
    Chombo4::makeIntersectingSpheres(a_argc, a_argv);

    Chombo4::pout() << "making intersecting spheres (smoothed)" << endl;
    Chombo4::makeSmoothedIntersectingSpheres(a_argc, a_argv);
    
#if DIM==3    
    Chombo4::pout() << "making cylindered sphere" << endl;
    Chombo4::makeCylinderedSphere(a_argc, a_argv);
#endif    
  }

  Chombo4::pout() << "printing time table " << endl;
  CH_TIMER_REPORT();
#ifdef CH_MPI
  Chombo4::pout() << "about to call MPI Finalize" << std::endl;
  MPI_Finalize();
#endif
  return 0;
}
