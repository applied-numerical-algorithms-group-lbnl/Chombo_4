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
#include "Chombo_EBLevelBoxData.H"
#include "Chombo_EBLevelFluxData.H"
#include "EBMultigrid.H"
#include "EBMultigridFunctions.H"

#include <iomanip>

#define MAX_ORDER 2

#define PI 3.141592653589793

using std::cout;
using std::endl;
using std::shared_ptr;

//Function to fill out, cell by cell, an EB data box
PROTO_KERNEL_START
unsigned int  setRhsPtF(int              a_pt[DIM],
                        Var<Real, 1>     a_phi,
                        Real             a_dx)
{
  Real x = (Real(a_pt[0]) + 0.5)*(a_dx);
  Real y = (Real(a_pt[1]) + 0.5)*(a_dx);
  Real val = -2.0*PI*PI*sin(PI*x)*sin(PI*y);
  a_phi(0) = val;
  return 0;
}
PROTO_KERNEL_END(setRhsPtF, setRhsPt)

//Function to fill out, cell by cell, an EB data box
PROTO_KERNEL_START
unsigned int  setPhiPtF(int              a_pt[DIM],
                        Var<Real, 1>     a_phi,
                        Real             a_dx)
{
  Real x = (Real(a_pt[0]) + 0.5)*(a_dx);
  Real y = (Real(a_pt[1]) + 0.5)*(a_dx);
  Real val = sin(PI*x)*sin(PI*y);
  a_phi(0) = val;
  return 0;
}
PROTO_KERNEL_END(setPhiPtF, setPhiPt)

//=================================================
void
fillKappa(shared_ptr<GeometryService<MAX_ORDER> >& a_geoserv,
          Chombo4::DisjointBoxLayout               m_grids,
          shared_ptr<LevelData<EBGraph> >&         m_graphs,
          Chombo4::Box&                            m_domain,
          Chombo4::Copier&                         m_exchangeCopier,
          EBLevelBoxData<CELL, 1>& m_kappa)
{
  CH_TIME("EBAdvection::fillKappa");
  Chombo4::DataIterator dit = m_grids.dataIterator();
  for(int ibox = 0; ibox < dit.size(); ++ibox)
  {
    Chombo4::Box grid =m_grids[dit[ibox]];
    const EBGraph  & graph = (*m_graphs)[dit[ibox]];
    EBBoxData <CELL, Real, 1> & devidat = m_kappa[dit[ibox]];
    EBHostData<CELL, Real, 1>   hostdat(devidat.box(), graph);
    //fill kappa on the host then copy to the device
    a_geoserv->fillKappa(hostdat, grid, dit[ibox], m_domain);
    // now copy to the device
    EBLevelBoxData<CELL, 1>::copyToDevice(m_kappa[dit[ibox]], hostdat);
  }
  m_kappa.exchange(m_exchangeCopier);
}

//==================================================
int
runConvergence(int a_argc, char* a_argv[])
{

  //Default values definition
  //=========================

  //Value for covered cells
  Real coveredval = -1;
  //Total number of cells in each direction
  int nx      = 32;
  //Max number of cells per block in each direction
  int maxGrid = 32;
  //Ellipse Center values
  Real x0 = 0.5;
  Real y0 = 0.5;
  Real z0 = 0.5;
  //Ellipse stretching values
  Real A = 1.0;
  Real B = 1.0;
  Real C = 1.0;
  //Ellipse radius
  Real R = 0.25;
  //Number of streams for CUDA-GPU
  int nStream    = 8;
  //Domain BC - 0 Neumann - 1 Dirichelt
  int dombc = 1;
  //EB BC - 0 Neumann - 1 Dirichelt
  int ebbc  = 1;
  //Solver parameters
  Real alpha = 0.0;
  Real beta = 1.0;
  Real tol = 1e-11;
  int  maxIter = 300;
  int numSmooth  = 2;

  //Reading parameters from the inputs file
  ParmParse pp;

  pp.get("domainBC"  , dombc);
  pp.get("EBBC"      , ebbc);
  pp.get("nx"     , nx);
  pp.get("nstream", nStream);
  pp.get("maxGrid", maxGrid);
  pp.get("x0"     , x0);
  pp.get("y0"     , y0);
  pp.get("z0"     , z0);
  pp.get("A"      , A);
  pp.get("B"      , B);
  pp.get("C"      , C);
  pp.get("R"      , R);
  pp.get("coveredval", coveredval);
  pp.get("alpha"     , alpha);
  pp.get("beta"      , beta);
  pp.get("tolerance" , tol);
  pp.get("maxIter"   , maxIter);
  pp.get("numSmooth" , numSmooth);

  pout() << "nx      = " << nx       << endl;
  pout() << "maxGrid = " << maxGrid  << endl;
  pout() << "x0      = " << x0       << endl;
  pout() << "y0      = " << y0       << endl;
  pout() << "z0      = " << z0       << endl;
  pout() << "A       = " << A        << endl;
  pout() << "B       = " << B        << endl;
  pout() << "C       = " << C        << endl;
  pout() << "R       = " << R        << endl;

  pout() << "nstream = " << nStream  << endl;

  //Setting CUDA streams
#ifdef PROTO_CUDA
  Proto::DisjointBoxLayout::setNumStreams(nStream);
#endif

  //Setting ellipse parameters
  RealVect ABC, X0;
  ABC[0] = A;
  ABC[1] = B;
  X0[0] = x0;
  X0[1] = y0;
#if DIM==3
  ABC[2] = C;
  X0[2] = z0;
#endif

  //Setting Lo and Hi domain vectors
  IntVect domLo = IntVect::Zero;
  IntVect domHi  = (nx - 1)*IntVect::Unit;

  //Setting problem domain
  //EB and periodic do not mix
  Chombo4::ProblemDomain domain(domLo, domHi);

  //Declaring box domain -- needed for many function calls
  Chombo4::Box domain_box = domain.domainBox();

  //Splitting the domain into boxes and doing load balance
  Vector<Chombo4::Box> boxes;
  unsigned int blockfactor = 8;
  domainSplit(domain, boxes, maxGrid, blockfactor);

  //Generating a grids vector for different levels -- needed for EBMultigrid
  Vector<Chombo4::DisjointBoxLayout> vecgrids;
  pout() << "making grids" << endl;
  GeometryService<2>::generateGrids(vecgrids, domain.domainBox(), maxGrid);
  Chombo4::DisjointBoxLayout grids = vecgrids[0];
  grids.printBalance();

  //Ghost cells definition for communications
  IntVect dataGhostIV =   4*IntVect::Unit;
  Point   dataGhostPt = ProtoCh::getPoint(dataGhostIV);
  int geomGhost = 4;
  RealVect origin = RealVect::Zero();
  //Grid size definition
  Real dx = 1.0/nx;
  //Implicit function for the ellipse -- false argument means the calculation happens outside the ellipse
  shared_ptr<BaseIF> impfunc(new Proto::SimpleEllipsoidIF(ABC, X0, R, false));

  pout() << "defining geometry" << endl;
  //Geometry definition 
  GeometryService<2>* geomptr = new GeometryService<2>(impfunc, origin, dx, domain.domainBox(), vecgrids, geomGhost);
  shared_ptr< GeometryService<2> >  geoserv(geomptr);

  pout() << "making dictionary" << endl;
  //EB dictionary creation
  //Definition of the grids and grid size for different levels -- needed for EBMultigrid
  vector<Chombo4::Box>    vecdomain(vecgrids.size(), domain.domainBox());
  vector<Real>   vecdx    (vecgrids.size(), dx);
  for(int ilev = 1; ilev < vecgrids.size(); ilev++)
  { 
    vecdomain[ilev] = coarsen(vecdomain[ilev-1], 2);
    vecdx    [ilev] =           2*vecdx[ilev-1];
  }
  shared_ptr<EBDictionary<2, Real, CELL, CELL> >
    dictionary(new EBDictionary<2, Real, CELL, CELL>(geoserv, vecgrids, vecdomain, vecdx, dataGhostPt));

  //Stencil type for EBMultigrid
  string stenname = StencilNames::Poisson2;
  //BC for domain and eb definitions
  string dombcname, ebbcname;
  if(dombc == 0)
  { 
    dombcname = StencilNames::Neumann;
    pout() << "using Neumann BCs at domain" << endl;
  }
  else
  { 
    dombcname = StencilNames::Dirichlet;
    pout() << "using Dirichlet BCs at domain" << endl;
  }

  if(ebbc == 0)
  { 
    ebbcname = StencilNames::Neumann;
    pout() << "using Neumann BCs at EB" << endl;
  }
  else
  { 
    ebbcname = StencilNames::Dirichlet;
    pout() << "using Dirichlet BCs at EB" << endl;
  }
  string alldombc[2*DIM];
  for(int iface = 0; iface < 2*DIM; iface++)
  { 
    alldombc[iface] = dombcname;
  }

  //EB graphs declaration 
  shared_ptr<LevelData<EBGraph> > graphs = geoserv->getGraphs(domain.domainBox());

  //EBMultigrid declaration
  EBMultigrid solver(dictionary, geoserv, alpha, beta, dx, grids, stenname, dombcname, ebbcname, domain.domainBox(), dataGhostIV);
  EBMultigrid::s_numSmoothUp   = numSmooth;
  EBMultigrid::s_numSmoothDown = numSmooth;

  pout() << "making data" << endl;
  //EB Data declaration
  EBLevelBoxData<CELL,   1>  phi(grids, dataGhostIV, graphs);
  EBLevelBoxData<CELL,   1>  rhs(grids, dataGhostIV, graphs);
  EBLevelBoxData<CELL,   1>  scalcell(grids, dataGhostIV, graphs);
  EBLevelBoxData<CELL,   1>  kappa(grids, dataGhostIV, graphs);

  //Ghost cells operations are defined by Copier
  Chombo4::Copier exchangeCopier;
  exchangeCopier.exchangeDefine(grids, dataGhostIV);
  //fillKappa is a function that will fill kappa with volume-fraction generated by the ellipse
  fillKappa(geoserv,grids,graphs,domain_box,exchangeCopier,kappa);

  //Data iterator definition
  Chombo4::DataIterator dit = grids.dataIterator();
  pout() << "setting values" << endl;

  //Iteration over boxes
//#pragma omp parallel for
  for(int ibox = 0; ibox < dit.size(); ibox++)
  {
    //Accessing EB box data
    EBBoxData<CELL, Real, 1>& phibd = phi[dit[ibox]];
    EBBoxData<CELL, Real, 1>& rhsbd = rhs[dit[ibox]];
    EBBoxData<CELL, Real, 1>& solbd = scalcell[dit[ibox]];
    //Setting values to the entire EB box 
    phibd.setVal(0.0);
    //Extracting the box from the EB box data
    Bx rhsbx = rhsbd.box();
    Bx solbx = solbd.box();
    size_t numflopspt = 0;
    //Function to fill out, cell by cell, an EB box data
    ebforallInPlace_i(numflopspt, "setRhsPt", setRhsPt, rhsbx, rhsbd, dx);
    ebforallInPlace_i(numflopspt, "setPhiPt", setPhiPt, solbx, solbd, dx);
  }

  //Calling solve 
  solver.solve(phi, rhs, tol, maxIter, false);

  //Writing EB box data into hdf5 files
  writeEBLevelHDF5<1>(string("phi.hdf5"), phi, kappa, domain.domainBox(), graphs, coveredval, dx, 1.0, 0.0);
  writeEBLevelHDF5<1>(string("rhs.hdf5"), rhs, kappa, domain.domainBox(), graphs, coveredval, dx, 1.0, 0.0);
  writeEBLevelHDF5<1>(string("sol.hdf5"), scalcell, kappa, domain.domainBox(), graphs, coveredval, dx, 1.0, 0.0);
  pout() << "exiting runConvergence" << endl;
  return 0;
}

//========================================================
int
main(int a_argc, char* a_argv[])
{
#ifdef CH_MPI
  MPI_Init(&a_argc,&a_argv);
  pout() << "MPI INIT called" << std::endl;
  //needs to be called after MPI_Init
  CH_TIMER_SETFILE("geomConv.time.table");
#endif
  {
    // Check for an input file
    char* inFile = NULL;

    if (a_argc > 1)
      {
        inFile = a_argv[1];
      }
    else
      {
        pout() << "Usage: <executable name> <inputfile>" << endl;
        pout() << "No input file specified" << endl;
        return -1;
      }
    ParmParse pp(a_argc-2,a_argv+2,NULL,inFile);
    // Run the main code for Convergence
    runConvergence(a_argc, a_argv);
  }
#ifdef CH_MPI
  CH_TIMER_REPORT();
  pout() << "about to call MPI Finalize" << std::endl;
  MPI_Finalize();
#endif
}

