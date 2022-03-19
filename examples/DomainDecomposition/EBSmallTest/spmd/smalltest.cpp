#include <cmath>
#include <cstdio>
#include <iostream>


#include "EBProto.H"
#include "Chombo_EBLevelBoxData.H"
#include "Chombo_EBLevelFluxData.H"
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
#include "EBAdvection.H"

#include "Chombo_EBDataChoreography.H"
#include "Chombo_ProtoFactories.H"

#include <iomanip>
#define GEOM_MAX_ORDER 4

using std::cout;
using std::endl;
using std::shared_ptr;
using Proto::Var;
using Proto::SimpleEllipsoidIF;
using CH4_Data_Choreography::DistributedData;
using Chombo4::EBLevelFluxData;
/***/
int getCorrectMooch(const Proto::Point& a_iv)
{
  int val = a_iv[0] + 10*a_iv[1];
  return val;
}
/***/
void
fillTheMooch(DistributedData<EBHostData<CELL, int, 1> >& a_mooch,
             shared_ptr<GeometryService<GEOM_MAX_ORDER>    >     a_geoserv,
             Chombo4::DisjointBoxLayout             a_grids,
             Chombo4::Box                           a_domain,
             Real a_dx)
{
  shared_ptr<LevelData<EBGraph> > graphs = a_geoserv->getGraphs(a_domain);
  for(Chombo4::DataIterator dit = a_grids.dataIterator(); dit.ok(); ++dit)
  {
    ChBx grid       = a_grids[dit()];
    Bx  grbx = ProtoCh::getProtoBox(grid);
    const auto& graph = (*graphs)[dit()];
    for(auto bit = grbx.begin(); bit != grbx.end(); ++bit)
    {
      vector<EBIndex<CELL> > vofs = graph.getVoFs(*bit);
      int val = getCorrectMooch(*bit);
      for(int ivof = 0; ivof < vofs.size(); ivof++)
      {
        a_mooch[dit()](vofs[ivof], 0) = val;
      }
    }
  }
}
/***/
int
checkTheMooch(DistributedData<EBHostData<CELL, int, 1> > & a_mooch,
              shared_ptr<GeometryService<GEOM_MAX_ORDER>    >           a_geoserv,
              Chombo4::DisjointBoxLayout                   a_grids,
              Chombo4::Box                                 a_domain,
              IntVect                                      a_ghost,
              Real a_dx)
{
  shared_ptr<LevelData<EBGraph> > graphs = a_geoserv->getGraphs(a_domain);
  for(ChDit dit = a_grids.dataIterator(); dit.ok(); ++dit)
  {
    ChBx grid       = a_grids[dit()];
    grid.grow(a_ghost);
    grid &= a_domain;
    Bx  grbx = ProtoCh::getProtoBox(grid);
    const auto& graph = (*graphs)[dit()];
    for(auto bit = grbx.begin(); bit != grbx.end(); ++bit)
    {
      vector<EBIndex<CELL> > vofs = graph.getVoFs(*bit);
      int correctVal= getCorrectMooch(*bit);
      for(int ivof = 0; ivof < vofs.size(); ivof++)
      {
        int exchangedVal =  a_mooch[dit()](vofs[ivof], 0);
        if(exchangedVal != correctVal)
        {
          Chombo4::pout() << "checkTheMooch: mismatch at " << *bit << ", correct = " << correctVal << ", actual = " << exchangedVal << endl;
          return -1;
        }
      }
    }
  }
  return 0;
}
/***/
int
testMinimalSPMD(shared_ptr<GeometryService<GEOM_MAX_ORDER>    >  a_geoserv,
                Chombo4::DisjointBoxLayout                       a_grids,
                Chombo4::Box                                     a_domain,
                Real a_dx, Point  a_ghost)
{
  int retval = 0;
  Chombo4::pout() << "entering testMinimalSPMD" << std::endl;
  IntVect numghost(a_ghost);
   GraphConstructorFactory<EBHostData<CELL, int, 1> > 
     factory(a_geoserv->getGraphs(a_domain));
   
  DistributedData<EBHostData<CELL, int, 1> > mooch(a_grids, a_ghost, factory);
  fillTheMooch(mooch, a_geoserv, a_grids, a_domain, a_dx);
  mooch.exchange(true);
  retval = checkTheMooch(mooch, a_geoserv, a_grids, a_domain, a_ghost, a_dx);
  Chombo4::pout() << "leaving testMinimalSPMD" << std::endl;
  return retval;
}
/***/
PROTO_KERNEL_START
unsigned int  exactSnoochF(int           a_pt[DIM],
                           Var<Real, 1>  a_phi)
{
  Real snoochval = a_pt[0] + 100*a_pt[1];
  a_phi(0) = snoochval;
  return 0;
}
PROTO_KERNEL_END(exactSnoochF, exactSnooch)
/***/
void
fillTheSnooch(EBLevelBoxData<CELL,  1>                       & a_snooch,
              shared_ptr<GeometryService<GEOM_MAX_ORDER> >     a_geoserv,
              Chombo4::DisjointBoxLayout                       a_grids,
              Chombo4::Box                                     a_domain,
              Real a_dx)
{
  shared_ptr<LevelData<EBGraph> > graphs = a_geoserv->getGraphs(a_domain);
  Chombo4::DataIterator dit = a_grids.dataIterator();
  for(int ibox = 0; ibox < dit.size(); ibox++)
  {
    ChBx grid       = a_grids[dit[ibox]];
    Bx   grbx = ProtoCh::getProtoBox(grid);
    auto& snoochfab = a_snooch[dit[ibox]];
    Bx  inputBx = snoochfab.inputBox();
    ebforall_i(inputBx, exactSnooch, grbx, snoochfab);
  }
}
/***/
PROTO_KERNEL_START
unsigned int  checkSnoochF(int           a_pt[DIM],
                           Var<Real, 1>  a_phi)
{
  Real snoochval = a_pt[0] + 100*a_pt[1];
  Real phival  = a_phi(0);
  Real tol = 0.0001;
  Real err = snoochval - phival;
  if((err > tol) || (err < -tol))
  {
    printf("error in checksnooch at point (%d, %d) \n", a_pt[0], a_pt[1]);
  }
  return 0;
}
PROTO_KERNEL_END(checkSnoochF, checkSnooch)
/***/
void
checkTheSnooch(EBLevelBoxData<CELL,  1>                       & a_snooch,
               shared_ptr<GeometryService<GEOM_MAX_ORDER> >     a_geoserv,
               Chombo4::DisjointBoxLayout                       a_grids,
               Chombo4::Box                                     a_domain,
               Real a_dx, IntVect a_ghost)
{
  shared_ptr<LevelData<EBGraph> > graphs = a_geoserv->getGraphs(a_domain);
  Chombo4::DataIterator dit = a_grids.dataIterator();
  for(int ibox = 0; ibox < dit.size(); ibox++)
  {
    ChBx grid       = a_grids[dit[ibox]];
    ChBx grown = grow(grid, a_ghost);
    grown &= a_domain;
    Bx   grbx = ProtoCh::getProtoBox(grown);
    
    auto& snoochfab = a_snooch[dit[ibox]];
    Bx  inputBx = snoochfab.inputBox();
    ebforall_i(inputBx, checkSnooch, grbx, snoochfab);
  }
}

/***/
int
testEBLevelBoxData(shared_ptr<GeometryService<GEOM_MAX_ORDER> >  a_geoserv,
                   Chombo4::DisjointBoxLayout                    a_grids,
                   Chombo4::Box                                  a_domain,
                   Real a_dx, Point  a_ghost)
{
  auto graphs = a_geoserv->getGraphs(a_domain);
  Chombo4::pout() << "starting testEBLevelBoxData" << std::endl;
  IntVect numghost(a_ghost);
  
  EBLevelBoxData<CELL, 1>  snooch(a_grids, numghost, graphs);
  fillTheSnooch(snooch, a_geoserv, a_grids, a_domain, a_dx);
  snooch.exchange(true);
  checkTheSnooch(snooch, a_geoserv, a_grids, a_domain, a_dx, numghost);
  Chombo4::pout() << "leaving testEBLevelBoxData" << std::endl;
  return 0;
}
/***/
PROTO_KERNEL_START
unsigned int  exactFloochF(int           a_pt[DIM],
                           Var<Real, 1>  a_phi,
                           int facedir)
{
  Real snoochval = a_pt[0] + 100*a_pt[1] + 1000*facedir;
  a_phi(0) = snoochval;
  return 0;
}
PROTO_KERNEL_END(exactFloochF, exactFlooch)
/***/
void
fillTheFlooch(EBLevelFluxxData<1>                            & a_flooch,
              shared_ptr<GeometryService<GEOM_MAX_ORDER> >     a_geoserv,
              Chombo4::DisjointBoxLayout                       a_grids,
              Chombo4::Box                                     a_domain,
              Real a_dx)
{
  shared_ptr<LevelData<EBGraph> > graphs = a_geoserv->getGraphs(a_domain);
  Chombo4::DataIterator dit = a_grids.dataIterator();
  for(int ibox = 0; ibox < dit.size(); ibox++)
  {
    ChBx grid       = a_grids[dit[ibox]];
    Bx   grbx = ProtoCh::getProtoBox(grid);
    auto& floochfab = a_flooch[dit[ibox]];
    {
      auto& xflooch = *floochfab.m_xflux;
      Bx  xinputBx  = xflooch.inputBox();
      ebforall_i(xinputBx, exactFlooch, grbx, xflooch, 0);
    }
    {
      auto& yflooch = *floochfab.m_yflux;
      Bx  yinputBx = yflooch.inputBox();
      ebforall_i(yinputBx, exactFlooch, grbx, yflooch, 1);
    }
#if DIM==3
    {
      auto& zflooch = *floochfab.m_zflux;
      Bx  zinputBx = zflooch.inputBox();
      ebforall_i(zinputBx, exactFlooch, grbx, zflooch, 2);
    }
#endif    
  }
}
/***/
PROTO_KERNEL_START
unsigned int  checkFloochF(int           a_pt[DIM],
                           Var<Real, 1>  a_phi,
                           int a_facedir)
{
  Real snoochval = a_pt[0] + 100*a_pt[1] + 1000*a_facedir;
  Real phival  = a_phi(0);
  Real tol = 0.0001;
  Real err = snoochval - phival;
  if((err > tol) || (err < -tol))
  {
    printf("error in checkflooch at point (%d, %d) \n", a_pt[0], a_pt[1]);
  }
  return 0;
}
PROTO_KERNEL_END(checkSnoochF, checkSnooch)
/***/
void
checkTheFlooch(EBLevelFluxData<1>                             & a_flooch,
               shared_ptr<GeometryService<GEOM_MAX_ORDER> >     a_geoserv,
               Chombo4::DisjointBoxLayout                       a_grids,
               Chombo4::Box                                     a_domain,
               Real a_dx, IntVect a_ghost)
{
  shared_ptr<LevelData<EBGraph> > graphs = a_geoserv->getGraphs(a_domain);
  Chombo4::DataIterator dit = a_grids.dataIterator();
  for(int ibox = 0; ibox < dit.size(); ibox++)
  {
    ChBx grid       = a_grids[dit[ibox]];
    ChBx grown = grow(grid, a_ghost);
    grown &= a_domain;
    Bx   grbx = ProtoCh::getProtoBox(grown);
    
    auto& floochfab = a_snooch[dit[ibox]];
    {
      auto& xfloochfab = *floochfab.m_xflux;
      Bx  xinputBx = xfloochfab.inputBox();
      ebforall_i(xinputBx, checkFlooch, grbx, xfloochfab, 0);
    }
    {
      auto& yfloochfab = *floochfab.m_yflux;
      Bx  yinputBx = yfloochfab.inputBox();
      ebforall_i(yinputBx, checkFlooch, grbx, yfloochfab, 1);
    }
#if DIM==3
    {
      auto& zfloochfab = *floochfab.m_zflux;
      Bx  zinputBx = zfloochfab.inputBox();
      ebforall_i(zinputBx, checkFlooch, grbx, zfloochfab, 2);
    }
#endif    
  }
}
int
testEBLevelFluxData(shared_ptr<GeometryService<GEOM_MAX_ORDER> >  a_geoserv,
                    Chombo4::DisjointBoxLayout                    a_grids,
                    Chombo4::Box                                  a_domain,
                    Real a_dx, Point  a_ghost)
{
  auto graphs = a_geoserv->getGraphs(a_domain);
  Chombo4::pout() << "starting testEBFluxBoxData" << std::endl;
  IntVect numghost(a_ghost);
  
  EBLevelFluxData<1>  flooch(a_grids, numghost, graphs);
  fillTheFlooch(snooch, a_geoserv, a_grids, a_domain, a_dx);
  flooch.exchange(true);
  checkTheFlooch(snooch, a_geoserv, a_grids, a_domain, a_dx, numghost);
  Chombo4::pout() << "leaving testEBLevelFluxData" << std::endl;
  return 0;
}

/***/
void makeGrids(Chombo4::DisjointBoxLayout& a_grids,
               Real             & a_dx,
               const int        & a_nx,
               const int        & a_maxGrid)
{
  Chombo4::pout() << "making grids" << endl;

  IntVect domLo = IntVect::Zero;
   IntVect domHi  = (a_nx - 1)*IntVect::Unit;

  // EB and periodic do not mix
  Chombo4::ProblemDomain domain(domLo, domHi);

  Vector<Chombo4::Box> boxes;
  unsigned int blockfactor = 8;
  domainSplit(domain, boxes, a_maxGrid, blockfactor);
  
  Vector<int> procs;

  a_dx = 1.0/a_nx;
  LoadBalance(procs, boxes);
  a_grids = Chombo4::DisjointBoxLayout(boxes, procs, domain);
  a_grids.printBalance();

}

///
shared_ptr<BaseIF>  getImplicitFunction(Real  & a_geomCen,
                                        Real  & a_geomRad,
                                        int   & a_whichGeom)

{
  using Proto::BaseIF;
  shared_ptr<BaseIF>  retval;
  ParmParse pp;
  
  a_geomCen = 0;
  a_geomRad = 1;
  pp.get("which_geom", a_whichGeom);
  if(a_whichGeom == -1)
  {
    using Proto::AllRegularIF;
    Chombo4::pout() << "all regular geometry" << endl;
    retval = shared_ptr<BaseIF>(new AllRegularIF());
  }
  else if(a_whichGeom == 0)
  {
    using Proto::SimpleEllipsoidIF;
    Chombo4::pout() << "sphere" << endl;

    pp.get("geom_cen", a_geomCen);
    pp.get("geom_rad", a_geomRad);
    Chombo4::pout() << "geom_cen = " << a_geomCen       << endl;
    Chombo4::pout() << "geom_rad = " << a_geomRad       << endl;

    RealVect ABC = RealVect::Unit(); //this is what it makes it a sphere instead of an ellipse
    RealVect  X0 = RealVect::Unit();
    X0 *= a_geomCen;

    retval = shared_ptr<BaseIF>(new SimpleEllipsoidIF(ABC, X0, a_geomRad, true));//true is for inside regular
  }
  else if(a_whichGeom ==  1)
  {
    using Proto::PlaneIF;
    Chombo4::pout() << "plane" << endl;
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
//=================================================
void defineGeometry(Chombo4::DisjointBoxLayout& a_grids,
                    Real             & a_dx,
                    Real             & a_geomCen,
                    Real             & a_geomRad,
                    int              & a_whichGeom,
                    int              & a_nx,
                    shared_ptr<GeometryService<GEOM_MAX_ORDER> >&  a_geoserv)
{
  Chombo4::pout() << "defining geometry" << endl;

  ParmParse pp;
  int maxGrid = 32;
    
  pp.get("nx"        , a_nx);
  pp.get("maxGrid", maxGrid);

  Chombo4::pout() << "nx       = " << a_nx     << endl;
  Chombo4::pout() << "maxGrid  = " << maxGrid  << endl;

  makeGrids(a_grids, a_dx, a_nx, maxGrid);
  Chombo4::Box domain = a_grids.physDomain().domainBox();
  int geomGhost = 4;
  RealVect origin = RealVect::Zero();

  Chombo4::pout() << "creating implicit function" << endl;
  shared_ptr<BaseIF>  impfunc = getImplicitFunction(a_geomCen, a_geomRad, a_whichGeom);

  Chombo4::pout() << "creating geometry service" << endl;
  a_geoserv  = shared_ptr<GeometryService<GEOM_MAX_ORDER> >(new GeometryService<GEOM_MAX_ORDER>(impfunc, origin, a_dx, domain, a_grids, geomGhost));
}

int
runTests(int a_argc, char* a_argv[])
{

  using Chombo4::ProblemDomain;
  using Chombo4::DisjointBoxLayout;
  using Chombo4::LevelBoxData;
  using Chombo4::Copier;
  using Chombo4::DataIterator;
  using Chombo4::MayDay;
  using Proto::BaseIF;
  
  Real dx;
  Chombo4::DisjointBoxLayout grids;

  shared_ptr<GeometryService<GEOM_MAX_ORDER> >  geoserv;

  Real geomCen;
  Real geomRad;
  int whichGeom;
  int nx;
  defineGeometry(grids, dx, geomCen, geomRad, whichGeom, nx,  geoserv);

  IntVect dataGhostIV =   4*IntVect::Unit;
  Point   dataGhostPt = ProtoCh::getPoint(dataGhostIV); 

  Chombo4::pout() << "making dictionary" << endl;
  Chombo4::Box domain = grids.physDomain().domainBox();
  int retval = testMinimalSPMD(geoserv, grids, domain, dx, dataGhostPt);
  if(retval != 0)
  {
    Chombo4::pout() << "problem in testMinimalSPMD" << endl;
  }
  //because this gets tested on the device , no return values are possible.
  //why yes, it is a lovely way to run a computer.
  testEBLevelBoxData(geoserv, grids, domain, dx, dataGhostPt);

  return retval;
}


int main(int a_argc, char* a_argv[])
{
#ifdef CH_MPI
  MPI_Init(&a_argc, &a_argv);
  Chombo4::pout() << "MPI INIT called" << std::endl;
#endif
  //needs to be called after MPI_Init
  CH_TIMER_SETFILE("ebadvect.time.table");
  int retval = 0;
  {
    if (a_argc < 2)
    {
      cerr<< " usage " << a_argv[0] << " <input_file_name> " << endl;
      exit(0);
    }
    char* in_file = a_argv[1];
    ParmParse  pp(a_argc-2,a_argv+2,NULL,in_file);
    retval = runTests(a_argc, a_argv);
    
  }

  Chombo4::pout() << "printing time table " << endl;
  CH_TIMER_REPORT();
#ifdef CH_MPI
  Chombo4::pout() << "about to call MPI Finalize" << std::endl;
  MPI_Finalize();
#endif
  return retval;
}


