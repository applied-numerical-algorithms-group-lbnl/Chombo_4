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
testMinimalSPMD(shared_ptr<GeometryService<GEOM_MAX_ORDER>    >     a_geoserv,
                Chombo4::DisjointBoxLayout             a_grids,
                Chombo4::Box                           a_domain,
                Real a_dx, Point  a_ghost)
{
  IntVect numghost(a_ghost);
   GraphConstructorFactory<EBHostData<CELL, int, 1> > 
     factory(a_geoserv->getGraphs(a_domain));
   
  DistributedData<EBHostData<CELL, int, 1> > mooch(a_grids, a_ghost, factory);
  fillTheMooch(mooch, a_geoserv, a_grids, a_domain, a_dx);
  mooch.exchange();
  int retval = checkTheMooch(mooch, a_geoserv, a_grids, a_domain, a_ghost, a_dx);
  return retval;
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

  Chombo4::pout() << "defining geometry" << endl;
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


