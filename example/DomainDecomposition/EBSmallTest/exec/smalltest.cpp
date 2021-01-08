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
#include "SetupFunctions.H"

#include <iomanip>

#define MAX_ORDER 2

using std::cout;
using std::endl;
using std::shared_ptr;
using Proto::Var;
using Proto::SimpleEllipsoidIF;

typedef Var<Real,DIM> Vec;
typedef Var<Real,  1> Sca;
#if DIM==2
void 
testDebugFunctions(const EBGraph& a_graph)
{
  Point lo(3, 5);
  Point hi(7, 9);
  Bx bx(lo, hi);
  EBBoxData<CELL, Real, 1> cellfab(bx, a_graph);
  unsigned int numflops = 0;
  pout() << "testing dumpEB1" << endl;
  ebforallInPlace_i(numflops, "InitCell", InitCell,  bx, cellfab);
  dumpEB1(&cellfab);

  EBFluxData<Real, 1> fluxdat(bx, a_graph);
  ebforallInPlace_i(numflops, "InitCell", InitCell,  bx, cellfab);
  
  pout() << "testing dumpXFace" << endl;
  int idir = 0;
  ebforallInPlace_i(numflops, "InitFace", InitFace, fluxdat.m_xflux->box(),
                    *fluxdat.m_xflux, idir);
  dumpXFace(&fluxdat.m_xflux);


  pout() << "testing dumpYFace" << endl;
  idir = 1;
  ebforallInPlace_i(numflops, "InitFace", InitFace, fluxdat.m_yflux->box(),
                    *fluxdat.m_yflux, idir);
  dumpYFace(&fluxdat.m_yflux);

}

/***/
void 
testHiStencil(shared_ptr<EBEncyclopedia<2, Real> >   a_brit,
              shared_ptr<GeometryService<2>    >   a_geoserv,
              Chombo4::DisjointBoxLayout                     a_grids,
              Chombo4::Box                                   a_domain,
              Real a_dx, Point  a_ghost)
{
  using Chombo4::ProblemDomain;
  using Chombo4::DisjointBoxLayout;
  using Chombo4::LevelBoxData;
  using Chombo4::Copier;
  using Chombo4::DataIterator;
  using Chombo4::MayDay;
  using Proto::BaseIF;
  
  pout() << "testing hi cell to face stencil" << endl;
  shared_ptr<LevelData<EBGraph> > graphs = a_geoserv->getGraphs(a_domain);

  a_brit->registerCellToFace(StencilNames::CellToFaceHi, 
                             StencilNames::NoBC,
                             StencilNames::NoBC,
                             a_domain, a_domain, false, Point::Ones());

  DataIterator dit = a_grids.dataIterator();
  for(int ibox = 0; ibox < dit.size(); ibox++)
  {
    Chombo4::Box grid = a_grids[dit[ibox]];
    Bx  bx = ProtoCh::getProtoBox(grid);
    Bx grown = bx.grow(a_ghost);
    EBGraph graph = (*graphs)[dit[ibox]];

    EBBoxData<CELL, Real, 1> cellfab(grown, graph);
    unsigned int numflops = 0;
    pout() << "cell data looks like:" << endl;
    ebforallInPlace_i(numflops, "InitCell", InitCell,  grown, cellfab);

    dumpEB1(&cellfab);

    EBFluxData<Real, 1> fluxdat(grown, graph);
    int idir = 0;
    pout() << "XFace data for Side::Lo looks like:" << endl;
    a_brit->applyCellToFace(StencilNames::CellToFaceHi, 
                            StencilNames::NoBC,
                            a_domain, fluxdat, cellfab, idir, ibox, true, 1.0);

    dumpXFace(&fluxdat.m_xflux);
  }
}
/***/
void 
testLoStencil(shared_ptr<EBEncyclopedia<2, Real> >       a_brit,
                   shared_ptr<GeometryService<2>    >    a_geoserv,
                   Chombo4::DisjointBoxLayout                     a_grids,
              Chombo4::Box                                   a_domain,
                   Real a_dx, Point  a_ghost)
{
  using Chombo4::ProblemDomain;
  using Chombo4::DisjointBoxLayout;
  using Chombo4::LevelBoxData;
  using Chombo4::Copier;
  using Chombo4::DataIterator;
  using Chombo4::MayDay;
  using Proto::BaseIF;
  
  pout() << "testing low cell to face stencil" << endl;
  shared_ptr<LevelData<EBGraph> > graphs = a_geoserv->getGraphs(a_domain);

  a_brit->registerCellToFace(StencilNames::CellToFaceLo, 
                             StencilNames::NoBC,
                             StencilNames::NoBC,
                             a_domain, a_domain, false, Point::Ones());

  DataIterator dit = a_grids.dataIterator();
  for(int ibox = 0; ibox < dit.size(); ibox++)
  {
    auto grid = a_grids[dit[ibox]];
    Bx  bx = ProtoCh::getProtoBox(grid);
    Bx grown = bx.grow(a_ghost);
    EBGraph graph = (*graphs)[dit[ibox]];

    EBBoxData<CELL, Real, 1> cellfab(grown, graph);
    unsigned int numflops = 0;
    pout() << "cell data looks like:" << endl;
    ebforallInPlace_i(numflops, "InitCell", InitCell,  grown, cellfab);
    dumpEB1(&cellfab);

    EBFluxData<Real, 1> fluxdat(grown, graph);
    int idir = 0;
    a_brit->applyCellToFace(StencilNames::CellToFaceLo, 
                            StencilNames::NoBC,
                            a_domain, fluxdat, cellfab, idir, ibox, true, 1.0);

    pout() << "XFace data for Side::Hi data looks like:" << endl;
    dumpXFace(&fluxdat.m_xflux);
  }
}
/***/
void 
testSimpleStencils(shared_ptr<EBEncyclopedia<2, Real> >   a_brit,
                   shared_ptr<GeometryService<2>    >     a_geoserv,
                   Chombo4::DisjointBoxLayout                      a_grids,
                   Chombo4::Box                                    a_domain,
                   Real a_dx, Point  a_ghost)
{
  pout() << "testing cell to face stencil" << endl;
  testHiStencil(a_brit, a_geoserv, a_grids, a_domain, a_dx, a_ghost);
  testLoStencil(a_brit, a_geoserv, a_grids, a_domain, a_dx, a_ghost);
}
/***/
#endif

//=================================================

void makeGrids(Chombo4::DisjointBoxLayout& a_grids,
               Real             & a_dx,
               const int        & a_nx,
               const int        & a_maxGrid)
{
  pout() << "making grids" << endl;

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
    pout() << "all regular geometry" << endl;
    retval = shared_ptr<BaseIF>(new AllRegularIF());
  }
  else if(a_whichGeom == 0)
  {
    using Proto::SimpleEllipsoidIF;
    pout() << "sphere" << endl;

    pp.get("geom_cen", a_geomCen);
    pp.get("geom_rad", a_geomRad);
    pout() << "geom_cen = " << a_geomCen       << endl;
    pout() << "geom_rad = " << a_geomRad       << endl;

    RealVect ABC = RealVect::Unit(); //this is what it makes it a sphere instead of an ellipse
    RealVect  X0 = RealVect::Unit();
    X0 *= a_geomCen;

    retval = shared_ptr<BaseIF>(new SimpleEllipsoidIF(ABC, X0, a_geomRad, true));//true is for inside regular
  }
  else if(a_whichGeom ==  1)
  {
    using Proto::PlaneIF;
    pout() << "plane" << endl;
    RealVect normal, startPt;
    vector<double> v_norm, v_start;
    pp.getarr("geom_normal", v_norm, 0, DIM);
    pp.getarr("geom_start_pt", v_start, 0, DIM);
    for(int idir = 0; idir < DIM; idir++)
    {
      normal[ idir] = v_norm[ idir];
      startPt[idir] = v_start[idir];
      pout() << "normal ["<< idir << "] = " << normal [idir]  << endl;
      pout() << "startPt["<< idir << "] = " << startPt[idir]  << endl;
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
                    shared_ptr<GeometryService<MAX_ORDER> >&  a_geoserv)
{
  pout() << "defining geometry" << endl;

  ParmParse pp;
  int maxGrid = 32;
    
  pp.get("nx"        , a_nx);
  pp.get("maxGrid", maxGrid);

  pout() << "nx       = " << a_nx     << endl;
  pout() << "maxGrid  = " << maxGrid  << endl;

  makeGrids(a_grids, a_dx, a_nx, maxGrid);
  Chombo4::Box domain = a_grids.physDomain().domainBox();
  int geomGhost = 4;
  RealVect origin = RealVect::Zero();

  pout() << "creating implicit function" << endl;
  shared_ptr<BaseIF>  impfunc = getImplicitFunction(a_geomCen, a_geomRad, a_whichGeom);

  pout() << "creating geometry service" << endl;
  a_geoserv  = shared_ptr<GeometryService<MAX_ORDER> >(new GeometryService<MAX_ORDER>(impfunc, origin, a_dx, domain, a_grids, geomGhost));
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

  pout() << "defining geometry" << endl;
  shared_ptr<GeometryService<MAX_ORDER> >  geoserv;

  Real geomCen;
  Real geomRad;
  int whichGeom;
  int nx;
  defineGeometry(grids, dx, geomCen, geomRad, whichGeom, nx,  geoserv);

  IntVect dataGhostIV =   4*IntVect::Unit;
  Point   dataGhostPt = ProtoCh::getPoint(dataGhostIV); 

  pout() << "making dictionary" << endl;
  Chombo4::Box domain = grids.physDomain().domainBox();
  shared_ptr<EBEncyclopedia<2, Real> > 
    brit(new EBEncyclopedia<2, Real>(geoserv, grids, domain, dx, dataGhostPt));

  testSimpleStencils(brit, geoserv, grids, domain, dx, dataGhostPt);

  return 0;
}


int main(int a_argc, char* a_argv[])
{
#ifdef CH_MPI
  MPI_Init(&a_argc, &a_argv);
  pout() << "MPI INIT called" << std::endl;
#endif
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
    runTests(a_argc, a_argv);
  }

  pout() << "printing time table " << endl;
  CH_TIMER_REPORT();
#ifdef CH_MPI
  pout() << "about to call MPI Finalize" << std::endl;
  MPI_Finalize();
#endif
  return 0;
}


