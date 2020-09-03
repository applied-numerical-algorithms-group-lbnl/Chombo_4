#include <cstdio>
#include <cstring>
#include <cassert>
#include <cmath>

#include <vector>
#include <memory>

#include <iostream>
#include <fstream>
#include <sstream>

#include "Proto.H"
#include "Chombo_LAPACKMatrix.H"

PROTO_KERNEL_START 
unsigned int  setTrigStuffF(int              a_pt[DIM],
                            Var<Real, 1>     a_lphCalc,
                            Var<Real, 1>     a_lphExac,
                            Var<Real, 1>     a_phiExac,
                            Real             a_rad,
                            Real             a_cen,
                            Real             a_dx)
{
  Real pi = 4.*atan(1.0);
  Real x = Real(a_pt[0])*(a_dx + 0.5) - a_cen;
  Real y = Real(a_pt[1])*(a_dx + 0.5) - a_cen;
  Real rsq = x*x + y*y;
  Real sinval = sin(pi*x)*sin(pi*y);
#if DIM==3  
  Real z = Real(a_pt[1])*(a_dx + 0.5);
  rsq += z*z;
  sinval *= sin(pi*z);
#endif
  Real cosval = cos(pi*(rsq -(R*R)));
  a_phiExac(0) = cosval;
  a_lphExac(0) = -pi*pi*cosval;
  a_lphCalc(0) = 0;
  return 0;
}
PROTO_KERNEL_END(setTrigStuffF, setTrigStuff)

vector<EBIndex<CELL> > getNeighbors(unsigned int          a_radius,
                                    const EBIndex<CELL> & a_start,
                                    const EBGraph       & a_graph)
{
  vector<EBIndex<CELL> > retval;
  Bx ptbox(a_start.m_pt, a_start.m_pt);
  Bx area = ptbox.grow(a_radius);
  area &=  a_graph.getDomain();
  for(auto bit = area.begin(); bit != bit.end(); ++bit)
  {
    vector<EBIndex<CELL> > cellvofs = a_graph.getVoFs(*bit);
    retval.insert(retval.end(), cellvofs.begin(), cellvofs.end());
  }
  return retval;
}
/***/
std::shared_ptr< AggStencil<CELL, CELL, Real> >
getAggSten(const EBGraph & a_graph,
           const Box     & a_valid,
           const Point   & a_ghost)
{
  vector< EBIndex<CELL> > dstVoFs;
  vector< LocalStencil< CELL, Real> > localstens;
  
  std::shared_ptr< AggStencil<CELL, CELL, Real> > retval(dstVoFs,
                                                         localstens,
                                                         graph, graph,
                                                         valid, valid,
                                                         ghost, ghost);
  return retval;
}
/***/
template< unsigned int order>
Real
getTruncationError(int a_nx)
{
  using Chombo4::Box;
  using Chombo4::DisjointBoxLayout;

  IntVect dataGhostIV =   3*IntVect::Unit;
  Point   dataGhostPt = ProtoCh::getPoint(dataGhostIV); 

  Box domain(IntVect::Zero, (a_nx-1)*IntVect::Unit);
  Vector<Box> boxes(1, domain);
  Vector<int> procs(1, 0);
  Real dx = 1.0/domain.size(0);
  std::array<bool, DIM> periodic;
  for(int idir = 0; idir < DIM; idir++) periodic[idir]=true;
  DisjointBoxLayout grids(boxes, procs);

  int geomGhost = 4;

  RealVect ABC = RealVect::Unit();
  RealVect X0  = RealVect::Unit();
  
  RealVect origin= RealVect::Zero();
  Real R = 0.45;
  Real C = 0.5;
  X0 *= C;
  shared_ptr<BaseIF>                  impfunc(new SimpleEllipsoidIF(ABC, X0, R, true));
  shared_ptr<GeometryService<order> > geoserv(new GeometryService<2>(impfunc, origin, dx, domain, grids, geomGhost));
  shared_ptr<Chombo4::LevelData<EBGraph>  > graphs = geoserv->getGraphs(domain);
  pout() << "making data" << endl;
  EBLevelBoxData<CELL,   1>  phiExacLD(grids, dataGhostIV, graphs);
  EBLevelBoxData<CELL,   1>  lphExacLD(grids, dataGhostIV, graphs);
  EBLevelBoxData<CELL,   1>  lphCalcLD(grids, dataGhostIV, graphs);
  
  DataIterator dit = grids.dataIterator();
  for(int ibox = 0; ibox < dit.size(); ibox++)
  {
    Box     valid =     grids[dit[ibox]];
    EBGraph graph =    graphs[dit[ibox]];
    auto& lphCalc = lphCalcLD[dit[ibox]];
    auto& lphExac = lphExacLD[dit[ibox]];
    auto& phiExac = phiExacLD[dit[ibox]];
    ebforall_i(inputBox, setTrigStuff, valid, lphCalc, lphExac, phiExac, R, C, dx);
    std::shared_ptr< AggStencil<CELL, CELL, Real> > aggsten = getAggSten(graph, valid, ghost);
    aggsten->apply(lphCalc, phiCalc, 1.0, false);
  }

    
}
/***/
int main(int argc, char* argv[])
{
  Real nxFine = 64;
  Real nxCoar = nxFine/2;
  
  Real errFine = getTruncationError<4>(nxFine);
  Real errCoar = getTruncationError<4>(nxCoar);
  int  iresult = 0;
  Real order = 0;
  if(errFine > 1.0e-6)
  {
    order = log(errCoar/errFine)/log(2);
    if(order < 0.5) //because I am using pointwise functions, 2nd order is the best you can do here
    {
      iresult =  1;
    }
  }
  if(iresult == 0)
  {
    printf("hoeb test PASSED with order = %f \n", order);
  }
  else
  {
    printf("hoeb test FAILED with order =  %f\n", order);
  }
  return 0;

}
