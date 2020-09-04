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
template< unsigned int order>
Real
getTruncationError(int a_nx)
{
  EBBoxData<CELL, 1> phiFine, phiMedi, diff;
  getPhi(phiFine, nx);
  getPhi(phiCoar, nx/2);

  average(diff, phiFine, phiCoar);

  Real   retval = maxnorm(diff);
  return retval;
    
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
