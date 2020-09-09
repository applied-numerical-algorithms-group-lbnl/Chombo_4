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
void
getPhi(EBBoxData<CELL, 1>& phi, int a_nx, DisjointBoxLayout a_grids, Box a_domain)
{
  IntVect dataGhostIV =   IntVect::Unit;
  Point   dataGhostPt = ProtoCh::getPoint(dataGhostIV); 
  int geomGhost = 2;
  RealVect ABC = RealVect::Unit();
  RealVect X0  = RealVect::Unit();
  RealVect origin= RealVect::Zero();
  Real R = 0.45;
  Real C = 0.5;
  X0 *= C;
  SineSphere<order> phigen(R, C);
  shared_ptr<BaseIF>                  impfunc(new SimpleEllipsoidIF(ABC, X0, R, true));
  shared_ptr<GeometryService<order> > geoserv(new GeometryService<2>(impfunc, origin, dx, a_domain, a_grids, geomGhost));
  shared_ptr<Chombo4::LevelData<EBGraph>  > graphs = geoserv->getGraphs(a_domain);
  pout() << "making data" << endl;
  a_phi.define(a_grids, dataGhostPt, graphs);

  DataIterator dit = a_grids.dataIterator();
  for(int ibox = 0; ibox < dit.size(); ibox++)
  {
    Box chgrid = dit[dit[ibox]];
    Bx  prgrid = ProtoCh::getProtoBox(chgrid);
    auto graph = graphs[dit[ibox]];
    EBHostData<CELL, 1> hostdat(prgrid, graph, 1);
    for(auto bit = prgrid.begin(); bit != prgrid.end(); ++bit)
    {
      auto vofs = graph.getVoFs(*bit);
      for(int ivof = 0; ivof < vofs.size(); ivof++)
      {
        hostdat(vofs[ivof], 0) = phigen(graph, dx, voldat, vofs[ivof]);
      }
    }
    EBLevelBoxData<CELL, 1>::copyToDevice(a_phi[dit[ibox]], hostdat);
  }
}

void
average(EBLevelBoxData<CELL, 1> & a_diff,
        EBLevelBoxData<CELL, 1> & a_phiFine,
        EBLevelBoxData<CELL, 1> & a_phiCoar)
{
}
/***/
template< unsigned int order>
Real
getTruncationError(int a_nx)
{
  Box domFine(IntVect::Zero, (a_nx-1)*IntVect::Unit);
  Box domCoar = coarsen(domFine, 2);
  Vector<Box> boxes(1, domain);
  Vector<int> procs(1, 0);
  Real dx = 1.0/domain.size(0);
  DisjointBoxLayout gridsFine(boxes, procs);
  DisjointBoxLayout gridsCoar;
  coarsen(gridsCoar, gridsFine, 2);
  
  EBLevelBoxData<CELL, 1> phiFine, phiMedi, diff;
  getPhi(phiFine, nx,   gridsFine, domFine);
  getPhi(phiCoar, nx/2, gridsCoar, domCoar);

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
