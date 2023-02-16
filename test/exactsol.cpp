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
#include "EBProto.H"
#include "Chombo_EBChombo.H"
#include "Chombo_GeometryService.H"
#include "Chombo_EBLevelBoxData.H"
#include "Chombo_LAPACKMatrix.H"
#include "Proto_EBExactSolutions.H"
#include "Chombo_EBDictionary.H"

/***/
template< unsigned int order>
void
getPhi(EBLevelBoxData<CELL, 1>            & a_phi,
       shared_ptr<GeometryService<order> >& a_geoserv ,
       int a_nx, DisjointBoxLayout a_grids, Box a_domain, Real a_dx,
       BaseExactSolution& a_phigen)
{
  typedef IndexedMoments<DIM  , order> IndMomDIM;
  typedef HostIrregData<CELL,      IndMomDIM, 1>  VoluData;

  shared_ptr<Chombo4::LevelData<EBGraph>  > graphs = a_geoserv->getGraphs(a_domain);
  shared_ptr<LevelData<VoluData> > volDatLD= a_geoserv->getVoluData(a_domain);

  DataIterator dit = a_grids.dataIterator();
  for(int ibox = 0; ibox < dit.size(); ibox++)
  {
    Box chgrid = a_grids[dit[ibox]];
    Bx  prgrid = ProtoCh::getProtoBox(chgrid);
    auto graph =(*graphs)[dit[ibox]];
    Bx inputBox = a_phi[dit[ibox]].inputBox();
    EBHostData<CELL, Real, 1> hostdat(inputBox, graph, 1);
    auto& voldat = (*volDatLD)[dit[ibox]];
    for(auto bit = prgrid.begin(); bit != prgrid.end(); ++bit)
    {
      auto vofs = graph.getVoFs(*bit);
      for(int ivof = 0; ivof < vofs.size(); ivof++)
      {
        hostdat(vofs[ivof], 0) = a_phigen(graph, a_dx, voldat, vofs[ivof]);
      }
    }
    EBLevelBoxData<CELL, 1>::copyToDevice(a_phi[dit[ibox]], hostdat);
  }
}

template< unsigned int order>
void
average(EBLevelBoxData<CELL, 1> & a_diff,
        EBLevelBoxData<CELL, 1> & a_phiFine,
        EBLevelBoxData<CELL, 1> & a_phiCoar,
        shared_ptr<GeometryService<order> >& a_geoserv ,
        Vector<DisjointBoxLayout> a_grids,
        Vector<Box> a_domains,
        Vector<Real> a_dxes)
{
  typedef EBStencil<order, Real, CELL, CELL> ebstencil_t;
  IntVect dataGhostIV =a_diff.ghostVect();
  Point   dataGhostPt = ProtoCh::getPoint(dataGhostIV);
  
  shared_ptr<EBDictionary<order, Real, CELL, CELL> > 
    dictionary(new EBDictionary<order, Real, CELL, CELL>(a_geoserv, a_grids, a_domains, a_dxes, dataGhostPt));
  string nobcName("no_bcs");
  string restrictName = string("Restriction");
  Box domFine =a_domains[0];
  Box domCoar =a_domains[1];
  dictionary->registerStencil(restrictName, nobcName, nobcName, domFine, domCoar, false);
  
  DataIterator dit = a_grids[0].dataIterator();

  for(int ibox = 0; ibox < dit.size(); ++ibox)
  {
    auto& difffab  =    a_diff[dit[ibox]];
    auto& finecalc = a_phiFine[dit[ibox]];
    auto& coarcalc = a_phiCoar[dit[ibox]];
    //finer level owns the stencil (and the operator)

    shared_ptr<ebstencil_t> stencil =
      dictionary->getEBStencil(restrictName, nobcName, domFine, domCoar, ibox);
    //set resc = Ave(resf) (true is initToZero)
    stencil->apply(difffab, finecalc,  true, 1.0);
    difffab -= coarcalc;
  }
}
/***/
template< unsigned int order>
Real
getTruncationError(int a_nx)
{
  Box domFine(IntVect::Zero, (a_nx-1)*IntVect::Unit);
  Box domCoar = coarsen(domFine, 2);
  Vector<Box> boxes(1, domFine);
  Vector<int> procs(1, 0);
  int nxFine = domFine.size(0);
  int nxCoar = nxFine/2;
  Real dxFine = 1.0/domFine.size(0);
  Real dxCoar = 2.*dxFine;
  DisjointBoxLayout gridsFine(boxes, procs);
  DisjointBoxLayout gridsCoar;
  coarsen(gridsCoar, gridsFine, 2);
  
  int geomGhost = 2;
  RealVect ABC = RealVect::Unit();
  RealVect X0  = RealVect::Unit();
  RealVect origin= RealVect::Zero();
  Real R = 0.45;
  Real C = 0.5;
  X0 *= C;
//  for(int idir = 0; idir < DIM; idir++)
//  {
//    pout() << "direction = " << idir << endl;
//    for(int ipower = 0; ipower < 3; ipower++)
//    {
//      pout() << "power  = " << ipower << endl;
//    
//      std::pair<Point, Real> mono;
//      //mono.first = Point::Ones(1);
//      mono.first = Point::Zeros();
//      mono.first[idir] = ipower;
//      mono.second = 1;
//
//      //SineSphereEF<order> phigen(R, C);
//      MonomialEF<order> phigen(mono);
//
//
//      RealVect loca = RealVect::Unit();
//      loca *= 4.0;
//      IndexTM<int, DIM> curDeriv= IndexTM<int, DIM>::Zero;
//      for(int xderiv = 0; xderiv < 3; xderiv++)
//      {
//        curDeriv[0] = xderiv;
//        for(int yderiv = 0; yderiv < 3; yderiv++)
//        {
//          curDeriv[1] = yderiv;
//          int zderiv = 0;
//#if DIM==3
//          for(zderiv = 0; zderiv < 5; zderiv++)
//          {
//            curDeriv[2] = zderiv;
//#endif
//            int sum = curDeriv.sum();
//            if(sum < 5)
//            {
//              Real derivval = phigen.getDerivative(curDeriv, loca);
//              pout() << "(xd yd zd) derivval: ("
//                     << xderiv << ", "<< yderiv << ", " << zderiv << ")= "
//                     << derivval << std::endl;
//            }
//#if DIM==3
//          }
//#endif    
//        }
//      }
//    }
//  }

  Vector<DisjointBoxLayout>           grids(2);
  grids[0] = gridsFine;
  grids[1] = gridsCoar;
  shared_ptr<BaseIF>                  impfunc(new SimpleEllipsoidIF(ABC, X0, R, true));
  shared_ptr<GeometryService<order> > geoserv(new GeometryService<order>(impfunc, origin, dxFine, domFine, grids, geomGhost));

  
  pout() << "making data holders" << endl;
  IntVect dataGhostIV = IntVect::Unit;
  shared_ptr<Chombo4::LevelData<EBGraph>  > graphFine = geoserv->getGraphs(domFine);
  shared_ptr<Chombo4::LevelData<EBGraph>  > graphCoar = geoserv->getGraphs(domCoar);
  EBLevelBoxData<CELL, 1> phiFine(gridsFine, dataGhostIV, graphFine);
  EBLevelBoxData<CELL, 1> phiCoar(gridsCoar, dataGhostIV, graphCoar);
  EBLevelBoxData<CELL, 1> phiDiff(gridsCoar, dataGhostIV, graphCoar);
  
  pout() << "filling data holders" << endl;
  std::pair<Point, Real> mono;
  //mono.first = Point::Ones(1);
  mono.first = Point::Zeros();
  mono.first[0] = 1;
  mono.second = 1;

  //SineSphereEF<order> phigen(R, C);
  MonomialEF<order> phigen(mono);
  getPhi<order>(phiFine, geoserv, nxFine, gridsFine, domFine, dxFine, phigen);
  getPhi<order>(phiCoar, geoserv, nxCoar, gridsCoar, domCoar, dxCoar, phigen);

  Vector<Box> domains(2);
  domains[0] = domFine;
  domains[1] = domCoar;

  Vector<Real> dxes(2);
  dxes[0] = dxFine;
  dxes[1] = dxCoar;
  average<order>(phiDiff, phiFine, phiCoar, geoserv, grids, domains, dxes);

  phiFine.writeToFileHDF5("phiFine.hdf5", 0.0);
  phiCoar.writeToFileHDF5("phiCoar.hdf5", 0.0);
  phiDiff.writeToFileHDF5("phiDiff.hdf5", 0.0);
  
  Real   retval = phiDiff.maxNorm(0);
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
