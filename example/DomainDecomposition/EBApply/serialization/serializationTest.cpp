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
#include "Chombo_EBLevelBoxData.H"
#include "Chombo_EBLevelFluxData.H"
#include "DebugFunctions.H"
#include "implem/Proto_HostIrregData.H"
#include <iomanip>

#include "Chombo_NamespaceHeader.H"

#define MAX_ORDER 2

using std::cout;
using std::endl;
using std::shared_ptr;

shared_ptr<vector<EBIndex<CELL> > > getIndicies(const Bx & a_bx)
{
  shared_ptr<vector<EBIndex<CELL> > > retval( new vector<EBIndex<CELL> >);
  for(auto bit = a_bx.begin(); bit != a_bx.end(); ++bit)
  {
    EBIndex<CELL> index;
    index.m_pt = *bit;
    index.m_vofIDMe = 0;
    index.m_vofIDLo = 0;
    index.m_isBoundary = false;
    retval->push_back(index);
  }
  return retval;
}

int
test1(int a_argc, char* a_argv[])
{
  cout << "entering test1" << endl;
  Bx unitBx(Point::Zeroes(), Point::Ones());
  shared_ptr<vector<EBIndex<CELL> > > indicies = getIndicies(unitBx);
  HostIrregData<CELL, IndexedMoments<2,2>, 1> holder;
  holder.define(unitBx, indicies);
  cout << "leaving  test1" << endl;
  return 0;
}

int
test2(int a_argc, char* a_argv[])
{
  cout << "entering test2" << endl;
  Bx srcBx(Point::Zeroes(), Point::Ones( ));
  Bx dstBx(Point::Ones(),   Point::Ones(2));
  Bx intBx = srcBx & dstBx;
  Real srcVal = 1.;
  Real dstVal = 2.;
  
  shared_ptr<vector<EBIndex<CELL> > > srcInd = getIndicies(srcBx);
  shared_ptr<vector<EBIndex<CELL> > > dstInd = getIndicies(dstBx);
  shared_ptr<vector<EBIndex<CELL> > > intInd = getIndicies(intBx);

  HostIrregData<CELL, Real, 1> src, dst;
  cout << "defining data" << endl;
  src.define(srcBx, srcInd);
  dst.define(dstBx, dstInd);
  src.setVal(srcVal);
  dst.setVal(dstVal);
  for(int iind = 0; iind < dstInd->size(); iind++)
  {
    auto intMom = dst((*dstInd)[iind], 0);
    intMom -= dstVal;
    auto diff = std::abs(intMom);
    if(diff > 1.0e-6)
    {
      MayDay::Abort("bad value in test 2: 1");
    }
  }
  for(int iind = 0; iind < srcInd->size(); iind++)
  {
    auto intMom = src((*srcInd)[iind], 0);
    intMom -= srcVal;
    auto diff = std::abs(intMom);
    if(diff > 1.0e-6)
    {
      MayDay::Abort("bad value in test 2: 2");
    }
  }

  cout << "doing serialization dance" << endl;
  int sizeSrc = src.charsize(intBx, 0, 1);
  int sizeDst = dst.charsize(intBx, 0, 1);
  if (sizeSrc != sizeDst)
  {
    MayDay::Abort("test2 LinearizationTest failure: dest and source have different sizes");
  }
  Vector<char> buffer(sizeSrc);
  void* buf = (void*)&(buffer[0]);
  src.linearOut(buf, intBx, 0, 1);
  dst.linearIn (buf, intBx, 0, 1);

  cout << "checking the answer" << endl;
  for(int iind = 0; iind < intInd->size(); iind++)
  {
    auto intMom = dst((*intInd)[iind], 0);
    intMom -= srcVal;
    auto diff = std::abs(intMom);
    if(diff > 1.0e-6)
    {
      MayDay::Abort("bad value in test 2: 3");
    }
  }
  cout << "leaving  test2" << endl;
  return 0;
}
int
test3(int a_argc, char* a_argv[])
{
  cout << "entering test3" << endl;
  Bx srcBx(Point::Zeroes(), Point::Ones( ));
  Bx dstBx(Point::Ones(),   Point::Ones(2));
  Bx intBx = srcBx & dstBx;
  IndexedMoments<2,2> srcVal, dstVal;
  srcVal.setVal(1.0);
  dstVal.setVal(2.0);
  
  shared_ptr<vector<EBIndex<CELL> > > srcInd = getIndicies(srcBx);
  shared_ptr<vector<EBIndex<CELL> > > dstInd = getIndicies(dstBx);
  shared_ptr<vector<EBIndex<CELL> > > intInd = getIndicies(intBx);

  HostIrregData<CELL, IndexedMoments<2,2>, 1> src, dst;
  cout << "defining data" << endl;
  src.define(srcBx, srcInd);
  dst.define(dstBx, dstInd);
  src.setVal(srcVal);
  dst.setVal(dstVal);

  cout << "doing serialization dance" << endl;
  int sizeSrc = src.charsize(intBx, 0, 1);
  int sizeDst = dst.charsize(intBx, 0, 1);
  if (sizeSrc != sizeDst)
  {
    MayDay::Abort("test 3 LinearizationTest failure: dest and source have different sizes");
  }
  Vector<char> buffer(sizeSrc);
  void* buf = (void*)&(buffer[0]);
  src.linearOut(buf, intBx, 0, 1);
  dst.linearIn (buf, intBx, 0, 1);

  cout << "checking the answer" << endl;
  for(int iind = 0; iind < intInd->size(); iind++)
  {
    auto intMom = dst((*intInd)[iind], 0);
    intMom *= -1.0;
    intMom += srcVal;
    for(int imom = 0; imom < intMom.size(); imom++)
    {
      auto diff = std::abs(intMom[imom]);
      if(diff > 1.0e-6)
      {
        MayDay::Abort("bad value in test 3");
      }
    }
  }
  cout << "leaving  test3" << endl;
  return 0;
}

#include "Chombo_NamespaceFooter.H"

int main(int a_argc, char* a_argv[])
{
#ifdef CH_MPI
  MPI_Init(&a_argc, &a_argv);
  pout() << "MPI INIT called" << std::endl;
#endif
  //needs to be called after MPI_Init
  CH_TIMER_SETFILE("ebapply.time.table");
  {
    if (a_argc < 2)
    {
      cerr<< " usage " << a_argv[0] << " <input_file_name> " << endl;
      exit(0);
    }
    char* in_file = a_argv[1];
    ParmParse  pp(a_argc-2,a_argv+2,NULL,in_file);
    Chombo4::test1(a_argc, a_argv);
    Chombo4::test2(a_argc, a_argv);
  }

  pout() << "printing time table " << endl;
  CH_TIMER_REPORT();
#ifdef CH_MPI
  pout() << "about to call MPI Finalize" << std::endl;
  MPI_Finalize();
#endif
  return 0;
}
