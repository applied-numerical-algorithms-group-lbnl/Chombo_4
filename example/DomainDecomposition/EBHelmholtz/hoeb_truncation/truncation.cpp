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
#include "Chombo_EBDictionary.H"
#include "Chombo_EBChombo.H"
#include "EBMultigrid.H"

#include "DebugFunctions.H"
#include "Hoeb_ExactSolutions.H"
//this one defines HOEB_MAX_ORDER
#include "Hoeb_Utilities.H"
#include "Hoeb_LAPACKMatrix.H"
#include <iomanip>

#ifdef CH_USE_DOUBLE
Real g_tol = 1.0e-12;
#else
Real g_tol = 1.0e-6;
#endif
using std::endl;
using hoeb::LAPACKMatrix;
using Chombo4::pout;
//This is the dirty work of getting functions into the symbol table for debugging.
unsigned int
printNeighborhood(hoeb::Neighborhood<CELL>*  a_blerg)
{
  if(a_blerg != NULL)
  {
    a_blerg->poutAll();
  }
  return 0;
}
unsigned int
printCompositeStencil(hoeb::CompositeStencil*  a_blerg)
{
  if(a_blerg != NULL)
  {
    a_blerg->poutAll();
  }
  return 0;
}

//shamelessly swiped from Hans' Chombo4 test
unsigned int solveLeastSquaresTest()
{
  using Chombo4::pout;
  int size = 3;
  /* matrix A */
  Real cA[3*3] = { 3.1, 1.0 , 3.4,   
                   1.3,-6.9 , 7.2,    
                   -5.7, 5.8 ,-8.8}   ;

  /* matrix B */
  Real cB[3*3] = {-1.3, -0.1, 1.8,   
                  -1.2, -0.3, 1.9,    
                  -1.2, -0.2,  1.8}  ; 

  /* least squares solution from octave, A \ B */
  Real cC[3*3] = {1.,                1.,                1.,
                  0.941236852109999, 0.973843619313142, 0.944531745025979 ,
                  0.938461538461539, 0.953846153846154, 0.938461538461539};


  LAPACKMatrix A(size, size, cA);
  LAPACKMatrix B(size, size, cB);
  LAPACKMatrix C(size, size, cC);
  pout() << "least squares test" << endl;
  pout() << "A = " << endl;
  A.poutAll();
  pout() << "B = " << endl;
  B.poutAll();
  solveLeastSquares(A, B);
  pout() << "Answer = " << endl;
  B.poutAll();
  pout() << "correct Answer = " << endl;
  C.poutAll();
  

  // Check the answer
  C -= B;
  for(int irow = 0; irow < size; irow++)
  {
    for(int icol = 0; icol < size; icol++)
    {
      if(std::abs(C(irow,icol)) > g_tol)
      {
        pout() << "at row = " << irow << ", icol = " << icol << endl;
        pout() << "least square test returned error of "<< C(irow,icol) << endl;
        return 7;
      }
    }
  }

  return 0;
}

unsigned int inverseTest()
{
  using Chombo4::pout;
  int n =4;
  LAPACKMatrix A(n, n);
  Real val = 1;

  for(int irow = 0; irow < n; irow++)
  {
    for(int icol = 0; icol < n; icol++)
    {
      bool odd = (irow%2 == 1);
      if(odd && (irow== icol))
        A(irow, icol) = -val;
      else
        A(irow, icol) =  val;
      val += 1;
    }
  }
  LAPACKMatrix Ainv = A;
  int test = Ainv.invert();
  if(test != 0)
  {
    pout() << "we started with a singular matrix" << endl;
    return 6;
  }

  LAPACKMatrix AAinv;
  multiply(AAinv, A, Ainv);

  pout() << "inverse test" << endl;
  pout() << "A = " << endl;
  A.poutAll();
  pout() << endl;

  pout() << "Ainv = " << endl;
  Ainv.poutAll();
  pout() << endl;

  pout() << "A * Ainv = " << endl;
  AAinv.poutAll();
  pout() << endl;

  unsigned int retval = 0;
  for(int irow = 0; irow < n; irow++)
  {
    for(int icol = 0; icol < n; icol++)
    {
      Real corrval = 0;
      if(irow == icol)
      {
        corrval = 1.0;
      }
      if(std::abs(AAinv(irow, icol)-corrval) > 1.0e-3) retval += 1;
    }
  }
  return retval;
}
/**/
unsigned int inverseWithSVDTest()
{
  using Chombo4::pout;
  int n =4;
  int m = n;
  LAPACKMatrix A(m, n);
  Real val = 1;
  A.setVal(0.);

  for(int irow = 0; irow < m; irow++)
  {
    for(int icol = 0; icol < n; icol++)
    {
      bool odd = (irow%2 == 1);
      if(odd && (irow== icol))
        A(irow, icol) = -val;
      else
        A(irow, icol) =  val;
      val += 1;
    }
  }
  LAPACKMatrix Ainv = A;
  int test = Ainv.invertUsingSVD(1, 1.0e-6);
  if(test != 0)
  {
    pout() << "we started with a singular matrix" << endl;
    return 6;
  }

//  A.truncate(n, n);
//  Ainv.truncate(n, n);

  LAPACKMatrix AAinv;
  Ainv.transpose();
  multiply(AAinv, A, Ainv);

  pout() << "inverse with SVD test" << endl;
  pout() << "A = " << endl;
  A.poutAll();
  pout() << endl;

  pout() << "Ainv = " << endl;
  Ainv.poutAll();
  pout() << endl;

  pout() << "A * Ainv = " << endl;
  AAinv.poutAll();
  pout() << endl;

  unsigned int retval = 0;
  for(int irow = 0; irow < n; irow++)
  {
    for(int icol = 0; icol < n; icol++)
    {
      Real corrval = 0;
      if(irow == icol)
      {
        corrval = 1.0;
      }
      if(std::abs(AAinv(irow, icol)-corrval) > 1.0e-3) retval += 1;
    }
  }
  return retval;
}

/**/
unsigned int inverseWithLeastSquaresTest()
{
  using Chombo4::pout;
  int n =4;
  LAPACKMatrix A(n, n);
  Real val = 1;
  A.setVal(0.);
  for(int irow = 0; irow < n; irow++)
  {
    for(int icol = 0; icol < n; icol++)
    {
      bool odd = (irow%2 == 1);
      if(odd && (irow== icol))
        A(irow, icol) = -val;
      else
        A(irow, icol) =  val;
      val += 1;
    }
  }
  LAPACKMatrix Ainv = A;
  int test = Ainv.invertUsingLeastSquares();
  if(test != 0)
  {
    pout() << "we started with a singular matrix" << endl;
    return 6;
  }

  A.truncate(n, n);
  Ainv.truncate(n, n);

  LAPACKMatrix AAinv;
  multiply(AAinv, A, Ainv);

  pout() << "inverse with Least Squares test" << endl;
  pout() << "A = " << endl;
  A.poutAll();
  pout() << endl;

  pout() << "Ainv = " << endl;
  Ainv.poutAll();
  pout() << endl;

  pout() << "A * Ainv = " << endl;
  AAinv.poutAll();
  pout() << endl;

  unsigned int retval = 0;
  for(int irow = 0; irow < n; irow++)
  {
    for(int icol = 0; icol < n; icol++)
    {
      Real corrval = 0;
      if(irow == icol)
      {
        corrval = 1.0;
      }
      if(std::abs(AAinv(irow, icol)-corrval) > 1.0e-3) retval += 1;
    }
  }
  return retval;
}

/**/
unsigned int transposeTest()
{
  using Chombo4::pout;
  Real cA[3*3] = 
    {
      1,2,3,
      4,5,6,
      7,8,9
    };
  Real cAtrancorrect [3*3] = 
    {
      1,4,7,
      2,5,8,
      3,6,9
    };

  pout() << "transpose test" << endl;
  pout() << "A = " << endl;
  LAPACKMatrix A(3, 3, cA); 
  A.poutAll();
  
  LAPACKMatrix Atran = A;
  Atran.transpose();
  pout() << "A transpose = " << endl;
  Atran.poutAll();

  Real eps = 1.0e-6;
  for(int irow = 0; irow < 3; irow++)
  {
    for(int icol = 0; icol < 3; icol++)
    {
      Real compval = Atran(irow,icol);
      int offset   = Atran.offset(irow,icol);
      Real corrval = cAtrancorrect[offset];
      if(std::abs(compval - corrval) > eps)
      {
        Chombo4::MayDay::Warning("tranposition error");
        return 7;
      }
    }
  }
  return 0;
}


unsigned int testReducedRankLS()
{
  using Chombo4::pout;
  int m = 4;
  int n = 3;

  // A, the overdetermined least squares matrix, is m x n, m >= n
  // (row major order, these are the columns
  Real cA[12] = { 1, 1, 0,
                  0, 1, 0,
                  1, 0, 1,
                  0, 0, 1 };
  Real cB[] = {1, 1, 1, 0};  //last zero implicit in old code
  Real cexactRRLS[] = {1, 0, 0};
  
  LAPACKMatrix A(m, n, cA);
  LAPACKMatrix B(m, 1, cB);
  LAPACKMatrix exac(n, 1, cexactRRLS);

  pout() << "reduced rank test" << endl;
  pout() << "A = " << endl;
  A.poutAll();
  pout() << "B = " << endl;
  B.poutAll();
  solveReducedRankLS(A, B);

  pout() << "answer = " << endl;
  B.poutAll();
  pout() << "right answer =" << endl;
  exac.poutAll()        ;
  for (int i=0; i < n; ++i) 
  {
    Real calc =    B(i,0);
    Real corr = exac(i,0);
    if(std::abs(calc-corr) > g_tol)
    {
      pout() << "reduced rank ls error = "<< calc -corr << endl;
      return 7;
    }
  }
  return 0;
}


///
unsigned int testLSTransposeResid()
{
  using Chombo4::pout;
  const int m = 4;
  const int n = 3;

  // A, the underdetermined least squares matrix, is m x n, m < n
  // (row major order, these are the columns
  Real cA[12] = { 1, 1, 0, 0,
                  1, 0, 0, 0,
                  1, 0, 0, 0};
  Real cB[] = {0, 0, 1, 0};

  Real cExactLS[] = {.5, -.5, 0, 0};

  LAPACKMatrix A(m, n, cA);
  LAPACKMatrix B(m, 1, cB);
  LAPACKMatrix exactLS(m, 1, cExactLS);

  int retval = solveLSTSVDOnce(A, B);
  if(retval != 0) return retval;

  for (int i=0; i < m; ++i) 
  {
    if(std::abs(B(i, 0) - exactLS(i,0)) > g_tol)
    {
      Real residerr = B(i, 0) - exactLS(i,0);
      pout() << "lst point error  = "<< residerr << endl;
      return -7;
    }
  }

  return 0;
}

/**/
unsigned int testLSTransposeCOF()
{
  using Chombo4::pout;
  int m = 4;
  int n = 3;

  // A, the underdetermined least squares matrix, is m x n, m < n

  Real cA[12] = { 1, 1, 0, 0,
                  1, 0, 1, 0,
                  1, 0, 0, 1};
  Real cB[] = {1, 1, 1, 0};

  Real cexactRRLS[] = {.75, .25, .25, .25};

  LAPACKMatrix A(m, n, cA);
  LAPACKMatrix B(m, 1, cB);
  LAPACKMatrix exactRRLS(m, 1, cexactRRLS);

  solveLSTSVDOnce(A, B);

  for (int i=0; i < m; ++i) 
  {
    Real calc = B(i,0);
    Real corr = exactRRLS(i,0);
    if(std::abs(calc-corr) > g_tol)
    {
      pout() << "lst transpose  cof error = " << calc-corr << endl;
      return 7;
    }
  }
  return 0;
}
/**/
unsigned int lapackTest()
{

  using Chombo4::pout;
  LAPACKMatrix::s_checkConditionNumber = true;
  int icode = 0; // to be returned, sum of return values (negative or zero)

  unsigned int retval = 0;
  retval  = testLSTransposeResid();
  icode += retval;
  if(retval != 0)
  {
    pout() << "Error: testLSTTransposeResid test returned with value = " << retval << endl;
  }
  else
  {
    pout() << "LSTTransposeResid test passed" << endl;
  }

  retval = testLSTransposeCOF();
  icode += retval;
  if(retval != 0)
  {
    pout() << "Error: testLSTransposecof test returned with value = " << retval << endl;
  }
  else
  {
    pout() << "lst transpose cof test passed" << endl;
  }


  /**/
  retval = inverseTest();
  icode += retval;
  if(retval != 0)
  {
    pout() << "Error: inverse test returned with value = " << retval<< endl;
  }
  else
  {
    pout() << "inverse test passed" << endl;
  }
  /**/

  retval = 0; 
  retval = solveLeastSquaresTest();
  icode += retval;
  if(retval != 0)
  {
    pout() << "Error: leastSquares test returned with value = " << retval << endl;
  }
  else
  {
    pout() << "least squares test passed" << endl;
  }
  /**/

  retval = testReducedRankLS();
  icode += retval;
  if(retval != 0)
  {
    pout() << "Error: reduced rank ls test returned with value = " << retval << endl;
  }
  else
  {
    pout() << "reduced rank ls test passed" << endl;
  }

  retval = inverseWithLeastSquaresTest();
  icode += retval;
  if(retval != 0)
  {
    pout() << "Error: inverse with least squares test returned with value = " << retval<< endl;
  }
  else
  {
    pout() << "inverse with least squares test passed" << endl;
  }
  /**/

  /**/
  retval = inverseWithSVDTest();
  icode += retval;
  if(retval != 0)
  {
    pout() << "Error: inverse with SVD test returned with value = " << retval<< endl;
  }
  else
  {
    pout() << "inverse with SVD test passed" << endl;
  }
  /**/

  /**/
  retval = transposeTest();
  icode += retval;
  if(retval != 0)
  {
    pout() << "Error: transpose test returned with value = " << retval << endl;
  }
  else
  {
    pout() << "transpose test passed" << endl;
  }
  /**/

  if(icode == 0)
  {
    pout() << "all tests passed" << endl;
  }
  else
  {
    pout() << "not all tests passed, returned total value " << icode << endl;
  }

  return icode;
}
/****/
void
getKappaLphiHomogeneous(EBLevelBoxData<CELL, 1>                                            &  a_klp,
                        const EBLevelBoxData<CELL, 1>                                      &  a_phi,
                        const shared_ptr<hoeb::graph_distrib_t>                              &  a_graphs,
                        const Chombo4::DisjointBoxLayout                                   &  a_grids,
                        const Chombo4::Box                                                 &  a_domain,
                        const Real                                                         &  a_dx,
                        const shared_ptr<EBDictionary<HOEB_MAX_ORDER, Real, CELL, CELL> >  &  a_dictionary,
                        const shared_ptr< GeometryService<HOEB_MAX_ORDER> >                &  a_geoserv,
                        string a_stencilname)
{

  Chombo4::DataIterator dit = a_grids.dataIterator();
  //register it for every box
  for(unsigned int ibox = 0; ibox < dit.size(); ++ibox)
  {
    string stencilName;
    string ebbcName;
    vector<     EBIndex<CELL>  >          dstVoFs;
    vector<LocalStencil<CELL, Real> >     stencil;
    Proto::Box                            srcValid;
    Proto::Box                            dstValid;
    Proto::Box                            srcDomain;
    Proto::Box                            dstDomain;
    Point                                 srcGhost;
    Point                                 dstGhost;
    bool                                  needDiagonalWeights;
    
    if(a_stencilname == string("Devendran"))
    {
      bool printStuff = false;
      Real alpha = 0;
      Real beta  = 1;
      int stenRad = 3;
      int maxRad  = 4;
      string dombcarray[2*DIM];
      for(int ivec = 0; ivec < 2*DIM; ivec++)
      {
        dombcarray[ivec] = string("Neumann");
      }
      hoeb::
        getHomogeneousDharshiStencil(stencilName,        
                                     dstVoFs,            
                                     stencil,            
                                     srcValid,           
                                     dstValid,           
                                     srcDomain,          
                                     dstDomain,          
                                     needDiagonalWeights,
                                     a_geoserv,
                                     a_grids,
                                     a_domain,
                                     a_dx,
                                     ibox,
                                     dombcarray,
                                     ebbcName,
                                     alpha, beta,
                                     stenRad, maxRad, printStuff);
    }
    else if(a_stencilname == string("Schwartz"))
    {
      hoeb::
        schwartzLaplStencil(stencilName,        
                            ebbcName,           
                            dstVoFs,            
                            stencil,            
                            srcValid,           
                            dstValid,           
                            srcDomain,          
                            dstDomain,          
                            srcGhost,           
                            dstGhost,
                            needDiagonalWeights,
                            a_geoserv,
                            a_grids,
                            a_domain,
                            a_dx,
                            ibox);
    }

    ///registering stencil
    a_dictionary->registerStencil(stencilName,        
                                  ebbcName,           
                                  dstVoFs,            
                                  stencil,            
                                  srcValid,           
                                  dstValid,           
                                  srcDomain,          
                                  dstDomain,          
                                  srcGhost,           
                                  dstGhost,
                                  needDiagonalWeights,
                                  ibox);
    //now apply it
    for(unsigned int ibox = 0; ibox < dit.size(); ++ibox)
    {
      auto      & lphfab = a_klp[dit[ibox]];
      const auto& phifab = a_phi[dit[ibox]];
      auto stencil = a_dictionary->getEBStencil(stencilName, ebbcName, a_domain, a_domain, ibox);
      //set resc = Ave(resf) (true is initToZero)
      stencil->apply(lphfab, phifab,  true, 1.0);
    }
  }
}
typedef CH4_Data_Choreography::DistributedData<EBGraph>   graph_distrib_t;
/*******/ 
void
getKappaLphi(EBLevelBoxData<CELL, 1>                                            &  a_klp,
             const EBLevelBoxData<CELL, 1>                                      &  a_phi,
             const shared_ptr<graph_distrib_t>                                  &  a_graphs,
             const Chombo4::DisjointBoxLayout                                   &  a_grids,
             const Chombo4::Box                                                 &  a_domain,
             const Real                                                         &  a_dx,
             const shared_ptr<EBDictionary<HOEB_MAX_ORDER, Real, CELL, CELL> >  &  a_dictionary,
             const shared_ptr< GeometryService<HOEB_MAX_ORDER> >                &  a_geoserv,
             string                                                                a_stencilname)
{

  typedef GraphConstructorFactory<EBHostData<CELL, Real, 1> > hostfactorycell_t;
  typedef CH4_Data_Choreography::DistributedData<EBHostData<CELL, Real, 1> > cell_distrib_t;
  
  IntVect dataGhostIV = a_klp.ghostVect();
  cell_distrib_t  hostklp(a_grids, dataGhostIV, hostfactorycell_t(a_graphs));


  Chombo4::DataIterator dit = a_grids.dataIterator();
  for(unsigned int ibox = 0; ibox < dit.size(); ++ibox)
  {
    auto chomboBox =a_grids[dit[ibox]];
    Proto::Box bx = getProtoBox(chomboBox);
    const auto & graph= (*a_graphs)[dit[ibox]];
    for(auto bit = bx.begin(); bit != bx.end(); ++bit)
    {
      auto pt = *bit;
      auto vofs = graph.getVoFs(pt);
      for(int ivof = 0; ivof < vofs.size(); ivof++)
      {
        const auto& vof = vofs[ivof];
        Real vofval = hoeb::getKappaLphiVoF(vof,
                                            dit[ibox],
                                            a_phi,        
                                            a_graphs,     
                                            a_grids,      
                                            a_domain,     
                                            a_dx,         
                                            a_dictionary, 
                                            a_geoserv,    
                                            a_stencilname);
        auto& klpfab = hostklp[dit[ibox]];
        klpfab(vof, 0) = vofval;
      }
    }
  }
  EBLevelBoxData<CELL,1>::copyToDevice(a_klp, hostklp); 
}
/*******/ 
PROTO_KERNEL_START 
void subtractionTractionF(Var<Real, 1>    a_error,
                          Var<Real, 1>    a_klpFtoC,
                          Var<Real, 1>    a_klpCoar)
{
  Real ftoc  = a_klpFtoC(0);
  Real coar  = a_klpCoar(0);
  Real errr  = ftoc - coar;
//  pout() << ftoc << " - " << coar << " = " << errr << endl;
  a_error(0) = errr;
}
PROTO_KERNEL_END(subtractionTractionF, subtractionTraction)
/****/
void
getKLPhiError(EBLevelBoxData<CELL,   1>                                           &  a_errCoar, 
              const shared_ptr<graph_distrib_t    >                               &  a_graphsFine,
              const Chombo4::DisjointBoxLayout                                    &  a_gridsFine,
              const Chombo4::Box                                                  &  a_domFine,
              const Real                                                          &  a_dxFine,
              const shared_ptr<graph_distrib_t    >                               &  a_graphsCoar,
              const Chombo4::DisjointBoxLayout                                    &  a_gridsCoar,
              const Chombo4::Box                                                  &  a_domCoar,
              const Real                                                          &  a_dxCoar,
              const shared_ptr<EBDictionary<HOEB_MAX_ORDER, Real, CELL, CELL> >   &  a_dictionary,
              const shared_ptr< GeometryService<HOEB_MAX_ORDER> >                 &  a_geoserv,
              string a_stencilname)
{
  ParmParse pp;
  int nghost = 6;
  IntVect dataGhostIV = nghost*IntVect::Unit;
  EBLevelBoxData<CELL,   1>  phiFine(a_gridsFine, dataGhostIV, a_graphsFine);
  EBLevelBoxData<CELL,   1>  phiCoar(a_gridsCoar, dataGhostIV, a_graphsCoar);
  EBLevelBoxData<CELL,   1>  klpFine(a_gridsFine, dataGhostIV, a_graphsFine);
  EBLevelBoxData<CELL,   1>  klpCoar(a_gridsCoar, dataGhostIV, a_graphsCoar);
  EBLevelBoxData<CELL,   1>  klpFtoC(a_gridsCoar, dataGhostIV, a_graphsCoar);
  
  hoeb::fillPhi(phiFine, a_graphsFine, a_gridsFine, a_domFine, a_dxFine, a_geoserv);
  hoeb::fillPhi(phiCoar, a_graphsCoar, a_gridsCoar, a_domCoar, a_dxCoar, a_geoserv);

  getKappaLphi(klpFine, phiFine, a_graphsFine, a_gridsFine, a_domFine, a_dxFine, a_dictionary, a_geoserv, a_stencilname);
  getKappaLphi(klpCoar, phiCoar, a_graphsCoar, a_gridsCoar, a_domCoar, a_dxCoar, a_dictionary, a_geoserv, a_stencilname);

  hoeb::restrictKappaLphi(klpFtoC, klpFine,
                          a_graphsFine, a_gridsFine, a_domFine, a_dxFine,                    
                          a_graphsCoar, a_gridsCoar, a_domCoar, a_dxCoar,
                          a_dictionary, a_geoserv);

  //error = Ave(klphifine) - klphicoar 
  Chombo4::DataIterator dit = a_gridsCoar.dataIterator();
  for(unsigned int ibox = 0; ibox < dit.size(); ibox++)
  {
    auto& ftocfab =   klpFtoC[dit[ibox]];
    auto& coarfab =   klpCoar[dit[ibox]];
    auto& errfab  = a_errCoar[dit[ibox]];
    auto  inputbx = ftocfab.inputBox();
    auto  validbx = (*a_graphsCoar)[dit[ibox]].validBox();
    ebforall(inputbx, subtractionTraction, validbx, errfab, ftocfab, coarfab);
  }
}
/****/
void
getICError(EBLevelBoxData<CELL,   1>                                           &  a_errCoar, 
           const shared_ptr<graph_distrib_t>                                   &  a_graphsFine,
           const Chombo4::DisjointBoxLayout                                    &  a_gridsFine,
           const Chombo4::Box                                                  &  a_domFine,
           const Real                                                          &  a_dxFine,
           const shared_ptr<graph_distrib_t>                                   &  a_graphsCoar,
           const Chombo4::DisjointBoxLayout                                    &  a_gridsCoar,
           const Chombo4::Box                                                  &  a_domCoar,
           const Real                                                          &  a_dxCoar,
           const shared_ptr<EBDictionary<HOEB_MAX_ORDER, Real, CELL, CELL> >   &  a_dictionary,
           const shared_ptr< GeometryService<HOEB_MAX_ORDER> >                 &  a_geoserv)
{
  ParmParse pp;
  int nghost = 6;

  IntVect dataGhostIV =   nghost*IntVect::Unit;
  EBLevelBoxData<CELL,   1>  phiFine(a_gridsFine, dataGhostIV, a_graphsFine);
  EBLevelBoxData<CELL,   1>  phiCoar(a_gridsCoar, dataGhostIV, a_graphsCoar);
  EBLevelBoxData<CELL,   1>  phiFtoC(a_gridsCoar, dataGhostIV, a_graphsCoar);
  
  hoeb::fillPhi(phiFine, a_graphsFine, a_gridsFine, a_domFine, a_dxFine, a_geoserv);
  hoeb::fillPhi(phiCoar, a_graphsCoar, a_gridsCoar, a_domCoar, a_dxCoar, a_geoserv);

  hoeb::restrictPhi(phiFtoC, phiFine,
                    a_graphsFine, a_gridsFine, a_domFine, a_dxFine,                    
                    a_graphsCoar, a_gridsCoar, a_domCoar, a_dxCoar,
                    a_dictionary, a_geoserv);

  //error = Ave(phifine) - phicoar 
  Chombo4::DataIterator dit = a_gridsCoar.dataIterator();
  for(unsigned int ibox = 0; ibox < dit.size(); ibox++)
  {
    auto& ftocfab =   phiFtoC[dit[ibox]];
    auto& coarfab =   phiCoar[dit[ibox]];
    auto& errfab  = a_errCoar[dit[ibox]];
    auto  inputbx = ftocfab.inputBox();
    auto  validbx = (*a_graphsCoar)[dit[ibox]].validBox();
    ebforall(inputbx, subtractionTraction, validbx, errfab, ftocfab, coarfab);
  }
}

/****/
unsigned int
runTruncationErrorTest(string a_stencilname)
{
  using Chombo4::pout;
  Real coveredval = -1;
  int nx      = 32;
  int maxGrid = 32;
    
  ParmParse pp;

  pp.get("nx"        , nx);
  pp.get("maxGrid"  , maxGrid);
  pp.get("coveredval", coveredval);         


  pout() << "nx"        << " = " <<  nx         << endl;
  pout() << "max_grid"  << " = " <<  maxGrid    << endl;
  pout() << "coveredval"<< " = " <<  coveredval << endl;         
  

  IntVect domLo = IntVect::Zero;
  IntVect domHi  = (nx - 1)*IntVect::Unit;

  Chombo4::ProblemDomain domain(domLo, domHi);

  vector<Chombo4::DisjointBoxLayout> vecgrids;
  pout() << "making grids" << endl;
  GeometryService<2>::generateGrids(vecgrids, domain.domainBox(), maxGrid);

  Chombo4::DisjointBoxLayout gridsFine = vecgrids[0];
  Chombo4::DisjointBoxLayout gridsMedi = vecgrids[1];
  Chombo4::DisjointBoxLayout gridsCoar = vecgrids[2];
  IntVect dataGhostIV =   6*IntVect::Unit;
  Point   dataGhostPt = ProtoCh::getPoint(dataGhostIV); 
  int geomGhost = 6;
  RealVect origin = RealVect::Zero();
  
  shared_ptr<BaseIF>    impfunc = hoeb::getImplicitFunction();
  pout() << "defining geometry" << endl;
  
  //Real dx = 1.0/(Real(nx));
  Real dx = 32.0/(Real(nx));
  vector<Chombo4::Box>    vecdomain(vecgrids.size(), domain.domainBox());
  vector<Real>   vecdx    (vecgrids.size(), dx);
  for(int ilev = 1; ilev < vecgrids.size(); ilev++)
  {
    vecdomain[ilev] = coarsen(vecdomain[ilev-1], 2);
    vecdx    [ilev] =           2*vecdx[ilev-1];
  }
  
  Real dxFine = vecdx[0];
  Real dxMedi = vecdx[1];
  Real dxCoar = vecdx[2];
  GeometryService<HOEB_MAX_ORDER>* geomptr =
    new GeometryService<HOEB_MAX_ORDER>
    (impfunc, origin, dxFine, domain.domainBox(), vecgrids, geomGhost);
  
  shared_ptr< GeometryService<HOEB_MAX_ORDER> >  geoserv(geomptr);

  pout() << "making dictionary" << endl;

  Chombo4::Box domFine = vecdomain[0];
  Chombo4::Box domMedi = vecdomain[1];
  Chombo4::Box domCoar = vecdomain[2];

  shared_ptr<graph_distrib_t> graphsFine = geoserv->getGraphs(domFine);
  shared_ptr<graph_distrib_t> graphsMedi = geoserv->getGraphs(domMedi);
  shared_ptr<graph_distrib_t> graphsCoar = geoserv->getGraphs(domCoar);

  shared_ptr<EBDictionary<HOEB_MAX_ORDER, Real, CELL, CELL> >  dictionary
    (new     EBDictionary<HOEB_MAX_ORDER, Real, CELL, CELL>
     (geoserv, vecgrids, vecdomain, vecdx, dataGhostPt));

  EBLevelBoxData<CELL,   1>  errMedi(gridsMedi, dataGhostIV, graphsMedi);
  EBLevelBoxData<CELL,   1>  errCoar(gridsCoar, dataGhostIV, graphsCoar);
  

  getKLPhiError(errMedi, 
                graphsFine, gridsFine, domFine, dxFine,
                graphsMedi, gridsMedi, domMedi, dxMedi,
                dictionary, geoserv, a_stencilname);

  getKLPhiError(errCoar, 
                graphsMedi, gridsMedi, domMedi, dxMedi,
                graphsCoar, gridsCoar, domCoar, dxCoar,
                dictionary, geoserv, a_stencilname);


  pout() << "writing to file" << endl;
  string fileCoar = a_stencilname + string("errCoar.hdf5");
  string fileMedi = a_stencilname + string("errMedi.hdf5");
  
  errCoar.writeToFileHDF5(fileCoar, 0.0);
  errMedi.writeToFileHDF5(fileMedi, 0.0);
  //Norm!
  Real normMedi = errMedi.maxNorm(0);
  Real normCoar = errCoar.maxNorm(0);

  Real tol = 1.0e-16;
  Real order = 0;
  if((normCoar > tol) && (normMedi > tol))
  {
    order = log(normCoar/normMedi)/log(2.0);
  }
  pout() << "||klphi errMedi||_max = " << normMedi << std::endl;
  pout() << "||klphi errCoar||_max = " << normCoar << std::endl;
  pout() << a_stencilname << " Richardson truncation error order for kappa(L(phi))= " << order << std::endl;
  return 0;
}

/****/
unsigned int
runInitialConditionTest()
{
  using Chombo4::pout;
  Real coveredval = -1;
  int nx      = 32;
  int maxGrid = 32;
    
  ParmParse pp;

  pp.get("nx"        , nx);
  pp.get("maxGrid"  , maxGrid);
  pp.get("coveredval", coveredval);         


  pout() << "nx"        << " = " <<  nx         << endl;
  pout() << "max_grid"  << " = " <<  maxGrid    << endl;
  pout() << "coveredval"<< " = " <<  coveredval << endl;         
  

  IntVect domLo = IntVect::Zero;
  IntVect domHi  = (nx - 1)*IntVect::Unit;

  Chombo4::ProblemDomain domain(domLo, domHi);

  vector<Chombo4::DisjointBoxLayout> vecgrids;
  pout() << "making grids" << endl;
  GeometryService<2>::generateGrids(vecgrids, domain.domainBox(), maxGrid);

  Chombo4::DisjointBoxLayout gridsFine = vecgrids[0];
  Chombo4::DisjointBoxLayout gridsMedi = vecgrids[1];
  Chombo4::DisjointBoxLayout gridsCoar = vecgrids[2];
  IntVect dataGhostIV =   6*IntVect::Unit;
  Point   dataGhostPt = ProtoCh::getPoint(dataGhostIV); 
  int geomGhost = 6;
  RealVect origin = RealVect::Zero();
  
  shared_ptr<BaseIF>    impfunc = hoeb::getImplicitFunction();
  pout() << "defining geometry" << endl;
  Real dx = 1.0/(Real(nx));
  vector<Chombo4::Box>    vecdomain(vecgrids.size(), domain.domainBox());
  vector<Real>   vecdx    (vecgrids.size(), dx);
  for(int ilev = 1; ilev < vecgrids.size(); ilev++)
  {
    vecdomain[ilev] = coarsen(vecdomain[ilev-1], 2);
    vecdx    [ilev] =           2*vecdx[ilev-1];
  }
  
  Real dxFine = vecdx[0];
  Real dxMedi = vecdx[1];
  Real dxCoar = vecdx[2];
  GeometryService<HOEB_MAX_ORDER>* geomptr =
    new GeometryService<HOEB_MAX_ORDER>
    (impfunc, origin, dxFine, domain.domainBox(), vecgrids, geomGhost);
  
  shared_ptr< GeometryService<HOEB_MAX_ORDER> >  geoserv(geomptr);

  pout() << "making dictionary" << endl;

  Chombo4::Box domFine = vecdomain[0];
  Chombo4::Box domMedi = vecdomain[1];
  Chombo4::Box domCoar = vecdomain[2];

  shared_ptr<graph_distrib_t> graphsFine = geoserv->getGraphs(domFine);
  shared_ptr<graph_distrib_t> graphsMedi = geoserv->getGraphs(domMedi);
  shared_ptr<graph_distrib_t> graphsCoar = geoserv->getGraphs(domCoar);

  
  shared_ptr<EBDictionary<HOEB_MAX_ORDER, Real, CELL, CELL> >  dictionary
    (new     EBDictionary<HOEB_MAX_ORDER, Real, CELL, CELL>
     (geoserv, vecgrids, vecdomain, vecdx, dataGhostPt));

  EBLevelBoxData<CELL,   1>  errMedi(gridsMedi, dataGhostIV, graphsMedi);
  EBLevelBoxData<CELL,   1>  errCoar(gridsCoar, dataGhostIV, graphsCoar);
  

  getICError(errMedi, 
             graphsFine, gridsFine, domFine, dxFine,
             graphsMedi, gridsMedi, domMedi, dxMedi,
             dictionary, geoserv);

  getICError(errCoar, 
             graphsMedi, gridsMedi, domMedi, dxMedi,
             graphsCoar, gridsCoar, domCoar, dxCoar,
             dictionary, geoserv);


  //Norm!
  Real normMedi = errMedi.maxNorm(0);
  Real normCoar = errCoar.maxNorm(0);

  Real tol = 1.0e-16;
  Real order = 0;
  if((normCoar > tol) && (normMedi > tol))
  {
    order = log(normCoar/normMedi)/log(2.0);
  }
  pout() << "writing to file" << endl;
  errCoar.writeToFileHDF5(string("IC.errCoar.hdf5"), 0.0);
  errMedi.writeToFileHDF5(string("IC.errMedi.hdf5"), 0.0);
  
  pout() << "||IC errMedi||_max = " << normMedi << std::endl;
  pout() << "||IC errCoar||_max = " << normCoar << std::endl;
  pout() << "Richardson truncation error order for IC= " << order << std::endl;

  return 0;
}

/****/

/****/
template <CENTERING cent>
void
getDevendranFlux(EBLevelBoxData<cent,   1>                                           &  a_deviflux, 
                 const shared_ptr<graph_distrib_t    >                               &  a_graphs,
                 const Chombo4::DisjointBoxLayout                                    &  a_grids,
                 const Chombo4::Box                                                  &  a_dom,
                 const Real                                                          &  a_dx,
                 const shared_ptr< GeometryService<HOEB_MAX_ORDER> >                 &  a_geoserv)
{
  using Chombo4::pout;
  //fill on the host and transfer over to the device rather than go  through the pain and
  //suffering necessary to squeeze this into forall
  IntVect ghost = a_deviflux.ghostVect();
  typedef GraphConstructorFactory<EBHostData<cent, Real, 1> > hostfactoryface_t;
  typedef GraphConstructorFactory<EBHostData<CELL, Real, 1> > hostfactorycell_t;
  typedef CH4_Data_Choreography::DistributedData<EBHostData<CELL, Real, 1> > cell_distrib_t;
  typedef CH4_Data_Choreography::DistributedData<EBHostData<cent, Real, 1> > cent_distrib_t;
  cent_distrib_t  hostflux(a_grids, 1, ghost, hostfactoryface_t(a_graphs));
  int nghost = 6;
  IntVect dataGhostIV = nghost*IntVect::Unit;
  EBLevelBoxData<CELL,   1>  deviphi(a_grids, nghost*IntVect::Unit, a_graphs);
  hoeb::fillPhi(deviphi, a_graphs, a_grids, a_dom, a_dx, a_geoserv);
  
  cell_distrib_t   hostphi(a_grids, 1, ghost, hostfactorycell_t(a_graphs));
  EBLevelBoxData<CELL,   1>::copyToHost(hostphi, deviphi);
  

  int facedir = 0;
  if(cent == XFACE)
  {
    facedir = 0;
  }
  else if(cent == YFACE)
  {
    facedir = 1;
  }
  else if(cent == ZFACE)
  {
    facedir = 2;
  }
  else
  {
    PROTO_ASSERT(false, "bad face type, no donut");
  }
  string dombcname[2*DIM];
  string ebbcname("Dirichlet");
  for(int ivec = 0; ivec < 2*DIM; ivec++)
  {
    dombcname[ivec] = string("Dirichlet");
  }
  
  Chombo4::DataIterator dit = a_grids.dataIterator();
  for(unsigned int ibox = 0; ibox < dit.size(); ibox++)
  {
    hostflux[dit[ibox]].setVal(0.);  //fill all those nasty  covered values with zero
    auto chomboBox =a_grids[dit[ibox]];
    Proto::Box bx = getProtoBox(chomboBox);
    const auto & graph= (*a_graphs)[dit[ibox]];
    for(auto bit = bx.begin(); bit != bx.end(); ++bit)
    {
      auto pt = *bit;
      auto vofs = graph.getVoFs(pt);
      for(int ivof = 0; ivof < vofs.size(); ivof++)
      {
        const auto& vof = vofs[ivof];
        for(SideIterator sit; sit.ok(); ++sit)
        {
          vector<EBIndex<cent> > faces = Proto::getFaces<cent>(vof, sit(), graph);
          for(int iface  = 0; iface < faces.size(); ++iface)
          {
            const auto& face = faces[iface];
            bool dividebyarea = true; //we want the average flux here
            Real fluxpt =
              hoeb::getDevendranFluxFace<cent>(hostphi[dit[ibox]], (*a_graphs)[dit[ibox]],
                                               face, vof,  dombcname, ebbcname, a_geoserv,
                                               a_dom, a_dx, facedir, sit(),
                                               dit[ibox], dividebyarea, false);
            hostflux[dit[ibox]](face, 0) = fluxpt;
          }
        } 
      } 
    }
  }
  EBLevelBoxData<cent, 1>::copyToDevice(a_deviflux, hostflux);
}
/****/
template <CENTERING cent>
void
getDevendranFluxError(EBLevelBoxData<cent,   1>                                           &  a_errCoar, 
                      const shared_ptr<graph_distrib_t>                                   &  a_graphsFine,
                      const Chombo4::DisjointBoxLayout                                    &  a_gridsFine,
                      const Chombo4::Box                                                  &  a_domFine,
                      const Real                                                          &  a_dxFine,
                      const shared_ptr<graph_distrib_t>                               &  a_graphsCoar,
                      const Chombo4::DisjointBoxLayout                                    &  a_gridsCoar,
                      const Chombo4::Box                                                  &  a_domCoar,
                      const Real                                                          &  a_dxCoar,
                      const shared_ptr<EBDictionary<HOEB_MAX_ORDER, Real, cent, cent> >   &  a_dictionary,
                      const shared_ptr< GeometryService<HOEB_MAX_ORDER> >                 &  a_geoserv,
                      bool a_finestLevel)
{
  ParmParse pp;
  int nghost = 6;

  IntVect dataGhostIV =   nghost*IntVect::Unit;
  EBLevelBoxData<cent,   1>  fluxFine(a_gridsFine, dataGhostIV, a_graphsFine);
  EBLevelBoxData<cent,   1>  fluxCoar(a_gridsCoar, dataGhostIV, a_graphsCoar);
  EBLevelBoxData<cent,   1>  fluxFtoC(a_gridsCoar, dataGhostIV, a_graphsCoar);
  
  getDevendranFlux(fluxFine, a_graphsFine, a_gridsFine, a_domFine, a_dxFine, a_geoserv);
  getDevendranFlux(fluxCoar, a_graphsCoar, a_gridsCoar, a_domCoar, a_dxCoar, a_geoserv);

  hoeb::restrictFlux<cent>(fluxFtoC, fluxFine,
                    a_graphsFine, a_gridsFine, a_domFine, a_dxFine,                    
                    a_graphsCoar, a_gridsCoar, a_domCoar, a_dxCoar,
                    a_dictionary, a_geoserv);

  //error = Ave(phifine) - phicoar 
  Chombo4::DataIterator dit = a_gridsCoar.dataIterator();
  for(unsigned int ibox = 0; ibox < dit.size(); ibox++)
  {
    auto& ftocfab =   fluxFtoC[dit[ibox]];
    auto& coarfab =   fluxCoar[dit[ibox]];
    auto& errfab  = a_errCoar[dit[ibox]];
    auto  inputbx = ftocfab.inputBox();
    auto  validbx = (*a_graphsCoar)[dit[ibox]].validBox();
    ebforall(inputbx, subtractionTraction, validbx, errfab, ftocfab, coarfab);
  }

  if(a_finestLevel)
  {
    if(cent == XFACE)
    {
      fluxFine.writeToFileHDF5(string("fluxFine_x.hdf5"), 0.0);
      fluxCoar.writeToFileHDF5(string("fluxCoar_x.hdf5"), 0.0);
      fluxFtoC.writeToFileHDF5(string("fluxFtoC_x.hdf5"), 0.0);
    }
    else if(cent == YFACE)
    {
      fluxFine.writeToFileHDF5(string("fluxFine_y.hdf5"), 0.0);
      fluxCoar.writeToFileHDF5(string("fluxCoar_y.hdf5"), 0.0);
      fluxFtoC.writeToFileHDF5(string("fluxFtoC_y.hdf5"), 0.0);
    }
  }
}
/****/
template<CENTERING cent>
unsigned int
devendranFluxTruncation(string a_prefix)
{
#if 0 
  using Chombo4::pout;
  Real coveredval = -1;
  int nx      = 32;
  int maxGrid = 32;
    
  ParmParse pp;

  pp.get("nx"        , nx);
  pp.get("maxGrid"  , maxGrid);
  pp.get("coveredval", coveredval);         


  pout() << "nx"        << " = " <<  nx         << endl;
  pout() << "max_grid"  << " = " <<  maxGrid    << endl;
  pout() << "coveredval"<< " = " <<  coveredval << endl;         
  

  IntVect domLo = IntVect::Zero;
  IntVect domHi  = (nx - 1)*IntVect::Unit;

  Chombo4::ProblemDomain domain(domLo, domHi);

  vector<Chombo4::DisjointBoxLayout> vecgrids;
  pout() << "making grids" << endl;
  GeometryService<2>::generateGrids(vecgrids, domain.domainBox(), maxGrid);

  Chombo4::DisjointBoxLayout gridsFine = vecgrids[0];
  Chombo4::DisjointBoxLayout gridsMedi = vecgrids[1];
  Chombo4::DisjointBoxLayout gridsCoar = vecgrids[2];
  IntVect dataGhostIV =   6*IntVect::Unit;
  Point   dataGhostPt = ProtoCh::getPoint(dataGhostIV); 
  int geomGhost = 6;
  RealVect origin = RealVect::Zero();
  
  shared_ptr<BaseIF>    impfunc = hoeb::getImplicitFunction();
  pout() << "defining geometry" << endl;
  Real dx = 1.0/(Real(nx));
  vector<Chombo4::Box>    vecdomain(vecgrids.size(), domain.domainBox());
  vector<Real>   vecdx    (vecgrids.size(), dx);
  for(int ilev = 1; ilev < vecgrids.size(); ilev++)
  {
    vecdomain[ilev] = coarsen(vecdomain[ilev-1], 2);
    vecdx    [ilev] =           2*vecdx[ilev-1];
  }
  
  Real dxFine = vecdx[0];
  Real dxMedi = vecdx[1];
  Real dxCoar = vecdx[2];
  GeometryService<HOEB_MAX_ORDER>* geomptr =
    new GeometryService<HOEB_MAX_ORDER>
    (impfunc, origin, dxFine, domain.domainBox(), vecgrids, geomGhost);
  
  shared_ptr< GeometryService<HOEB_MAX_ORDER> >  geoserv(geomptr);

  pout() << "making dictionary" << endl;

  Chombo4::Box domFine = vecdomain[0];
  Chombo4::Box domMedi = vecdomain[1];
  Chombo4::Box domCoar = vecdomain[2];

  shared_ptr<graph_distrib_t> graphsFine = geoserv->getGraphs(domFine);
  shared_ptr<graph_distrib_t> graphsMedi = geoserv->getGraphs(domMedi);
  shared_ptr<graph_distrib_t> graphsCoar = geoserv->getGraphs(domCoar);

  
  shared_ptr<EBDictionary<HOEB_MAX_ORDER, Real, cent, cent>  >  dictionary
    (new     EBDictionary<HOEB_MAX_ORDER, Real, cent, cent>
     (geoserv, vecgrids, vecdomain, vecdx, dataGhostPt));

  EBLevelBoxData<cent,   1>  errMedi(gridsMedi, dataGhostIV, graphsMedi);
  EBLevelBoxData<cent,   1>  errCoar(gridsCoar, dataGhostIV, graphsCoar);
  

  getDevendranFluxError<cent>(errMedi, 
                              graphsFine, gridsFine, domFine, dxFine,
                              graphsMedi, gridsMedi, domMedi, dxMedi,
                              dictionary, geoserv, true);

  getDevendranFluxError<cent>(errCoar, 
                              graphsMedi, gridsMedi, domMedi, dxMedi,
                              graphsCoar, gridsCoar, domCoar, dxCoar,
                              dictionary, geoserv, false);


  //Norm!
  Real normMedi = errMedi.maxNorm(0);
  Real normCoar = errCoar.maxNorm(0);

  Real tol = 1.0e-16;
  Real order = 0;
  if((normCoar > tol) && (normMedi > tol))
  {
    order = log(normCoar/normMedi)/log(2.0);
  }
  pout() << "writing to file" << endl;
  string fileCoar = a_prefix + string(".errCoar.hdf5");
  string fileMedi = a_prefix + string(".errMedi.hdf5");
  errCoar.writeToFileHDF5(fileCoar, 0.0);
  errMedi.writeToFileHDF5(fileMedi, 0.0);
  
  pout() << "|| " << a_prefix<< "  errMedi ||_max = " << normMedi << std::endl;
  pout() << "|| " << a_prefix<< "  errCoar ||_max = " << normCoar << std::endl;
  pout() << "Richardson truncation error order for Devendran flux = " << order << std::endl;

#endif
  
  return 0;
}
/****/
unsigned int
runTests()
{

  using Chombo4::pout;
  bool runIC, runDCT, runSCT, runLAP, runDevFluxTE;
  ParmParse pp;
  pp.get("do_lapack_test"           , runLAP);
  pp.get("do_initial_condition_test", runIC);
  pp.get("do_devendran_test"        , runDCT);
  pp.get("do_schwartz_test"         , runSCT);
  pp.get("do_devendran_flux_test"   , runDevFluxTE);
  unsigned int errcode = 0;
  if(runDevFluxTE)
  {
    unsigned int errcode1 = devendranFluxTruncation<XFACE>(string("devfluxX")); hoeb::checkError(errcode1, "devendran_flux_truncation_error_x");
    errcode += 10*errcode1;
    unsigned int errcode0 = devendranFluxTruncation<YFACE>(string("devfluxY")); hoeb::checkError(errcode0, "devendran_flux_truncation_error_y");
    errcode += errcode0;
#if DIM==3
    unsigned int errcode2 = devendranFluxTruncation<ZFACE>(string("devfluxZ")); hoeb::checkError(errcode2, "devendran_flux_truncation_error_z");
    errcode += 100*errcode2;
#endif    
  }
  if(runLAP)
  {
    unsigned int errcode3 = lapackTest()                 ; hoeb::checkError(errcode3, "lapack");
    errcode += 1000*errcode3;
  }
  if(runIC)
  {
    unsigned int errcode4 = runInitialConditionTest()     ; hoeb::checkError(errcode4, "initial_condition");
    errcode += 10000*errcode4;
  }
  if(runSCT)
  {
    unsigned int errcode5 = runTruncationErrorTest(string("Schwartz")) ; hoeb::checkError(errcode5, "Schwartz_stencil");
    errcode += 100000*errcode5;
  }
  if(runDCT)
  {
    unsigned int errcode6 = runTruncationErrorTest(string("Devendran")); hoeb::checkError(errcode6, "Devendran_stencil");
    errcode += 1000000*errcode6;
  }
  return errcode;
}
/****/

int main(int a_argc, char* a_argv[])
{
#ifdef CH_USE_PETSC  
  //because of some kind of solipsistic madness, PetscInitialize calls MPI_INIT
  PetscInt ierr = PetscInitialize(&a_argc, &a_argv, "./.petscrc",PETSC_NULL); CHKERRQ(ierr);
#else  
#ifdef CH_MPI
  MPI_Init(&a_argc, &a_argv);
  pout() << "MPI INIT called" << std::endl;
#endif
#endif
  using Chombo4::pout;
  //needs to be called after MPI_Init
  CH_TIMER_SETFILE("trunc.time.table");
  {
    printNeighborhood(NULL);
    printCompositeStencil(NULL);
    
    if (a_argc < 2)
    {
      cerr<< " usage " << a_argv[0] << " <input_file_name> " << endl;
      exit(0);
    }
    char* in_file = a_argv[1];
    ParmParse  pp(a_argc-2,a_argv+2,NULL,in_file);
    runTests();
  }

  pout() << "printing time table " << endl;
  CH_TIMER_REPORT();
#ifdef CH_MPI
#ifdef CH_USE_PETSC
  pout() << "about to call petsc Finalize" << std::endl;
  PetscFinalize();
#else  
  pout() << "about to call MPI Finalize" << std::endl;
  MPI_Finalize();
#endif
#endif
  return 0;
}
