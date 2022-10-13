#ifndef __hoeb_utilities__
#define __hoeb_utilities__
#include "Proto.H"
#include "Chombo_EBChombo.H"
#include "Chombo_EBLevelBoxData.H"
#include "Chombo_GeometryService.H"
#include "Chombo_EBDictionary.H"
#include "Hoeb_LAPACKMatrix.H"
#include "Hoeb_Neighborhood.H"
#include "Hoeb_CompositeStencil.H"
#define HOEB_MAX_ORDER 4

namespace hoeb
{
  typedef CH4_Data_Choreography::DistributedData<EBGraph>   graph_distrib_t;
  typedef GraphConstructorFactory<EBHostData<CELL, Real, 1> > hostfactorycell_t;
  typedef CH4_Data_Choreography::DistributedData<EBHostData<CELL, Real, 1> > cell_distrib_t;
  
  template <CENTERING cent>
  Real
  getNeumannBoundaryFlux(const shared_ptr<hoeb::BaseExactSolution<HOEB_MAX_ORDER> >& a_exactSolution,
                         const EBHostData<CELL, Real, 1>                     & a_phi,
                         const EBGraph                                       & a_graph,
                         const EBIndex<cent>                                 & a_face,
                         const EBIndex<CELL>                                 & a_startingVoF,
                         const std::string                                     a_dombcname[2*DIM],
                         const std::string                                   & a_ebbcname,
                         const shared_ptr< GeometryService<HOEB_MAX_ORDER> > & a_geoserv,
                         const Chombo4::Box                                  & a_domain,
                         Real a_dx, int a_facedir, Side::LoHiSide a_vofside,
                         const Chombo4::DataIndex& a_dit, bool a_divideByArea) 
  {
    Chombo4::MayDay::Error("we have a problem");
    return 0;
  }

  
/***/
  template < >
  inline Real
  getNeumannBoundaryFlux(const shared_ptr<hoeb::BaseExactSolution<HOEB_MAX_ORDER> >& a_exactSolution,
                         const EBHostData<CELL, Real, 1>                     & a_phi,
                         const EBGraph                                       & a_graph,
                         const EBIndex<XFACE>                                & a_face,
                         const EBIndex<CELL>                                 & a_startingVoF,
                         const std::string                                     a_dombcname[2*DIM],
                         const std::string                                   & a_ebbcname,
                         const shared_ptr< GeometryService<HOEB_MAX_ORDER> > & a_geoserv,
                         const Chombo4::Box                                  & a_domain,
                         Real a_dx, int a_facedir,
                         Side::LoHiSide a_vofside, const Chombo4::DataIndex& a_dit,
                         bool a_divideByArea) 
  {
    const auto& facld = a_geoserv->getXFaceData(a_domain);
    const auto& dat   = (*facld)[a_dit];
    Real gradval = a_exactSolution->
      evaluateGradDMO(a_graph, a_dx, dat, a_face, 0, a_divideByArea);
    return gradval;
    
  }
/***/
  template < >
  inline Real
  getNeumannBoundaryFlux(const shared_ptr<hoeb::BaseExactSolution<HOEB_MAX_ORDER> >& a_exactSolution,
                         const EBHostData<CELL, Real, 1>                     & a_phi,
                         const EBGraph                                       & a_graph,
                         const EBIndex<YFACE>                                & a_face,
                         const EBIndex<CELL>                                 & a_startingVoF,
                         const std::string                                     a_dombcname[2*DIM],
                         const std::string                                   & a_ebbcname,
                         const shared_ptr< GeometryService<HOEB_MAX_ORDER> > & a_geoserv,
                         const Chombo4::Box                                  & a_domain,
                         Real a_dx, int a_facedir,
                         Side::LoHiSide a_vofside,
                         const Chombo4::DataIndex& a_dit, bool a_divideByArea) 
  {
    const auto& facld = a_geoserv->getYFaceData(a_domain);
    const auto& dat   = (*facld)[a_dit];
    Real gradval = a_exactSolution->
      evaluateGradDMO(a_graph, a_dx, dat, a_face, 1, a_divideByArea);
    return gradval;
    
  }
/***/
  template < >
  inline Real
  getNeumannBoundaryFlux(const shared_ptr<hoeb::BaseExactSolution<HOEB_MAX_ORDER> >& a_exactSolution,
                         const EBHostData<CELL, Real, 1>                     & a_phi,
                         const EBGraph                                       & a_graph,
                         const EBIndex<ZFACE>                                & a_face,
                         const EBIndex<CELL>                                 & a_startingVoF,
                         const std::string                                     a_dombcname[2*DIM],
                         const std::string                                   & a_ebbcname,
                         const shared_ptr< GeometryService<HOEB_MAX_ORDER> > & a_geoserv,
                         const Chombo4::Box                                  & a_domain,
                         Real a_dx, int a_facedir,
                         Side::LoHiSide a_vofside,
                         const Chombo4::DataIndex& a_dit, bool a_divideByArea) 
  {
    const auto& facld = a_geoserv->getZFaceData(a_domain);
    const auto& dat   = (*facld)[a_dit];
    Real gradval = a_exactSolution->
      evaluateGradDMO(a_graph, a_dx, dat, a_face, 2, a_divideByArea);
    return gradval;
    
  }
/***/
  template < >
  inline Real
  getNeumannBoundaryFlux(const shared_ptr<hoeb::BaseExactSolution<HOEB_MAX_ORDER> >& a_exactSolution,
                         const EBHostData<CELL, Real, 1>                     & a_phi,
                         const EBGraph                                       & a_graph,
                         const EBIndex<BOUNDARY>                             & a_face,
                         const EBIndex<CELL>                                 & a_startingVoF,
                         const std::string                                     a_dombcname[2*DIM],
                         const std::string                                   & a_ebbcname,
                         const shared_ptr< GeometryService<HOEB_MAX_ORDER> > & a_geoserv,
                         const Chombo4::Box                                  & a_domain,
                         Real a_dx, int a_facedir,
                         Side::LoHiSide a_vofside,
                         const Chombo4::DataIndex& a_dit, bool a_divideByArea) 
  {
    const auto& ebfacldx = a_geoserv->getEBNormalData_x(a_domain);
    const auto& ebfacldy = a_geoserv->getEBNormalData_y(a_domain);
    const auto& ebfacldz = a_geoserv->getEBNormalData_z(a_domain);
    const auto& ebmom    = a_geoserv->getEBFaceData(a_domain);
    EBIndex<CELL> vof  = a_face.getVoF();
    Real gradval = a_exactSolution->evaluateGradDIM(a_graph, a_dx,
                                                    ebfacldx, ebfacldy, ebfacldz, ebmom,
                                                    a_dit, vof, a_divideByArea);
    return gradval;
  }
  ////
  inline void
  shiftMomentAndFillRow(LAPACKMatrix                            & a_mat,
                        IndexedMoments<    DIM, HOEB_MAX_ORDER> & a_bcmom,
                        const IndexTM<Real, DIM>                & a_xbardistance,
                        const unsigned int                      & a_currentRow)
  {
                          
    //shift so it has the same xbar
    a_bcmom.shift(a_xbardistance);
    Real volume = a_bcmom[IndexTM<int, DIM>::Zero];
    bool divide = (volume > 1.0e-16);
    if(!divide)
    {
      a_mat.setSmallCellRow(a_currentRow);
    }
    else
    {
      for(MomentIterator<DIM, HOEB_MAX_ORDER> momit; momit.ok(); ++momit)
      {
        auto momind = momit();
        unsigned int currentCol = IndexedMoments<DIM, HOEB_MAX_ORDER>::indexOf(momind);
        a_mat(a_currentRow, currentCol) = a_bcmom[momind]/volume;
      }
    }
  }
  /////
  inline void
  fillMatrixRowForEBBC(LAPACKMatrix                                        & a_mat,
                       const shared_ptr< GeometryService<HOEB_MAX_ORDER> > & a_geoserv,
                       const Chombo4::Box                                  & a_srcDomain,
                       const EBIndex<BOUNDARY>                             & a_face,
                       const EBGraph                                       & a_graph,
                       const IndexTM<Real, DIM>                            & a_xbardistance,
                       const Chombo4::DataIndex                            & a_datind,
                       const unsigned int                                  & a_currentRow,
                       int a_facedir, Real a_dx, string a_bcname)
  {
    IndexedMoments<    DIM, HOEB_MAX_ORDER> bcmom;
    if(a_bcname == string("Dirichlet"))
    {
      bcmom =  hoeb_basics::getEBIndexedMoment(a_face, a_geoserv, a_srcDomain,
                                  a_datind, a_graph, a_dx);
    }
    else if(a_bcname == string("Neumann"))
    {
      bcmom.setVal(0.);
      for(unsigned int normDir = 0; normDir < DIM; normDir++)
      {
        for(MomentIterator<DIM, HOEB_MAX_ORDER> momit; momit.ok(); ++momit)
        {
          auto momind = momit();
          //Real bcval = 0;
          
          IndexedMoments<DIM, HOEB_MAX_ORDER> normMom=
            hoeb_basics::getEBNormalMoment(a_face, a_geoserv, a_srcDomain,
                              a_datind, a_graph, normDir, a_dx);

          if(momind[normDir] > 0)
          {
            int pd = momind[normDir];
            IndexTM<int, DIM> derind =momind - IndexTM<int, DIM>::BASISV(normDir);
            Real derval = normMom[derind];
            bcmom[momind] += pd*derval;
          }
        }
      }
    }
    shiftMomentAndFillRow(a_mat, bcmom, a_xbardistance, a_currentRow);
  }

  ///
  template < CENTERING cent>
  void
  fillMatrixRowForDomBC
  (LAPACKMatrix                                                         & a_mat,
   const HostIrregData<cent, IndexedMoments<DIM-1, HOEB_MAX_ORDER>, 1 > & a_facdatfab,
   const EBIndex<cent>                                                  & a_face,
   const EBGraph                                                        & a_graph,
   const IndexTM<Real, DIM>                                             & a_xbardistance,
   const unsigned int                                                   & a_currentRow,
   int a_facedir, Real  a_dx, string a_bcname)
  {
                       
    IndexedMoments<DIM-1, HOEB_MAX_ORDER>
      facmom = hoeb_basics::getFaceIndexedMoment<DIM-1,cent>(a_face, a_facdatfab, a_graph, a_dx, a_facedir);
    //need to translate this into higher dim because that is what the matrix sees
    int facedir = 0;
    IndexedMoments<    DIM, HOEB_MAX_ORDER> volmom = hoeb_basics::getVolMomFromFacMom<DIM>(facmom, facedir);
    IndexedMoments<    DIM, HOEB_MAX_ORDER> bcmom;
    if(a_bcname == string("Dirichlet"))
    {
      bcmom = volmom;
    }
    else if(a_bcname == string("Neumann"))
    {
      for(MomentIterator<DIM, HOEB_MAX_ORDER> momit; momit.ok(); ++momit)
      {
        auto momind = momit();
        Real bcval = 0;
        if(momind[a_facedir] > 0)
        {
          int pd = momind[a_facedir];
          IndexTM<int, DIM> derind =momind - IndexTM<int, DIM>::BASISV(a_facedir);
          Real derval = volmom[derind];
          bcval = pd*derval;
        }
        bcmom[momind] = bcval;
      }
    }
    else
    {
      PROTO_ASSERT(false, "fillmatrixrow: bogus bc name");
    }
    shiftMomentAndFillRow(a_mat, bcmom, a_xbardistance, a_currentRow);
  }


  /******/  
  template <CENTERING cent>
  LAPACKMatrix
  getMomentMatrix(LAPACKMatrix                                        & a_weights,
                  const Neighborhood<cent>                            & a_block,
                  const shared_ptr< GeometryService<HOEB_MAX_ORDER> > & a_geoserv,
                  const Chombo4::Box                                  & a_srcDomain,
                  const RealVect                                      & a_xbar,
                  const Chombo4::DataIndex                            & a_dit,
                  Real                                                  a_dx)
  {
    typedef IndexedMoments<DIM  , HOEB_MAX_ORDER> IndMomDIM;
    typedef IndexedMoments<DIM-1, HOEB_MAX_ORDER> IndMomSDMinOne;
    unsigned int Pcols = IndMomDIM::size();
    unsigned int Nrows = a_block.size();
    if(Nrows < Pcols)
    {
      Chombo4::MayDay::Warning("not enough degrees of freedom");
    }
    LAPACKMatrix retval(Nrows, Pcols);
    a_weights.define(Nrows, Nrows);
    a_weights.setVal(0.);
    auto dit = a_geoserv->getDBL(a_srcDomain).dataIterator();
    
    const auto & graphsldptr =  a_geoserv->getGraphs(  a_srcDomain);
    const auto & volmomldptr =  a_geoserv->getVoluData(a_srcDomain);
    
    const auto & graph     = (*graphsldptr)[a_dit];
    const auto & volmomfab = (*volmomldptr)[a_dit];

    unsigned int currentRow = 0;
    {
      //volumes
      for(unsigned int ivol = 0; ivol < a_block.m_volumes.size(); ivol++)
      {
        EBIndex<cent> vof    = a_block.m_volumes[ivol];
        RealVect vofloc = hoeb_basics::getVoFLocation(vof, a_dx);
        IndexTM<Real, DIM> distance = hoeb_basics::getIndexDistance(vofloc, a_xbar);
        IndexedMoments<DIM  , HOEB_MAX_ORDER>
          volmom = hoeb_basics::getIndexedMoment<DIM, cent>(vof, volmomfab, graph, a_dx);

        shiftMomentAndFillRow(retval, volmom, distance, currentRow);
        
        a_weights(currentRow, currentRow) = a_block.m_volumesWgt[currentRow];
        currentRow++;
      }
    }
    {
      //now for the  domain x faces
      const auto & xfacldptr  =  a_geoserv->getXFaceData(a_srcDomain);
      const auto & xfacdatfab = (*xfacldptr)[a_dit];

      for(unsigned int ivol = 0; ivol < a_block.m_xfacesDiri.size(); ivol++)
      {
        const auto & face = a_block.m_xfacesDiri[ivol];
        RealVect vofloc = hoeb_basics::getFaceLocation<XFACE>(face, a_dx, 0);
        IndexTM<Real, DIM> distance = hoeb_basics::getIndexDistance(vofloc, a_xbar);
        fillMatrixRowForDomBC<XFACE>(retval, xfacdatfab, face, graph, distance,
                                     currentRow, 0, a_dx, string("Dirichlet"));

        a_weights(currentRow,currentRow) = a_block.m_xfacesDiriWgt[ivol];
        currentRow++;
      }
      for(unsigned int ivol = 0; ivol < a_block.m_xfacesNeum.size(); ivol++)
      {
        const auto & face = a_block.m_xfacesNeum[ivol];
        RealVect vofloc = hoeb_basics::getFaceLocation<XFACE>(face, a_dx, 0);
        IndexTM<Real, DIM> distance = hoeb_basics::getIndexDistance(vofloc, a_xbar);
        fillMatrixRowForDomBC<XFACE>(retval, xfacdatfab, face, graph, distance,
                                     currentRow, 0, a_dx, string("Neumann"));

        a_weights(currentRow,currentRow) = a_block.m_xfacesNeumWgt[ivol];
        currentRow++;
      }
    }
    {
      //now for the  domain y faces
      const auto & yfacldptr  =  a_geoserv->getYFaceData(a_srcDomain);
      const auto & yfacdatfab = (*yfacldptr)[a_dit];
      for(unsigned int ivol = 0; ivol < a_block.m_yfacesDiri.size(); ivol++)
      {
        const auto & face = a_block.m_yfacesDiri[ivol];
        RealVect vofloc = hoeb_basics::getFaceLocation<YFACE>(face, a_dx, 1);
        IndexTM<Real, DIM> distance = hoeb_basics::getIndexDistance(vofloc, a_xbar);
        fillMatrixRowForDomBC<YFACE>(retval, yfacdatfab, face, graph, distance,
                                     currentRow, 1, a_dx, string("Dirichlet"));
      
        a_weights(currentRow, currentRow) = a_block.m_yfacesDiriWgt[ivol];
        currentRow++;
      }
      for(unsigned int ivol = 0; ivol < a_block.m_yfacesNeum.size(); ivol++)
      {
        const auto & face = a_block.m_yfacesNeum[ivol];
        RealVect vofloc             = hoeb_basics::getFaceLocation<YFACE>(face, a_dx, 1);
        IndexTM<Real, DIM> distance = hoeb_basics::getIndexDistance(vofloc, a_xbar);
        fillMatrixRowForDomBC<YFACE>(retval, yfacdatfab, face, graph, distance,
                                     currentRow, 1, a_dx, string("Neumann"));
      
        a_weights(currentRow, currentRow) = a_block.m_yfacesNeumWgt[ivol];
        currentRow++;
      }
    }
#if DIM==3
    {
      //now for the  domain z faces
      const auto & zfacldptr  =  a_geoserv->getZFaceData(a_srcDomain);
      const auto & zfacdatfab = (*zfacldptr)[a_dit];
      for(unsigned int ivol = 0; ivol < a_block.m_zfacesDiri.size(); ivol++)
      {
        const auto & face = a_block.m_zfacesDiri[ivol];
        RealVect vofloc             = hoeb_basics::getFaceLocation<ZFACE>(face, a_dx, 2);
        IndexTM<Real, DIM> distance = hoeb_basics::getIndexDistance(vofloc, a_xbar);
        fillMatrixRowForDomBC<ZFACE>(retval, zfacdatfab, face, graph, distance,
                                     currentRow, 2, a_dx, string("Dirichlet"));
      
        a_weights(currentRow,currentRow) = a_block.m_zfacesDiriWgt[ivol];
        currentRow++;
      }
      for(unsigned int ivol = 0; ivol < a_block.m_zfacesNeum.size(); ivol++)
      {
        const auto & face = a_block.m_zfacesNeum[ivol];
        RealVect vofloc             = hoeb_basics::getFaceLocation<ZFACE>(face, a_dx, 2);
        IndexTM<Real, DIM> distance = hoeb_basics::getIndexDistance(vofloc, a_xbar);
        fillMatrixRowForDomBC<ZFACE>(retval, zfacdatfab, face, graph, distance,
                                     currentRow, 2, a_dx, string("Neumann"));
      
        a_weights(currentRow,currentRow) = a_block.m_zfacesNeumWgt[ivol];
        currentRow++;
      }
    }
#endif
    {
      //Cut faces
      for(unsigned int ivol = 0; ivol < a_block.m_cutfacesDiri.size(); ivol++)
      {
        const auto & face = a_block.m_cutfacesDiri[ivol];
        //-1 gets the cell centered location
        RealVect vofloc             = hoeb_basics::getFaceLocation<BOUNDARY>(face, a_dx, -1);
        IndexTM<Real, DIM> distance = hoeb_basics::getIndexDistance(vofloc, a_xbar);
        fillMatrixRowForEBBC(retval, a_geoserv, a_srcDomain, face, graph, distance,
                             a_dit, currentRow, 2, a_dx, 
                             string("Dirichlet"));
      
        a_weights(currentRow,currentRow) = a_block.m_cutfacesDiriWgt[ivol];
        currentRow++;
      }
      for(unsigned int ivol = 0; ivol < a_block.m_cutfacesNeum.size(); ivol++)
      {
        const auto & face = a_block.m_cutfacesNeum[ivol];
        //-1 gets the cell centered location
        RealVect vofloc             = hoeb_basics::getFaceLocation<BOUNDARY>(face, a_dx, -1);
        IndexTM<Real, DIM> distance = hoeb_basics::getIndexDistance(vofloc, a_xbar);
        fillMatrixRowForEBBC(retval, a_geoserv, a_srcDomain, face, graph, distance,
                             a_dit, currentRow, 2, a_dx, 
                             string("Neumann"));
      
        a_weights(currentRow,currentRow) = a_block.m_cutfacesNeumWgt[ivol];
        currentRow++;
      }
    }

    return retval;
  }  


  //Moore-Penrose inverse of the moment matrix is calulcated using weighted least squares
  inline LAPACKMatrix
  getMoorePenroseInverse(const LAPACKMatrix & a_M, bool a_printStuff)
  {


    LAPACKMatrix AT = a_M;
    AT.transpose();
    LAPACKMatrix ATA;
    multiply(ATA, AT, a_M);
    LAPACKMatrix ATAinv = ATA;
    ATAinv.invert();
    
    if(a_printStuff)
    {
      using Chombo4::pout;
      LAPACKMatrix Amat = a_M;
      LAPACKMatrix checkMat;
      multiply(checkMat, ATAinv, ATA);
      pout() << "A*Ainv = " << endl;
      checkMat.poutAll();
      checkMat.poutMaxMins();

      static int filenumber = 0;
      string filemat = string("Amat.") + std::to_string(filenumber) + string(".matrix");
      string fileata = string("ATA.") + std::to_string(filenumber) + string(".matrix");
      string fileinv = string("ATAinv.") + std::to_string(filenumber) + string(".matrix");
      string filechk = string("ATAATAinv.") + std::to_string(filenumber) + string(".matrix");
      checkMat.writeToFile(filechk);
      ATA.writeToFile(     fileata);
      Amat.writeToFile(    filemat);
      ATAinv.writeToFile(  fileinv);
      filenumber++;
    }

    return ATAinv;
  }
  ///
  inline      LAPACKMatrix
  getStencilMatrix(const LAPACKMatrix & a_Qmat,
                   const LAPACKMatrix & a_Gmat,
                   const LAPACKMatrix & a_AvInvmat,
                   const LAPACKMatrix & a_weight)
  {
    LAPACKMatrix GTmat = a_Gmat;
    GTmat.transpose();
    //this is GT W
    LAPACKMatrix GTW;
    multiply(GTW, GTmat, a_weight);
    
    //this is Av-1 GT W
    LAPACKMatrix AvinvGTW;
    multiply(AvinvGTW, a_AvInvmat, GTW);
    
    LAPACKMatrix Smat;
    multiply(Smat, a_Qmat, AvinvGTW);
    return Smat;
  }

   
  ///
  template <CENTERING cent>
  inline LAPACKMatrix
  getDharshiQMatrix(const EBIndex<cent>                                 & a_face,
                    const Neighborhood<CELL>                            & a_block,
                    const shared_ptr< GeometryService<HOEB_MAX_ORDER> > & a_geoserv,
                    const Chombo4::Box                                  & a_srcDomain,
                    const RealVect                                      & a_xbar,
                    Real                                                  a_dx,
                    int                                                   a_facedir,
                    const Chombo4::DataIndex                            & a_datind)
  {
    IndexedMoments<DIM-1, HOEB_MAX_ORDER> 
      facmom = hoeb_basics::getIndexedMomentFace<cent>(a_face, a_geoserv, a_srcDomain, a_datind, a_dx);
    IndexedMoments<    DIM, HOEB_MAX_ORDER> volmom = hoeb_basics::getVolMomFromFacMom<DIM>(facmom, a_facedir);

    typedef       IndexTM<int, DIM-1>          IvDMO;
    IvDMO zeroind = IvDMO::Zero;

    Real momzero = facmom[zeroind];
    RealVect facloc             = hoeb_basics::getFaceLocation<cent>(a_face, a_dx, a_facedir);
    IndexTM<Real, DIM> distance = hoeb_basics::getIndexDistance(facloc, a_xbar);
    volmom.shift(distance);

    typedef IndexedMoments<DIM  , HOEB_MAX_ORDER> IndMomDIM;
    unsigned int Pcols = IndMomDIM::size();
    LAPACKMatrix Qmat(1, Pcols);


    for(MomentIterator<DIM, HOEB_MAX_ORDER> momit; momit.ok(); ++momit)
    {
      auto momind = momit();
      Real matval = 0;
      if(momind[a_facedir] > 0)
      {
        int pd = momind[a_facedir];
        IndexTM<int, DIM> derind =momind - IndexTM<int, DIM>::BASISV(a_facedir);
        Real derval = volmom[derind];
        if(std::abs(momzero) > 1.0e-10)
        {
          matval = pd*derval/momzero;
        }
      }
      unsigned int currentCol = IndexedMoments<DIM, HOEB_MAX_ORDER>::indexOf(momind);
      Qmat(0,currentCol) = matval;
    }
    
    return Qmat;
  }

  ///template specialization to get the eb flux q matrix
  template < >
  inline LAPACKMatrix
  getDharshiQMatrix(const EBIndex<BOUNDARY>                             & a_face,
                    const Neighborhood<CELL>                            & a_block,
                    const shared_ptr< GeometryService<HOEB_MAX_ORDER> > & a_geoserv,
                    const Chombo4::Box                                  & a_srcDomain,
                    const RealVect                                      & a_xbar,
                    Real                                                  a_dx,
                    int                                                   a_facedir,
                    const Chombo4::DataIndex                            & a_datind)
  {
    typedef IndexedMoments<DIM  , HOEB_MAX_ORDER> IndMomDIM;
    unsigned int Pcols = IndMomDIM::size();
    LAPACKMatrix Qmat(1, Pcols);

    Qmat.setVal(0.);
    const auto & graphsldptr = a_geoserv->getGraphs(    a_srcDomain);
    const auto & graph = (*graphsldptr)[a_datind];
    for(unsigned int normDir = 0; normDir < DIM; normDir++)
    {
      IndexedMoments<DIM, HOEB_MAX_ORDER> normMom=
        hoeb_basics::getEBNormalMoment(a_face, a_geoserv, a_srcDomain,
                                       a_datind, graph, normDir, a_dx);

      //begin debug
      //normMom.print();
      //end debug
      for(MomentIterator<DIM, HOEB_MAX_ORDER> momit; momit.ok(); ++momit)
      {
        auto momind = momit();
        if(momind[normDir] > 0)
        {
          int pd = momind[normDir];
          IndexTM<int, DIM> derind =momind - IndexTM<int, DIM>::BASISV(normDir);
          Real derval = normMom[derind];
          unsigned int currentCol = IndexedMoments<DIM, HOEB_MAX_ORDER>::indexOf(momind);
          Qmat(0,currentCol) += pd*derval;
        }
      }
    }
    return Qmat;
  }
  
  ///  This function very much relies upon the ordering  imposed in getMomentMatrix
  inline 
  CompositeStencil
  getStencilFromMatrix(const LAPACKMatrix      & a_stenmat,
                       const Neighborhood<CELL>& a_block)
  {
    CompositeStencil retval;
    //(set in getMomentMatrix)  first bits are cells
    int imat = 0;
    for(unsigned int ivol = 0; ivol < a_block.m_volumes.size(); ivol++)
    {
      const auto& vof = a_block.m_volumes[ivol];
      Real wgt = a_stenmat(0, imat);
      retval.m_cellSten.add(vof, wgt);
      imat ++;
    }
    //(set in getMomentMatrix)  second is x domain faces diri then neum
    for(unsigned int ifac = 0; ifac < a_block.m_xfacesDiri.size(); ifac++)
    {
      const auto& domFace = a_block.m_xfacesDiri[ifac];
      Real wgt = a_stenmat(0, imat);
      retval.m_xfaceDiriSten.add(domFace, wgt);
      imat ++;
    }
    for(unsigned int ifac = 0; ifac < a_block.m_xfacesNeum.size(); ifac++)
    {
      const auto& domFace = a_block.m_xfacesNeum[ifac];
      Real wgt = a_stenmat(0, imat);
      retval.m_xfaceNeumSten.add(domFace, wgt);
      imat ++;
    }

    //(set in getMomentMatrix)  third is y domain faces diri then neum
    for(unsigned int ifac = 0; ifac < a_block.m_yfacesDiri.size(); ifac++)
    {
      const auto& domFace = a_block.m_yfacesDiri[ifac];
      Real wgt = a_stenmat(0, imat);
      retval.m_yfaceDiriSten.add(domFace, wgt);
      imat ++;
    }
    for(unsigned int ifac = 0; ifac < a_block.m_yfacesNeum.size(); ifac++)
    {
      const auto& domFace = a_block.m_yfacesNeum[ifac];
      Real wgt = a_stenmat(0, imat);
      retval.m_yfaceNeumSten.add(domFace, wgt);
      imat ++;
    }

#if DIM==3    
    //(set  in getMomentMatrix)  fourth is z domain faces diri then neum
    for(unsigned int ifac = 0; ifac < a_block.m_zfacesDiri.size(); ifac++)
    {
      const auto& domFace = a_block.m_zfacesDiri[ifac];
      Real wgt = a_stenmat(0, imat);
      retval.m_zfaceDiriSten.add(domFace, wgt);
      imat ++;
    }
    for(unsigned int ifac = 0; ifac < a_block.m_zfacesNeum.size(); ifac++)
    {
      const auto& domFace = a_block.m_zfacesNeum[ifac];
      Real wgt = a_stenmat(0, imat);
      retval.m_zfaceNeumSten.add(domFace, wgt);
      imat ++;
    }
#endif
    
    //(set  in getMomentMatrix)  last are the cut faces --- diri then neum
    for(unsigned int ifac = 0; ifac < a_block.m_cutfacesDiri.size(); ifac++)
    {
      const auto& domFace = a_block.m_cutfacesDiri[ifac];
      Real wgt = a_stenmat(0, imat);
      retval.m_cutfaceDiriSten.add(domFace, wgt);
      imat ++;
    }
    for(unsigned int ifac = 0; ifac < a_block.m_cutfacesNeum.size(); ifac++)
    {
      const auto& domFace = a_block.m_cutfacesNeum[ifac];
      Real wgt = a_stenmat(0, imat);
      retval.m_cutfaceNeumSten.add(domFace, wgt);
      imat ++;
    }
    return retval;
  }
                       
  /****/
  template <CENTERING cent>
  void
  restrictFlux(EBLevelBoxData<cent, 1>                                           &  a_phiFToC,
               const EBLevelBoxData<cent, 1>                                     &  a_phiFine,
               const shared_ptr<LevelData<EBGraph> >                             &  a_graphsFine,
               const Chombo4::DisjointBoxLayout                                  &  a_gridsFine,
               const Chombo4::Box                                                &  a_domFine,
               const Real                                                        &  a_dxFine,
               const shared_ptr<LevelData<EBGraph> >                             &  a_graphsCoar,
               const Chombo4::DisjointBoxLayout                                  &  a_gridsCoar,
               const Chombo4::Box                                                &  a_domCoar,
               const Real                                                        &  a_dxCoar,
               const shared_ptr<EBDictionary<HOEB_MAX_ORDER, Real, cent, cent> > &  a_dictionary,
               const shared_ptr< GeometryService<HOEB_MAX_ORDER> >               &  a_geoserv)
  {
    std::string  nobcname        = string("no_bcs");
    std::string  restrictionName = string("face_restriction");
    string dombc[2*DIM];
    for(unsigned int idom = 0; idom < 2*DIM;  idom++)
    {
      dombc[idom] = nobcname;
    }
    a_dictionary->registerStencil(restrictionName, dombc, nobcname, a_domFine, a_domCoar, false);
    Chombo4::DataIterator dit = a_gridsCoar.dataIterator();
    for(unsigned int ibox = 0; ibox < dit.size(); ++ibox)
    {
      auto      & coarfab = a_phiFToC[dit[ibox]];
      const auto& finefab = a_phiFine[dit[ibox]];
      auto stencil = a_dictionary->getEBStencil(restrictionName, nobcname, a_domFine, a_domCoar, ibox);
      //set resc = Ave(resf) (true is initToZero)
      stencil->apply(coarfab, finefab,  true, 1.0);
    }
  }


  
  ///
  
  template <CENTERING cent>
  inline CompositeStencil
  getDevendranIntFluxDACompositeStencil(const EBIndex<cent>                                 & a_face,
                                        const EBIndex<CELL>                                 & a_startingVoF,
                                        const std::string                                     a_dombcname[2*DIM],
                                        const std::string                                   & a_ebbcname,
                                        const shared_ptr< GeometryService<HOEB_MAX_ORDER> > & a_geoserv,
                                        const Chombo4::Box                                  & a_srcDomain,
                                        const Chombo4::DataIndex                            & a_dit,
                                        Real a_dx, int a_facedir,
                                        Side::LoHiSide a_vofside,
                                        int a_stenRadius, int a_maxStenRadius,
                                        bool a_printStuff)
  {
    return
      getDevendranIntFluxDACompositeStencil(a_face,
                                            a_dombcname[2*DIM],
                                            a_ebbcname,
                                            a_geoserv,
                                            a_srcDomain,
                                            a_dit,
                                            a_dx, a_facedir,a_vofside,
                                            a_stenRadius, a_maxStenRadius,
                                            a_printStuff);
  }

  ///
  /**
     Return integral F dA  (do NOT divide out the area).
     For EB Fluxes, it is Integral F-n dA with an outward-facing normal.
     starting vof never used
  */
  template <CENTERING cent>
  inline CompositeStencil
  getDevendranIntFluxDACompositeStencil(const EBIndex<cent>                                 & a_face,
                                        const std::string                                     a_dombcname[2*DIM],
                                        const std::string                                   & a_ebbcname,
                                        const shared_ptr< GeometryService<HOEB_MAX_ORDER> > & a_geoserv,
                                        const Chombo4::Box                                  & a_srcDomain,
                                        const Chombo4::DataIndex& a_dit,
                                        Real a_dx, int a_facedir,
                                        Side::LoHiSide a_vofside,
                                        int a_stenRadius, int a_maxStenRadius,
                                        bool a_printStuff)
  {
    //facedir == -1 means the EB
    const auto & graphsldptr = a_geoserv->getGraphs(    a_srcDomain);

    const auto & datind = a_dit;
    const auto & graph = (*graphsldptr)[datind];
    CompositeStencil  retval;

    //when I know the answer, short circuit all the other stuff and return it
    int iindex = -1;
    if(a_facedir >= 0)
    {
      iindex = ebp_index(a_facedir, a_vofside);
    }
    if((a_facedir==-1) && (a_ebbcname  == string("Neumann")))
    {
      return retval;
    }
    else if ((a_facedir >= 0) && a_face.m_isBoundary && (a_dombcname[iindex] == string("Neumann")))
    {
      return retval;
    }
    else
    {
      bool useBoundaries = true;

      //get the neighborhood near the target face
      Neighborhood<CELL>  block =
        getCellNeighbors<cent>(a_face, graph, a_dombcname, a_ebbcname, a_facedir,
                               a_stenRadius, a_maxStenRadius, useBoundaries,
                               a_printStuff);
      if(a_printStuff)
      {
        Chombo4::pout() << "in getDevendranIntFluxDACompositeStencil " << endl;
        Chombo4::pout() << "face = "<< a_face << endl;
        block.poutAll();
      }
      const EBIndex<cent>& face = a_face;
      Real                 dx = a_dx;
      int                  facedir = a_facedir;
      RealVect xbar =   hoeb_basics::getFaceLocation(face, dx, facedir);

      if(a_printStuff)
      {
        Chombo4::pout() << "xbar = " << xbar <<  endl;
        
      }

      
      //this is the diagonal weighting matrix.
      LAPACKMatrix weight;
      //this matrix contains all the moments and boundary condition constraints
      LAPACKMatrix Mmat = getMomentMatrix<CELL>(weight, block, a_geoserv,
                                                a_srcDomain, xbar, a_dit, a_dx);

      if(a_printStuff)
      {
        Chombo4::pout() << "moment matrix = " << endl;
        Mmat.poutAll();
        Chombo4::pout() << "weight matrix = " << endl;
        weight.poutDiag();
      }
      
      LAPACKMatrix Gmat; //G = WM
      multiply(Gmat, weight, Mmat);
      
      if(a_printStuff)
      {
        Chombo4::pout() << "Gmat = " << endl;
        Gmat.poutAll();
      }

      if(a_printStuff)
      {
        using Chombo4::pout;
        static int filenumber = 0;
        string filegmat = string("Gmat.") + std::to_string(filenumber) + string(".matrix");
        string filemmat = string("Mmat.") + std::to_string(filenumber) + string(".matrix");
        string filewgt  = string("weight") + std::to_string(filenumber) + string(".matrix");
        Gmat.writeToFile(filegmat);
        Mmat.writeToFile(filemmat);
        weight.writeToFile(filewgt);

        filenumber++;
      }
      //Moore-Penrose inverse of the moment matrix
      //is calulcated using weighted least squares
      //because G is weighted version of M
      LAPACKMatrix GTGinvMat = getMoorePenroseInverse(Gmat, a_printStuff);
      
      if(a_printStuff)
      {
        using Chombo4::pout;
        Chombo4::pout() << "GTGinvMat = " << endl;
        GTGinvMat.poutAll();
      }

      //this is the matrix that defines the operator
      LAPACKMatrix Qmat =
        getDharshiQMatrix<cent>
        (a_face, block, a_geoserv, a_srcDomain, xbar,  a_dx, a_facedir, datind);

      if(a_printStuff)
      {
        Chombo4::pout() << "Qmat = " << endl;
        Qmat.poutAll();
      }
      Real qmatzero  =  Qmat(0, 0);
      int ideb = 0;
      if(std::abs(qmatzero) > 1.0e-6)
      {
        ideb = 1;
      }
      LAPACKMatrix Smat = getStencilMatrix(Qmat, Gmat, GTGinvMat, weight);
      
      retval = hoeb::getStencilFromMatrix(Smat, block);
    //begin debug
   // 
   //   retval.m_cellSten.clear();
   //   Real coeff = 1.; //flux *da, in 2D takes out dx 
   //   EBIndex<CELL> vofhi = a_face.getVoF(Side::Hi); 
   //   EBIndex<CELL> voflo = a_face.getVoF(Side::Lo); 
   //   retval.m_cellSten.add(vofhi, coeff);
   //   retval.m_cellSten.add(voflo,-coeff);
    //end debug
      
      if(a_printStuff)
      {
        Real sumStenc = retval.m_cellSten.sum();
        Chombo4::pout() << "composite stencil: "<<  endl;
        retval.poutAll();
        Chombo4::pout() << "composite stencil sum of homogeneous bits =" << sumStenc <<  endl;
      }
    }
    return retval;
  }
  ///
  /**
     Return integral F dA  (do NOT divide out the area).
     For EB Fluxes, it is Integral F-n dA with an outward-facing normal.
     this only returns the bit that acts on cell data.   The other stuff is assumed to be zero
     (homogeneous operator--for use with EBStencil).
  */
  
  template <CENTERING cent>
  inline LocalStencil<CELL, Real> 
  getDharshiIntFluxDAStencil(const EBIndex<cent>                                 & a_face,
                             const EBIndex<CELL>                                 & a_startingVoF,
                             const std::string                                     a_dombcname[2*DIM],
                             const std::string                                   & a_ebbcname,
                             const shared_ptr< GeometryService<HOEB_MAX_ORDER> > & a_geoserv,
                             const Chombo4::Box                                  & a_srcDomain,
                             const Chombo4::DataIndex                            & a_dit,
                             Real a_dx, int a_facedir, Side::LoHiSide a_vofside,
                             int a_stenRadius, int a_maxStenRadius,
                             bool a_printStuff)
  {
    CompositeStencil composite
      = getDevendranIntFluxDACompositeStencil(a_face,     
                                              a_dombcname,
                                              a_ebbcname, 
                                              a_geoserv,  
                                              a_srcDomain,
                                              a_dit, a_dx, a_facedir, a_vofside,
                                              a_stenRadius, a_maxStenRadius,
                                              a_printStuff);
    return composite.m_cellSten;
  }



  /******/  
  inline LocalStencil<CELL, Real> 
  getHomogeneousDharshiStencil(const EBIndex<CELL>                                 & a_vof,
                               const std::string                                     a_dombcname[2*DIM],
                               const std::string                                   & a_ebbcname,
                               const shared_ptr< GeometryService<HOEB_MAX_ORDER> > & a_geoserv,
                               Proto::Box                                          & a_srcDomain,
                               unsigned int a_ibox, Real a_alpha, Real a_beta, Real a_dx,
                               int a_stenRadius, int a_maxStenRadius,
                               bool a_printStuff)
  
  {
    LocalStencil<CELL, Real> vofsten;
    //use the age-old trick of building up the stencil by construction
    auto dbl = a_geoserv->getDBL(a_srcDomain);
    auto dit = dbl.dataIterator();
    const auto & graphsldptr = a_geoserv->getGraphs(    a_srcDomain);
    const auto & graph = (*graphsldptr)[dit[a_ibox]];

    if(a_printStuff)
    {
      Chombo4::pout() << "GetHomogeneousDharshiStencil for vof = :: " << a_vof.m_pt << endl;
    };
    
    for(SideIterator sit; sit.ok(); ++sit)
    {
      int isign = sign(sit());
      {
        auto xfaces = graph.getXFaces(a_vof, sit());
        for(unsigned int iface = 0; iface < xfaces.size(); iface++)
        {
          // only really want to print one face's worth of stuff
          bool printstuff = a_printStuff && (isign == 1);
          auto face = xfaces[iface];
          if(printstuff)
          {
            Chombo4::pout() << "GetHomogeneousDharshiStencil:: ";
            Chombo4::pout() << "X face " << iface << " = " << face.m_pt << endl;
          }
          LocalStencil<CELL, Real> fluxsten =
            getDharshiIntFluxDAStencil<XFACE>(face,
                                              a_vof,a_dombcname,a_ebbcname,
                                              a_geoserv, a_srcDomain, dit[a_ibox],
                                              a_dx, 0, sit(), 
                                              a_stenRadius, a_maxStenRadius, printstuff);
          fluxsten *= Real(isign);
          vofsten += fluxsten;
        }
      }
      {                                                // 
        auto yfaces = graph.getYFaces(a_vof, sit());
        for(unsigned int iface = 0; iface < yfaces.size(); iface++)
        {
          // only really want to print one face's worth of stuff
          bool printstuff = a_printStuff && (isign == 1);
          //comment if you want to print y stuff
          printstuff = false;
          
          auto face = yfaces[iface];
          if(printstuff)
          {
            Chombo4::pout() << "Y face " << iface << " = " << face.m_pt << endl;
          }
          LocalStencil<CELL, Real> fluxsten =
            getDharshiIntFluxDAStencil<YFACE>(face,
                                              a_vof,a_dombcname,a_ebbcname,
                                              a_geoserv, a_srcDomain, dit[a_ibox],
                                              a_dx, 1, sit(), 
                                              a_stenRadius, a_maxStenRadius, printstuff);
          fluxsten *= Real(isign);
          vofsten += fluxsten;
        }
      }
#if DIM==3      
      {
        auto zfaces = graph.getZFaces(a_vof, sit());
        for(unsigned int iface = 0; iface < zfaces.size(); iface++)
        {
          auto face = zfaces[iface];

          LocalStencil<CELL, Real> fluxsten =
            getDharshiIntFluxDAStencil<ZFACE>(face,
                                              a_vof, a_dombcname, a_ebbcname,
                                              a_geoserv, a_srcDomain, dit[a_ibox],
                                              a_dx, 2, sit(), 
                                              a_stenRadius, a_maxStenRadius, false);
                                              
          fluxsten *= Real(isign);
          vofsten += fluxsten;
        }
      }
#endif
    }
    {
      EBIndex<BOUNDARY> face = a_vof.getCutFace();
      LocalStencil<CELL, Real>
        fluxsten =
        getDharshiIntFluxDAStencil<BOUNDARY>(face,
                                             a_vof,a_dombcname,a_ebbcname,
                                             a_geoserv, a_srcDomain, dit[a_ibox],
                                             a_dx, -1, Side::Invalid, 
                                             a_stenRadius, a_maxStenRadius, false);
                                             
      vofsten += fluxsten;
    }
    //need to divide by dx^d to get kappa*lapl(phi)
    Real cellvolinv = 1;
    for(int idir = 0; idir < DIM; idir++)
    {
      cellvolinv /= a_dx;
    }
    vofsten *= cellvolinv;
    
    return vofsten;
  }

  /******/
  inline void
  getHomogeneousDharshiStencil(string                                              & a_stencilName,
                               vector<EBIndex<CELL> >                              & a_dstVoFs,
                               vector<LocalStencil<CELL, Real> >                   & a_stencil,
                               Proto::Box                                          & a_srcValid,
                               Proto::Box                                          & a_dstValid,
                               Proto::Box                                          & a_srcDomain,
                               Proto::Box                                          & a_dstDomain,
                               bool                                                & a_needDiagonalWeights,
                               const shared_ptr< GeometryService<HOEB_MAX_ORDER> > & a_geoserv,
                               const Chombo4::DisjointBoxLayout                    & a_grids,
                               const Chombo4::Box                                  & a_domain,
                               const Real                                          & a_dx,
                               unsigned int                                          a_ibox,
                               const string                                          a_dombcarray[2*DIM],
                               const string                                        & a_ebbcName,
                               Real a_alpha, Real a_beta,
                               int a_stenRadius, int a_maxStenRadius,
                               bool a_printStuff)

  {
    /* geometric multigrid is not a great idea here*/
    a_needDiagonalWeights = false;
    a_srcDomain = ProtoCh::getProtoBox(a_domain);
    a_dstDomain = ProtoCh::getProtoBox(a_domain);
    Chombo4::ParmParse pp;
    using Chombo4::pout;
    Real alpha = a_alpha;
    Real beta =  a_beta;

    a_stencilName = string("Dharshi_Laplacian");

    const auto & graphsldptr = a_geoserv->getGraphs(    a_srcDomain);
    const auto & voldatldptr = a_geoserv->getVoluData(  a_srcDomain);
    const auto & ebfdatldptr = a_geoserv->getEBFaceData(a_srcDomain);
    const auto & xfadatldptr = a_geoserv->getXFaceData( a_srcDomain);
    const auto & yfadatldptr = a_geoserv->getYFaceData( a_srcDomain);
    const auto & zfadatldptr = a_geoserv->getZFaceData( a_srcDomain);
    const auto & dbl         = a_geoserv->getDBL(       a_srcDomain);
    auto dit = dbl.dataIterator();

    const auto & graph = (*graphsldptr)[dit[a_ibox]];
    a_srcValid  = graph.validBox();
    a_dstValid  = graph.validBox();

    int idst = 0;
    for(auto bit = a_dstValid.begin(); bit != a_dstValid.end(); ++bit)
    {
      auto vofs = graph.getVoFs(*bit);
      for(unsigned int ivof = 0; ivof < vofs.size(); ivof++)
      {
        bool printStuff = a_printStuff;
        //begin debug
        auto vofpt = vofs[ivof].m_pt;
        Point ptd1(D_DECL(1,  1 , 1 ));
        Point ptd2(D_DECL(15, 15, 15));
        int ideb = 0;
        if((vofpt == ptd1) || (vofpt == ptd2))
        {
          Chombo4::pout() << "*********************************************"  << endl;
          Chombo4::pout() << "turning  up verbosity for point " << vofpt << endl;
          Chombo4::pout() << "*********************************************"  << endl;
          ideb = 1;
          printStuff = true;
        }
        //end debug
        LocalStencil<CELL, Real> vofsten =
          getHomogeneousDharshiStencil(vofs[ivof],
                                       a_dombcarray, a_ebbcName,
                                       a_geoserv, a_srcDomain, a_ibox,
                                       alpha, beta, a_dx,
                                       a_stenRadius, a_maxStenRadius, printStuff);
        a_dstVoFs.push_back(vofs[ivof]);
        idst++;
        a_stencil.push_back(vofsten);
      }
    }
  }


  ///get the stencil for the laplacian described in Schwartz, et al.
  void
  schwartzLaplStencil(string                                              & a_stencilName,
                      string                                              & a_ebbcName,
                      vector<EBIndex<CELL> >                              & a_dstVoFs,
                      vector<LocalStencil<CELL, Real> >                   & a_stencil,
                      Proto::Box                                          & a_srcValid,
                      Proto::Box                                          & a_dstValid,
                      Proto::Box                                          & a_srcDomain,
                      Proto::Box                                          & a_dstDomain,
                      Proto::Point                                        & a_srcGhost,
                      Proto::Point                                        & a_dstGhost,
                      bool                                                & a_needDiagonalWeights,
                      const shared_ptr< GeometryService<HOEB_MAX_ORDER> > & a_geoserv,
                      const Chombo4::DisjointBoxLayout                    & a_grids,
                      const Chombo4::Box                                  & a_domain,
                      const Real                                          & a_dx,
                      unsigned int                                          a_ibox);
  
  shared_ptr<hoeb::BaseExactSolution<HOEB_MAX_ORDER> >   getBaseExactSoltuion();
  
  ///return an implicit function based on stuff in the input file
  shared_ptr<BaseIF> getImplicitFunction();

  ///
  /**
     Fill in scalar from exactsolutions class
  */
  void
  fillPhi(EBLevelBoxData<CELL, 1>                                &  a_phi,
          const shared_ptr<graph_distrib_t    >                  &  a_graphs,
          const Chombo4::DisjointBoxLayout                       &  a_grids,
          const Chombo4::Box                                     &  a_domFine,
          const Real                                             &  a_dx,
          const shared_ptr< GeometryService<HOEB_MAX_ORDER> >    &  a_geoserv);

  ///
  /**
   **/
  void
  restrictKappaLphi(EBLevelBoxData<CELL, 1>                                           &  a_klpFToC,
                    const EBLevelBoxData<CELL, 1>                                     &  a_klpFine,
                    const shared_ptr<graph_distrib_t>                                 &  a_graphsFine,
                    const Chombo4::DisjointBoxLayout                                  &  a_gridsFine,
                    const Chombo4::Box                                                &  a_domFine,
                    const Real                                                        &  a_dxFine,
                    const shared_ptr<graph_distrib_t>                                 &  a_graphsCoar,
                    const Chombo4::DisjointBoxLayout                                  &  a_gridsCoar,
                    const Chombo4::Box                                                &  a_domCoar,
                    const Real                                                        &  a_dxCoar,
                    const shared_ptr<EBDictionary<HOEB_MAX_ORDER, Real, CELL, CELL> > &  a_dictionary,
                    const shared_ptr< GeometryService<HOEB_MAX_ORDER> >               &  a_geoserv);

  ///
  /**
   **/
  void
  restrictPhi(EBLevelBoxData<CELL, 1>                                           &  a_phiFToC,
              const EBLevelBoxData<CELL, 1>                                     &  a_phiFine,
              const shared_ptr<graph_distrib_t    >                             &  a_graphsFine,
              const Chombo4::DisjointBoxLayout                                  &  a_gridsFine,
              const Chombo4::Box                                                &  a_domFine,
              const Real                                                        &  a_dxFine,
              const shared_ptr<graph_distrib_t    >                             &  a_graphsCoar,
              const Chombo4::DisjointBoxLayout                                  &  a_gridsCoar,
              const Chombo4::Box                                                &  a_domCoar,
              const Real                                                        &  a_dxCoar,
              const shared_ptr<EBDictionary<HOEB_MAX_ORDER, Real, CELL, CELL> > &  a_dictionary,
              const shared_ptr< GeometryService<HOEB_MAX_ORDER> >               &  a_geoserv);

  ///
  /**
   **/
  template<CENTERING cent>
  void
  restrictFlux(EBLevelBoxData<cent, 1>                                          &  a_phiFToC,
              const EBLevelBoxData<cent, 1>                                     &  a_phiFine,
              const shared_ptr<graph_distrib_t    >                             &  a_graphsFine,
              const Chombo4::DisjointBoxLayout                                  &  a_gridsFine,
              const Chombo4::Box                                                &  a_domFine,
              const Real                                                        &  a_dxFine,
              const shared_ptr<graph_distrib_t    >                             &  a_graphsCoar,
              const Chombo4::DisjointBoxLayout                                  &  a_gridsCoar,
              const Chombo4::Box                                                &  a_domCoar,
              const Real                                                        &  a_dxCoar,
              const shared_ptr<EBDictionary<HOEB_MAX_ORDER, Real, cent, cent> > &  a_dictionary,
               const shared_ptr< GeometryService<HOEB_MAX_ORDER> >               &  a_geoserv);
  
  
}

#endif
