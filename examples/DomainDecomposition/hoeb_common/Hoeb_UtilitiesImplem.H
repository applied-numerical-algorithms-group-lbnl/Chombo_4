#ifndef __hoeb_utilities_implem_
#define __hoeb_utilities_implem_
#include "Hoeb_Neighborhood.H"
namespace hoeb
{
/***/
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
/***/
  template <CENTERING cent>
  Real
  getDevendranFluxFace(const EBHostData<CELL, Real, 1>                     & a_phi,
                       const EBGraph                                       & a_graph,
                       const EBIndex<cent>                                 & a_face,
                       const EBIndex<CELL>                                 & a_startingVoF,
                       const std::string                                     a_dombcname[2*DIM],
                       const std::string                                   & a_ebbcname,
                       const shared_ptr< GeometryService<HOEB_MAX_ORDER> > & a_geoserv,
                       const Chombo4::Box                                  & a_domain,
                       Real a_dx, int a_facedir,
                       Side::LoHiSide a_vofside, const Chombo4::DataIndex& a_dit, bool a_divideByArea,                       bool a_printStuff)
  {
    Real retval = 0;
    shared_ptr<hoeb::BaseExactSolution<HOEB_MAX_ORDER> > exactSol =  hoeb::getBaseExactSoltuion();
    unsigned int bindex = ebp_index(a_facedir, a_vofside);
    bool isNeumannBoundaryFace = (a_face.m_isBoundary && (a_dombcname[bindex] == string("Neumann")));
    if(isNeumannBoundaryFace)
    {
      retval =  getNeumannBoundaryFlux<cent>(exactSol,
                                             a_phi,
                                             a_graph,
                                             a_face,
                                             a_startingVoF,
                                             a_dombcname,
                                             a_ebbcname,
                                             a_geoserv,
                                             a_domain,
                                             a_dx, a_facedir,
                                             a_vofside, a_dit, a_divideByArea);
    }
    else
    {
      hoeb::CompositeStencil composite =
        hoeb::getDevendranIntFluxDACompositeStencil(a_face,
                                                    a_startingVoF,
                                                    a_dombcname,
                                                    a_ebbcname,
                                                    a_geoserv,
                                                    a_domain,
                                                    a_dit,                              
                                                    a_dx, a_facedir,
                                                    a_vofside,
                                                    a_printStuff);



      retval= composite.evaluate(a_phi, *exactSol, a_graph, a_geoserv, a_domain, a_dx, a_dit);
      if(a_divideByArea)
      {
        //the stencil evaluates to int(flux da) so we have to divide out the area
        IndexedMoments<DIM-1, HOEB_MAX_ORDER> facemom;
        if(a_graph.isFaceRegular(a_face))
        {
          facemom.setToRegular(a_dx);
        }
        else
        {
          
          EBIndex<cent> dummy;
          typedef CH4_Data_Choreography::DistributedData<HostIrregData<cent, IndexedMoments<DIM-1, HOEB_MAX_ORDER>, 1 > > face_distrib_t;
          shared_ptr< face_distrib_t > geodataptr = getFaceData<cent>(a_geoserv, a_domain, dummy);
          facemom = (*geodataptr)[a_dit](a_face, 0);
        }
        Real area = facemom[IndexTM<int, DIM-1>::Zero];
        if(area > 1.0e-10)
        {
          retval *= (1.0/area);
        }
      }
    }
    
    return retval;
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
        matval = pd*derval;
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
  /**
     Return integral F dA  (do NOT divide out the area).
     For EB Fluxes, it is Integral F-n dA with an outward-facing normal.
     starting vof never used
  */
  template <CENTERING cent>
  CompositeStencil
  getDevendranIntFluxDACompositeStencil(const EBIndex<cent>                                 & a_face,
                                        const EBIndex<CELL>                                 & a_startingVoF,
                                        const std::string                                     a_dombcname[2*DIM],
                                        const std::string                                   & a_ebbcname,
                                        const shared_ptr< GeometryService<HOEB_MAX_ORDER> > & a_geoserv,
                                        const Chombo4::Box                                  & a_srcDomain,
                                        const Chombo4::DataIndex& a_dit,
                                        Real a_dx, int a_facedir,
                                        Side::LoHiSide a_vofside,
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
      //get the neighborhood near the target face
      Neighborhood<CELL>  block =
        getCellNeighbors<cent>(a_face, graph, a_dombcname, a_ebbcname, a_facedir);
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
      
      LAPACKMatrix Smat = getStencilMatrix(Qmat, Gmat, GTGinvMat, weight);
      
      retval = hoeb::getStencilFromMatrix(Smat, block);
      
      if(a_printStuff)
      {
        Chombo4::pout() << "composite stencil = " << endl;
        retval.poutAll();
      }
    }
    return retval;
  }

}
#endif
