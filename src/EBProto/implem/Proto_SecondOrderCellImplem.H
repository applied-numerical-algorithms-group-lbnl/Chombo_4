                                                                                                                                  
#ifndef _Proto_SecondOrderCellImplem_H_
#define _Proto_SecondOrderCellImplem_H_
#include "Proto_Normalizor.H"

namespace Proto
{
  
  /*****************************/  
  template<CENTERING cent, int order>
  void
  Base2ndOrderFlux<cent, order>::
  getFaceInterpStencil(LocalStencil<cent, double>                   &  a_sten,
                       const EBIndex<cent>                          &  a_face,
                       const EBGraph                                &  a_graph,
                       const HostIrregData<cent, IndMomSDMinOne, 1> &  a_facdata,
                       double                                          a_dx) const
  {
    a_sten.clear();
    double areaFrac;
    RealVect centroid;
    int  facedir = a_face.direction();
    int tandirs[DIM-1];
    int itan = 0;
    for(int idir = 0; idir < DIM; idir++)
    {
      if(idir != facedir)
      {
        tandirs[itan] = idir;
        itan++; 
      }
    }

    Normalizor<order>::normalizedFaceMoments(areaFrac, centroid, a_facdata, a_face, a_graph, a_dx);

    Box domain = a_graph.getDomain();
    bool dropOrder = false;
    Side::LoHiSide sides[DIM-1];
    for(int itan = 0; itan < DIM-1; itan++)
    {
      if(centroid[tandirs[itan]] > 0)
      {
        sides[itan] = Side::Hi;
      }
      else
      {
        sides[itan] = Side::Lo;
      }
    }

#if DIM==2    
    EBIndex<cent> faces[2];
    faces[0] = a_face;
    bool hasFace = a_graph.getAdjacentFace<cent>(faces[1], a_face,tandirs[0], sides[0]);      
    if(!hasFace) dropOrder = true;
//debug
//    dropOrder = true;
//end debug
    if(dropOrder)
    {
      a_sten.add(a_face, 1.0);
    }
    else
    {
      double xbar = std::abs(centroid[tandirs[0]]);
      a_sten.add(faces[0], 1.0-xbar);
      a_sten.add(faces[1],     xbar);
    }
#else
    EBIndex<cent> faces[4];
    faces[0] = a_face;
    double xbar = std::abs(centroid[tandirs[0]]);
    double ybar = std::abs(centroid[tandirs[1]]);
    bool hasFacex = a_graph.getAdjacentFace<cent>(faces[1], a_face, tandirs[0], sides[0]);      
    bool hasFacey = a_graph.getAdjacentFace<cent>(faces[2], a_face, tandirs[1], sides[1]);      
    if((!hasFacex)  || (!hasFacey))
    {
      dropOrder = true;
    }
    if(dropOrder)
    {
      a_sten.add(a_face, 1.0);
    }
    else
    {
      EBIndex<cent> facexy, faceyx;
      bool hasFaceyx = a_graph.getAdjacentFace<cent>(faceyx, faces[1], tandirs[1], sides[1]);      
      bool hasFacexy = a_graph.getAdjacentFace<cent>(facexy, faces[2], tandirs[0], sides[0]);      
      if((!hasFaceyx) || (!hasFacexy) || (facexy != faceyx))
      {
        dropOrder = true;
        a_sten.add(a_face, 1.0);
      }
      else
      {
        faces[3] = facexy;
        double  f00coef = 1.0 - xbar - ybar + xbar*ybar;
        double  f10coef = xbar - xbar*ybar;
        double  f01coef = ybar - xbar*ybar;
        double  f11coef = xbar*ybar;
        a_sten.add(faces[0], f00coef);
        a_sten.add(faces[1], f10coef);
        a_sten.add(faces[2], f01coef);
        a_sten.add(faces[3], f11coef);
      }
    }
#endif
  }
/*******************************/
  template<CENTERING cent, int order>
  void
  Base2ndOrderFlux<cent, order>::
  getCentroidStencil(LocalStencil<CELL, double>                   &  a_sten,
                     const EBIndex<cent>                          &  a_face,
                     const EBGraph                                &  a_graph,
                     const HostIrregData<cent, IndMomSDMinOne, 1> &  a_facdata,
                     string                                          a_dombc,
                     double                                          a_dx) const
  {
#if 1
    LocalStencil<cent, double> facesten;
    getFaceInterpStencil(facesten, a_face,
                         a_graph,    
                         a_facdata, 
                         a_dx);   
#else
    LocalStencil<cent, double> facesten;
    facesten.add(a_face, 1.0);
#endif 
    a_sten.clear();
    for(int isten = 0; isten < facesten.size(); isten++)
    {
      EBIndex<cent>& face   = facesten.m_entries[isten].m_vof;
      double       & weight = facesten.m_entries[isten].m_weight;
      LocalStencil<CELL, double> centSten;
      getFaceCenteredFluxStencil(centSten, face, a_graph, a_facdata, a_dombc, a_dx);
      centSten *= weight;
      a_sten   += centSten;
    }
  }

           
///
  template <CENTERING cent, int order>
  void
  Poisson2ndOrder<cent,order>::
  getFaceCenteredFluxStencil(LocalStencil<CELL, double>                   &  a_facesten,
                             const EBIndex<cent>                          &  a_face,
                             const EBGraph                                &  a_graph,
                             const HostIrregData<cent, IndMomSDMinOne, 1> &  a_facdata,
                             string                                          a_dombc,
                             double                                          a_dx) const
  {
    EBIndex<CELL> voflo = a_face.getVoF(Side::Lo);
    EBIndex<CELL> vofhi = a_face.getVoF(Side::Hi);
    a_facesten.clear();
    int ideb = 0;
    //if(true)
    if(!a_face.m_isBoundary)
    {
      a_facesten.add(vofhi,  1./a_dx);
      a_facesten.add(voflo, -1./a_dx);
    }
    else if(a_dombc == string("Neumann"))
    {
      //zero flux here
      ideb++;
    }
    else if(a_dombc == string("Dirichlet"))
    {
      //figure out low or high side
      if(a_face.m_vofIDLo < 0)
      {
        a_facesten.add(vofhi,  2./a_dx);
      }
      else if(a_face.m_vofIDMe < 0)
      {
        a_facesten.add(voflo,  -2./a_dx);
      }
      else
      {
        PROTO_ASSERT(false, "logic error in getfacecenteredfluxstencil");
      }
      
    }
    else
    {
      PROTO_ASSERT(false, "bc not found in getfacecenteredfluxstencil");
    }
    
            
  }

  template < int order>
  void
  SecondOrderStencil<order>::
  get2ndOrderDivFStencil(LocalStencil< CELL,  double>          & a_vofStencil,
                         const EBIndex<CELL>                   & a_vof,
                         const EBGraph                         & a_graph,
                         const VoluData                        & a_voludata,
                         const EBFaData                        & a_ebfadata,
                         const XFacData                        & a_xfacdata,
                         const YFacData                        & a_yfacdata,
                         const ZFacData                        & a_zfacdata,
                         const Base2ndOrderFlux<XFACE, order>  & a_xfaceflux,
                         const Base2ndOrderFlux<YFACE, order>  & a_yfaceflux,
                         const Base2ndOrderFlux<ZFACE, order>  & a_zfaceflux,
                         string                                  a_dombc[2*DIM],
                         string                                  a_ebbcName,
                         double                                  a_dx)
  {
    PR_TIME("get2ndOrderDivFStencil");
    a_vofStencil.clear();
#if DIM==2
//    int ideb = 0;
//    Box domain = a_graph.getDomain();
//    Point ivdeb(2, 9);
//    if((a_vof.m_pt == ivdeb) && (domain.size(0) == 32))
//    {
//      ideb = 1;
//    }
#endif

    for (SideIterator sit; sit.ok(); ++sit)
    {
      vector< EBIndex<XFACE> > xfaces = a_graph.getXFaces(a_vof, sit());
      int isign = sign(sit());
      for (int iface = 0; iface < xfaces.size(); iface++)
      {
        LocalStencil<CELL, double> fluxStencil;
        double   areaFrac;
        RealVect centroid;
        Normalizor<order>::normalizedFaceMoments(areaFrac, centroid, a_xfacdata, xfaces[iface], a_graph, a_dx);

        int iindex = ebp_index(0, sit());
        a_xfaceflux.getCentroidStencil(fluxStencil, xfaces[iface],
                                       a_graph,    
                                       a_xfacdata,
                                       a_dombc[iindex],
                                       a_dx);       
        double factor = double(isign)*areaFrac/a_dx;
        fluxStencil *= factor;
        a_vofStencil += fluxStencil;
      }
    }

    for (SideIterator sit; sit.ok(); ++sit)
    {
      vector< EBIndex<YFACE> > yfaces = a_graph.getYFaces(a_vof, sit());
      int isign = sign(sit());
      for (int iface = 0; iface < yfaces.size(); iface++)
      {
        LocalStencil<CELL, double> fluxStencil;
        double   areaFrac;
        RealVect centroid;
        Normalizor<order>::normalizedFaceMoments(areaFrac, centroid, a_yfacdata, yfaces[iface], a_graph, a_dx);
        int iindex = ebp_index(1, sit());
        a_yfaceflux.getCentroidStencil(fluxStencil, yfaces[iface],
                                       a_graph,    
                                       a_yfacdata,
                                       a_dombc[iindex],
                                       a_dx);       

        fluxStencil *= double(isign)*areaFrac/a_dx;
        a_vofStencil += fluxStencil;
      }
    }
#if DIM==3
    for (SideIterator sit; sit.ok(); ++sit)
    {
      vector< EBIndex<ZFACE> > zfaces = a_graph.getZFaces(a_vof, sit());
      int isign = sign(sit());
      for (int iface = 0; iface < zfaces.size(); iface++)
      {
        LocalStencil<CELL, double> fluxStencil;
        double   areaFrac;
        RealVect centroid;
        Normalizor<order>::normalizedFaceMoments(areaFrac, centroid, a_zfacdata, zfaces[iface], a_graph, a_dx);

        int iindex = ebp_index(1, sit());
        a_zfaceflux.getCentroidStencil(fluxStencil, zfaces[iface],
                                       a_graph,    
                                       a_zfacdata, 
                                       a_dombc[iindex],
                                       a_dx);       


        fluxStencil *= double(isign)*areaFrac/a_dx;
        a_vofStencil += fluxStencil;
      }
    }
#endif      

//    if(ideb == 1)
//    {
//      cout << "before ebflux, ";
//      cout << "stencil for " << a_vof.m_pt << ":" << endl;
//      a_vofStencil.print();
//    }

    LocalStencil<CELL, double> ebflux;
    EBIndex<BOUNDARY> face = a_vof.getCutFace();
    double weight; //for inhomogeneous bcs
    if(a_ebbcName == string("Neumann"))
    {
      NeumannEBFlux<order>::
        getEBFluxStencil(ebflux, weight,
                         face,
                         a_graph,
                         a_voludata,
                         a_ebfadata,
                         a_xfacdata,
                         a_yfacdata,
                         a_zfacdata,
                         a_dx);
    }
    else if(a_ebbcName == string("Dirichlet"))
    {
      DirichletEBFlux<order>::
        getEBFluxStencil(ebflux, weight,
                         face,
                         a_graph,
                         a_voludata,
                         a_ebfadata,
                         a_xfacdata,
                         a_yfacdata,
                         a_zfacdata,
                         a_dx);
//      if(ideb==1)
//      {
//        cout << "ebflux stencil before multiply =" << endl;
//        ebflux.print();
//      }
    }
    else
    {
      PROTO_ASSERT(false, "could not find EBBC string");
    }
      
    RealVect centroid, normal;
    double bndryArea, alpha;
    EBIndex<BOUNDARY> cutface = a_vof.getCutFace();
    Normalizor<order>::normalizedBndryMoments(alpha, centroid, a_ebfadata,  cutface, a_graph, a_dx);
    Normalizor<order>::getNormalAndBoundaryArea(bndryArea, normal,   
                                                a_voludata,
                                                a_ebfadata,
                                                a_xfacdata,
                                                a_yfacdata,
                                                a_zfacdata, 
                                                a_vof, a_graph, a_dx);
    //need the boundary area that preserves free stream instead of the moments version
    alpha = bndryArea;
    ebflux *= alpha/a_dx;

//    if(ideb == 1)
//    {
//      cout << "ebflux after multiply " << endl;
//      ebflux.print();
//    }

    //for some reason, EBAMRPoissonOp calculates beta weight without the EB flux
    double diagweight;
    bool foundVoF = a_vofStencil.getValueAtVoF(diagweight, a_vof);
    PR_assert(foundVoF); //if this fails, the stencil has no diagonal value
    a_vofStencil.setBetaWeight(diagweight);

//    if(ideb == 1)
//    {
//      cout << "betaweight = " << diagweight << endl;
//    }
    
    a_vofStencil += ebflux;

//    if(ideb == 1)
//    {
//      cout << "after ebflux, ";
//      cout << "stencil for " << a_vof.m_pt << ":" << endl;
//      a_vofStencil.print();
//    }

  }

  template < int order>
  void
  DirichletEBFlux<order>::
  calculateWeightingMatrix(RealVect           x0,
                           vector<RealVect>&  xp,
                           vector<RealVect>&  weightMatrix,
                           bool&              detZero)
  {
    int stenSize = xp.size();

    vector<RealVect> deltaX = xp;
    for (int isten = 0; isten < stenSize; isten++)
    {
      deltaX[isten] -= x0;
    }

    vector<RealVect>    aTransA(DIM, RealVect::Zero());
    vector<RealVect> invATransA(DIM, RealVect::Zero());
    for (int idir = 0; idir < DIM; idir++)
    {
      for (int jdir = 0; jdir < DIM; jdir++)
      {
        for (int isten = 0; isten < stenSize; isten++)
        {
          aTransA[idir][jdir] = aTransA[idir][jdir]
            + deltaX[isten][idir]*deltaX[isten][jdir];
        }
      }
    }

    double det;
#if DIM == 2
    det = aTransA[0][0] * aTransA[1][1] - aTransA[0][1] * aTransA[1][0];
    if (det < 1.e-15 && det > -1.e-15)
    {
      detZero = true;
    }
    else
    {
      invATransA[0][0] =  aTransA[1][1] / det;
      invATransA[0][1] = -aTransA[0][1] / det;
      invATransA[1][0] = -aTransA[1][0] / det;
      invATransA[1][1] =  aTransA[0][0] / det;
    }
#elif DIM == 3
    det = aTransA[0][0] * ( aTransA[1][1] * aTransA[2][2]
                            - aTransA[1][2] * aTransA[2][1])
      + aTransA[0][1] * ( aTransA[1][2] * aTransA[2][0]
                          - aTransA[1][0] * aTransA[2][2])
      + aTransA[0][2] * ( aTransA[1][0] * aTransA[2][1]
                          - aTransA[1][1] * aTransA[2][0]);

    if (det < 1.e-15 && det > -1.e-15)
    {
      detZero = true;
    }
    else
    {
      invATransA[0][0] = ( aTransA[1][1] * aTransA[2][2]
                           - aTransA[1][2] * aTransA[2][1]) / det;
      invATransA[0][1] = ( aTransA[1][2] * aTransA[2][0]
                           - aTransA[1][0] * aTransA[2][2]) / det;
      invATransA[0][2] = ( aTransA[1][0] * aTransA[2][1]
                           - aTransA[1][1] * aTransA[2][0]) / det;
      invATransA[1][0] = ( aTransA[2][1] * aTransA[0][2]
                           - aTransA[2][2] * aTransA[0][1]) / det;
      invATransA[1][1] = ( aTransA[2][2] * aTransA[0][0]
                           - aTransA[2][0] * aTransA[0][2]) / det;
      invATransA[1][2] = ( aTransA[2][0] * aTransA[0][1]
                           - aTransA[2][1] * aTransA[0][0]) / det;
      invATransA[2][0] = ( aTransA[0][1] * aTransA[1][2]
                           - aTransA[0][2] * aTransA[1][1]) / det;
      invATransA[2][1] = ( aTransA[0][2] * aTransA[1][0]
                           - aTransA[0][0] * aTransA[1][2]) / det;
      invATransA[2][2] = ( aTransA[0][0] * aTransA[1][1]
                           - aTransA[0][1] * aTransA[1][0]) / det;
    }
#else
    THIS_IS_AN_ERROR_MESSAGE__THIS_WILL_ONLY_COMPILE_WHEN_DIM_IS_2_OR_3;
#endif

    //if (!detZero)
    {
      weightMatrix = vector<RealVect>(stenSize,RealVect::Zero());
      for (int idir = 0; idir < DIM; idir++)
      {
        for (int isten = 0; isten < stenSize; isten++)
        {
          for (int jdir = 0; jdir < DIM; jdir++)
          {
            weightMatrix[isten][idir] += invATransA[idir][jdir] * deltaX[isten][jdir];
          }
        }
      }
    }
  }
  
  /***/
  template < int order>
  void
  DirichletEBFlux<order>::
  getEBFluxStencil(LocalStencil<CELL, double>    &  a_stencil,
                   double                        &  a_weight, //for inhomogeneous bcs
                   const EBIndex<BOUNDARY>       &  a_face,
                   const EBGraph                 & a_graph,
                   const VoluData                & a_voludata,
                   const EBFaData                & a_ebfadata,
                   const XFacData                & a_xfacdata,
                   const YFacData                & a_yfacdata,
                   const ZFacData                & a_zfacdata,
                   const double                  & a_dx)
  {
    EBIndex<CELL> vof = a_face.getVoF();
    unsigned int radius = 2;
    vector< EBIndex<CELL> > volsinclusive=  a_graph.getAllVoFsWithinRadius(vof.m_pt, radius);
    vector<EBIndex<CELL> >  vols;
    for(int ivof = 0; ivof < volsinclusive.size(); ivof++)
    {
      if(volsinclusive[ivof] != vof)
      {
        vols.push_back(volsinclusive[ivof]);
      }
    }
    RealVect centroid, normal;
    double bndryArea, alpha;
    Normalizor<order>::normalizedBndryMoments(alpha, centroid, a_ebfadata, a_face, a_graph, a_dx);
    Normalizor<order>::getNormalAndBoundaryArea(bndryArea, normal,   
                                                a_voludata,
                                                a_ebfadata,
                                                a_xfacdata,
                                                a_yfacdata,
                                                a_zfacdata, 
                                                vof, a_graph, a_dx);

#if DIM==2
    unsigned int minStenSize = 3;
#else 
    unsigned int minStenSize = 7;
#endif

    bool dropOrder= (vols.size() < minStenSize);

    
    if (!dropOrder)
    {
////begin debug
//#if DIM==2
//      int ideb = 0;
//      int domlen = a_graph.getDomain().size(0);
//      Point ivdeb(0, 9);
//      if((vof.m_pt == ivdeb) && (domlen == 16))
//      {
//        ideb = 1;
//      }
//#endif      
////end debug      
      RealVect x0;
      Point iv0 = vof.m_pt;
      for (int idir = 0; idir < DIM; idir++)
      {
        x0[idir] = a_dx * (0.5 + centroid[idir] + iv0[idir]);
//        x0[idir] = a_dx * (0.5 + iv0[idir]);
      }

      vector<RealVect> xp(vols.size());
      for (int isten = 0; isten < vols.size(); isten++)
      {
        for (int idir = 0; idir < DIM; idir++)
        {
          xp[isten][idir] = a_dx * (0.5 + vols[isten].m_pt[idir]);
        }
      }

      unsigned int volsize = vols.size();
      vector<RealVect> invATransAdeltaX(volsize, RealVect::Zero());


      bool detZero = false;
      calculateWeightingMatrix(x0, xp, invATransAdeltaX, detZero);

      a_stencil.clear();
      a_weight = 0.0;

      for (int isten = 0; isten < vols.size(); isten++)
      {
        double dphidnWeight = 0.0;
        for (int idir = 0; idir < DIM; idir++)
        {
          dphidnWeight -= invATransAdeltaX[isten][idir] * normal[idir];
        }
              
        a_stencil.add(vols[isten],dphidnWeight);
        a_weight -= dphidnWeight;
      }
    }
    else
    {
      a_stencil.clear();
      a_weight = 0.0;
    }
    
  }
/**/
  inline double getStenSign(string a_domainBCName)
  {
    double stenSign;
    if(a_domainBCName == string("Dirichlet"))
    {
      stenSign = -1;
    }
    else if(a_domainBCName == string("Neumann"))
    {
      stenSign =  1;
    }
    else
    {
      PROTO_ASSERT(false, "unknown bc name");
    }
    return stenSign;
  }
  /***/
  template < int order>
  void
  SecondOrderStencil<order>::
  getDomainBCs(vector<Stencil<double> >    a_regStencilBC[2*DIM],       
               vector<Box>                 a_BCApplyBoxes[2*DIM],       
               string                      a_domainBCName[2*DIM],               
               const Box                 & a_validBox,
               const Box                 & a_domain,
               const Point               & a_srcGhost)
    
  {
    Point lobox = a_validBox.low();
    Point hibox = a_validBox.high();
    Point lodom = a_domain.low();
    Point hidom = a_domain.high();
    Box srcBox = a_validBox.grow(a_srcGhost);

    for(int idir = 0; idir < DIM; idir++)
    {
      int numghost = a_srcGhost[idir]; //need a bc stencil for each ghost cell
      Point loapp  = srcBox.low();
      Point hiapp  = srcBox.high();
      loapp[idir]  = lodom[idir]-1;
      hiapp[idir]  = lodom[idir]-1;

      if(lobox[idir] == lodom[idir])
      {
        int loindex = ebp_index(idir, Side::Lo);
        a_regStencilBC[loindex].resize(numghost);
        a_BCApplyBoxes[loindex].resize(numghost);
        string bcname = a_domainBCName[loindex];
        if((bcname == string("Dirichlet")) || (bcname == string("Neumann")))
        {
          double  stenSign = getStenSign(bcname);
          Box appBox0(loapp, hiapp);

          for(int igh = 0; igh < numghost; igh++)
          {
            a_BCApplyBoxes[loindex][igh] = appBox0.shift(idir, -igh);
            a_regStencilBC[loindex][igh] = (stenSign)*Shift(Point::Basis(idir, 2*igh+1));
          }
        }
        else if(bcname == string("LinearExtrapolation"))
        {
          //just going to do one ghost cell as extrapolation.  the rest are piecewise constant
          Box appBox0(loapp, hiapp);
          a_BCApplyBoxes[loindex][0] = appBox0;
          a_regStencilBC[loindex][0] = (2.)*Shift(Point::Basis(idir, 1)) + (-1.)*Shift(Point::Basis(idir, 2));

          for(int igh = 1; igh < numghost; igh++)
          {
            a_BCApplyBoxes[loindex][igh] = appBox0.shift(idir, -igh);
            a_regStencilBC[loindex][igh] = (1.)*Shift(Point::Basis(idir, 1));
          }
        }
        else if(bcname == string("no_bcs"))
        {
          //do nothing
        }
        else
        {
          std::cout << "bogus bc name = " << bcname << std::endl;
          PROTO_ASSERT(false, "getsecondorderbcs: bogus bc string lo");
        }
      }
      if(hibox[idir] == hidom[idir])
      {
        int hiindex = ebp_index(idir, Side::Hi);
        a_regStencilBC[hiindex].resize(numghost);
        a_BCApplyBoxes[hiindex].resize(numghost);
        string bcname = a_domainBCName[hiindex];
        Point loapp = srcBox.low();
        Point hiapp = srcBox.high();
        loapp[idir] = hidom[idir]+1;
        hiapp[idir] = hidom[idir]+1;

        if((bcname == string("Dirichlet")) || (bcname == string("Neumann")))
        {
          double  stenSign = getStenSign(bcname);
          Box appBox0(loapp, hiapp);
        
          for(int igh = 0; igh < numghost; igh++)
          {
            a_BCApplyBoxes[hiindex][igh] = appBox0.shift(idir,  igh);
            a_regStencilBC[hiindex][igh] = (stenSign)*Shift(Point::Basis(idir, -2*igh-1));
          }
        }
        else if(bcname == string("LinearExtrapolation"))
        {
          //just going to do one ghost cell as extrapolation.  the rest are piecewise constant
          Box appBox0(loapp, hiapp);
          a_BCApplyBoxes[hiindex][0] = appBox0;
          a_regStencilBC[hiindex][0] = (2.)*Shift(Point::Basis(idir, -1)) + (-1.)*Shift(Point::Basis(idir, -2));

          for(int igh = 1; igh < numghost; igh++)
          {
            a_BCApplyBoxes[hiindex][igh] = appBox0.shift(idir, igh);
            a_regStencilBC[hiindex][igh] = (1.)*Shift(Point::Basis(idir, -1));
          }
        }
        else if(bcname == string("no_bcs"))
        {
          //do nothing
        }
        else
        {
          std::cout << "bogus bc name = " << bcname << std::endl;
          PROTO_ASSERT(false, "getsecondorderbcs: bogus bc string hi");
        }
      }
    }
  }
}
#endif
