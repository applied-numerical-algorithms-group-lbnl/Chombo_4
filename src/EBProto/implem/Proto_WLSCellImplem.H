                                                                                                                                   
#ifndef _Proto_WLSCellImplem_H_
#define _Proto_WLSCellImplem_H_

namespace Proto
{
  
  template<CENTERING cent, int order>
  Real
  BaseWLSFlux<cent, order>::
  computeInhomogTerm(const Chombo4::LAPACKMatrix    & a_stencil,
                     const Vector<IndMomSDMinOne>   & a_facmom,
                     const Vector<IndMomDIM>        & a_ebfmom,
                     const int                        a_nvolneigh,
		     const Vector<RealVect>         & a_facneighLoc,
		     const Vector<int>              & a_facneighDir,
                     const Vector<Side::LoHiSide>   & a_facneighSid,
		     const Vector<RealVect>         & a_ebfneighLoc,
	             const Vector<Vector<IndMomDIM> >& a_ebfneighNorm,
		     string                           a_ebbcName,
                     string                           a_dombcName[2*DIM],
                     double                           a_dx) const
  {
    Vector<Real> domBCValues, ebBCValues;
    int nrows = a_facmom.size();
    domBCValues.resize(nrows);
    for(int iface = 0; iface < nrows; iface++)
    {
      RealVect a_faceLoc = a_facneighLoc[iface];
      int a_faceDir = a_facneighDir[iface];
      Side::LoHiSide a_faceSid = a_facneighSid[iface];
      int index = ebp_index(a_faceDir, a_faceSid);
      if(a_dombcName[index] == string("Neumann"))
      {
        domBCValues[iface] = NeumannExactSoln<order>::
                               getFaceIntegral(a_facmom[iface],
					       a_faceLoc,
					       a_faceDir,
					       a_dx);
      }
      else if(a_dombcName[index] == string("Dirichlet"))
      {
        domBCValues[iface] = DirichletExactSoln<order>::
                               getFaceIntegral(a_facmom[iface],
					       a_faceLoc,
					       a_faceDir,
					       a_dx);
      }
      else
      {
        PROTO_ASSERT(false, "could not find domBC string");
      }

      //domBCValues[iface] = getFaceIntegral(a_faces[iface], a_ebisBox, a_dx);
      //bool divideable = divideByArea(a_aveValues[iface], a_faces[iface], a_ebisBox, a_dx);
      //if(!divideable)
      //{
      //  domBCValues[iface] = 0.;
      //}
    }

    nrows = a_ebfmom.size();
    ebBCValues.resize(nrows);
    for(int ivof = 0; ivof < nrows; ivof++)
    {
      RealVect a_ebfLoc = a_ebfneighLoc[ivof];
      if(a_ebbcName == string("Neumann"))
      {
        ebBCValues[ivof] = NeumannExactSoln<order>::
                               getEBIntegral(a_ebfneighNorm[ivof],
				             a_ebfLoc,
				             a_dx);
      }
      else if(a_ebbcName == string("Dirichlet"))
      {
        ebBCValues[ivof] = DirichletExactSoln<order>::
                               getEBIntegral(a_ebfmom[ivof],
				             a_ebfLoc,
				             a_dx);
      }
      else
      {
        PROTO_ASSERT(false, "could not find domBC string");
      }

      //bool divideable = divideByArea(a_aveValues[iface], a_faces[iface], a_ebisBox, a_dx);
      //if(!divideable)
      //{
      //  ebBCValues[iface] = 0.;
      //}
    }
  
    int nfac = a_facmom.size();
    int nebf = a_ebfmom.size();
  
    Real inhomogTerm = 0.0;
    for(int iface = 0; iface < nfac; iface++)
    {
      inhomogTerm += domBCValues[iface]*a_stencil(0, a_nvolneigh + iface);
    }
    for(int iface = 0; iface < nebf; iface++)
    {
      inhomogTerm += ebBCValues[iface]*a_stencil(0, a_nvolneigh + nfac + iface);
    }

    return inhomogTerm;
  }

  template<CENTERING cent, int order>
  void
  BaseWLSFlux<cent, order>::
  getEBFluxMatrix(Chombo4::LAPACKMatrix    & a_fluxMatrix,
                  const EBIndex<cent>      & a_face,
		  RealVect                 & a_faceLoc,
		  vector<IndMomDIM>        & a_ebnorm,
		  double                   & a_bndryArea,
                  double                     a_dx) const
  {
    Real areaScaleFactor  = (D_TERM(1.0,  *a_dx, *a_dx));
    MomentIterator<DIM, order> momit;
    a_fluxMatrix.define(1, momit.size());
    a_fluxMatrix.setVal(0.);

    RealVect startLoc;
    for (int idir = 0; idir < DIM; idir++)
    {
      startLoc[idir] = a_dx * (0.5 + a_face.m_pt[idir]);
      if (a_face.direction() == idir) startLoc[idir] = a_dx * (a_face.m_pt[idir]);
    }

    for(int idir = 0; idir < DIM; idir++)
    {

      RealVect diff = startLoc;
      diff -= a_faceLoc;
      IndexTM<Real, DIM> diffrv;
      for(int idir = 0; idir < DIM; idir++)
      {
        diffrv[idir] = diff[idir];
      }
      a_ebnorm[idir].shift(diffrv);

      // the contribution to the flux matrix from the normal in this particular direction
      Chombo4::LAPACKMatrix fluxIncr(1, momit.size());
      fluxIncr.setVal(0.);
      for(momit.reset(); momit.ok(); ++momit)
      {
        IndexTM<int, DIM> p = momit();
        Real pd = p[idir];
        //get column number associated with this particular moment
        int jmom = IndexedMoments<DIM, order>::indexOf(p);
        if(pd > 0)
        {
          IndexTM<int, DIM> pshift = p - BASISV_TM<int, DIM>(idir);
          fluxIncr(0, jmom) = pd*a_ebnorm[idir][pshift];///(a_bndryArea*areaScaleFactor*POW(a_dx, momit().sum()));
	  cout<<" bndryArea "<<a_bndryArea<<endl;
        }
        else
        {
          fluxIncr(0, jmom) = 0;
        }
      }
      a_fluxMatrix += fluxIncr;
    }
    for(int idir = 0; idir < DIM; idir++)
    {
      for(momit.reset(); momit.ok(); ++momit)
      {
        IndexTM<int, DIM> p = momit();
        Real pd = p[idir];
        //get column number associated with this particular moment
        int jmom = IndexedMoments<DIM, order>::indexOf(p);
        ::pout() << ",icol = "  << jmom <<  ", Q(0,icol) = " << a_fluxMatrix(0, jmom) << endl;
      }
    }
   
  }

  template<CENTERING cent, int order>
  void
  BaseWLSFlux<cent, order>::
  getFaceFluxMatrix(Chombo4::LAPACKMatrix    & a_fluxMatrix,
                    const EBIndex<cent>      & a_face,
		    RealVect                 & a_faceLoc,
	            IndMomSDMinOne           & a_faceMom,
                    double                     a_dx) const
  {
    Real areaScaleFactor  = (D_TERM(1.0,  *a_dx, *a_dx));

    RealVect startLoc;
    for (int idir = 0; idir < DIM; idir++)
    {
      startLoc[idir] = a_dx * (0.5 + a_face.m_pt[idir]);
      if (a_face.direction() == idir) startLoc[idir] = a_dx * (a_face.m_pt[idir]);
    }
    // shift center of moments to a_faceLoc.
    // the moments are centered at the face center initially
    RealVect diff = startLoc;
    diff -= a_faceLoc;
    IndexTM<Real, DIM-1> diffrv;
    int index = 0;
    for(int idir = 0; idir < DIM; idir++)
    {
      if(idir != a_face.direction())
      {
        diffrv[index] = diff[idir];
        index++;
      }
    }
    a_faceMom.shift(diffrv);

    MomentIterator<DIM, order> momit;
    a_fluxMatrix.define(1, momit.size());
    const int& faceDir = a_face.direction();
    for(momit.reset(); momit.ok(); ++momit)
    {
      IndexTM<int, DIM> p = momit();
      int pd = p[faceDir];
    
      int jmom = IndexedMoments<DIM, order>::indexOf(p);

      //all other integrals are zero since along face
      //(y-ybar) =0  so (y-ybar)^py = 0 unless py == 0.
      if(pd == 1)
      {
        // pfaceshift are all components of p that are not in the faceDir direction
        IndexTM<int, DIM-1> pfaceshift;
        int index = 0;
        for(int idir = 0; idir < DIM; idir++)
        {
           if(idir  != faceDir)
           {
             pfaceshift[index] =  p[idir];
             index++;
           }
        }
        // for readability, multiply by pd even though pd = 1
        //also pfaceshift is really p-e^d because we have set p
        a_fluxMatrix(0, jmom) = pd*a_faceMom[pfaceshift];///areaScaleFactor;
      }
      else
      {
        a_fluxMatrix(0, jmom) = 0;
      }
      ::pout() << ",icol = "  << jmom <<  ", Q(0,icol) = " << a_fluxMatrix(0, jmom) << endl;
    }
  }

  template<CENTERING cent, int order>
  void
  BaseWLSFlux<cent, order>::
  getWeightsMatrixFace(Chombo4::LAPACKMatrix  & a_weights,
                       Vector<int>            & a_weightDisVol,
                       Vector<int>            & a_weightDisFac,
                       Vector<int>            & a_weightDisEB) const
  {
    int nrows = a_weightDisVol.size()+a_weightDisFac.size()+a_weightDisEB.size();
    a_weights.define(nrows, nrows);
    a_weights.setVal(1.0);
    a_weights.setToIdentity(); //this sets all the off diagonals to zero
    int powerE = -3;
    // same code as least squares solvers
    for(int ineigh = 0; ineigh < nrows; ineigh++)
    {
      Real blerg;
      if (ineigh<a_weightDisVol.size())
      {
        int dist = a_weightDisVol[ineigh];
        blerg = POW(dist,powerE);
      }
      else if (ineigh>=a_weightDisVol.size() && ineigh<a_weightDisVol.size()+a_weightDisFac.size())
      {
        int dist = a_weightDisFac[ineigh-a_weightDisVol.size()];
        blerg = POW(dist,powerE);
      }
      if (ineigh>=a_weightDisVol.size()+a_weightDisFac.size())
      {
        int dist = a_weightDisEB[ineigh-a_weightDisVol.size()-a_weightDisFac.size()];
        blerg = POW(dist,powerE);
      }
      a_weights(ineigh, ineigh) = blerg;
      ::pout() << ",icol,irow = "  << ineigh <<  ", W(irow,icol) = " << a_weights(ineigh, ineigh) << endl;
    }
  }

  template<CENTERING cent, int order>
  void
  BaseWLSFlux<cent, order>::
  getStencilWeights(Chombo4::LAPACKMatrix  & a_phiCoef,
                    Chombo4::LAPACKMatrix  & a_Qshift,
                    Chombo4::LAPACKMatrix  & a_Mvol,
                    Chombo4::LAPACKMatrix  & a_weights,
                    double                   a_dx) const
  {
    Chombo4::LAPACKMatrix MT, WT, MTWT;
    MT = a_Mvol;
    MT.transpose();
    WT = a_weights;
    WT.transpose();
    multiply(MTWT, MT, WT);

    Chombo4::LAPACKMatrix QT;
    QT = a_Qshift;
    QT.transpose();

    Chombo4::LAPACKMatrix X, sT;
    int info = solveLSTSVD(X, MTWT, QT,
                           10,
                           1.0e-13);

    if(info != 0)
    {
      Chombo4::MayDay::Warning("non zero info returned");
    }

    multiply(sT, WT, X);
    a_phiCoef = sT;
    a_phiCoef.transpose();
  }

  template<CENTERING cent, int order>
  vector<EBIndex<CELL> >
  BaseWLSFlux<cent, order>::
  getFlapVoFs(const EBGraph& a_graph,  EBIndex<CELL>& a_start) const
  {
    vector<EBIndex<CELL> > retval;

    for(int idir  = 0; idir <  DIM; idir++)
    {
      for(SideIterator sit; sit.ok(); ++sit)
      {
        vector<EBIndex<CELL> > neighbors = a_graph.getNeighborVoFs(a_start, idir, sit());
        retval.insert(retval.end(), neighbors.begin(), neighbors.end());
      }
    }
    return retval;
  } 

  template<CENTERING cent, int order>
  vector<EBIndex<CELL> >
  BaseWLSFlux<cent, order>::
  getVoFsNeighbors(const EBIndex<cent> & a_face,
                   unsigned int          a_radius,
		   Vector<int>         & a_weightDisVol,
                   const EBGraph       & a_graph) const
  {
    vector<EBIndex<CELL> > retval;
    if (cent == BOUNDARY) 
    {
      EBIndex<CELL> voflo = a_face.getVoF();
      retval.push_back(voflo);
      a_weightDisVol.push_back(1);
    }
    else
    {
      EBIndex<CELL> voflo = a_face.getVoF(Side::Lo);
      EBIndex<CELL> vofhi = a_face.getVoF(Side::Hi);
      if (!a_face.m_isBoundary)
      {
        retval.push_back(voflo);
        a_weightDisVol.push_back(1);      
        retval.push_back(vofhi);	   
        a_weightDisVol.push_back(1);      
      } 
      else 
      { 
        if (a_face.m_vofIDMe>=0) 
        {
          retval.push_back(vofhi);
          a_weightDisVol.push_back(1);
        }
        if (a_face.m_vofIDLo>=0) 
        {
          retval.push_back(voflo);
          a_weightDisVol.push_back(1);
        }
      }
    }
    int startV = 0;
    for (int irad = 1; irad<=a_radius; irad++){
      //vector< EBIndex<CELL> > volsinclusive=  a_graph.getAllVoFsWithinRadius(voflo.m_pt, radius);
      int neighSize = retval.size();
      for (int ineigh = startV; ineigh<neighSize; ineigh++){
        vector< EBIndex<CELL> > volsinclusive=  getFlapVoFs(a_graph, retval[ineigh]);
        for(int ivof = 0; ivof < volsinclusive.size(); ivof++)
        {
          bool incl = true;
          for(int ivof2 = 0; ivof2 < retval.size(); ivof2++)
          {
            if (volsinclusive[ivof]==retval[ivof2])
            {
              incl=false;
              break;
            }
          }
          if (incl) 
          {
	    retval.push_back(volsinclusive[ivof]);
	    a_weightDisVol.push_back(irad+1);
	  }
        }
      }
      startV = neighSize;
    }
    return retval;
  }

  template<CENTERING cent, int order>
  void
  BaseWLSFlux<cent, order>::
  getMomentsInfoFromNeighbors(vector<EBIndex<CELL> >  & a_vols,
		              const EBGraph           & a_graph,
                              const VoluData          & a_voludata,
                              const EBFaData          & a_ebfadata,
                              const XFacData          & a_xfacdata,
                              const YFacData          & a_yfacdata,
                              const ZFacData          & a_zfacdata,
			      const EBNormalData      & a_ebnormxdata,
                              const EBNormalData      & a_ebnormydata,
#if DIM==3
                              const EBNormalData      & a_ebnormzdata,
#endif
                              Vector<IndMomDIM>       & a_volmom,
                              Vector<IndMomDIM>       & a_ebfmom,
                              Vector<IndMomSDMinOne>  & a_facmom,
                              Vector<RealVect>        & a_volneighLoc,
                              Vector<RealVect>        & a_ebfneighLoc,
                              Vector<RealVect>        & a_facneighLoc,
                              Vector<int>             & a_weightDisVol,
                              Vector<int>             & a_weightDisEB,
                              Vector<int>             & a_weightDisFac,
                              Vector<int>             & a_facneighDir,
			      Vector<Side::LoHiSide>  & a_facneighSid,
			      Vector<Vector<IndMomDIM> >& a_ebfneighNorm,
			      Vector<double>         & a_volneighFrac,
                              Vector<double>         & a_ebfneighFrac,
                              Vector<double>         & a_facneighFrac,
			      double                    a_dx) const
  {
    Vector< EBIndex<XFACE> > Xfacneigh;
    Vector< EBIndex<YFACE> > Yfacneigh;
    IndMomDIM regVolMom;
    regVolMom.setToRegular(a_dx);
    IndMomSDMinOne regFaceMom;
    regFaceMom.setToRegular(a_dx);
    for(int ivof = 0; ivof < a_vols.size(); ivof++)
    {
      RealVect volLoc;
      EBIndex<CELL> volInd = a_vols[ivof];
      for (int idir = 0; idir < DIM; idir++)
      {
        volLoc[idir] = a_dx * (0.5 + volInd.m_pt[idir]);
      }
      cout<< " vol order "<<volInd.m_pt<<endl;
      a_volneighLoc.push_back(volLoc);

      if(a_graph.isIrregular(volInd.m_pt))
      {
	a_volmom.push_back(a_voludata(volInd, 0));
	RealVect centroid;
        double kappa;
        Normalizor<order>::normalizedVoluMoments(kappa, centroid, a_voludata, volInd, a_graph, a_dx);
	a_volneighFrac.push_back(kappa);
	const EBIndex<BOUNDARY> ebf = volInd.getCutFace();
	a_ebfmom.push_back(a_ebfadata(ebf, 0));
        a_ebfneighLoc.push_back(volLoc);
        a_weightDisEB.push_back(a_weightDisVol[ivof]);
	RealVect normal;
        double bndryArea;
        // Get normal and boundary area
        Normalizor<order>::getNormalAndBoundaryArea(bndryArea, normal,
            a_voludata, a_ebfadata, a_xfacdata, a_yfacdata, a_zfacdata, volInd, a_graph, a_dx);
	Vector<IndMomDIM> ebnorm (DIM);
        ebnorm[0] = a_ebnormxdata(volInd, 0);
        ebnorm[1] = a_ebnormydata(volInd, 0);
#if DIM==3
        ebnorm[2] = a_ebnormzdata(volInd, 0);
#endif
	a_ebfneighNorm.push_back(ebnorm);
	a_ebfneighFrac.push_back(bndryArea);
        cout<< " eb order "<<volInd.m_pt<<endl;
      }
      else
      {
	a_volmom.push_back(regVolMom);
	a_volneighFrac.push_back(1.0);
      }
      int side = 0;
      for(SideIterator sit; sit.ok(); ++sit)
      {
        vector<EBIndex<XFACE> > xfaces = a_graph.getXFaces(volInd, sit());
        vector<EBIndex<YFACE> > yfaces = a_graph.getYFaces(volInd, sit());
#if DIM==3
        vector<EBIndex<ZFACE> > zfaces = a_graph.getZFaces(volInd, sit());
#endif
        if(xfaces.size() > 0)
        {
          for(int iface = 0; iface < xfaces.size(); iface++)
          {
            const EBIndex<XFACE> x_face = xfaces[iface];
            // check you have not added face to the list of face neigbors
            bool found = false;
            for(int jface = 0; jface < Xfacneigh.size(); jface++)
            {
              if(Xfacneigh[jface] == x_face)
              {
                found = true;
              }
            }
            if(!found && x_face.m_isBoundary)
            {
              RealVect XfaceLoc;
              for (int idir = 0; idir < DIM; idir++)
              {
                XfaceLoc[idir] = a_dx * (0.5 + x_face.m_pt[idir]);
                if (idir == 0) XfaceLoc[idir] = a_dx * (x_face.m_pt[idir]);
              }
    	      if(a_graph.isIrregular(volInd.m_pt))
              {
                a_facmom.push_back(a_xfacdata(x_face, 0));
		double   areaFrac;
                RealVect centroid;
                Normalizor<order>::normalizedFaceMoments(areaFrac, centroid, a_xfacdata, x_face, a_graph, a_dx);
	        a_facneighFrac.push_back(areaFrac);
    	      }
    	      else
              {
                a_facmom.push_back(regFaceMom);
	        a_facneighFrac.push_back(1.0);
    	      }
              Xfacneigh.push_back(x_face);
              a_facneighLoc.push_back(XfaceLoc);
              a_facneighDir.push_back(x_face.direction());
	      a_facneighSid.push_back(sit());
              a_weightDisFac.push_back(a_weightDisVol[ivof]);
            }
          }
        }
        if(yfaces.size() > 0)
        {
          for(int iface = 0; iface < yfaces.size(); iface++)
          {
            const EBIndex<YFACE> y_face = yfaces[iface];
            // check you have not added face to the list of face neigbors
            bool found = false;
            for(int jface = 0; jface < Yfacneigh.size(); jface++)
            { 
              if(Yfacneigh[jface] == y_face)
              {
                found = true;
              }
            } 
            if(!found && y_face.m_isBoundary)
            {
              RealVect YfaceLoc;
              for (int idir = 0; idir < DIM; idir++)
              {
                YfaceLoc[idir] = a_dx * (0.5 + y_face.m_pt[idir]);
                if (idir == 1) YfaceLoc[idir] = a_dx * (y_face.m_pt[idir]);
              }     
              if(a_graph.isIrregular(volInd.m_pt))
              {
                a_facmom.push_back(a_yfacdata(y_face, 0));
		double   areaFrac;
                RealVect centroid;
                Normalizor<order>::normalizedFaceMoments(areaFrac, centroid, a_yfacdata, y_face, a_graph, a_dx);
                a_facneighFrac.push_back(areaFrac);
              }
              else
              {
                a_facmom.push_back(regFaceMom);
                a_facneighFrac.push_back(1.0);
              }
              Yfacneigh.push_back(y_face);
              a_facneighLoc.push_back(YfaceLoc);
              a_facneighDir.push_back(y_face.direction());
	      a_facneighSid.push_back(sit());
              a_weightDisFac.push_back(a_weightDisVol[ivof]);
            }
          }
        }
#if DIM==3
	if(zfaces.size() > 0)
        {
          for(int iface = 0; iface < zfaces.size(); iface++)
          {
            const EBIndex<ZFACE> z_face = zfaces[iface];
            // check you have not added face to the list of face neigbors
            bool found = false;
            for(int jface = 0; jface < Zfacneigh.size(); jface++)
            {
              if(Zfacneigh[jface] == z_face)
              {
                found = true;
              }
            }
            if(!found && z_face.m_isBoundary)
            {
              RealVect ZfaceLoc;
              for (int idir = 0; idir < DIM; idir++)
              {
                ZfaceLoc[idir] = a_dx * (0.5 + z_face.m_pt[idir]);
                if (idir == 2) ZfaceLoc[idir] = a_dx * (z_face.m_pt[idir]);
              }
              if(a_graph.isIrregular(volInd.m_pt))
              {
                a_facmom.push_back(a_zfacdata(z_face, 0));
		double   areaFrac;
                RealVect centroid;
                Normalizor<order>::normalizedFaceMoments(areaFrac, centroid, a_zfacdata, z_face, a_graph, a_dx);
                a_facneighFrac.push_back(areaFrac);
              }
              else
              {
                a_facmom.push_back(regFaceMom);
                a_facneighFrac.push_back(1.0);
              }
              Zfacneigh.push_back(z_face);
              a_facneighLoc.push_back(ZfaceLoc);
              a_facneighDir.push_back(z_face.direction());
	      a_facneighSid.push_back(sit());
              a_weightDisFac.push_back(a_weightDisVol[ivof]);
            }
          }
        }
#endif
      }//End SideIterator loop
    }//End vol neighbors loop
  }

  template<CENTERING cent, int order>
  void
  BaseWLSFlux<cent, order>::
  getShiftedMomentMatrix(Chombo4::LAPACKMatrix    & a_Mvol,
                           RealVect               & a_faceLoc,
                           double                   a_dx,
                           Vector<IndMomDIM>      & a_volmom,
                           Vector<IndMomDIM>      & a_ebfmom,
                           Vector<IndMomSDMinOne> & a_facmom,
                           Vector<RealVect>       & a_volneighLoc,
                           Vector<RealVect>       & a_ebfneighLoc,
                           Vector<RealVect>       & a_facneighLoc,
			   Vector<int>            & a_facneighDir,
			   Vector<Side::LoHiSide> & a_facneighSid,
			   Vector<Vector<IndMomDIM> > & a_ebfneighNorm,
			   Vector<double>         & a_volneighFrac,
                           Vector<double>         & a_ebfneighFrac,
                           Vector<double>         & a_facneighFrac,
			   string                   a_ebbcName,
                           string                   a_dombcName[2*DIM]) const
  {
    Vector<IndMomDIM> fac2volmom;
    int nvol = a_volmom.size();
    int nfac = a_facmom.size();
    int nebf = a_ebfmom.size();
    fac2volmom.resize(nfac);

    int nrows = nvol + nfac + nebf;
    cout<< "nvol + nfac + nebf "<<nvol<<" "<<nfac<<" "<<nebf<<std::endl;
    MomentIterator<DIM, order> momit;
    int M = momit.size();
    int N = nrows;
    a_Mvol.define(N, M);
    cout<< "MVOL "<<N<<" "<<M<<std::endl;
  
    Real  volScaleFactor  = (D_TERM(a_dx, *a_dx, *a_dx));
    Real areaScaleFactor  = (D_TERM(1.0,  *a_dx, *a_dx));

    //Adding volume moments
    for(int irow = 0; irow < nvol; irow++)
    {
      //the moments start at the cell CENTER
      RealVect diff = a_volneighLoc[irow];
      diff -= a_faceLoc;
      IndexTM<Real, DIM> diffrv;
      for(int idir = 0; idir < DIM; idir++)
        {
          diffrv[idir] = diff[idir];
        }
      a_volmom[irow].shift(diffrv);

      if(Chombo4::LAPACKMatrix::s_verbose)
      {
        ::pout() << "irow = " << irow << endl;
      }
      for(momit.reset(); momit.ok(); ++momit)
      {
        //need to use the lower order moment declaration to
        //get the correct column number
        int icol = IndexedMoments<DIM, order>::indexOf(momit());
        //a_Mvol(irow, icol) = a_volmom[irow][momit()]/(a_volneighFrac[irow]*volScaleFactor*POW(a_dx, momit().sum()));
        a_Mvol(irow, icol) = a_volmom[irow][momit()]/(a_volneighFrac[irow]*volScaleFactor);
        if(Chombo4::LAPACKMatrix::s_verbose)
        {
          ::pout() << ",icol = "  << icol << ", momind = " << momit() << ", mat(irow,icol) = " << a_Mvol(irow, icol) << endl;
        }
      }
    }

    //Adding face moments
    for(int iface = 0; iface < nfac; iface++)
    {
      int irow = nvol + iface;
      RealVect diff = a_facneighLoc[iface];
      diff -= a_faceLoc;
      IndexTM<Real, DIM> diffrv;
      for(int idir = 0; idir < DIM; idir++)
      {
        diffrv[idir] = diff[idir];
      }

      int a_faceDir = a_facneighDir[iface];
      Side::LoHiSide a_faceSid = a_facneighSid[iface];
      int index = ebp_index(a_faceDir, a_faceSid);
      if(a_dombcName[index] == string("Neumann"))
      {
        NeumannBC<order>::
          getFaceMoments(fac2volmom[iface],
                         a_faceDir,
                         a_facmom[iface]);
      }
      else if(a_dombcName[index] == string("Dirichlet"))
      {
        DirichletBC<order>::
          getFaceMoments(fac2volmom[iface],
                         a_faceDir,
                         a_facmom[iface]);
      }
      else
      {
        PROTO_ASSERT(false, "could not find domBC string");
      }

      fac2volmom[iface].shift(diffrv);

      if(Chombo4::LAPACKMatrix::s_verbose)
      {
        ::pout() << "irow = " << irow << endl;
      }
      for(momit.reset(); momit.ok(); ++momit)
      {
        //need to use the lower order moment declaration to
        //get the correct column number
        int icol = IndexedMoments<DIM, order>::indexOf(momit());
        a_Mvol(irow, icol) = fac2volmom[iface][momit()]/(a_facneighFrac[iface]*areaScaleFactor*POW(a_dx, momit().sum()));
        //a_Mvol(irow, icol) = fac2volmom[iface][momit()]/(a_facneighFrac[iface]*areaScaleFactor);
        if(Chombo4::LAPACKMatrix::s_verbose)
        {
          ::pout() << ",icol = "  << icol << ", momind = " << momit() << ", mat(irow,icol) = " << a_Mvol(irow, icol) << endl;
        }
      }
    }

    //Adding eb moments
    for(int iebf = 0; iebf < nebf; iebf++)
    {
      int irow = nvol + nfac + iebf;
      RealVect diff = a_ebfneighLoc[iebf];
      diff -= a_faceLoc;
      IndexTM<Real, DIM> diffrv;
      for(int idir = 0; idir < DIM; idir++)
      {
        diffrv[idir] = diff[idir];
      }

      if(a_ebbcName == string("Neumann"))
      {
        NeumannBC<order>::
          getEBMoments(a_ebfmom[iebf],
		       a_ebfneighNorm[iebf]);
      }
      else if(a_ebbcName == string("Dirichlet"))
      {
      }
      else
      {
        PROTO_ASSERT(false, "could not find ebBC string");
      }

      a_ebfmom[iebf].shift(diffrv);

      if(Chombo4::LAPACKMatrix::s_verbose)
      {
        ::pout() << "irow = " << irow << endl;
      }
      for(momit.reset(); momit.ok(); ++momit)
      {
        //need to use the lower order moment declaration to
        //get the correct column number
        int icol = IndexedMoments<DIM, order>::indexOf(momit());
        //a_Mvol(irow, icol) = a_ebfmom[iebf][momit()]/(a_ebfneighFrac[iebf]*areaScaleFactor*POW(a_dx, momit().sum()));
        a_Mvol(irow, icol) = a_ebfmom[iebf][momit()]/(a_ebfneighFrac[iebf]*areaScaleFactor);
        if(Chombo4::LAPACKMatrix::s_verbose)
        {
          ::pout() << ",icol = "  << icol << ", momind = " << momit() << ", mat(irow,icol) = " << a_Mvol(irow, icol) << endl;
        }
      }
    }

 }

/*******************************/
  template<CENTERING cent, int order>
  void
  BaseWLSFlux<cent, order>::
  getFaceCenteredStencil(LocalStencil<CELL, double> & a_sten,
                         const EBIndex<cent>        & a_face,
                         const EBGraph              & a_graph,
                         const VoluData             & a_voludata, 
                         const EBFaData             & a_ebfadata,
                         const XFacData             & a_xfacdata,
                         const YFacData             & a_yfacdata,
                         const ZFacData             & a_zfacdata,
			 const EBNormalData         & a_ebnormxdata,
                         const EBNormalData         & a_ebnormydata,
#if DIM==3
                         const EBNormalData         & a_ebnormzdata,
#endif
			 string                       a_ebbcName,
                         string                       a_dombcName[2*DIM],
                         const int                  & a_dir,
			 IndMomSDMinOne             & a_faceMom,
                         double                       a_dx) const
  { 
    Chombo4::LAPACKMatrix::s_verbose = true; 

    RealVect faceLoc;
    for (int idir = 0; idir < DIM; idir++)
    {
      faceLoc[idir] = a_dx * (0.5 + a_face.m_pt[idir]);
      if (a_dir == idir) faceLoc[idir] = a_dx * (a_face.m_pt[idir]);
    }

    unsigned int radius = 2;
    Vector<int> weightDisVol;
    Vector<int> weightDisFac;
    Vector<int> weightDisEB;
    vector<EBIndex<CELL> >  vols = getVoFsNeighbors(a_face, radius, weightDisVol, a_graph);

    Vector<IndMomDIM> volmom, ebfmom;
    Vector<IndMomSDMinOne> facmom;
    Vector<RealVect> volneighLoc;
    Vector<RealVect> ebfneighLoc;
    Vector<RealVect> facneighLoc;
    Vector<double> volneighFrac;
    Vector<double> ebfneighFrac;
    Vector<double> facneighFrac;
    Vector<Vector<IndMomDIM> > ebfneighNorm;
    Vector<int> facneighDir;
    Vector<Side::LoHiSide> facneighSid;

    getMomentsInfoFromNeighbors(vols,
		                a_graph, a_voludata, a_ebfadata,
				a_xfacdata, a_yfacdata, a_zfacdata,
                                a_ebnormxdata, a_ebnormydata,
#if DIM==3
                                a_ebnormzdata,
#endif				
		                volmom, ebfmom, facmom,
				volneighLoc, ebfneighLoc, facneighLoc,
				weightDisVol, weightDisEB, weightDisFac,
				facneighDir, facneighSid, ebfneighNorm, 
				volneighFrac, ebfneighFrac, facneighFrac, a_dx);

    //first fill M with moments  shifted to the center of the face
    Chombo4::LAPACKMatrix  Mvol;
    getShiftedMomentMatrix(Mvol, faceLoc, a_dx, 
		           volmom, ebfmom, facmom,
		           volneighLoc, ebfneighLoc, facneighLoc,
			   facneighDir, facneighSid, ebfneighNorm,
			   volneighFrac, ebfneighFrac, facneighFrac, 
			   a_ebbcName, a_dombcName);

    // RHS of stencil equation
    Chombo4::LAPACKMatrix Qshift;
    getFaceFluxMatrix(Qshift, a_face, faceLoc, a_faceMom, a_dx);

    // Weights matrix
    Chombo4::LAPACKMatrix weights;
    getWeightsMatrixFace(weights, weightDisVol, weightDisFac, weightDisEB);

    // Compute phiCoef
    Chombo4::LAPACKMatrix phiCoef;
    getStencilWeights(phiCoef, Qshift, Mvol, weights, a_dx);

    //Filling stencil
    a_sten.clear();
    for(int ivof = 0; ivof < vols.size(); ivof++)
    {
      Real thisCo = phiCoef(0, ivof);
      int ivar = 0; //scalar
      a_sten.add(vols[ivof], thisCo);
      cout<<"coef "<<vols[ivof].m_pt<<"  "<<thisCo<<endl;
    }

  }

  template < int order>
  void
  WLSStencil<order>::
  getWLSDivFStencil(LocalStencil< CELL,  double>          & a_vofStencil,
                         const EBIndex<CELL>                   & a_vof,
                         const EBGraph                         & a_graph,
                         const VoluData                        & a_voludata,
                         const EBFaData                        & a_ebfadata,
                         const XFacData                        & a_xfacdata,
                         const YFacData                        & a_yfacdata,
                         const ZFacData                        & a_zfacdata,
			 const EBNormalData                    & a_ebnormxdata,
                         const EBNormalData                    & a_ebnormydata,
#if DIM==3
                         const EBNormalData                    & a_ebnormzdata,
#endif
                         const BaseWLSFlux<XFACE, order>  & a_xfaceflux,
                         const BaseWLSFlux<YFACE, order>  & a_yfaceflux,
                         const BaseWLSFlux<ZFACE, order>  & a_zfaceflux,
                         const BaseWLSFlux<BOUNDARY, order>  & a_ebflux,
                         string                             a_ebbcName,
                         string                             a_dombcName[2*DIM],
                         double                                  a_dx)
  {
    PR_TIME("getWLSDivFStencil");
    a_vofStencil.clear();

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

        IndMomSDMinOne faceMom;
	if (areaFrac < 1.0) 
	{
          faceMom = a_xfacdata(xfaces[iface], 0); 
	}
	else
	{
          faceMom.setToRegular(a_dx);
	}
        cout<<" here "<<xfaces[iface].m_pt<<endl;
        a_xfaceflux.getFaceCenteredStencil(fluxStencil, xfaces[iface],
                                       a_graph,
				       a_voludata,
				       a_ebfadata,    
                                       a_xfacdata,
                                       a_yfacdata,
                                       a_zfacdata,
				       a_ebnormxdata,
				       a_ebnormydata,
#if DIM==3
				       a_ebnormzdata,
#endif				  
				       a_ebbcName,
				       a_dombcName,
				       0,
				       faceMom, 
                                       a_dx);       
        cout<<" stencilX face "<<xfaces[iface].m_pt<<endl;
        cout<<" fluxStencil "<<endl;
        fluxStencil.print();
        fluxStencil *= double(isign);
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
        IndMomSDMinOne faceMom;
	if (areaFrac < 1.0) 
	{
          faceMom = a_yfacdata(yfaces[iface], 0); 
	}
	else
	{
          faceMom.setToRegular(a_dx);
	}

	a_yfaceflux.getFaceCenteredStencil(fluxStencil, yfaces[iface],
                                       a_graph,
                                       a_voludata,
                                       a_ebfadata,
                                       a_xfacdata,
                                       a_yfacdata,
                                       a_zfacdata,
				       a_ebnormxdata,
                                       a_ebnormydata,
#if DIM==3
                                       a_ebnormzdata,
#endif
				       a_ebbcName,
                                       a_dombcName,
                                       1,
				       faceMom, 
                                       a_dx);
	cout<<" stencilY face "<<yfaces[iface].m_pt<<endl;
        cout<<" fluxStencil "<<endl;
        fluxStencil.print();
        fluxStencil *= double(isign);
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
        IndMomSDMinOne faceMom;
	if (areaFrac < 1.0) 
	{
          faceMom = a_zfacdata(zfaces[iface], 0); 
	}
	else
	{
          faceMom.setToRegular(a_dx);
	}

        a_zfaceflux.getFaceCenteredStencil(fluxStencil, zfaces[iface],
                                           a_graph,    
                                           a_voludata,
                                           a_ebfadata,
                                           a_xfacdata,
                                           a_yfacdata,
                                           a_zfacdata,
					   a_ebnormxdata,
                                           a_ebnormydata,
#if DIM==3
                                           a_ebnormzdata,
#endif
                                           a_ebbcName,
                                           a_dombcName,
                                           2,
				           faceMom,
                                           a_dx);       


        fluxStencil *= double(isign)*areaFrac;
        a_vofStencil += fluxStencil;
      }
    }
#endif      

//    if(ideb == 1)
//    {
      cout << "before ebflux, ";
      cout << "stencil for " << a_vof.m_pt << ":" << endl;
      a_vofStencil.print();

    LocalStencil<CELL, double> ebflux;
    double bndryArea, alpha;
    if(a_graph.isIrregular(a_vof.m_pt))
    {
      EBIndex<BOUNDARY> cutface = a_vof.getCutFace();
      double weight; //for inhomogeneous bcs
      RealVect centroid, normal;
      IndMomDIM ebmom = a_ebfadata(cutface, 0);
      vector<IndMomDIM> ebnorm (DIM);
      ebnorm[0] = a_ebnormxdata(a_vof, 0);
      ebnorm[1] = a_ebnormydata(a_vof, 0);
#if DIM==3
      ebnorm[2] = a_ebnormzdata(a_vof, 0);
#endif
      Normalizor<order>::normalizedBndryMoments(alpha, centroid, a_ebfadata,  cutface, a_graph, a_dx);
      Normalizor<order>::getNormalAndBoundaryArea(bndryArea, normal,
                                                  a_voludata,
                                                  a_ebfadata,
                                                  a_xfacdata,
                                                  a_yfacdata,
                                                  a_zfacdata,
                                                  a_vof, a_graph, a_dx);
      if(a_ebbcName == string("Neumann"))
      {
        //neumann has no stencil
      }
      else if(a_ebbcName == string("Dirichlet"))
      {
          a_ebflux.getEBFluxStencil(ebflux, weight,
                           cutface,
                           a_graph,
                           a_voludata,
                           a_ebfadata,
                           a_xfacdata,
                           a_yfacdata,
                           a_zfacdata,
			   a_ebnormxdata,
                           a_ebnormydata,
#if DIM==3
                           a_ebnormzdata,
#endif
			   ebnorm,
  			   a_ebbcName,
                           a_dombcName,
			   ebmom,
			   normal,
			   bndryArea,
                           a_dx);
      }
      else
      {
        PROTO_ASSERT(false, "could not find EBBC string");
      }
    } 
    //for some reason, EBAMRPoissonOp calculates beta weight without the EB flux
    double diagweight;
    bool foundVoF = a_vofStencil.getValueAtVoF(diagweight, a_vof);
    PR_assert(foundVoF); //if this fails, the stencil has no diagonal value
    a_vofStencil.setBetaWeight(diagweight);
 
    if(a_graph.isIrregular(a_vof.m_pt)) 
    {
      cout<<" stencil ebface "<<a_vof.m_pt<<endl;
      cout<<" ebfluxStencil "<<endl;
      ebflux.print();
      //ebflux *= bndryArea;
      a_vofStencil += ebflux;
    }
#if 0
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
#endif
  }

  template < int order>
  void
  DirichletBC<order>::
  getFaceMoments(IndMomDIM      & a_fac2volmom,
                 int              a_faceDir,
                 IndMomSDMinOne & a_facmom)
  {
    MomentIterator<DIM, order> momit;

    for(momit.reset(); momit.ok(); ++momit)
    {
      if(momit()[a_faceDir] == 0)
      {
        IndexTM<int,DIM-1> areaMomIndex;
        int index = 0;
        for(int idir = 0; idir < DIM; idir++)
        {
          if(idir  != a_faceDir)
          {
            areaMomIndex[index] =  momit()[idir];
            index++;
          }
        }
        IndexTM<int, DIM> p = momit();
        a_fac2volmom[momit()] = a_facmom[areaMomIndex];
      }
      else
      {
        a_fac2volmom[momit()] = 0.0;
      }
    }
  }

  template < int order>
  void
  NeumannBC<order>::
  getFaceMoments(IndMomDIM      & a_fac2volmom,
                 int              a_faceDir,
                 IndMomSDMinOne & a_facmom)
  {
    MomentIterator<DIM, order> momit;
    IndMomDIM volmom;

    for(momit.reset(); momit.ok(); ++momit)
    {
      if(momit()[a_faceDir] == 0)
      {
        IndexTM<int,DIM-1> areaMomIndex;
        int index = 0;
        for(int idir = 0; idir < DIM; idir++)
        {
          if(idir  != a_faceDir)
          {
            areaMomIndex[index] =  momit()[idir];
            index++;
          }
        }
        IndexTM<int, DIM> p = momit();
        volmom[momit()] = a_facmom[areaMomIndex];
      }
      else
      {
        volmom[momit()] = 0.0;
      }
    }

    // shift the volume moments and multiply by p[facedir]
    for(momit.reset(); momit.ok(); ++momit)
    {
      IndexTM<int, DIM> p = momit();
      if(p[a_faceDir] == 0)
      {
        a_fac2volmom[p] = 0.0;
      }
      else
      {
        // Multiplying by p[facedir] is technically redundant because the below
        // expression evaluates to 0 except when p[facedir] = 1, but we write
        // the code this way for readability.
        a_fac2volmom[p] = p[a_faceDir]*volmom[p-BASISV_TM<int,DIM>(a_faceDir)];
      }
    }
  }

  template < int order>
  void
  NeumannBC<order>::
  getEBMoments(IndMomDIM         & a_ebfmom,
	       Vector<IndMomDIM> & a_ebfNorm)
  {
    // important because we will sum over each of the dimensions
    MomentIterator<DIM, order> momit;
    a_ebfmom.setToZero();
    for(int idir = 0; idir < DIM; idir++)
    {
      // add in pd*int((x-x0)^p-ed*nd*dA)
      for(momit.reset(); momit.ok(); ++momit)
      {
        IndexTM<int, DIM> jmom = momit();
        if(jmom[idir] > 0)
        {
          a_ebfmom[jmom] += jmom[idir]*a_ebfNorm[idir][jmom-BASISV_TM<int,DIM>(idir)];
        }
      }
    }
  }

  /***/
  template<CENTERING cent, int order>
  void
  BaseWLSFlux<cent, order>::
  getEBFluxStencil(LocalStencil<CELL, double>    & a_stencil,
                   double                        & a_weight, //for inhomogeneous bcs
                   const EBIndex<cent>           & a_face,
                   const EBGraph                 & a_graph,
                   const VoluData                & a_voludata,
                   const EBFaData                & a_ebfadata,
                   const XFacData                & a_xfacdata,
                   const YFacData                & a_yfacdata,
                   const ZFacData                & a_zfacdata,
		   const EBNormalData            & a_ebnormxdata,
                   const EBNormalData            & a_ebnormydata,
#if DIM==3
                   const EBNormalData            & a_ebnormzdata,
#endif
		   vector<IndMomDIM>             & a_ebnorm,
		   string                          a_ebbcName,
                   string                          a_dombcName[2*DIM],
		   IndMomDIM                     & a_ebmom,
		   RealVect                      & a_normal,
		   double                        & a_bndryArea,
                   const double                  & a_dx) const
  {
    Chombo4::LAPACKMatrix::s_verbose = true;

    RealVect faceLoc;
    for (int idir = 0; idir < DIM; idir++)
    {
      faceLoc[idir] = a_dx * (0.5 + a_face.m_pt[idir]);
    }
    unsigned int radius = 2;
    Vector<int> weightDisVol;
    Vector<int> weightDisFac;
    Vector<int> weightDisEB;
    vector<EBIndex<CELL> >  vols = getVoFsNeighbors(a_face, radius, weightDisVol, a_graph);

    Vector<IndMomDIM> volmom, ebfmom;
    Vector<IndMomSDMinOne> facmom;
    Vector<RealVect> volneighLoc;
    Vector<RealVect> ebfneighLoc;
    Vector<RealVect> facneighLoc;
    Vector<double> volneighFrac;
    Vector<double> ebfneighFrac;
    Vector<double> facneighFrac;
    Vector<Vector<IndMomDIM> > ebfneighNorm;
    Vector<int> facneighDir;
    Vector<Side::LoHiSide> facneighSid;

    getMomentsInfoFromNeighbors(vols,
                                a_graph, a_voludata, a_ebfadata,
                                a_xfacdata, a_yfacdata, a_zfacdata,
                                a_ebnormxdata, a_ebnormydata,
#if DIM==3
                                a_ebnormzdata,
#endif				
                                volmom, ebfmom, facmom,
                                volneighLoc, ebfneighLoc, facneighLoc,
                                weightDisVol, weightDisEB, weightDisFac,
                                facneighDir, facneighSid, ebfneighNorm, 
				volneighFrac, ebfneighFrac, facneighFrac, a_dx);

    //first fill M with moments  shifted to the center of the face
    Chombo4::LAPACKMatrix  Mvol;
    getShiftedMomentMatrix(Mvol, faceLoc, a_dx, 
		           volmom, ebfmom, facmom,
		           volneighLoc, ebfneighLoc, facneighLoc,
			   facneighDir, facneighSid, ebfneighNorm,
			   volneighFrac, ebfneighFrac, facneighFrac, 
			   a_ebbcName, a_dombcName);

    // RHS of stencil equation
    Chombo4::LAPACKMatrix Qshift;
    getEBFluxMatrix(Qshift, a_face, faceLoc, a_ebnorm, a_bndryArea, a_dx);

    // Weights matrix
    Chombo4::LAPACKMatrix weights;
    getWeightsMatrixFace(weights, weightDisVol, weightDisFac, weightDisEB);

    // Compute phiCoef
    Chombo4::LAPACKMatrix phiCoef;
    getStencilWeights(phiCoef, Qshift, Mvol, weights, a_dx);

    //Filling stencil
    a_stencil.clear();
    for(int ivof = 0; ivof < vols.size(); ivof++)
    {
      Real thisCo = phiCoef(0, ivof);
      int ivar = 0; //scalar
      a_stencil.add(vols[ivof], thisCo);
      cout<<"coef "<<vols[ivof].m_pt<<"  "<<thisCo<<endl;
    }

  }

  template < int order>
  Real
  NeumannExactSoln<order>::
  getFaceIntegral(const IndMomSDMinOne  & a_faceMom,
		  const RealVect        & a_faceLoc,
                  int                     a_faceDir,
	          const Real              a_dx) const 
  {
    Real retval = 0.0;
    MomentIterator<DIM, order> momit;
    for(momit.reset(); momit.ok(); ++momit)
    {
      IndexTM<int,DIM-1>  pSDMinOne = momit();
      IndexTM<int,DIM>  pSpaceDim;
      // face moment is zero unless p in face direction is 1 (because we take a derivative in 
      // the face direction)
      int index = 0;
      for(int idir = 0; idir < DIM; idir++)
      {
        if(idir == a_faceDir)
        {
          pSpaceDim[idir] = 1;
        }
        else
        {
          pSpaceDim[idir] = pSDMinOne[index];
          index++;
        }
      }

      Real moment = a_faceMom[pSDMinOne];
 
      IntVect power;
      RealVect point;
      Vector<int> powerVec;
      Vector<Real> pointVec;

      ParmParse pp;
      pp.getarr("phi_power",powerVec,0,DIM);
      pp.getarr("phi_point",pointVec,0,DIM);
      for(int idir = 0; idir < DIM; idir++)
      {
        point[idir] = pointVec[idir];
        power[idir] = powerVec[idir];
      }

      MonomialPhi<order> mp(power,point);

      // pSpaceDim[a_face.direction()] is 1 so no need to multiply
      retval += moment*mp.getTaylorCoefficient(a_faceLoc,pSpaceDim);
    }
  
    return retval;
  }

  template < int order>
  Real
  DirichletExactSoln<order>::
  getFaceIntegral(const IndMomSDMinOne  & a_faceMom,
		  const RealVect        & a_faceLoc,
                  int                     a_faceDir,
	          const Real              a_dx) const 
  {
    Real retval = 0.0;
    MomentIterator<DIM, order> momit;
    for(momit.reset(); momit.ok(); ++momit)
    {
      IndexTM<int,DIM> p = momit();
      Real moment = a_faceMom[p];

      IntVect power;
      RealVect point;
      Vector<int> powerVec;
      Vector<Real> pointVec;

      ParmParse pp;
      pp.getarr("phi_power",powerVec,0,DIM);
      pp.getarr("phi_point",pointVec,0,DIM);
      for(int idir = 0; idir < DIM; idir++)
      {
        point[idir] = pointVec[idir];
        power[idir] = powerVec[idir];
      }

      MonomialPhi<order> mp(power,point);

      retval += moment*mp.getTaylorCoefficient(a_faceLoc,p);

    }

    return retval;
  }

  template < int order>
  Real
  DirichletExactSoln<order>::
  getEBIntegral(const IndMomDIM  & a_ebfMom,
		const RealVect   & a_ebfLoc,
	        const Real         a_dx) const 
  {

    Real retval = 0.0;
    MomentIterator<DIM, order> momit;
    for(momit.reset(); momit.ok(); ++momit)
    {
      IndexTM<int,DIM> p = momit();
      Real moment = a_ebfMom[p];

      IntVect power;
      RealVect point;
      Vector<int> powerVec;
      Vector<Real> pointVec;

      ParmParse pp;
      pp.getarr("phi_power",powerVec,0,DIM);
      pp.getarr("phi_point",pointVec,0,DIM);
      for(int idir = 0; idir < DIM; idir++)
      {
        point[idir] = pointVec[idir];
        power[idir] = powerVec[idir];
      }

      MonomialPhi<order> mp(power,point);

      retval += moment*mp.getTaylorCoefficient(a_ebfLoc,p);

    }

    return retval;
  }

  template < int order>
  Real
  NeumannExactSoln<order>::
  getEBIntegral(const Vector<IndMomDIM>  & a_ebfNormalMom,
		const RealVect           & a_ebfLoc,
	        const Real                 a_dx) const 
  {

    Real retval = 0.0;
    for(int ndir = 0; ndir < DIM; ndir++)
    {
      MomentIterator<DIM, order> momit;
      for(momit.reset(); momit.ok(); ++momit)
      {
        IndexTM<int,DIM> p = momit();
        Real pd = p[ndir];

        if(pd > 0)
	{
	  IndexTM<int,DIM> pShift = p - BASISV_TM<int, DIM>(ndir);
	  Real moment = a_ebfNormalMom[ndir][pShift];
  
          IntVect power;
          RealVect point;
          Vector<int> powerVec;
          Vector<Real> pointVec;
    
          ParmParse pp;
          pp.getarr("phi_power",powerVec,0,DIM);
          pp.getarr("phi_point",pointVec,0,DIM);
          for(int idir = 0; idir < DIM; idir++)
          {
            point[idir] = pointVec[idir];
            power[idir] = powerVec[idir];
          }
    
          MonomialPhi<order> mp(power,point);
    
          retval += pd*moment*mp.getTaylorCoefficient(a_ebfLoc,p);
	}
      }
    }

    return retval;

  }


}
#endif
