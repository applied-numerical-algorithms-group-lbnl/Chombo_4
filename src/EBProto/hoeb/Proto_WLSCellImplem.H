#ifndef _Proto_WLSCellImplem_H_
#define _Proto_WLSCellImplem_H_

namespace Proto
{

  //Get neighbors VoFs through the faces of a specific cell.
  template<CENTERING cent, int order>
  vector<EBIndex<CELL> >
  BaseWLSFlux<cent, order>::
  getFlapVoFs(const EBGraph& a_graph, const EBIndex<CELL>& a_start) const
  {
    PR_TIME("getFlapVoFs");
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

  //Get neighbors VoFs of a given face within a specific radius.
  //Radius starts at one for the cells sharing the face, and 
  //increases by one for the neighbor cells sharing a face with the original
  //cells of radius one.
  template<CENTERING cent, int order>
  Vector<EBIndex<CELL> >
  BaseWLSFlux<cent, order>::
  getVoFsNeighbors(const EBIndex<cent> & a_face,
                   const EBIndex<CELL> & a_vof,
                   const int             a_radius,
                   Neighbors<order>    & a_neigh,
                   const EBGraph       & a_graph,
                   int                   a_neighType) const
  {
    PR_TIME("getVoFsNeighbors");
    Vector<int>  weightDisVol;
    Vector<EBIndex<CELL> > retval;
    MomentIterator<DIM,order> momit;
    if (a_neighType == 1)
    {
      retval.push_back(a_vof);
      weightDisVol.push_back(1);
    }
    else
    {
      if (cent == BOUNDARY) 
      {
        EBIndex<CELL> voflo = a_face.getVoF();
        retval.push_back(voflo);
        weightDisVol.push_back(1);
      }
      else
      {
        EBIndex<CELL> voflo = a_face.getVoF(Side::Lo);
        EBIndex<CELL> vofhi = a_face.getVoF(Side::Hi);
        if (!a_face.m_isBoundary)
        {
          retval.push_back(voflo);
          weightDisVol.push_back(1);      
          retval.push_back(vofhi);         
          weightDisVol.push_back(1);      
        } 
        else 
        { 
          if (a_face.m_vofIDMe>=0) 
          {
            retval.push_back(vofhi);
            weightDisVol.push_back(1);
          }
          if (a_face.m_vofIDLo>=0) 
          {
            retval.push_back(voflo);
            weightDisVol.push_back(1);
          }
        }
      }
    }
    int startV = 0;
    int new_rad = a_radius;
    for (int irad = 1; irad<=new_rad; irad++)
    {
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
            weightDisVol.push_back(irad+1);
          }
        }
      }
      startV = neighSize;
      if (irad==new_rad && retval.size()<=1.0*momit.size()) new_rad++;
    }
    a_neigh.setWeightDisVol(weightDisVol);
    return retval;
  }

  //Gets all the moments and geometric info for the neighbors cells, boundary faces
  //and EB faces. Stores all this info in the Neighbors object.
  template<CENTERING cent, int order>
  void
  BaseWLSFlux<cent, order>::
  getMomentsInfoFromNeighbors(Vector<EBIndex<CELL> >  & a_vols,
                              const EBGraph           & a_graph,
                              const VoluData          & a_voludata,
                              const EBFaData          & a_ebfadata,
                              const XFacData          & a_xfacdata,
                              const YFacData          & a_yfacdata,
                              const ZFacData          & a_zfacdata,
                              const EBNormalData      & a_ebnormxdata,
                              const EBNormalData      & a_ebnormydata,
                              const EBNormalData      & a_ebnormzdata,
                              Neighbors<order>  & a_neigh,
                              int                       a_neighType,
                              string                    a_ebbcName,
                              string                    a_dombcName[2*DIM],
                              double                    a_dx) const
  {
    PR_TIME("getMomentsInfoFromNeighbors");
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
    Vector<int> weightDisFac;
    Vector<int> weightDisEB;

    Vector< EBIndex<XFACE> > Xfacneigh;
    Vector< EBIndex<YFACE> > Yfacneigh;
    Vector< EBIndex<ZFACE> > Zfacneigh;
    IndMomDIM regVolMom;
    regVolMom.setToRegular(a_dx);
    IndMomSDMinOne regFaceMom;
    regFaceMom.setToRegular(a_dx);
    Vector<int> weightDisVol = a_neigh.m_weightDisVol;

    Vector<Point> volp; 
    Vector<Point> facp; 
    Vector<Point> ebfp; 

    for(int ivof = 0; ivof < a_vols.size(); ivof++)
    {
      RealVect volLoc;
      EBIndex<CELL> volInd = a_vols[ivof];
      for (int idir = 0; idir < DIM; idir++)
      {
        volLoc[idir] = a_dx * (0.5 + volInd.m_pt[idir]);
      }
      volp.push_back(volInd.m_pt);
      volneighLoc.push_back(volLoc);

      if(a_graph.isIrregular(volInd.m_pt))
      {
        volmom.push_back(a_voludata(volInd, 0));
        volneighFrac.push_back(a_voludata(volInd, 0)[IndexTM<int,DIM>::Zero]/POW(a_dx,DIM));
        if (this->m_fluxType==0){
          const EBIndex<BOUNDARY> ebf = volInd.getCutFace();
          Vector<IndMomDIM> ebnorm (DIM);
          ebnorm[0] = a_ebnormxdata(volInd, 0);
          ebnorm[1] = a_ebnormydata(volInd, 0);
#if DIM==3
          ebnorm[2] = a_ebnormzdata(volInd, 0);
#endif
          ebfmom.push_back(a_ebfadata(ebf, 0));
          ebfneighLoc.push_back(volLoc);
          weightDisEB.push_back(weightDisVol[ivof]);
          ebfneighNorm.push_back(ebnorm);
          ebfneighFrac.push_back(a_ebfadata(ebf, 0)[IndexTM<int,DIM>::Zero]/POW(a_dx,DIM-1));
          ebfp.push_back(volInd.m_pt);
        }
      }
      else
      {
        volmom.push_back(regVolMom);
        volneighFrac.push_back(1.0);
      }
      if (this->m_fluxType==0){
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
              int index = ebp_index(0, sit());
              if(!found && x_face.m_isBoundary
                 && (a_dombcName[index] != string("Extrapolated")))
              {
                RealVect XfaceLoc;
                for (int idir = 0; idir < DIM; idir++)
                {
                  XfaceLoc[idir] = a_dx * (0.5 + x_face.m_pt[idir]);
                  if (idir == 0) XfaceLoc[idir] = a_dx * (x_face.m_pt[idir]);
                }
                if(a_graph.isIrregular(volInd.m_pt))
                {
                  facmom.push_back(a_xfacdata(x_face, 0));
                  facneighFrac.push_back(a_xfacdata(x_face, 0)[IndexTM<int,DIM-1>::Zero]/POW(a_dx,DIM-1));
                }
                else
                {
                  facmom.push_back(regFaceMom);
                  facneighFrac.push_back(1.0);
                }
                Xfacneigh.push_back(x_face);
                facneighLoc.push_back(XfaceLoc);
                facneighDir.push_back(x_face.direction());
                facneighSid.push_back(sit());
                weightDisFac.push_back(weightDisVol[ivof]);
                facp.push_back(x_face.m_pt);
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
                  facmom.push_back(a_yfacdata(y_face, 0));
                  facneighFrac.push_back(a_yfacdata(y_face, 0)[IndexTM<int,DIM-1>::Zero]/POW(a_dx,DIM-1));
                }
                else
                {
                  facmom.push_back(regFaceMom);
                  facneighFrac.push_back(1.0);
                }
                Yfacneigh.push_back(y_face);
                facneighLoc.push_back(YfaceLoc);
                facneighDir.push_back(y_face.direction());
                facneighSid.push_back(sit());
                weightDisFac.push_back(weightDisVol[ivof]);
                facp.push_back(y_face.m_pt);
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
                  facmom.push_back(a_zfacdata(z_face, 0));
                  facneighFrac.push_back(a_zfacdata(z_face, 0)[IndexTM<int,DIM-1>::Zero]/POW(a_dx,DIM-1));
                }
                else
                {
                  facmom.push_back(regFaceMom);
                  facneighFrac.push_back(1.0);
                }
                Zfacneigh.push_back(z_face);
                facneighLoc.push_back(ZfaceLoc);
                facneighDir.push_back(z_face.direction());
                facneighSid.push_back(sit());
                weightDisFac.push_back(weightDisVol[ivof]);
              }
            }
          }
#endif
        }//End SideIterator loop
      }
    }//End vol neighbors loop

    if(a_ebbcName == string("no_bcs"))
      {
        ebfmom.resize(0);
        ebfneighNorm.resize(0);
        ebfneighLoc.resize(0);
        ebfneighFrac.resize(0);
        weightDisEB.resize(0);
      }
    a_neigh.setMomentsInfo(volmom, facmom,
                           ebfmom, ebfneighNorm,
                           a_dx);
    
    a_neigh.setLocInfo(volneighLoc, facneighLoc,
                       ebfneighLoc);
      
    a_neigh.setWeightDisFacAndEB(weightDisFac,
                                 weightDisEB);

    a_neigh.setFracInfo(volneighFrac,  
                        ebfneighFrac,
                        facneighFrac);

    a_neigh.setAdditionalFaceInfo(facneighDir,
                                  facneighSid);

    if (m_verbose)
    {
      pout()<<" Index order for the matrices " << endl;
      for (int i=0; i<volp.size(); i++)
      {
        pout()<<" Vol "<<volp[i]<<endl;
      }
      for (int i=0; i<facp.size(); i++)
      {
        pout()<<" Face "<<facp[i]<<endl;
      }
      for (int i=0; i<ebfp.size(); i++)
      {
        pout()<<" EBFace "<<ebfp[i]<<endl;
      }
    }
  }

  //Compute the Moments Matrix
  template<CENTERING cent, int order>
  void
  BaseWLSFlux<cent, order>::
  getShiftedMomentMatrix(Chombo4::LAPACKMatrix    & a_Mvol,
                           RealVect               & a_faceLoc,
                           double                   a_dx,
                           const Neighbors<order> & a_neigh,
                           string                   a_ebbcName,
                           string                   a_dombcName[2*DIM]) const
  {
    PR_TIME("getShiftedMomentMatrix");
    Vector<IndMomDIM> volmom = a_neigh.m_volMom;
    Vector<IndMomDIM> ebfmom = a_neigh.m_ebfMom;
    Vector<IndMomSDMinOne> facmom = a_neigh.m_facMom;
    Vector<RealVect> volneighLoc = a_neigh.m_volLoc;
    Vector<RealVect> ebfneighLoc = a_neigh.m_ebfLoc;
    Vector<RealVect> facneighLoc = a_neigh.m_facLoc;
    Vector<int> facneighDir = a_neigh.m_facneighDir;
    Vector<Side::LoHiSide> facneighSid = a_neigh.m_facneighSid;
    Vector<Vector<IndMomDIM> > ebfneighNorm = a_neigh.getEBNormMoments();
    Vector<double> volneighFrac = a_neigh.m_volneighFrac;
    Vector<double> ebfneighFrac = a_neigh.m_ebfneighFrac;
    Vector<double> facneighFrac = a_neigh.m_facneighFrac;
    Vector<IndMomDIM> fac2volmom;

    int nvol = volmom.size();
    int nfac = facmom.size();
    int nebf = ebfmom.size();
    fac2volmom.resize(nfac);

    int nrows = nvol + nfac + nebf;
    MomentIterator<DIM, order> momit;
    int M = momit.size();
    int N = nrows;
    a_Mvol.define(N, M);
    if (m_verbose) 
    {
      pout()<< "Size Moments matrix "<<N<<" "<<M<<std::endl;
      pout()<< "With nvol + nfac + nebf "<<nvol<<" "<<nfac<<" "<<nebf<<std::endl;
      pout()<< "Moments matrix" <<endl;
    }
  
    Real  volScaleFactor  = (D_TERM(a_dx, *a_dx, *a_dx));
    Real areaScaleFactor  = (D_TERM(1.0,  *a_dx, *a_dx));
    double eps = 1e-9;

    //Adding volume moments
    for(int irow = 0; irow < nvol; irow++)
    {
      //the moments start at the cell CENTER
      RealVect diff = volneighLoc[irow];
      diff -= a_faceLoc;
      IndexTM<Real, DIM> diffrv;
      for(int idir = 0; idir < DIM; idir++)
        {
          diffrv[idir] = diff[idir]/a_dx;
        }
      volmom[irow].shift(diffrv);

      if(Chombo4::LAPACKMatrix::s_verbose)
      {
        ::pout() << "irow = " << irow << endl;
      }
      for(momit.reset(); momit.ok(); ++momit)
      {
        //need to use the lower order moment declaration to
        //get the correct column number
        int icol = IndexedMoments<DIM, order>::indexOf(momit());
        if (volneighFrac[irow] < eps)
        {
          if (icol == 0)
          {
            a_Mvol(irow, icol) = 1.0;
          }
          else
          {
            a_Mvol(irow, icol) = 0.0;
          }
        }
        else
        {
          a_Mvol(irow, icol) = volmom[irow][momit()]/(volneighFrac[irow]);
        }
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
      RealVect diff = facneighLoc[iface];
      diff -= a_faceLoc;
      IndexTM<Real, DIM> diffrv;
      for(int idir = 0; idir < DIM; idir++)
      {
        diffrv[idir] = diff[idir]/a_dx;
      }

      int a_faceDir = facneighDir[iface];
      Side::LoHiSide faceSid = facneighSid[iface];
      int index = ebp_index(a_faceDir, faceSid);
      if(a_dombcName[index] == string("Neumann"))
      {
        NeumannBC<order>::
          getFaceMoments(fac2volmom[iface],
                         a_faceDir,
                         facmom[iface]);
      }
      else if(a_dombcName[index] == string("Dirichlet"))
      {
        DirichletBC<order>::
          getFaceMoments(fac2volmom[iface],
                         a_faceDir,
                         facmom[iface]);
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
        if (facneighFrac[iface] < eps)
        {
          if (icol == 0)
          {
            a_Mvol(irow, icol) = 1.0;
          }
          else
          {
            a_Mvol(irow, icol) = 0.0;
          }
        }
        else
        {
          a_Mvol(irow, icol) = fac2volmom[iface][momit()]/(facneighFrac[iface]);
        }
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
      RealVect diff = ebfneighLoc[iebf];
      diff -= a_faceLoc;
      IndexTM<Real, DIM> diffrv;
      for(int idir = 0; idir < DIM; idir++)
      {
        diffrv[idir] = diff[idir]/a_dx;
      }

      if(a_ebbcName == string("Neumann"))
      {
        NeumannBC<order>::
          getEBMoments(ebfmom[iebf],
                       ebfneighNorm[iebf]);
      }
      else if(a_ebbcName == string("Dirichlet"))
      {
      }
      else
      {
        PROTO_ASSERT(false, "could not find ebBC string");
      }

      ebfmom[iebf].shift(diffrv);

      if(Chombo4::LAPACKMatrix::s_verbose)
      {
        ::pout() << "irow = " << irow << endl;
      }
      for(momit.reset(); momit.ok(); ++momit)
      {
        //need to use the lower order moment declaration to
        //get the correct column number
        int icol = IndexedMoments<DIM, order>::indexOf(momit());
        if (ebfneighFrac[iebf] < eps)
        {
          if (icol == 0)
          {
            a_Mvol(irow, icol) = 1.0;
          }
          else
          {
            a_Mvol(irow, icol) = 0.0;
          }
        }
        else
        {
          a_Mvol(irow, icol) = ebfmom[iebf][momit()]/(ebfneighFrac[iebf]);
        }
        if(Chombo4::LAPACKMatrix::s_verbose)
        {
          ::pout() << ",icol = "  << icol << ", momind = " << momit() << ", mat(irow,icol) = " << a_Mvol(irow, icol) << endl;
        }
      }
    }
  }

  //Compute rhs matrix for a face
  template<CENTERING cent, int order>
  void
  BaseWLSFlux<cent, order>::
  getFaceFluxMatrix(Chombo4::LAPACKMatrix    & a_fluxMatrix,
                    const EBIndex<cent>      & a_face,
                    RealVect                 & a_faceLoc,
                    IndMomSDMinOne           & a_faceMom,
                    double                     a_signNormal,
                    double                     a_dx) const
  {
    PR_TIME("getFaceFluxMatrix");

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
        diffrv[index] = diff[idir]/a_dx;
        index++;
      }
    }
    a_faceMom.shift(diffrv);

    MomentIterator<DIM, order> momit;
    a_fluxMatrix.define(1, momit.size());
    a_fluxMatrix.setVal(0.);
    const int& faceDir = a_face.direction();
    if (m_verbose) pout() << " RHS - Face Flux Matrix " << endl;
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
        a_fluxMatrix(0, jmom) = pd*a_faceMom[pfaceshift];
      }
      else
      {
        a_fluxMatrix(0, jmom) = 0;
      }
      if(Chombo4::LAPACKMatrix::s_verbose)
      {
        ::pout() << ",icol = "  << jmom <<  ", Q(0,icol) = " << a_fluxMatrix(0, jmom) << endl;
      }
    }
  }

  //Compute rhs matrix for a EB face
  template<CENTERING cent, int order>
  void
  BaseWLSFlux<cent, order>::
  getEBFluxMatrix(Chombo4::LAPACKMatrix    & a_fluxMatrix,
                  const EBIndex<cent>      & a_face,
                  RealVect                 & a_faceLoc,
                  Vector<IndMomDIM>        & a_ebnorm,
                  double                     a_dx) const
  {
    PR_TIME("getEBFluxMatrix");
    MomentIterator<DIM, order> momit;
    a_fluxMatrix.define(1, momit.size());
    a_fluxMatrix.setVal(0.);

    RealVect startLoc;
    for (int idir = 0; idir < DIM; idir++)
    {
      startLoc[idir] = a_dx * (0.5 + a_face.m_pt[idir]);
    }

    for(int idir = 0; idir < DIM; idir++)
    {

      RealVect diff = startLoc;
      diff -= a_faceLoc;
      IndexTM<Real, DIM> diffrv;
      for(int idir = 0; idir < DIM; idir++)
      {
        diffrv[idir] = diff[idir]/a_dx;
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
          fluxIncr(0, jmom) = pd*a_ebnorm[idir][pshift];
        }
        else
        {
          fluxIncr(0, jmom) = 0;
        }
      }
      a_fluxMatrix += fluxIncr;
    }
    if (m_verbose)  pout() << " RHS - EB Flux Matrix " << endl;
    for(int idir = 0; idir < DIM; idir++)
    {
      for(momit.reset(); momit.ok(); ++momit)
      {
        IndexTM<int, DIM> p = momit();
        //Real pd = p[idir];
        //get column number associated with this particular moment
        int jmom = IndexedMoments<DIM, order>::indexOf(p);
        if(Chombo4::LAPACKMatrix::s_verbose)
        {
          ::pout() << ",icol = "  << jmom <<  ", Q(0,icol) = " << a_fluxMatrix(0, jmom) << endl;
        }
      }
    }
   
  }

  //Compute weight matrix
  template<CENTERING cent, int order>
  void
  BaseWLSFlux<cent, order>::
  getWeightsMatrixFace(Chombo4::LAPACKMatrix  & a_weights,
                       const Neighbors<order> & a_neigh,
                       RealVect               & a_faceLoc) const
  { 
    PR_TIME("getWeightsMatrixFace");
    Vector<int> weightDisVol = a_neigh.m_weightDisVol;
    Vector<int> weightDisFac = a_neigh.m_weightDisFac;
    Vector<int> weightDisEB = a_neigh.m_weightDisEB;

    int nrows = weightDisVol.size()+weightDisFac.size()+weightDisEB.size();
    a_weights.define(nrows, nrows);
    a_weights.setVal(1.0);
    a_weights.setToIdentity(); //this sets all the off diagonals to zero

    int max_rad = 0.0;
    for(int ineigh = 0; ineigh < weightDisVol.size(); ineigh++)
    {
      if (max_rad<weightDisVol[ineigh])
      {
        max_rad = weightDisVol[ineigh];
      }
    }

    //int powerE;
    //if (max_rad-1 == m_radius) powerE = m_powerE;
    //else powerE = m_powerE - (2.0*(max_rad-m_radius-1));
    int powerE = m_powerE;
    if (m_verbose)
      {
        pout() << " Weights Matrix " << endl;
      }
    for(int ineigh = 0; ineigh < nrows; ineigh++)
    {
      Real blerg;
      if (ineigh<weightDisVol.size())
      {
        int dist = weightDisVol[ineigh];
        blerg = POW(dist,powerE);
      }
      else if (ineigh>=weightDisVol.size() && ineigh<weightDisVol.size()+weightDisFac.size())
      {
        int dist = weightDisFac[ineigh-weightDisVol.size()];
        blerg = POW(dist,powerE);
      }
      if (ineigh>=weightDisVol.size()+weightDisFac.size())
      {
        int dist = weightDisEB[ineigh-weightDisVol.size()-weightDisFac.size()];
        blerg = POW(dist,powerE);
      }
      a_weights(ineigh, ineigh) = blerg;
      // sane bounds for weights - not a hard rule, just the convention we use
      CH_assert(a_weights(ineigh, ineigh) > 0.);
      CH_assert(a_weights(ineigh, ineigh) <= 1.);
      if(Chombo4::LAPACKMatrix::s_verbose)
      {
        pout() << ",icol,irow = " << ineigh <<  ", W(irow,icol) = " << a_weights(ineigh, ineigh) << endl;
      }
    }
  }

  //Compute stencil coefficients matrix
  template<CENTERING cent, int order>
  void
  BaseWLSFlux<cent, order>::
  getStencilWeights(Chombo4::LAPACKMatrix        & a_phiCoef,
                    const Chombo4::LAPACKMatrix  & a_Qshift,
                    const Chombo4::LAPACKMatrix  & a_Mvol,
                    const Chombo4::LAPACKMatrix  & a_weights,
                    double                         a_dx) const
  {
    PR_TIME("getStencilWeights");
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

  //Compute inhomogeneous term
  template<CENTERING cent, int order>
  Real
  BaseWLSFlux<cent, order>::
  computeInhomogTerm(const Chombo4::LAPACKMatrix    & a_stencil,
                     const Neighbors<order>         & a_neigh,
                     const int                        a_nvolneigh,
                     string                           a_ebbcName,
                     string                           a_dombcName[2*DIM],
                     double                           a_dx) const
  {
    PR_TIME("computeInhomogTerm");
    Vector<IndMomSDMinOne> a_facmom = a_neigh.m_facMom;
    Vector<IndMomDIM> a_ebfmom= a_neigh.m_ebfMom;
    Vector<Vector<IndMomDIM> > a_ebfneighNorm = a_neigh.m_ebfneighNorm;
    Vector<RealVect> a_facneighLoc = a_neigh.m_facLoc;
    Vector<int> a_facneighDir = a_neigh.m_facneighDir;
    Vector<Side::LoHiSide> a_facneighSid = a_neigh.m_facneighSid;
    Vector<RealVect> a_ebfneighLoc = a_neigh.m_ebfLoc;
    Vector<double> ebfneighFrac = a_neigh.m_ebfneighFrac;
    Vector<double> facneighFrac = a_neigh.m_facneighFrac;
    Real areaScaleFactor = (D_TERM(1.0,  *a_dx, *a_dx));

    Vector<Real> domBCValues, ebBCValues;
    int nrows = a_facmom.size();
    domBCValues.resize(nrows);
    for(int iface = 0; iface < nrows; iface++)
    {
      RealVect a_faceLocDiff = a_facneighLoc[iface];
      int a_faceDir = a_facneighDir[iface];
      Side::LoHiSide a_faceSid = a_facneighSid[iface];
      int index = ebp_index(a_faceDir, a_faceSid);
      if(a_dombcName[index] == string("Neumann"))
      {
        domBCValues[iface] = NeumannExactSoln<order>::
                               getFaceIntegral(a_facmom[iface],
                                               a_faceLocDiff,
                                               a_faceDir,
                                               m_exactSol,
                                               a_dx);
      }
      else if(a_dombcName[index] == string("Dirichlet"))
      {
        domBCValues[iface] = DirichletExactSoln<order>::
                               getFaceIntegral(a_facmom[iface],
                                               a_faceLocDiff,
                                               a_faceDir,
                                               m_exactSol,
                                               a_dx);
      }
      else
      {
        PROTO_ASSERT(false, "could not find domBC string");
      }

      if (facneighFrac[iface]<1.0e-11)
      {
        domBCValues[iface] = 0.0;
      }
      else
      {
        domBCValues[iface] = domBCValues[iface]/(facneighFrac[iface]*areaScaleFactor);
      }
    }

    nrows = a_ebfmom.size();
    ebBCValues.resize(nrows);
    for(int ivof = 0; ivof < nrows; ivof++)
    {
      RealVect a_ebfLocDiff = a_ebfneighLoc[ivof]; 
      if(a_ebbcName == string("Neumann"))
      {
        ebBCValues[ivof] = NeumannExactSoln<order>::
                               getEBIntegral(a_ebfneighNorm[ivof],
                                             a_ebfLocDiff,
                                             m_exactSol,
                                             a_dx);
      }
      else if(a_ebbcName == string("Dirichlet"))
      {
        ebBCValues[ivof] = DirichletExactSoln<order>::
                               getEBIntegral(a_ebfmom[ivof],
                                             a_ebfLocDiff,
                                             m_exactSol,
                                             a_dx);
      }
      else
      {
        PROTO_ASSERT(false, "could not find domBC string");
      }

      if (ebfneighFrac[ivof]<1.0e-11)
      {
        ebBCValues[ivof] = 0.0;
      }
      else
      {
        ebBCValues[ivof] = ebBCValues[ivof]/(ebfneighFrac[ivof]*areaScaleFactor);
      }
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
      if (m_verbose) pout()<<" ebbc val "<<ebBCValues[iface]<<"  coef "<<a_stencil(0, a_nvolneigh + nfac + iface)<<endl;
    }

    return inhomogTerm;
  }

  //Compute face centered flux stencil
  template<CENTERING cent, int order>
  void
  BaseWLSFlux<cent, order>::
  getFaceCenteredStencil(LocalStencil<CELL, double> & a_sten,
                         Real                       & a_inhomogTerm, //for inhomogeneous bcs
                         const EBIndex<CELL>        & a_vof,
                         const EBIndex<cent>        & a_face,
                         const EBGraph              & a_graph,
                         const VoluData             & a_voludata, 
                         const EBFaData             & a_ebfadata,
                         const XFacData             & a_xfacdata,
                         const YFacData             & a_yfacdata,
                         const ZFacData             & a_zfacdata,
                         const EBNormalData         & a_ebnormxdata,
                         const EBNormalData         & a_ebnormydata,
                         const EBNormalData         & a_ebnormzdata,
                         string                       a_ebbcName,
                         string                       a_dombcName[2*DIM],
                         const int                    a_dir,
                         IndMomSDMinOne             & a_faceMom,
                         int                        & a_neighType,
                         double                       a_signNormal, 
                         double                       a_dx) const
  { 
    PR_TIME("getFaceCenteredStencil");
    Chombo4::LAPACKMatrix::s_verbose = m_verbose; 

    Neighbors<order> neighs;

    RealVect faceLoc;
    for (int idir = 0; idir < DIM; idir++)
    {
      faceLoc[idir] = a_dx * (0.5 + a_face.m_pt[idir]);
      if (a_dir == idir) faceLoc[idir] = a_dx * (a_face.m_pt[idir]);
    }

    int radius = m_radius;
    Vector<EBIndex<CELL> >  vols = getVoFsNeighbors(a_face, a_vof, radius, neighs, a_graph, a_neighType);
    MomentIterator<DIM,order> momit;
    if (vols.size()<=momit.size()) pout()<<"DANGER: NUMBER OF NEIGHS "<<vols.size()<<" and NUMBER OF coef "<<momit.size()<<endl;
    getMomentsInfoFromNeighbors(vols,
                                a_graph, a_voludata, a_ebfadata,
                                a_xfacdata, a_yfacdata, a_zfacdata,
                                a_ebnormxdata, a_ebnormydata,
                                a_ebnormzdata,
                                neighs, a_neighType,
                                a_ebbcName, a_dombcName, a_dx);

    //first fill M with moments shifted to the center of the face
    Chombo4::LAPACKMatrix Mvol;
    getShiftedMomentMatrix(Mvol, faceLoc, a_dx,
                           neighs, 
                           a_ebbcName, a_dombcName);

    // RHS of stencil equation
    Chombo4::LAPACKMatrix Qshift;
    getFaceFluxMatrix(Qshift, a_face, faceLoc, a_faceMom, a_signNormal, a_dx);

    // Weights matrix
    Chombo4::LAPACKMatrix weights;
    RealVect centroid_pos;
    centroid_pos[0] = a_faceMom[1]/a_faceMom[0];
#if DIM==3
    centroid_pos[1] = a_faceMom[1+order]/a_faceMom[0];
#endif
    int centr_iter = 0;
    for (int idir = 0; idir < DIM; idir++)
    {
      if (a_dir != idir) 
      {
        faceLoc[idir] += centroid_pos[centr_iter]*a_dx;
        centr_iter++;
      }
    }
    if (m_verbose) pout()<<" centroid comp "<<faceLoc<<endl;
    getWeightsMatrixFace(weights, neighs, faceLoc);

    // Compute phiCoef
    Chombo4::LAPACKMatrix phiCoef;
    getStencilWeights(phiCoef, Qshift, Mvol, weights, a_dx);

    //Filling stencil
    a_sten.clear();
    Real coef_acc = 0.0;
    if (m_verbose)  pout() << " Coefs Matrix " << endl;
    for(int ivof = 0; ivof < vols.size(); ivof++)
    {
      Real thisCo = phiCoef(0, ivof);
      a_sten.add(vols[ivof], thisCo);
      if (m_verbose){
        pout()<<"stencil coefs for Vof: "<<vols[ivof].m_pt<<" Coef: "<<thisCo<<endl;
      }
      coef_acc += thisCo;
    }

    int nvolneigh = vols.size();
    a_inhomogTerm = computeInhomogTerm(phiCoef, neighs,
                                       nvolneigh, 
                                       a_ebbcName, a_dombcName,
                                       a_dx);
    if (m_verbose) pout()<<"For phi=1 "<<coef_acc<<" "<<a_inhomogTerm<<" "<<coef_acc+a_inhomogTerm<<endl;
  }

  //Compute EB face stencil
  template<CENTERING cent, int order>
  void
  BaseWLSFlux<cent, order>::
  getEBFluxStencil(LocalStencil<CELL, double>    & a_stencil,
                   Real                          & a_inhomogTerm, //for inhomogeneous bcs
                   const EBIndex<CELL>           & a_vof,
                   const EBIndex<cent>           & a_face,
                   const EBGraph                 & a_graph,
                   const VoluData                & a_voludata,
                   const EBFaData                & a_ebfadata,
                   const XFacData                & a_xfacdata,
                   const YFacData                & a_yfacdata,
                   const ZFacData                & a_zfacdata,
                   const EBNormalData            & a_ebnormxdata,
                   const EBNormalData            & a_ebnormydata,
                   const EBNormalData            & a_ebnormzdata,
                   Vector<IndMomDIM>             & a_ebnorm,
                   string                          a_ebbcName,
                   string                          a_dombcName[2*DIM],
                   int                             a_neighType,
                   const double                    a_dx) const
  {
    PR_TIME("getEBFluxStencil");
    Chombo4::LAPACKMatrix::s_verbose = m_verbose;

    Neighbors<order> neighs;

    RealVect faceLoc;
    for (int idir = 0; idir < DIM; idir++)
    {
      faceLoc[idir] = a_dx * (0.5 + a_face.m_pt[idir]);
    }
    int radius = m_radius;
    Vector<EBIndex<CELL> >  vols = getVoFsNeighbors(a_face, a_vof, radius, neighs, a_graph, 0);

    getMomentsInfoFromNeighbors(vols,
                                a_graph, a_voludata, a_ebfadata,
                                a_xfacdata, a_yfacdata, a_zfacdata,
                                a_ebnormxdata, a_ebnormydata,
                                a_ebnormzdata,
                                neighs, a_neighType,
                                a_ebbcName, a_dombcName, a_dx);

    //first fill M with moments  shifted to the center of the face
    Chombo4::LAPACKMatrix  Mvol;
    getShiftedMomentMatrix(Mvol, faceLoc, a_dx,
                           neighs, 
                           a_ebbcName, a_dombcName);

    // RHS of stencil equation
    Chombo4::LAPACKMatrix Qshift;
    getEBFluxMatrix(Qshift, a_face, faceLoc, a_ebnorm, a_dx);

    // Weights matrix
    Chombo4::LAPACKMatrix weights;
    getWeightsMatrixFace(weights, neighs, faceLoc);

    // Compute phiCoef
    Chombo4::LAPACKMatrix phiCoef;
    getStencilWeights(phiCoef, Qshift, Mvol, weights, a_dx);

    //Filling stencil
    a_stencil.clear();
    Real coef_acc = 0.0;
    if (m_verbose)  pout() << " Coefs Matrix " << endl;
    for(int ivof = 0; ivof < vols.size(); ivof++)
    {
      Real thisCo = phiCoef(0, ivof);
      a_stencil.add(vols[ivof], thisCo);
      if (m_verbose){
        pout()<<"stencil coefs for Vof: "<<vols[ivof].m_pt<<" Coef: "<<thisCo<<endl;
      }
      coef_acc += thisCo;
    }

    int nvolneigh = vols.size();
    a_inhomogTerm = computeInhomogTerm(phiCoef, neighs,
                                       nvolneigh,
                                       a_ebbcName, a_dombcName,
                                       a_dx);
    if (m_verbose) pout()<<"For phi=1 "<<coef_acc<<" "<<a_inhomogTerm<<" "<<coef_acc+a_inhomogTerm<<endl;
  }

//Compute rhs matrix for a face
  template<CENTERING cent, int order>
  void
  DivWLSFlux<cent, order>::
  getFaceFluxMatrix(Chombo4::LAPACKMatrix    & a_fluxMatrix,
                    const EBIndex<cent>      & a_face,
                    RealVect                 & a_faceLoc,
                    IndMomSDMinOne           & a_faceMom,
                    double                     a_signNormal,
                    double                     a_dx) const
  {
    PR_TIME("getFaceFluxMatrix");

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
        diffrv[index] = diff[idir]/a_dx;
        index++;
      }
    }
    a_faceMom.shift(diffrv);

    MomentIterator<DIM, order> momit;
    a_fluxMatrix.define(1, momit.size());
    a_fluxMatrix.setVal(0.);
    const int& faceDir = a_face.direction();
    if (this->m_verbose) pout() << " RHS - Face Flux Matrix " << endl;
    for(momit.reset(); momit.ok(); ++momit)
    {
      IndexTM<int, DIM> p = momit();
      int pd = p[faceDir];

      int jmom = IndexedMoments<DIM, order>::indexOf(p);

      //all other integrals are zero since along face
      //(y-ybar) =0  so (y-ybar)^py = 0 unless py == 0.
      if(pd == 0)
      {
        // pfaceshift are all components of p that are not in the faceDir direction
        IndexTM<int, DIM-1> pfaceshift;
        int index = 0;
        for(int idir = 0; idir < DIM; idir++)
        {
           if(idir != faceDir)
           {
             pfaceshift[index] = p[idir];
             index++;
           }
        }
        ::pout() << pfaceshift;
        // for readability, multiply by pd even though pd = 1
        //also pfaceshift is really p-e^d because we have set p
        a_fluxMatrix(0, jmom) = a_faceMom[pfaceshift];
      }
      else
      {
        a_fluxMatrix(0, jmom) = 0;
      }
      if(Chombo4::LAPACKMatrix::s_verbose)
      {
        ::pout() << ",icol = "  << jmom <<  ", Q(0,icol) = " << a_fluxMatrix(0, jmom) << endl;
      }
    }
  }

  //Compute rhs matrix for a EB face
  template<CENTERING cent, int order>
  void
  DivWLSFlux<cent, order>::
  getEBFluxMatrix(Chombo4::LAPACKMatrix    & a_fluxMatrix,
                  const EBIndex<cent>      & a_face,
                  RealVect                 & a_faceLoc,
                  Vector<IndMomDIM>        & a_ebnorm,
                  double                     a_dx) const
  {
    PR_TIME("getEBFluxMatrix");
    MomentIterator<DIM, order> momit;
    a_fluxMatrix.define(1, momit.size());
    a_fluxMatrix.setVal(0.);

    RealVect startLoc;
    for (int idir = 0; idir < DIM; idir++)
    {
      startLoc[idir] = a_dx * (0.5 + a_face.m_pt[idir]);
    }

    for(int idir = 0; idir < DIM; idir++)
    {

      RealVect diff = startLoc;
      diff -= a_faceLoc;
      IndexTM<Real, DIM> diffrv;
      for(int idir = 0; idir < DIM; idir++)
      {
        diffrv[idir] = diff[idir]/a_dx;
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
        if(pd >= 0)
        {
          IndexTM<int, DIM> pshift = p - BASISV_TM<int, DIM>(idir);
          fluxIncr(0, jmom) = a_ebnorm[idir][pshift];
        }
        else
        {
          fluxIncr(0, jmom) = 0;
        }
      }
      a_fluxMatrix += fluxIncr;
    }
    if (this->m_verbose)  pout() << " RHS - EB Flux Matrix " << endl;
    for(int idir = 0; idir < DIM; idir++)
    {
      for(momit.reset(); momit.ok(); ++momit)
      {
        IndexTM<int, DIM> p = momit();
        //Real pd = p[idir];
        //get column number associated with this particular moment
        int jmom = IndexedMoments<DIM, order>::indexOf(p);
        if(Chombo4::LAPACKMatrix::s_verbose)
        {
          ::pout() << ",icol = "  << jmom <<  ", Q(0,icol) = " << a_fluxMatrix(0, jmom) << endl;
        }
      }
    }

  }

  //Compute rhs matrix for a face
  template<int order>
  void
  BaseWLSCell<order>::
  getEvalCellMatrix(Chombo4::LAPACKMatrix    & a_cellMatrix,
                    const EBIndex<CELL>      & a_cell,
                    RealVect                 & a_cellLoc,
                    IndMomDIM                & a_cellMom,
                    double                     a_dx,
                    const IndexTM<int, DIM>    a_deriv) const
  {
    PR_TIME("getEvalCellMatrix");

    RealVect startLoc;
    for (int idir = 0; idir < DIM; idir++)
      {
        startLoc[idir] = a_dx * (0.5 + a_cell.m_pt[idir]);
      }
    // shift center of moments to a_cellLoc
    // the moments are centered at the cell center initially.
    // Move them to the centroid
    RealVect diff = startLoc - a_cellLoc;
    IndexTM<Real, DIM> diffrv;
    for(int idir = 0; idir < DIM; idir++)
      {
        diffrv[idir] = diff[idir]/a_dx;
      }
    a_cellMom.shift(diffrv);

    MomentIterator<DIM, order> momit;
    a_cellMatrix.define(1, momit.size());
    a_cellMatrix.setVal(0.);
    if (this->m_verbose) pout() << " RHS - Cell Eval Matrix " << endl;
    for(momit.reset(); momit.ok(); ++momit)
      {
        IndexTM<int, DIM> p = momit();
        int jmom = IndexedMoments<DIM, order>::indexOf(p);
        IndexTM<int, DIM> pmq = p - a_deriv; // p minus q
        if (pmq.min() < 0) // dropping bellow a zero moment has value of zero
          {
            a_cellMatrix(0, jmom) = 0.;
          }
        else
          {
            Real kappa = a_cellMom[0];
            Real derivCoef = p.factorial()/pmq.factorial();
            a_cellMatrix(0, jmom) = derivCoef*a_cellMom[p-a_deriv]/kappa;
          }
        if(Chombo4::LAPACKMatrix::s_verbose)
          {
            pout() << ",icol = "  << jmom <<  ", Q(0,icol) = " << a_cellMatrix(0, jmom) << endl;
          }
      }
  }

  //Compute face centered flux stencil
  template<int order>
  void
  BaseWLSCell<order>::
  getCellStencil(LocalStencil<CELL, double> & a_sten,
                 Real                       & a_inhomogTerm, //for inhomogeneous bcs
                 const EBIndex<CELL>        & a_vof,
                 const EBGraph              & a_graph,
                 const VoluData             & a_voludata,
                 const EBFaData             & a_ebfadata,
                 const XFacData             & a_xfacdata,
                 const YFacData             & a_yfacdata,
                 const ZFacData             & a_zfacdata,
                 const EBNormalData         & a_ebnormxdata,
                 const EBNormalData         & a_ebnormydata,
                 const EBNormalData         & a_ebnormzdata,
                 string                       a_ebbcName,
                 string                       a_dombcName[2*DIM],
                 IndMomDIM                 & a_cellMom,
                 double                       a_dx,
                 const IndexTM<int, DIM> a_deriv) const
  {
    PR_TIME("getCellStencil");
    Chombo4::LAPACKMatrix::s_verbose = this->m_verbose;

    Neighbors<order> neighs;
    EBIndex<CELL> cell = a_vof; // a dummy, since we want cell bases neighbors
    RealVect cellLoc;
    int neighType=1; // about the cell
    for (int idir = 0; idir < DIM; idir++)
      {
        cellLoc[idir] = a_dx * (0.5 + cell.m_pt[idir]);
      }

    Vector<EBIndex<CELL> >  vols = this->getVoFsNeighbors(cell, a_vof, this->m_radius, neighs, a_graph, neighType);
    MomentIterator<DIM,order> momit;
    if (vols.size()<=momit.size()) pout()<<"DANGER: NUMBER OF NEIGHS "<<vols.size()<<" and NUMBER OF coef "<<momit.size()<<endl;

    this->getMomentsInfoFromNeighbors(vols,
                                      a_graph, a_voludata, a_ebfadata,
                                      a_xfacdata, a_yfacdata, a_zfacdata,
                                      a_ebnormxdata, a_ebnormydata,
                                      a_ebnormzdata,
                                      neighs, neighType,
                                      a_ebbcName, a_dombcName, a_dx);

    //first fill M with moments,  shifted to the cell center
    Chombo4::LAPACKMatrix  Mvol;
    this->getShiftedMomentMatrix(Mvol, cellLoc, a_dx,
                                 neighs,
                                 a_ebbcName, a_dombcName);

    // RHS of stencil equation
    Chombo4::LAPACKMatrix Qshift;
    this->getEvalCellMatrix(Qshift, cell, cellLoc, a_cellMom, a_dx, a_deriv);

    // Weights matrix, actually weighted about the cell center
    Chombo4::LAPACKMatrix weights;
    this->getWeightsMatrixFace(weights, neighs, cellLoc);

    // Compute phiCoef
    Chombo4::LAPACKMatrix phiCoef;
    this->getStencilWeights(phiCoef, Qshift, Mvol, weights, a_dx);

    //Filling stencil
    a_sten.clear();
    Real coef_acc = 0.0;
    if (this->m_verbose)  pout() << " Cell Coefs Matrix " << endl;
    for(int ivof = 0; ivof < vols.size(); ivof++)
      {
        Real thisCo = phiCoef(0, ivof);
        a_sten.add(vols[ivof], thisCo);
        if (this->m_verbose){
          pout()<<"stencil coefs for Vof: "<<vols[ivof].m_pt<<" Coef: "<<thisCo<<endl;
        }
        coef_acc += thisCo;
      }

    int nvolneigh = vols.size();
    a_inhomogTerm = this->computeInhomogTerm(phiCoef, neighs,
                                             nvolneigh,
                                             a_ebbcName, a_dombcName,
                                             a_dx);
    if (this->m_verbose) pout()<<"For phi=1 "<<coef_acc<<" "<<a_inhomogTerm<<" "<<coef_acc+a_inhomogTerm<<endl;
  }

  //Getting flux stencil for laplacian operator at a face by a weighted average between 
  //two cell-centered flux stencils states computed by the neighbor cells sharing the face
  template < int order>
  LocalStencil< CELL,  double>
  WLSStencil<order>::
  getLapFluxStencil(LocalStencil< CELL,  double> & a_fluxMStencil,
                    LocalStencil< CELL,  double> & a_fluxPStencil, 
                    RealVect                     & a_centroid,
                    const EBIndex<CELL>          & a_voflo, 
                    const EBIndex<CELL>          & a_vofhi,
                    int                            a_dir,
                    const VoluData               & a_voludata,
                    const EBGraph                & a_graph,
                    double                         a_dx)
  {
    LocalStencil< CELL,  double> fluxStencil;
    RealVect centroidM;
    double kappaM;
    Normalizor<order>::normalizedVoluMoments(kappaM, centroidM, a_voludata, a_voflo, a_graph, a_dx);
    RealVect centroidP;
    double kappaP;
    Normalizor<order>::normalizedVoluMoments(kappaP, centroidP, a_voludata, a_vofhi, a_graph, a_dx);
    Real expD = 1.0;
    Real distM = 1.0;
    Real distP = 1.0;
    if (a_dir == 0 ){
      distM = pow(pow(pow(abs(centroidM[0]-0.5),2)+pow(abs(centroidM[1]-a_centroid[1]),2),0.5),expD);
      distP = pow(pow(pow(abs(centroidP[0]+0.5),2)+pow(abs(centroidP[1]-a_centroid[1]),2),0.5),expD);
    }else if (a_dir == 1){
      distM = pow(pow(pow(abs(centroidM[1]-0.5),2)+pow(abs(centroidM[0]-a_centroid[0]),2),0.5),expD);
      distP = pow(pow(pow(abs(centroidP[1]+0.5),2)+pow(abs(centroidP[0]-a_centroid[0]),2),0.5),expD);
    }
    //distM = pow(abs(centroidM[a_dir]-0.5),expD);
    //distP = pow(abs(centroidP[a_dir]+0.5),expD);
    a_fluxMStencil *= distP/(distM+distP);
    a_fluxPStencil *= distM/(distM+distP);
    fluxStencil += a_fluxMStencil;
    fluxStencil += a_fluxPStencil;
    return fluxStencil;
  }

  template < int order>
  double
  WLSStencil<order>::
  getLapinhomogTerm(double                       & a_inhomogTermM,
                    double                       & a_inhomogTermP, 
                    RealVect                     & a_centroid,
                    const EBIndex<CELL>          & a_voflo, 
                    const EBIndex<CELL>          & a_vofhi,
                    int                            a_dir,
                    const VoluData               & a_voludata,
                    const EBGraph                & a_graph,
                    double                         a_dx)
  {
    double inhomogTerm;
    RealVect centroidM;
    double kappaM;
    Normalizor<order>::normalizedVoluMoments(kappaM, centroidM, a_voludata, a_voflo, a_graph, a_dx);
    RealVect centroidP;
    double kappaP;
    Normalizor<order>::normalizedVoluMoments(kappaP, centroidP, a_voludata, a_vofhi, a_graph, a_dx);
    Real expD = 1.0;
    Real distM = 1.0;
    Real distP = 1.0;
    if (a_dir == 0 ){
      distM = pow(pow(pow(abs(centroidM[0]-0.5),2)+pow(abs(centroidM[1]-a_centroid[1]),2),0.5),expD);
      distP = pow(pow(pow(abs(centroidP[0]+0.5),2)+pow(abs(centroidP[1]-a_centroid[1]),2),0.5),expD);
    }else if (a_dir == 1){
      distM = pow(pow(pow(abs(centroidM[1]-0.5),2)+pow(abs(centroidM[0]-a_centroid[0]),2),0.5),expD);
      distP = pow(pow(pow(abs(centroidP[1]+0.5),2)+pow(abs(centroidP[0]-a_centroid[0]),2),0.5),expD);
    }
    //Real distM = pow(abs(centroidM[a_dir]-0.5),expD);
    //Real distP = pow(abs(centroidP[a_dir]+0.5),expD);
    a_inhomogTermM *= distP/(distM+distP);
    a_inhomogTermP *= distM/(distM+distP);
    inhomogTerm = a_inhomogTermM + a_inhomogTermP;
    return inhomogTerm;
  }

  //Gets Laplacian WLS flux divergence stencil
  template < int order>
  void
  WLSStencil<order>::
  getLaplaWLSDivFStencil(LocalStencil< CELL,  double>          & a_vofStencil,
                         const EBIndex<CELL>                   & a_vof,
                         const EBGraph                         & a_graph,
                         const VoluData                        & a_voludata,
                         const EBFaData                        & a_ebfadata,
                         const XFacData                        & a_xfacdata,
                         const YFacData                        & a_yfacdata,
                         const ZFacData                        & a_zfacdata,
                         const EBNormalData                    & a_ebnormxdata,
                         const EBNormalData                    & a_ebnormydata,
                         const EBNormalData                    & a_ebnormzdata,
                         const BaseWLSFlux<XFACE, order>  & a_xfaceflux,
                         const BaseWLSFlux<YFACE, order>  & a_yfaceflux,
                         const BaseWLSFlux<ZFACE, order>  & a_zfaceflux,
                         const BaseWLSFlux<BOUNDARY, order>  & a_ebflux,
                         string                             a_ebbcName,
                         string                             a_dombcName[2*DIM],
                         Real                               & a_inhoTerm,
                         int                                & a_neighType,
                         double                                  a_dx)
  {
    PR_TIME("getLaplaWLSDivFStencil");
    a_vofStencil.clear();
    Real volScaleFactor  = (D_TERM(a_dx, *a_dx, *a_dx));
    Real areaScaleFactor = (D_TERM(1.0,  *a_dx, *a_dx));

    for (SideIterator sit; sit.ok(); ++sit)
    {
      vector< EBIndex<XFACE> > xfaces = a_graph.getXFaces(a_vof, sit());
      int isign = sign(sit());
      for (int iface = 0; iface < xfaces.size(); iface++)
      {
        LocalStencil<CELL, double> fluxStencil;
        Real inhomogTerm=0.0; //for inhomogeneous bcs
        double factor = double(isign);
        std::map<EBIndex<XFACE>,LocalStencil<CELL, double> >::iterator it;
        it = xfacesFluxStencil.find(xfaces[iface]);
        if (it != xfacesFluxStencil.end())
        {
          fluxStencil = it->second;
          inhomogTerm = xfacesInhTerm[xfaces[iface]];
          fluxStencil *= factor;
          a_inhoTerm += inhomogTerm*factor;
          a_vofStencil += fluxStencil;
        }
        else
        {
          double   areaFrac;
          RealVect centroid;
          Normalizor<order>::normalizedFaceMoments(areaFrac, centroid, a_xfacdata, xfaces[iface], a_graph, a_dx);
  
          IndMomSDMinOne faceMom;
          if (areaFrac < 1.0) 
          {
            faceMom = a_xfacdata(xfaces[iface], 0); 
            for(MomentIterator<DIM-1, order> iter; iter.ok(); ++iter)
            {
              const IndexTM<int, DIM-1>& p = iter();
              faceMom[iter()] = faceMom[iter()]/(areaScaleFactor*POW(a_dx,p));
            }
          }
          else
          {
            faceMom.setToRegular(a_dx);
            for(MomentIterator<DIM-1, order> iter; iter.ok(); ++iter)
            {
              const IndexTM<int, DIM-1>& p = iter();
              faceMom[iter()] = faceMom[iter()]/(areaScaleFactor*POW(a_dx,p));
            }
          }
  
          if (a_xfaceflux.m_verbose) pout()<<" stencil Xface "<<xfaces[iface].m_pt<<endl;
  
          if (!xfaces[iface].m_isBoundary && a_neighType){
            LocalStencil<CELL, double> fluxMStencil;
            LocalStencil<CELL, double> fluxPStencil;
            double inhomogTermM;
            double inhomogTermP;
            EBIndex<CELL> voflo = xfaces[iface].getVoF(Side::Lo);
            EBIndex<CELL> vofhi = xfaces[iface].getVoF(Side::Hi);
  
            a_xfaceflux.getFaceCenteredStencil(fluxMStencil, inhomogTermM, voflo, xfaces[iface],
                                               a_graph,
                                               a_voludata,
                                               a_ebfadata,    
                                               a_xfacdata,
                                               a_yfacdata,
                                               a_zfacdata,
                                               a_ebnormxdata,
                                               a_ebnormydata,
                                               a_ebnormzdata,
                                               a_ebbcName,
                                               a_dombcName,
                                               0,
                                               faceMom,
                                               a_neighType, 
                                               factor, 
                                               a_dx); 
    
            a_xfaceflux.getFaceCenteredStencil(fluxPStencil, inhomogTermP, vofhi, xfaces[iface],
                                               a_graph,
                                               a_voludata,
                                               a_ebfadata,
                                               a_xfacdata,
                                               a_yfacdata,
                                               a_zfacdata,
                                               a_ebnormxdata,
                                               a_ebnormydata,
                                               a_ebnormzdata,
                                               a_ebbcName,
                                               a_dombcName,
                                               0,
                                               faceMom,
                                               a_neighType,
                                               factor, 
                                               a_dx);
  
            fluxStencil = getLapFluxStencil(fluxMStencil, fluxPStencil, centroid, voflo, vofhi, 0, a_voludata, a_graph, a_dx);
            inhomogTerm = getLapinhomogTerm(inhomogTermM, inhomogTermP, centroid, voflo, vofhi, 0, a_voludata, a_graph, a_dx);
            //if (a_xfaceflux.m_verbose)
            //{
            //  pout()<<" xflux"<<endl;
            //  pout()<<" weight for M "<<distP/(distM+distP)<<endl;
            //  pout()<<" weight for P "<<distM/(distM+distP)<<endl;
            //}
  
          }else{
            a_xfaceflux.getFaceCenteredStencil(fluxStencil, inhomogTerm, a_vof, xfaces[iface],
                                               a_graph,
                                               a_voludata,
                                               a_ebfadata,
                                               a_xfacdata,
                                               a_yfacdata,
                                               a_zfacdata,
                                               a_ebnormxdata,
                                               a_ebnormydata,
                                               a_ebnormzdata,
                                               a_ebbcName,
                                               a_dombcName,
                                               0,
                                               faceMom,
                                               a_neighType,
                                               factor, 
                                               a_dx);
          }
  
          if (a_xfaceflux.m_verbose)
          {       
            pout()<<" inhomogTerm Xface "<<inhomogTerm<<endl;
            pout()<<" fluxStencil "<<endl;
            fluxStencil.print();
          }
          xfacesFluxStencil[xfaces[iface]] = fluxStencil;
          xfacesInhTerm[xfaces[iface]] = inhomogTerm; 
          fluxStencil *= factor;
          a_inhoTerm += inhomogTerm*factor;
          a_vofStencil += fluxStencil;
        }
      }
    }

    for (SideIterator sit; sit.ok(); ++sit)
    {
      vector< EBIndex<YFACE> > yfaces = a_graph.getYFaces(a_vof, sit());
      int isign = sign(sit());
      for (int iface = 0; iface < yfaces.size(); iface++)
      {
        LocalStencil<CELL, double> fluxStencil;
        Real inhomogTerm = 0.0; //for inhomogeneous bcs
        double factor = double(isign);
        std::map<EBIndex<YFACE>,LocalStencil<CELL, double> >::iterator it;
        it = yfacesFluxStencil.find(yfaces[iface]);
        if (it != yfacesFluxStencil.end())
        {
          fluxStencil = it->second;
          inhomogTerm = yfacesInhTerm[yfaces[iface]];
          fluxStencil *= factor;
          a_inhoTerm += inhomogTerm*factor;
          a_vofStencil += fluxStencil;
        }
        else
        {       
          double   areaFrac;
          RealVect centroid;
          Normalizor<order>::normalizedFaceMoments(areaFrac, centroid, a_yfacdata, yfaces[iface], a_graph, a_dx);
          IndMomSDMinOne faceMom;
          if (areaFrac < 1.0) 
          {
            faceMom = a_yfacdata(yfaces[iface], 0);
            for(MomentIterator<DIM-1, order> iter; iter.ok(); ++iter)
            {
              const IndexTM<int, DIM-1>& p = iter();
              faceMom[iter()] = faceMom[iter()]/(areaScaleFactor*POW(a_dx,p));
            }
          }
          else
          {
            faceMom.setToRegular(a_dx);
            for(MomentIterator<DIM-1, order> iter; iter.ok(); ++iter)
            {
              const IndexTM<int, DIM-1>& p = iter();
              faceMom[iter()] = faceMom[iter()]/(areaScaleFactor*POW(a_dx,p));
            }
          }
  
          if (a_yfaceflux.m_verbose) pout()<<" stencilY face "<<yfaces[iface].m_pt<<endl;
  
          if (!yfaces[iface].m_isBoundary && a_neighType){
            LocalStencil<CELL, double> fluxMStencil;
            LocalStencil<CELL, double> fluxPStencil;
            double inhomogTermM;
            double inhomogTermP;
            EBIndex<CELL> voflo = yfaces[iface].getVoF(Side::Lo);
            EBIndex<CELL> vofhi = yfaces[iface].getVoF(Side::Hi);
  
            a_yfaceflux.getFaceCenteredStencil(fluxMStencil, inhomogTermM, voflo, yfaces[iface],
                                           a_graph,
                                           a_voludata,
                                           a_ebfadata,
                                           a_xfacdata,
                                           a_yfacdata,
                                           a_zfacdata,
                                           a_ebnormxdata,
                                           a_ebnormydata,
                                           a_ebnormzdata,
                                           a_ebbcName,
                                           a_dombcName,
                                           1,
                                           faceMom,
                                           a_neighType,
                                           factor, 
                                           a_dx);
    
            a_yfaceflux.getFaceCenteredStencil(fluxPStencil, inhomogTermP, vofhi, yfaces[iface],
                                           a_graph,
                                           a_voludata,
                                           a_ebfadata,
                                           a_xfacdata,
                                           a_yfacdata,
                                           a_zfacdata,
                                           a_ebnormxdata,
                                           a_ebnormydata,
                                           a_ebnormzdata,
                                           a_ebbcName,
                                           a_dombcName,
                                           1,
                                           faceMom,
                                           a_neighType,
                                           factor, 
                                           a_dx);
            fluxStencil = getLapFluxStencil(fluxMStencil, fluxPStencil, centroid, voflo, vofhi, 1, a_voludata, a_graph, a_dx);
            inhomogTerm = getLapinhomogTerm(inhomogTermM, inhomogTermP, centroid, voflo, vofhi, 1, a_voludata, a_graph, a_dx);
            //if (a_yfaceflux.m_verbose)
            //{
            //  pout()<<" yflux"<<endl;
            //  pout()<<" weight for M "<<distP/(distM+distP)<<endl;
            //  pout()<<" weight for P "<<distM/(distM+distP)<<endl;
            //}
          }else{
            a_yfaceflux.getFaceCenteredStencil(fluxStencil, inhomogTerm, a_vof, yfaces[iface],
                                               a_graph,
                                               a_voludata,
                                               a_ebfadata,
                                               a_xfacdata,
                                               a_yfacdata,
                                               a_zfacdata,
                                               a_ebnormxdata,
                                               a_ebnormydata,
                                               a_ebnormzdata,
                                               a_ebbcName,
                                               a_dombcName,
                                               1,
                                               faceMom, 
                                               a_neighType, 
                                               factor, 
                                               a_dx);
          }
   
          if (a_yfaceflux.m_verbose)
          {
            pout()<<" inhomogTerm Yface "<<inhomogTerm<<endl;
            pout()<<" fluxStencil "<<endl;
            fluxStencil.print();
          }          
          yfacesFluxStencil[yfaces[iface]] = fluxStencil;
          yfacesInhTerm[yfaces[iface]] = inhomogTerm;
          fluxStencil *= factor;
          a_inhoTerm += inhomogTerm*factor;
          a_vofStencil += fluxStencil;
        }
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
        Real inhomogTerm; //for inhomogeneous bcs
        double factor = double(isign);
        std::map<EBIndex<ZFACE>,LocalStencil<CELL, double> >::iterator it;
        it = zfacesFluxStencil.find(zfaces[iface]);
        if (it != zfacesFluxStencil.end())
        {
          fluxStencil = it->second;
          inhomogTerm = zfacesInhTerm[zfaces[iface]];
          fluxStencil *= factor;
          a_inhoTerm += inhomogTerm*factor;
          a_vofStencil += fluxStencil;
        }
        else
        {
          double   areaFrac;
          RealVect centroid;
          Normalizor<order>::normalizedFaceMoments(areaFrac, centroid, a_zfacdata, zfaces[iface], a_graph, a_dx);
          IndMomSDMinOne faceMom;
          if (areaFrac < 1.0) 
          {
            faceMom = a_zfacdata(zfaces[iface], 0);
            for(MomentIterator<DIM-1, order> iter; iter.ok(); ++iter)
            {
              const IndexTM<int, DIM-1>& p = iter();
              faceMom[iter()] = faceMom[iter()]/(areaScaleFactor*POW(a_dx,p));
            } 
          }
          else
          {
            faceMom.setToRegular(a_dx);
            for(MomentIterator<DIM-1, order> iter; iter.ok(); ++iter)
            {
              const IndexTM<int, DIM-1>& p = iter();
              faceMom[iter()] = faceMom[iter()]/(areaScaleFactor*POW(a_dx,p));
            }
          }
  
          if (!zfaces[iface].m_isBoundary && a_neighType){
            LocalStencil<CELL, double> fluxMStencil;
            LocalStencil<CELL, double> fluxPStencil;
            double inhomogTermM;
            double inhomogTermP;
            EBIndex<CELL> voflo = zfaces[iface].getVoF(Side::Lo);
            EBIndex<CELL> vofhi = zfaces[iface].getVoF(Side::Hi);
  
            a_zfaceflux.getFaceCenteredStencil(fluxMStencil, inhomogTermM, voflo, zfaces[iface],
                                           a_graph,
                                           a_voludata,
                                           a_ebfadata,
                                           a_xfacdata,
                                           a_yfacdata,
                                           a_zfacdata,
                                           a_ebnormxdata,
                                           a_ebnormydata,
                                           a_ebnormzdata,
                                           a_ebbcName,
                                           a_dombcName,
                                           2,
                                           faceMom,
                                           a_neighType,
                                           factor, 
                                           a_dx);
    
            a_zfaceflux.getFaceCenteredStencil(fluxPStencil, inhomogTermP, vofhi, zfaces[iface],
                                           a_graph,
                                           a_voludata,
                                           a_ebfadata,
                                           a_xfacdata,
                                           a_yfacdata,
                                           a_zfacdata,
                                           a_ebnormxdata,
                                           a_ebnormydata,
                                           a_ebnormzdata,
                                           a_ebbcName,
                                           a_dombcName,
                                           2,
                                           faceMom,
                                           a_neighType,
                                           factor, 
                                           a_dx);
            fluxStencil = getLapFluxStencil(fluxMStencil, fluxPStencil, centroid, voflo, vofhi, 2, a_voludata, a_graph, a_dx);
            inhomogTerm = getLapinhomogTerm(inhomogTermM, inhomogTermP, centroid, voflo, vofhi, 2, a_voludata, a_graph, a_dx);
          }else{
            a_zfaceflux.getFaceCenteredStencil(fluxStencil, inhomogTerm, a_vof, zfaces[iface],
                                               a_graph,    
                                               a_voludata,
                                               a_ebfadata,
                                               a_xfacdata,
                                               a_yfacdata,
                                               a_zfacdata,
                                               a_ebnormxdata,
                                               a_ebnormydata,
                                               a_ebnormzdata,
                                               a_ebbcName,
                                               a_dombcName,
                                               2,
                                               faceMom,
                                               a_neighType, 
                                               factor, 
                                               a_dx);       
  
          }
          zfacesFluxStencil[zfaces[iface]] = fluxStencil;
          zfacesInhTerm[zfaces[iface]] = inhomogTerm;
          fluxStencil *= factor;
          a_inhoTerm += inhomogTerm*factor;
          a_vofStencil += fluxStencil;
        }
      }
    }
#endif      

    LocalStencil<CELL, double> ebflux;
    double bndryArea, alpha;
    Real inhomogTerm = 0.0; //for inhomogeneous bcs
    if(a_graph.isIrregular(a_vof.m_pt))
    {
      EBIndex<BOUNDARY> cutface = a_vof.getCutFace();
      RealVect centroid, normal;
      Vector<IndMomDIM> ebnorm (DIM);
      ebnorm[0] = a_ebnormxdata(a_vof, 0);
      for(MomentIterator<DIM, order> iter; iter.ok(); ++iter)
      {
        const IndexTM<int, DIM>& p = iter();
        ebnorm[0][iter()] = ebnorm[0][iter()]/(areaScaleFactor*POW(a_dx,p));
      }
      ebnorm[1] = a_ebnormydata(a_vof, 0);
      for(MomentIterator<DIM, order> iter; iter.ok(); ++iter)
      {
        const IndexTM<int, DIM>& p = iter();
        ebnorm[1][iter()] = ebnorm[1][iter()]/(areaScaleFactor*POW(a_dx,p));
      }
#if DIM==3
      ebnorm[2] = a_ebnormzdata(a_vof, 0);
      for(MomentIterator<DIM, order> iter; iter.ok(); ++iter)
      {
        const IndexTM<int, DIM>& p = iter();
        ebnorm[2][iter()] = ebnorm[2][iter()]/(areaScaleFactor*POW(a_dx,p));
      }
#endif
      Normalizor<order>::normalizedBndryMoments(alpha, centroid, a_ebfadata,  cutface, a_graph, a_dx);
      Normalizor<order>::getNormalAndBoundaryArea(bndryArea, normal,
                                                  a_voludata,
                                                  a_ebfadata,
                                                  a_xfacdata,
                                                  a_yfacdata,
                                                  a_zfacdata,
                                                  a_vof, a_graph, a_dx);

      if (a_ebflux.m_verbose) pout()<<" stencil ebface "<<a_vof.m_pt<<endl;

      if(a_ebbcName == string("Neumann"))
      {
        //neumann has no stencil
        RealVect faceLoc;
        for (int idir = 0; idir < DIM; idir++)
        {
          faceLoc[idir] = a_dx * (0.5 + cutface.m_pt[idir]);
        }
        inhomogTerm = NeumannExactSoln<order>::
                               getEBIntegral(ebnorm,
                                             faceLoc,
                                             a_ebflux.m_exactSol,
                                             a_dx);
        inhomogTerm = inhomogTerm/(areaScaleFactor);

      }
      else if(a_ebbcName == string("Dirichlet"))
      {
        a_ebflux.getEBFluxStencil(ebflux, inhomogTerm, a_vof,
                                  cutface,
                                  a_graph,
                                  a_voludata,
                                  a_ebfadata,
                                  a_xfacdata,
                                  a_yfacdata,
                                  a_zfacdata,
                                  a_ebnormxdata,
                                  a_ebnormydata,
                                  a_ebnormzdata,
                                  ebnorm,
                                  a_ebbcName,
                                  a_dombcName,
                                  a_neighType, 
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
      if (a_ebflux.m_verbose)
      {
        pout()<<" inhomogTerm ebface "<<inhomogTerm<<endl;
        pout()<<" alpha "<<alpha<<" bndryArea "<<bndryArea<<endl;
        pout()<<" ebfluxStencil "<<endl;
        ebflux.print();
      }
      a_inhoTerm += inhomogTerm;
      a_vofStencil += ebflux;
    }
  }

  //Getting flux stencil for advection operator at a face by a weighted average between 
  //two cell-centered flux stencils states computed by the neighbor cells sharing the face
  template < int order>
  LocalStencil< CELL,  double>
  WLSStencil<order>::
  getAdvFluxStencil(LocalStencil< CELL,  double> & a_fluxMStencil,
                    LocalStencil< CELL,  double> & a_fluxPStencil, 
                    RealVect                     & a_centroid,
                    const EBIndex<CELL>          & a_voflo, 
                    const EBIndex<CELL>          & a_vofhi,
                    const RealVect               & a_v_faceLoc,
                    int                            a_dir,
                    const VoluData               & a_voludata,
                    const EBGraph                & a_graph,
                    double                         a_dx)
  {
    LocalStencil< CELL,  double> fluxStencil;
    RealVect centroidM;
    double kappaM;
    Normalizor<order>::normalizedVoluMoments(kappaM, centroidM, a_voludata, a_voflo, a_graph, a_dx);
    RealVect centroidP;
    double kappaP;
    Normalizor<order>::normalizedVoluMoments(kappaP, centroidP, a_voludata, a_vofhi, a_graph, a_dx);
    Real expD = 1.0;
    Real distM = 1.0;
    Real distP = 1.0;
#if DIM==3
    if (a_dir == 0 ){
      distM = pow(pow(pow(abs(centroidM[0]-0.5),2)+pow(abs(centroidM[1]-a_centroid[1]),2)+pow(abs(centroidM[2]-a_centroid[2]),2),0.5),expD);
      distP = pow(pow(pow(abs(centroidP[0]+0.5),2)+pow(abs(centroidP[1]-a_centroid[1]),2)+pow(abs(centroidP[2]-a_centroid[2]),2),0.5),expD);
    }else if (a_dir == 1){
      distM = pow(pow(pow(abs(centroidM[1]-0.5),2)+pow(abs(centroidM[0]-a_centroid[0]),2)+pow(abs(centroidM[2]-a_centroid[2]),2),0.5),expD);
      distP = pow(pow(pow(abs(centroidP[1]+0.5),2)+pow(abs(centroidP[0]-a_centroid[0]),2)+pow(abs(centroidP[2]-a_centroid[2]),2),0.5),expD);
    }else if (a_dir == 2){
      distM = pow(pow(pow(abs(centroidM[2]-0.5),2)+pow(abs(centroidM[0]-a_centroid[0]),2)+pow(abs(centroidM[1]-a_centroid[1]),2),0.5),expD);
      distP = pow(pow(pow(abs(centroidP[2]+0.5),2)+pow(abs(centroidP[0]-a_centroid[0]),2)+pow(abs(centroidP[1]-a_centroid[1]),2),0.5),expD);
    }
#elif DIM==2
    if (a_dir == 0 ){
      distM = pow(pow(pow(abs(centroidM[0]-0.5),2)+pow(abs(centroidM[1]-a_centroid[1]),2),0.5),expD);
      distP = pow(pow(pow(abs(centroidP[0]+0.5),2)+pow(abs(centroidP[1]-a_centroid[1]),2),0.5),expD);
    }else if (a_dir == 1){
      distM = pow(pow(pow(abs(centroidM[1]-0.5),2)+pow(abs(centroidM[0]-a_centroid[0]),2),0.5),expD);
      distP = pow(pow(pow(abs(centroidP[1]+0.5),2)+pow(abs(centroidP[0]-a_centroid[0]),2),0.5),expD);
    }
#endif
    bool upwind_stencil = false;
    //This will take the upwind stencil between stencilM and stencilP
    if (upwind_stencil)
    {
      if (a_v_faceLoc[a_dir]>1e-9){
        fluxStencil += a_fluxMStencil;
      }else if (a_v_faceLoc[a_dir]<-1e-9){
        fluxStencil += a_fluxPStencil;
      }else{
        a_fluxMStencil *= 0.5;
        a_fluxPStencil *= 0.5;
        fluxStencil += a_fluxMStencil;
        fluxStencil += a_fluxPStencil;
      }
    }
    //This will take the centroid average stencil between stencilM and stencilP
    else
    {
      a_fluxMStencil *= distP/(distM+distP);
      a_fluxPStencil *= distM/(distM+distP);
      fluxStencil += a_fluxMStencil;
      fluxStencil += a_fluxPStencil;
    } 
    return fluxStencil;
  }

  //Gets Advection Laplacian WLS flux divergence stencil
  template < int order>
  void
  WLSStencil<order>::
  getAdvWLSDivFStencil(LocalStencil< CELL,  double>          & a_vofStencil,
                         const EBIndex<CELL>                   & a_vof,
                         const EBGraph                         & a_graph,
                         const VoluData                        & a_voludata,
                         const EBFaData                        & a_ebfadata,
                         const XFacData                        & a_xfacdata,
                         const YFacData                        & a_yfacdata,
                         const ZFacData                        & a_zfacdata,
                         const EBNormalData                    & a_ebnormxdata,
                         const EBNormalData                    & a_ebnormydata,
                         const EBNormalData                    & a_ebnormzdata,
                         const ExactVelAdvectionWLSFlux<XFACE, order>  & a_xfaceflux,
                         const ExactVelAdvectionWLSFlux<YFACE, order>  & a_yfaceflux,
                         const ExactVelAdvectionWLSFlux<ZFACE, order>  & a_zfaceflux,
                         const ExactVelAdvectionWLSFlux<BOUNDARY, order>  & a_ebflux,
                         string                             a_ebbcName,
                         string                             a_dombcName[2*DIM],
                         Real                               & a_inhoTerm,
                         int                                & a_neighType,
                         double                                  a_dx)
  {
    PR_TIME("getAdvWLSDivFStencil");
    a_vofStencil.clear();
    Real volScaleFactor  = (D_TERM(a_dx, *a_dx, *a_dx));
    Real areaScaleFactor = (D_TERM(1.0,  *a_dx, *a_dx));

    for (SideIterator sit; sit.ok(); ++sit)
    {
      vector< EBIndex<XFACE> > xfaces = a_graph.getXFaces(a_vof, sit());
      int isign = sign(sit());
      for (int iface = 0; iface < xfaces.size(); iface++)
      {
        LocalStencil<CELL, double> fluxMStencil;
        LocalStencil<CELL, double> fluxPStencil;
        LocalStencil<CELL, double> fluxStencil;
        Real inhomogTerm=0.0; //for inhomogeneous bcs
        double factor = double(isign);
        double   areaFrac;
        RealVect centroid;
        Normalizor<order>::normalizedFaceMoments(areaFrac, centroid, a_xfacdata, xfaces[iface], a_graph, a_dx);
        std::map<EBIndex<XFACE>,LocalStencil<CELL, double> >::iterator it;
        it = xfacesFluxStencil.find(xfaces[iface]);
        if (it != xfacesFluxStencil.end())
        {
          fluxStencil = it->second;
          inhomogTerm = xfacesInhTerm[xfaces[iface]];
          fluxStencil *= factor;
          a_inhoTerm += inhomogTerm*factor;
          a_vofStencil += fluxStencil;
        }
        else
        {
          EBIndex<CELL> voflo = xfaces[iface].getVoF(Side::Lo);
          EBIndex<CELL> vofhi = xfaces[iface].getVoF(Side::Hi);
  
          RealVect faceLoc;
          for (int idir = 0; idir < DIM; idir++)
          {
            faceLoc[idir] = a_dx * (0.5 + xfaces[iface].m_pt[idir]);
            if (idir == 0) faceLoc[idir] = a_dx * (xfaces[iface].m_pt[idir]);
          }
          RealVect v_faceLoc = a_xfaceflux.m_exactVel->getVelocityAtPoint(faceLoc);
  
          IndMomSDMinOne faceMom;
          if (areaFrac < 1.0) 
          {
            faceMom = a_xfacdata(xfaces[iface], 0); 
            for(MomentIterator<DIM-1, order> iter; iter.ok(); ++iter)
            {
              const IndexTM<int, DIM-1>& p = iter();
              faceMom[iter()] = faceMom[iter()]/(areaScaleFactor*POW(a_dx,p));
            }
          }
          else
          {
            faceMom.setToRegular(a_dx);
            for(MomentIterator<DIM-1, order> iter; iter.ok(); ++iter)
            {
              const IndexTM<int, DIM-1>& p = iter();
              faceMom[iter()] = faceMom[iter()]/(areaScaleFactor*POW(a_dx,p));
            }
          }
  
          if (a_xfaceflux.m_verbose) pout()<<" stencil Xface "<<xfaces[iface].m_pt<<" centroid "<<(a_dx*centroid)+faceLoc<<endl;
  
          if (!xfaces[iface].m_isBoundary){

            a_xfaceflux.getFaceCenteredStencil(fluxMStencil, inhomogTerm, voflo, xfaces[iface],
                                               a_graph,
                                               a_voludata,
                                               a_ebfadata,    
                                               a_xfacdata,
                                               a_yfacdata,
                                               a_zfacdata,
                                               a_ebnormxdata,
                                               a_ebnormydata,
                                               a_ebnormzdata,
                                               a_ebbcName,
                                               a_dombcName,
                                               0,
                                               faceMom,
                                               a_neighType,
                                               factor, 
                                               a_dx); 

            a_xfaceflux.getFaceCenteredStencil(fluxPStencil, inhomogTerm, vofhi, xfaces[iface],
                                               a_graph,
                                               a_voludata,
                                               a_ebfadata,
                                               a_xfacdata,
                                               a_yfacdata,
                                               a_zfacdata,
                                               a_ebnormxdata,
                                               a_ebnormydata,
                                               a_ebnormzdata,
                                               a_ebbcName,
                                               a_dombcName,
                                               0,
                                               faceMom,
                                               a_neighType,
                                               factor, 
                                               a_dx);

            fluxStencil = getAdvFluxStencil(fluxMStencil, fluxPStencil, centroid, voflo, vofhi, v_faceLoc, 0, a_voludata, a_graph, a_dx);
  
          }else{
            a_xfaceflux.getFaceCenteredStencil(fluxStencil, inhomogTerm, a_vof, xfaces[iface],
                                               a_graph,
                                               a_voludata,
                                               a_ebfadata,
                                               a_xfacdata,
                                               a_yfacdata,
                                               a_zfacdata,
                                               a_ebnormxdata,
                                               a_ebnormydata,
                                               a_ebnormzdata,
                                               a_ebbcName,
                                               a_dombcName,
                                               0,
                                               faceMom,
                                               a_neighType,
                                               factor, 
                                               a_dx);
          }
          if (a_xfaceflux.m_verbose)
          {       
            pout()<<" inhomogTerm Xface "<<inhomogTerm<<endl;
            pout()<<" fluxStencil "<<endl;
            fluxStencil.print();
          }
          //xfacesFluxStencil[xfaces[iface]] = fluxStencil;
          //xfacesInhTerm[xfaces[iface]] = inhomogTerm;
          //fluxStencil *= factor;
          a_inhoTerm += inhomogTerm;
          a_vofStencil += fluxStencil;
       
        }
      }
    }

    for (SideIterator sit; sit.ok(); ++sit)
    {
      vector< EBIndex<YFACE> > yfaces = a_graph.getYFaces(a_vof, sit());
      int isign = sign(sit());
      for (int iface = 0; iface < yfaces.size(); iface++)
      {
        LocalStencil<CELL, double> fluxMStencil;
        LocalStencil<CELL, double> fluxPStencil;
        LocalStencil<CELL, double> fluxStencil;
        Real inhomogTerm=0.0; //for inhomogeneous bcs
        double factor = double(isign);
        double   areaFrac;
        RealVect centroid;
        Normalizor<order>::normalizedFaceMoments(areaFrac, centroid, a_yfacdata, yfaces[iface], a_graph, a_dx);
        std::map<EBIndex<YFACE>,LocalStencil<CELL, double> >::iterator it;
        it = yfacesFluxStencil.find(yfaces[iface]);
        if (it != yfacesFluxStencil.end())
        {
          fluxStencil = it->second;
          inhomogTerm = yfacesInhTerm[yfaces[iface]];
          fluxStencil *= factor;
          a_inhoTerm += inhomogTerm*factor;
          a_vofStencil += fluxStencil;
        }
        else
        {
          EBIndex<CELL> voflo = yfaces[iface].getVoF(Side::Lo);
          EBIndex<CELL> vofhi = yfaces[iface].getVoF(Side::Hi);
  
          RealVect faceLoc;
          for (int idir = 0; idir < DIM; idir++)
          {
            faceLoc[idir] = a_dx * (0.5 + yfaces[iface].m_pt[idir]);
            if (idir == 1) faceLoc[idir] = a_dx * (yfaces[iface].m_pt[idir]);
          }
          RealVect v_faceLoc = a_yfaceflux.m_exactVel->getVelocityAtPoint(faceLoc);
  
          IndMomSDMinOne faceMom;
          if (areaFrac < 1.0) 
          {
            faceMom = a_yfacdata(yfaces[iface], 0);
            for(MomentIterator<DIM-1, order> iter; iter.ok(); ++iter)
            {
              const IndexTM<int, DIM-1>& p = iter();
              faceMom[iter()] = faceMom[iter()]/(areaScaleFactor*POW(a_dx,p));
            }
          }
          else
          {
            faceMom.setToRegular(a_dx);
            for(MomentIterator<DIM-1, order> iter; iter.ok(); ++iter)
            {
              const IndexTM<int, DIM-1>& p = iter();
              faceMom[iter()] = faceMom[iter()]/(areaScaleFactor*POW(a_dx,p));
            }
          }
  
          if (a_yfaceflux.m_verbose) pout()<<" stencil Yface "<<yfaces[iface].m_pt<<" centroid "<<(a_dx*centroid)+faceLoc<<endl;
  
          if (!yfaces[iface].m_isBoundary){
            a_yfaceflux.getFaceCenteredStencil(fluxMStencil, inhomogTerm, voflo, yfaces[iface],
                                           a_graph,
                                           a_voludata,
                                           a_ebfadata,
                                           a_xfacdata,
                                           a_yfacdata,
                                           a_zfacdata,
                                           a_ebnormxdata,
                                           a_ebnormydata,
                                           a_ebnormzdata,
                                           a_ebbcName,
                                           a_dombcName,
                                           1,
                                           faceMom,
                                           a_neighType,
                                           factor, 
                                           a_dx);

            a_yfaceflux.getFaceCenteredStencil(fluxPStencil, inhomogTerm, vofhi, yfaces[iface],
                                           a_graph,
                                           a_voludata,
                                           a_ebfadata,
                                           a_xfacdata,
                                           a_yfacdata,
                                           a_zfacdata,
                                           a_ebnormxdata,
                                           a_ebnormydata,
                                           a_ebnormzdata,
                                           a_ebbcName,
                                           a_dombcName,
                                           1,
                                           faceMom,
                                           a_neighType,
                                           factor, 
                                           a_dx);

            fluxStencil = getAdvFluxStencil(fluxMStencil, fluxPStencil, centroid, voflo, vofhi, v_faceLoc, 1, a_voludata, a_graph, a_dx);
          }else{
            a_yfaceflux.getFaceCenteredStencil(fluxStencil, inhomogTerm, a_vof, yfaces[iface],
                                               a_graph,
                                               a_voludata,
                                               a_ebfadata,
                                               a_xfacdata,
                                               a_yfacdata,
                                               a_zfacdata,
                                               a_ebnormxdata,
                                               a_ebnormydata,
                                               a_ebnormzdata,
                                               a_ebbcName,
                                               a_dombcName,
                                               1,
                                               faceMom, 
                                               a_neighType, 
                                               factor, 
                                               a_dx);
          }
          if (a_yfaceflux.m_verbose)
          {
            pout()<<" inhomogTerm Yface "<<inhomogTerm<<endl;
            pout()<<" fluxStencil "<<endl;
            fluxStencil.print();
          }
          //yfacesFluxStencil[yfaces[iface]] = fluxStencil;
          //yfacesInhTerm[yfaces[iface]] = inhomogTerm;
          //fluxStencil *= factor;
          a_inhoTerm += inhomogTerm;
          a_vofStencil += fluxStencil;
        }
      }
    }
#if DIM==3
    for (SideIterator sit; sit.ok(); ++sit)
    {
      vector< EBIndex<ZFACE> > zfaces = a_graph.getZFaces(a_vof, sit());
      int isign = sign(sit());
      for (int iface = 0; iface < zfaces.size(); iface++)
      {
        LocalStencil<CELL, double> fluxMStencil;
        LocalStencil<CELL, double> fluxPStencil;
        LocalStencil<CELL, double> fluxStencil;
        Real inhomogTerm=0.0; //for inhomogeneous bcs
        double factor = double(isign);
        std::map<EBIndex<ZFACE>,LocalStencil<CELL, double> >::iterator it;
        it = zfacesFluxStencil.find(zfaces[iface]);
        if (it != zfacesFluxStencil.end())
        {
          fluxStencil = it->second;
          inhomogTerm = zfacesInhTerm[zfaces[iface]];
          fluxStencil *= factor;
          a_inhoTerm += inhomogTerm*factor;
          a_vofStencil += fluxStencil;
        }
        else
        {
          EBIndex<CELL> voflo = zfaces[iface].getVoF(Side::Lo);
          EBIndex<CELL> vofhi = zfaces[iface].getVoF(Side::Hi); 
          double   areaFrac;
          RealVect centroid;
          Normalizor<order>::normalizedFaceMoments(areaFrac, centroid, a_zfacdata, zfaces[iface], a_graph, a_dx);
          RealVect faceLoc;
          for (int idir = 0; idir < DIM; idir++)
          {
            faceLoc[idir] = a_dx * (0.5 + zfaces[iface].m_pt[idir]);
            if (idir == 2) faceLoc[idir] = a_dx * (zfaces[iface].m_pt[idir]);
          }
          RealVect v_faceLoc = a_zfaceflux.m_exactVel->getVelocityAtPoint(faceLoc);
  
          IndMomSDMinOne faceMom;
          if (areaFrac < 1.0) 
          {
            faceMom = a_zfacdata(zfaces[iface], 0);
            for(MomentIterator<DIM-1, order> iter; iter.ok(); ++iter)
            {
              const IndexTM<int, DIM-1>& p = iter();
              faceMom[iter()] = faceMom[iter()]/(areaScaleFactor*POW(a_dx,p));
            } 
          }
          else
          {
            faceMom.setToRegular(a_dx);
            for(MomentIterator<DIM-1, order> iter; iter.ok(); ++iter)
            {
              const IndexTM<int, DIM-1>& p = iter();
              faceMom[iter()] = faceMom[iter()]/(areaScaleFactor*POW(a_dx,p));
            }
          }
  
          if (a_zfaceflux.m_verbose) pout()<<" stencil Zface "<<zfaces[iface].m_pt<<" centroid "<<(a_dx*centroid)+faceLoc<<endl;

          if (!zfaces[iface].m_isBoundary){
            a_zfaceflux.getFaceCenteredStencil(fluxMStencil, inhomogTerm, voflo, zfaces[iface],
                                           a_graph,
                                           a_voludata,
                                           a_ebfadata,
                                           a_xfacdata,
                                           a_yfacdata,
                                           a_zfacdata,
                                           a_ebnormxdata,
                                           a_ebnormydata,
                                           a_ebnormzdata,
                                           a_ebbcName,
                                           a_dombcName,
                                           2,
                                           faceMom,
                                           a_neighType,
                                           factor, 
                                           a_dx);
    
            a_zfaceflux.getFaceCenteredStencil(fluxPStencil, inhomogTerm, vofhi, zfaces[iface],
                                           a_graph,
                                           a_voludata,
                                           a_ebfadata,
                                           a_xfacdata,
                                           a_yfacdata,
                                           a_zfacdata,
                                           a_ebnormxdata,
                                           a_ebnormydata,
                                           a_ebnormzdata,
                                           a_ebbcName,
                                           a_dombcName,
                                           2,
                                           faceMom,
                                           a_neighType,
                                           factor, 
                                           a_dx);
            fluxStencil = getAdvFluxStencil(fluxMStencil, fluxPStencil, centroid, voflo, vofhi, v_faceLoc, 2, a_voludata, a_graph, a_dx);
          }else{
            a_zfaceflux.getFaceCenteredStencil(fluxStencil, inhomogTerm, a_vof, zfaces[iface],
                                               a_graph,    
                                               a_voludata,
                                               a_ebfadata,
                                               a_xfacdata,
                                               a_yfacdata,
                                               a_zfacdata,
                                               a_ebnormxdata,
                                               a_ebnormydata,
                                               a_ebnormzdata,
                                               a_ebbcName,
                                               a_dombcName,
                                               2,
                                               faceMom,
                                               a_neighType, 
                                               factor, 
                                               a_dx);       
  
          }
          //zfacesFluxStencil[zfaces[iface]] = fluxStencil;
          //zfacesInhTerm[zfaces[iface]] = inhomogTerm;
          //fluxStencil *= factor;
          a_inhoTerm += inhomogTerm;
          a_vofStencil += fluxStencil;
        }
      }
    }
#endif      

    //for some reason, EBAMRPoissonOp calculates beta weight without the EB flux
    double diagweight;
    bool foundVoF = a_vofStencil.getValueAtVoF(diagweight, a_vof);
    PR_assert(foundVoF); //if this fails, the stencil has no diagonal value
    a_vofStencil.setBetaWeight(diagweight);
  }

  //Gets WLS flux for the divergence stencil
  template <int order>
  void
  WLSStencil<order>::
  getWLSDivFCompStencil(LocalStencil<CELL, double>    & a_vofStencilx,
                        LocalStencil<CELL, double>    & a_vofStencily,
                        LocalStencil<CELL, double>    & a_vofStencilz,
                         const EBIndex<CELL>                   & a_vof,
                         const EBGraph                         & a_graph,
                         const VoluData                        & a_voludata,
                         const EBFaData                        & a_ebfadata,
                         const XFacData                        & a_xfacdata,
                         const YFacData                        & a_yfacdata,
                         const ZFacData                        & a_zfacdata,
                         const EBNormalData                    & a_ebnormxdata,
                         const EBNormalData                    & a_ebnormydata,
                         const EBNormalData                    & a_ebnormzdata,
                         const DivWLSFlux<XFACE, order>  & a_xfaceflux,
                         const DivWLSFlux<YFACE, order>  & a_yfaceflux,
                         const DivWLSFlux<ZFACE, order>  & a_zfaceflux,
                         const DivWLSFlux<BOUNDARY, order>  & a_ebflux,
                         string                             a_ebbcName,
                         string                             a_dombcName[2*DIM],
                         Real                               & a_inhoTerm,
                         int                                & a_neighType,
                        double                                  a_dx)
  {
    PR_TIME("getWLSDivFCompStencil");
    // a_vofStencil.clear();
    a_vofStencilx.clear();
    a_vofStencily.clear();
#if DIM==3
    a_vofStencilz.clear();
#endif
    a_inhoTerm = 0;
    Real volScaleFactor  = (D_TERM(a_dx, *a_dx, *a_dx));
    Real areaScaleFactor = (D_TERM(1.0,  *a_dx, *a_dx));

    for (SideIterator sit; sit.ok(); ++sit)
    {
      vector< EBIndex<XFACE> > xfaces = a_graph.getXFaces(a_vof, sit());
      int isign = sign(sit());
      for (int iface = 0; iface < xfaces.size(); iface++)
      {
        LocalStencil<CELL, double> fluxMStencil;
        LocalStencil<CELL, double> fluxPStencil;
        LocalStencil<CELL, double> fluxStencil;
        Real inhomogTerm=0.0; //for inhomogeneous bcs
        double factor = double(isign);
        std::map<EBIndex<XFACE>,LocalStencil<CELL, double> >::iterator it;
        it = xfacesFluxStencil.find(xfaces[iface]);
        if (it != xfacesFluxStencil.end())
        {
          fluxStencil = it->second;
          inhomogTerm = xfacesInhTerm[xfaces[iface]];
          fluxStencil *= factor;
          a_inhoTerm += inhomogTerm*factor;
          a_vofStencilx += fluxStencil;
        }
        else
        {
          EBIndex<CELL> voflo = xfaces[iface].getVoF(Side::Lo);
          EBIndex<CELL> vofhi = xfaces[iface].getVoF(Side::Hi);
          double   areaFrac;
          RealVect centroid;
          Normalizor<order>::normalizedFaceMoments(areaFrac, centroid, a_xfacdata, xfaces[iface], a_graph, a_dx);

          RealVect faceLoc;
          for (int idir = 0; idir < DIM; idir++)
          {
            faceLoc[idir] = a_dx * (0.5 + xfaces[iface].m_pt[idir]);
            if (idir == 0) faceLoc[idir] = a_dx * (xfaces[iface].m_pt[idir]);
          }

          IndMomSDMinOne faceMom;
          if (areaFrac < 1.0)
          {
            faceMom = a_xfacdata(xfaces[iface], 0);
            for(MomentIterator<DIM-1, order> iter; iter.ok(); ++iter)
            {
              const IndexTM<int, DIM-1>& p = iter();
              faceMom[iter()] = faceMom[iter()]/(areaScaleFactor*POW(a_dx,p));
            }
          }
          else
          {
            faceMom.setToRegular(a_dx);
            for(MomentIterator<DIM-1, order> iter; iter.ok(); ++iter)
            {
              const IndexTM<int, DIM-1>& p = iter();
              faceMom[iter()] = faceMom[iter()]/(areaScaleFactor*POW(a_dx,p));
            }
          }

          if (a_xfaceflux.m_verbose) pout()<<" stencil Xface "<<xfaces[iface].m_pt<<endl;

          if (!xfaces[iface].m_isBoundary){
            double inhomogTermM = 0;
            double inhomogTermP = 0;
            a_xfaceflux.getFaceCenteredStencil(fluxMStencil, inhomogTermM, voflo, xfaces[iface],
                                               a_graph,
                                               a_voludata,
                                               a_ebfadata,
                                               a_xfacdata,
                                               a_yfacdata,
                                               a_zfacdata,
                                               a_ebnormxdata,
                                               a_ebnormydata,
                                               a_ebnormzdata,
                                               a_ebbcName,
                                               a_dombcName,
                                               0,
                                               faceMom,
                                               a_neighType,
                                               factor,
                                               a_dx);

            a_xfaceflux.getFaceCenteredStencil(fluxPStencil, inhomogTermP, vofhi, xfaces[iface],
                                               a_graph,
                                               a_voludata,
                                               a_ebfadata,
                                               a_xfacdata,
                                               a_yfacdata,
                                               a_zfacdata,
                                               a_ebnormxdata,
                                               a_ebnormydata,
                                               a_ebnormzdata,
                                               a_ebbcName,
                                               a_dombcName,
                                               0,
                                               faceMom,
                                               a_neighType,
                                               factor,
                                               a_dx);

            // fluxStencil = getAdvFluxStencil(fluxMStencil, fluxPStencil, centroid, voflo, vofhi, v_faceLoc, 0, a_voludata, a_graph, a_dx);
            // Really just the face average
            fluxStencil = getLapFluxStencil(fluxMStencil, fluxPStencil, centroid, voflo, vofhi, 0, a_voludata, a_graph, a_dx);
            inhomogTerm = getLapinhomogTerm(inhomogTermM, inhomogTermP, centroid, voflo, vofhi, 0, a_voludata, a_graph, a_dx);
            //if (a_xfaceflux.m_verbose)
            //{
            //  pout()<<" xflux"<<endl;
            //  pout()<<" weight for M "<<distP/(distM+distP)<<endl;
            //  pout()<<" weight for P "<<distM/(distM+distP)<<endl;
            //}

          }else{
            a_xfaceflux.getFaceCenteredStencil(fluxStencil, inhomogTerm, a_vof, xfaces[iface],
                                               a_graph,
                                               a_voludata,
                                               a_ebfadata,
                                               a_xfacdata,
                                               a_yfacdata,
                                               a_zfacdata,
                                               a_ebnormxdata,
                                               a_ebnormydata,
                                               a_ebnormzdata,
                                               a_ebbcName,
                                               a_dombcName,
                                               0,
                                               faceMom,
                                               a_neighType,
                                               factor,
                                               a_dx);
          }
          if (a_xfaceflux.m_verbose)
          {
            pout()<<" inhomogTerm Xface "<<inhomogTerm<<endl;
            pout()<<" divFluxStencil "<<endl;
            fluxStencil.print();
          }
          xfacesFluxStencil[xfaces[iface]] = fluxStencil;
          xfacesInhTerm[xfaces[iface]] = inhomogTerm;
          double factor = double(isign);
          fluxStencil *= factor;
          a_inhoTerm += inhomogTerm*factor;
          a_vofStencilx += fluxStencil;

        }
      }
    }

    for (SideIterator sit; sit.ok(); ++sit)
    {
      vector< EBIndex<YFACE> > yfaces = a_graph.getYFaces(a_vof, sit());
      int isign = sign(sit());
      for (int iface = 0; iface < yfaces.size(); iface++)
      {
        LocalStencil<CELL, double> fluxMStencil;
        LocalStencil<CELL, double> fluxPStencil;
        LocalStencil<CELL, double> fluxStencil;
        Real inhomogTerm=0.0; //for inhomogeneous bcs
        double factor = double(isign);
        std::map<EBIndex<YFACE>,LocalStencil<CELL, double> >::iterator it;
        it = yfacesFluxStencil.find(yfaces[iface]);
        if (it != yfacesFluxStencil.end())
        {
          fluxStencil = it->second;
          inhomogTerm = yfacesInhTerm[yfaces[iface]];
          fluxStencil *= factor;
          a_inhoTerm += inhomogTerm*factor;
          a_vofStencily += fluxStencil;
        }
        else
        {
          EBIndex<CELL> voflo = yfaces[iface].getVoF(Side::Lo);
          EBIndex<CELL> vofhi = yfaces[iface].getVoF(Side::Hi);
          double   areaFrac;
          RealVect centroid;
          Normalizor<order>::normalizedFaceMoments(areaFrac, centroid, a_yfacdata, yfaces[iface], a_graph, a_dx);

          RealVect faceLoc;
          for (int idir = 0; idir < DIM; idir++)
          {
            faceLoc[idir] = a_dx * (0.5 + yfaces[iface].m_pt[idir]);
            if (idir == 1) faceLoc[idir] = a_dx * (yfaces[iface].m_pt[idir]);
          }
          //RealVect v_faceLoc = a_yfaceflux.m_exactVel->getVelocityAtPoint(faceLoc);

          IndMomSDMinOne faceMom;
          if (areaFrac < 1.0)
          {
            faceMom = a_yfacdata(yfaces[iface], 0);
            for(MomentIterator<DIM-1, order> iter; iter.ok(); ++iter)
            {
              const IndexTM<int, DIM-1>& p = iter();
              faceMom[iter()] = faceMom[iter()]/(areaScaleFactor*POW(a_dx,p));
            }
          }
          else
          {
            faceMom.setToRegular(a_dx);
            for(MomentIterator<DIM-1, order> iter; iter.ok(); ++iter)
            {
              const IndexTM<int, DIM-1>& p = iter();
              faceMom[iter()] = faceMom[iter()]/(areaScaleFactor*POW(a_dx,p));
            }
          }

          if (a_yfaceflux.m_verbose) pout()<<" stencilY face "<<yfaces[iface].m_pt<<endl;

          if (!yfaces[iface].m_isBoundary){
            double inhomogTermM = 0;
            double inhomogTermP = 0;
            a_yfaceflux.getFaceCenteredStencil(fluxMStencil, inhomogTermM, voflo, yfaces[iface],
                                               a_graph,
                                               a_voludata,
                                               a_ebfadata,
                                               a_xfacdata,
                                               a_yfacdata,
                                               a_zfacdata,
                                               a_ebnormxdata,
                                               a_ebnormydata,
                                               a_ebnormzdata,
                                               a_ebbcName,
                                               a_dombcName,
                                               1,
                                               faceMom,
                                               a_neighType,
                                               factor,
                                               a_dx);

            a_yfaceflux.getFaceCenteredStencil(fluxPStencil, inhomogTermP, vofhi, yfaces[iface],
                                               a_graph,
                                               a_voludata,
                                               a_ebfadata,
                                               a_xfacdata,
                                               a_yfacdata,
                                               a_zfacdata,
                                               a_ebnormxdata,
                                               a_ebnormydata,
                                               a_ebnormzdata,
                                               a_ebbcName,
                                               a_dombcName,
                                               1,
                                               faceMom,
                                               a_neighType,
                                               factor,
                                               a_dx);
            // fluxStencil = getAdvFluxStencil(fluxMStencil, fluxPStencil, centroid, voflo, vofhi, v_faceLoc, 1, a_voludata, a_graph, a_dx);
                        // Really just the face average
            fluxStencil = getLapFluxStencil(fluxMStencil, fluxPStencil, centroid, voflo, vofhi, 1, a_voludata, a_graph, a_dx);
            inhomogTerm = getLapinhomogTerm(inhomogTermM, inhomogTermP, centroid, voflo, vofhi, 1, a_voludata, a_graph, a_dx);
            //if (a_yfaceflux.m_verbose)
            //{
            //  pout()<<" yflux"<<endl;
            //  pout()<<" weight for M "<<distP/(distM+distP)<<endl;
            //  pout()<<" weight for P "<<distM/(distM+distP)<<endl;
            //}
          }
          else
          {
            a_yfaceflux.getFaceCenteredStencil(fluxStencil, inhomogTerm, a_vof, yfaces[iface],
                                               a_graph,
                                               a_voludata,
                                               a_ebfadata,
                                               a_xfacdata,
                                               a_yfacdata,
                                               a_zfacdata,
                                               a_ebnormxdata,
                                               a_ebnormydata,
                                               a_ebnormzdata,
                                               a_ebbcName,
                                               a_dombcName,
                                               1,
                                               faceMom,
                                               a_neighType,
                                               factor,
                                               a_dx);
          }
          if (a_yfaceflux.m_verbose)
          {
            pout()<<" inhomogTerm Yface "<<inhomogTerm<<endl;
            pout()<<" divFluxStencil "<<endl;
            fluxStencil.print();
          }
          yfacesFluxStencil[yfaces[iface]] = fluxStencil;
          yfacesInhTerm[yfaces[iface]] = inhomogTerm;
          double factor = double(isign);
          fluxStencil *= factor;
          a_inhoTerm += inhomogTerm*factor;
          a_vofStencily += fluxStencil;
        }
      }
    }
#if DIM==3
    for (SideIterator sit; sit.ok(); ++sit)
    {
      vector< EBIndex<ZFACE> > zfaces = a_graph.getZFaces(a_vof, sit());
      int isign = sign(sit());
      for (int iface = 0; iface < zfaces.size(); iface++)
      {
        LocalStencil<CELL, double> fluxMStencil;
        LocalStencil<CELL, double> fluxPStencil;
        LocalStencil<CELL, double> fluxStencil;
        Real inhomogTerm=0.0; //for inhomogeneous bcs
        double factor = double(isign);
        std::map<EBIndex<ZFACE>,LocalStencil<CELL, double> >::iterator it;
        it = zfacesFluxStencil.find(zfaces[iface]);
        if (it != zfacesFluxStencil.end())
        {
          fluxStencil = it->second;
          inhomogTerm = zfacesInhTerm[zfaces[iface]];
          fluxStencil *= factor;
          a_inhoTerm += inhomogTerm*factor;
          a_vofStencilz += fluxStencil;
        }
        else
        {
          EBIndex<CELL> voflo = zfaces[iface].getVoF(Side::Lo);
          EBIndex<CELL> vofhi = zfaces[iface].getVoF(Side::Hi);
          double   areaFrac;
          RealVect centroid;
          Normalizor<order>::normalizedFaceMoments(areaFrac, centroid, a_zfacdata, zfaces[iface], a_graph, a_dx);
          RealVect faceLoc;
          for (int idir = 0; idir < DIM; idir++)
          {
            faceLoc[idir] = a_dx * (0.5 + zfaces[iface].m_pt[idir]);
            if (idir == 2) faceLoc[idir] = a_dx * (zfaces[iface].m_pt[idir]);
          }
          //RealVect v_faceLoc = a_zfaceflux.m_exactVel->getVelocityAtPoint(faceLoc);

          IndMomSDMinOne faceMom;
          if (areaFrac < 1.0)
          {
            faceMom = a_zfacdata(zfaces[iface], 0);
            for(MomentIterator<DIM-1, order> iter; iter.ok(); ++iter)
            {
              const IndexTM<int, DIM-1>& p = iter();
              faceMom[iter()] = faceMom[iter()]/(areaScaleFactor*POW(a_dx,p));
            }
          }
          else
          {
            faceMom.setToRegular(a_dx);
            for(MomentIterator<DIM-1, order> iter; iter.ok(); ++iter)
            {
              const IndexTM<int, DIM-1>& p = iter();
              faceMom[iter()] = faceMom[iter()]/(areaScaleFactor*POW(a_dx,p));
            }
          }

          if (!zfaces[iface].m_isBoundary){
            double inhomogTermM = 0;
            double inhomogTermP = 0;
            a_zfaceflux.getFaceCenteredStencil(fluxMStencil, inhomogTermM, voflo, zfaces[iface],
                                               a_graph,
                                               a_voludata,
                                               a_ebfadata,
                                               a_xfacdata,
                                               a_yfacdata,
                                               a_zfacdata,
                                               a_ebnormxdata,
                                               a_ebnormydata,
                                               a_ebnormzdata,
                                               a_ebbcName,
                                               a_dombcName,
                                               2,
                                               faceMom,
                                               a_neighType,
                                               factor,
                                               a_dx);

            a_zfaceflux.getFaceCenteredStencil(fluxPStencil, inhomogTermP, vofhi, zfaces[iface],
                                               a_graph,
                                               a_voludata,
                                               a_ebfadata,
                                               a_xfacdata,
                                               a_yfacdata,
                                               a_zfacdata,
                                               a_ebnormxdata,
                                               a_ebnormydata,
                                               a_ebnormzdata,
                                               a_ebbcName,
                                               a_dombcName,
                                               2,
                                               faceMom,
                                               a_neighType,
                                               factor,
                                               a_dx);
            // fluxStencil = getAdvFluxStencil(fluxMStencil, fluxPStencil, centroid, voflo, vofhi, v_faceLoc, 2, a_voludata, a_graph, a_dx);
            // Really just the face average
            fluxStencil = getLapFluxStencil(fluxMStencil, fluxPStencil, centroid, voflo, vofhi, 2, a_voludata, a_graph, a_dx);
            inhomogTerm = getLapinhomogTerm(inhomogTermM, inhomogTermP, centroid, voflo, vofhi, 2, a_voludata, a_graph, a_dx);
          }else{
            a_zfaceflux.getFaceCenteredStencil(fluxStencil, inhomogTerm, a_vof, zfaces[iface],
                                               a_graph,
                                               a_voludata,
                                               a_ebfadata,
                                               a_xfacdata,
                                               a_yfacdata,
                                               a_zfacdata,
                                               a_ebnormxdata,
                                               a_ebnormydata,
                                               a_ebnormzdata,
                                               a_ebbcName,
                                               a_dombcName,
                                               2,
                                               faceMom,
                                               a_neighType,
                                               factor,
                                               a_dx);

          }
          zfacesFluxStencil[zfaces[iface]] = fluxStencil;
          zfacesInhTerm[zfaces[iface]] = inhomogTerm;
          double factor = double(isign);
          fluxStencil *= factor;
          a_inhoTerm += inhomogTerm*factor;
          a_vofStencilz += fluxStencil;
        }
      }
    }
#endif

    //for some reason, EBAMRPoissonOp calculates beta weight without the EB flux
      {
        double diagweight;
        // a_vofStencilx.print();
        bool foundVoF = a_vofStencilx.getValueAtVoF(diagweight, a_vof);
        PR_assert(foundVoF); //if this fails, the stencil has no diagonal value
        a_vofStencilx.setBetaWeight(diagweight);

        foundVoF = a_vofStencily.getValueAtVoF(diagweight, a_vof);
        PR_assert(foundVoF); //if this fails, the stencil has no diagonal value
        a_vofStencilx.setBetaWeight(diagweight);
#if DIM==3
        foundVoF = a_vofStencilz.getValueAtVoF(diagweight, a_vof);
        PR_assert(foundVoF); //if this fails, the stencil has no diagonal value
        a_vofStencilz.setBetaWeight(diagweight);
#endif
      }
  }

  //Gets WLS the cell gradient stencil
  template <int order>
  void
  WLSStencil<order>::
  getWLSGradStencil(LocalStencil<CELL, double>    & a_vofStencilx,
                    LocalStencil<CELL, double>    & a_vofStencily,
                    LocalStencil<CELL, double>    & a_vofStencilz,
                    const EBIndex<CELL>                   & a_vof,
                    const EBGraph                         & a_graph,
                    const VoluData                        & a_voludata,
                    const EBFaData                        & a_ebfadata,
                    const XFacData                        & a_xfacdata,
                    const YFacData                        & a_yfacdata,
                    const ZFacData                        & a_zfacdata,
                    const EBNormalData                    & a_ebnormxdata,
                    const EBNormalData                    & a_ebnormydata,
                    const EBNormalData                    & a_ebnormzdata,
                    const BaseWLSCell<order>              & a_cellConstruct,
                    string                             a_ebbcName,
                    string                             a_dombcName[2*DIM],
                    vector<Real>                               & a_inhoTerm,
                    int                                & a_neighType,
                    double                                  a_dx)
  {
    PR_TIME("getWLSGradStencil");
    // a_vofStencil.clear();
    a_vofStencilx.clear();
    a_vofStencily.clear();
#if DIM==3
    a_vofStencilz.clear();
#endif
    CH_assert(a_inhoTerm.size() >= DIM);
    for (int dir=0; dir!=DIM; dir++)
      {
        a_inhoTerm[dir] = 0;
      }

    Real volScaleFactor  = (D_TERM(a_dx, *a_dx, *a_dx));
    Real areaScaleFactor = (D_TERM(1.0,  *a_dx, *a_dx));

    double   voluFrac;
    RealVect centroid;
    Normalizor<order>::normalizedVoluMoments(voluFrac, centroid, a_voludata, a_vof, a_graph, a_dx);

    // x direction
    {
      IndMomDIM voluMom;
      if (voluFrac < 1.0)
        {
          voluMom = a_voludata(a_vof, 0);
          for(MomentIterator<DIM, order> iter; iter.ok(); ++iter)
            {
              const auto& p = iter();
              voluMom[iter()] = voluMom[iter()]/(volScaleFactor*POW(a_dx, p));
            }
        }
      else
        {
          voluMom.setToRegular(a_dx);
          for(MomentIterator<DIM, order> iter; iter.ok(); ++iter)
            {
              const auto& p = iter();
              voluMom[iter()] = voluMom[iter()]/(volScaleFactor*POW(a_dx, p));
            }
        }

      if (a_cellConstruct.m_verbose) pout() << " stencil xdir "<< a_vof.m_pt << endl;

      IndexTM<int, DIM> deriv = IndexTM<int, DIM>::Zero;
      deriv[0] = 1; // x-derivative stencil
      a_cellConstruct.getCellStencil(a_vofStencilx,
                                     a_inhoTerm[0],
                                     a_vof,
                                     a_graph,
                                     a_voludata,
                                     a_ebfadata,
                                     a_xfacdata,
                                     a_yfacdata,
                                     a_zfacdata,
                                     a_ebnormxdata,
                                     a_ebnormydata,
                                     a_ebnormzdata,
                                     a_ebbcName,
                                     a_dombcName,
                                     voluMom,
                                     a_dx,
                                     deriv);
      if (a_cellConstruct.m_verbose)
        {
          pout() << " inhomogTerm xdir " << a_inhoTerm[0] << endl;
          pout() << " stencil " << endl;
          a_vofStencilx.print();
        }
    }

    {
      IndMomDIM voluMom;
      if (voluFrac < 1.0)
        {
          voluMom = a_voludata(a_vof, 0);
          for(MomentIterator<DIM, order> iter; iter.ok(); ++iter)
            {
              const auto& p = iter();
              voluMom[iter()] = voluMom[iter()]/(volScaleFactor*POW(a_dx, p));
            }
        }
      else
        {
          voluMom.setToRegular(a_dx);
          for(MomentIterator<DIM, order> iter; iter.ok(); ++iter)
            {
              const auto& p = iter();
              voluMom[iter()] = voluMom[iter()]/(volScaleFactor*POW(a_dx, p));
            }
        }

      if (a_cellConstruct.m_verbose) pout() << " stencil ydir "<< a_vof.m_pt << endl;

      IndexTM<int, DIM> deriv = IndexTM<int, DIM>::Zero;
      deriv[1] = 1; // x-derivative stencil
      a_cellConstruct.getCellStencil(a_vofStencily,
                                     a_inhoTerm[1],
                                     a_vof,
                                     a_graph,
                                     a_voludata,
                                     a_ebfadata,
                                     a_xfacdata,
                                     a_yfacdata,
                                     a_zfacdata,
                                     a_ebnormxdata,
                                     a_ebnormydata,
                                     a_ebnormzdata,
                                     a_ebbcName,
                                     a_dombcName,
                                     voluMom,
                                     a_dx,
                                     deriv);
      if (a_cellConstruct.m_verbose)
        {
          pout() << " inhomogTerm ydir " << a_inhoTerm[1] << endl;
          pout() << " stencil " << endl;
          a_vofStencily.print();
        }
    }

    //for some reason, EBAMRPoissonOp calculates beta weight without the EB flux
      {
        double diagweight;
        // a_vofStencilx.print();
        bool foundVoF = a_vofStencilx.getValueAtVoF(diagweight, a_vof);
        PR_assert(foundVoF); //if this fails, the stencil has no diagonal value
        a_vofStencilx.setBetaWeight(diagweight);

        foundVoF = a_vofStencily.getValueAtVoF(diagweight, a_vof);
        PR_assert(foundVoF); //if this fails, the stencil has no diagonal value
        a_vofStencilx.setBetaWeight(diagweight);
#if DIM==3
        foundVoF = a_vofStencilz.getValueAtVoF(diagweight, a_vof);
        PR_assert(foundVoF); //if this fails, the stencil has no diagonal value
        a_vofStencilz.setBetaWeight(diagweight);
#endif
      }
  }

  //Compute Dirichlet BCs at faces for the Moment Matrix.
  template < int order>
  void
  DirichletBC<order>::
  getFaceMoments(IndMomDIM      & a_fac2volmom,
                 int              a_faceDir,
                 IndMomSDMinOne & a_facmom)
  {
    PR_TIME("getFaceMomentsDirichletBC");
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

  //Compute Neumann BCs at faces for the Moment Matrix.
  template < int order>
  void
  NeumannBC<order>::
  getFaceMoments(IndMomDIM      & a_fac2volmom,
                 int              a_faceDir,
                 IndMomSDMinOne & a_facmom)
  {
    PR_TIME("getFaceMomentsNeumannBC");
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

  //Compute Neumann BCs at EB faces for the Moment Matrix.
  template < int order>
  void
  NeumannBC<order>::
  getEBMoments(IndMomDIM         & a_ebfmom,
               Vector<IndMomDIM> & a_ebfNorm)
  {
    PR_TIME("getEBMomentsNeumannBC");
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

  //Compute Neumann BCs at faces and EB Faces. Required for inhomogeneous conditions.
  //getIntegralNormalDerivative
  template < int order>
  Real
  NeumannExactSoln<order>::
  getFaceIntegral(const IndMomSDMinOne  & a_faceMom,
                  const RealVect        & a_faceLoc,
                  int                     a_faceDir,
                  const shared_ptr<BaseExactSolution<order> >  & m_exactSol,
                  const Real              a_dx)  
  {
    PR_TIME("getFaceIntegralNeumannExactSoln");
    Real retval = 0.0;
    Real areaScaleFactor = (D_TERM(1.0,  *a_dx, *a_dx));
    MomentIterator<DIM-1, order> momit;
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
 
      // pSpaceDim[a_face.direction()] is 1 so no need to multiply
      retval += moment*(areaScaleFactor*POW(a_dx,pSpaceDim))*m_exactSol->getDerivative(pSpaceDim,a_faceLoc)/(pSpaceDim.factorial());
    }
  
    return retval;
  }

  //Compute Dirichlet BCs at faces and EB Faces. Required for inhomogeneous conditions.
  //getIntegralPhi
  template < int order>
  Real
  DirichletExactSoln<order>::
  getFaceIntegral(const IndMomSDMinOne  & a_faceMom,
                  const RealVect        & a_faceLoc,
                  int                     a_faceDir,
                  const shared_ptr<BaseExactSolution<order> >  & m_exactSol,
                  const Real              a_dx) 
  {
    PR_TIME("getFaceIntegralDirichletExactSoln");
    Real retval = 0.0;
    Real areaScaleFactor = (D_TERM(1.0,  *a_dx, *a_dx));
    MomentIterator<DIM-1, order> momit;
    for(momit.reset(); momit.ok(); ++momit)
    {
      IndexTM<int,DIM-1> pSDMinOne = momit();
      IndexTM<int,DIM> pSpaceDim;

      // face moment m_p is zero except when p[a_face.direction()] = 0 
      int index = 0;
      for(int idir = 0; idir < DIM; idir++)
      {
        if(idir == a_faceDir)
        {
          pSpaceDim[idir] = 0;
        }
        else
        {
          pSpaceDim[idir] = pSDMinOne[index];
          index++;
        }
      }

      Real moment = a_faceMom[pSDMinOne];

      retval += moment*(areaScaleFactor*POW(a_dx,pSpaceDim))*m_exactSol->getDerivative(pSpaceDim,a_faceLoc)/(pSpaceDim.factorial());

    }

    return retval;
  }

  //Compute Dirichlet BCs at faces and EB Faces. Required for inhomogeneous conditions.
  //getIntegralPhiAlongEB
  template < int order>
  Real
  DirichletExactSoln<order>::
  getEBIntegral(const IndMomDIM  & a_ebfMom,
                const RealVect   & a_ebfLoc,
                const shared_ptr<BaseExactSolution<order> >  & m_exactSol,
                const Real         a_dx) 
  {
    PR_TIME("getEBIntegralDirichletExactSoln");

    Real retval = 0.0;
    Real areaScaleFactor = (D_TERM(1.0,  *a_dx, *a_dx));

    MomentIterator<DIM, order> momit;
    for(momit.reset(); momit.ok(); ++momit)
    {
      IndexTM<int,DIM> p = momit();
      Real moment = a_ebfMom[p];

      retval += moment*(areaScaleFactor*POW(a_dx,p))*m_exactSol->getDerivative(p,a_ebfLoc)/(p.factorial());

    }

    return retval;
  }

  //Compute Neumann BCs at faces and EB Faces. Required for inhomogeneous conditions.
  //getIntegralGradPhiDotNormal
  template < int order>
  Real
  NeumannExactSoln<order>::
  getEBIntegral(const Vector<IndMomDIM>  & a_ebfNormalMom,
                const RealVect           & a_ebfLoc,
                const shared_ptr<BaseExactSolution<order> >  & m_exactSol,
                const Real                 a_dx) 
  {
    PR_TIME("getEBIntegralNeumannExactSoln");

    Real retval = 0.0;
    Real areaScaleFactor = (D_TERM(1.0,  *a_dx, *a_dx));
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
  
          retval += pd*moment*(areaScaleFactor*POW(a_dx,p))*m_exactSol->getDerivative(p,a_ebfLoc)/(p.factorial());
        }
      }
    }

    return retval;

  }

  //Compute rhs matrix for a face
  template<CENTERING cent, int order>
  void
  ExactVelAdvectionWLSFlux<cent, order>::
  getFaceFluxMatrix(Chombo4::LAPACKMatrix    & a_fluxMatrix,
                    const EBIndex<cent>      & a_face,
                    RealVect                 & a_faceLoc,
                    IndMomSDMinOne           & a_faceMom,
                    double                     a_signNormal,
                    double                     a_dx) const
  {
    PR_TIME("getFaceFluxMatrix");
    Real areaScaleFactor  = (D_TERM(1.0,  *a_dx, *a_dx));
    int a_faceDir = a_face.direction();

    MomentIterator<DIM, order> momit;
    IndMomDIM facemomDIM;

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
        diffrv[index] = diff[idir]/a_dx;
        index++;
      }
    }
    a_faceMom.shift(diffrv);

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
        facemomDIM[momit()] = a_faceMom[areaMomIndex];
      }
      else
      {
        facemomDIM[momit()] = 0.0;
      }
    }

    m_exactVel->setDx(a_dx);
    IndMomDIM velCoeff = m_exactVel->getVelCompTaylorCoeffs(a_faceLoc, a_face.direction());
    a_fluxMatrix.define(1, momit.size());
    a_fluxMatrix.setVal(0.);
    const int& faceDir = a_face.direction();
    if (this->m_verbose) pout() << " RHS - Adv Face Flux Matrix " << endl;

    for(momit.reset(); momit.ok(); ++momit)
    {   
      IndexTM<int, DIM> p = momit();
      int icol = IndexedMoments<DIM, order>::indexOf(p);
      a_fluxMatrix(0, icol) = 0.0;
   
      MomentIterator<DIM, order> velMomit;
      for(velMomit.reset(); velMomit.ok(); ++velMomit)
      {
        IndexTM<int, DIM> q = velMomit();
        IndexTM<int, DIM> pPlusQ = p+q;
 
        if(pPlusQ.sum() <= order)
        {
          a_fluxMatrix(0,icol) += velCoeff[q]*facemomDIM[pPlusQ]*a_signNormal;
        }
      }
      if(Chombo4::LAPACKMatrix::s_verbose)
      {
        pout() << ",icol = "  << icol <<  ", Q(0,icol) = " << a_fluxMatrix(0, icol) << " a_signNormal "<<a_signNormal<<endl;
      }
    }
  }
  /***/

  //Compute rhs matrix for a EB face
  template<CENTERING cent, int order>
  void
  ExactVelAdvectionWLSFlux<cent, order>::
  getEBFluxMatrix(Chombo4::LAPACKMatrix    & a_fluxMatrix,
                  const EBIndex<cent>      & a_face,
                  RealVect                 & a_faceLoc,
                  Vector<IndMomDIM>        & a_ebnorm,
                  double                     a_dx) const
  {
    PR_TIME("getEBFluxMatrix");
    Chombo4::MayDay::Error("should never get here  (eb flux always = 0)"); 
  }

  //Compute weight matrix
  template<CENTERING cent, int order>
  void
  ExactVelAdvectionWLSFlux<cent, order>::
  getWeightsMatrixFace(Chombo4::LAPACKMatrix  & a_weights,
                       const Neighbors<order> & a_neigh,
                       RealVect               & a_faceLoc) const
  { 
    PR_TIME("getWeightsMatrixFace");
    Vector<int> weightDisVol = a_neigh.m_weightDisVol;
    Vector<int> weightDisFac = a_neigh.m_weightDisFac;
    Vector<int> weightDisEB = a_neigh.m_weightDisEB;
    Vector<RealVect> volneighLoc = a_neigh.m_volLoc;
    Vector<RealVect> ebfneighLoc = a_neigh.m_ebfLoc;
    Vector<IndMomDIM> volMom = a_neigh.m_volMom; 
    Vector<IndMomDIM> ebfMom = a_neigh.m_ebfMom; 
    RealVect v_faceLoc = m_exactVel->getVelocityAtPoint(a_faceLoc);
    Real dx = m_exactVel->getDx();

    int nrows = weightDisVol.size()+weightDisFac.size()+weightDisEB.size();
    a_weights.define(nrows, nrows);
    a_weights.setVal(1.0);
    a_weights.setToIdentity(); //this sets all the off diagonals to zero

    int max_rad = 0.0;
    for(int ineigh = 0; ineigh < weightDisVol.size(); ineigh++)
    {
      if (max_rad<weightDisVol[ineigh])
      {
        max_rad = weightDisVol[ineigh];
      }
    }

    Real coeff = 0.5;
    RealVect x_upw = (a_faceLoc/dx) - coeff*(v_faceLoc);

    int powerE = this->m_powerE;;
    //if (max_rad-1 == this->m_radius) powerE = this->m_powerE;
    //else powerE = this->m_powerE - (2.0*(max_rad-this->m_radius-1));
    if (this->m_verbose)  pout() << " Weights Matrix Adv " << endl;
    for(int ineigh = 0; ineigh < nrows; ineigh++)
    {
      Real blerg;
      if (ineigh<weightDisVol.size())
      {
        Real centroid_x = volMom[ineigh][1]/(volMom[ineigh][IndexTM<int,DIM>::Zero]);
        Real centroid_y = volMom[ineigh][1+order]/(volMom[ineigh][IndexTM<int,DIM>::Zero]);
#if DIM==3
        Real centroid_z = volMom[ineigh][2*(1+order)]/(volMom[ineigh][IndexTM<int,DIM>::Zero]);
#endif
        RealVect volLoc = volneighLoc[ineigh];
        volLoc[0] += centroid_x*dx;
        volLoc[1] += centroid_y*dx;
#if DIM==3
        volLoc[2] += centroid_z*dx;
#endif
        RealVect weight = (volLoc/dx) - x_upw;
#if DIM==2
        Real weight_mag = pow(pow(weight[0],2)+pow(weight[1],2),0.5);
#elif DIM==3
        Real weight_mag = pow(pow(weight[0],2)+pow(weight[1],2)+pow(weight[2],2),0.5);
#endif
        Real max_weight = std::max(2*weight_mag,1.0);
        blerg = pow(max_weight,-order-1);
      }
      else if (ineigh>=weightDisVol.size() && ineigh<weightDisVol.size()+weightDisFac.size())
      {
        int dist = weightDisFac[ineigh-weightDisVol.size()];
        //blerg = POW(dist,powerE);
      }
      if (ineigh>=weightDisVol.size()+weightDisFac.size())
      {
        int dist = weightDisEB[ineigh-weightDisVol.size()-weightDisFac.size()];
        //blerg = POW(dist,powerE);
      }
      a_weights(ineigh, ineigh) = blerg;
      if(Chombo4::LAPACKMatrix::s_verbose)
      {
        pout() << ",icol,irow = "  << ineigh <<  ", W(irow,icol) = " << a_weights(ineigh, ineigh) << endl;
      }
    }
  }

  template<int order>
  RealVect
  ConstantVectorVel<order>::
  getVelocityAtPoint(const RealVect& a_x) const
  { 
    return m_value;
  } 
      
  /***/
  template<int order>
  IndexedMoments<DIM  , order>
  ConstantVectorVel<order>::
  getVelCompTaylorCoeffs(const RealVect   & a_x0,
                         const int        & a_velcomp) const
  {
    IndMomDIM retval;
    for(MomentIterator<DIM, order> momit; momit.ok(); ++momit)
      {
        const IndexTM<int, DIM>& p = momit();
        int orderP = p.sum();
        if(orderP == 0)
          {
            retval[p] = m_value[a_velcomp];
          }
        else
          {
            retval[p] =  0;
          }
        retval[p] = retval[p]*(POW(this->m_dx,p));
      }
  
    //taylor coeff so divide by p! (no need because all derivs zero)
    //retval.divideByFactorial();
    return retval;
  }

/***/
template<int order>
RealVect
SolidBodyRotVel<order>::
getVelocityAtPoint(const RealVect& a_x) const
{
  RealVect rvec = a_x - m_center;
  RealVect retval;
  retval[0] = -m_omega*(rvec[1]);
  retval[1] = m_omega*(rvec[0]);
  for(int idir = 2; idir < DIM; idir++)
    {
      retval[idir] = 0.0;
    }

  return retval;
}

/***/
template<int order>
IndexedMoments<DIM  , order>
SolidBodyRotVel<order>::
getVelCompTaylorCoeffs(const RealVect   & a_x0,
                       const int        & a_velcomp) const
{
  IndMomDIM retval;
  RealVect rvec = a_x0 - m_center;
  for(MomentIterator<DIM, order> momit; momit.ok(); ++momit)
    {
      const IndexTM<int, DIM>& p = momit();
      int orderP = p.sum();
      //vel = C*(-y, x)
      if(orderP == 0)
        {
          if(a_velcomp == 0)
            {
              retval[p] = -m_omega*(rvec[1]);
            }
          else if(a_velcomp == 1)
            {
              retval[p] =  m_omega*(rvec[0]);
            }
          else
            {
              retval[p] =  0;
            }
        }
      else if (orderP == 1)
        {
          if((a_velcomp == 0) &&  (p[1]== 1))
            {
              retval[p] = -m_omega;
            }
          else if((a_velcomp == 1) &&  (p[0]== 1))
            {
              retval[p] =  m_omega;
            }
          else
            {
              retval[p] = 0;
            }
        }
      else //no higher orders 
        {
          retval[p] =  0;
        }
        retval[p] = retval[p]*(pow(this->m_dx,orderP));
    }

  //taylor coeff so divide by p! (no need because only first derivs non zero)
  //retval.divideByFactorial();
  return retval;
}


/***/
template<int order>
RealVect
PoiseuilleVectorVel<order>::
getVelocityAtPoint(const RealVect& a_x) const
{
  RealVect rvec = a_x - m_center;
  Real rvec_mag = pow(rvec[0],2);
  RealVect retval;
  retval[0] = 0;
  retval[1] = m_wallValue + (m_centerValue*(1.0-(rvec_mag/pow(m_radius,2.0))));

  return retval;
}

/***/
template<int order>
IndexedMoments<DIM  , order>
PoiseuilleVectorVel<order>::
getVelCompTaylorCoeffs(const RealVect   & a_x0,
                       const int        & a_velcomp) const
{
  IndMomDIM retval;
  RealVect rvec = a_x0 - m_center;
  Real rvec_mag = abs(rvec[0]);
  for(MomentIterator<DIM, order> momit; momit.ok(); ++momit)
    {
      const IndexTM<int, DIM>& p = momit();
      int orderP = p.sum();
      if(orderP == 0)
        {
          if(a_velcomp == 0)
            {
              retval[p] = 0;
            }
          else if(a_velcomp == 1)
            {
              retval[p] = m_wallValue + (m_centerValue*(1.0-pow(rvec_mag/m_radius,2.0)));
            }
          else
            {
              retval[p] =  0;
            }
        }
      else if (orderP == 1)
        {
          if(a_velcomp == 0)
            {
              retval[p] = 0;
            }
          else if((a_velcomp == 1) &&  (p[0]== 1))
            {
              retval[p] =  -2*m_centerValue*(rvec[0])/pow(m_radius,2.0);
            }
          else if((a_velcomp == 1) &&  (p[1]== 1))
            {
              retval[p] =  0;
            }
          else
            {
              retval[p] = 0;
            }
        }
      else if (orderP == 2)
        {
          if(a_velcomp == 0)
            {
              retval[p] = 0;
            }
          else if((a_velcomp == 1) &&  (p[0]== 2))
            {
              retval[p] =  -2*m_centerValue/pow(m_radius,2.0);
            }
          else if((a_velcomp == 1) &&  (p[1]== 2))
            {
              retval[p] =  0;
            }
          else
            {
              retval[p] = 0;
            }
        }
      else //no higher orders 
        {
          retval[p] =  0;
        }
        retval[p] = retval[p]*(pow(this->m_dx,orderP))/p.factorial();
    }

  //taylor coeff so divide by p! (no need because only first derivs non zero)
  //retval.divideByFactorial();
  return retval;
}



}
#endif
