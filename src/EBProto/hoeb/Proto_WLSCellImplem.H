#ifndef _Proto_WLSCellImplem_H_
#define _Proto_WLSCellImplem_H_

namespace Proto
{
  //Compute the Moments Matrix
  template<CENTERING cent, int order>
  void
  BaseWLSFlux<cent, order>::
  getShiftedMomentMatrix(Chombo4::LAPACKMatrix    & a_Mvol,
                         RealVect               & a_faceLoc,
                         Real                     a_dx,
                         const Neighbors<order> & a_neigh,
                         string                   a_ebbcName,
                         string                   a_dombcName[2*DIM]) const
  {
    PR_TIME("getShiftedMomentMatrix");
    using Chombo4::pout;

    int nvol = a_neigh.m_volMom.size();
    int nfac = a_neigh.m_facMom.size();
    int nebf = a_neigh.m_ebfMom.size();

    int nrows = nvol + nfac + nebf;
    MomentIterator<DIM, order> momit;
    a_Mvol.define(nrows, momit.size());

    if (m_verbose) 
    {
      pout() << "Size Moments matrix "<< nrows << " " << momit.size() << std::endl;
      pout() << "With nvol + nfac + nebf " << nvol << " " << nfac << " " << nebf << std::endl;
      pout() << "Moments matrix" << endl;
    }

    // Start of the volume block matrix (row, column)
    IndexTM<int, 2> start = IndexTM<int, 2>{0, 0};
    // fill in the sub-matrix
    this->getShiftedMomentMatrixVolumes(a_Mvol,
                                        a_faceLoc,
                                        a_dx,
                                        a_neigh,
                                        start);

    // Start of the face block matrix
    start = IndexTM<int, 2>{nvol, 0};
    // fill in the sub-matrix
    this->getShiftedMomentMatrixRegFaces(a_Mvol,
                                         a_faceLoc,
                                         a_dx,
                                         a_neigh,
                                         a_dombcName,
                                         start);

    // Start of the EB block matrix
    start = IndexTM<int, 2>{nvol+nfac, 0};
    // fill in the sub-matrix
    this->getShiftedMomentMatrixEBFaces(a_Mvol,
                                        a_faceLoc,
                                        a_dx,
                                        a_neigh,
                                        a_ebbcName,
                                        start);

    if (m_verbose)
      {
        pout() << "Moment Matrix" << std::endl;
        a_Mvol.poutMatlab();
      }
  }

      //Compute the moments matrix portion for cells in the neighborhood
  template<CENTERING cent, int order>
  void
  BaseWLSFlux<cent, order>::
  getShiftedMomentMatrixVolumes(Chombo4::LAPACKMatrix  & a_Mblock,
                                RealVect                 a_faceLoc,
                                Real                     a_dx,
                                const Neighbors<order> & a_neigh,
                                IndexTM<int, 2>         a_start) const
  {
    PR_TIME("getShiftedMomentMatrixVolumes");
    using Chombo4::pout;
    Vector<IndMomDIM> volmom = a_neigh.m_volMom;
    Vector<RealVect> volneighLoc = a_neigh.m_volLoc;
    const Vector<Real> volneighFrac = a_neigh.m_volneighFrac;

    int nvol = volmom.size();
    MomentIterator<DIM, order> momit;

    //Adding volume moments
    for(int irow = 0; irow < nvol; irow++)
      {
        //the moments start at the cell CENTER
        shiftMomentToLoc(volmom[irow], volneighLoc[irow], a_faceLoc, a_dx);

        if(Chombo4::LAPACKMatrix::s_verbose)
          {
            pout() << "irow = " << irow << endl;
          }
        int mrow = irow + a_start[0];
        for(momit.reset(); momit.ok(); ++momit)
          {
            //need to use the lower order moment declaration to
            //get the correct column number
            int icol = IndexedMoments<DIM, order>::indexOf(momit());
            int mcol = icol + a_start[1];
            if (volneighFrac[irow] < m_eps)
              {
                Chombo4::MayDay::Warning("Zeroing all but the first element of a volume moment row");
                if (icol == 0)
                  {
                    a_Mblock(mrow, mcol) = 1.0;
                  }
                else
                  {
                    a_Mblock(mrow, mcol) = 0.0;
                  }
              }
            else
              {
                a_Mblock(mrow, mcol) = volmom[irow][momit()]/(volneighFrac[irow]);
              }
            if(Chombo4::LAPACKMatrix::s_verbose)
              {
                pout() << ",icol = "  << icol << ", momind = " << momit() << ", mat(irow,icol) = " << a_Mblock(irow, icol) << endl;
              }
          }
      }
  }

  //Compute the moments matrix portion for regular faces in the neighborhood
  template<CENTERING cent, int order>
  void
  BaseWLSFlux<cent, order>::
  getShiftedMomentMatrixRegFaces(Chombo4::LAPACKMatrix  & a_Mblock,
                                 RealVect                 a_faceLoc,
                                 Real                     a_dx,
                                 const Neighbors<order> & a_neigh,
                                 string                   a_dombcName[2*DIM],
                                 IndexTM<int, 2>         a_start) const
  {
    PR_TIME("getShiftedMomentMatrixFaces");
    using Chombo4::pout;
    Vector<IndMomSDMinOne> facmom = a_neigh.m_facMom;
    Vector<RealVect> facneighLoc = a_neigh.m_facLoc;
    Vector<int> facneighDir = a_neigh.m_facneighDir;
    Vector<Side::LoHiSide> facneighSid = a_neigh.m_facneighSid;
    const Vector<Real> facneighFrac = a_neigh.m_facneighFrac;
    Vector<IndMomDIM> fac2volmom;
    MomentIterator<DIM, order> momit;
    int nfac = facmom.size();
    fac2volmom.resize(nfac);

    //Adding face moments
    for(int iface = 0; iface < nfac; iface++)
    {
      int irow = iface;

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

      //fac2volmom[iface].shift(diffrv);
      shiftMomentToLoc(fac2volmom[iface], facneighLoc[iface], a_faceLoc, a_dx);

      if(Chombo4::LAPACKMatrix::s_verbose)
      {
        pout() << "irow = " << irow << endl;
      }
      int mrow = irow + a_start[0];
      for(momit.reset(); momit.ok(); ++momit)
      {
        //need to use the lower order moment declaration to
        //get the correct column number
        int icol = IndexedMoments<DIM, order>::indexOf(momit());
        int mcol = icol + a_start[1];
        if (facneighFrac[iface] < m_eps)
        {
          Chombo4::MayDay::Warning("Zeroing all but the first element of a face moment row");
          if (icol == 0)
          {
            a_Mblock(mrow, mcol) = 1.0;
          }
          else
          {
            a_Mblock(mrow, mcol) = 0.0;
          }
        }
        else
        {
          a_Mblock(mrow, mcol) = fac2volmom[iface][momit()]/(facneighFrac[iface]);
        }
        if(Chombo4::LAPACKMatrix::s_verbose)
        {
          pout() << ",icol = "  << icol << ", momind = " << momit() << ", mat(irow,icol) = " << a_Mblock(mrow, mcol) << endl;
        }
      }
    }
  }

  //Compute the moments matrix portion for cut faces in the neighborhood
  template<CENTERING cent, int order>
  void
  BaseWLSFlux<cent, order>::
  getShiftedMomentMatrixEBFaces(Chombo4::LAPACKMatrix  & a_Mblock,
                                RealVect                 a_faceLoc,
                                Real                     a_dx,
                                const Neighbors<order> & a_neigh,
                                string                   a_ebbcName,
                                IndexTM<int, 2>         a_start) const
  {
    PR_TIME("getShiftedMomentMatrixEBFaces");
    using Chombo4::pout;
    Vector<IndMomDIM> ebfmom = a_neigh.m_ebfMom;
    const Vector<RealVect> ebfneighLoc = a_neigh.m_ebfLoc;
    Vector<Vector<IndMomDIM> > ebfneighNorm = a_neigh.getEBNormMoments();
    const Vector<Real> ebfneighFrac = a_neigh.m_ebfneighFrac;
    MomentIterator<DIM, order> momit;
    int nebf = ebfmom.size();

    //Adding eb moments
    for(int iebf = 0; iebf < nebf; iebf++)
    {
      int irow = iebf;
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
        pout() << "EB boundary (" << a_ebbcName << ") unrecognized" << std::endl;
        PROTO_ASSERT(false, "could not find ebBC string ");
      }

      shiftMomentToLoc(ebfmom[iebf], ebfneighLoc[iebf], a_faceLoc, a_dx);

      if(Chombo4::LAPACKMatrix::s_verbose)
      {
        pout() << "irow = " << irow << endl;
      }
      int mrow = irow + a_start[0];
      for(momit.reset(); momit.ok(); ++momit)
      {
        //need to use the lower order moment declaration to
        //get the correct column number
        int icol = IndexedMoments<DIM, order>::indexOf(momit());
        int mcol = icol + a_start[1];
        if (ebfneighFrac[iebf] < m_eps)
        {
          Chombo4::MayDay::Warning("Zeroing all but the first element of a eb face moment row");
          if (icol == 0)
          {
            a_Mblock(mrow, mcol) = 1.0;
          }
          else
          {
            a_Mblock(mrow, mcol) = 0.0;
          }
        }
        else
        {
          a_Mblock(mrow, mcol) = ebfmom[iebf][momit()]/(ebfneighFrac[iebf]);
        }
        if(Chombo4::LAPACKMatrix::s_verbose)
        {
          pout() << ",icol = "  << icol << ", momind = " << momit() << ", mat(irow,icol) = " << a_Mblock(mrow, mcol) << endl;
        }
      }
    }
  }

  //Compute the moments matrix portion for cut faces in the neighborhood
  template<CENTERING cent, int order>
  void
  BaseWLSFlux<cent, order>::
  getShiftedMomentMatrixEBNormalFaces(Chombo4::LAPACKMatrix  & a_Mblock,
                                      RealVect                 a_faceLoc,
                                      Real                     a_dx,
                                      const Neighbors<order> & a_neigh,
                                      int                      a_dir,
                                      string                   a_ebbcName,
                                      IndexTM<int, 2>         a_start) const
  {
    PR_TIME("getShiftedMomentMatrixEBNormalFaces");
    using Chombo4::pout;
    Vector<IndMomDIM> ebfmom = a_neigh.m_ebfMom;
    Vector<RealVect> ebfneighLoc = a_neigh.m_ebfLoc;
    Vector<Vector<IndMomDIM> > ebfneighNorm = a_neigh.getEBNormMoments();
    const Vector<Real> ebfneighFrac = a_neigh.m_ebfneighFrac;
    MomentIterator<DIM, order> momit;
    int nebf = ebfmom.size();

    //Adding eb moments
    for(int iebf = 0; iebf < nebf; iebf++)
    {
      int irow = iebf;
      shiftMomentToLoc(ebfneighNorm[iebf][a_dir], ebfneighLoc[iebf], a_faceLoc, a_dx);

      if(Chombo4::LAPACKMatrix::s_verbose)
      {
        pout() << "irow = " << irow << endl;
      }
      int mrow = irow + a_start[0];
      for(momit.reset(); momit.ok(); ++momit)
      {
        //need to use the lower order moment declaration to
        //get the correct column number
        int icol = IndexedMoments<DIM, order>::indexOf(momit());
        int mcol = icol + a_start[1];
        if ((ebfneighFrac[iebf] < m_eps) && (icol != 0))
        {
          Chombo4::MayDay::Warning("Zeroing all but the first element of a eb normal moment row");
          a_Mblock(mrow, mcol) = 0.0;
        }
        else
        {
          a_Mblock(mrow, mcol) = -ebfneighNorm[iebf][a_dir][momit()]/(ebfneighFrac[iebf]); // FIXME sign
        }

        if(Chombo4::LAPACKMatrix::s_verbose)
        {
          pout() << ",icol = "  << icol << ", momind = " << momit() << ", mat(irow,icol) = " << a_Mblock(mrow, mcol) << endl;
        }
      }
    }
  }

  template<CENTERING cent, int order>
  inline void
  BaseWLSFlux<cent, order>::
  checkSubMatrixSize(Chombo4::LAPACKMatrix  & a_Mblock,
                     const IndexTM<int, 2>   a_blockSize,
                     const IndexTM<int, 2>   a_start) const
  {
    // check the rows
    CH_assert(a_start[0] > 0);
    CH_assert(a_blockSize[0] > 0);
    CH_assert(a_start[0] + a_blockSize[0] < a_Mblock.numRow());
    // check the columns
    CH_assert(a_start[1] > 0);
    CH_assert(a_blockSize[1] > 0);
    CH_assert(a_start[1] + a_blockSize[1] < a_Mblock.numCol());
  }


  template<CENTERING cent, int order>
  inline void
  BaseWLSFlux<cent, order>::
  shiftMomentToLoc(IndMomDIM & a_moment,
                   RealVect    a_origLoc,
                   RealVect    a_newLoc,
                   Real        a_dx) const
  {
    RealVect diff = a_origLoc - a_newLoc;
    IndexTM<Real, DIM> diffrv;
    for(int dir = 0; dir != DIM; dir++)
      {
        diffrv[dir] = diff[dir]/a_dx;
      }
    a_moment.shift(diffrv);
  }

  //Compute rhs matrix for a face
  template<CENTERING cent, int order>
  void
  LapWLSFlux<cent, order>::
  getFaceFluxMatrix(Chombo4::LAPACKMatrix    & a_fluxMatrix,
                    const EBIndex<cent>      & a_face,
                    RealVect                 & a_faceLoc,
                    IndMomSDMinOne           & a_faceMom,
                    Real                     a_signNormal,
                    Real                     a_dx) const
  {
    PR_TIME("getFaceFluxMatrix");
    using Chombo4::pout;
    //Real areaScaleFactor  = getRegularArea(a_dx);

    RealVect startLoc;
    for (int dir = 0; dir < DIM; dir++)
    {
      startLoc[dir] = a_dx * (0.5 + a_face.m_pt[dir]);
      if (a_face.direction() == dir) startLoc[dir] = a_dx * (a_face.m_pt[dir]);
    }
    // shift center of moments to a_faceLoc.
    // the moments are centered at the face center initially
    RealVect diff = startLoc;
    diff -= a_faceLoc;
    IndexTM<Real, DIM-1> diffrv;
    int index = 0;
    for(int dir = 0; dir < DIM; dir++)
    {
      if(dir != a_face.direction())
      {
        diffrv[index] = diff[dir]/a_dx;
        index++;
      }
    }
    a_faceMom.shift(diffrv);

    MomentIterator<DIM, order> momit;
    a_fluxMatrix.define(1, momit.size());
    a_fluxMatrix.setVal(0.);
    const int& faceDir = a_face.direction();
    if (m_verbose) {pout() << " RHS - Face Flux Matrix " << endl;}
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
        for(int dir = 0; dir < DIM; dir++)
        {
           if(dir  != faceDir)
           {
             pfaceshift[index] =  p[dir];
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
        pout() << ",icol = "  << jmom <<  ", Q(0,icol) = " << a_fluxMatrix(0, jmom) << endl;
      }
    }
  }

  //Compute rhs matrix for a EB face
  template<CENTERING cent, int order>
  void
  LapWLSFlux<cent, order>::
  getEBFluxMatrix(Chombo4::LAPACKMatrix    & a_fluxMatrix,
                  const EBIndex<cent>      & a_face,
                  RealVect                 & a_faceLoc,
                  Vector<IndMomDIM>        & a_ebnorm,
                  Real                     a_dx) const
  {
    PR_TIME("getEBFluxMatrix");
    MomentIterator<DIM, order> momit;
    a_fluxMatrix.define(1, momit.size());
    a_fluxMatrix.setVal(0.);

    RealVect startLoc;
    for (int dir = 0; dir < DIM; dir++)
    {
      startLoc[dir] = a_dx * (0.5 + a_face.m_pt[dir]);
    }

    for(int dir = 0; dir < DIM; dir++)
    {

      RealVect diff = startLoc;
      diff -= a_faceLoc;
      IndexTM<Real, DIM> diffrv;
      for(int dir = 0; dir < DIM; dir++)
      {
        diffrv[dir] = diff[dir]/a_dx;
      }
      a_ebnorm[dir].shift(diffrv);

      // the contribution to the flux matrix from the normal in this particular direction
      Chombo4::LAPACKMatrix fluxIncr(1, momit.size());
      fluxIncr.setVal(0.);
      for(momit.reset(); momit.ok(); ++momit)
      {
        IndexTM<int, DIM> p = momit();
        Real pd = p[dir];
        //get column number associated with this particular moment
        int jmom = IndexedMoments<DIM, order>::indexOf(p);
        if(pd > 0)
        {
          IndexTM<int, DIM> pshift = p - BASISV_TM<int, DIM>(dir);
          fluxIncr(0, jmom) = pd*a_ebnorm[dir][pshift];
        }
        else
        {
          fluxIncr(0, jmom) = 0;
        }
      }
      a_fluxMatrix += fluxIncr;
    }
    if (m_verbose) { pout() << " RHS - EB Flux Matrix " << endl;}
    for(int dir = 0; dir < DIM; dir++)
    {
      for(momit.reset(); momit.ok(); ++momit)
      {
        IndexTM<int, DIM> p = momit();
        //Real pd = p[dir];
        //get column number associated with this particular moment
        int jmom = IndexedMoments<DIM, order>::indexOf(p);
        if(Chombo4::LAPACKMatrix::s_verbose)
        {
          pout() << ",icol = "  << jmom <<  ", Q(0,icol) = " << a_fluxMatrix(0, jmom) << endl;
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
    using Chombo4::pout;
    Vector<int> weightDisVol = a_neigh.m_weightDisVol;
    Vector<int> weightDisFac = a_neigh.m_weightDisFac;
    Vector<int> weightDisEB = a_neigh.m_weightDisEB;

    int nrows = weightDisVol.size()+weightDisFac.size()+weightDisEB.size();
    a_weights.define(nrows, nrows);
    a_weights.setToIdentity(); //this sets all the off diagonals to zero

    // int max_rad = 0.0;
    // for(int ineigh = 0; ineigh < weightDisVol.size(); ineigh++)
    // {
    //   if (max_rad<weightDisVol[ineigh])
    //   {
    //     max_rad = weightDisVol[ineigh];
    //   }
    // }
    //int powerE;
    //if (max_rad-1 == m_radius) powerE = m_powerE;
    //else powerE = m_powerE - (2.0*(max_rad-m_radius-1));

    // add the volume weights
    IndexTM<int, 2> start = 0*IndexTM<int, 2>::Unit;
    getWeightsMatrixComponent(a_weights, weightDisVol, start);
    // add the regular face weights
    start = (weightDisVol.size())*IndexTM<int, 2>::Unit;
    getWeightsMatrixComponent(a_weights, weightDisFac, start);
    // add the eb face weights
    start = (weightDisVol.size() + weightDisFac.size())*IndexTM<int, 2>::Unit;
    getWeightsMatrixComponent(a_weights, weightDisEB, start);

    if(Chombo4::LAPACKMatrix::s_verbose)
      {
        pout() << " Weights Matrix " << endl;
        a_weights.poutDiagMatlab();
      }
  }

  //Compute weight matrix for volume contributions
  template<CENTERING cent, int order>
  void
  BaseWLSFlux<cent, order>::
  getWeightsMatrixComponent(Chombo4::LAPACKMatrix  & a_weights,
                            const Vector<int>      & a_weightDist,
                            IndexTM<int, 2>         a_start) const
  {
    PR_TIME("getWeightsMatrixFace");

    int powerE = m_stenAttb.weightPower();

    for(int ineigh = 0; ineigh != a_weightDist.size(); ineigh++)
    {
      int dist = a_weightDist[ineigh];
      Real weight = POW(dist, powerE);
      a_weights(ineigh+a_start[0], ineigh+a_start[1]) = weight;
      // sane bounds for weights - not a hard rule, just the convention we use
      CH_assert((weight > 0.) && (weight <= 1.));
    }
  }

  //   //Compute stencil coefficients vector
  // // Solves for C = G*(W*M)^-1*W
  // template<CENTERING cent, int order>
  // void
  // BaseWLSFlux<cent, order>::
  // buildCompositeSystem(Chombo4::LAPACKMatrix        & a_C,
  //                      const Chombo4::LAPACKMatrix  & a_G,
  //                      const Chombo4::LAPACKMatrix  & a_M,
  //                      const Chombo4::LAPACKMatrix  & a_W) const
  // {
  //   Chombo4::LAPACKMatrix C, G, M, W;

  // }


  //Compute stencil coefficients vector
  // Solves for C = G*(W*M)^-1*W
  template<CENTERING cent, int order>
  void
  BaseWLSFlux<cent, order>::
  computeWLSStencilValues(Chombo4::LAPACKMatrix        & a_stencilVals,
                          const Chombo4::LAPACKMatrix  & a_evalOp,
                          const Chombo4::LAPACKMatrix  & a_Moments,
                          const Chombo4::LAPACKMatrix  & a_weights) const
  {
    PR_TIME("computeWLSStencilValues");
    Chombo4::LAPACKMatrix MT, WT, MTWT;
    MT = a_Moments;
    MT.transpose();
    WT = a_weights;
    WT.transpose(); // FIXME diagonal so this shouldn't matter. The bigger issue is we store the entire W
    multiply(MTWT, MT, WT);

    //if (m_verbose)
      // {
      //   // invert and check
      //   Chombo4::LAPACKMatrix Minv, I_M, MW, MWinv, I_MW;
      //   // Minv = a_Moments;
      //   // Minv.pseudoInvertUsingQR();
      //   // multiply(I_M, Minv, a_Moments);
      //   // a_Moments.checkConditionNumber();
      //   // pout() << "Show M*M^-1 " << std::endl;
      //   // I_M.poutAll();

      //   MW = MTWT;
      //   MW.transpose();
      //   MWinv = MW;
      //   MW.pseudoInvertUsingSVD(10,
      //                             1.0e-13);
      //   multiply(I_MW, MW, MWinv);
      //   MW.checkConditionNumber();
      //   pout() << "Show WM*WM^-1 " << std::endl;
      //   I_MW.poutAll();
      // }

    Chombo4::LAPACKMatrix QT;
    QT = a_evalOp;
    QT.transpose();

    Chombo4::LAPACKMatrix X, sT;

    int info = solveLSTSVD(X, MTWT, QT,
                           10,
                           1.0e-13); // FIXME why this tolerance?

    if(info != 0)
    {
      Chombo4::MayDay::Warning("non zero info returned");
    }

    multiply(sT, WT, X);
    a_stencilVals = sT;
    a_stencilVals.transpose();
  }

  //Compute inhomogeneous term
  template<CENTERING cent, int order>
  Real
  LapWLSFlux<cent, order>::
  computeInhomogTerm(const Chombo4::LAPACKMatrix    & a_stencil,
                     const Neighbors<order>         & a_neigh,
                     const int                        a_nvolneigh,
                     const EBDomainBC<1>     & a_domainBC,
                     Real                             a_dx) const
  {
    PR_TIME("computeInhomogTerm");
    using Chombo4::pout;
    Real inhomogTermReg = computeInhomogTermReg(a_stencil,
                                                a_neigh,
                                                a_nvolneigh,
                                                a_domainBC,
                                                a_dx);

    Real inhomogTermEB = computeInhomogTermEB(a_stencil,
                                                a_neigh,
                                                a_nvolneigh,
                                                a_domainBC,
                                                a_dx);

    return inhomogTermReg + inhomogTermEB;
  }


  template<CENTERING cent, int order>
  Real
  LapWLSFlux<cent, order>::
  computeInhomogTermReg(const Chombo4::LAPACKMatrix     & a_stencil,
                        const Neighbors<order>          & a_neigh,
                        const int                         a_nvolneigh,
                        const EBDomainBC<1>      & a_domainBC,
                        Real                            a_dx) const
  {
    PR_TIME("computeInhomogTerm");
    using Chombo4::pout;
    Vector<IndMomSDMinOne> a_facmom = a_neigh.m_facMom;
    Vector<RealVect> a_facneighLoc = a_neigh.m_facLoc;
    Vector<int> a_facneighDir = a_neigh.m_facneighDir;
    Vector<Side::LoHiSide> a_facneighSid = a_neigh.m_facneighSid;
    Vector<Real> facneighFrac = a_neigh.m_facneighFrac;
    Real areaScaleFactor = getRegularArea(a_dx);

    Real inhomogTerm = 0.0;

    for(int iface = 0; iface < a_neigh.sizeFaces(); iface++)
    {
      Real bcValue = 0;
      RealVect a_faceLocDiff = a_facneighLoc[iface];
      int faceDir = a_facneighDir[iface];
      Side::LoHiSide faceSide = a_facneighSid[iface];
      string bcName = m_domainBC.getBCtypes().getRegComp(faceDir, faceSide);
      if(bcName == string("Neumann"))
      {
        bcValue = NeumannExactSoln<order>::
                               getFaceIntegral(a_facmom[iface],
                                               a_faceLocDiff,
                                               faceDir,
                                               m_domainBC.getBCfuncs().getRegComp(faceDir, faceSide),
                                               a_dx);
      }
      else if(bcName == string("Dirichlet"))
      {
        bcValue = DirichletExactSoln<order>::
                               getFaceIntegral(a_facmom[iface],
                                               a_faceLocDiff,
                                               faceDir,
                                               m_domainBC.getBCfuncs().getRegComp(faceDir, faceSide),
                                               a_dx);
      }
      else
      {
        PROTO_ASSERT(false, "could not find domBC string");
      }

      if (facneighFrac[iface] < this->m_eps)
      {
        bcValue = 0.0;
      }
      else
      {
        bcValue /= (facneighFrac[iface]*areaScaleFactor);
      }

      inhomogTerm += bcValue*a_stencil(0, a_neigh.indexFromFace(iface));
    }

    return inhomogTerm;
  }

  template<CENTERING cent, int order>
  Real
  LapWLSFlux<cent, order>::
  computeInhomogTermEB(const Chombo4::LAPACKMatrix     & a_stencil,
                       const Neighbors<order>          & a_neigh,
                       const int                         a_nvolneigh,
                       const EBDomainBC<1>      & a_domainBC,
                       Real                            a_dx) const
  {
    PR_TIME("computeInhomogTermEB");
    using Chombo4::pout;
    Vector<IndMomDIM> a_ebfmom= a_neigh.m_ebfMom;
    Vector<Vector<IndMomDIM> > a_ebfneighNorm = a_neigh.m_ebfneighNorm;
    Vector<RealVect> a_ebfneighLoc = a_neigh.m_ebfLoc;
    Vector<Real> ebfneighFrac = a_neigh.m_ebfneighFrac;
    Real areaScaleFactor = getRegularArea(a_dx);

    Real inhomogTerm = 0.0;

    for(int iface = 0; iface < a_neigh.sizeEB(); iface++)
    {
      Real bcValue = 0;
      RealVect a_ebfLocDiff = a_ebfneighLoc[iface];
      string a_ebbcName = a_domainBC.getBCtypes().getEBComp();
      if(a_ebbcName == string("Neumann"))
      {
        bcValue = NeumannExactSoln<order>::
                               getEBIntegral(a_ebfneighNorm[iface],
                                             a_ebfLocDiff,
                                             a_domainBC.getBCfuncs().getEBComp(),
                                             a_dx);
      }
      else if(a_ebbcName == string("Dirichlet"))
      {
        bcValue = DirichletExactSoln<order>::
                               getEBIntegral(a_ebfmom[iface],
                                             a_ebfLocDiff,
                                             a_domainBC.getBCfuncs().getEBComp(),
                                             a_dx);
      }
      else
      {
        PROTO_ASSERT(false, "could not find domBC string");
      }

      if (ebfneighFrac[iface] < this->m_eps)
      {
        bcValue = 0.0;
      }
      else
      {
        bcValue /= (ebfneighFrac[iface]*areaScaleFactor);
      }

      inhomogTerm += bcValue*a_stencil(0, a_neigh.indexFromEB(iface));
    }

    return inhomogTerm;
  }
  
  // //Compute inhomogeneous term
  // template<CENTERING cent, int order>
  // Real
  // BaseWLSFlux<cent, order>::
  // computeInhomogTerm(const Chombo4::LAPACKMatrix    & a_stencil,
  //                    const Neighbors<order>         & a_neigh,
  //                    const int                        a_nvolneigh,
  //                    const EBDomainBC<DIM>   & a_domainBC,
  //                    int                              a_bcComp,
  //                    Real                           a_dx) const
  // {
  //   PR_TIME("computeInhomogTerm");
  //   using Chombo4::pout;
  //   Vector<IndMomSDMinOne> a_facmom = a_neigh.m_facMom;
  //   Vector<IndMomDIM> a_ebfmom= a_neigh.m_ebfMom;
  //   Vector<Vector<IndMomDIM> > a_ebfneighNorm = a_neigh.m_ebfneighNorm;
  //   Vector<RealVect> a_facneighLoc = a_neigh.m_facLoc;
  //   Vector<int> a_facneighDir = a_neigh.m_facneighDir;
  //   Vector<Side::LoHiSide> a_facneighSid = a_neigh.m_facneighSid;
  //   Vector<RealVect> a_ebfneighLoc = a_neigh.m_ebfLoc;
  //   Vector<Real> ebfneighFrac = a_neigh.m_ebfneighFrac;
  //   Vector<Real> facneighFrac = a_neigh.m_facneighFrac;
  //   Real areaScaleFactor = getRegularArea(a_dx);

  //   Vector<Real> domBCValues, ebBCValues;
  //   int nrows = a_facmom.size();
  //   domBCValues.resize(nrows);
  //   for(int iface = 0; iface < nrows; iface++)
  //   {
  //     RealVect a_faceLocDiff = a_facneighLoc[iface];
  //     int a_faceDir = a_facneighDir[iface];
  //     Side::LoHiSide a_faceSid = a_facneighSid[iface];
  //     auto exactSolSide = a_domainBC.getBCfuncs().getRegComp(a_faceDir, a_faceSid, a_bcComp);
  //     int index = ebp_index(a_faceDir, a_faceSid);
  //     if(a_dombcName[index] == string("Neumann"))
  //     {
  //       domBCValues[iface] = NeumannExactSoln<order>::
  //                              getFaceIntegral(a_facmom[iface],
  //                                              a_faceLocDiff,
  //                                              a_faceDir,
  //                                              exactSolSide,
  //                                              a_dx);
  //     }
  //     else if(a_dombcName[index] == string("Dirichlet"))
  //     {
  //       domBCValues[iface] = DirichletExactSoln<order>::
  //                              getFaceIntegral(a_facmom[iface],
  //                                              a_faceLocDiff,
  //                                              a_faceDir,
  //                                              exactSolSide,
  //                                              a_dx);
  //     }
  //     else
  //     {
  //       PROTO_ASSERT(false, "could not find domBC string");
  //     }

  //     if (facneighFrac[iface]<1.0e-11)
  //     {
  //       domBCValues[iface] = 0.0;
  //     }
  //     else
  //     {
  //       domBCValues[iface] = domBCValues[iface]/(facneighFrac[iface]*areaScaleFactor);
  //     }
  //   }

  //   nrows = a_ebfmom.size();
  //   ebBCValues.resize(nrows);
  //   auto exactSolEB = a_domainBC.getBCfuncs().getEBComp(a_bcComp);
  //   for(int ivof = 0; ivof < nrows; ivof++)
  //   {
  //     RealVect a_ebfLocDiff = a_ebfneighLoc[ivof];
  //     if(a_ebbcName == string("Neumann"))
  //     {
  //       ebBCValues[ivof] = NeumannExactSoln<order>::
  //                              getEBIntegral(a_ebfneighNorm[ivof],
  //                                            a_ebfLocDiff,
  //                                            exactSolEB,
  //                                            a_dx);
  //     }
  //     else if(a_ebbcName == string("Dirichlet"))
  //     {
  //       ebBCValues[ivof] = DirichletExactSoln<order>::
  //                              getEBIntegral(a_ebfmom[ivof],
  //                                            a_ebfLocDiff,
  //                                            exactSolEB,
  //                                            a_dx);
  //     }
  //     else
  //     {
  //       PROTO_ASSERT(false, "could not find domBC string");
  //     }

  //     if (ebfneighFrac[ivof]<1.0e-11)
  //     {
  //       ebBCValues[ivof] = 0.0;
  //     }
  //     else
  //     {
  //       ebBCValues[ivof] = ebBCValues[ivof]/(ebfneighFrac[ivof]*areaScaleFactor);
  //     }
  //   }

  //   int nfac = a_facmom.size();
  //   int nebf = a_ebfmom.size();

  //   Real inhomogTerm = 0.0;
  //   for(int iface = 0; iface < nfac; iface++)
  //   {
  //     inhomogTerm += domBCValues[iface]*a_stencil(0, a_nvolneigh + iface);
  //   }
  //   for(int iface = 0; iface < nebf; iface++)
  //   {
  //     inhomogTerm += ebBCValues[iface]*a_stencil(0, a_nvolneigh + nfac + iface);
  //     if (m_verbose)
  //       {
  //         pout() << " ebbc val " << ebBCValues[iface] << "  coef "
  //                << a_stencil(0, a_nvolneigh + nfac + iface) << endl;
  //       }
  //   }

  //   return inhomogTerm;
  // }

  //Compute face centered flux stencil
  template<CENTERING cent, int order>
  void
  LapWLSFlux<cent, order>::
  getFaceCenteredStencil(LocalStencil<CELL, Real> & a_sten,
                         Real                       & a_inhomogTerm, //for inhomogeneous bcs
                         const EBIndex<CELL>        & a_vof,
                         const EBIndex<cent>        & a_face,
                         const EBGraph              & a_graph,
                         const MomentHolder<order>  & a_momentData,
                         const int                    a_dir,
                         IndMomSDMinOne             & a_faceMom,
                         int                        & a_neighType,
                         Real                       a_signNormal,
                         Real                       a_dx) const
  {
    PR_TIME("getFaceCenteredStencil");
    using Chombo4::pout;
    Chombo4::LAPACKMatrix::s_verbose = m_verbose;
    //Chombo4::LAPACKMatrix::s_verbose = true;

    // FIXME, no reason to unpack these here
    // further refactoring needed
    string a_ebbcName = m_domainBC.getBCtypes().getEBComp();
    string a_dombcName[2*DIM];
    m_domainBC.getBCtypes().fillFlatArray(a_dombcName);

    Neighbors<order> neighs;

    RealVect faceLoc;
    for (int dir = 0; dir < DIM; dir++)
    {
      faceLoc[dir] = a_dx * (0.5 + a_face.m_pt[dir]);
      if (a_dir == dir) faceLoc[dir] = a_dx * (a_face.m_pt[dir]);
    }

    int radius = m_stenAttb.radius();
    Vector<EBIndex<CELL> > vols;
    MomentIterator<DIM,order> momit;
    int grows = 0;
    // HACK - this grows outlet stencils to make sure they are well posed, but there should be a better way to check this. FIXME Hard coded for high side boundaries
    //if ((a_vof.m_pt[a_dir] >= a_graph.getDomain().high()[a_dir]) && (a_dombcName[ebp_index(a_dir, Side::Hi)] == string("no_bcs")))
      // {
      //   radius += 2;
      //   pout() << "growing no_bc boundary at " << a_vof.m_pt << std::endl;
      // }
    // REVIEW Make sure stencils are big enough. Size is not a sufficient measure, this should check the stencil is solvable
    while (neighs.size() <= momit.size()+DIM)
      {
        vols = neighs.template getVoFsNeighbors<cent>(a_face, a_vof, radius, a_graph, a_neighType);

        neighs.getMomentsInfoFromNeighbors(vols, a_graph, a_momentData, m_domainBC, a_dx);

        // if (grows>0) pout() << "neigh set tight for " << a_vof.m_pt << ", growing radius by 1" << endl;
        radius += 1;
        grows += 1;
      }
    CH_assert(vols.size() > momit.size());
    CH_assert(neighs.size() > momit.size());

    // first fill M with moments shifted to the center of the face
    Chombo4::LAPACKMatrix Mvol;
    this->getShiftedMomentMatrix(Mvol, faceLoc, a_dx,
                           neighs, 
                           a_ebbcName, a_dombcName);
    // Real cn = getInverseOfConditionNumber(Mvol);
    // if (cn < 1e-5)
    //   pout() << "condition number small in " << a_vof.m_pt << a_dir << std::endl;
    // RHS of stencil equation
    // The G vector
    Chombo4::LAPACKMatrix Qshift;
    getFaceFluxMatrix(Qshift, a_face, faceLoc, a_faceMom, a_signNormal, a_dx);

    // Weights matrix W
    Chombo4::LAPACKMatrix weights;
    RealVect centroid_pos;
    centroid_pos[0] = a_faceMom[1]/a_faceMom[0];
#if DIM==3
    centroid_pos[1] = a_faceMom[1+order]/a_faceMom[0];
#endif
    int centr_iter = 0;
    for (int dir = 0; dir < DIM; dir++)
    {
      if (a_dir != dir)
      {
        faceLoc[dir] += centroid_pos[centr_iter]*a_dx;
        centr_iter++;
      }
    }
    if (m_verbose) {pout()<<" centroid comp "<<faceLoc<<endl;}
    this->getWeightsMatrixFace(weights, neighs, faceLoc);

    // Solve for the stencil vector S
    Chombo4::LAPACKMatrix phiCoef;
    this->computeWLSStencilValues(phiCoef, Qshift, Mvol, weights);

    //Filling stencil, prune out insignificant values
    a_sten.clear();
    Real stenMaxAbs = -1e300;
    for(int ivof = 0; ivof < vols.size(); ivof++)
      {
        stenMaxAbs = std::max(std::abs(phiCoef(0, ivof)), stenMaxAbs);
      }
    Real coef_acc = 0.0;
    if (m_verbose) { pout() << " Coefs Matrix " << endl;}
    for(int ivof = 0; ivof < vols.size(); ivof++)
    {
      Real thisCo = phiCoef(0, ivof);
      if (std::abs(thisCo) > stenMaxAbs*1e-15) { a_sten.add(vols[ivof], thisCo); }
      //else {pout()<< "ignoring stencil coef of "<< thisCo <<" Coef: "<<thisCo<<endl;}
      if (m_verbose){
        pout()<<"stencil coefs for Vof: "<<vols[ivof].m_pt<<" Coef: "<<thisCo<<endl;
      }
      coef_acc += thisCo;
    }

    int nvolneigh = vols.size();
    a_inhomogTerm = computeInhomogTerm(phiCoef, neighs,
                                       nvolneigh, 
                                       m_domainBC,
                                       a_dx);
    if (m_verbose)
      {
        pout() << "For phi=1 " << coef_acc << " " << a_inhomogTerm << " "
               << coef_acc+a_inhomogTerm<<endl;
      }
  }

  //Compute EB face stencil
  template<CENTERING cent, int order>
  void
  LapWLSFlux<cent, order>::
  getEBFluxStencil(LocalStencil<CELL, Real>      & a_stencil,
                   Real                          & a_inhomogTerm, //for inhomogeneous bcs
                   const EBIndex<CELL>           & a_vof,
                   const EBIndex<cent>           & a_face,
                   const EBGraph                 & a_graph,
                   const MomentHolder<order>     & a_momentData,
                   Vector<IndMomDIM>             & a_ebnorm,
                   int                             a_neighType,
                   const Real                    a_dx) const
  {
    PR_TIME("getEBFluxStencil");
    using Chombo4::pout;
    Chombo4::LAPACKMatrix::s_verbose = m_verbose;
    // Chombo4::LAPACKMatrix::s_verbose = true;

    // FIXME, no reason to unpack these here
    // further refactoring needed
    string a_ebbcName = m_domainBC.getBCtypes().getEBComp();
    string a_dombcName[2*DIM];
    m_domainBC.getBCtypes().fillFlatArray(a_dombcName);

    if(a_ebbcName == string("Neumann"))
      {
        //neumann has no stencil, just an eval on the face
        Real areaScaleFactor = getRegularArea(a_dx);
        RealVect faceLoc;
        for (int idir = 0; idir < DIM; idir++)
          {
            faceLoc[idir] = a_dx * (0.5 + a_face.m_pt[idir]);
          }
        a_inhomogTerm = NeumannExactSoln<order>::
          getEBIntegral(a_ebnorm,
                        faceLoc,
                        m_domainBC.getBCfuncs().getEBComp(),
                        a_dx);
        a_inhomogTerm /= areaScaleFactor;

      }
    else if(a_ebbcName == string("Dirichlet"))
      {
        // Dirichlet does require a stencil, so do a similar dance as in the grid aligned faces
        Neighbors<order> neighs;

        RealVect faceLoc;
        for (int dir = 0; dir < DIM; dir++)
          {
            faceLoc[dir] = a_dx * (0.5 + a_face.m_pt[dir]);
          }
        int radius = m_stenAttb.radius();

        Vector<EBIndex<CELL> > vols;
        MomentIterator<DIM,order> momit;
        int grows = 0;
        // HACK - this grows outlet stencils to make sure they are well posed, but there should be a better way to check this. FIXME Hard coded for high side boundaries
        //if ((a_vof.m_pt[a_dir] >= a_graph.getDomain().high()[a_dir]) && (a_dombcName[ebp_index(a_dir, Side::Hi)] == string("no_bcs")))
        // {
        //   radius += 2;
        //   pout() << "growing no_bc boundary at " << a_vof.m_pt << std::endl;
        // }
        // REVIEW Make sure stencils are big enough. Size is not a sufficient measure, this should check the stencil is solvable
        while (neighs.size() <= momit.size()+DIM)
          {
            vols = neighs.template getVoFsNeighbors<cent>(a_face, a_vof, radius, a_graph, a_neighType);

            neighs.getMomentsInfoFromNeighbors(vols, a_graph, a_momentData, m_domainBC, a_dx);

            // if (grows>0) pout() << "neigh set tight for " << a_vof.m_pt << ", growing radius by 1" << endl;
            radius += 1;
            grows += 1;
          }
        CH_assert(vols.size() > momit.size());
        CH_assert(neighs.size() > momit.size());

        //first fill M with moments  shifted to the center of the face
        Chombo4::LAPACKMatrix  Mvol;
        this->getShiftedMomentMatrix(Mvol, faceLoc, a_dx,
                               neighs,
                               a_ebbcName, a_dombcName);

        // RHS of stencil equation
        Chombo4::LAPACKMatrix Qshift;
        getEBFluxMatrix(Qshift, a_face, faceLoc, a_ebnorm, a_dx);

        // Weights matrix
        Chombo4::LAPACKMatrix weights;
        this->getWeightsMatrixFace(weights, neighs, faceLoc);

        // Compute phiCoef
        Chombo4::LAPACKMatrix phiCoef;
        this->computeWLSStencilValues(phiCoef, Qshift, Mvol, weights);

        //Filling stencil
        a_stencil.clear();
        Real coef_acc = 0.0;
        if (m_verbose) { pout() << " Coefs Matrix " << endl;}
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
                                           m_domainBC,
                                           a_dx);
        if (m_verbose)
          {
            pout() << "For phi=1 " << coef_acc << " " << a_inhomogTerm << " "
                   << coef_acc+a_inhomogTerm<<endl;
          }
      }
    else
      {
        PROTO_ASSERT(false, "could not find EBBC string");
      }
  }


  template<CENTERING cent, int order>
  void
  DivWLSFlux<cent, order>::
  getFaceCenteredStencilFull(std::array<LocalStencil<CELL, Real>, DIM >& a_sten,
                             Real                       & a_inhomogTerm, //for inhomogeneous bcs
                             const EBIndex<CELL>        & a_vof,
                             const EBIndex<cent>        & a_face,
                             const EBGraph              & a_graph,
                             const MomentHolder<order>  & a_momentData,
                             const int                    a_dir,
                             IndMomSDMinOne             & a_faceMom,
                             int                        & a_neighType,
                             Real                       a_signNormal,
                             Real                       a_dx) const
  {
    PR_TIME("getFaceCenteredStencil");
    using Chombo4::pout;
    Chombo4::LAPACKMatrix::s_verbose = m_verbose;
    //Chombo4::LAPACKMatrix::s_verbose = true;

    // FIXME, no reason to unpack these here
    // further refactoring needed
    string a_ebbcName = m_domainBC.getBCtypes().getEBComp(a_dir);
    string a_dombcName[2*DIM];
    m_domainBC.getBCtypes().fillFlatArray(a_dombcName, a_dir);

    EBDomainBC<1> domainBCdir = m_domainBC.makeSingleComponent(a_dir);

    Neighbors<order> neighs;

    RealVect faceLoc;
    for (int dir = 0; dir < DIM; dir++)
    {
      faceLoc[dir] = a_dx * (0.5 + a_face.m_pt[dir]);
      if (a_dir == dir) faceLoc[dir] = a_dx * (a_face.m_pt[dir]);
    }

    int radius = m_stenAttb.radius();
    Vector<EBIndex<CELL> > vols;
    MomentIterator<DIM,order> momit;
    int grows = 0;
    // HACK - this grows outlet stencils to make sure they are well posed, but there should be a better way to check this. FIXME Hard coded for high side boundaries
    // if ((a_vof.m_pt[a_dir] >= a_graph.getDomain().high()[a_dir]) && (a_dombcName[ebp_index(a_dir, Side::Hi)] == string("no_bcs")))
    //   {
    //     radius += 2;
    //     pout() << "growing no_bc boundary at " << a_vof.m_pt << std::endl;
    //   }
    // REVIEW Make sure stencils are big enough. Size is not a sufficient measure, this should check the stencil is solvable
    while (neighs.size() <= momit.size()+DIM)
      {
        vols = neighs.getVoFsNeighbors(a_face, a_vof, radius, a_graph, a_neighType);

        neighs.getMomentsInfoFromNeighbors(vols, a_graph, a_momentData, domainBCdir, a_dx, m_verbose);

        // if (grows>0) pout() << "neigh set tight for " << a_vof.m_pt << ", growing radius by 1" << endl;
        radius += 1;
        grows += 1;
      }
    CH_assert(vols.size() > momit.size());
    CH_assert(neighs.size() > momit.size());


    // Real cn = getInverseOfConditionNumber(Mvol);
    // if (cn < 1e-5)
    //   pout() << "condition number small in " << a_vof.m_pt << a_dir << std::endl;
    // RHS of stencil equation
    // The G vector
    Chombo4::LAPACKMatrix Qshift;
    getFaceFluxMatrix(Qshift, a_face, faceLoc, a_faceMom, a_signNormal, a_dx);

    // // Solve for the stencil vector S // ORIGINAL VERSION
    // Chombo4::LAPACKMatrix phiCoef;
    // computeWLSStencilValues(phiCoef, Qshift, Mvol, weights);

    int nmoms = momit.size();
    int nvolneigh = neighs.sizeVols();
    int nfac = neighs.sizeFaces();
    int nebf = neighs.sizeEB();
    int ntotal = neighs.size();

    // build a new composite system
    Chombo4::LAPACKMatrix M;
    Chombo4::LAPACKMatrix W;
    Chombo4::LAPACKMatrix G;
    Chombo4::LAPACKMatrix S;

    // build the composite system
    // in 2D:
    // [ M_vols     0  ]
    // [ 0      M_vols ]
    // [ M_xface    0  ]
    // [ 0      M_yface]
    // [ M_xeb  M_yeb  ]
    {
      int ntotal_coup = DIM*(nvolneigh + nfac) + nebf;
      int nmoms_coup = DIM*nmoms;
      M.define(ntotal_coup, nmoms_coup);
      M.setVal(0.);
      W.define(ntotal_coup, ntotal_coup);
      W.setVal(0.);
      // fill the volume and face moments and weights
      for (int dir=0; dir!=DIM; dir++)
        {
          // build each block of the volume moments
          IndexTM<int, 2> start = IndexTM<int, 2>{dir*nvolneigh, dir*nmoms}; // row, col block beginning
          this->getShiftedMomentMatrixVolumes(M, faceLoc, a_dx, neighs, start);
          IndexTM<int, 2> startW = start[0]*IndexTM<int, 2>::Unit; // diagonal entry, shared row ordering with moments
          this->getWeightsMatrixComponent(W, neighs.m_weightDisVol, startW);

          // build each block of the regular face moments
          start = IndexTM<int, 2>{dir*nfac + DIM*(nvolneigh), dir*nmoms}; // row, col block beginning
          this->getShiftedMomentMatrixRegFaces(M, faceLoc, a_dx, neighs, a_dombcName, start);
          startW = start[0]*IndexTM<int, 2>::Unit; // diagonal entry, shared row ordering with moments
          this->getWeightsMatrixComponent(W, neighs.m_weightDisFac, startW);

          // build each block of the eb moments
          start = IndexTM<int, 2>{DIM*(nvolneigh + nfac), dir*nmoms}; // row, col block beginning
          this->getShiftedMomentMatrixEBNormalFaces(M, faceLoc, a_dx, neighs, dir, a_ebbcName, start);
        }
      // weights for EB section only get added once because they are coupled
      IndexTM<int, 2> startW = (DIM*(nvolneigh + nfac))*IndexTM<int, 2>::Unit; // diagonal entry, shared row ordering with moments
      this->getWeightsMatrixComponent(W, neighs.m_weightDisEB, startW);

      // The G matrix, is Qshift for this dir, and zero otherwise
      G.define(1, DIM*momit.size());
      G.setVal(0.);
      //const int& faceDir = a_face.direction();
      for(momit.reset(); momit.ok(); ++momit)
        {
          IndexTM<int, DIM> p = momit();
          // int pd = p[faceDir];
          int jmom = IndexedMoments<DIM, order>::indexOf(p);
          int offset = a_dir*momit.size();
          G(0, jmom+offset) = Qshift(0, jmom);
        }
    }
    // comment out to get old decoupled version
    this->computeWLSStencilValues(S, G, M, W);

    if (m_verbose)
      {
        pout() << "W matrix" << std::endl;
        W.poutDiagMatlab();
        pout() << "M matrix" << std::endl;
        M.poutMatlab();
        pout() << "G matrix" << std::endl;
        G.poutMatlab();
        pout() << "S matrix" << std::endl;
        S.poutMatlab();
      }

    for(int dir=0; dir!=DIM; dir++)
      {
        a_sten[dir].clear();
      }
    a_inhomogTerm = 0;
    //Filling stencil, prune out insignificant values
    //Real stenMaxAbs = S.maxNorm();

    Real coef_acc = 0.0;
    if (m_verbose) { pout() << " Coefs Matrix " << endl;}
    for(int dir=0; dir!=DIM; dir++)
      {
        // grab a single components of the stencils to evaluate
        //Chombo4::LAPACKMatrix Scomp(1, ntotal, (S.dataPtr()+dir*S.offset(0, ntotal)));
        Chombo4::LAPACKMatrix Scomp(1, ntotal);
        Scomp.setVal(0); // don't add eb part, do it later
        {
          for(int row_vol=0; row_vol!=nvolneigh; row_vol++)
            {
              int row = row_vol;
              int row_offset = row_vol + dir*nvolneigh; //ntotal;
              Scomp(0, row) = S(0, row_offset);
            }
          for(int row_fac=0; row_fac!=nfac; row_fac++) // do the regular face moments
            {
              int row = nvolneigh + row_fac;
              int row_offset = row_fac + dir*nfac + DIM*(nvolneigh);
              Scomp(0, row) = S(0, row_offset);
            }
        }

        // Add the stencils components for the U-dir component
        for(int ivof = 0; ivof < vols.size(); ivof++)
          {
            Real thisCo = Scomp(0, ivof);
            //if (std::abs(thisCo) > stenMaxAbs*1e-15) // prune small stencil values
              {
                a_sten[dir].add(vols[ivof], thisCo);
              }
            //else {pout()<< "ignoring stencil coef of "<< thisCo << endl;}
            if (m_verbose){
              pout() << "stencil coefs for Vof: " << vols[ivof].m_pt
                     << " dir: " << dir
                     <<" Coef: " << thisCo<<endl;
            }
            coef_acc += thisCo;
          }

        // add the regular face contributions of U-dir
        a_inhomogTerm += computeInhomogTermReg(Scomp,
                                               neighs,
                                               dir,
                                               m_domainBC,
                                               a_dx);
      }

    // add the EB face contribution
    //for(int dir=0; dir!=DIM; dir++)
      {
        Chombo4::LAPACKMatrix Scomp(1, ntotal);
        Scomp.setVal(0); // now do only the eb parts
        {
          for(int row_ebf=0; row_ebf!=nebf; row_ebf++) // do the eb moments
            {
              int row = nvolneigh + nfac + row_ebf;
              int row_offset = row_ebf + DIM*(nvolneigh + nfac);
              Scomp(0, row) = S(0, row_offset);
            }
        }
        a_inhomogTerm += computeInhomogTermEB(Scomp,
                                              neighs,
                                              m_domainBC,
                                              a_dx);
      }
    if (m_verbose)
      {
        pout() << "For phi=1 " << coef_acc << " " << a_inhomogTerm << " "
               << coef_acc+a_inhomogTerm<<endl;
      }
  }

  //Compute EB face stencil
  template<CENTERING cent, int order>
  void
  DivWLSFlux<cent, order>::
  getEBFluxStencilFull(std::array<LocalStencil<CELL, Real>, SpaceDim >& a_sten,
                       Real                          & a_inhomogTerm, //for inhomogeneous bcs
                       const EBIndex<CELL>           & a_vof,
                       const EBIndex<cent>           & a_face,
                       const EBGraph                 & a_graph,
                       const MomentHolder<order>     & a_momentData,
                       int                             a_neighType,
                       const Real                    a_dx) const
  {
    PR_TIME("getEBFluxStencil");
    using Chombo4::pout;
    Chombo4::LAPACKMatrix::s_verbose = m_verbose;
    // Chombo4::LAPACKMatrix::s_verbose = true;
    std::array<IndMomDIM, DIM> ebnorm;

    //auto& a_voludata = a_momentData.volumeMoments();
    auto& a_ebfadata = a_momentData.ebFaceMoments();
#if DIM >= 1
    //auto& a_xfacdata = a_momentData.template faceMoments<XFACE>();
    auto& a_ebnormxdata = a_momentData.ebNormal(0);
    ebnorm[0] = a_ebnormxdata(a_vof, 0);
#endif
#if DIM >= 2
    //auto& a_yfacdata = a_momentData.template faceMoments<YFACE>();
    auto& a_ebnormydata = a_momentData.ebNormal(1);
    ebnorm[1] = a_ebnormydata(a_vof, 0);
#endif
#if DIM >= 3
    //auto& a_zfacdata = a_momentData.template faceMoments<ZFACE>();
    auto& a_ebnormzdata = a_momentData.ebNormal(2);
    ebnorm[2] = a_ebnormzdata(a_vof, 0);
#endif

    // FIXME, I sure hope all these are the same
    string a_ebbcName = m_domainBC.getBCtypes().getEBComp(0);
    PROTO_ASSERT(m_domainBC.getBCtypes().getEBComp(0) == m_domainBC.getBCtypes().getEBComp(1),
                 "All EB BC's in all directions should be of the same type for divergence op");

    Real inhomogTerm = 0.0; //for inhomogeneous bcs

    IndMomDIM ebfmom = a_ebfadata(a_face, 0);
    Normalizor<order>::normalizeFaceMoment(ebfmom, a_dx);

    if(a_ebbcName == string("Dirichlet"))
      {
        // no stencil required, just a BC contribution
        RealVect faceLoc;
        for (int dir = 0; dir < DIM; dir++)
          {
            faceLoc[dir] = a_dx * (0.5 + a_face.m_pt[dir]);
          }
        inhomogTerm = 0;
        for (int dir=0; dir!=DIM; dir++)
          {
            Normalizor<order>::normalizeFaceMoment(ebnorm[dir], a_dx);
            inhomogTerm += DirichletExactSoln<order>::
              getEBIntegral(ebnorm[dir],
                            faceLoc,
                            m_domainBC.getBCfuncs().getEBComp(dir),
                            a_dx);
          }
        a_inhomogTerm = inhomogTerm/getRegularArea(a_dx);
      }
    else
      {
        PROTO_ASSERT(false, "could not find EBBC string");
      }

  }

//Compute rhs matrix for a face
  template<CENTERING cent, int order>
  void
  DivWLSFlux<cent, order>::
  getFaceFluxMatrix(Chombo4::LAPACKMatrix    & a_fluxMatrix,
                    const EBIndex<cent>      & a_face,
                    RealVect                 & a_faceLoc,
                    IndMomSDMinOne           & a_faceMom,
                    Real                     a_signNormal,
                    Real                     a_dx) const
  {
    PR_TIME("getFaceFluxMatrix");
    using Chombo4::pout;
    RealVect startLoc;
    for (int dir = 0; dir < DIM; dir++)
    {
      startLoc[dir] = a_dx * (0.5 + a_face.m_pt[dir]);
      if (a_face.direction() == dir) startLoc[dir] = a_dx * (a_face.m_pt[dir]);
    }
    // shift center of moments to a_faceLoc.
    // the moments are centered at the face center initially
    RealVect diff = startLoc;
    diff -= a_faceLoc;
    IndexTM<Real, DIM-1> diffrv;
    int index = 0;
    for(int dir = 0; dir < DIM; dir++)
    {
      if(dir != a_face.direction())
      {
        diffrv[index] = diff[dir]/a_dx;
        index++;
      }
    }
    a_faceMom.shift(diffrv);

    MomentIterator<DIM, order> momit;
    a_fluxMatrix.define(1, momit.size());
    a_fluxMatrix.setVal(0.);
    const int& faceDir = a_face.direction();
    if (this->m_verbose){ pout() << " RHS - Face Flux Matrix Div" << endl;}
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
        for(int dir = 0; dir < DIM; dir++)
        {
           if(dir != faceDir)
           {
             pfaceshift[index] = p[dir];
             index++;
           }
        }
        a_fluxMatrix(0, jmom) = a_faceMom[pfaceshift];
      }
      else
      {
        a_fluxMatrix(0, jmom) = 0;
      }
      if(Chombo4::LAPACKMatrix::s_verbose)
      {
        pout() << ",icol = "  << jmom <<  ", Q(0,icol) = " << a_fluxMatrix(0, jmom) << endl;
      }
    }
  }

  template<CENTERING cent, int order>
  Real
  DivWLSFlux<cent, order>::
  computeInhomogTermReg(const Chombo4::LAPACKMatrix     & a_stencil,
                        const Neighbors<order>          & a_neigh,
                        const int                         a_dir,
                        const EBDomainBC<DIM>    & a_domainBC,
                        Real                            a_dx) const
  {
    PR_TIME("computeInhomogTerm");
    using Chombo4::pout;
    Vector<IndMomSDMinOne> a_facmom = a_neigh.m_facMom;
    Vector<RealVect> a_facneighLoc = a_neigh.m_facLoc;
    Vector<int> a_facneighDir = a_neigh.m_facneighDir;
    Vector<Side::LoHiSide> a_facneighSid = a_neigh.m_facneighSid;
    Vector<Real> facneighFrac = a_neigh.m_facneighFrac;
    Real areaScaleFactor = getRegularArea(a_dx);

    Real inhomogTerm = 0.0;

    for(int iface = 0; iface < a_neigh.sizeFaces(); iface++)
    {
      Real bcValue = 0;
      RealVect faceLocDiff = a_facneighLoc[iface];
      int faceDir = a_facneighDir[iface];
      Side::LoHiSide faceSide = a_facneighSid[iface];
      if(m_domainBC.getBCtypes().getRegComp(faceDir, faceSide, a_dir) == string("Neumann"))
      {
        bcValue = NeumannExactSoln<order>::
                               getFaceIntegral(a_facmom[iface],
                                               faceLocDiff,
                                               faceDir,
                                               m_domainBC.getBCfuncs().getRegComp(faceDir, faceSide, a_dir),
                                               a_dx);
      }
      else if(m_domainBC.getBCtypes().getRegComp(faceDir, faceSide, a_dir) == string("Dirichlet"))
      {
        bcValue = DirichletExactSoln<order>::
                               getFaceIntegral(a_facmom[iface],
                                               faceLocDiff,
                                               faceDir,
                                               m_domainBC.getBCfuncs().getRegComp(faceDir, faceSide, a_dir),
                                               a_dx);
      }
      else
      {
        PROTO_ASSERT(false, "could not find domBC string");
      }

      if (facneighFrac[iface] < this->m_eps)
      {
        bcValue = 0.0;
      }
      else
      {
        bcValue /= (facneighFrac[iface]*areaScaleFactor);
      }

      inhomogTerm += bcValue*a_stencil(0, a_neigh.indexFromFace(iface));

      if (m_verbose)
        {
          pout() << "face loc" << faceLocDiff <<
            " boundary value " << bcValue << std::endl;
        }
    }
    return inhomogTerm;
  }

  template<CENTERING cent, int order>
  Real
  DivWLSFlux<cent, order>::
  computeInhomogTermEB(const Chombo4::LAPACKMatrix     & a_stencil,
                       const Neighbors<order>          & a_neigh,
                       const EBDomainBC<DIM>    & a_domainBC,
                       Real                            a_dx) const
  {
    PR_TIME("computeInhomogTermEB");
    using Chombo4::pout;
    Vector<IndMomDIM> a_ebfmom= a_neigh.m_ebfMom;
    Vector<Vector<IndMomDIM> > a_ebfneighNorm = a_neigh.m_ebfneighNorm;
    Vector<RealVect> a_ebfneighLoc = a_neigh.m_ebfLoc;
    Vector<Real> ebfneighFrac = a_neigh.m_ebfneighFrac;
    Real areaScaleFactor = getRegularArea(a_dx);

    Real inhomogTerm = 0.0;

    for(int iface = 0; iface < a_neigh.sizeEB(); iface++)
    {
      Real bcValue = 0;
      RealVect a_ebfLocDiff = a_ebfneighLoc[iface];
      if(m_domainBC.getBCtypes().getEBComp(0) == string("Neumann"))
      {
        // bcValue = NeumannExactSoln<order>::
        //                        getEBIntegral(a_ebfneighNorm[iface],
        //                                      a_ebfLocDiff,
        //                                      m_exactSolEB,
        //                                      a_dx);
      }
      else if(m_domainBC.getBCtypes().getEBComp(0) == string("Dirichlet"))
      {
        // bcValue is U dot n, sum up each component
        for (int dir=0; dir!=DIM; dir++)
          {
            //Normalizor<order>::normalizeFaceMoment(a_ebfneighNorm[dir][iface], a_dx);
            bcValue -= DirichletExactSoln<order>::
              getEBIntegral(a_ebfneighNorm[iface][dir],
                            a_ebfLocDiff,
                            m_domainBC.getBCfuncs().getEBComp(dir),
                            a_dx);
          }
      }
      else
      {
        PROTO_ASSERT(false, "could not find domBC string");
      }

      if (ebfneighFrac[iface] < this->m_eps)
      {
        bcValue = 0.0;
      }
      else
      {
        bcValue /= (ebfneighFrac[iface]*areaScaleFactor);
      }

      inhomogTerm += bcValue*a_stencil(0, a_neigh.indexFromEB(iface));
    }
    return inhomogTerm;
  }

  //Compute rhs matrix for a face
  template<int order>
  void
  GradWLSCell<order>::
  getEvalCellMatrix(Chombo4::LAPACKMatrix    & a_cellMatrix,
                    const EBIndex<CELL>      & a_cell,
                    RealVect                 & a_cellLoc,
                    IndMomDIM                & a_cellMom,
                    Real                     a_dx,
                    const IndexTM<int, DIM>    a_deriv) const
  {
    PR_TIME("getEvalCellMatrix");
    using Chombo4::pout;
    RealVect startLoc;
    for (int dir = 0; dir < DIM; dir++)
      {
        startLoc[dir] = a_dx * (0.5 + a_cell.m_pt[dir]);
      }
    // shift center of moments to a_cellLoc
    // the moments are centered at the cell center initially.
    // Move them to the centroid
    RealVect diff = startLoc - a_cellLoc;
    IndexTM<Real, DIM> diffrv;
    for(int dir = 0; dir < DIM; dir++)
      {
        diffrv[dir] = diff[dir]/a_dx;
      }
    a_cellMom.shift(diffrv);

    MomentIterator<DIM, order> momit;
    a_cellMatrix.define(1, momit.size());
    a_cellMatrix.setVal(0.);
    if (this->m_verbose) {pout() << " RHS - Cell Eval Matrix " << endl;}
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
  GradWLSCell<order>::
  getCellStencil(LocalStencil<CELL, Real> & a_sten,
                 Real                       & a_inhomogTerm, //for inhomogeneous bcs
                 const EBIndex<CELL>        & a_vof,
                 const EBGraph              & a_graph,
                 const MomentHolder<order>  & a_momentData,
                 IndMomDIM                  & a_cellMom,
                 Real                       a_dx,
                 const IndexTM<int, DIM>    a_deriv) const
  {
    PR_TIME("getCellStencil");
    using Chombo4::pout;
    Chombo4::LAPACKMatrix::s_verbose = this->m_verbose;
    // Chombo4::LAPACKMatrix::s_verbose = true;

    string a_ebbcName = m_domainBC.getBCtypes().getEBComp();
    string a_dombcName[2*DIM];
    m_domainBC.getBCtypes().fillFlatArray(a_dombcName);

    Neighbors<order> neighs;
    EBIndex<CELL> cell = a_vof; // a dummy, since we want cell bases neighbors
    RealVect cellLoc;
    int neighType=1; // about the cell
    for (int dir = 0; dir < DIM; dir++)
      {
        cellLoc[dir] = a_dx * (0.5 + cell.m_pt[dir]);
      }
    int radius = m_stenAttb.radius();
    Vector<EBIndex<CELL> > vols;
    MomentIterator<DIM,order> momit;
    int grows = 0;

    // REVIEW Make sure stencils are big enough. Size is not a sufficient measure, this should check the stencil is solvable
    while (neighs.size() <= momit.size()+DIM)
      {
        vols = neighs.template getVoFsNeighbors<CELL>(a_vof, a_vof, radius, a_graph, neighType);

        neighs.getMomentsInfoFromNeighbors(vols, a_graph, a_momentData, m_domainBC, a_dx);

        // if (grows>0) pout() << "neigh set tight for " << a_vof.m_pt << ", growing radius by 1" << endl;
        radius += 1;
        grows += 1;
      }
    CH_assert(vols.size() > momit.size());
    CH_assert(neighs.size() > momit.size());

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
    this->computeWLSStencilValues(phiCoef, Qshift, Mvol, weights);

    //Filling stencil
    a_sten.clear();
    Real coef_acc = 0.0;
    if (this->m_verbose) { pout() << " Cell Coefs Matrix " << endl;}
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
                                             m_domainBC,
                                             a_dx);
    if (this->m_verbose) {pout()<<"For phi=1 "<<coef_acc<<" "<<a_inhomogTerm<<" "<<coef_acc+a_inhomogTerm<<endl;}
  }

 //Compute inhomogeneous term
  template<int order>
  Real
  GradWLSCell<order>::
  computeInhomogTerm(const Chombo4::LAPACKMatrix    & a_stencil,
                     const Neighbors<order>         & a_neigh,
                     const int                        a_nvolneigh,
                     const EBDomainBC<1>     & a_domainBC,
                     Real                             a_dx) const
  {
    PR_TIME("computeInhomogTerm");
    using Chombo4::pout;
    Real inhomogTermReg = computeInhomogTermReg(a_stencil,
                                                a_neigh,
                                                a_nvolneigh,
                                                a_domainBC,
                                                a_dx);

    Real inhomogTermEB = computeInhomogTermEB(a_stencil,
                                                a_neigh,
                                                a_nvolneigh,
                                                a_domainBC,
                                                a_dx);

    return inhomogTermReg + inhomogTermEB;
  }


  template<int order>
  Real
  GradWLSCell<order>::
  computeInhomogTermReg(const Chombo4::LAPACKMatrix     & a_stencil,
                        const Neighbors<order>          & a_neigh,
                        const int                         a_nvolneigh,
                        const EBDomainBC<1>      & a_domainBC,
                        Real                            a_dx) const
  {
    PR_TIME("computeInhomogTerm");
    using Chombo4::pout;
    Vector<IndMomSDMinOne> a_facmom = a_neigh.m_facMom;
    Vector<RealVect> a_facneighLoc = a_neigh.m_facLoc;
    Vector<int> a_facneighDir = a_neigh.m_facneighDir;
    Vector<Side::LoHiSide> a_facneighSid = a_neigh.m_facneighSid;
    Vector<Real> facneighFrac = a_neigh.m_facneighFrac;
    Real areaScaleFactor = getRegularArea(a_dx);

    Real inhomogTerm = 0.0;

    for(int iface = 0; iface < a_neigh.sizeFaces(); iface++)
    {
      Real bcValue = 0;
      RealVect a_faceLocDiff = a_facneighLoc[iface];
      int faceDir = a_facneighDir[iface];
      Side::LoHiSide faceSide = a_facneighSid[iface];
      string bcName = m_domainBC.getBCtypes().getRegComp(faceDir, faceSide);
      if(bcName == string("Neumann"))
      {
        bcValue = NeumannExactSoln<order>::
                               getFaceIntegral(a_facmom[iface],
                                               a_faceLocDiff,
                                               faceDir,
                                               m_domainBC.getBCfuncs().getRegComp(faceDir, faceSide),
                                               a_dx);
      }
      else if(bcName == string("Dirichlet"))
      {
        bcValue = DirichletExactSoln<order>::
                               getFaceIntegral(a_facmom[iface],
                                               a_faceLocDiff,
                                               faceDir,
                                               m_domainBC.getBCfuncs().getRegComp(faceDir, faceSide),
                                               a_dx);
      }
      else
      {
        PROTO_ASSERT(false, "could not find domBC string");
      }

      if (facneighFrac[iface] < this->m_eps)
      {
        bcValue = 0.0;
      }
      else
      {
        bcValue /= (facneighFrac[iface]*areaScaleFactor);
      }

      inhomogTerm += bcValue*a_stencil(0, a_neigh.indexFromFace(iface));
    }

    return inhomogTerm;
  }

  template<int order>
  Real
  GradWLSCell<order>::
  computeInhomogTermEB(const Chombo4::LAPACKMatrix     & a_stencil,
                       const Neighbors<order>          & a_neigh,
                       const int                         a_nvolneigh,
                       const EBDomainBC<1>      & a_domainBC,
                       Real                            a_dx) const
  {
    PR_TIME("computeInhomogTermEB");
    using Chombo4::pout;
    Vector<IndMomDIM> a_ebfmom= a_neigh.m_ebfMom;
    Vector<Vector<IndMomDIM> > a_ebfneighNorm = a_neigh.m_ebfneighNorm;
    Vector<RealVect> a_ebfneighLoc = a_neigh.m_ebfLoc;
    Vector<Real> ebfneighFrac = a_neigh.m_ebfneighFrac;
    Real areaScaleFactor = getRegularArea(a_dx);

    Real inhomogTerm = 0.0;

    for(int iface = 0; iface < a_neigh.sizeEB(); iface++)
    {
      Real bcValue = 0;
      RealVect a_ebfLocDiff = a_ebfneighLoc[iface];
      string a_ebbcName = a_domainBC.getBCtypes().getEBComp();
      if(a_ebbcName == string("Neumann"))
      {
        bcValue = NeumannExactSoln<order>::
                               getEBIntegral(a_ebfneighNorm[iface],
                                             a_ebfLocDiff,
                                             a_domainBC.getBCfuncs().getEBComp(),
                                             a_dx);
      }
      else if(a_ebbcName == string("Dirichlet"))
      {
        bcValue = DirichletExactSoln<order>::
                               getEBIntegral(a_ebfmom[iface],
                                             a_ebfLocDiff,
                                             a_domainBC.getBCfuncs().getEBComp(),
                                             a_dx);
      }
      else
      {
        PROTO_ASSERT(false, "could not find domBC string");
      }

      if (ebfneighFrac[iface] < this->m_eps)
      {
        bcValue = 0.0;
      }
      else
      {
        bcValue /= (ebfneighFrac[iface]*areaScaleFactor);
      }

      inhomogTerm += bcValue*a_stencil(0, a_neigh.indexFromEB(iface));
    }

    return inhomogTerm;
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
        for(int dir = 0; dir < DIM; dir++)
        {
          if(dir  != a_faceDir)
          {
            areaMomIndex[index] =  momit()[dir];
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
        for(int dir = 0; dir < DIM; dir++)
        {
          if(dir  != a_faceDir)
          {
            areaMomIndex[index] =  momit()[dir];
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
    for(int dir = 0; dir < DIM; dir++)
    {
      // add in pd*int((x-x0)^p-ed*nd*dA)
      for(momit.reset(); momit.ok(); ++momit)
      {
        IndexTM<int, DIM> jmom = momit();
        if(jmom[dir] > 0)
        {
          a_ebfmom[jmom] += jmom[dir]*a_ebfNorm[dir][jmom-BASISV_TM<int,DIM>(dir)];
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
                  const shared_ptr<BaseExactSolution >  & a_exactSol,
                  const Real              a_dx)  
  {
    PR_TIME("getFaceIntegralNeumannExactSoln");
    Real retval = 0.0;
    Real areaScaleFactor = getRegularArea(a_dx);
    MomentIterator<DIM-1, order> momit;
    for(momit.reset(); momit.ok(); ++momit)
    {
      IndexTM<int,DIM-1>  pSDMinOne = momit();
      IndexTM<int,DIM>  pSpaceDim;
      // face moment is zero unless p in face direction is 1 (because we take a derivative in 
      // the face direction)
      int index = 0;
      for(int dir = 0; dir < DIM; dir++)
      {
        if(dir == a_faceDir)
        {
          pSpaceDim[dir] = 1;
        }
        else
        {
          pSpaceDim[dir] = pSDMinOne[index];
          index++;
        }
      }

      Real moment = a_faceMom[pSDMinOne];
 
      // pSpaceDim[a_face.direction()] is 1 so no need to multiply
      retval += moment*(areaScaleFactor*POW(a_dx,pSpaceDim))*a_exactSol->getDerivative(pSpaceDim,a_faceLoc)/(pSpaceDim.factorial());
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
                  const shared_ptr<BaseExactSolution >  & a_exactSol,
                  const Real              a_dx) 
  {
    PR_TIME("getFaceIntegralDirichletExactSoln");
    Real retval = 0.0;
    Real areaScaleFactor = getRegularArea(a_dx);
    MomentIterator<DIM-1, order> momit;
    for(momit.reset(); momit.ok(); ++momit)
    {
      IndexTM<int,DIM-1> pSDMinOne = momit();
      IndexTM<int,DIM> pSpaceDim;

      // face moment m_p is zero except when p[a_face.direction()] = 0 
      int index = 0;
      for(int dir = 0; dir < DIM; dir++)
      {
        if(dir == a_faceDir)
        {
          pSpaceDim[dir] = 0;
        }
        else
        {
          pSpaceDim[dir] = pSDMinOne[index];
          index++;
        }
      }

      Real moment = a_faceMom[pSDMinOne];

      retval += moment*(areaScaleFactor*POW(a_dx,pSpaceDim))*a_exactSol->getDerivative(pSpaceDim,a_faceLoc)/(pSpaceDim.factorial());

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
                const shared_ptr<BaseExactSolution >  & a_exactSol,
                const Real         a_dx) 
  {
    PR_TIME("getEBIntegralDirichletExactSoln");

    Real retval = 0.0;
    Real areaScaleFactor = getRegularArea(a_dx);

    MomentIterator<DIM, order> momit;
    for(momit.reset(); momit.ok(); ++momit)
    {
      IndexTM<int,DIM> p = momit();
      Real moment = a_ebfMom[p];

      retval += moment*(areaScaleFactor*POW(a_dx,p))*a_exactSol->getDerivative(p,a_ebfLoc)/(p.factorial());

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
                const shared_ptr<BaseExactSolution >  & a_exactSol,
                const Real                 a_dx) 
  {
    PR_TIME("getEBIntegralNeumannExactSoln");

    Real retval = 0.0;
    Real areaScaleFactor = getRegularArea(a_dx);
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
  
          retval += pd*moment*(areaScaleFactor*POW(a_dx,p))*a_exactSol->getDerivative(p,a_ebfLoc)/(p.factorial());
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
                    Real                     a_signNormal,
                    Real                     a_dx) const
  {
    PR_TIME("getFaceFluxMatrix");
    using Chombo4::pout;
    //Real areaScaleFactor  = getRegularArea(a_dx);
    int a_faceDir = a_face.direction();

    MomentIterator<DIM, order> momit;
    IndMomDIM facemomDIM;

    RealVect startLoc;
    for (int dir = 0; dir < DIM; dir++)
    {
      startLoc[dir] = a_dx * (0.5 + a_face.m_pt[dir]);
      if (a_face.direction() == dir) startLoc[dir] = a_dx * (a_face.m_pt[dir]);
    }
    // shift center of moments to a_faceLoc.
    // the moments are centered at the face center initially
    RealVect diff = startLoc;
    diff -= a_faceLoc;
    IndexTM<Real, DIM-1> diffrv;
    int index = 0;
    for(int dir = 0; dir < DIM; dir++)
    {
      if(dir != a_face.direction())
      {
        diffrv[index] = diff[dir]/a_dx;
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
        for(int dir = 0; dir < DIM; dir++)
        {
          if(dir  != a_faceDir)
          {
            areaMomIndex[index] =  momit()[dir];
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
    //const int& faceDir = a_face.direction();
    if (this->m_verbose) {pout() << " RHS - Adv Face Flux Matrix " << endl;}

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
                  Real                     a_dx) const
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
    using Chombo4::pout;
    Vector<int> weightDisVol = a_neigh.m_weightDisVol;
    Vector<int> weightDisFac = a_neigh.m_weightDisFac;
    Vector<int> weightDisEB = a_neigh.m_weightDisEB;
    Vector<RealVect> volneighLoc = a_neigh.m_volLoc;
    Vector<RealVect> ebfneighLoc = a_neigh.m_ebfLoc;
    Vector<IndMomDIM> volMom = a_neigh.m_volMom; 
    Vector<IndMomDIM> ebfMom = a_neigh.m_ebfMom; 
    RealVect v_faceLoc = m_exactVel->getVelocityAtPoint(a_faceLoc);
#if DIM==2
    Real v_faceLoc_mag = pow(pow(v_faceLoc[0],2)+pow(v_faceLoc[1],2),0.5);
#elif DIM==3
    Real v_faceLoc_mag = pow(pow(v_faceLoc[0],2)+pow(v_faceLoc[1],2)+pow(v_faceLoc[2],2),0.5);
#endif
    v_faceLoc[0] = v_faceLoc[0]/v_faceLoc_mag;
    v_faceLoc[1] = v_faceLoc[1]/v_faceLoc_mag;
#if DIM==3
    v_faceLoc[2] = v_faceLoc[2]/v_faceLoc_mag;
#endif
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

    //int powerE = this->m_powerE;
    //if (max_rad-1 == m_stenAttb.radius()) powerE = this->m_powerE;
    //else powerE = this->m_powerE - (2.0*(max_rad-m_stenAttb.radius()-1));
    if (this->m_verbose) { pout() << " Weights Matrix Adv " << endl;}
    for(int ineigh = 0; ineigh < nrows; ineigh++)
    {
      Real blerg = 0.;
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
        //int dist = weightDisFac[ineigh-weightDisVol.size()];
        //blerg = POW(dist,powerE);
        Chombo4::MayDay::Error("how did I get here? getWeightsMatrixFace");
      }
      if (ineigh>=weightDisVol.size()+weightDisFac.size())
      {
        //int dist = weightDisEB[ineigh-weightDisVol.size()-weightDisFac.size()];
        //blerg = POW(dist,powerE);
        Chombo4::MayDay::Error("how did I get here? getWeightsMatrixFace");
      }
      if (blerg == 0.)
        {
          Chombo4::MayDay::Error("how do I have zero weight? getWeightsMatrixFace");
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
  for(int dir = 2; dir < DIM; dir++)
    {
      retval[dir] = 0.0;
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

} // end namespace proto
#endif
