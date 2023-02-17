#ifndef _Proto_WLSStencilImplem_H_
#define _Proto_WLSStencilImplem_H_

namespace Proto
{
  // Get a single face flux stencil at a face by a distance weighted average between
  // two cell-centered flux stencils states computed by the neighbor cells sharing the face
  // This was originally designed for the Laplacian stencil, but is generally applicable
  template <int order>
  void
  WLSStencil<order>::
  getDistanceWeightedFluxAvgCoeffs(
    Real&    a_loSideWeight,
    Real&    a_hiSideWeight,
    const RealVect        a_faceCentroid,
    const EBIndex<CELL> & a_voflo,
    const EBIndex<CELL> & a_vofhi,
    const int             a_dir,
    const VoluData       & a_voludata,
    const EBGraph        & a_graph,
    const Real             a_dx)
  {
    Real expD = 1.0; // fixed power?
    // get the normalized centroids
    RealVect centroidM, centroidP;
    Real kappaM, kappaP;
    Normalizor<order>::normalizedVoluMoments(kappaM, centroidM, a_voludata, a_voflo, a_graph, a_dx);
    Normalizor<order>::normalizedVoluMoments(kappaP, centroidP, a_voludata, a_vofhi, a_graph, a_dx);
    // shift centroid to be relative to the shared face between them
    centroidM -= 0.5*BASISREALV(a_dir);
    centroidP += 0.5*BASISREALV(a_dir);
    // get the distance from cell centroid to face centroid
    Real distM = pow(((centroidM-a_faceCentroid)).vectorLength(), expD);
    Real distP = pow(((centroidP-a_faceCentroid)).vectorLength(), expD);
    // eval the weighted distances
    a_loSideWeight = distP/(distM+distP);
    a_hiSideWeight = distM/(distM+distP);
    // a_loSideWeight = distM/(distM+distP); // this does better?
    // a_hiSideWeight = distP/(distM+distP);
  }

    // Get a single face flux stencil at a face by a distance weighted average between
  // two cell-centered flux stencils states computed by the neighbor cells sharing the face
  // This was originally designed for the Laplacian stencil, but is generally applicable
  template <int order>
  LocalStencil<CELL,  Real>
  WLSStencil<order>::
  getDistanceWeightedFluxAvgStencil(
    LocalStencil< CELL,  Real> & a_fluxMStencil,
    LocalStencil< CELL,  Real> & a_fluxPStencil,
    const RealVect                 a_centroid,
    const EBIndex<CELL>          & a_voflo,
    const EBIndex<CELL>          & a_vofhi,
    const int                      a_dir,
    const VoluData               & a_voludata,
    const EBGraph                & a_graph,
    const Real                     a_dx)
  {
    LocalStencil< CELL,  Real> fluxStencil;

    Real weightM, weightP;
    getDistanceWeightedFluxAvgCoeffs(
      weightM,
      weightP,
      a_centroid,
      a_voflo,
      a_vofhi,
      a_dir,
      a_voludata,
      a_graph,
      a_dx);

    a_fluxMStencil *= weightM;
    a_fluxPStencil *= weightP;
    fluxStencil += a_fluxMStencil;
    fluxStencil += a_fluxPStencil;
    return fluxStencil;
  }

/**
 * This function does stuff
 */
  template < int order>
  Real
  WLSStencil<order>::
  getDistanceWeightedFluxAvgStencil(
    Real                           a_inhomogTermM,
    Real                           a_inhomogTermP,
    const RealVect                 a_centroid,
    const EBIndex<CELL>          & a_voflo,
    const EBIndex<CELL>          & a_vofhi,
    const int                      a_dir,
    const VoluData               & a_voludata,
    const EBGraph                & a_graph,
    const Real                     a_dx)
  {
    Real inhomogTerm;

    Real weightM, weightP;
    getDistanceWeightedFluxAvgCoeffs(
      weightM,
      weightP,
      a_centroid,
      a_voflo,
      a_vofhi,
      a_dir,
      a_voludata,
      a_graph,
      a_dx);

    a_inhomogTermM *= weightM;
    a_inhomogTermP *= weightP;
    inhomogTerm = a_inhomogTermM + a_inhomogTermP;
    return inhomogTerm;
  }


  //Gets Laplacian WLS flux divergence stencil
  //FIXME Short rant here.
  //      This class is bad and there are few things to fix. First every dimension is unlooped, and asking for mistakes. If 3D breaks look here first.
  //      Also the flux opjects really don't need to live outside of this function.
  //      So we pretty much do the same thing for every operator, with a little templating the sum over each flux part could just be done once and reused
  //      I don't really know why these functions all live in a class together.
  template < int order>
  void
  WLSStencil<order>::
  getLapWLSDivFStencil(LocalStencil< CELL,  Real>          & a_vofStencil,
                       Real                                & a_inhoTerm,
                       const EBIndex<CELL>                 & a_vof,
                       const EBGraph                       & a_graph,
                       const MomentHolder<order>           & a_momentData,
                       const LapWLSFlux<XFACE   , order>   & a_xfaceflux,
                       const LapWLSFlux<YFACE   , order>   & a_yfaceflux,
                       const LapWLSFlux<ZFACE   , order>   & a_zfaceflux,
                       const LapWLSFlux<BOUNDARY, order>   & a_ebflux,
                       int                                 a_neighType,
                       Real                                a_dx)
  {
    PR_TIME("getLaplaWLSDivFStencil");
    using Chombo4::pout;
    a_vofStencil.clear();
    Real areaScaleFactor = getRegularArea(a_dx);

    // unpack all the things, but ideally they just need passed
    auto& a_voludata = a_momentData.volumeMoments();
    auto& a_ebfadata = a_momentData.ebFaceMoments();
#if DIM >= 1
    auto& a_xfacdata = a_momentData.template faceMoments<XFACE>();
    auto& a_ebnormxdata = a_momentData.ebNormal(0);
#endif
#if DIM >= 2
    auto& a_yfacdata = a_momentData.template faceMoments<YFACE>();
    auto& a_ebnormydata = a_momentData.ebNormal(1);
#endif
#if DIM >= 3
    auto& a_zfacdata = a_momentData.template faceMoments<ZFACE>();
    auto& a_ebnormzdata = a_momentData.ebNormal(2);
#endif

    for (SideIterator sit; sit.ok(); ++sit)
    {
      vector< EBIndex<XFACE> > xfaces = a_graph.getXFaces(a_vof, sit());
      int isign = sign(sit());
      for (int iface = 0; iface < xfaces.size(); iface++)
      {
        LocalStencil<CELL, Real> fluxStencil;
        Real inhomogTerm=0.0; //for inhomogeneous bcs
        Real factor = Real(isign);
        std::map<EBIndex<XFACE>,LocalStencil<CELL, Real> >::iterator it;
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
          Real   areaFrac;
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

          if (a_xfaceflux.m_verbose) {pout()<<" stencil Xface "<<xfaces[iface].m_pt<<endl;}

          if (!xfaces[iface].m_isBoundary && a_neighType){
            LocalStencil<CELL, Real> fluxMStencil;
            LocalStencil<CELL, Real> fluxPStencil;
            Real inhomogTermM;
            Real inhomogTermP;
            EBIndex<CELL> voflo = xfaces[iface].getVoF(Side::Lo);
            EBIndex<CELL> vofhi = xfaces[iface].getVoF(Side::Hi);

            a_xfaceflux.getFaceCenteredStencil(fluxMStencil, inhomogTermM, voflo, xfaces[iface],
                                               a_graph,
                                               a_momentData,
                                               0,
                                               faceMom,
                                               a_neighType,
                                               factor,
                                               a_dx);

            a_xfaceflux.getFaceCenteredStencil(fluxPStencil, inhomogTermP, vofhi, xfaces[iface],
                                               a_graph,
                                               a_momentData,
                                               0,
                                               faceMom,
                                               a_neighType,
                                               factor,
                                               a_dx);

            fluxStencil = getDistanceWeightedFluxAvgStencil(fluxMStencil, fluxPStencil, centroid, voflo, vofhi, 0, a_voludata, a_graph, a_dx);
            inhomogTerm = getDistanceWeightedFluxAvgStencil(inhomogTermM, inhomogTermP, centroid, voflo, vofhi, 0, a_voludata, a_graph, a_dx);


          }else{
            a_xfaceflux.getFaceCenteredStencil(fluxStencil, inhomogTerm, a_vof, xfaces[iface],
                                               a_graph,
                                               a_momentData,
                                               0,
                                               faceMom,
                                               a_neighType,
                                               factor,
                                               a_dx);
          }

          if (a_xfaceflux.m_verbose)
          {
            pout()<<" inhomogTerm Xface "<<inhomogTerm<<endl;
            pout()<<" lapfluxStencil "<<endl;
            fluxStencil.print();
            fluxStencil.printSum();
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
        LocalStencil<CELL, Real> fluxStencil;
        Real inhomogTerm = 0.0; //for inhomogeneous bcs
        Real factor = Real(isign);
        std::map<EBIndex<YFACE>,LocalStencil<CELL, Real> >::iterator it;
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
          Real   areaFrac;
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

          if (a_yfaceflux.m_verbose) {pout()<<" stencilY face "<<yfaces[iface].m_pt<<endl;}

          if (!yfaces[iface].m_isBoundary && a_neighType){
            LocalStencil<CELL, Real> fluxMStencil;
            LocalStencil<CELL, Real> fluxPStencil;
            Real inhomogTermM;
            Real inhomogTermP;
            EBIndex<CELL> voflo = yfaces[iface].getVoF(Side::Lo);
            EBIndex<CELL> vofhi = yfaces[iface].getVoF(Side::Hi);

            a_yfaceflux.getFaceCenteredStencil(fluxMStencil, inhomogTermM, voflo, yfaces[iface],
                                           a_graph,
                                           a_momentData,
                                           1,
                                           faceMom,
                                           a_neighType,
                                           factor,
                                           a_dx);

            a_yfaceflux.getFaceCenteredStencil(fluxPStencil, inhomogTermP, vofhi, yfaces[iface],
                                           a_graph,
                                           a_momentData,
                                           1,
                                           faceMom,
                                           a_neighType,
                                           factor,
                                           a_dx);
            fluxStencil = getDistanceWeightedFluxAvgStencil(fluxMStencil, fluxPStencil, centroid, voflo, vofhi, 1, a_voludata, a_graph, a_dx);
            inhomogTerm = getDistanceWeightedFluxAvgStencil(inhomogTermM, inhomogTermP, centroid, voflo, vofhi, 1, a_voludata, a_graph, a_dx);
            //if (a_yfaceflux.m_verbose)
            //{
            //  pout()<<" yflux"<<endl;
            //  pout()<<" weight for M "<<distP/(distM+distP)<<endl;
            //  pout()<<" weight for P "<<distM/(distM+distP)<<endl;
            //}
          }else{
            a_yfaceflux.getFaceCenteredStencil(fluxStencil, inhomogTerm, a_vof, yfaces[iface],
                                               a_graph,
                                               a_momentData,
                                               1,
                                               faceMom,
                                               a_neighType,
                                               factor,
                                               a_dx);
          }

          if (a_yfaceflux.m_verbose)
          {
            pout()<<" inhomogTerm Yface "<<inhomogTerm<<endl;
            pout()<<" lapfluxStencil "<<endl;
            fluxStencil.print();
            fluxStencil.printSum();
          }
          yfacesFluxStencil[yfaces[iface]] = fluxStencil;
          yfacesInhTerm[yfaces[iface]] = inhomogTerm;
          fluxStencil *= factor;
          a_inhoTerm += inhomogTerm*factor;
          a_vofStencil += fluxStencil;
        }
      }
    }
#if DIM>=3
    for (SideIterator sit; sit.ok(); ++sit)
    {
      vector< EBIndex<ZFACE> > zfaces = a_graph.getZFaces(a_vof, sit());
      int isign = sign(sit());
      for (int iface = 0; iface < zfaces.size(); iface++)
      {
        LocalStencil<CELL, Real> fluxStencil;
        Real inhomogTerm; //for inhomogeneous bcs
        Real factor = Real(isign);
        std::map<EBIndex<ZFACE>,LocalStencil<CELL, Real> >::iterator it;
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
          Real   areaFrac;
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
            LocalStencil<CELL, Real> fluxMStencil;
            LocalStencil<CELL, Real> fluxPStencil;
            Real inhomogTermM;
            Real inhomogTermP;
            EBIndex<CELL> voflo = zfaces[iface].getVoF(Side::Lo);
            EBIndex<CELL> vofhi = zfaces[iface].getVoF(Side::Hi);

            a_zfaceflux.getFaceCenteredStencil(fluxMStencil, inhomogTermM, voflo, zfaces[iface],
                                           a_graph,
                                           a_momentData,
                                           2,
                                           faceMom,
                                           a_neighType,
                                           factor,
                                           a_dx);

            a_zfaceflux.getFaceCenteredStencil(fluxPStencil, inhomogTermP, vofhi, zfaces[iface],
                                           a_graph,
                                           a_momentData,
                                           2,
                                           faceMom,
                                           a_neighType,
                                           factor,
                                           a_dx);
            fluxStencil = getDistanceWeightedFluxAvgStencil(fluxMStencil, fluxPStencil, centroid, voflo, vofhi, 2, a_voludata, a_graph, a_dx);
            inhomogTerm = getDistanceWeightedFluxAvgStencil(inhomogTermM, inhomogTermP, centroid, voflo, vofhi, 2, a_voludata, a_graph, a_dx);
          }else{
            a_zfaceflux.getFaceCenteredStencil(fluxStencil, inhomogTerm, a_vof, zfaces[iface],
                                               a_graph,
                                               a_momentData,
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

    LocalStencil<CELL, Real> ebflux;
    Real alpha;
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
#if DIM>=3
      ebnorm[2] = a_ebnormzdata(a_vof, 0);
      for(MomentIterator<DIM, order> iter; iter.ok(); ++iter)
      {
        const IndexTM<int, DIM>& p = iter();
        ebnorm[2][iter()] = ebnorm[2][iter()]/(areaScaleFactor*POW(a_dx,p));
      }
#endif
      Normalizor<order>::normalizedBndryMoments(alpha, centroid, a_ebfadata,  cutface, a_graph, a_dx);

      if (a_ebflux.m_verbose) {pout()<<" stencil ebface "<<a_vof.m_pt<<endl;}

      a_ebflux.getEBFluxStencil(ebflux, inhomogTerm, a_vof,
                                cutface,
                                a_graph,
                                a_momentData,
                                ebnorm,
                                a_neighType,
                                a_dx);

      if (a_ebflux.m_verbose)
      {
        pout() << " inhomogTerm ebface " << inhomogTerm << endl;
        pout() << " alpha " << alpha << endl;
        pout() << " ebfluxStencil " << endl;
        ebflux.print();
      }
      a_inhoTerm += inhomogTerm;
      a_vofStencil += ebflux;
    }

    //for some reason, EBAMRPoissonOp calculates beta weight without the EB flux
    Real diagweight;
    bool foundVoF = a_vofStencil.getValueAtVoF(diagweight, a_vof);
    PR_assert(foundVoF); //if this fails, the stencil has no diagonal value
    a_vofStencil.setBetaWeight(diagweight);

    // offset the stencil to respect constants
    // Real stenSum, inhoSum;
    // for(int ivec =  0; ivec < a_vofStencil.size(); ivec++)
    //   {
    //     stenSum += a_vofStencil.m_entries[ivec].m_weight;
    //   }
    // stenSum /= a_vofStencil.size();
    // pout() << "Row sum adj = " << stenSum << std::endl;
    // for(int ivec =  0; ivec < a_vofStencil.size(); ivec++)
    //   {
    //     a_vofStencil.m_entries[ivec].m_weight -= stenSum;
    //   }

  }


  //Getting flux stencil for advection operator at a face by a weighted average between
  //two cell-centered flux stencils states computed by the neighbor cells sharing the face
  template < int order>
  LocalStencil< CELL,  Real>
  WLSStencil<order>::
  getAdvFluxStencil(LocalStencil< CELL,  Real> & a_fluxMStencil,
                    LocalStencil< CELL,  Real> & a_fluxPStencil,
                    RealVect                     & a_centroid,
                    const EBIndex<CELL>          & a_voflo,
                    const EBIndex<CELL>          & a_vofhi,
                    const RealVect               & a_v_faceLoc,
                    int                            a_dir,
                    const VoluData               & a_voludata,
                    const EBGraph                & a_graph,
                    Real                         a_dx)
  {
    LocalStencil< CELL,  Real> fluxStencil;
    RealVect centroidM;
    Real kappaM;
    Normalizor<order>::normalizedVoluMoments(kappaM, centroidM, a_voludata, a_voflo, a_graph, a_dx);
    RealVect centroidP;
    Real kappaP;
    Normalizor<order>::normalizedVoluMoments(kappaP, centroidP, a_voludata, a_vofhi, a_graph, a_dx);
    Real expD = 1.0;
    Real distM = 1.0;
    Real distP = 1.0;
#if DIM>=3
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
  getAdvWLSDivFStencil(LocalStencil< CELL,  Real>          & a_vofStencil,
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
                         Real                                  a_dx)
  {
    PR_TIME("getAdvWLSDivFStencil");
    using Chombo4::pout;
    a_vofStencil.clear();
    //Real volScaleFactor  = (D_TERM(a_dx, *a_dx, *a_dx));
    Real areaScaleFactor = getRegularArea(a_dx);

    for (SideIterator sit; sit.ok(); ++sit)
    {
      vector< EBIndex<XFACE> > xfaces = a_graph.getXFaces(a_vof, sit());
      int isign = sign(sit());
      for (int iface = 0; iface < xfaces.size(); iface++)
      {
        LocalStencil<CELL, Real> fluxMStencil;
        LocalStencil<CELL, Real> fluxPStencil;
        LocalStencil<CELL, Real> fluxStencil;
        Real inhomogTerm=0.0; //for inhomogeneous bcs
        Real factor = Real(isign);
        Real   areaFrac;
        RealVect centroid;
        Normalizor<order>::normalizedFaceMoments(areaFrac, centroid, a_xfacdata, xfaces[iface], a_graph, a_dx);
        std::map<EBIndex<XFACE>,LocalStencil<CELL, Real> >::iterator it;
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

          if (a_xfaceflux.m_verbose) {pout()<<" stencil Xface "<<xfaces[iface].m_pt<<" centroid "<<(a_dx*centroid)+faceLoc<<endl;}

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
        LocalStencil<CELL, Real> fluxMStencil;
        LocalStencil<CELL, Real> fluxPStencil;
        LocalStencil<CELL, Real> fluxStencil;
        Real inhomogTerm=0.0; //for inhomogeneous bcs
        Real factor = Real(isign);
        Real   areaFrac;
        RealVect centroid;
        Normalizor<order>::normalizedFaceMoments(areaFrac, centroid, a_yfacdata, yfaces[iface], a_graph, a_dx);
        std::map<EBIndex<YFACE>,LocalStencil<CELL, Real> >::iterator it;
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

          if (a_yfaceflux.m_verbose) {pout()<<" stencil Yface "<<yfaces[iface].m_pt<<" centroid "<<(a_dx*centroid)+faceLoc<<endl;}

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
#if DIM>=3
    for (SideIterator sit; sit.ok(); ++sit)
    {
      vector< EBIndex<ZFACE> > zfaces = a_graph.getZFaces(a_vof, sit());
      int isign = sign(sit());
      for (int iface = 0; iface < zfaces.size(); iface++)
      {
        LocalStencil<CELL, Real> fluxMStencil;
        LocalStencil<CELL, Real> fluxPStencil;
        LocalStencil<CELL, Real> fluxStencil;
        Real inhomogTerm=0.0; //for inhomogeneous bcs
        Real factor = Real(isign);
        std::map<EBIndex<ZFACE>,LocalStencil<CELL, Real> >::iterator it;
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
          Real   areaFrac;
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

          if (a_zfaceflux.m_verbose) {pout()<<" stencil Zface "<<zfaces[iface].m_pt<<" centroid "<<(a_dx*centroid)+faceLoc<<endl;}

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
    Real diagweight;
    bool foundVoF = a_vofStencil.getValueAtVoF(diagweight, a_vof);
    PR_assert(foundVoF); //if this fails, the stencil has no diagonal value
    a_vofStencil.setBetaWeight(diagweight);
  }


  // template <int order>
  // template<int dir>
  // void
  // WLSStencil<order>::
  // getWLSDivFCompStencilDir(
  //   std::array<LocalStencil<CELL, Real>, SpaceDim >& a_vofStencil,
  //   const EBIndex<CELL>                 & a_vof,
  //   const EBGraph                       & a_graph,
  //   const MomentHolder<order>           & a_momentData,
  //   const DivWLSFlux<XFACE   , order>  & a_xfaceflux, // Why build these separate?
  //   const DivWLSFlux<YFACE   , order>  & a_yfaceflux,
  //   const DivWLSFlux<ZFACE   , order>  & a_zfaceflux,
  //   const DivWLSFlux<BOUNDARY, order>  & a_ebflux,
  //   string                                a_ebbcName,
  //   vector<string>                        a_dombcName[2*DIM],
  //   Real                                & a_inhoTerm,
  //   int                                 a_neighType,
  //   Real                                a_dx)
  // {
  //   PR_TIME("getWLSDivFCompStencil");
  //   using Chombo4::pout;

  //   CH_assert((dir >= 0) && (dir < DIM));
  //   constexpr CENTERING FACE = REG_FACES[dir];

  //   auto& voludata = a_momentData.volumeMoments();
  //   auto& ebfadata = a_momentData.ebFaceMoments();
  //   auto& dirFacData = a_momentData.template faceMoments<FACE>();
  //   auto& ebnormdata = a_momentData.ebNormal(dir);

  //   // auto& a_faceflux = D_SELECT(a_xfaceflux, a_yfaceflux, a_zfaceflux);
  //   auto& a_faceflux = D_SELECT(a_xfaceflux, a_yfaceflux, a_zfaceflux);

  //   auto& vofStencil = a_vofStencil[dir];
  //   vofStencil.clear();
  //   Real areaScaleFactor = getRegularArea(a_dx);

  //   for (SideIterator sit; sit.ok(); ++sit)
  //   {
  //     // vector< EBIndex<FACE> > faces = a_graph.getFaces<FACE>(a_vof, sit());
  //     vector< EBIndex<FACE> > faces = getFaces<FACE>(a_vof, sit(), a_graph);
  //     //vector< EBIndex<XFACE> > xfaces = a_graph.getXFaces(a_vof, sit());
  //     int isign = sign(sit());
  //     for (int iface = 0; iface < faces.size(); iface++)
  //     {
  //       auto thisFace = faces[iface];
  //       LocalStencil<CELL, Real> fluxMStencil;
  //       LocalStencil<CELL, Real> fluxPStencil;
  //       LocalStencil<CELL, Real> fluxStencil;
  //       Real inhomogTerm = 0.; //for inhomogeneous bcs
  //       Real factor = Real(isign);

  //       // solve new stencil
  //       {
  //         EBIndex<CELL> voflo = faces[iface].getVoF(Side::Lo);
  //         EBIndex<CELL> vofhi = faces[iface].getVoF(Side::Hi);
  //         Real   areaFrac;
  //         RealVect centroid;
  //         Normalizor<order>::normalizedFaceMoments(areaFrac, centroid, dirFacData, faces[iface], a_graph, a_dx);
  //         RealVect faceLoc;
  //         for (int idir = 0; idir < DIM; idir++)
  //         {
  //           faceLoc[idir] = a_dx * (0.5 + faces[iface].m_pt[idir]);
  //           if (idir == 0) faceLoc[idir] = a_dx * (faces[iface].m_pt[idir]);
  //         }

  //         IndMomSDMinOne faceMom;
  //         if (areaFrac < 1.0)
  //           {
  //             faceMom = dirFacData(faces[iface], 0);
  //             for(MomentIterator<DIM-1, order> iter; iter.ok(); ++iter)
  //               {
  //                 const IndexTM<int, DIM-1>& p = iter();
  //                 faceMom[iter()] = faceMom[iter()]/(areaScaleFactor*POW(a_dx,p));
  //               }
  //           }
  //         else
  //           {
  //             faceMom.setToRegular(a_dx);
  //             for(MomentIterator<DIM-1, order> iter; iter.ok(); ++iter)
  //               {
  //                 const IndexTM<int, DIM-1>& p = iter();
  //                 faceMom[iter()] = faceMom[iter()]/(areaScaleFactor*POW(a_dx,p));
  //               }
  //           }

  //         if (a_faceflux.m_verbose)
  //           {
  //             pout() << std::endl <<" stencil Face "<<faces[iface].m_pt<<endl;
  //           }

  //         string xdombcName[2*DIM];
  //         for(int ivec = 0; ivec < 2*DIM; ivec++)
  //           {
  //             xdombcName[ivec] = a_dombcName[ivec][0];
  //           }
  //         if (!faces[iface].m_isBoundary){
  //           Real inhomogTermM = 0;
  //           Real inhomogTermP = 0;
  //           // the low side stencil (M)inus
  //           a_faceflux.getFaceCenteredStencil(fluxMStencil, inhomogTermM, voflo, faces[iface],
  //                                              a_graph,
  //                                              a_momentData,
  //                                              a_ebbcName,
  //                                              xdombcName,
  //                                              0,
  //                                              faceMom,
  //                                              a_neighType,
  //                                              factor,
  //                                              a_dx);
  //           // the high side stencil (P)lus
  //           a_faceflux.getFaceCenteredStencil(fluxPStencil, inhomogTermP, vofhi, faces[iface],
  //                                              a_graph,
  //                                              a_momentData,
  //                                              a_ebbcName,
  //                                              xdombcName,
  //                                              0,
  //                                              faceMom,
  //                                              a_neighType,
  //                                              factor,
  //                                              a_dx);

  //           // Really just get the face average
  //           fluxStencil = getDistanceWeightedFluxAvgStencil(fluxMStencil, fluxPStencil, centroid, voflo, vofhi, 0, voludata, a_graph, a_dx);
  //           inhomogTerm = getDistanceWeightedFluxAvgStencil(inhomogTermM, inhomogTermP, centroid, voflo, vofhi, 0, voludata, a_graph, a_dx);
  //         }else{
  //           a_faceflux.getFaceCenteredStencil(fluxStencil, inhomogTerm, a_vof, faces[iface],
  //                                              a_graph,
  //                                              a_momentData,
  //                                              a_ebbcName,
  //                                              xdombcName,
  //                                              0,
  //                                              faceMom,
  //                                              a_neighType,
  //                                              factor,
  //                                              a_dx);
  //         }
  //         if (a_faceflux.m_verbose)
  //           {
  //             pout()<<" inhomogTerm Face "<<inhomogTerm<<endl;
  //             pout()<<" divFluxStencil "<<endl;
  //             fluxStencil.print();
  //             fluxStencil.printSum();
  //           }
  //         // facesFluxStencil[faces[iface]] = fluxStencil;
  //         // facesInhTerm[faces[iface]] = inhomogTerm;
  //         Real factor = Real(isign);
  //         fluxStencil *= factor;
  //         a_inhoTerm += inhomogTerm*factor;
  //         a_vofStencil[dir] += fluxStencil;
  //       }
  //     }
  //   }
  // }

  //Gets WLS flux for the divergence stencil
  template <int order>
  void
  WLSStencil<order>::
  getWLSDivFCompStencil(std::array<LocalStencil<CELL, Real>, SpaceDim >& a_vofStencil,
                        Real                               & a_inhoTerm,
                        const EBIndex<CELL>                   & a_vof,
                        const EBGraph                         & a_graph,
                        const MomentHolder<order>             & a_momentData,
                        const DivWLSFlux<XFACE, order>  & a_xfaceflux,
                        const DivWLSFlux<YFACE, order>  & a_yfaceflux,
                        const DivWLSFlux<ZFACE, order>  & a_zfaceflux,
                        const DivWLSFlux<BOUNDARY, order>  & a_ebflux,
                        int                                a_neighType,
                        Real                                  a_dx)
  {
    PR_TIME("getWLSDivFCompStencil");
    using Chombo4::pout;

    // One day this may work and I can "loop" over each direction, but putting on hold for now
    // D_CALLS(getWLSDivFCompStencilDir,
    //         a_vofStencil,
    //         a_vof,
    //         a_graph,
    //         a_momentData,
    //         // a_xfaceflux,
    //         // a_yfaceflux,
    //         // a_zfaceflux,
    //         // a_ebflux,
    //         a_ebbcName,
    //         a_dombcName,
    //         a_inhoTerm,
    //         a_neighType,
    //         a_dx);

    auto& a_voludata = a_momentData.volumeMoments();
    //auto& a_ebfadata = a_momentData.ebFaceMoments();
#if DIM >= 1
    auto& a_xfacdata = a_momentData.template faceMoments<XFACE>();
    //auto& a_ebnormxdata = a_momentData.ebNormal(0);
#endif
#if DIM >= 2
    auto& a_yfacdata = a_momentData.template faceMoments<YFACE>();
    //auto& a_ebnormydata = a_momentData.ebNormal(1);
#endif
#if DIM >= 3
    auto& a_zfacdata = a_momentData.template faceMoments<ZFACE>();
    //auto& a_ebnormzdata = a_momentData.ebNormal(2);
#endif

    // a_vofStencil.clear();
    auto& a_vofStencilx = a_vofStencil[0];
    a_vofStencilx.clear();
#if DIM>=2
    auto& a_vofStencily = a_vofStencil[1];
    a_vofStencily.clear();
#elif DIM>=3
    auto& a_vofStencilz = a_vofStencil[2];
    a_vofStencilz.clear();
#endif
    a_inhoTerm = 0.;
    Real areaScaleFactor = getRegularArea(a_dx);

    /// First do all the x-face stuff
    for (SideIterator sit; sit.ok(); ++sit)
    {
      //vector< EBIndex<FACE> > xfaces = a_graph.getFaces<FACE>(a_vof, sit());
      vector< EBIndex<XFACE> > xfaces = a_graph.getXFaces(a_vof, sit());
      int isign = sign(sit());
      for (int iface = 0; iface < xfaces.size(); iface++)
      {
        std::array<LocalStencil<CELL, Real>, SpaceDim > fullFluxStencilM, fullFluxStencilP, fullFluxStencil;

        auto thisFace = xfaces[iface];
        LocalStencil<CELL, Real> fluxMStencil;
        LocalStencil<CELL, Real> fluxPStencil;
        LocalStencil<CELL, Real> fluxStencil;
        Real inhomogTerm = 0.; //for inhomogeneous bcs
        Real factor = Real(isign);
        // Cache solved stencils
        // std::map<EBIndex<XFACE>,LocalStencil<CELL, Real> >::iterator it;
        // it = xfacesFluxStencil.find(xfaces[iface]);
        // if (it != xfacesFluxStencil.end())
        // {
        //   fluxStencil = it->second;
        //   inhomogTerm = xfacesInhTerm[xfaces[iface]];
        //   fluxStencil *= factor;
        //   a_inhoTerm += inhomogTerm*factor;
        //   a_vofStencilx += fluxStencil;
        // }
        // // solve new stencil
        // else
        {
          EBIndex<CELL> voflo = xfaces[iface].getVoF(Side::Lo);
          EBIndex<CELL> vofhi = xfaces[iface].getVoF(Side::Hi);
          Real   areaFrac;
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
            }
          else
            {
              faceMom.setToRegular(a_dx);
            }
          Normalizor<order>::normalizeFaceMoment(faceMom, a_dx);

          if (a_xfaceflux.m_verbose)
            {
              pout() << std::endl <<" stencil Xface "<<xfaces[iface].m_pt<<endl;
            }

          if (!xfaces[iface].m_isBoundary){
            Real inhomogTermM = 0;
            Real inhomogTermP = 0;
            // the low side stencil (M)inus
            a_xfaceflux.getFaceCenteredStencilFull(fullFluxStencilM,
                                               inhomogTermM, voflo, xfaces[iface],
                                               a_graph,
                                               a_momentData,
                                               0,
                                               faceMom,
                                               a_neighType,
                                               factor,
                                               a_dx);
            // the high side stencil (P)lus
            a_xfaceflux.getFaceCenteredStencilFull(fullFluxStencilP,
                                               inhomogTermP, vofhi, xfaces[iface],
                                               a_graph,
                                               a_momentData,
                                               0,
                                               faceMom,
                                               a_neighType,
                                               factor,
                                               a_dx);

            // average the fluxes on the face
            for (int dir = 0; dir!=DIM; dir++)
            {
              fullFluxStencil[dir] = getDistanceWeightedFluxAvgStencil(fullFluxStencilM[dir], fullFluxStencilP[dir], centroid, voflo, vofhi, 0, a_voludata, a_graph, a_dx);
            }
            inhomogTerm = getDistanceWeightedFluxAvgStencil(inhomogTermM, inhomogTermP, centroid, voflo, vofhi, 0, a_voludata, a_graph, a_dx);
          }else{
            a_xfaceflux.getFaceCenteredStencilFull(fullFluxStencil,
                                               inhomogTerm, a_vof, xfaces[iface],
                                               a_graph,
                                               a_momentData,
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
              fluxStencil.printSum();
            }
          // xfacesFluxStencil[xfaces[iface]] = fluxStencil;
          // xfacesInhTerm[xfaces[iface]] = inhomogTerm;
          Real factor = Real(isign);
          fluxStencil *= factor;
          a_inhoTerm += inhomogTerm*factor;
          for (int dir = 0; dir!=DIM; dir++)
            {
              fullFluxStencil[dir] *= factor;
              a_vofStencil[dir] += fullFluxStencil[dir];
            }
        }
      }
    }

    for (SideIterator sit; sit.ok(); ++sit)
    {
      vector< EBIndex<YFACE> > yfaces = a_graph.getYFaces(a_vof, sit());
      int isign = sign(sit());
      for (int iface = 0; iface < yfaces.size(); iface++)
      {
        std::array<LocalStencil<CELL, Real>, SpaceDim > fullFluxStencilM, fullFluxStencilP, fullFluxStencil;

        LocalStencil<CELL, Real> fluxMStencil;
        LocalStencil<CELL, Real> fluxPStencil;
        LocalStencil<CELL, Real> fluxStencil;
        Real inhomogTerm=0.0; //for inhomogeneous bcs
        Real factor = Real(isign);
        // std::map<EBIndex<YFACE>,LocalStencil<CELL, Real> >::iterator it;
        // it = yfacesFluxStencil.find(yfaces[iface]);
        // if (it != yfacesFluxStencil.end())
        // {
        //   fluxStencil = it->second;
        //   inhomogTerm = yfacesInhTerm[yfaces[iface]];
        //   fluxStencil *= factor;
        //   a_inhoTerm += inhomogTerm*factor;
        //   a_vofStencily += fluxStencil;
        // }
        // else
        {
          EBIndex<CELL> voflo = yfaces[iface].getVoF(Side::Lo);
          EBIndex<CELL> vofhi = yfaces[iface].getVoF(Side::Hi);
          Real   areaFrac;
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

          if (a_yfaceflux.m_verbose) {pout()<<" stencilY face "<<yfaces[iface].m_pt<<endl;}

          if (!yfaces[iface].m_isBoundary){
            Real inhomogTermM = 0;
            Real inhomogTermP = 0;
            a_yfaceflux.getFaceCenteredStencilFull(fullFluxStencilM,
                                               inhomogTermM, voflo, yfaces[iface],
                                               a_graph,
                                               a_momentData,
                                               1,
                                               faceMom,
                                               a_neighType,
                                               factor,
                                               a_dx);

            a_yfaceflux.getFaceCenteredStencilFull(fullFluxStencilP,
                                               inhomogTermP, vofhi, yfaces[iface],
                                               a_graph,
                                               a_momentData,
                                               1,
                                               faceMom,
                                               a_neighType,
                                               factor,
                                               a_dx);
            // fluxStencil = getAdvFluxStencil(fluxMStencil, fluxPStencil, centroid, voflo, vofhi, v_faceLoc, 1, a_voludata, a_graph, a_dx);
                        // Really just the face average
            // fluxStencil = getDistanceWeightedFluxAvgStencil(fluxMStencil, fluxPStencil, centroid, voflo, vofhi, 1, a_voludata, a_graph, a_dx);
            for (int dir = 0; dir!=DIM; dir++)
              {
                fullFluxStencil[dir] = getDistanceWeightedFluxAvgStencil(fullFluxStencilM[dir], fullFluxStencilP[dir], centroid, voflo, vofhi, 0, a_voludata, a_graph, a_dx);
              }
            inhomogTerm = getDistanceWeightedFluxAvgStencil(inhomogTermM, inhomogTermP, centroid, voflo, vofhi, 1, a_voludata, a_graph, a_dx);
          }
          else
          {
            a_yfaceflux.getFaceCenteredStencilFull(fullFluxStencil,
                                               inhomogTerm, a_vof, yfaces[iface],
                                               a_graph,
                                               a_momentData,
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
            fluxStencil.printSum();
          }
          // yfacesFluxStencil[yfaces[iface]] = fluxStencil;
          // yfacesInhTerm[yfaces[iface]] = inhomogTerm;
          Real factor = Real(isign);
          fluxStencil *= factor;
          a_inhoTerm += inhomogTerm*factor;
          for (int dir = 0; dir!=DIM; dir++)
            {
              fullFluxStencil[dir] *= factor;
              a_vofStencil[dir] += fullFluxStencil[dir];
            }
        }
      }
    }
#if DIM>=3
    for (SideIterator sit; sit.ok(); ++sit)
    {
      vector< EBIndex<ZFACE> > zfaces = a_graph.getZFaces(a_vof, sit());
      int isign = sign(sit());
      for (int iface = 0; iface < zfaces.size(); iface++)
      {
        LocalStencil<CELL, Real> fluxMStencil;
        LocalStencil<CELL, Real> fluxPStencil;
        LocalStencil<CELL, Real> fluxStencil;
        Real inhomogTerm=0.0; //for inhomogeneous bcs
        Real factor = Real(isign);
        std::map<EBIndex<ZFACE>,LocalStencil<CELL, Real> >::iterator it;
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
          Real   areaFrac;
          RealVect centroid;
          Normalizor<order>::normalizedFaceMoments(areaFrac, centroid, a_zfacdata, zfaces[iface], a_graph, a_dx);
          RealVect faceLoc;
          for (int idir = 0; idir < DIM; idir++)
          {
            faceLoc[idir] = a_dx * (0.5 + zfaces[iface].m_pt[idir]);
            if (idir == 2) faceLoc[idir] = a_dx * (zfaces[iface].m_pt[idir]);
          }

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

          if (a_zfaceflux.m_verbose)
            {
              pout() << std::endl <<" stencil Zface "<< zfaces[iface].m_pt << endl;
            }

          string zdombcName[2*DIM];
          for(int ivec = 0; ivec < 2*DIM; ivec++)
            {
              zdombcName[ivec] = a_dombcName[ivec][2];
            }
          if (!zfaces[iface].m_isBoundary){
            Real inhomogTermM = 0;
            Real inhomogTermP = 0;
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
                                               zdombcName,
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
                                               zdombcName,
                                               2,
                                               faceMom,
                                               a_neighType,
                                               factor,
                                               a_dx);

            // Really just the face average
            fluxStencil = getDistanceWeightedFluxAvgStencil(fluxMStencil, fluxPStencil, centroid, voflo, vofhi, 2, a_voludata, a_graph, a_dx);
            inhomogTerm = getDistanceWeightedFluxAvgStencil(inhomogTermM, inhomogTermP, centroid, voflo, vofhi, 2, a_voludata, a_graph, a_dx);
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
                                               zdombcName,
                                               2,
                                               faceMom,
                                               a_neighType,
                                               factor,
                                               a_dx);

          }
          if (a_zfaceflux.m_verbose)
            {
              pout()<<" inhomogTerm Zface "<<inhomogTerm<<endl;
              pout()<<" divFluxStencil "<<endl;
              fluxStencil.print();
              fluxStencil.printSum();
            }
          zfacesFluxStencil[zfaces[iface]] = fluxStencil;
          zfacesInhTerm[zfaces[iface]] = inhomogTerm;
          Real factor = Real(isign);
          fluxStencil *= factor;
          a_inhoTerm += inhomogTerm*factor;
          a_vofStencilz += fluxStencil;
        }
      }
    }
#endif

    {
      if(a_graph.isIrregular(a_vof.m_pt))
        {
          LocalStencil<CELL, Real> ebflux;
          Real inhomogTerm; //for inhomogeneous bcs
          // no stencil on EB face, yet
          std::array<LocalStencil<CELL, Real>, SpaceDim > fullFluxStencil;
          a_ebflux.getEBFluxStencilFull(fullFluxStencil,
                                        inhomogTerm,
                                        a_vof,
                                        a_vof.getCutFace(),
                                        a_graph,
                                        a_momentData,
                                        a_neighType,
                                        a_dx);
          a_inhoTerm += inhomogTerm;
        }
    }

    //for some reason, EBAMRPoissonOp calculates beta weight without the EB flux
      {
        Real diagweight;
        // a_vofStencilx.print();
        bool foundVoF = a_vofStencilx.getValueAtVoF(diagweight, a_vof);
        PR_assert(foundVoF); //if this fails, the stencil has no diagonal value
        a_vofStencilx.setBetaWeight(diagweight);
#if DIM>=2
        foundVoF = a_vofStencily.getValueAtVoF(diagweight, a_vof);
        PR_assert(foundVoF); //if this fails, the stencil has no diagonal value
        a_vofStencily.setBetaWeight(diagweight);
#endif
#if DIM>=3
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
  getWLSGradStencil(std::array<LocalStencil<CELL, Real>, SpaceDim >& a_vofStencil,
                    vector<Real>                                & a_inhoTerm,
                    const EBIndex<CELL>                 & a_vof,
                    const EBGraph                       & a_graph,
                    const MomentHolder<order>           & a_momentData,
                    const GradWLSCell<order>            & a_cellConstruct,
                    int                                 a_neighType,
                    Real                                a_dx)
  {
    PR_TIME("getWLSGradStencil");
    using Chombo4::pout;

    a_vofStencil[0].clear();
#if DIM>=2
    a_vofStencil[1].clear();
#endif
#if DIM>=3
    a_vofStencil[2].clear();
#endif
    CH_assert(a_inhoTerm.size() >= DIM);
    for (int dir=0; dir!=DIM; dir++)
      {
        a_inhoTerm[dir] = 0;
      }

    Real volScaleFactor  = getRegularVolume(a_dx);

    auto& a_voludata = a_momentData.volumeMoments();
    Real   voluFrac;
    RealVect centroid;
    Normalizor<order>::normalizedVoluMoments(voluFrac, centroid, a_voludata, a_vof, a_graph, a_dx);

    for (int dir=0; dir!=DIM; dir++)
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

        if (a_cellConstruct.m_verbose) {pout() << " stencil dir " << dir << " " << a_vof.m_pt << endl;}

        IndexTM<int, DIM> deriv = IndexTM<int, DIM>::Zero;
        deriv[dir] = 1; // dir-derivative stencil
        a_cellConstruct.getCellStencil(a_vofStencil[dir],
                                       a_inhoTerm[dir],
                                       a_vof,
                                       a_graph,
                                       a_momentData,
                                       voluMom,
                                       a_dx,
                                       deriv);
        if (a_cellConstruct.m_verbose)
          {
            pout() << " inhomogTerm dir " << dir << " " << a_inhoTerm[dir] << endl;
            pout() << " grad stencil " << dir << " " << endl;
            a_vofStencil[dir].print();
          }
      }

   //  // x direction
//     {
//       IndMomDIM voluMom;
//       if (voluFrac < 1.0)
//         {
//           voluMom = a_voludata(a_vof, 0);
//           for(MomentIterator<DIM, order> iter; iter.ok(); ++iter)
//             {
//               const auto& p = iter();
//               voluMom[iter()] = voluMom[iter()]/(volScaleFactor*POW(a_dx, p));
//             }
//         }
//       else
//         {
//           voluMom.setToRegular(a_dx);
//           for(MomentIterator<DIM, order> iter; iter.ok(); ++iter)
//             {
//               const auto& p = iter();
//               voluMom[iter()] = voluMom[iter()]/(volScaleFactor*POW(a_dx, p));
//             }
//         }

//       if (a_cellConstruct.m_verbose) {pout() << " stencil xdir "<< a_vof.m_pt << endl;}

//       IndexTM<int, DIM> deriv = IndexTM<int, DIM>::Zero;
//       deriv[0] = 1; // x-derivative stencil
//       a_cellConstruct.getCellStencil(a_vofStencil[0],
//                                      a_inhoTerm[0],
//                                      a_vof,
//                                      a_graph,
//                                      a_momentData,
//                                      voluMom,
//                                      a_dx,
//                                      deriv);
//       if (a_cellConstruct.m_verbose)
//         {
//           pout() << " inhomogTerm xdir " << a_inhoTerm[0] << endl;
//           pout() << " grad stencil " << endl;
//           a_vofStencil[0].print();
//         }
//     }
//     // y-direction
//     {
//       IndMomDIM voluMom;
//       if (voluFrac < 1.0)
//         {
//           voluMom = a_voludata(a_vof, 0);
//           for(MomentIterator<DIM, order> iter; iter.ok(); ++iter)
//             {
//               const auto& p = iter();
//               voluMom[iter()] = voluMom[iter()]/(volScaleFactor*POW(a_dx, p));
//             }
//         }
//       else
//         {
//           voluMom.setToRegular(a_dx);
//           for(MomentIterator<DIM, order> iter; iter.ok(); ++iter)
//             {
//               const auto& p = iter();
//               voluMom[iter()] = voluMom[iter()]/(volScaleFactor*POW(a_dx, p));
//             }
//         }

//       if (a_cellConstruct.m_verbose) {pout() << " stencil ydir "<< a_vof.m_pt << endl;}

//       IndexTM<int, DIM> deriv = IndexTM<int, DIM>::Zero;
//       deriv[1] = 1; // y-derivative stencil
//       a_cellConstruct.getCellStencil(a_vofStencil[1],
//                                      a_inhoTerm[1],
//                                      a_vof,
//                                      a_graph,
//                                      a_momentData,
//                                      voluMom,
//                                      a_dx,
//                                      deriv);
//       if (a_cellConstruct.m_verbose)
//         {
//           pout() << " inhomogTerm ydir " << a_inhoTerm[1] << endl;
//           pout() << " stencil " << endl;
//           a_vofStencil[1].print();
//         }
//     }
// #if DIM>=3
//     // z direction
//     {
//       IndMomDIM voluMom;
//       if (voluFrac < 1.0)
//         {
//           voluMom = a_voludata(a_vof, 0);
//           for(MomentIterator<DIM, order> iter; iter.ok(); ++iter)
//             {
//               const auto& p = iter();
//               voluMom[iter()] = voluMom[iter()]/(volScaleFactor*POW(a_dx, p));
//             }
//         }
//       else
//         {
//           voluMom.setToRegular(a_dx);
//           for(MomentIterator<DIM, order> iter; iter.ok(); ++iter)
//             {
//               const auto& p = iter();
//               voluMom[iter()] = voluMom[iter()]/(volScaleFactor*POW(a_dx, p));
//             }
//         }

//       if (a_cellConstruct.m_verbose) {pout() << " stencil zdir "<< a_vof.m_pt << endl;}

//       IndexTM<int, DIM> deriv = IndexTM<int, DIM>::Zero;
//       deriv[2] = 1; // z-derivative stencil
//       a_cellConstruct.getCellStencil(a_vofStencil[2],
//                                      a_inhoTerm[2],
//                                      a_vof,
//                                      a_graph,
//                                      a_momentData,
//                                      voluMom,
//                                      a_dx,
//                                      deriv);
//       if (a_cellConstruct.m_verbose)
//         {
//           pout() << " inhomogTerm zdir " << a_inhoTerm[2] << endl;
//           pout() << " grad stencil " << endl;
//           a_vofStencil[2].print();
//         }
//     }
// #endif

    // REVIEW does this matter?
    // //for some reason, EBAMRPoissonOp calculates beta weight without the EB flux
//       {
//         Real diagweight;
//         // a_vofStencilx.print();
//         bool foundVoF = a_vofStencil[0].getValueAtVoF(diagweight, a_vof);
//         PR_assert(foundVoF); //if this fails, the stencil has no diagonal value
//         a_vofStencil[0].setBetaWeight(diagweight);

//         foundVoF = a_vofStencil[1].getValueAtVoF(diagweight, a_vof);
//         PR_assert(foundVoF); //if this fails, the stencil has no diagonal value
//         a_vofStencil[1].setBetaWeight(diagweight);
// #if DIM>=3
//         foundVoF = a_vofStencil[2].getValueAtVoF(diagweight, a_vof);
//         PR_assert(foundVoF); //if this fails, the stencil has no diagonal value
//         a_vofStencil[2].setBetaWeight(diagweight);
// #endif
//       }
  }

}
#endif
