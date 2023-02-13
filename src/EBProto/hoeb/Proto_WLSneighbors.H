#ifndef _Proto_WLSneighbors_H_
#define _Proto_WLSneighbors_H_

#include "Chombo_ProtoInterface.H"
#include "EBProto.H"
#include "Proto_MomentHolder.H"

namespace Proto
{

  //Class that will contain all the information related to neighbors VoFs,
  //boundary faces and eb faces in a certain radius.
  //FIXME this class is too transparent about how data is stored,
  //hiding away vector storage would be helpful
  template<int order>
  class Neighbors
  {
  public:
    typedef IndexedMoments<DIM  , order> IndMomDIM;
    typedef IndexedMoments<DIM-1, order> IndMomSDMinOne;
    typedef MomentIterator<DIM  , order> MomItDIM;
    typedef MomentIterator<DIM-1, order> MomItSDMinOne;

    Neighbors()
      {
        clear();
      };

    inline void clear()
      {
        m_cells.resize(0);
        m_volMom.resize(0);
        m_facMom.resize(0);
        m_ebfMom.resize(0);
        m_volLoc.resize(0);
        m_facLoc.resize(0);
        m_ebfLoc.resize(0);
        m_weightDisVol.resize(0);
        m_weightDisFac.resize(0);
        m_weightDisEB.resize(0);
        m_volneighFrac.resize(0);
        m_ebfneighFrac.resize(0);
        m_facneighFrac.resize(0);
        m_ebfneighNorm.resize(0);
        m_facneighDir.resize(0);
        m_facneighSid.resize(0);
      };

    //Get neighbors VoFs through the faces of a specific cell.
    vector<EBIndex<CELL> >
    getFlapVoFs(const EBGraph       & a_graph,
                const EBIndex<CELL> & a_start) const;

    //Get neighbors VoFs of a given face within a specific radius.
    //Radius starts at one for the cells sharing the face, and
    //increases by one for the neighbor cells sharing a face with the original
    //cells of radius one.
    template<CENTERING cent>
    Vector<EBIndex<CELL> >
    getVoFsNeighbors(const EBIndex<cent> & a_face,
                     const EBIndex<CELL> & a_vof,
                     const int             a_radius,
                     const EBGraph       & a_graph,
                     int                   a_neighType);

    template<CENTERING cent>
    Vector<EBIndex<CELL> >
    getNeighborsFromFace(const EBIndex<cent> & a_face,
                         const int             a_baseRadius,
                         const EBGraph       & a_graph) const;

    Vector<EBIndex<CELL> >
    getNeighborsFromCell(const EBIndex<CELL> & a_vof,
                         const int             a_baseRadius,
                         const EBGraph       & a_graph) const;

    void
    growNeighbors(Vector<EBIndex<CELL> > &a_base,
                  const EBIndex<CELL> & a_vof,
                  const EBGraph       & a_graph);

    // New version, use this
    void
    getMomentsInfoFromNeighbors(Vector<EBIndex<CELL> > & a_vols,
                                const EBGraph          & a_graph,
                                const MomentHolder<order>  & a_momentData,
                                string                   a_ebbcName,
                                string                   a_dombcName[2*DIM], // order is is not obvious
                                Real                   a_dx,
                                bool a_verbose = false,
                                int a_fluxType = 0);

    // This is aiming to generate neighbors and moments, and have some way of checking the neighborhood is big enough
    // FIXME this is a work in progress, and uses some adhoc checks at the moment
    void
    getNeighborsAndMoments(const EBGraph              & a_graph,
                           const MomentHolder<order>  & a_momentData,
                           string                       a_ebbcName,
                           string                       a_dombcName[2*DIM], // order is is not obvious
                           Real                         a_dx,
                           bool a_verbose = false)
      {
        MomentIterator<DIM,order> momit;
        Vector<EBIndex<CELL> > vols;
        // REVIEW Make sure stencils are big enough. Size is not a sufficient measure, this should check the stencil is solvable
        // while (this->size() <= momit.size()+DIM)
        //   {
        //     vols = neighs.getVoFsNeighbors(a_face, a_vof, radius, a_graph, a_neighType);

        //     getMomentsInfoFromNeighbors(vols,
        //                                 a_graph,
        //                                 a_momentData,
        //                                 a_ebbcName, a_dombcName, a_dx, a_verbose);

        //     // if (grows>0) pout() << "neigh set tight for " << a_vof.m_pt << ", growing radius by 1" << endl;
        //     radius += 1;
        //     grows += 1;
        //   }
        // CH_assert(vols.size() > momit.size());
        // CH_assert(neighs.size() > momit.size());
      }

    inline void
    setWeightDisVol (Vector<int> & a_weightDisVol)
      {
        m_weightDisVol = a_weightDisVol;
      };

    inline void
    setMomentsInfo(Vector<IndMomDIM>      & a_volMom,
                   Vector<IndMomSDMinOne> & a_facMom,
                   Vector<IndMomDIM>      & a_ebfMom,
                   Vector<Vector<IndMomDIM> > & a_ebfneighNorm,
                   Real                     a_dx)
      {
        m_volMom.resize(a_volMom.size());
        Real volScaleFactor  = (D_TERM(a_dx, *a_dx, *a_dx));
        Real areaScaleFactor = (D_TERM(1.0,  *a_dx, *a_dx));
        for(int i = 0; i < a_volMom.size(); ++i)
          {
            IndMomDIM momspt = a_volMom[i];
            for(MomentIterator<DIM, order> iter; iter.ok(); ++iter)
              {
                const IndexTM<int, DIM>& p = iter();
                m_volMom[i][iter()] = momspt[iter()]/(volScaleFactor*POW(a_dx,p));
              }
          }
        m_facMom.resize(a_facMom.size());
        for(int i = 0; i < a_facMom.size(); ++i)
          {
            IndMomSDMinOne momspt = a_facMom[i];
            for(MomentIterator<DIM-1, order> iter; iter.ok(); ++iter)
              {
                const IndexTM<int, DIM-1>& p = iter();
                m_facMom[i][iter()] = momspt[iter()]/(areaScaleFactor*POW(a_dx,p));
              }
          }
        m_ebfMom.resize(a_ebfMom.size());
        for(int i = 0; i < a_ebfMom.size(); ++i)
          {
            IndMomDIM momspt = a_ebfMom[i];
            for(MomentIterator<DIM, order> iter; iter.ok(); ++iter)
              {
                const IndexTM<int, DIM>& p = iter();
                m_ebfMom[i][iter()] = momspt[iter()]/(areaScaleFactor*POW(a_dx,p));
              }
          }
        m_ebfneighNorm.resize(a_ebfneighNorm.size());
        for (int i = 0; i < a_ebfneighNorm.size(); ++i)
          {
            m_ebfneighNorm[i].resize(a_ebfneighNorm[i].size());
            for(int j = 0; j < a_ebfneighNorm[i].size(); ++j)
              {
                IndMomDIM momspt = a_ebfneighNorm[i][j];
                for(MomentIterator<DIM, order> iter; iter.ok(); ++iter)
                  {
                    const IndexTM<int, DIM>& p = iter();
                    m_ebfneighNorm[i][j][iter()] = momspt[iter()]/(areaScaleFactor*POW(a_dx,p));
                  }
              }
          }
      };

    inline Vector<Vector<IndMomDIM> > getEBNormMoments () const
      {
        Vector<Vector<IndMomDIM> > ebfneighNorm;
        ebfneighNorm.resize(m_ebfneighNorm.size());
        for (int i = 0; i < m_ebfneighNorm.size(); ++i)
          {
            ebfneighNorm[i] = m_ebfneighNorm[i];
          }
        return ebfneighNorm;
      };

    inline void setLocInfo (Vector<RealVect> & a_volLoc,
                            Vector<RealVect> & a_facLoc,
                            Vector<RealVect> & a_ebfLoc)
      {
        m_volLoc = a_volLoc;
        m_facLoc = a_facLoc;
        m_ebfLoc = a_ebfLoc;
      };

    inline void setWeightDisFacAndEB (Vector<int> & a_weightDisFac,
                                      Vector<int> & a_weightDisEB)
      {
        m_weightDisFac = a_weightDisFac;
        m_weightDisEB = a_weightDisEB;
      };

    inline void setFracInfo (Vector<Real> & a_volneighFrac,
                             Vector<Real> & a_ebfneighFrac,
                             Vector<Real> & a_facneighFrac)
      {
        m_volneighFrac = a_volneighFrac;
        m_ebfneighFrac = a_ebfneighFrac;
        m_facneighFrac = a_facneighFrac;
      };

    inline void setAdditionalFaceInfo (Vector<int>            & a_facneighDir,
                                       Vector<Side::LoHiSide> & a_facneighSid)
      {
        m_facneighDir = a_facneighDir;
        m_facneighSid = a_facneighSid;
      };

    // How many neighbors are there
    inline unsigned int
    sizeVols() const
      {
        return m_volLoc.size();
      }

    inline unsigned int
    sizeFaces() const
      {
        return m_facLoc.size();
      }

    inline unsigned int
    sizeEB() const
      {
        return m_ebfLoc.size();
      }

    inline unsigned int
    size() const
      {
        return sizeVols() + sizeFaces() + sizeEB();
        // m_volMom.size() + m_facMom.size() + m_ebfMom.size();
      }

    // get a global neighbor index from the neighbor collection
    inline unsigned int
    indexFromVol(unsigned int a_volIndex) const
      {
        CH_assert(a_volIndex < sizeVols());
        return a_volIndex;
      }

    inline unsigned int
    indexFromFace(unsigned int a_faceIndex) const
      {
        CH_assert(a_faceIndex < sizeFaces());
        return a_faceIndex + sizeVols();
      }

    inline unsigned int
    indexFromEB(unsigned int a_ebIndex) const
      {
        CH_assert(a_ebIndex < sizeEB());
        return a_ebIndex + sizeVols() + sizeFaces();
      }

    inline bool
    isRegular() const
      {
        return (sizeEB() == 0);
      }

    inline bool
    isIterior() const
      {
        return ((sizeEB() == 0) & (sizeFaces() == 0));
      }

    // Make sure neighbor storage matches
    inline void
    check() const
    {
      CH_assert(sizeVols() == m_volMom.size());
      CH_assert(sizeVols() == m_volLoc.size());
      CH_assert(sizeVols() == m_weightDisVol.size());
      CH_assert(sizeVols() == m_volneighFrac.size());

      CH_assert(sizeFaces() == m_facMom.size());
      CH_assert(sizeFaces() == m_facLoc.size());
      CH_assert(sizeFaces() == m_weightDisFac.size());
      CH_assert(sizeFaces() == m_facneighFrac.size());
      CH_assert(sizeFaces() == m_facneighDir.size());
      CH_assert(sizeFaces() == m_facneighSid.size());

      CH_assert(sizeEB() == m_ebfMom.size());
      CH_assert(sizeEB() == m_ebfLoc.size());
      CH_assert(sizeEB() == m_weightDisEB.size());
      CH_assert(sizeEB() == m_ebfneighFrac.size());
      CH_assert(sizeEB() == m_ebfneighNorm.size());
    }

    // the cells used in the neighborhood
    Vector<EBIndex<CELL> >  m_cells;
    //Volume moments
    Vector<IndMomDIM>       m_volMom;
    //Face moments
    Vector<IndMomSDMinOne>  m_facMom;
    //EB Faces moments
    Vector<IndMomDIM>       m_ebfMom;
    //EB Normal moments
    Vector<Vector<IndMomDIM> > m_ebfneighNorm;
    //Cell location, will be the cell center
    Vector<RealVect>        m_volLoc;
    //Boundary Face location, will be the face center
    Vector<RealVect>        m_facLoc;
    //EB Face location, will be the cell center containing the EB Face
    Vector<RealVect>        m_ebfLoc;
    //Contain the manhattan distance from the face
    //or EB face being evaluated to the cell
    Vector<int>             m_weightDisVol;
    //Contain the manhattan distance from the face
    //or EB face being evaluated to the boundary face
    Vector<int>             m_weightDisFac;
    //Contain the manhattan distance from the face
    //or EB face being evaluated to the EB Face
    Vector<int>             m_weightDisEB;
    //Volume fraction
    Vector<Real>          m_volneighFrac;
    //EB Face area fraction
    Vector<Real>          m_ebfneighFrac;
    //Face area fraction
    Vector<Real>          m_facneighFrac;
    //Face direction
    Vector<int>             m_facneighDir;
    //Face Side -- lo or hi
    Vector<Side::LoHiSide>  m_facneighSid;
  };


  //Get neighbors VoFs through the faces of a specific cell.
  template<int order>
  vector<EBIndex<CELL> >
  Neighbors<order>::
  getFlapVoFs(const EBGraph& a_graph,
              const EBIndex<CELL>& a_start) const
  {
    PR_TIME("getFlapVoFs");
    vector<EBIndex<CELL> > retval;

    for(int idir  = 0; idir <  DIM; idir++)
      {
        for(SideIterator sit; sit.ok(); ++sit)
          {
            vector<EBIndex<CELL> > neighbors = a_graph.getNeighborVoFs(a_start, idir, sit());
            // append neighbors to retval
            retval.insert(retval.end(), neighbors.begin(), neighbors.end());
          }
      }
    return retval;
  }

  //Get neighbors VoFs of a given face within a specific radius.
  //Radius starts at one for the cells sharing the face, and
  //increases by one for the neighbor cells sharing a face with the original
  //cells of radius one.
  template<int order>
  template<CENTERING cent>
  Vector<EBIndex<CELL> >
  Neighbors<order>::
  getVoFsNeighbors(const EBIndex<cent> & a_face,
                   const EBIndex<CELL> & a_vof,
                   const int             a_radius,
                   const EBGraph       & a_graph,
                   int                   a_neighType)
  {
    PR_TIME("getVoFsNeighbors");
    m_weightDisVol.resize(0);
    m_cells.resize(0);
    MomentIterator<DIM,order> momit;
    if (a_neighType == 1)
      {
        m_cells.push_back(a_vof);
        m_weightDisVol.push_back(1);
      }
    else
      {
        if (cent == BOUNDARY)
          {
            EBIndex<CELL> voflo = a_face.getVoF();
            m_cells.push_back(voflo);
            m_weightDisVol.push_back(1);
          }
        else
          {
            EBIndex<CELL> voflo = a_face.getVoF(Side::Lo);
            EBIndex<CELL> vofhi = a_face.getVoF(Side::Hi);
            if (!a_face.m_isBoundary)
              {
                m_cells.push_back(voflo);
                m_weightDisVol.push_back(1);
                m_cells.push_back(vofhi);
                m_weightDisVol.push_back(1);
              }
            else
              {
                if (a_face.m_vofIDMe>=0)
                  {
                    m_cells.push_back(vofhi);
                    m_weightDisVol.push_back(1);
                  }
                if (a_face.m_vofIDLo>=0)
                  {
                    m_cells.push_back(voflo);
                    m_weightDisVol.push_back(1);
                  }
              }
          }
      }
    int startV = 0;
    int new_rad = a_radius;
    for (int irad = 1; irad<=new_rad; irad++)
      {
        int neighSize = m_cells.size();
        for (int ineigh = startV; ineigh<neighSize; ineigh++){
          vector< EBIndex<CELL> > volsinclusive = getFlapVoFs(a_graph, m_cells[ineigh]);
          for(int ivof = 0; ivof < volsinclusive.size(); ivof++)
            {
              bool incl = true;
              for(int ivof2 = 0; ivof2 < m_cells.size(); ivof2++)
                {
                  if (volsinclusive[ivof]==m_cells[ivof2])
                    {
                      incl=false;
                      break;
                    }
                }
              if (incl)
                {
                  m_cells.push_back(volsinclusive[ivof]);
                  m_weightDisVol.push_back(irad+1);
                }
            }
        }
        startV = neighSize;
        if (irad==new_rad && m_cells.size()<=1.0*momit.size()) new_rad++;
      }
    CH_assert(m_cells.size() == m_weightDisVol.size());
    return m_cells;
  }

  //Gets all the moments and geometric info for the neighbors cells, boundary faces
  //and EB faces. Stores all this info in the Neighbors object.
  template<int order>
  void
  Neighbors<order>::
  getMomentsInfoFromNeighbors(Vector<EBIndex<CELL> >  & a_vols,
                              const EBGraph           & a_graph,
                              const MomentHolder<order>  & a_momentData,
                              string                    a_ebbcName,
                              string                    a_dombcName[2*DIM],
                              Real                    a_dx,
                              bool a_verbose,
                              int a_fluxType)
  {
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

    PR_TIME("getMomentsInfoFromNeighbors");
    using Chombo4::pout;
    Vector<IndMomDIM> volmom, ebfmom;
    Vector<IndMomSDMinOne> facmom;
    Vector<RealVect> volneighLoc;
    Vector<RealVect> ebfneighLoc;
    Vector<RealVect> facneighLoc;
    Vector<Real> volneighFrac;
    Vector<Real> ebfneighFrac;
    Vector<Real> facneighFrac;
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
    Vector<int> weightDisVol = m_weightDisVol;

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
            if (a_fluxType==0){
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
        if (a_fluxType==0){
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
                         && (a_dombcName[index] != string("no_bcs")))
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
                      int index = ebp_index(1, sit());
                      if(!found && y_face.m_isBoundary
                         && (a_dombcName[index] != string("no_bcs")))
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
                      int index = ebp_index(2, sit());
                      if(!found && z_face.m_isBoundary
                         && (a_dombcName[index] != string("no_bcs")))
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
    setMomentsInfo(volmom, facmom,
                           ebfmom, ebfneighNorm,
                           a_dx);

    setLocInfo(volneighLoc, facneighLoc,
                       ebfneighLoc);

    setWeightDisFacAndEB(weightDisFac,
                                 weightDisEB);

    setFracInfo(volneighFrac,
                        ebfneighFrac,
                        facneighFrac);

    setAdditionalFaceInfo(facneighDir,
                                  facneighSid);

    if (a_verbose)
      {
        pout()<<" Index order of moments for the matrices " << endl;
        for (int i=0; i!=sizeVols(); i++)
          {
            pout() << " Vol " << i << ", neighbor " << indexFromVol(i) << endl;
          }
        for (int i=0; i!=sizeFaces(); i++)
          {
            pout() << " Face " << i << ", neighbor " << indexFromFace(i) << endl;
          }
        for (int i=0; i!=sizeEB(); i++)
          {
            pout() << " EB " << i << ", neighbor " << indexFromEB(i) << endl;
          }
      }
    check();
  }

} // end proto namespace

#endif //_Proto_WLSneighbors_H_
