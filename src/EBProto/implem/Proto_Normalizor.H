#ifndef _Proto_Normalizor__H_
#define _Proto_Normalizor__H_


#include "EBProto.H"

namespace Proto
{

  // Utilities related to moments
  // get the volume of a regular cell, given dx. Dimensionality is accounted for
  inline Real
  getRegularVolume(const Real a_dx)
  {
    return D_TERM(a_dx, *a_dx, *a_dx);
  }

  // get the volume of a regular face, given dx. Dimensionality is accounted for
  inline Real
  getRegularArea(const Real a_dx)
  {
    return D_TERM(1, *a_dx, *a_dx);
  }


  /// This class is used to provide nomalized moments (removes dependence on grid spacing)
  template <int order>
  class Normalizor
  {
  public:
    typedef IndexTM<int  , DIM>                  IvDIM;
    typedef IndexTM<int  , DIM-1>                IvSDMinOne;
    typedef IndexedMoments<DIM  , order> IndMomDIM;
    typedef IndexedMoments<DIM-1, order> IndMomSDMinOne;
    typedef MomentIterator<DIM  , order> MomItDIM;
    typedef MomentIterator<DIM-1, order> MomItSDMinOne;
    typedef HostIrregData<CELL    ,  IndMomDIM , 1>  VoluData;
    typedef HostIrregData<BOUNDARY,  IndMomDIM , 1>  EBFaData;
    typedef HostIrregData<XFACE, IndMomSDMinOne, 1>  XFacData;
    typedef HostIrregData<YFACE, IndMomSDMinOne, 1>  YFacData;
    typedef HostIrregData<ZFACE, IndMomSDMinOne, 1>  ZFacData;

    ///
    // the centroid is relative to the cartesian face center with unit length
    template <CENTERING cent>
    inline static void
    normalizedFaceMoments(double                                                      & a_areaFrac,
                          RealVect                                                    & a_centroid,
                          const HostIrregData<cent, IndexedMoments<DIM-1, order>, 1 > & a_facedat,
                          const EBIndex<cent>                                         & a_face,
                          const EBGraph                                               & a_graph,
                          double                                                        a_dx   )
    {


      if(a_graph.isFaceRegular(a_face))
      {
        a_areaFrac = 1;
        a_centroid  = RealVect::Zero();
      }
      else if(a_graph.isFaceCovered(a_face))
      {
        a_areaFrac = 0;
        a_centroid  = RealVect::Zero();
      }
      else
      {
        double fullfacearea = 1.;
        for(int idir = 0; idir < DIM-1; idir++) fullfacearea *= a_dx;
        const IndexedMoments<DIM-1, order> & facemom = a_facedat(a_face, 0); 
        double area = facemom[IvSDMinOne::Zero];
        a_areaFrac = area/fullfacearea;
        a_centroid = RealVect::Zero();
        if(area > 0)
        {
          int iindex = 0;
          for (int idir = 0; idir < DIM; ++idir)
          {
            if(idir != a_face.direction())
            {
              IvSDMinOne mono = BASISV_TM<int,DIM-1>(iindex);
              a_centroid[idir] = facemom[mono];
              //normalize  by area
              a_centroid[idir] /= area;
              //normalize by dx
              a_centroid[idir] /= a_dx;
              iindex++;
            }
          }
        }
      }
    }


    ///
    template <CENTERING cent>
    inline static void
    normalizedFaceMoments(double                                                       & a_areaFrac,
                          RealVect                                                     & a_centroid,
                          const HostIrregData<XFACE, IndexedMoments<DIM-1, order>, 1 > & a_xfacedat,
                          const HostIrregData<YFACE, IndexedMoments<DIM-1, order>, 1 > & a_yfacedat,
                          const HostIrregData<ZFACE, IndexedMoments<DIM-1, order>, 1 > & a_zfacedat,
                          const EBIndex<cent>                                          & a_face,
                          const EBGraph                                                & a_graph,
                          double                                                         a_dx,  
                          int a_idir)
    {
      if(a_idir == 0)
      {
        EBIndex<XFACE> face;
        face.m_pt         =  a_face.m_pt;
        face.m_vofIDMe    =  a_face.m_vofIDMe;
        face.m_vofIDLo    =  a_face.m_vofIDLo;
        face.m_isBoundary =  a_face.m_isBoundary;
        normalizedFaceMoments(a_areaFrac, a_centroid, a_xfacedat, face, a_graph, a_dx);
      }
      else if(a_idir == 1)
      {
        EBIndex<YFACE> face;
        face.m_pt         =  a_face.m_pt;
        face.m_vofIDMe    =  a_face.m_vofIDMe;
        face.m_vofIDLo    =  a_face.m_vofIDLo;
        face.m_isBoundary =  a_face.m_isBoundary;
        normalizedFaceMoments(a_areaFrac, a_centroid, a_yfacedat, face, a_graph, a_dx);
      }
#if DIM==3
      else if(a_idir == 2)
      {
        EBIndex<ZFACE> face;
        face.m_pt         =  a_face.m_pt;
        face.m_vofIDMe    =  a_face.m_vofIDMe;
        face.m_vofIDLo    =  a_face.m_vofIDLo;
        face.m_isBoundary =  a_face.m_isBoundary;
        normalizedFaceMoments(a_areaFrac, a_centroid, a_zfacedat, face, a_graph, a_dx);
      }
#endif
      else
      {
        PROTO_ASSERT(false, "bogus direction");
      }

    }
    
    // the centroid is relative to the cartesian cell center with unit length
    inline static void
    normalizedVoluMoments(double                                                      & a_volFrac,
                          RealVect                                                    & a_centroid,
                          const HostIrregData<CELL, IndexedMoments<DIM  , order> , 1> & a_celldat,
                          const EBIndex<CELL>                                         & a_vof,
                          const EBGraph                                               & a_graph,
                          double                                                        a_dx   )
    {

      if(a_graph.isRegular(a_vof.m_pt))
      {
        a_volFrac = 1;
        a_centroid  = RealVect::Zero();
      }
      else if(a_graph.isCovered(a_vof.m_pt))
      {
        a_volFrac = 0;
        a_centroid  = RealVect::Zero();
      }
      else
      {
        double fullcellvol = 1.;
        for(int idir = 0; idir < DIM; idir++) fullcellvol *= a_dx;
        const IndexedMoments<DIM  , order> & moments = a_celldat(a_vof, 0); 
        double volume = moments[IvDIM::Zero];
        a_volFrac = volume/fullcellvol;
        a_centroid = RealVect::Zero();
        if(volume > 0)
        {
          for (int idir = 0; idir < DIM; ++idir)
          {
            IvDIM mono = BASISV_TM<int,DIM>(idir);
            a_centroid[idir] = moments[mono];
            //normalize  by volume
            a_centroid[idir] /= volume;
            //normalize by dx
            a_centroid[idir] /= a_dx;
          }
        }
      }
    }


    ///
    inline static void
    normalizedBndryMoments(double                                                         & a_areaFrac,
                           RealVect                                                       & a_centroid,
                           const EBFaData                                                 & a_ebfadat,
                           const EBIndex<BOUNDARY>                                        & a_vof,
                           const EBGraph                                                  & a_graph,
                           double                                                           a_dx   )
    {

      if(a_graph.isRegular(a_vof.m_pt))
      {
        a_centroid  = RealVect::Zero();
      }
      else if(a_graph.isCovered(a_vof.m_pt))
      {
        a_centroid  = RealVect::Zero();
      }
      else
      {
//        int ideb = 0;
#if DIM==2
//        Point ivdeb(1,19);
//        if(ivdeb == a_vof.m_pt)
//        {
//          ideb = 1;
//        }
#endif
        double fullcellvol = 1.;
        for(int idir = 0; idir < DIM; idir++) fullcellvol *= a_dx;
        double fullcellarea = fullcellvol/a_dx;
        //double fullarea = fullcellvol/a_dx;
        const IndexedMoments<DIM,order> & moments = a_ebfadat(a_vof, 0); 
        double area = moments[IvDIM::Zero];
        a_areaFrac = area/fullcellarea;
        a_centroid = RealVect::Zero();
        if(area > 0)
        {
          for (int idir = 0; idir < DIM; ++idir)
          {
            IvDIM mono = BASISV_TM<int,DIM>(idir);
            a_centroid[idir] = moments[mono];
            //normalize  by area
            a_centroid[idir] /= area;
            //normalize by dx
            a_centroid[idir] /= a_dx;
          }
        }
      }
    }


    //sum of area fractions on a particular side of a vof
    inline static double
    sumArea(const VoluData                & a_voludata,
            const EBFaData                & a_ebfadata,
            const XFacData                & a_xfacdata,
            const YFacData                & a_yfacdata,
            const ZFacData                & a_zfacdata,
            const EBGraph                 & a_graph,
            const EBIndex<CELL>           & a_vof,
            const int                     & a_idir,
            const Side::LoHiSide          & a_sd,
            const double                  & a_dx)
    {
      
      double retval = 0;
      if(a_graph.isRegular(a_vof.m_pt))
      {
        retval = 1;
      }
      else if(a_graph.isCovered(a_vof.m_pt))
      {
        retval = 0;
      }
      else
      {
        if(a_idir == 0)
        {
          vector< EBIndex<XFACE> >faces = a_graph.getXFaces(a_vof, a_sd);
          for(int iface = 0; iface < faces.size(); iface++)
          {
            double   areaFrac;
            RealVect centroid;
            normalizedFaceMoments<XFACE>(areaFrac,
                                         centroid,
                                         a_xfacdata,
                                         faces[iface],
                                         a_graph,
                                         a_dx   );
            retval += areaFrac;
          }
        }
        else if(a_idir == 1)
        {
          vector< EBIndex<YFACE> >faces = a_graph.getYFaces(a_vof, a_sd);

          for(int iface = 0; iface < faces.size(); iface++)
          {
            double   areaFrac;
            RealVect centroid;
            normalizedFaceMoments<YFACE>(areaFrac,
                                         centroid,
                                         a_yfacdata,
                                         faces[iface],
                                         a_graph,
                                         a_dx   );
            retval += areaFrac;
          }
        }
#if DIM==3
        else if(a_idir == 2)
        {
          vector< EBIndex<ZFACE> >faces = a_graph.getZFaces(a_vof, a_sd);
          for(int iface = 0; iface < faces.size(); iface++)
          {
            double areaFrac;
            RealVect centroid;
            normalizedFaceMoments<ZFACE>(areaFrac,
                                         centroid,
                                         a_zfacdata,
                                         faces[iface],
                                         a_graph,
                                         a_dx   );
            retval += areaFrac;
          }
        }
#endif
        else
        {
          PROTO_ASSERT(false, "bogus idir");
        }
      } //else (is irregular)

      return retval;
    }

    template <typename data_t>
    inline static void 
    getNormalAndBoundaryArea(double                        & a_bndryArea,
                             RealVect                      & a_normal,
                             const VoluData                & a_voludata,
                             const EBFaData                & a_ebfadata,
                             const XFacData                & a_xfacdata,
                             const YFacData                & a_yfacdata,
                             const ZFacData                & a_zfacdata,
                             const EBIndex<CELL>           & a_vof,
                             const EBGraph                 & a_graph,
                             data_t                          a_dx   )
    {
      a_normal = BASISREALV(0); //give it something valid in case there is no area
      a_bndryArea = 0;
      if(a_graph.isIrregular(a_vof.m_pt))
      {

        data_t irregArea = 0.0;
        data_t hiArea[DIM];
        data_t loArea[DIM];
        for(int idir = 0; idir < DIM; idir++)
        {
          hiArea[idir] =  sumArea(a_voludata,
                                  a_ebfadata,
                                  a_xfacdata,
                                  a_yfacdata,
                                  a_zfacdata,
                                  a_graph,
                                  a_vof, idir,
                                  Side::Hi,
                                  a_dx);
          loArea[idir] =  sumArea(a_voludata,
                                  a_ebfadata,
                                  a_xfacdata,
                                  a_yfacdata,
                                  a_zfacdata,
                                  a_graph,
                                  a_vof, idir,
                                  Side::Lo,
                                  a_dx);
          irregArea += (hiArea[idir]-loArea[idir])*(hiArea[idir]-loArea[idir]);
        }
        a_bndryArea = sqrt(irregArea);
        if(a_bndryArea > 0)
        {
          for(int idir = 0; idir < DIM; idir++)
          {
            a_normal[idir] = (hiArea[idir]-loArea[idir])/a_bndryArea;
          }
        }
      }
    }

    // normalize a volume moment by removing dependence on dx
    inline static void
    normalizeFaceMoment(IndexedMoments<DIM-1, order> & a_faceMoment,
                        const Real a_dx)
      {
        Real areaScaleFactor = getRegularArea(a_dx);
        for(MomentIterator<DIM-1, order> iter; iter.ok(); ++iter)
          {
            const IndexTM<int, DIM-1>& p = iter();
            a_faceMoment[iter()] /= (areaScaleFactor*POW(a_dx, p));
          }
      }

    inline static void
    normalizeFaceMoment(IndexedMoments<DIM, order> & a_faceMoment,
                        const Real a_dx)
      {
        Real areaScaleFactor = getRegularArea(a_dx);
        for(MomentIterator<DIM, order> iter; iter.ok(); ++iter)
          {
            const IndexTM<int, DIM>& p = iter();
            a_faceMoment[iter()] /= (areaScaleFactor*POW(a_dx, p));
          }
      }

    // normalize a face moment by removing dependence on dx
    inline static void
    normalizeVolumeMoment(IndexedMoments<DIM, order> & a_volMoment,
                          const Real a_dx)
      {
        Real volumeScaleFactor = getRegularVolume(a_dx);
        for(MomentIterator<DIM, order> iter; iter.ok(); ++iter)
          {
            const IndexTM<int, DIM>& p = iter();
            a_volMoment[iter()] /= (volumeScaleFactor*POW(a_dx, p));
          }
      }

    template <CENTERING cent>
    inline static RealVect
    faceCenter(EBIndex<cent> a_face,
               const Real a_dx)
      {
        RealVect faceLoc;
        for (int idir = 0; idir != DIM; idir++)
          {
            faceLoc[idir] = a_dx * (0.5 + a_face.m_pt[idir]);
          }
        return faceLoc;
      }
  };

}
#include "Proto_EBStencilArchiveImplem.H"

#endif
  
