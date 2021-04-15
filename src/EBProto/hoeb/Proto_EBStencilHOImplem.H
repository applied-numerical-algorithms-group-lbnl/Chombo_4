
#ifndef _Proto_EBStencilHOImplem_H_
#define _Proto_EBStencilHOImplem_H_
#include "../implem/Proto_Normalizor.H"
#include "Proto_WLSCell.H"

namespace Proto
{
  
  ///this deep down in the bowels, everything is a template specialization
  ///returns whether there are domain boundary conditions
  template <int order>
  bool 
  EBStencilHO<CELL, CELL,  order, double>::
  getStencil(vector<EBIndex<CELL> >                    & a_dstVoFs,
             vector<LocalStencil<CELL, double> >       & a_stencil,
             Stencil<double>                           & a_regStencil,
             vector<Stencil<double> >                    a_regStencilBCS[2*DIM],
             vector<Box>                                 a_BCApplyBoxes[2*DIM],       
             Box                                       & a_regApplyBox,
             bool                                      & a_bcOnly,
             string                                      a_stencilName,
             string                                      a_domainBCName[2*DIM],
             string                                      a_ebbcName,
             const Box                                 & a_srcValid,
             const Box                                 & a_dstValid,
             const Box                                 & a_srcDomain,
             const Box                                 & a_dstDomain,
             const Point                               & a_srcGhost,
             const Point                               & a_dstGhost,
             const EBGraph                             & a_srcGraph,
             const EBGraph                             & a_dstGraph,
             const VoluData                            & a_voludata,
             const EBFaData                            & a_ebfadata,
             const XFacData                            & a_xfacdata,
             const YFacData                            & a_yfacdata,
             const ZFacData                            & a_zfacdata,
	     const EBNormalData                        & a_ebnormxdata,
             const EBNormalData                        & a_ebnormydata,
             const EBNormalData                        & a_ebnormzdata,
             double                                      a_dx,
             bool                                        a_periodiconly, 
             Point                                       a_dstGrown)               
  {
    a_bcOnly = false;
    if(a_stencilName == string("Weighted_Least_Squares_Poisson"))
    {
      bool forceNeumann = false;
      return
        getWLSStencil(a_dstVoFs,
                      a_stencil,
                      a_regStencil,
                      a_regStencilBCS,
                      a_BCApplyBoxes,
                      a_regApplyBox,
                      a_stencilName,
                      a_domainBCName,
                      a_ebbcName,
                      a_srcValid,
                      a_srcDomain,
                      a_srcGhost,
                      a_dstGhost,
                      a_srcGraph,
                      a_voludata,
                      a_ebfadata,
                      a_xfacdata,
                      a_yfacdata,
                      a_zfacdata,
		      a_ebnormxdata,
                      a_ebnormydata,
                      a_ebnormzdata,
                      a_dx,
                      a_periodiconly,
                      forceNeumann);
    }
    else if (a_stencilName == string("Weighted_Least_Squares_Poisson_All_Neumann"))
    {
      bool forceNeumann = true;
      return
        getWLSStencil(a_dstVoFs,
                      a_stencil,
                      a_regStencil,
                      a_regStencilBCS,
                      a_BCApplyBoxes,
                      a_regApplyBox,
                      a_stencilName,
                      a_domainBCName,
                      a_ebbcName,
                      a_srcValid,
                      a_srcDomain,
                      a_srcGhost,
                      a_dstGhost,
                      a_srcGraph,
                      a_voludata,
                      a_ebfadata,
                      a_xfacdata,
                      a_yfacdata,
                      a_zfacdata,
                      a_ebnormxdata,
                      a_ebnormydata,
                      a_ebnormzdata,
                      a_dx,
                      a_periodiconly,
                      forceNeumann);
    }
    else
    {
      PR_error("stencil_archive_cell_cell::stencil string not found");
    }
    return false;//never gets here but this shuts up some compilers;
  }
 
  template <int order>
  bool
  EBStencilHO<CELL, CELL,  order, double>::
  getWLSStencil(vector<EBIndex<CELL> >                    & a_dstVoFs,                    
                vector<LocalStencil<CELL, double> >       & a_stencil,                    
                Stencil<double>                           & a_regStencil,         
                vector<Stencil<double> >                    a_regStencilBC[2*DIM],       
                vector<Box>                                 a_BCApplyBoxes[2*DIM],       
                Box                                       & a_regApplyBox,                
                string                                      a_stencilName,                
                string                                      a_domainBCName[2*DIM],               
                string                                      a_ebbcName,                   
                const Box                                 & a_validBox,                   
                const Box                                 & a_domain,                     
                const Point                               & a_srcGhost,                   
                const Point                               & a_dstGhost,                   
                const EBGraph                             & a_graph,                      
                const VoluData                            & a_voludata,                   
                const EBFaData                            & a_ebfadata,                   
                const XFacData                            & a_xfacdata,                   
                const YFacData                            & a_yfacdata,                   
                const ZFacData                            & a_zfacdata,
                const EBNormalData                        & a_ebnormxdata,
                const EBNormalData                        & a_ebnormydata,
                const EBNormalData                        & a_ebnormzdata,
                double                                      a_dx,                         
                bool                                        a_periodiconly,               
                bool  a_forceNeumann)               
  {
//    using std::cout;
//    using std::endl;
//    a_regStencil = Stencil<double>::Laplacian();
//    cout << "done with regular sten init" << endl;
//    a_regStencil *= 1.0/(a_dx*a_dx);
//    cout << "done with regular sten scaling" << endl;
    a_regApplyBox = a_validBox;
    //this stencil has a span so we have to evaluate the stencil on a grown set of the irregular cells so the regular and irregular data actually talk
    //getIrregLocations(a_dstVoFs, a_regStencil, a_validBox, a_graph);
    //cout << "done with getIrregLocations " <<a_dstVoFs.size()<< endl;

    a_dstVoFs = a_graph.getAllVoFs(a_validBox);

    BaseWLSFlux<XFACE, order> xflux;
    BaseWLSFlux<YFACE, order> yflux;
    BaseWLSFlux<ZFACE, order> zflux;
    BaseWLSFlux<BOUNDARY, order> ebflux;
    a_stencil.resize(a_dstVoFs.size());
    string ebbc = a_ebbcName;
    string dombc[2*DIM];
    if(a_forceNeumann)
    {
      ebbc = string("Neumann");
    }
    if(!a_periodiconly)
    {
      for(int ivec = 0; ivec < 2*DIM; ivec++)
      {
        if(a_forceNeumann)
        {
          dombc[ivec] = string("Neumann");
        }
        else
        {
          dombc[ivec] = a_domainBCName[ivec];
        }
      }
    }
    for(int ivof = 0; ivof < a_dstVoFs.size(); ivof++)
    {
      //cout << "should not get here with no irregular cells" << endl;
      Point pt = a_dstVoFs[ivof].m_pt;
      cout << "point1 "<<pt<<endl;
      WLSStencil<order>::
        getWLSDivFStencil(a_stencil[ivof],
                               a_dstVoFs[ivof],
                               a_graph,
                               a_voludata,
                               a_ebfadata,
                               a_xfacdata, a_yfacdata, a_zfacdata,
			       a_ebnormxdata, a_ebnormydata, 
			       a_ebnormzdata,
                               xflux, yflux, zflux, ebflux, ebbc, dombc,
                               a_dx);
      RealVect centroid;
      double kappa;
      Normalizor<order>::normalizedVoluMoments(kappa, centroid, a_voludata, a_dstVoFs[ivof], a_graph, a_dx);
      a_stencil[ivof]*= 1.0/(a_dx*a_dx);
      cout << "divF, ";
      cout << "stencil for " << a_dstVoFs[ivof].m_pt << ":" << endl;
      a_stencil[ivof].print();
    }
    return true;
  }
  /***/
/////////////end cell to cell

}
#endif
