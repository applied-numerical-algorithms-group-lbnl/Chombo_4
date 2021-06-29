
#ifndef _Proto_EBStencilHOImplem_H_
#define _Proto_EBStencilHOImplem_H_
#include "Chombo_ParmParse.H"
#include "../implem/Proto_Normalizor.H"
#include "Proto_WLSCell.H"

namespace Proto
{
  //Gets high-order stencil. First implementation gives you the stencil for Poisson operator
  //All Neumann stencil is needed for EBMultigrid, it still needs to be checked.
  //Returns whether there are domain boundary conditions
  template <int order>
  bool 
  EBStencilHO<CELL, CELL,  order, double>::
  getStencil(vector<EBIndex<CELL> >                    & a_dstVoFs,
             vector<LocalStencil<CELL, double> >       & a_stencil,
             vector<double>                            & a_inhoContr,
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
        getWLSPoissonStencil(a_dstVoFs,
                             a_stencil,
                             a_inhoContr,
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
        getWLSPoissonStencil(a_dstVoFs,
                             a_stencil,
                             a_inhoContr,
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
    else if (a_stencilName == string("Weighted_Least_Squares_Advection"))
    {
      bool forceNeumann = false;
      return
        getWLSAdvectionStencil(a_dstVoFs,
                               a_stencil,
                               a_inhoContr,
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
  getWLSPoissonStencil(vector<EBIndex<CELL> >                    & a_dstVoFs,                    
                       vector<LocalStencil<CELL, double> >       & a_stencil,                    
                       vector<double>                            & a_inhoContr,
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
//    using std::pout();
//    using std::endl;
//    a_regStencil = Stencil<double>::Laplacian();
//    pout() << "done with regular sten init" << endl;
//    a_regStencil *= 1.0/(a_dx*a_dx);
//    pout() << "done with regular sten scaling" << endl;
    //this stencil has a span so we have to evaluate the stencil on a grown set of the irregular cells so the regular and irregular data actually talk
    //getIrregLocations(a_dstVoFs, a_regStencil, a_validBox, a_graph);
    //pout() << "done with getIrregLocations " <<a_dstVoFs.size()<< endl;

    //a_stencil will compute its corresponding stencil for all the cells in the domain
    a_regApplyBox = a_validBox;
    a_dstVoFs = a_graph.getAllVoFs(a_validBox);
    a_stencil.resize(a_dstVoFs.size());
    a_inhoContr.resize(a_dstVoFs.size());

    int radius = 2;
    int powerE = -3;
    bool verbose = 0;
    int exactSolution = 1;
    Vector<int> ptVec;

    //Reading parameters from the inputs file for WLS stencils
    ParmParse pp;
    pp.get("radius"   , radius);
    pp.get("powerE"   , powerE);
    pp.query("verbose"   , verbose);

    //Exact solution
    pp.get("exact_solution",exactSolution);

    if (verbose)
    {
      pp.getarr("pt_verbose",ptVec,0,DIM);
    } 

    string ebbc = a_ebbcName;
    string dombc[2*DIM];

    //Assigning domain and EB boundary conditions
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

    //Looping over all the valid VoFs
    for(int ivof = 0; ivof < a_dstVoFs.size(); ivof++)
    {
      Point pt = a_dstVoFs[ivof].m_pt;
      bool pt_verbose = false;
      if (verbose)
      {
#if DIM==3
	if (ptVec[0]==pt[0] && ptVec[1]==pt[1] && ptVec[2]==pt[2])
#else
	if (ptVec[0]==pt[0] && ptVec[1]==pt[1])
#endif
        {
	  pt_verbose = true;
          pout() << "Debugging Point "<<pt<<endl;
        }
	else
        {
          pt_verbose = false;
        }
      }

      //Creating the WLS fluxes needed for flux diveregence
      BaseWLSFlux<XFACE, order> xflux(radius, powerE, exactSolution, pt_verbose);
      BaseWLSFlux<YFACE, order> yflux(radius, powerE, exactSolution, pt_verbose);
      BaseWLSFlux<ZFACE, order> zflux(radius, powerE, exactSolution, pt_verbose);
      BaseWLSFlux<BOUNDARY, order> ebflux(radius, powerE, exactSolution, pt_verbose);
      Real inhoTerm = 0.0;
      //Neighbor list type. 0-FaceBased or 1-CellBased
      int neighType = 0;

      //Computing the flux divergence stencil
      WLSStencil<order>::
        getLaplaWLSDivFStencil(a_stencil[ivof],
                          a_dstVoFs[ivof],
                          a_graph,
                          a_voludata,
                          a_ebfadata,
                          a_xfacdata, a_yfacdata, a_zfacdata,
                          a_ebnormxdata, a_ebnormydata, 
                          a_ebnormzdata,
                          xflux, yflux, zflux, ebflux, ebbc, dombc,
			  inhoTerm,
                          neighType,
                          a_dx);

      RealVect centroid;
      double kappa;
      Normalizor<order>::normalizedVoluMoments(kappa, centroid, a_voludata, a_dstVoFs[ivof], a_graph, a_dx);
      a_stencil[ivof]*= 1.0/(a_dx*a_dx);
      inhoTerm*= 1.0/(a_dx*a_dx);
      a_inhoContr[ivof] = inhoTerm;
      if (pt_verbose) 
      {
        pout() << "divF stencil for " << a_dstVoFs[ivof].m_pt << ":" << endl;
        a_stencil[ivof].print();
	Real accDivF = 0.0;
	for(int ivec =  0; ivec < a_stencil[ivof].size(); ivec++)
        {
          accDivF += a_stencil[ivof].m_entries[ivec].m_weight;
        }

        pout() << "For phi=1 acc " << accDivF <<" inhoTerm "<<inhoTerm<< endl;
      }
    }
    return true;
  }
  /***/
/////////////end cell to cell

  template <int order>
  bool
  EBStencilHO<CELL, CELL,  order, double>::
  getWLSAdvectionStencil(vector<EBIndex<CELL> >                    & a_dstVoFs,                    
                         vector<LocalStencil<CELL, double> >       & a_stencil,                    
                         vector<double>                            & a_inhoContr,
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
    vector<LocalStencil<CELL, double> > diff_stencil;
    //a_stencil will compute its corresponding stencil for all the cells in the domain
    pout()<<setprecision(15)<<scientific;
    a_regApplyBox = a_validBox;
    a_dstVoFs = a_graph.getAllVoFs(a_validBox);
    a_stencil.resize(a_dstVoFs.size());
    diff_stencil.resize(a_dstVoFs.size());
    a_inhoContr.resize(a_dstVoFs.size());

    int radius = 2;
    int powerE = -3;
    bool verbose = 0;
    int exactSolution = 1;
    Vector<int> ptVec;

    //Reading parameters from the inputs file for WLS stencils
    ParmParse pp;
    pp.get("radius"   , radius);
    pp.get("powerE"   , powerE);
    pp.query("verbose"   , verbose);

    //Exact solution
    pp.get("exact_solution",exactSolution);

    if (verbose)
    {
      pp.getarr("pt_verbose",ptVec,0,DIM);
    } 

    string ebbc = a_ebbcName;
    string dombc[2*DIM];

    //Assigning domain and EB boundary conditions
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

    //Looping over all the valid VoFs
    for(int ivof = 0; ivof < a_dstVoFs.size(); ivof++)
    {
      Point pt = a_dstVoFs[ivof].m_pt;
      bool pt_verbose = false;
      if (verbose)
      {
#if DIM==3
	if (ptVec[0]==pt[0] && ptVec[1]==pt[1] && ptVec[2]==pt[2])
#else
	if (ptVec[0]==pt[0] && ptVec[1]==pt[1])
#endif
        {
	  pt_verbose = true;
          pout() << "Debugging Point "<<pt<<endl;
        }
	else
        {
          pt_verbose = false;
        }
      }

      //Creating the WLS fluxes needed for flux diveregence
      ExactVelAdvectionWLSFlux<XFACE, order> xflux(radius, powerE, exactSolution, pt_verbose);
      ExactVelAdvectionWLSFlux<YFACE, order> yflux(radius, powerE, exactSolution, pt_verbose);
      ExactVelAdvectionWLSFlux<ZFACE, order> zflux(radius, powerE, exactSolution, pt_verbose);
      ExactVelAdvectionWLSFlux<BOUNDARY, order> ebflux(radius, powerE, exactSolution, pt_verbose);
      Real inhoTerm = 0.0;
      //Neighbor list type. 0-FaceBased or 1-CellBased
      int neighType = 1;

      //Computing the flux divergence stencil
      WLSStencil<order>::
        getAdvWLSDivFStencil(a_stencil[ivof],
                          a_dstVoFs[ivof],
                          a_graph,
                          a_voludata,
                          a_ebfadata,
                          a_xfacdata, a_yfacdata, a_zfacdata,
                          a_ebnormxdata, a_ebnormydata, 
                          a_ebnormzdata,
                          xflux, yflux, zflux, ebflux, ebbc, dombc,
			  inhoTerm,
                          neighType,
                          a_dx);

#if 0
      //Creating the WLS fluxes needed for flux diveregence
      BaseWLSFlux<XFACE, order> xfluxL(radius, powerE, exactSolution, pt_verbose);
      BaseWLSFlux<YFACE, order> yfluxL(radius, powerE, exactSolution, pt_verbose);
      BaseWLSFlux<ZFACE, order> zfluxL(radius, powerE, exactSolution, pt_verbose);
      BaseWLSFlux<BOUNDARY, order> ebfluxL(radius, powerE, exactSolution, pt_verbose);

      neighType = 1;
      WLSStencil<order>::
        getLaplaWLSDivFStencil(diff_stencil[ivof],
                          a_dstVoFs[ivof],
                          a_graph,
                          a_voludata,
                          a_ebfadata,
                          a_xfacdata, a_yfacdata, a_zfacdata,
                          a_ebnormxdata, a_ebnormydata,
                          a_ebnormzdata,
                          xfluxL, yfluxL, zfluxL, ebfluxL, ebbc, dombc,
                          inhoTerm,
                          neighType,
                          a_dx);
#endif
      RealVect centroid;
      double kappa;
      Normalizor<order>::normalizedVoluMoments(kappa, centroid, a_voludata, a_dstVoFs[ivof], a_graph, a_dx);
      a_stencil[ivof]*= 1.0/(a_dx);
      //diff_stencil[ivof]*= 1.0/(a_dx*a_dx);
      //Real factor= POW(1.0, 1-order);
      //factor *= 1.0;
      //Real diffco = -factor*POW(a_dx, order);
      //Real diffco = -factor;
      //diff_stencil[ivof] *= diffco;
      inhoTerm*= 1.0/(a_dx);
      //a_stencil[ivof]+=diff_stencil[ivof];
      a_inhoContr[ivof] = inhoTerm;
      if (pt_verbose) 
      {
        pout() << "divF stencil for " << a_dstVoFs[ivof].m_pt << ":" << endl;
        a_stencil[ivof].print();
	Real accDivF = 0.0;
	for(int ivec =  0; ivec < a_stencil[ivof].size(); ivec++)
        {
          accDivF += a_stencil[ivof].m_entries[ivec].m_weight;
        }

        pout() << "For phi=1 acc " << accDivF <<" inhoTerm "<<inhoTerm<< endl;
      }
    }
    return true;
  }
  /***/
/////////////end cell to cell


}
#endif
