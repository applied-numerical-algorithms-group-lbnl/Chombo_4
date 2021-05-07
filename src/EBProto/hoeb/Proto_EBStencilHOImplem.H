
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
//    using std::cout;
//    using std::endl;
//    a_regStencil = Stencil<double>::Laplacian();
//    cout << "done with regular sten init" << endl;
//    a_regStencil *= 1.0/(a_dx*a_dx);
//    cout << "done with regular sten scaling" << endl;
    //this stencil has a span so we have to evaluate the stencil on a grown set of the irregular cells so the regular and irregular data actually talk
    //getIrregLocations(a_dstVoFs, a_regStencil, a_validBox, a_graph);
    //cout << "done with getIrregLocations " <<a_dstVoFs.size()<< endl;

    //a_stencil will compute its corresponding stencil for all the cells in the domain
    a_regApplyBox = a_validBox;
    a_dstVoFs = a_graph.getAllVoFs(a_validBox);
    a_stencil.resize(a_dstVoFs.size());
    a_inhoContr.resize(a_dstVoFs.size());

    int radius = 2;
    int powerE = -3;
    int num_mono = 0;
    bool verbose = 0;
    Vector<int> ptVec;
    Vector<int> phiPowerX;
    Vector<int> phiPowerY;
    Vector<int> phiPowerZ;
    Vector<Real> phiCoef;

    //Reading parameters from the inputs file for WLS stencils
    ParmParse pp;
    pp.get("radius"   , radius);
    pp.get("powerE"   , powerE);
    pp.query("verbose"   , verbose);
    pp.get("phi_num_mono"   , num_mono);
    pp.getarr("phi_power_x",phiPowerX,0,num_mono);
    pp.getarr("phi_power_y",phiPowerY,0,num_mono);
    pp.getarr("phi_power_z",phiPowerZ,0,num_mono);
    pp.getarr("phi_coef",phiCoef,0,num_mono);

    vector<pair<Point, Real> > entries;
    for (int i=0; i<num_mono; i++)
    {
      std::pair<Point, Real> mono;
      mono.first = Point::Zeros();
      mono.first[0] = phiPowerX[i];
      mono.first[1] = phiPowerY[i];
#if DIM==3
      mono.first[2] = phiPowerZ[i];
#endif
      mono.second = phiCoef[i];
      entries.push_back(mono);
    }

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
	if (ptVec[0]==pt[0] && ptVec[1]==pt[1])
        {
	  pt_verbose = true;
          cout << "Debugging Point "<<pt<<endl;
        }
	else
        {
          pt_verbose = false;
        }
      }

      //Creating the WLS fluxes needed for flux diveregence
      BaseWLSFlux<XFACE, order> xflux(radius, powerE, entries, pt_verbose);
      BaseWLSFlux<YFACE, order> yflux(radius, powerE, entries, pt_verbose);
      BaseWLSFlux<ZFACE, order> zflux(radius, powerE, entries, pt_verbose);
      BaseWLSFlux<BOUNDARY, order> ebflux(radius, powerE, entries, pt_verbose);
      Real inhoTerm = 0.0;

      //Computing the flux divergence stencil
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
			  inhoTerm,
                          a_dx);

      RealVect centroid;
      double kappa;
      Normalizor<order>::normalizedVoluMoments(kappa, centroid, a_voludata, a_dstVoFs[ivof], a_graph, a_dx);
      a_stencil[ivof]*= 1.0/(a_dx*a_dx);
      inhoTerm*= 1.0/(a_dx*a_dx);
      a_inhoContr[ivof] = inhoTerm;
      if (pt_verbose) 
      {
        cout << "divF stencil for " << a_dstVoFs[ivof].m_pt << ":" << endl;
        a_stencil[ivof].print();
	Real accDivF = 0.0;
	for(int ivec =  0; ivec < a_stencil[ivof].size(); ivec++)
        {
          accDivF += a_stencil[ivof].m_entries[ivec].m_weight;
        }

        cout << "For phi=1 acc " << accDivF <<" inhoTerm "<<inhoTerm<< endl;
      }
    }
    return true;
  }
  /***/
/////////////end cell to cell

}
#endif
