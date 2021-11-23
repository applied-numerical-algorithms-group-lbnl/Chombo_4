
#ifndef _Proto_EBStencilHOImplem_H_
#define _Proto_EBStencilHOImplem_H_
#include "Chombo_ParmParse.H"
#include "../implem/Proto_Normalizor.H"
#include "Proto_WLSCell.H"

namespace Proto
{

  template <int order>
  vector<EBIndex<CELL> >
  EBStencilHO<CELL, CELL,  order, double>::
  getFlapVoFs(const EBGraph& a_graph,  EBIndex<CELL>& a_start)
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

  template <int order>
  bool
  EBStencilHO<CELL, CELL,  order, double>::
  isRegularStencil(const EBIndex<CELL> & a_vof,
                   const int             a_radius,
                   const EBGraph       & a_graph,
                   const VoluData      & a_voludata,
                   const double          a_dx)
  {
    PR_TIME("isRegularStencil");
    bool retval;
    Vector<EBIndex<CELL> > vof_cont;
    vof_cont.push_back(a_vof);
    RealVect centroid;
    double kappa;
    Normalizor<order>::normalizedVoluMoments(kappa, centroid, a_voludata, a_vof, a_graph, a_dx);
    if (kappa < 1.0) retval=false;  
    int startV = 0;
    for (int irad = 1; irad<=a_radius; irad++)
    {
      int neighSize = vof_cont.size();
      for (int ineigh = startV; ineigh<neighSize; ineigh++){
        vector< EBIndex<CELL> > volsinclusive=  getFlapVoFs(a_graph, vof_cont[ineigh]);
        for(int ivof = 0; ivof < volsinclusive.size(); ivof++)
        {
          bool incl = true;
          for(int ivof2 = 0; ivof2 < vof_cont.size(); ivof2++)
          {
            if (volsinclusive[ivof]==vof_cont[ivof2])
            {
              incl=false;
              break;
            }
          }
          if (incl) 
          {
            vof_cont.push_back(volsinclusive[ivof]);
            Normalizor<order>::normalizedVoluMoments(kappa, centroid, a_voludata, volsinclusive[ivof], a_graph, a_dx);
            if (kappa < 1.0) retval=false;  
          }
        }
      }
      startV = neighSize;
    }
    return retval;
  }

  template <int order>
  void
  EBStencilHO<CELL, CELL,  order, double>::
  printWLSStencil(int                                  radius,
                  const LocalStencil< CELL,  double> & a_stencil,
                  const Point                        & a_vofPt) 
  {
    using Chombo4::pout;
    pout()<<setprecision(16)<<scientific;
    for (int j = radius+1; j>=-(radius+1); j--)
    {
      for (int i = -(radius+1); i<=radius+1; i++)
      {
        bool act = false;
        for(int ivec =  0; ivec < a_stencil.size(); ivec++)
        {
          if (a_vofPt[0]+i==a_stencil.m_entries[ivec].m_vof.m_pt[0] &&
              a_vofPt[1]+j==a_stencil.m_entries[ivec].m_vof.m_pt[1])
          {
            pout()<<"("<<a_stencil.m_entries[ivec].m_vof.m_pt[0]<<","<<a_stencil.m_entries[ivec].m_vof.m_pt[1]<<") ";                
            act = true;
          }
        }
        if (!act) pout()<<"("<<a_vofPt[0]+i<<","<<a_vofPt[1]+j<<") ";
      }
      pout()<<endl;
    }
    //pout()    << std::setprecision(1);
    for (int j = radius+1; j>=-(radius+1); j--)
    {
      for (int i = -(radius+1); i<=radius+1; i++)
      { 
        bool act = false;
        for(int ivec =  0; ivec < a_stencil.size(); ivec++)
        {
          if (a_vofPt[0]+i==a_stencil.m_entries[ivec].m_vof.m_pt[0] &&
              a_vofPt[1]+j==a_stencil.m_entries[ivec].m_vof.m_pt[1])
          {
            if(a_stencil.m_entries[ivec].m_weight<0.0)
            {
              pout()<<"-"<<abs(a_stencil.m_entries[ivec].m_weight)<<"  ";
            }
            else
            {                
              pout()<<" "<<abs(a_stencil.m_entries[ivec].m_weight)<<"  ";
            }
            act = true;
          } 
        }
        if (!act) pout()<<"          N/A        "<<"  ";
      }
      pout()<<endl;
    }
    pout()    << std::setprecision(16);
  }

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
             Point                                       a_dstGrown,
             const shared_ptr<BaseExactSolution<order> > & a_boundarySolutions)
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
                             forceNeumann,
                             a_boundarySolutions);
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
                             forceNeumann,
                             a_boundarySolutions);
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
                       bool  a_forceNeumann,
                       const shared_ptr<BaseExactSolution<order> > & a_boundarySolutions)
  {
    using Chombo4::pout;
    pout()<<setprecision(15)<<scientific;
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
    //Neighbor list type. 0-FaceBased or 1-CellBased
    int neighListType = 0;

    Vector<int> ptVec;

    //Reading parameters from the inputs file for WLS stencils
    ParmParse pp;
    pp.get("radius"   , radius);
    pp.get("powerE"   , powerE);
    pp.query("verbose"   , verbose);
    pp.query("neighListType"   , neighListType);

    //Exact solution
    if (a_boundarySolutions == nullptr)
      {
        pp.get("exact_solution",exactSolution);
      }
    else
      {
        pout() << "WLSPoissionStencil using externally specified BC" << std::endl;
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

    int ncells = 32;
    int nrepRows = -1;
    //pp.get("nx"   , ncells);
    Point lodom = a_domain.low();
    Point hidom = a_domain.high();
    ncells = hidom[1] - lodom[1] + 1;
    pp.query("nrepRows"   , nrepRows);
    vector<map<int,int > > map1;
    map<pair<int,int>,EBIndex<CELL> > map2;
    map1.resize(ncells);

    for(int k=0; k<ncells; k++)
    {
      //Looping over all the valid VoFs
      for(int ivof = 0; ivof < a_dstVoFs.size(); ivof++)
      {
        Point pt = a_dstVoFs[ivof].m_pt;
        pair<int,int> ptPair (pt[0],pt[1]);
        map2[ptPair] = a_dstVoFs[ivof];
        if (pt[1]==k)
        {
          map1[k][pt[0]] = ivof;
        }
      }
    }

    map<int,int >::iterator it;
    map<pair<int,int>,EBIndex<CELL> >::iterator it2;
    vector<LocalStencil<CELL, double> > row_stencil;
    vector<double> row_inhoContr;
    vector<EBIndex<CELL> > row_index;

    int initRep = nrepRows;
    while(initRep<2*radius)
    {
      initRep += nrepRows;
    }

    if (nrepRows == -1)
    {
      initRep = 0;
      nrepRows = ncells;
    }

    //Looping over all the valid VoFs
    for(int k = initRep; k<initRep + nrepRows; k++)
    {
      for (it = map1[k].begin(); it != map1[k].end(); it++)
      {
      //for(int ivof = 0; ivof < a_dstVoFs.size(); ivof++)
      //{
        //EBIndex<CELL> vofInd = a_dstVoFs[ivof];
        int ivof = it->second;
        EBIndex<CELL> vofInd = a_dstVoFs[ivof];
        Point pt = vofInd.m_pt;
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
        // BaseWLSFlux<XFACE, order> xflux(radius, powerE, exactSolution, pt_verbose);
        // BaseWLSFlux<YFACE, order> yflux(radius, powerE, exactSolution, pt_verbose);
        // BaseWLSFlux<ZFACE, order> zflux(radius, powerE, exactSolution, pt_verbose);
        // BaseWLSFlux<BOUNDARY, order> ebflux(radius, powerE, exactSolution, pt_verbose);
        //Creating the WLS fluxes needed for flux divergence - U on a face
        BaseWLSFlux<XFACE, order> xflux;
        BaseWLSFlux<YFACE, order> yflux;
        BaseWLSFlux<ZFACE, order> zflux;
        BaseWLSFlux<BOUNDARY, order> ebflux;
        if (a_boundarySolutions == nullptr)
          {
            //pout() << "Using uniform div boundary values" << std::endl;
            xflux = BaseWLSFlux<XFACE, order>(radius, powerE, exactSolution, pt_verbose);
            yflux = BaseWLSFlux<YFACE, order>(radius, powerE, exactSolution, pt_verbose);
            zflux = BaseWLSFlux<ZFACE, order>(radius, powerE, exactSolution, pt_verbose);
            ebflux = BaseWLSFlux<BOUNDARY, order>(radius, powerE, exactSolution, pt_verbose);
          }
        else
          {
            //pout() << "Using specified div boundary values" << std::endl;
            xflux = BaseWLSFlux<XFACE, order>(radius, powerE, a_boundarySolutions, 0, pt_verbose);
            yflux = BaseWLSFlux<YFACE, order>(radius, powerE, a_boundarySolutions, 0, pt_verbose);
            zflux = BaseWLSFlux<ZFACE, order>(radius, powerE, a_boundarySolutions, 0, pt_verbose);
            ebflux = BaseWLSFlux<BOUNDARY, order>(radius, powerE, a_boundarySolutions, 0, pt_verbose);
          }

        Real inhoTerm = 0.0;
  
        //Computing the flux divergence stencil
        WLSStencil<order> lapSten;
        lapSten.
          getLaplaWLSDivFStencil(a_stencil[ivof],
                            vofInd,
                            a_graph,
                            a_voludata,
                            a_ebfadata,
                            a_xfacdata, a_yfacdata, a_zfacdata,
                            a_ebnormxdata, a_ebnormydata, 
                            a_ebnormzdata,
                            xflux, yflux, zflux, ebflux, ebbc, dombc,
                            inhoTerm,
                            neighListType,
                            a_dx);
  
        RealVect centroid;
        double kappa;
        Normalizor<order>::normalizedVoluMoments(kappa, centroid, a_voludata, vofInd, a_graph, a_dx);
        a_stencil[ivof]*= 1.0/(a_dx*a_dx);
        inhoTerm*= 1.0/(a_dx*a_dx);
        a_inhoContr[ivof] = inhoTerm;
        row_stencil.push_back(a_stencil[ivof]);
        row_inhoContr.push_back(a_inhoContr[ivof]);
        row_index.push_back(vofInd);
        if (pt_verbose) 
        {
          pout() << "divF stencil for " << vofInd.m_pt << ":" << endl;
          //a_stencil[ivof].print();
          printWLSStencil(radius, a_stencil[ivof], vofInd.m_pt);
          Real accDivF = 0.0;
          for(int ivec =  0; ivec < a_stencil[ivof].size(); ivec++)
          {
            accDivF += a_stencil[ivof].m_entries[ivec].m_weight;
          }
  
          //pout() << "For phi=1 acc " << accDivF <<" inhoTerm "<<inhoTerm<< endl;
          pout() << "Sum of all stencil weights " << accDivF <<" stencil size "<<a_stencil[ivof].size()<< endl;
        }
      }
    }

    if (nrepRows != ncells)
    {
      int kk =0;
      for(int k = 0; k<map1.size(); k++)
      {
        for (it = map1[k].begin(); it != map1[k].end(); it++)
        {
          if (k>=initRep && k<initRep + nrepRows) continue;
          int ivof = it->second;
          if (kk==row_stencil.size()) kk=0;
          EBIndex<CELL> vofInd = a_dstVoFs[ivof];
          int difi = vofInd.m_pt[0] - row_index[kk].m_pt[0];
          int difj = vofInd.m_pt[1] - row_index[kk].m_pt[1];
          LocalStencil<CELL, double> new_stencil;
          for(int ivec =  0; ivec < row_stencil[kk].size(); ivec++)
          {
              int pt_x = row_stencil[kk].m_entries[ivec].m_vof.m_pt[0] + difi;
              int pt_y = row_stencil[kk].m_entries[ivec].m_vof.m_pt[1] + difj;
              if (pt_y < 0) {
                pt_y += ncells;
                pt_x += ncells/nrepRows;
                if (difi==0) pt_x -= ncells/nrepRows;
              }
              if (pt_y > ncells-1) {
                pt_y -= ncells;
                pt_x -= ncells/nrepRows;
                if (difi==0) pt_x += ncells/nrepRows;
              }
              pair<int,int> pairPt (pt_x,pt_y);
              it2 = map2.find(pairPt);
              if (it2 != map2.end())
              {
                new_stencil.add(it2->second, row_stencil[kk].m_entries[ivec].m_weight);
              }
          }
          Point pt = vofInd.m_pt;
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
 
          a_stencil[ivof] = new_stencil; 
          a_inhoContr[ivof] = row_inhoContr[kk]; 
          kk++;

          if (pt_verbose) 
          {
            pout() << "divF stencil for " << vofInd.m_pt << ":" << endl;
            //a_stencil[ivof].print();
            printWLSStencil(radius, a_stencil[ivof], vofInd.m_pt);
            Real accDivF = 0.0;
            for(int ivec =  0; ivec < a_stencil[ivof].size(); ivec++)
            {
              accDivF += a_stencil[ivof].m_entries[ivec].m_weight;
            }
    
            pout() << "Original pt " << row_index[kk-1].m_pt[0] <<" , "<<row_index[kk-1].m_pt[1]<< endl;
            pout() << "Sum of all stencil weights " << accDivF <<" stencil size "<<a_stencil[ivof].size()<< endl;
          }
 

        }
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
    //a_stencil will compute its corresponding stencil for all the cells in the domain
    using Chombo4::pout;
    pout()<<setprecision(15)<<scientific;
    a_regApplyBox = a_validBox;
    a_dstVoFs = a_graph.getAllVoFs(a_validBox);
    a_stencil.resize(a_dstVoFs.size());
    a_inhoContr.resize(a_dstVoFs.size());

    int radius = 2;
    int powerE = -3;
    bool verbose = 0;
    int exactSolution = 1;
    //Neighbor list type. 0-FaceBased or 1-CellBased
    int neighListType = 1;

    Vector<int> ptVec;

    //Reading parameters from the inputs file for WLS stencils
    ParmParse pp;
    pp.get("radius"   , radius);
    pp.get("powerE"   , powerE);
    pp.query("verbose"   , verbose);
    pp.query("neighListType"   , neighListType);

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

    int ncells = 32;
    int nrepRows = -1;
    //pp.get("nx"   , ncells);
    Point lodom = a_domain.low();
    Point hidom = a_domain.high();
    ncells = hidom[1] - lodom[1] + 1;
    pp.query("nrepRows"   , nrepRows);
    //vector<map<int,int > > map1;
    vector<map<pair<int,int>,int > > map1;
    map<tuple<int,int,int>,EBIndex<CELL> > map2;
    map1.resize(ncells);

    for(int k=0; k<ncells; k++)
    {
      //Looping over all the valid VoFs
      for(int ivof = 0; ivof < a_dstVoFs.size(); ivof++)
      {
        Point pt = a_dstVoFs[ivof].m_pt;
#if DIM==3
        tuple<int,int,int> ptPair (pt[0],pt[1],pt[2]);
#elif DIM==2
        tuple<int,int,int> ptPair (pt[0],pt[1],0);
#endif
        map2[ptPair] = a_dstVoFs[ivof];
        if (pt[1]==k)
        {
#if DIM==3
          pair<int,int> ptPair2 (pt[0],pt[2]);
#elif DIM==2
          pair<int,int> ptPair2 (pt[0],0);
#endif
          map1[k][ptPair2] = ivof;
        }
      }
    }

    map<pair<int,int>,int >::iterator it;
    map<tuple<int,int,int>,EBIndex<CELL> >::iterator it2;
    vector<LocalStencil<CELL, double> > row_stencil;
    vector<double> row_inhoContr;
    vector<EBIndex<CELL> > row_index;

    int initRep = nrepRows;
    while(initRep<2*radius)
    {
      initRep += nrepRows;
    }

    if (nrepRows == -1)
    {
      initRep = 0;
      nrepRows = ncells;
    }

    bool storeRegSten = false;
    LocalStencil<CELL, double> regSten;
    Real regInhContr;
    //Looping over all the valid VoFs
    for(int k = initRep; k<initRep + nrepRows; k++)
    {
      for (it = map1[k].begin(); it != map1[k].end(); it++)
      {
      //for(int ivof = 0; ivof < a_dstVoFs.size(); ivof++)
      //{
        //EBIndex<CELL> vofInd = a_dstVoFs[ivof];
        int ivof = it->second;
        EBIndex<CELL> vofInd = a_dstVoFs[ivof];
        bool isReg = isRegularStencil(vofInd, radius, a_graph, a_voludata, a_dx);
        Point pt = vofInd.m_pt;
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

        //if (!isReg || (isReg && !storeRegSten))
        if (1)
        {
        //Creating the WLS fluxes needed for flux diveregence
        ExactVelAdvectionWLSFlux<XFACE, order> xflux(radius, powerE, exactSolution, pt_verbose);
        ExactVelAdvectionWLSFlux<YFACE, order> yflux(radius, powerE, exactSolution, pt_verbose);
        ExactVelAdvectionWLSFlux<ZFACE, order> zflux(radius, powerE, exactSolution, pt_verbose);
        ExactVelAdvectionWLSFlux<BOUNDARY, order> ebflux(radius, powerE, exactSolution, pt_verbose);
        Real inhoTerm = 0.0;
  
        //Computing the flux divergence stencil
        WLSStencil<order> advSten;
        advSten.
          getAdvWLSDivFStencil(a_stencil[ivof],
                            vofInd,
                            a_graph,
                            a_voludata,
                            a_ebfadata,
                            a_xfacdata, a_yfacdata, a_zfacdata,
                            a_ebnormxdata, a_ebnormydata, 
                            a_ebnormzdata,
                            xflux, yflux, zflux, ebflux, ebbc, dombc,
                            inhoTerm,
                            neighListType,
                            a_dx);
  
        RealVect centroid;
        double kappa;
        Normalizor<order>::normalizedVoluMoments(kappa, centroid, a_voludata, vofInd, a_graph, a_dx);
        a_stencil[ivof]*= 1.0/(kappa*a_dx);
        inhoTerm*= 1.0/(a_dx);
        //a_inhoContr[ivof] = inhoTerm;
        a_inhoContr[ivof] = 0.0;
        if (isReg)
        {
          regSten = a_stencil[ivof];
          regInhContr = a_inhoContr[ivof];
          storeRegSten = true;  
        }
        }
        else
        {
          a_stencil[ivof] = regSten;
          a_inhoContr[ivof] = regInhContr;
          int difi = vofInd.m_pt[0] - a_stencil[ivof].m_entries[0].m_vof.m_pt[0];
          int difj = vofInd.m_pt[1] - a_stencil[ivof].m_entries[0].m_vof.m_pt[1];
#if DIM==3
          int difk = vofInd.m_pt[2] - a_stencil[ivof].m_entries[0].m_vof.m_pt[2];
#endif
   	  for(int ivec =  0; ivec < a_stencil[ivof].size(); ivec++)
          {
            a_stencil[ivof].m_entries[ivec].m_vof.m_pt[0] += difi;
            a_stencil[ivof].m_entries[ivec].m_vof.m_pt[1] += difj;
#if DIM==3
            a_stencil[ivof].m_entries[ivec].m_vof.m_pt[2] += difk;
#endif
          }
        }
        row_stencil.push_back(a_stencil[ivof]);
        row_inhoContr.push_back(a_inhoContr[ivof]);
        row_index.push_back(vofInd);
        if (pt_verbose) 
        {
          pout() << "divF stencil for " << vofInd.m_pt << ":" << endl;
#if DIM==3
          a_stencil[ivof].print();
#elif DIM==2
          printWLSStencil(radius, a_stencil[ivof], vofInd.m_pt);
#endif
          Real accDivF = 0.0;
  	  for(int ivec =  0; ivec < a_stencil[ivof].size(); ivec++)
          {
            accDivF += a_stencil[ivof].m_entries[ivec].m_weight;
          }
  
          //pout() << "For phi=1 acc " << accDivF <<" inhoTerm "<<inhoTerm<< endl;
          pout() << "Sum of all stencil weights " << accDivF <<" stencil size "<<a_stencil[ivof].size()<< endl;
        }
      }
    }

    if (nrepRows != ncells)
    {
      int kk =0;
      for(int k = 0; k<map1.size(); k++)
      {
        for (it = map1[k].begin(); it != map1[k].end(); it++)
        {
          if (k>=initRep && k<initRep + nrepRows) continue;
          int ivof = it->second;
          if (kk==row_stencil.size()) kk=0;
          EBIndex<CELL> vofInd = a_dstVoFs[ivof];
          int difi = vofInd.m_pt[0] - row_index[kk].m_pt[0];
          int difj = vofInd.m_pt[1] - row_index[kk].m_pt[1];
#if DIM==3
          int difk = vofInd.m_pt[2] - row_index[kk].m_pt[2];
#endif
          LocalStencil<CELL, double> new_stencil;
          for(int ivec =  0; ivec < row_stencil[kk].size(); ivec++)
          {
              int pt_x = row_stencil[kk].m_entries[ivec].m_vof.m_pt[0] + difi;
              int pt_y = row_stencil[kk].m_entries[ivec].m_vof.m_pt[1] + difj;
#if DIM==3
              int pt_z = row_stencil[kk].m_entries[ivec].m_vof.m_pt[2] + difk;
#endif
              if (pt_y < 0) {
                pt_y += ncells;
                pt_x += ncells/nrepRows;
                if (difi==0) pt_x -= ncells/nrepRows;
              }
              if (pt_y > ncells-1) {
                pt_y -= ncells;
                pt_x -= ncells/nrepRows;
                if (difi==0) pt_x += ncells/nrepRows;
              }
#if DIM==3
              tuple<int,int,int> pairPt (pt_x,pt_y,pt_z);
#elif DIM==2
              tuple<int,int,int> pairPt (pt_x,pt_y,0);
#endif
              it2 = map2.find(pairPt);
              if (it2 != map2.end())
              {
                new_stencil.add(it2->second, row_stencil[kk].m_entries[ivec].m_weight);
              }
          }
          Point pt = vofInd.m_pt;
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
 
          a_stencil[ivof] = new_stencil; 
          a_inhoContr[ivof] = row_inhoContr[kk]; 
          kk++;

          if (pt_verbose) 
          {
            pout() << "divF stencil for " << vofInd.m_pt << ":" << endl;
#if DIM==3
            a_stencil[ivof].print();
#elif DIM==2
            printWLSStencil(radius, a_stencil[ivof], vofInd.m_pt);
#endif
            Real accDivF = 0.0;
            for(int ivec =  0; ivec < a_stencil[ivof].size(); ivec++)
            {
              accDivF += a_stencil[ivof].m_entries[ivec].m_weight;
            }
    
            pout() << "Original pt " << row_index[kk-1].m_pt[0] <<" , "<<row_index[kk-1].m_pt[1]<< endl;
            pout() << "Sum of all stencil weights " << accDivF <<" stencil size "<<a_stencil[ivof].size()<< endl;
          }
 

        }
      }
    }


    return true;
  }
  /***/
/////////////end cell to cell


  template <int order>
  bool
  EBStencilHO<CELL, CELL,  order, double>::
  getWLSDivergenceCompStencil(vector<EBIndex<CELL> >                    & a_dstVoFs,
                              vector<vector<LocalStencil<CELL, double> > >      & a_stencil,
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
                              bool  a_forceNeumann,
                              const vector<shared_ptr<BaseExactSolution<order> > >& a_boundarySolutions)
  {
    //a_stencil will compute its corresponding stencil for all the cells in the domain
    pout()<<setprecision(15)<<scientific;
    a_regApplyBox = a_validBox;
    a_dstVoFs = a_graph.getAllVoFs(a_validBox);
    a_stencil.resize(SpaceDim, vector<LocalStencil<CELL, double> >(a_dstVoFs.size()));
    a_inhoContr.resize(a_dstVoFs.size());

    int radius = 2;
    int powerE = -3;
    bool verbose = 0;
    int exactSolution = 1;
    //Neighbor list type. 0-FaceBased or 1-CellBased
    int neighListType = 1;

    Vector<int> ptVec;

    //Reading parameters from the inputs file for WLS stencils
    ParmParse pp;
    pp.get("radius"   , radius);
    pp.get("powerE"   , powerE);
    pp.query("verbose"   , verbose);
    pp.query("neighListType"   , neighListType);

    //Exact solution
    if (a_boundarySolutions.size() == 0)
      {
        pp.get("exact_solution",exactSolution);
      }
    else
      {
        pout() << "WLSDivergenceCompStencil using externally specified BC" << std::endl;
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

    // Setup periodic boundaries
    int ncells = 32;
    int nrepRows = -1;
    //pp.get("nx"   , ncells);
    ncells = 1.0/a_dx;
    pp.query("nrepRows", nrepRows);
    vector<map<int, int> > map1;
    map<pair<int,int>, EBIndex<CELL> > map2;
    map1.resize(ncells);
    for(int k=0; k<ncells; k++)
      {
        //Looping over all the valid VoFs
        for(int ivof = 0; ivof < a_dstVoFs.size(); ivof++)
          {
            Point pt = a_dstVoFs[ivof].m_pt;
            pair<int,int> ptPair (pt[0],pt[1]);
            map2[ptPair] = a_dstVoFs[ivof];
            if (pt[1]==k)
              {
                map1[k][pt[0]] = ivof;
              }
          }
      }
    map<int,int >::iterator it;
    map<pair<int,int>,EBIndex<CELL> >::iterator it2;
    vector<LocalStencil<CELL, double> > row_stencil;
    vector<double> row_inhoContr;
    vector<EBIndex<CELL> > row_index;
    int initRep = nrepRows;
    while(initRep<2*radius)
      {
        initRep += nrepRows;
      }

    if (nrepRows == -1)
      {
        initRep = 0;
        nrepRows = ncells;
      }

    //Looping over all the valid VoFs and setting up their stencils
    for(int k = initRep; k<initRep + nrepRows; k++)
      {
        for (it = map1[k].begin(); it != map1[k].end(); it++)
          {
            int ivof = it->second;
            EBIndex<CELL> vofInd = a_dstVoFs[ivof];
            Point pt = vofInd.m_pt;
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

            //Creating the WLS fluxes needed for flux divergence - U on a face
            BaseWLSFlux<XFACE, order> xflux;
            BaseWLSFlux<YFACE, order> yflux;
            BaseWLSFlux<ZFACE, order> zflux;
            BaseWLSFlux<BOUNDARY, order> ebflux;
            if (a_boundarySolutions.size() == 0)
              {
                //pout() << "Using uniform div boundary values" << std::endl;
                xflux = BaseWLSFlux<XFACE, order>(radius, powerE, exactSolution, pt_verbose);
                yflux = BaseWLSFlux<YFACE, order>(radius, powerE, exactSolution, pt_verbose);
                zflux = BaseWLSFlux<ZFACE, order>(radius, powerE, exactSolution, pt_verbose);
                ebflux = BaseWLSFlux<BOUNDARY, order>(radius, powerE, exactSolution, pt_verbose);
              }
            else
              {
                //pout() << "Using specified div boundary values" << std::endl;
                CH_assert(a_boundarySolutions.size() == 4);
                xflux = BaseWLSFlux<XFACE, order>(radius, powerE, a_boundarySolutions[0], 0, pt_verbose);
                yflux = BaseWLSFlux<YFACE, order>(radius, powerE, a_boundarySolutions[1], 0, pt_verbose);
                zflux = BaseWLSFlux<ZFACE, order>(radius, powerE, a_boundarySolutions[2], 0, pt_verbose);
                ebflux = BaseWLSFlux<BOUNDARY, order>(radius, powerE, a_boundarySolutions[3], 0, pt_verbose);
              }

            //Computing the flux divergence stencil
            Real inhoTerm = 0.0;
            WLSStencil<order> sten;
            sten.getWLSDivFCompStencil(a_stencil[0][ivof],
                                       a_stencil[1][ivof],
                                       a_stencil[1][ivof], // dummy for the z
                                       vofInd,
                                       a_graph,
                                       a_voludata,
                                       a_ebfadata,
                                       a_xfacdata, a_yfacdata, a_zfacdata,
                                       a_ebnormxdata, a_ebnormydata,
                                       a_ebnormzdata,
                                       xflux, yflux, zflux, ebflux, ebbc, dombc,
                                       inhoTerm,
                                       neighListType,
                                       a_dx);

            RealVect centroid;
            double kappa;
            Normalizor<order>::normalizedVoluMoments(kappa, centroid, a_voludata, vofInd, a_graph, a_dx);
            // a_stencil[ivof] *= 1.0/(a_dx);
            for (int dir=0; dir!=SpaceDim; dir++)
              {
                a_stencil[dir][ivof] *= 1.0/(a_dx);
              }
            a_inhoContr[ivof] = inhoTerm/(a_dx);
            row_stencil.push_back(a_stencil[0][ivof]); // FIXME for multi-comp
            row_inhoContr.push_back(a_inhoContr[ivof]);
            row_index.push_back(vofInd);

            if (pt_verbose)
              {
                pout() << "divF stencil for " << vofInd.m_pt << ":" << endl;
                pout() << " inhomg term " << a_inhoContr[ivof] << endl;
                printWLSStencil(radius, a_stencil[0][ivof], vofInd.m_pt);
                Real accDivF = 0.0;
                for(int ivec =  0; ivec < a_stencil[0][ivof].size(); ivec++)
                  {
                    accDivF += a_stencil[0][ivof].m_entries[ivec].m_weight;
                  }

                //pout() << "For phi=1 acc " << accDivF <<" inhoTerm "<<inhoTerm<< endl;
                pout() << "Sum of all stencil weights " << accDivF <<" stencil size "<<a_stencil[0][ivof].size()<< endl;
              }
          }
      }

//     // update stencils for periodic boundaries
//     if (nrepRows != ncells)
//       {
//         int kk =0;
//         for(int k = 0; k<map1.size(); k++)
//           {
//             for (it = map1[k].begin(); it != map1[k].end(); it++)
//               {
//                 if (k>=initRep && k<initRep + nrepRows) continue;
//                 int ivof = it->second;
//                 if (kk==row_stencil.size()) kk=0;
//                 EBIndex<CELL> vofInd = a_dstVoFs[ivof];
//                 int difi = vofInd.m_pt[0] - row_index[kk].m_pt[0];
//                 int difj = vofInd.m_pt[1] - row_index[kk].m_pt[1];
//                 LocalStencil<CELL, double> new_stencil;
//                 for(int ivec =  0; ivec < row_stencil[kk].size(); ivec++)
//                   {
//                     int pt_x = row_stencil[kk].m_entries[ivec].m_vof.m_pt[0] + difi;
//                     int pt_y = row_stencil[kk].m_entries[ivec].m_vof.m_pt[1] + difj;
//                     if (pt_y < 0) {
//                       pt_y += ncells;
//                       pt_x += ncells/nrepRows;
//                       if (difi==0) pt_x -= ncells/nrepRows;
//                     }
//                     if (pt_y > ncells-1) {
//                       pt_y -= ncells;
//                       pt_x -= ncells/nrepRows;
//                       if (difi==0) pt_x += ncells/nrepRows;
//                     }
//                     pair<int,int> pairPt (pt_x,pt_y);
//                     it2 = map2.find(pairPt);
//                     if (it2 != map2.end())
//                       {
//                         new_stencil.add(it2->second, row_stencil[kk].m_entries[ivec].m_weight);
//                       }
//                   }
//                 Point pt = vofInd.m_pt;
//                 bool pt_verbose = false;
//                 if (verbose)
//                   {
// #if DIM==3
//                     if (ptVec[0]==pt[0] && ptVec[1]==pt[1] && ptVec[2]==pt[2])
// #else
//                       if (ptVec[0]==pt[0] && ptVec[1]==pt[1])
// #endif
//                         {
//                           pt_verbose = true;
//                           pout() << "Debugging Point "<<pt<<endl;
//                         }
//                       else
//                         {
//                           pt_verbose = false;
//                         }
//                   }

//                 a_stencil[ivof] = new_stencil;
//                 a_inhoContr[ivof] = row_inhoContr[kk];
//                 kk++;

//                 if (pt_verbose)
//                   {
//                     pout() << "divF stencil for " << vofInd.m_pt << ":" << endl;
//                     printWLSStencil(radius, a_stencil[ivof][0], vofInd.m_pt);
//                     Real accDivF = 0.0;
//                     for(int ivec =  0; ivec < a_stencil[ivof][0].size(); ivec++)
//                       {
//                         accDivF += a_stencil[ivof][0].m_entries[ivec].m_weight;
//                       }

//                     pout() << "Original pt " << row_index[kk-1].m_pt[0] <<" , "<<row_index[kk-1].m_pt[1]<< endl;
//                     pout() << "Sum of all stencil weights " << accDivF <<" stencil size "<<a_stencil[ivof][0].size()<< endl;
//                   }


//               }
//           }
//       }


    return true;
  }


  template <int order>
  bool
  EBStencilHO<CELL, CELL,  order, double>::
  getWLSGradientStencil(vector<EBIndex<CELL> >                    & a_dstVoFs,
                              vector<vector<LocalStencil<CELL, double> > >      & a_stencil,
                              vector<vector<double> >                          & a_inhoContr,
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
                              bool  a_forceNeumann,
                              const shared_ptr<BaseExactSolution<order> >& a_boundarySolutions)
  {
    //a_stencil will compute its corresponding stencil for all the cells in the domain
    pout()<<setprecision(15)<<scientific;
    a_regApplyBox = a_validBox;
    a_dstVoFs = a_graph.getAllVoFs(a_validBox);
    a_stencil.resize(SpaceDim, vector<LocalStencil<CELL, double> >(a_dstVoFs.size()));
    a_inhoContr.resize(SpaceDim, vector<double>(a_dstVoFs.size()));

    int radius = 2;
    int powerE = -3;
    bool verbose = 0;
    int exactSolution = 1;
    //Neighbor list type. 0-FaceBased or 1-CellBased
    int neighListType = 1;

    Vector<int> ptVec;

    //Reading parameters from the inputs file for WLS stencils
    ParmParse pp;
    pp.get("radius"   , radius);
    pp.get("powerE"   , powerE);
    pp.query("verbose"   , verbose);
    pp.query("neighListType"   , neighListType);

    //Exact solution
    if (a_boundarySolutions == nullptr)
      {
        pp.get("exact_solution",exactSolution);
      }
    else
      {
        pout() << "WLSGradientStencil using externally specified BC" << std::endl;
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

    // Setup periodic boundaries
    int ncells = 32;
    int nrepRows = -1;
    //pp.get("nx"   , ncells);
    ncells = 1.0/a_dx;
    pp.query("nrepRows", nrepRows);
    vector<map<int, int> > map1;
    map<pair<int,int>, EBIndex<CELL> > map2;
    map1.resize(ncells);
    for(int k=0; k<ncells; k++)
      {
        //Looping over all the valid VoFs
        for(int ivof = 0; ivof < a_dstVoFs.size(); ivof++)
          {
            Point pt = a_dstVoFs[ivof].m_pt;
            pair<int,int> ptPair (pt[0],pt[1]);
            map2[ptPair] = a_dstVoFs[ivof];
            if (pt[1]==k)
              {
                map1[k][pt[0]] = ivof;
              }
          }
      }
    map<int,int >::iterator it;
    map<pair<int,int>,EBIndex<CELL> >::iterator it2;
    vector<LocalStencil<CELL, double> > row_stencil;
    vector<double> row_inhoContr;
    vector<EBIndex<CELL> > row_index;
    int initRep = nrepRows;
    while(initRep<2*radius)
      {
        initRep += nrepRows;
      }

    if (nrepRows == -1)
      {
        initRep = 0;
        nrepRows = ncells;
      }

    //Looping over all the valid VoFs and setting up their stencils
    for(int k = initRep; k<initRep + nrepRows; k++)
      {
        for (it = map1[k].begin(); it != map1[k].end(); it++)
          {
            int ivof = it->second;
            EBIndex<CELL> vofInd = a_dstVoFs[ivof];
            Point pt = vofInd.m_pt;
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

            //Creating the WLS fluxes needed for flux divergence - U on a face
            BaseWLSCell<order> cellConstr;
            if (a_boundarySolutions == nullptr)
              {
                // pout() << "Using uniform div boundary values" << std::endl;
                cellConstr = BaseWLSCell<order>(radius, powerE, exactSolution, pt_verbose);
              }
            else
              {
                // pout() << "Using specified div boundary values" << std::endl;
                cellConstr = BaseWLSCell<order>(radius, powerE, a_boundarySolutions, pt_verbose);
              }

            //Computing the flux divergence stencil
            vector<Real> inhoTerm(DIM, 0.);
            WLSStencil<order> sten;
            sten.getWLSGradStencil(a_stencil[0][ivof],
                                   a_stencil[1][ivof],
                                   a_stencil[1][ivof], // dummy for the z
                                   vofInd,
                                   a_graph,
                                   a_voludata,
                                   a_ebfadata,
                                   a_xfacdata, a_yfacdata, a_zfacdata,
                                   a_ebnormxdata, a_ebnormydata,
                                   a_ebnormzdata,
                                   cellConstr,
                                   ebbc, dombc,
                                   inhoTerm,
                                   neighListType,
                                   a_dx);

            RealVect centroid;
            double kappa;
            Normalizor<order>::normalizedVoluMoments(kappa, centroid, a_voludata, vofInd, a_graph, a_dx);
            // a_stencil[ivof] *= 1.0/(a_dx);
            // for (int dir=0; dir!=SpaceDim; dir++)
            //   {
            //     a_stencil[dir][ivof] *= 1.0/(a_dx);
            //     a_inhoContr[dir][ivof] = inhoTerm[dir]/(a_dx);
            //   }
            row_stencil.push_back(a_stencil[0][ivof]); // FIXME for multi-comp
            row_inhoContr.push_back(a_inhoContr[0][ivof]);
            row_index.push_back(vofInd);

            if (pt_verbose)
              {
                pout() << "grad stencil for " << vofInd.m_pt << ":" << endl;
                pout() << " inhomg term " << a_inhoContr[0][ivof] << endl;
                printWLSStencil(radius, a_stencil[0][ivof], vofInd.m_pt);
                Real accDivF = 0.0;
                for(int ivec =  0; ivec < a_stencil[0][ivof].size(); ivec++)
                  {
                    accDivF += a_stencil[0][ivof].m_entries[ivec].m_weight;
                  }

                //pout() << "For phi=1 acc " << accDivF <<" inhoTerm "<<inhoTerm<< endl;
                pout() << "Sum of all stencil weights " << accDivF <<" stencil size "<<a_stencil[0][ivof].size()<< endl;
              }
          }
      }

//     // update stencils for periodic boundaries
//     if (nrepRows != ncells)
//       {
//         int kk =0;
//         for(int k = 0; k<map1.size(); k++)
//           {
//             for (it = map1[k].begin(); it != map1[k].end(); it++)
//               {
//                 if (k>=initRep && k<initRep + nrepRows) continue;
//                 int ivof = it->second;
//                 if (kk==row_stencil.size()) kk=0;
//                 EBIndex<CELL> vofInd = a_dstVoFs[ivof];
//                 int difi = vofInd.m_pt[0] - row_index[kk].m_pt[0];
//                 int difj = vofInd.m_pt[1] - row_index[kk].m_pt[1];
//                 LocalStencil<CELL, double> new_stencil;
//                 for(int ivec =  0; ivec < row_stencil[kk].size(); ivec++)
//                   {
//                     int pt_x = row_stencil[kk].m_entries[ivec].m_vof.m_pt[0] + difi;
//                     int pt_y = row_stencil[kk].m_entries[ivec].m_vof.m_pt[1] + difj;
//                     if (pt_y < 0) {
//                       pt_y += ncells;
//                       pt_x += ncells/nrepRows;
//                       if (difi==0) pt_x -= ncells/nrepRows;
//                     }
//                     if (pt_y > ncells-1) {
//                       pt_y -= ncells;
//                       pt_x -= ncells/nrepRows;
//                       if (difi==0) pt_x += ncells/nrepRows;
//                     }
//                     pair<int,int> pairPt (pt_x,pt_y);
//                     it2 = map2.find(pairPt);
//                     if (it2 != map2.end())
//                       {
//                         new_stencil.add(it2->second, row_stencil[kk].m_entries[ivec].m_weight);
//                       }
//                   }
//                 Point pt = vofInd.m_pt;
//                 bool pt_verbose = false;
//                 if (verbose)
//                   {
// #if DIM==3
//                     if (ptVec[0]==pt[0] && ptVec[1]==pt[1] && ptVec[2]==pt[2])
// #else
//                       if (ptVec[0]==pt[0] && ptVec[1]==pt[1])
// #endif
//                         {
//                           pt_verbose = true;
//                           pout() << "Debugging Point "<<pt<<endl;
//                         }
//                       else
//                         {
//                           pt_verbose = false;
//                         }
//                   }

//                 a_stencil[ivof] = new_stencil;
//                 a_inhoContr[ivof] = row_inhoContr[kk];
//                 kk++;

//                 if (pt_verbose)
//                   {
//                     pout() << "divF stencil for " << vofInd.m_pt << ":" << endl;
//                     printWLSStencil(radius, a_stencil[ivof][0], vofInd.m_pt);
//                     Real accDivF = 0.0;
//                     for(int ivec =  0; ivec < a_stencil[ivof][0].size(); ivec++)
//                       {
//                         accDivF += a_stencil[ivof][0].m_entries[ivec].m_weight;
//                       }

//                     pout() << "Original pt " << row_index[kk-1].m_pt[0] <<" , "<<row_index[kk-1].m_pt[1]<< endl;
//                     pout() << "Sum of all stencil weights " << accDivF <<" stencil size "<<a_stencil[ivof][0].size()<< endl;
//                   }


//               }
//           }
//       }


    return true;
  }
}
#endif
