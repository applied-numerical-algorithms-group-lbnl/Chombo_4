
#ifndef _Proto_EBStencilArchiveImplem_H_
#define _Proto_EBStencilArchiveImplem_H_
#include "Proto_Normalizor.H"
#include "Proto_SecondOrderCell.H"

namespace Proto
{
  

  ///this deep down in the bowels, everything is a template specialization
  ///returns whether there are domain boundary conditions
  template <int order>
  bool 
  EBStencilArchive<CELL, CELL,  order, double>::
  getStencil(vector<EBIndex<CELL> >                    & a_dstVoFs,
             vector<LocalStencil<CELL, double> >       & a_stencil,
             Stencil<double>                           & a_regStencil,
             vector<Stencil<double> >                    a_regStencilBCS[2*DIM],
             vector<Box>                                 a_BCApplyBoxes[2*DIM],       
             Box                                       & a_regApplyBox,
             bool                                      & a_bcOnly,
               bool                                      & a_irregOnly,
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
             double                                      a_dx,
             bool                                        a_periodiconly, 
             Point                                       a_dstGrown)               
  {
    a_bcOnly = false;
    a_irregOnly = false;
    if(a_stencilName == string("Domain_BC_Only"))
    {
      a_bcOnly = true;
      //for this one srcdomain == dstdomain
      //             srcgraph  == dstgraph
      return
        getDomainBCStencil(a_regStencilBCS,
                           a_BCApplyBoxes,
                           a_bcOnly,
                           a_domainBCName,
                           a_srcValid,
                           a_srcDomain,                  
                           a_srcGhost);                

    }
    else if(a_stencilName == string("Second_Order_Poisson"))
    {
      bool forceNeumann = false;
      //for this one srcdomain == dstdomain
      //             srcgraph  == dstgraph
      return
        getSecondOrderPoissonStencil(a_dstVoFs,                 
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
                                     a_dx,                      
                                     a_periodiconly,
                                     forceNeumann);            

    }
    else if(a_stencilName == string("Second_Order_Poisson_All_Neumann"))
    {
      bool forceNeumann = true;
      //for this one srcdomain == dstdomain
      //             srcgraph  == dstgraph
      return
        getSecondOrderPoissonStencil(a_dstVoFs,                 
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
                                     a_dx,                      
                                     a_periodiconly,
                                     forceNeumann);            

    }
    else if( (a_stencilName.find("Volume_Weighted_Averaging_rad_")) != std::string::npos)
    {

      //figure out the radius
      auto pos = a_stencilName.find_first_of("0123456789");
      unsigned long radius = 999;
      if (std::string::npos != pos)
      {
        radius = std::stoul(a_stencilName.substr(pos));
      }
      //for this one srcdomain == dstdomain
      //             srcgraph  == dstgraph
      return
        getNormalizationStencil(a_dstVoFs,                 
                                a_stencil,                 
                                a_regStencil,      
                                a_regApplyBox,             
                                a_srcValid,                
                                a_srcDomain,                  
                                a_srcGhost,                
                                a_dstGhost,                
                                a_srcGraph,               
                                a_voludata,
                                radius,
                                a_dx);                      
    }
    else if( (a_stencilName.find("Restriction")) != std::string::npos)
    {
      bool kappaWeighted = true;
      if( (a_stencilName.find("Multigrid")) != std::string::npos)
      {
        kappaWeighted = false;
      }

      return
        getRestrictionStencil(a_dstVoFs,                 
                              a_stencil,                 
                              a_regStencil,      
                              a_regApplyBox,             
                              a_srcValid,                
                              a_dstValid,                
                              a_srcGhost,                
                              a_dstGhost,                
                              a_srcGraph,                   
                              a_dstGraph,                   
                              a_voludata,
                              a_dx,
                              kappaWeighted);            

    }
    else if( (a_stencilName.find("PWC_Prolongation_")) != std::string::npos)
    {

      //figure out the color
      auto pos = a_stencilName.find_first_of("0123456789");
      unsigned long color = 999;
      if (std::string::npos != pos)
      {
       color = std::stoul(a_stencilName.substr(pos));
      }

      return
        getProlongationStencil(a_dstVoFs,                 
                               a_stencil,                 
                               a_regStencil,      
                               a_regApplyBox,             
                               a_srcValid,
                               a_dstValid,
                               a_srcGhost,                
                               a_dstGhost,                
                               a_srcGraph,                   
                               a_dstGraph,
                               color);
    }
    else if( (a_stencilName.find("Volume_Weighted_Redistribution_rad_")) != std::string::npos)
    {
      //figure out the radius
      auto pos = a_stencilName.find_first_of("0123456789");
      unsigned long radius = 999;
      if (std::string::npos != pos)
      {
        radius = std::stoul(a_stencilName.substr(pos));
      }
      //for this one srcdomain == dstdomain
      //             srcgraph  == dstgraph
      return
        getRedistributionStencil(a_dstVoFs,                 
                                 a_stencil,                 
                                 a_regStencil,      
                                 a_regApplyBox,             
                                 a_srcValid,                
                                 a_srcDomain,                  
                                 a_srcGhost,                
                                 a_dstGhost,                
                                 a_srcGraph,               
                                 a_voludata,
                                 radius,
                                 a_dx);                      

    }
    else if( (a_stencilName.find("Slope_Low_")) != std::string::npos)
    {
      //figure out the radius
      auto pos = a_stencilName.find_first_of("0123456789");
      unsigned long dir = 999;
      if (std::string::npos != pos)
      {
        dir = std::stoul(a_stencilName.substr(pos));
      }
      unsigned   int idir = dir;

      //for this one srcdomain == dstdomain
      //             srcgraph  == dstgraph
      return
        getSlopeStencil(a_dstVoFs,                 
                        a_stencil,                 
                        a_regStencil,      
                        a_regStencilBCS,
                        a_BCApplyBoxes,
                        a_regApplyBox,  
                        a_domainBCName,
                        a_srcValid,                
                        a_srcDomain,                  
                        a_srcGhost,                
                        a_dstGhost,                
                        a_srcGraph,               
                        a_voludata,
                        a_dstGrown,
                        idir,
                        Side::Lo);       //because it  is low

    }
    else if( (a_stencilName.find("Slope_High_")) != std::string::npos)
    {
      //figure out the radius
      auto pos = a_stencilName.find_first_of("0123456789");
      unsigned long dir = 999;
      if (std::string::npos != pos)
      {
        dir = std::stoul(a_stencilName.substr(pos));
      }
      unsigned   int idir = dir;
      

      //for this one srcdomain == dstdomain
      //             srcgraph  == dstgraph
      return
        getSlopeStencil(a_dstVoFs,                 
                        a_stencil,                 
                        a_regStencil,      
                        a_regStencilBCS,
                        a_BCApplyBoxes,
                        a_regApplyBox,  
                        a_domainBCName,
                        a_srcValid,                
                        a_srcDomain,                  
                        a_srcGhost,                
                        a_dstGhost,                
                        a_srcGraph,               
                        a_voludata,
                        a_dstGrown,
                        idir,
                        Side::Hi);//because it  is high

    }
    else
    {
      PR_error("stencil_archive_cell_cell::stencil string not found");
    }
    return false;//never gets here but this shuts up some compilers;
  }
  ///
  template <int order>
  bool
  EBStencilArchive<CELL, CELL,  order, double>::
  getDomainBCStencil(vector<Stencil<double> >                    a_regStencilBC[2*DIM],       
                     vector<Box>                                 a_BCApplyBoxes[2*DIM],
                     bool                                        a_bcOnly,
                     string                                      a_domainBCName[2*DIM],
                     const  Box                                & a_validBox,
                     const  Box                                & a_domain,
                     const Point                               & a_srcGhost)                   
  {
    //a_bcOnly = true;
    SecondOrderStencil<order>::
      getDomainBCs(a_regStencilBC, a_BCApplyBoxes, a_domainBCName, a_validBox, a_domain, a_srcGhost);

    //bcs are needed here
    return true;
  }
  ///
  template <int order>
  bool
  EBStencilArchive<CELL, CELL,  order, double>::
  getSlopeStencil(vector<EBIndex<CELL> >                    & a_dstVoFs,                    
                  vector<LocalStencil<CELL, double> >       & a_stencil,                    
                  Stencil<double>                           & a_regStencil,
                  vector<Stencil<double> >                    a_regStencilBC[2*DIM],       
                  vector<Box>                                 a_BCApplyBoxes[2*DIM],       
                  Box                                       & a_regApplyBox,                
                  string                                      a_domainBCName[2*DIM],
                  const Box                                 & a_validBox,                   
                  const Box                                 & a_domain,                   
                  const Point                               & a_srcGhost,                   
                  const Point                               & a_dstGhost,                   
                  const EBGraph                             & a_graph,                      
                  const VoluData                            & a_voludata,
                  const Point                               & a_dstGrown,
                  unsigned int  a_slopeDir, Side::LoHiSide a_slopeSign)
  {
    Box grownBox = a_validBox.grow(a_dstGrown) & a_domain;

    a_regApplyBox = a_validBox.grow(a_dstGrown); //need to apply outside domain 
    if(a_slopeSign == Side::Hi)
    {
      a_regStencil = (1.0)*Shift(Point::Basis(a_slopeDir)) + (-1.0)*Shift(Point::Zeroes());
    }
    else
    {
      a_regStencil = (1.0)*Shift(Point::Zeroes())          +  (-1.0)*Shift(-Point::Basis(a_slopeDir));  
    }
    
    //this stencil has a span so we have to evaluate the stencil on a grown set of the irregular cells so the regular and irregular data actually talk
    getIrregLocations(a_dstVoFs, a_regStencil, grownBox, a_graph);

    a_stencil.resize(a_dstVoFs.size());
    for(unsigned int idst = 0; idst < a_dstVoFs.size(); idst++)
    {
      a_stencil[idst].clear();
      const EBIndex<CELL>&     dstVoF = a_dstVoFs[idst];

      int isign = sign(a_slopeSign);

      vector< EBIndex<CELL> >otherVoFs = a_graph.getNeighborVoFs(dstVoF, a_slopeDir, a_slopeSign);
      if(otherVoFs.size() > 0)
      {
        double myweight    = -double(isign);
        double otherweight =  double(isign)/otherVoFs.size();
        a_stencil[idst].add(dstVoF, myweight);
        for(int jvof = 0; jvof < otherVoFs.size(); jvof++)
        {
          a_stencil[idst].add(otherVoFs[jvof], otherweight);
        }
      }
      else
      {
        vector< EBIndex<CELL> >flipVoFs = a_graph.getNeighborVoFs(dstVoF, a_slopeDir, flip(a_slopeSign));
        if(flipVoFs.size() > 0)
        {
          double myweight    =  double(isign);
          double flipweight =  -double(isign)/flipVoFs.size();
          a_stencil[idst].add(dstVoF, myweight);
          for(int jvof = 0; jvof < flipVoFs.size(); jvof++)
          {
            a_stencil[idst].add(flipVoFs[jvof], flipweight);
          }
        }
      }
    }    
    
    SecondOrderStencil<order>::
      getDomainBCs(a_regStencilBC, a_BCApplyBoxes, a_domainBCName, a_validBox, a_graph.getDomain(), a_srcGhost);

    //bcs are needed here
    return true;
  }
/******/
  template <int order>
  bool
  EBStencilArchive<CELL, CELL,  order, double>::
  getRestrictionStencil(vector<EBIndex<CELL> >                    & a_dstVoFs,                    
                        vector<LocalStencil<CELL, double> >       & a_stencil,                    
                        Stencil<double>                           & a_regStencil,         
                        Box                                       & a_regApplyBox,                
                        const Box                                 & a_srcValid,                   
                        const Box                                 & a_dstValid,                   
                        const Point                               & a_srcGhost,                   
                        const Point                               & a_dstGhost,                   
                        const EBGraph                             & a_srcGraph,                      
                        const EBGraph                             & a_dstGraph,                      
                        const VoluData                            & a_voludata,
                        double                                      a_dx,
                        bool                                        a_volumeWeighted)               
  {
#if DIM==2
    unsigned int ncolor = 4;
#else
    unsigned int ncolor = 8;
#endif    
    //this stencil has no span that I can understand so getIrrregLocations not appropriate
    a_dstVoFs = a_dstGraph.getIrregVoFs(a_dstValid);
    a_regApplyBox = a_dstValid;

    a_regStencil = Stencil<double>();
    double dnumpts = double(ncolor);
    double rweight = 1.0/dnumpts;
    for(unsigned long icolor = 0; icolor < ncolor; icolor++)
    {
      Point colorpt= getColor(icolor);
      a_regStencil += rweight*Shift(colorpt);
    }    
    a_regStencil.srcRatio() = Point::Ones(2);

    a_stencil.resize(a_dstVoFs.size());
    for(int ivof = 0; ivof < a_dstVoFs.size(); ivof++)
    {
      const EBIndex<CELL>&   coarVoF  = a_dstVoFs[ivof];
      vector<EBIndex<CELL> > fineVoFs = a_dstGraph.refine(coarVoF);
      
      if(a_volumeWeighted)
      {
        double sumvol = 0;
        for(int ifine = 0; ifine <fineVoFs.size(); ifine++)
        {
          Point pt = fineVoFs[ifine].m_pt;
          double kappa = 1;
          if(a_srcGraph.isCovered(pt))
          {
            kappa = 0;
          }
          else if(a_voludata.contains(fineVoFs[ifine]))
          {
            RealVect centroid;
            Normalizor<order>::normalizedVoluMoments(kappa, centroid, a_voludata, fineVoFs[ifine], a_srcGraph, a_dx);
          }
          sumvol += kappa;
          a_stencil[ivof].add(fineVoFs[ifine], kappa);
        }
        if(sumvol > 0)
        {
          a_stencil[ivof] *= (1.0/sumvol);
        }
      }
      else
      {
        //in multigrid, the rhs is already volume weighted so it is just
        //coarse = (1/ncolors)*(sum(fine))
        for(int ifine = 0; ifine <fineVoFs.size(); ifine++)
        {
          a_stencil[ivof].add(fineVoFs[ifine], rweight);
        }
      }

    }
    return false;//no domain boundary conditions needed here
  }

/******/
  template <int order>
  bool
  EBStencilArchive<CELL, CELL,  order, double>::
  getProlongationStencil(vector<EBIndex<CELL> >                    & a_dstVoFs,                    
                         vector<LocalStencil<CELL, double> >       & a_stencil,                    
                         Stencil<double>                           & a_regStencil,         
                         Box                                       & a_regApplyBox,                
                         const Box                                 & a_srcValid,                   
                         const Box                                 & a_dstValid,                   
                         const Point                               & a_srcGhost,                   
                         const Point                               & a_dstGhost,                   
                         const EBGraph                             & a_srcGraph,                      
                         const EBGraph                             & a_dstGraph,                      
                         unsigned long                               a_icolor)               
  {
    a_regApplyBox = a_srcValid;
    a_regStencil  = Stencil<double>();
    Point colorpt = getColor(a_icolor);
    a_regStencil  =  (0.0)*Shift(Point::Zeros());
    a_regStencil.destRatio() = Point::Ones(2);
//    a_regStencil.destShift() = colorpt;
//    //since the regular sencil has no span
//    // BUT it has  a shift and so on so I will make it grown.
    PointSet irregset(a_dstValid);
    for(PointSetIterator ivsit(irregset); ivsit.ok(); ++ivsit)
    {
      vector<EBIndex<CELL> > vofs = a_dstGraph.getVoFs(ivsit());
      a_dstVoFs.insert(a_dstVoFs.end(), vofs.begin(), vofs.end());
    }
//    a_dstVoFs = a_dstGraph.getIrregVoFs(a_dstValid);

    
    a_stencil.resize(a_dstVoFs.size());
    for(int ivof = 0; ivof < a_dstVoFs.size(); ivof++)
    {
      const EBIndex<CELL>&   fineVoF = a_dstVoFs[ivof];
      EBIndex<CELL>          coarVoF = a_dstGraph.coarsen(fineVoF);
      double weight = 1;

      Point vofpt = fineVoF.m_pt;
      for(int idir = 0; idir < DIM; idir++)
      {
        int mod = (vofpt[idir])%2;
        if(mod != colorpt[idir])
        {
          weight  = 0;
        }
      }
      
      a_stencil[ivof].add(coarVoF, weight);
    }

    return false;//no domain boundary conditions needed here
  }

  ///
  template <int order>
  bool
  EBStencilArchive<CELL, CELL,  order, double>::
  getRedistributionStencil(vector<EBIndex<CELL> >                    & a_dstVoFs,                    
                           vector<LocalStencil<CELL, double> >       & a_stencil,                    
                           Stencil<double>                           & a_regStencil,
                           Box                                       & a_regApplyBox,                
                           const Box                                 & a_validBox,                   
                           const Box                                 & a_domain,                   
                           const Point                               & a_srcGhost,                   
                           const Point                               & a_dstGhost,                   
                           const EBGraph                             & a_graph,                      
                           const VoluData                            & a_voludata,
                           int  a_redistRad,
                           double a_dx)                   
  {
    Box grownBox = a_validBox.grow(a_redistRad) & a_domain;
    //there is no reg stencil.  This should be called with  incrementonly = true (initToZero = false)
    a_regStencil = (0.0)*Shift(Point::Zeroes());
    a_regApplyBox = a_validBox;
    
    //since the regular sencil has no span, we can just get the vofs
    vector<EBIndex<CELL> > srcVoFs = a_graph.getIrregVoFs(grownBox);
    for(unsigned int ivof = 0; ivof < srcVoFs.size(); ivof++)
    {
      const EBIndex<CELL>     srcVoF = srcVoFs[ivof];
      vector<EBIndex<CELL > > neighbors = a_graph.getAllVoFsWithinRadius(srcVoF.m_pt, a_redistRad);
      
      LocalStencil<CELL, double> redistSten;
      double totvol = 0;
      for(unsigned int inei = 0; inei < neighbors.size(); inei++)
      {
        double kappa;    RealVect centroid;
        Normalizor<order>::normalizedVoluMoments(kappa, centroid, a_voludata, neighbors[inei], a_graph, a_dx);
        totvol += kappa;
        redistSten.add(neighbors[inei], kappa);
      }
      if(totvol > 0)
      {
        redistSten *= (1./totvol);
      }
      //now incrementally to a_stencil the redistSten.  Only do so if the destination is within the valid box

      for(int isten = 0; isten < redistSten.size(); isten++)
      {
        EBIndex<CELL> stenpt = redistSten.m_entries[isten].m_vof;
        double        weight = redistSten.m_entries[isten].m_weight;
        if(a_validBox.contains(stenpt.m_pt))
        {
          bool found = false;
          int index  = -1;
          for(int idst = 0; idst < a_dstVoFs.size(); idst++)
          {
            if(a_dstVoFs[idst] == stenpt)
            {
              index = idst;
              found = true;
            }
          }
          if(found)
          {
            a_stencil[index].add(stenpt, weight);
          }
          else
          {
            a_dstVoFs.push_back(stenpt);
            LocalStencil<CELL, double> dststen;
            dststen.add(stenpt, weight);
            a_stencil.push_back(dststen);
          }
        }
      }
    }    
    //no bcs to apply here
    return false;
  }
  ///
  template <int order>
  bool
  EBStencilArchive<CELL, CELL,  order, double>::
  getNormalizationStencil(vector<EBIndex<CELL> >                    & a_dstVoFs,                    
                          vector<LocalStencil<CELL, double> >       & a_stencil,                    
                          Stencil<double>                           & a_regStencil,         
                          Box                                       & a_regApplyBox,                
                          const Box                                 & a_validBox,                   
                          const Box                                 & a_domain,                     
                          const Point                               & a_srcGhost,                   
                          const Point                               & a_dstGhost,                   
                          const EBGraph                             & a_graph,                      
                          const VoluData                            & a_voludata,                   
                          unsigned long                             & a_radius,
                          double                                      a_dx)
  {
    //regular stencil is the identity matrix
    a_regStencil = (1.0)*Shift(Point::Zeroes());
    a_regApplyBox = a_validBox;
    //since the regular sencil has no span, we can just get the vofs
    a_dstVoFs = a_graph.getIrregVoFs(a_validBox);
    a_stencil.resize(a_dstVoFs.size());

    for(unsigned int ivof = 0; ivof < a_dstVoFs.size(); ivof++)
    {
      vector<EBIndex<CELL> > neighbors;
      Point pt = a_dstVoFs[ivof].m_pt;
      Box ivbox(pt, pt);
      Box neighborhood = ivbox.grow(a_radius) & a_graph.getDomain();
      for(auto bit = neighborhood.begin(); bit != neighborhood.end(); ++bit)
      {
        vector<EBIndex<CELL> > cellVoFs = a_graph.getVoFs(*bit);
        neighbors.insert(neighbors.end(), cellVoFs.begin(), cellVoFs.end());
      }

      //because stuff is kappa weighted already, all stencil weights are 1/sum(kappa)
      double totvol = 0;
      for(unsigned int inei = 0; inei < neighbors.size(); inei++)
      {
        double kappa;    RealVect centroid;
        Normalizor<order>::normalizedVoluMoments(kappa, centroid, a_voludata, neighbors[inei], a_graph, a_dx);
        totvol += kappa;
      }
      double weight = 0;
      if(totvol > 0)
      {
        weight = (1./totvol);
      }
      for(unsigned int inei = 0; inei < neighbors.size(); inei++)
      {
        a_stencil[ivof].add(neighbors[inei], weight);
      }
    }

    //no bcs to apply here
    return false;
  }              
/***/           
  template <int order>
  void
  EBStencilArchive<CELL, CELL,  order, double>::    
  getIrregLocations(vector<EBIndex<CELL> >        & a_dstVoFs,                    
                    const Stencil<double>         & a_regStencil,         
                    const Box                     & a_boundingBox,
                    const EBGraph                 & a_dstGraph)
  {
    a_dstVoFs.clear();
    PointSet irregset = a_dstGraph.getIrregCells(a_boundingBox);
    Box span = a_regStencil.span();
    Point lo = span.low();
    Point hi = span.high();
    Point grow;
    for(int idir = 0; idir < DIM; idir++)
    {
      grow[idir] = std::max(std::abs(lo[idir]), std::abs(hi[idir]));
      irregset.grow(idir, grow[idir]);
    }

    irregset &= a_dstGraph.getDomain();
    irregset &= a_boundingBox;
    for(PointSetIterator ivsit(irregset); ivsit.ok(); ++ivsit)
    {
      vector<EBIndex<CELL> > vofs = a_dstGraph.getVoFs(ivsit());
      a_dstVoFs.insert(a_dstVoFs.end(), vofs.begin(), vofs.end());
    }
  }
  ///this deep down in the bowels, everything is a template specialization
  template <int order>
  bool
  EBStencilArchive<CELL, CELL,  order, double>::
  getSecondOrderPoissonStencil(vector<EBIndex<CELL> >                    & a_dstVoFs,                    
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
                               double                                      a_dx,                         
                               bool                                        a_periodiconly,               
                               bool  a_forceNeumann)               
  {

    a_regStencil = Stencil<double>::Laplacian();

    a_regStencil *= 1.0/(a_dx*a_dx);

    Box grownBox =a_validBox.grow(1);
    grownBox &= a_domain;
    Box interiorBox = grownBox.grow(-1);
    a_regApplyBox = interiorBox;
    
    //this stencil has a span so we have to evaluate the stencil on a grown set of the irregular cells so the regular and irregular data actually talk
    getIrregLocations(a_dstVoFs, a_regStencil, a_validBox, a_graph);
//    cout << "done with getIrregLocations" << endl;

    Poisson2ndOrder<XFACE, order> xflux;
    Poisson2ndOrder<YFACE, order> yflux;
    Poisson2ndOrder<ZFACE, order> zflux;
    Poisson2ndOrder<BOUNDARY, order> ebflux;
    a_stencil.resize(a_dstVoFs.size());
    string ebbc = a_ebbcName;
    if(a_forceNeumann)
    {
      ebbc = string("Neumann");
    }
    string dombc[2*DIM];
    for(int index = 0; index < 2*DIM; index++)
    {
      if(a_forceNeumann)
        {
          dombc[index] = string("Neumann");
        }
        else
        {
          dombc[index] = a_domainBCName[index];
        }
    }
    for(int ivof = 0; ivof < a_dstVoFs.size(); ivof++)
    {
//      cout << "should not get here with no irregular cells" << endl;
      Point pt = a_dstVoFs[ivof].m_pt;
      SecondOrderStencil<order>::
        get2ndOrderDivFStencil(a_stencil[ivof],
                               a_dstVoFs[ivof],
                               a_graph,
                               a_voludata,
                               a_ebfadata,
                               a_xfacdata, a_yfacdata, a_zfacdata,
                               xflux, yflux, zflux, dombc, ebbc, 
                               a_dx);
    }

    ///add to the irregular set all regular volumes that are
    ///adjacent to the domain boundary and create their stencils

    vector<LocalStencil<CELL, double> > boundarySten;
    vector<EBIndex<     CELL        > > boundaryVolu;
    //get the volumes first
    for(auto bit = a_validBox.begin(); bit != a_validBox.end(); ++bit)
    {
      auto pt = *bit;
      bool addThisOne = false;
      //cut cells are already accounted for
      if(a_graph.isRegular(pt))
      {
        for(unsigned int idir = 0; idir < DIM; idir++)
        {
          for(SideIterator sit; sit.ok(); ++sit)
          {
            int isign = sign(sit());
            Point shift = Point::Basis(idir, isign);
            Point neigh = pt+shift;
            if(!a_domain.contains(neigh))
            {
              addThisOne = true;
            }
          }
        }
      }
      if(addThisOne)
      {
        EBIndex<CELL> vof;
        vof.defineAsRegular(pt);
        LocalStencil<CELL, double> vofstencil;
        getRegularPoissonBoundaryStencil(vofstencil, vof,
                                         a_domain, a_dx, dombc);
        boundaryVolu.push_back(vof);
        boundarySten.push_back(vofstencil);
      }
    }

    //insert the newly-treated-as-irregular stencils and volumes into the output
    a_stencil.insert(a_stencil.end(), boundarySten.begin(), boundarySten.end());
    a_dstVoFs.insert(a_dstVoFs.end(), boundaryVolu.begin(), boundaryVolu.end());
    return false;  //we are now doing domain boundary conditions internally
  }
  ////
  template <int order>
  void
  EBStencilArchive<CELL, CELL,  order, double>::
  getRegularPoissonBoundaryStencil(LocalStencil<CELL, double> & a_stencil,
                                   EBIndex<CELL>              & a_startvof,
                                   const Box                  & a_domain,
                                   double                       a_dx,
                                   string                       a_dombc[2*DIM])
  {
    a_stencil.clear();
    Point pt = a_startvof.m_pt;
    for(unsigned int idir = 0; idir < DIM; idir++)
    {
      for(SideIterator sit; sit.ok(); ++sit)
      {
        LocalStencil<CELL, double> fluxStencil;  
        int iindex = ebp_index(idir, sit());
        int isign = sign(sit());
        Point shift = Point::Basis(idir, isign);
        Point neigh = pt+shift;
        EBIndex<CELL> neighvof;
        neighvof.defineAsRegular(neigh);
        double dxsq = a_dx*a_dx;
        if(a_domain.contains(neigh))
        {
          fluxStencil.add(neighvof  ,  1./dxsq);
          fluxStencil.add(a_startvof, -1./dxsq);
        }
        else if(a_dombc[iindex] == string("Neumann"))
        {
          //no flux from a neumann boundary
        }
        else if(a_dombc[iindex] == string("Dirichlet"))
        {
          //with dirichlet,  we use the second order approximation
          //that phineigh = - phistart
          fluxStencil.add(a_startvof, -2./dxsq);
        }
        else
        {
          PROTO_ASSERT(false, "boundary condition string not  found");
        }
        a_stencil += fluxStencil;
      }
    }
    //We need this for relaxation (ignore bc as regular code did)
    double diagweight = -2.0*DIM/(a_dx*a_dx);
    a_stencil.setBetaWeight(diagweight);
  }
                                     
                                     
  /***/
/////////////end cell to cell

  ///begin CELL to fcent
  ///this deep down in the bowels, everything is a template specialization
  ///returns whether there are domain boundary conditions
  template <CENTERING fcent, int order>
  bool 
  EBStencilArchive<CELL,fcent,  order, double>::
  getStencil(vector<EBIndex<fcent> >                   & a_dstFaces,
             vector<LocalStencil<CELL, double> >       & a_stencil,
             Stencil<double>                           & a_regStencil,
             vector<Stencil<double> >                    a_regStencilBCS[2*DIM],
             vector<Box>                                 a_BCApplyBoxes[2*DIM],       
             Box                                       & a_regApplyBox,
             bool                                      & a_bcOnly,
             bool                                      & a_irregOnly,
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
             double                                      a_dx,
             bool                                        a_periodiconly, 
             Point                                       a_dstGrown)               
    {
      a_bcOnly = false;
      a_irregOnly = false;
      if(a_stencilName == string("AverageCellToFace"))
      {
        return
          getCellToFaceAveStencil(a_dstFaces,
                                  a_stencil,
                                  a_regStencil,
                                  a_regStencilBCS,
                                  a_BCApplyBoxes,       
                                  a_regApplyBox,
                                  a_domainBCName,
                                  a_srcValid,
                                  a_srcDomain,
                                  a_srcGhost,
                                  a_dstGhost,
                                  a_dstGraph,
                                  a_dx,
                                  a_dstGrown);

      }
      else if(a_stencilName == string("MAC_Gradient"))
      {
        return
          getMACGradientStencil(a_dstFaces,
                                a_stencil,
                                a_regStencil,
                                a_regStencilBCS,
                                a_BCApplyBoxes,       
                                a_regApplyBox,
                                a_srcValid,
                                a_srcDomain,
                                a_srcGhost,
                                a_dstGhost,
                                a_dstGraph,
                                a_dx,
                                a_dstGrown);

      }
      else if(a_stencilName == string("Cell_To_Face_Low"))
      {
        return
          getCellToFaceStencil(a_dstFaces,
                               a_stencil,
                               a_regStencil,
                               a_regApplyBox,
                               a_srcValid,
                               a_srcDomain,
                               a_srcGhost,
                               a_dstGhost,
                               a_dstGraph,
                               a_dx,
                               a_dstGrown,
                               Side::Lo);
                               
      }
      else if(a_stencilName == string("Cell_To_Face_High"))
      {
                               
        return
          getCellToFaceStencil(a_dstFaces,
                               a_stencil,
                               a_regStencil,
                               a_regApplyBox,
                               a_srcValid,
                               a_srcDomain,
                               a_srcGhost,
                               a_dstGhost,
                               a_dstGraph,
                               a_dx,
                               a_dstGrown,
                               Side::Hi);
      }
      else
      {
        PROTO_ASSERT(false, "string not  found");
      }

      return true;
    }

///
  template <CENTERING fcent, int order>
  bool 
  EBStencilArchive<CELL,fcent,  order, double>::
  getCellToFaceStencil(vector<EBIndex<fcent> >                   & a_dstFaces,
                       vector<LocalStencil<CELL, double> >       & a_stencil,
                       Stencil<double>                           & a_regStencil,
                       Box                                       & a_regApplyBox,
                       const Box                                 & a_validBox,
                       const Box                                 & a_srcDomain,
                       const Point                               & a_srcGhost,
                       const Point                               & a_dstGhost,
                       const EBGraph                             & a_graph,
                       double                                      a_dx,
                       Point                                       a_dstGrown,
                       Side::LoHiSide                              a_side)               
  {

    a_regStencil = (0.0)*Shift(Point::Zeros());
    a_regApplyBox = a_validBox;

    Box grown = a_validBox.grow(a_dstGrown);
    grown &= a_graph.getDomain();
    a_dstFaces = getAllFaces<fcent>(grown, a_graph);
    a_stencil.resize(a_dstFaces.size());
    for(unsigned int iface = 0; iface < a_dstFaces.size(); iface++)
    {
      EBIndex<fcent>& dstface = a_dstFaces[iface];
      EBIndex<CELL> vof = dstface.getVoF(a_side);
      a_stencil[iface].add(vof, 1.0);
    }

    return false;
  }
///
  template <CENTERING fcent, int order>
  bool 
  EBStencilArchive<CELL,fcent,  order, double>::
  getCellToFaceAveStencil(vector<EBIndex<fcent> >                   & a_dstFaces,
                          vector<LocalStencil<CELL, double> >       & a_stencil,
                          Stencil<double>                           & a_regStencil,
                          vector<Stencil<double> >                    a_regStencilBCS[2*DIM],
                          vector<Box>                                 a_BCApplyBoxes[2*DIM],       
                          Box                                       & a_regApplyBox,
                          string                                      a_domainBCName[2*DIM],
                          const Box                                 & a_validBox,
                          const Box                                 & a_srcDomain,
                          const Point                               & a_srcGhost,
                          const Point                               & a_dstGhost,
                          const EBGraph                             & a_graph,
                          double                                      a_dx,
                          Point                                       a_dstGrown)
  {
    int facedir;
    if(fcent == XFACE)
    {
      facedir = 0;
    }
    else if(fcent == YFACE)
    {
      facedir = 1;
    }
    else if(fcent == ZFACE)
    {
      facedir = 2;
    }
    else
    {
      PROTO_ASSERT(false, "bogus face direction");
    }
    a_regStencil = (1.0)*Shift(Point::Zeros());
    a_regApplyBox = a_validBox;

    Box grown = a_validBox.grow(a_dstGrown);
    a_dstFaces = getAllFaces<fcent>(grown, a_graph);
    a_stencil.resize(a_dstFaces.size());
    for(unsigned int iface = 0; iface < a_dstFaces.size(); iface++)
    {
      EBIndex<CELL> loVoF = a_dstFaces[iface].getVoF(Side::Lo);
      EBIndex<CELL> hiVoF = a_dstFaces[iface].getVoF(Side::Hi);
      a_stencil[iface].add(loVoF, 0.5);
      a_stencil[iface].add(hiVoF, 0.5);
    }

   
    SecondOrderStencil<order>::
      getDomainBCs(a_regStencilBCS, a_BCApplyBoxes, a_domainBCName, a_validBox, a_graph.getDomain(), a_srcGhost);
    return true;
  }

///
  template <CENTERING fcent, int order>
  bool 
  EBStencilArchive<CELL,fcent,  order, double>::
  getMACGradientStencil(vector<EBIndex<fcent> >                   & a_dstFaces,
                        vector<LocalStencil<CELL, double> >       & a_stencil,
                        Stencil<double>                           & a_regStencil,
                        vector<Stencil<double> >                    a_regStencilBCS[2*DIM],
                        vector<Box>                                 a_BCApplyBoxes[2*DIM],       
                        Box                                       & a_regApplyBox,
                        const Box                                 & a_validBox,
                        const Box                                 & a_srcDomain,
                        const Point                               & a_srcGhost,
                        const Point                               & a_dstGhost,
                        const EBGraph                             & a_graph,
                        double                                      a_dx,
                        Point                                       a_dstGrown)
  {
    a_regStencil = (1.0)*Shift(Point::Zeros());
    a_regApplyBox = a_validBox;

    Box grown = a_validBox.grow(a_dstGrown);
    a_dstFaces = getAllFaces<fcent>(grown, a_graph);
    a_stencil.resize(a_dstFaces.size());
    for(unsigned int iface = 0; iface < a_dstFaces.size(); iface++)
    {
      EBIndex<CELL> loVoF = a_dstFaces[iface].getVoF(Side::Lo);
      EBIndex<CELL> hiVoF = a_dstFaces[iface].getVoF(Side::Hi);
      a_stencil[iface].add(loVoF, -1.0/a_dx);
      a_stencil[iface].add(hiVoF,  1.0/a_dx);
    }


    using std::cout;
    using std::endl;
    string domainBCName[2*DIM];
    for(int ivec = 0; ivec < 2*DIM; ivec++)
    {
      domainBCName[ivec]= string("LinearExtrapolation");
    }
    
    SecondOrderStencil<order>::
      getDomainBCs(a_regStencilBCS, a_BCApplyBoxes, domainBCName, a_validBox, a_graph.getDomain(), a_srcGhost);
    return true;
  }

////end cell to face
  ///begin  fcent to cell

  ///this deep down in the bowels, everything is a template specialization
  ///returns whether there are domain boundary conditions
  template <CENTERING fcent, int order>
  bool 
  EBStencilArchive<fcent,  CELL, order, double>::
  getStencil(vector<EBIndex<CELL> >                   & a_dstVoFs,
             vector<LocalStencil<fcent, double> >       & a_stencil,
             Stencil<double>                           & a_regStencil,
             vector<Stencil<double> >                    a_regStencilBCS[2*DIM],
             vector<Box>                                 a_BCApplyBoxes[2*DIM],       
             Box                                       & a_regApplyBox,
             bool                                      & a_bcOnly,
             bool                                      & a_irregOnly,
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
             double                                      a_dx,
             bool                                        a_periodiconly, 
             Point                                       a_dstGrown)               
    {

      a_bcOnly = false;
      a_irregOnly = false;
      if(a_stencilName == string("Divergence"))
      {
        return
          getDivergenceStencil(a_dstVoFs,
                               a_stencil,
                               a_regStencil,
                               a_regApplyBox,
                               a_srcValid,
                               a_srcDomain,
                               a_srcGhost,
                               a_dstGhost,
                               a_voludata,
                               a_ebfadata,
                               a_xfacdata,
                               a_yfacdata,
                               a_zfacdata,
                               a_dstGraph,
                               a_dx);
                               
      }
      else if(a_stencilName == string("AverageFaceToCell"))
      {
        return
          getFaceToCellAveStencil(a_dstVoFs,
                                  a_stencil,
                                  a_regStencil,
                                  a_regApplyBox,
                                  a_srcValid,
                                  a_srcDomain,
                                  a_srcGhost,
                                  a_dstGhost,
                                  a_voludata,
                                  a_ebfadata,
                                  a_xfacdata,
                                  a_yfacdata,
                                  a_zfacdata,
                                  a_dstGraph,
                                  a_dx);
                               
      }
      else
      {
        PROTO_ASSERT(false, "string not  found");
      }

      return true;
    }

  template <CENTERING fcent, int order>
  void
  EBStencilArchive<fcent, CELL,  order, double>::    
  getIrregLocations(vector<EBIndex<CELL> >        & a_dstVoFs,                    
                    const Stencil<double>         & a_regStencil,         
                    const Box                     & a_boundingBox,
                    const EBGraph                 & a_dstGraph,
                    const int                     & a_facedir)
  {
    a_dstVoFs.clear();
    PointSet irregset = a_dstGraph.getIrregCells(a_boundingBox);
    Box span = a_regStencil.span();
    Point lo = span.low();
    Point hi = span.high();
    Point grow = Point::Zero();
    int idir = a_facedir;
    grow[idir] = std::max(std::abs(lo[idir]), std::abs(hi[idir]));
    irregset.grow(idir, grow[idir]);


    irregset &= a_dstGraph.getDomain();
    irregset &= a_boundingBox;
    for(PointSetIterator ivsit(irregset); ivsit.ok(); ++ivsit)
    {
      vector<EBIndex<CELL> > vofs = a_dstGraph.getVoFs(ivsit());
      a_dstVoFs.insert(a_dstVoFs.end(), vofs.begin(), vofs.end());
    }
  }
///
  template <CENTERING fcent, int order>
  bool 
  EBStencilArchive<fcent, CELL, order, double>::
  getDivergenceStencil(vector<EBIndex<CELL > >                   & a_dstVoFs,
                       vector<LocalStencil<fcent, double> >      & a_stencil,
                       Stencil<double>                           & a_regStencil,
                       Box                                       & a_regApplyBox,
                       const Box                                 & a_validBox,
                       const Box                                 & a_srcDomain,
                       const Point                               & a_srcGhost,
                       const Point                               & a_dstGhost,
                       const VoluData                            & a_voludata,
                       const EBFaData                            & a_ebfadata,
                       const XFacData                            & a_xfacdata,
                       const YFacData                            & a_yfacdata,
                       const ZFacData                            & a_zfacdata,
                       const EBGraph                             & a_graph,
                       double                                      a_dx)
  {
    int facedir;
    if(fcent == XFACE)
    {
      facedir = 0;
    }
    else if(fcent == YFACE)
    {
      facedir = 1;
    }
    else if(fcent == ZFACE)
    {
      facedir = 2;
    }
    else
    {
      PROTO_ASSERT(false, "bogus face direction");
    }
    
//    a_regStencil = (1.0/a_dx)*Shift(Point::Basis(facedir, 1)) +  (-1.0/a_dx)*Shift(Point::Zeros());
    //this stencil has a span so we have to evaluate the stencil on a grown set of the irregular cells so the regular and irregular data actually talk
//    getIrregLocations(a_dstVoFs, a_regStencil, a_validBox, a_graph, facedir);
    a_regApplyBox = a_validBox;
  
    a_regStencil = (1.0)*Shift(Point::Zeros());

    a_dstVoFs = a_graph.getAllVoFs(a_validBox);
    a_stencil.resize(a_dstVoFs.size());
    for(unsigned int ivof = 0; ivof < a_dstVoFs.size(); ivof++)
    {
      EBIndex<CELL> vof = a_dstVoFs[ivof];
      RealVect centroid;
      double kappa;
      RealVect normal;
      double bndryArea;
      Normalizor<order>::getNormalAndBoundaryArea(bndryArea, normal,   
                                                a_voludata,
                                                a_ebfadata,
                                                a_xfacdata,
                                                a_yfacdata,
                                                a_zfacdata, 
                                                vof, a_graph, a_dx);
      Normalizor<order>::normalizedVoluMoments(kappa, centroid, a_voludata, vof, a_graph, a_dx);
      
      for(SideIterator sit; sit.ok(); ++sit)
      {
        vector<EBIndex<fcent> > faces = getSideFaces<fcent>(vof, sit(), a_graph);
        for(int iface = 0; iface < faces.size(); iface++)
        {
          const EBIndex<fcent>& face = faces[iface];
          double areaFrac;
          RealVect centroid;
          Normalizor<order>::normalizedFaceMoments(areaFrac, centroid, a_xfacdata, a_yfacdata, a_zfacdata, face, a_graph, a_dx, facedir);
          double weight = sign(sit())*areaFrac/a_dx;

          a_stencil[ivof].add(face, weight);
        }
      }
    }
    
    //no bcs here.
    return false;
  }
  template <CENTERING fcent, int order>
  bool 
  EBStencilArchive<fcent, CELL, order, double>::
  getFaceToCellAveStencil(vector<EBIndex<CELL > >                   & a_dstVoFs,
                          vector<LocalStencil<fcent, double> >      & a_stencil,
                          Stencil<double>                           & a_regStencil,
                          Box                                       & a_regApplyBox,
                          const Box                                 & a_validBox,
                          const Box                                 & a_srcDomain,
                          const Point                               & a_srcGhost,
                          const Point                               & a_dstGhost,
                          const VoluData                            & a_voludata,
                          const EBFaData                            & a_ebfadata,
                          const XFacData                            & a_xfacdata,
                          const YFacData                            & a_yfacdata,
                          const ZFacData                            & a_zfacdata,
                          const EBGraph                             & a_graph,
                          double                                      a_dx)
  {
    int facedir;
    if(fcent == XFACE)
    {
      facedir = 0;
    }
    else if(fcent == YFACE)
    {
      facedir = 1;
    }
    else if(fcent == ZFACE)
    {
      facedir = 2;
    }
    else
    {
      PROTO_ASSERT(false, "bogus face direction");
    }

    //reg stencil did not work here.
    a_regStencil = (1.0)*Shift(Point::Zeros());
    a_regApplyBox = a_validBox;

    a_dstVoFs = a_graph.getAllVoFs(a_validBox);
    a_stencil.resize(a_dstVoFs.size());

    for(unsigned int ivof = 0; ivof < a_dstVoFs.size(); ivof++)
    {
      EBIndex<CELL> vof = a_dstVoFs[ivof];
      vector<EBIndex<fcent> > faceslo = getSideFaces<fcent>(vof, Side::Lo, a_graph);
      vector<EBIndex<fcent> > faceshi = getSideFaces<fcent>(vof, Side::Hi, a_graph);
      vector<EBIndex<fcent> > allfaces;
      allfaces.insert(allfaces.end(), faceslo.begin(), faceslo.end());
      allfaces.insert(allfaces.end(), faceshi.begin(), faceshi.end());
      int nfacestot = allfaces.size();
      double weight = 0;
      if(nfacestot > 0)
      {
        weight = 1.0/double(nfacestot);
      }
      for(int iface = 0; iface < allfaces.size(); iface++)
      {
        a_stencil[ivof].add(allfaces[iface], weight);
      }
    }
    
    //no bcs here.
    return false;
  }

///end face to cell



  ///begin  fcent to fcent
  ///this deep down in the bowels, everything is a template specialization
  ///returns whether there are domain boundary conditions
  template <CENTERING fcent, int order>
  bool 
  EBStencilArchive<fcent,  fcent, order, double>::
  getStencil(vector<EBIndex<fcent> >                   & a_dstVoFs,
             vector<LocalStencil<fcent, double> >       & a_stencil,
             Stencil<double>                           & a_regStencil,
             vector<Stencil<double> >                    a_regStencilBCS[2*DIM],
             vector<Box>                                 a_BCApplyBoxes[2*DIM],       
             Box                                       & a_regApplyBox,
             bool                                      & a_bcOnly,
             bool                                      & a_irregOnly,
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
             double                                      a_dx,
             bool                                        a_periodiconly, 
             Point                                       a_dstGrown)
    {

      a_bcOnly = false;
      if(a_stencilName == string("InterpolateToFaceCentroid"))
      {
        return
          getFaceCentroidInterpStencil(a_dstVoFs,
                                       a_stencil,
                                       a_regStencil,
                                       a_regApplyBox,
                                       a_irregOnly,
                                       a_srcValid,
                                       a_srcDomain,
                                       a_srcGhost,
                                       a_dstGhost,
                                       a_xfacdata,
                                       a_yfacdata,
                                       a_zfacdata,
                                       a_dstGraph,
                                       a_dx);
      }
      else if(a_stencilName == string("extrap_to_domain_face"))
      {
        return
          getDomainExtrapolationStencil(a_dstVoFs,
                                       a_stencil,
                                       a_regStencil,
                                       a_regApplyBox,
                                       a_irregOnly,
                                       a_srcValid,
                                       a_srcDomain,
                                       a_srcGhost,
                                       a_dstGhost,
                                       a_xfacdata,
                                       a_yfacdata,
                                       a_zfacdata,
                                       a_dstGraph,
                                       a_dx);
      }
      else
      {
        PROTO_ASSERT(false, "string not  found");
      }

      return true;
    }

  template <CENTERING fcent, int order>
  inline LocalStencil<fcent, double>
  getFaceCentroidInterpStencil(const EBIndex<fcent>                           & a_face, 
                               const HostIrregData<XFACE, IndexedMoments<DIM-1, order>, 1>  & a_xfacdata,
                               const HostIrregData<YFACE, IndexedMoments<DIM-1, order>, 1>  & a_yfacdata,
                               const HostIrregData<ZFACE, IndexedMoments<DIM-1, order>, 1>  & a_zfacdata,
                               const EBGraph                                  & a_graph,
                               const double                                   & a_dx)
  {
    PROTO_ASSERT(false, "not implemented");
  }


  template <int order>
  inline LocalStencil<XFACE, double>
  getFaceCentroidInterpStencil(const EBIndex<XFACE>                           & a_face, 
                               const HostIrregData<XFACE, IndexedMoments<DIM-1, order>, 1>  & a_xfacdata,
                               const HostIrregData<YFACE, IndexedMoments<DIM-1, order>, 1>  & a_yfacdata,
                               const HostIrregData<ZFACE, IndexedMoments<DIM-1, order>, 1>  & a_zfacdata,
                               const EBGraph                                  & a_graph,
                               const double                                   & a_dx)
  {
    LocalStencil<XFACE, double> retval;
    Base2ndOrderFlux<XFACE, order> xflux;
    xflux.getFaceInterpStencil(retval, a_face, a_graph, a_xfacdata, a_dx);

    return retval;
  }
  

  template <int order>
  inline LocalStencil<YFACE, double>
  getFaceCentroidInterpStencil(const EBIndex<YFACE>                           & a_face, 
                               const HostIrregData<XFACE, IndexedMoments<DIM-1, order>, 1>  & a_xfacdata,
                               const HostIrregData<YFACE, IndexedMoments<DIM-1, order>, 1>  & a_yfacdata,
                               const HostIrregData<ZFACE, IndexedMoments<DIM-1, order>, 1>  & a_zfacdata,
                               const EBGraph                                  & a_graph,
                               const double                                   & a_dx)
  {
    LocalStencil<YFACE, double> retval;
    Base2ndOrderFlux<YFACE, order> yflux;
    yflux.getFaceInterpStencil(retval, a_face, a_graph, a_yfacdata, a_dx);
    return retval;
  }
  

  template <int order>
  inline LocalStencil<ZFACE, double>
  getFaceCentroidInterpStencil(const EBIndex<ZFACE>                           & a_face, 
                               const HostIrregData<XFACE, IndexedMoments<DIM-1, order>, 1>  & a_xfacdata,
                               const HostIrregData<YFACE, IndexedMoments<DIM-1, order>, 1>  & a_yfacdata,
                               const HostIrregData<ZFACE, IndexedMoments<DIM-1, order>, 1>  & a_zfacdata,
                               const EBGraph                                  & a_graph,
                               const double                                   & a_dx)
  {
    LocalStencil<ZFACE, double> retval;
    Base2ndOrderFlux<ZFACE, order> zflux;
    zflux.getFaceInterpStencil(retval, a_face, a_graph, a_zfacdata, a_dx);
    return retval;
  }
  
  template <CENTERING fcent, int order>
  bool 
  EBStencilArchive<fcent, fcent, order, double>::
  getFaceCentroidInterpStencil(vector<EBIndex<fcent > >                  & a_dstFaces,
                               vector<LocalStencil<fcent, double> >      & a_stencil,
                               Stencil<double>                           & a_regStencil,
                               Box                                       & a_regApplyBox,
                               bool                                      & a_irregOnly,
                               const Box                                 & a_validBox,
                               const Box                                 & a_srcDomain,
                               const Point                               & a_srcGhost,
                               const Point                               & a_dstGhost,
                               const XFacData                            & a_xfacdata,
                               const YFacData                            & a_yfacdata,
                               const ZFacData                            & a_zfacdata,
                               const EBGraph                             & a_graph,
                               const double                              & a_dx)
  {
    //well defined in regular faces
    a_irregOnly = false;
    int facedir;
    if(fcent == XFACE)
    {
      facedir = 0;
    }
    else if(fcent == YFACE)
    {
      facedir = 1;
    }
    else if(fcent == ZFACE)
    {
      facedir = 2;
    }
    else
    {
      PROTO_ASSERT(false, "bogus face direction");
    }
    //no regular stencil span so we don't need to do anything fancy with dstFaces
    a_regStencil = (1.0)*Shift(Point::Zeros());
    a_regApplyBox = a_validBox.growHi(facedir, 1);

    Box grown = a_validBox.grow(1);
    a_dstFaces = getIrregFaces<fcent>(grown, a_graph);

    a_stencil.resize(a_dstFaces.size());
    for(unsigned int iface = 0; iface < a_dstFaces.size(); iface++)
    {
      const EBIndex<fcent>& face = a_dstFaces[iface];
      a_stencil[iface] =
        Proto::getFaceCentroidInterpStencil(face, a_xfacdata, a_yfacdata, a_zfacdata, a_graph, a_dx);
    }
    //no bcs here.
    return false;
  }
  template <CENTERING fcent, int order>
  bool 
  EBStencilArchive<fcent, fcent, order, double>::
  getDomainExtrapolationStencil(vector<EBIndex<fcent > >                  & a_dstFaces,
                                vector<LocalStencil<fcent, double> >      & a_stencil,
                                Stencil<double>                           & a_regStencil,
                                Box                                       & a_regApplyBox,
                                bool                                      & a_irregOnly,
                                const Box                                 & a_validBox,
                                const Box                                 & a_srcDomain,
                                const Point                               & a_srcGhost,
                                const Point                               & a_dstGhost,
                                const XFacData                            & a_xfacdata,
                                const YFacData                            & a_yfacdata,
                                const ZFacData                            & a_zfacdata,
                                const EBGraph                             & a_graph,
                                const double                              & a_dx)
  {
    int facedir;
    if(fcent == XFACE)
    {
      facedir = 0;
    }
    else if(fcent == YFACE)
    {
      facedir = 1;
    }
    else if(fcent == ZFACE)
    {
      facedir = 2;
    }
    else
    {
      PROTO_ASSERT(false, "bogus face direction");
    }
    //
    a_irregOnly = true;
    //no regular stencil span so we don't need to do anything fancy with dstFaces
    //just copies where it is not extrapolating
    a_regStencil = (1.0)*Shift(Point::Zeros());
    a_regApplyBox = a_validBox.growHi(facedir, 1);
    for(SideIterator sit; sit.ok(); ++sit)
    {
      Side::LoHiSide flipsd = flip(sit());
      int ishift   = sign(flipsd);
      Point pshift = Point::Basis(facedir, ishift);
      
      Box adjcell = a_validBox.adjCellSide(facedir, 1, sit());
      Box bndry   = adjcell.shift(pshift);
      for(auto bit = bndry.begin(); bit != bndry.end(); ++bit)
      {
        vector<EBIndex<CELL> > vofs = a_graph.getVoFs(*bit);
        for(unsigned int ivof = 0; ivof < vofs.size(); ivof++)
        {
          auto vof = vofs[ivof];
          vector<EBIndex<fcent> > bfaces = getFaces<fcent>(vof, sit(), a_graph);
          for(unsigned int iface = 0; iface < bfaces.size(); iface++)
          {
            if(bfaces[iface].m_isBoundary)
            {
              auto dstface = bfaces[iface];
              LocalStencil<fcent, double> locStencil;
              EBIndex<fcent> nearface, farface, farfarface;
              bool hasnear, hasfar, hasfarfar;
              vector<EBIndex<fcent> > nearfaces = getFaces<fcent>(vof, flipsd, a_graph);
              hasnear = (nearfaces.size() == 1) && (!nearfaces[0].isBoundary());
              if(hasnear)
              {
                nearface = nearfaces[0];
                EBIndex<CELL> nextvof = nearface.getVoF(flipsd);
                vector<EBIndex<fcent> > farfaces = getFaces<fcent>(nextvof, flipsd, a_graph);
                hasfar = ((farfaces.size() == 1) && (!farfaces[0].isBoundary()));
                if(hasfar)
                {
                  farface = farfaces[0];
                  EBIndex<CELL> nextnextvof = farface.getVoF(flipsd);
                  vector<EBIndex<fcent> > farfarfaces = getFaces<fcent>(nextnextvof, flipsd, a_graph);
                  hasfarfar = ((farfarfaces.size() == 1) && (!farfarfaces[0].isBoundary()));
                  if(hasfarfar)
                  {
                    farfarface = farfarfaces[0];
                  }
                }
              }

              if(hasfarfar)
              {
                locStencil.add(  nearface, 3.0);
                locStencil.add(   farface,-3.0);
                locStencil.add(farfarface, 1.0);
              }
              else if(hasfar)
              {
                locStencil.add(nearface, 2.0);
                locStencil.add( farface,-1.0);
              }
              else if(hasnear)
              {
                locStencil.add(nearface, 1.0);
              }
              if(hasnear || hasfar)
              {
                a_dstFaces.push_back(dstface);
                a_stencil.push_back(locStencil);
              }
            }
          } //end loop over faces
        }//end loop over vofs at edge
      }//end loop over cells at edge
    }//end loop over sides of the input box

    //no bcs here in the sense that we do not touch the source
    return false;
  }
///end face to cell
//now for the cross face stuff
//
  template <CENTERING srcCenter, CENTERING dstCenter, int order>
  class EBCrossStencilArchive
  {
  public:
    typedef IndexedMoments<DIM  , order> IndMomDIM;
    typedef IndexedMoments<DIM-1, order> IndMomSDMinOne;
    typedef MomentIterator<DIM  , order> MomItDIM;
    typedef MomentIterator<DIM-1, order> MomItSDMinOne;
    typedef HostIrregData<CELL,      IndMomDIM , 1>  VoluData;
    typedef HostIrregData<BOUNDARY,  IndMomDIM , 1>  EBFaData;
    typedef HostIrregData<XFACE, IndMomSDMinOne, 1>  XFacData;
    typedef HostIrregData<YFACE, IndMomSDMinOne, 1>  YFacData;
    typedef HostIrregData<ZFACE, IndMomSDMinOne, 1>  ZFacData;

    inline static bool
    getStencil(vector<EBIndex<dstCenter> >               & a_dstFaces,
               vector<LocalStencil<srcCenter,  double> > & a_stencil,
               Stencil<double>                           & a_regStencilInterior,
               vector<Stencil<double> >                    a_regStencilBCS[2*DIM],
               vector<Box>                                 a_BCApplyBoxes[2*DIM],       
               Box                                       & a_regApplyBox,
               bool                                      & a_bcOnly,
               bool                                      & a_irregOnly,
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
               double                                      a_dx,
               bool                                        a_periodiconly,
               Point                                       a_dstGrown,
               unsigned int a_srcDir, unsigned int a_dstDir)    
    {
      a_bcOnly = false;
      a_irregOnly = false;
      PR_assert( a_srcGraph.getDomain() == a_dstGraph.getDomain());
      Box domain = a_srcGraph.getDomain();
      
      if(a_stencilName == string("Tangential_Velocity_Correction"))
      {
        a_dstFaces = getAllFaces<dstCenter>(a_dstValid, a_dstGraph);
        a_stencil.resize(a_dstFaces.size());
        for(unsigned int idst = 0; idst < a_dstFaces.size(); idst++)
        {
          const EBIndex<dstCenter>& dstf = a_dstFaces[idst];

          vector<EBIndex<srcCenter> > allSrcFaces;
          for(SideIterator sito; sito.ok(); ++sito)
          {
            EBIndex<CELL> vof = dstf.getVoF(sito());
            if(domain.contains(vof.m_pt))
            {
              for(SideIterator siti; siti.ok(); ++siti)
              {   
                vector<EBIndex<srcCenter> > srcFaces= getSideFaces<srcCenter>(vof, siti(), a_dstGraph);
                allSrcFaces.insert(allSrcFaces.end(), srcFaces.begin(), srcFaces.end());
              }
            }
          }
          double weight = 0;
          if(allSrcFaces.size() > 0)
          {
            weight = 1.0/double(allSrcFaces.size());
          }
          for(unsigned int isten = 0; isten < allSrcFaces.size(); isten++)
          {
            a_stencil[idst].add(allSrcFaces[isten], weight);
          }
        }
      }
      else
      {
        PROTO_ASSERT(false, "could not find cross stencil name");
      }
      return false;

    }           
  };

  ///
  template <int order>
  bool 
  EBStencilArchive<YFACE, XFACE, order, double>::
  getStencil(vector<EBIndex<XFACE> >                   & a_dstFaces,
             vector<LocalStencil<YFACE, double> >      & a_stencil,
             Stencil<double>                           & a_regStencil,
             vector<Stencil<double> >                    a_regStencilBCS[2*DIM],
             vector<Box>                                 a_BCApplyBoxes[2*DIM],       
             Box                                       & a_regApplyBox,
             bool                                      & a_bcOnly,
               bool                                      & a_irregOnly,
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
             double                                      a_dx,
             bool                                        a_periodiconly, 
             Point                                       a_dstGrown)              
  {
    a_bcOnly = false;
    a_irregOnly = false;
    return EBCrossStencilArchive<YFACE,XFACE, order>::
      getStencil(a_dstFaces,
                 a_stencil,               
                 a_regStencil,            
                 a_regStencilBCS,
                 a_BCApplyBoxes,
                 a_regApplyBox,
                 a_bcOnly,
                 a_irregOnly,
                 a_stencilName,           
                 a_domainBCName,          
                 a_ebbcName,              
                 a_srcValid,              
                 a_dstValid,              
                 a_srcDomain,             
                 a_dstDomain,             
                 a_srcGhost,              
                 a_dstGhost,              
                 a_srcGraph,              
                 a_dstGraph,              
                 a_voludata,              
                 a_ebfadata,              
                 a_xfacdata,              
                 a_yfacdata,              
                 a_zfacdata,              
                 a_dx,                    
                 a_periodiconly,          
                 a_dstGrown,
                 1, 0);
  }



  template <int order>
  bool 
  EBStencilArchive<XFACE, YFACE, order, double>::
  getStencil(vector<EBIndex<YFACE> >                   & a_dstFaces,
             vector<LocalStencil<XFACE, double> >      & a_stencil,
             Stencil<double>                           & a_regStencil,
             vector<Stencil<double> >                    a_regStencilBCS[2*DIM],
             vector<Box>                                 a_BCApplyBoxes[2*DIM],       
             Box                                       & a_regApplyBox,
             bool                                      & a_bcOnly,
               bool                                      & a_irregOnly,
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
             double                                      a_dx,
             bool                                        a_periodiconly, 
             Point                                       a_dstGrown)              
  {
    a_bcOnly = false;
    a_irregOnly = false;
    return EBCrossStencilArchive<XFACE,YFACE, order>::
      getStencil(a_dstFaces,
                 a_stencil,               
                 a_regStencil,            
                 a_regStencilBCS,
                 a_BCApplyBoxes,
                 a_regApplyBox,           
                 a_bcOnly,
                 a_irregOnly,
                 a_stencilName,           
                 a_domainBCName,          
                 a_ebbcName,              
                 a_srcValid,              
                 a_dstValid,              
                 a_srcDomain,             
                 a_dstDomain,             
                 a_srcGhost,              
                 a_dstGhost,              
                 a_srcGraph,              
                 a_dstGraph,              
                 a_voludata,              
                 a_ebfadata,              
                 a_xfacdata,              
                 a_yfacdata,              
                 a_zfacdata,              
                 a_dx,                    
                 a_periodiconly,          
                 a_dstGrown,
                 0, 1);
  }


  ///
  template <int order>
  bool 
  EBStencilArchive<ZFACE, XFACE, order, double>::
  getStencil(vector<EBIndex<XFACE> >                   & a_dstFaces,
             vector<LocalStencil<ZFACE, double> >      & a_stencil,
             Stencil<double>                           & a_regStencil,
             vector<Stencil<double> >                    a_regStencilBCS[2*DIM],
             vector<Box>                                 a_BCApplyBoxes[2*DIM],       
             Box                                       & a_regApplyBox,
             bool                                      & a_bcOnly,
             bool                                      & a_irregOnly,
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
             double                                      a_dx,
             bool                                        a_periodiconly, 
             Point                                       a_dstGrown)              
  {
    a_bcOnly = false;
    a_irregOnly = false;
    return EBCrossStencilArchive<ZFACE,XFACE, order>::
      getStencil(a_dstFaces,
                 a_stencil,               
                 a_regStencil,            
                 a_regStencilBCS,
                 a_BCApplyBoxes,
                 a_regApplyBox,           
                 a_bcOnly,
                 a_irregOnly,
                 a_stencilName,           
                 a_domainBCName,          
                 a_ebbcName,              
                 a_srcValid,              
                 a_dstValid,              
                 a_srcDomain,             
                 a_dstDomain,             
                 a_srcGhost,              
                 a_dstGhost,              
                 a_srcGraph,              
                 a_dstGraph,              
                 a_voludata,              
                 a_ebfadata,              
                 a_xfacdata,              
                 a_yfacdata,              
                 a_zfacdata,              
                 a_dx,                    
                 a_periodiconly,          
                 a_dstGrown,
                 2, 0);
  }



  template <int order>
  bool 
  EBStencilArchive<XFACE, ZFACE, order, double>::
  getStencil(vector<EBIndex<ZFACE> >                   & a_dstFaces,
             vector<LocalStencil<XFACE, double> >      & a_stencil,
             Stencil<double>                           & a_regStencil,
             vector<Stencil<double> >                    a_regStencilBCS[2*DIM],
             vector<Box>                                 a_BCApplyBoxes[2*DIM],       
             Box                                       & a_regApplyBox,
             bool                                      & a_bcOnly,
               bool                                      & a_irregOnly,
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
             double                                      a_dx,
             bool                                        a_periodiconly, 
             Point                                       a_dstGrown)              
  {
    a_bcOnly = false;
    a_irregOnly = false;
    return EBCrossStencilArchive<XFACE,ZFACE, order>::
      getStencil(a_dstFaces,               
                 a_stencil,               
                 a_regStencil,            
                 a_regStencilBCS,
                 a_BCApplyBoxes,
                 a_regApplyBox,           
                 a_bcOnly,
                 a_irregOnly,
                 a_stencilName,           
                 a_domainBCName,          
                 a_ebbcName,              
                 a_srcValid,              
                 a_dstValid,              
                 a_srcDomain,             
                 a_dstDomain,             
                 a_srcGhost,              
                 a_dstGhost,              
                 a_srcGraph,              
                 a_dstGraph,              
                 a_voludata,              
                 a_ebfadata,              
                 a_xfacdata,              
                 a_yfacdata,              
                 a_zfacdata,              
                 a_dx,                    
                 a_periodiconly,          
                 a_dstGrown,
                 0, 2);
  }



//////////

  template <int order>
  bool 
  EBStencilArchive<ZFACE, YFACE, order, double>::
  getStencil(vector<EBIndex<YFACE> >                   & a_dstFaces,
             vector<LocalStencil<ZFACE, double> >      & a_stencil,
             Stencil<double>                           & a_regStencil,
             vector<Stencil<double> >                    a_regStencilBCS[2*DIM],
             vector<Box>                                 a_BCApplyBoxes[2*DIM],       
             Box                                       & a_regApplyBox,
             bool                                      & a_bcOnly,
             bool                                      & a_irregOnly,
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
             double                                      a_dx,
             bool                                        a_periodiconly, 
             Point                                       a_dstGrown)              
  {
    a_bcOnly = false;
    a_irregOnly = false;
    return EBCrossStencilArchive<ZFACE,YFACE, order>::
      getStencil(a_dstFaces,               
                 a_stencil,               
                 a_regStencil,            
                 a_regStencilBCS,
                 a_BCApplyBoxes,
                 a_regApplyBox,           
                 a_bcOnly,
                 a_irregOnly,
                 a_stencilName,           
                 a_domainBCName,          
                 a_ebbcName,              
                 a_srcValid,              
                 a_dstValid,              
                 a_srcDomain,             
                 a_dstDomain,             
                 a_srcGhost,              
                 a_dstGhost,              
                 a_srcGraph,              
                 a_dstGraph,              
                 a_voludata,              
                 a_ebfadata,              
                 a_xfacdata,              
                 a_yfacdata,              
                 a_zfacdata,              
                 a_dx,                    
                 a_periodiconly,          
                 a_dstGrown,
                 2, 1);
  }



  template <int order>
  bool 
  EBStencilArchive<YFACE, ZFACE, order, double>::
  getStencil(vector<EBIndex<ZFACE> >                   & a_dstFaces,
             vector<LocalStencil<YFACE, double> >      & a_stencil,
             Stencil<double>                           & a_regStencil,
             vector<Stencil<double> >                    a_regStencilBCS[2*DIM],
             vector<Box>                                 a_BCApplyBoxes[2*DIM],       
             Box                                       & a_regApplyBox,
             bool                                      & a_bcOnly,
             bool                                      & a_irregOnly,
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
             double                                      a_dx,
             bool                                        a_periodiconly, 
             Point                                       a_dstGrown)              
  {
    a_bcOnly = false;
    a_irregOnly = false;
    return EBCrossStencilArchive<YFACE,ZFACE, order>::
      getStencil(a_dstFaces, 
                 a_stencil,               
                 a_regStencil,            
                 a_regStencilBCS,
                 a_BCApplyBoxes,
                 a_regApplyBox,           
                 a_bcOnly,
                 a_irregOnly,
                 a_stencilName,           
                 a_domainBCName,          
                 a_ebbcName,              
                 a_srcValid,              
                 a_dstValid,              
                 a_srcDomain,             
                 a_dstDomain,             
                 a_srcGhost,              
                 a_dstGhost,              
                 a_srcGraph,              
                 a_dstGraph,              
                 a_voludata,              
                 a_ebfadata,              
                 a_xfacdata,              
                 a_yfacdata,              
                 a_zfacdata,              
                 a_dx,                    
                 a_periodiconly,          
                 a_dstGrown,
                 1, 2);
  }

#if 1
  //cell to boundary
  template <int order>
  bool
  EBStencilArchive<CELL, BOUNDARY, order, double>::
  getStencil(vector<EBIndex<BOUNDARY> >                & a_dstVoFs,
             vector<LocalStencil<CELL , double> >      & a_stencil,
             Stencil<double>                           & a_regStencil,
             vector<Stencil<double> >                    a_regStencilBCS[2*DIM],
             vector<Box>                                 a_BCApplyBoxes[2*DIM],       
             Box                                       & a_regApplyBox,
             bool                                      & a_bcOnly,
             bool                                      & a_irregOnly,
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
             double                                      a_dx,
             bool                                        a_periodiconly, 
             Point                                       a_dstGrown)
  {
    a_irregOnly = true;
    a_bcOnly = false;

    PR_assert( a_srcGraph.getDomain() == a_dstGraph.getDomain());
    Box domain = a_srcGraph.getDomain();

    if(a_stencilName == string("copy_cell_value_to_cut_face"))
    {
      a_dstVoFs = a_srcGraph.getCutFaces(a_srcValid);
      a_stencil.resize(a_dstVoFs.size());
      for(unsigned int idst = 0; idst < a_dstVoFs.size(); idst++)
      {
        auto& dstv = a_dstVoFs[idst];

        EBIndex<CELL> srcVoF = getEBIndex<CELL, BOUNDARY>(dstv);
          
        double weight = 1.;
        a_stencil[idst].clear();
        a_stencil[idst].add(srcVoF, weight);
      }
    }
    else
    {
      PROTO_ASSERT(false, "could not find  stencil name cell to eb");
    }
    return false;
  }

  //boundary to cell
  template <int order>
  bool
  EBStencilArchive<BOUNDARY, CELL, order, double>::
  getStencil(vector<EBIndex<CELL> >                    & a_dstVoFs,
             vector<LocalStencil<BOUNDARY, double> >   & a_stencil,
             Stencil<double>                           & a_regStencil,
             vector<Stencil<double> >                    a_regStencilBCS[2*DIM],
             vector<Box>                                 a_BCApplyBoxes[2*DIM],       
             Box                                       & a_regApplyBox,
             bool                                      & a_bcOnly,
             bool                                      & a_irregOnly,
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
             double                                      a_dx,
             bool                                        a_periodiconly, 
             Point                                       a_dstGrown)
  {
    a_irregOnly = true;
    a_bcOnly = false;

    PR_assert( a_srcGraph.getDomain() == a_dstGraph.getDomain());
    Box domain = a_srcGraph.getDomain();
      
    if(a_stencilName.find("cut_face_increment_to_kappa_div_dir_") != std::string::npos)
    {

      //figure out the radius
      auto pos = a_stencilName.find_first_of("0123456789");
      unsigned long longdir = 999;
      if (std::string::npos != pos)
      {
        longdir = std::stoul(a_stencilName.substr(pos));
      }
      unsigned int graddir = longdir;

      a_dstVoFs = a_srcGraph.getIrregVoFs(a_srcValid);
      a_stencil.resize(a_dstVoFs.size());
      for(unsigned int idst = 0; idst < a_dstVoFs.size(); idst++)
      {
        auto& dstv = a_dstVoFs[idst];

        EBIndex<BOUNDARY> srcCutFace=getEBIndex<BOUNDARY, CELL>(dstv);
          
        RealVect normal;
        double bndryArea;
        Normalizor<order>::getNormalAndBoundaryArea(bndryArea, normal,   
                                                    a_voludata,
                                                    a_ebfadata,
                                                    a_xfacdata,
                                                    a_yfacdata,
                                                    a_zfacdata, 
                                                    dstv, a_srcGraph, a_dx);
        double normcomp = normal[graddir];
        //negative because of our very weird normal standard (inward facing)
        double weight = -bndryArea*normcomp/a_dx; 
        a_stencil[idst].clear();
        a_stencil[idst].add(srcCutFace, weight);
      }
    }
    else
    {
      PROTO_ASSERT(false, "could not find  stencil name eb to cell");
    }
    return false;

  }
#endif  
}
#endif
