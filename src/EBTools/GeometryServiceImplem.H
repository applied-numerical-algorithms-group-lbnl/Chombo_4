#ifndef __GeometryServiceImplem__
#define __GeometryServiceImplem__

/*********************************************/
template <int order>
LevelGeometryService<order>::
LevelGeometryService(const shared_ptr<BaseIF>     &      a_baseIF,
                     const RealVect               &      a_origin,
                     const double                 &      a_dx,
                     const Box                    &      a_domain,
                     const DisjointBoxLayout      &      a_grids,
                     const int                    &      a_ghost)
{
  PR_TIME("LevelGeometryService::LevelGeometryService");

  m_baseIF  = a_baseIF;
  m_domain  = a_domain;
  m_origin  = a_origin;
  m_dx = a_dx;
  m_geomgrow = a_ghost*IntVect::Unit;

  m_degreeP = order + 1;


  m_grids  = a_grids;
  m_domain = a_domain;

  defineData();
//finest level of moment data comes from the geometry service
//typically, you need graphs to be a bit bigger than the actual stencil

  defineFinestLevel(a_origin, a_dx);
  
}
/*********************************************/
template<int order>
void LevelGeometryService<order>::
defineData()
{
  //the 2* is because of coarsening means ghost cells need to be more numerous on finer levels than coarser levels
  m_exchangeCopier.exchangeDefine(m_grids, 2*m_geomgrow);
  m_graphs   = shared_ptr<LevelData<EBGraph > > ( new LevelData<EBGraph> (m_grids, 1, 2*m_geomgrow, NullConstructorDataFactory<EBGraph >()));
  m_voludata = shared_ptr<LevelData<VoluData> > ( new LevelData<VoluData>(m_grids, 1, 2*m_geomgrow, NullConstructorDataFactory<VoluData>()));
  m_ebfadata = shared_ptr<LevelData<EBFaData> > ( new LevelData<EBFaData>(m_grids, 1, 2*m_geomgrow, NullConstructorDataFactory<EBFaData>()));
  m_xfacdata = shared_ptr<LevelData<XFacData> > ( new LevelData<XFacData>(m_grids, 1, 2*m_geomgrow, NullConstructorDataFactory<XFacData>()));
  m_yfacdata = shared_ptr<LevelData<YFacData> > ( new LevelData<YFacData>(m_grids, 1, 2*m_geomgrow, NullConstructorDataFactory<YFacData>()));
  m_zfacdata = shared_ptr<LevelData<ZFacData> > ( new LevelData<ZFacData>(m_grids, 1, 2*m_geomgrow, NullConstructorDataFactory<ZFacData>()));
}
/*********************************************/
template<int order>
void LevelGeometryService<order>::
defineFinestLevel(const RealVect               & a_origin,
                  const double                 & a_dx )
{
  using Proto::getIrregGraphFromNodes;
  DataIterator dit = m_grids.dataIterator();
  for(int ibox = 0; ibox < dit.size(); ibox++)
  {
    Box validBox = m_grids[dit[ibox]];
    Box grownBox = validBox.grow(2*m_geomgrow);
    Box interBoxCh = grownBox & m_domain;
    Bx interBox = getProtoBox(interBoxCh);
    Bx domain   = getProtoBox(m_domain);
    if(m_baseIF->entireBoxRegular(interBox, a_origin, a_dx))
    {
      (*m_graphs)[dit[ibox]].defineAsAllRegular(interBox, domain);
    }
    else if(m_baseIF->entireBoxCovered(interBox, a_origin, a_dx))
    {
      (*m_graphs)[dit[ibox]].defineAsAllCovered(interBox, domain);
    }
    else
    {
      HostBoxData<int> regIrregCovered;
      vector< IrregNode<order> > nodes;
      fillGraph(regIrregCovered, nodes, validBox, interBoxCh, m_domain);
      vector<IrregGraph> irrgraph;
      getIrregGraphFromNodes<order>(irrgraph, nodes);

      (*m_graphs  )[dit[ibox]].defineFromPrimitives(regIrregCovered, irrgraph, interBox, domain);
//now fill in the moditment information
      (*m_voludata)[dit[ibox]].define(interBox,  (*m_graphs)[dit[ibox]]);
      (*m_ebfadata)[dit[ibox]].define(interBox,  (*m_graphs)[dit[ibox]]);
      (*m_xfacdata)[dit[ibox]].define(interBox,  (*m_graphs)[dit[ibox]]);
      (*m_yfacdata)[dit[ibox]].define(interBox,  (*m_graphs)[dit[ibox]]);
#if DIM==3
      (*m_zfacdata)[dit[ibox]].define(interBox,  (*m_graphs)[dit[ibox]]);
#endif
      for(int inode = 0; inode < nodes.size(); inode++)
      {
        const IrregNode<order>& node = nodes[inode];
        EBIndex<CELL>     vof(node.m_cell, node.m_cellIndex);
        EBIndex<BOUNDARY> ebf = vof.getCutFace();
        (*m_voludata)[dit[ibox]](vof, 0)  = node.m_volumeMoments;
        (*m_ebfadata)[dit[ibox]](ebf, 0)  = node.m_EBMoments;
        for(SideIterator sit; sit.ok(); ++sit)
        {
          vector<EBIndex<XFACE> > xfaces = (*m_graphs)[dit[ibox]].getXFaces(vof, sit());
          vector<EBIndex<YFACE> > yfaces = (*m_graphs)[dit[ibox]].getYFaces(vof, sit());
#if DIM==3
          int zfacindex = ebp_index(2, sit());
          vector<EBIndex<ZFACE> > zfaces = (*m_graphs)[dit[ibox]].getZFaces(vof, sit());
#endif
          int yfacindex = ebp_index(1, sit());
          int xfacindex = ebp_index(0, sit());

          PR_assert((xfaces.size() <= 1) && (yfaces.size() <= 1));
          if(xfaces.size() > 0)
            (*m_xfacdata)[dit[ibox]](xfaces[0], 0)  = node.m_faceMoments[xfacindex];
          if(yfaces.size() > 0)
            (*m_yfacdata)[dit[ibox]](yfaces[0], 0)  = node.m_faceMoments[yfacindex];
#if DIM==3
          if(zfaces.size() > 0)
            (*m_zfacdata)[dit[ibox]](zfaces[0], 0)  = node.m_faceMoments[zfacindex];
#endif
        }
      }
    }
  }
  exchange();
}
/*********************************************/
template <int order>
void
LevelGeometryService<order>::
exchange()
{
  Interval interv(0,0);
  m_graphs  ->exchange(m_exchangeCopier);       
  m_voludata->exchange(m_exchangeCopier);       
  m_ebfadata->exchange(m_exchangeCopier);       
  m_xfacdata->exchange(m_exchangeCopier);       
  m_yfacdata->exchange(m_exchangeCopier);       
#if DIM==3
  m_zfacdata->exchange(m_exchangeCopier);
#endif
}
/*********************************************/
/// Funny private constructor for when grids change and we have to copy data over
/// to do coarsening
/*********************************************/
template <int order>
LevelGeometryService<order>::
LevelGeometryService(const Box                               &  a_domain,
                     const DisjointBoxLayout                 &  a_grids,          
                     const shared_ptr<LevelData<EBGraph> >   &  a_graphs,
                     const shared_ptr<LevelData<VoluData> >  &  a_voludata,
                     const shared_ptr<LevelData<EBFaData> >  &  a_ebfadata,
                     const shared_ptr<LevelData<XFacData> >  &  a_xfacdata,       
                     const shared_ptr<LevelData<YFacData> >  &  a_yfacdata,       
                     const shared_ptr<LevelData<ZFacData> >  &  a_zfacdata,       
                     const RealVect                          &  a_origin,         
                     const int                               &  a_order,          
                     const int                               &  a_degreeP,
                     const double                            &  a_dx,             
                     const IntVect                           &  a_geomgrow,       
                     const shared_ptr<BaseIF>                &  a_baseIF)              
{

  m_domain    =  a_domain;         
  m_grids     =  a_grids;          
  m_origin    =  a_origin;         
  m_order     =  a_order;          
  m_degreeP   =  a_degreeP;
  m_dx        =  a_dx;             
  m_geomgrow  =  a_geomgrow;
  m_baseIF    =  a_baseIF;         

  defineData();
  
  Interval interv(0,0);
  a_graphs  ->copyTo(interv, (*m_graphs  ), interv);       
  a_voludata->copyTo(interv, (*m_voludata), interv);       
  a_ebfadata->copyTo(interv, (*m_ebfadata), interv);       
  a_xfacdata->copyTo(interv, (*m_xfacdata), interv);       
  a_yfacdata->copyTo(interv, (*m_yfacdata), interv);       
  a_zfacdata->copyTo(interv, (*m_zfacdata), interv);
  
}
  
/*********************************************/
template <int order>
LevelGeometryService<order>::
LevelGeometryService(const LevelGeometryService   &      a_finerLevel,
                     const DisjointBoxLayout      &      a_grids)
{
  m_baseIF   = a_finerLevel.m_baseIF;
  m_domain   = coarsen(a_finerLevel.m_domain,2);
  m_origin   =         a_finerLevel.m_origin;
  m_dx       =      2.*a_finerLevel.m_dx;
  m_geomgrow =         a_finerLevel.m_geomgrow;
  m_degreeP =          a_finerLevel.m_degreeP;

  m_grids  = a_grids;


  defineData();


  LevelGeometryService* finerPtr = const_cast<LevelGeometryService*>(&a_finerLevel);
  bool deleteptr = false;
  if(!a_grids.compatible(a_finerLevel.m_grids))
  {
    DisjointBoxLayout refCoarGrids;
    refine(refCoarGrids, m_grids, 2);
    finerPtr = new LevelGeometryService(a_finerLevel.m_domain,          
                                        refCoarGrids,           
                                        a_finerLevel.m_graphs,          
                                        a_finerLevel.m_voludata,        
                                        a_finerLevel.m_ebfadata,        
                                        a_finerLevel.m_xfacdata,        
                                        a_finerLevel.m_yfacdata,        
                                        a_finerLevel.m_zfacdata,        
                                        a_finerLevel.m_origin,          
                                        a_finerLevel.m_order,           
                                        a_finerLevel.m_degreeP,         
                                        a_finerLevel.m_dx,              
                                        a_finerLevel.m_geomgrow,        
                                        a_finerLevel.m_baseIF);         

      
    deleteptr = true;
  }
  
  coarsenGraph(*finerPtr);

  coarsenData( *finerPtr);
  if(deleteptr)
  {
    //if we had to copy the fine stuff to a refcoar layout, we need to copy it back
    //because the fine to coarse thing changed the fine graph
    Interval interv(0,0);
    finerPtr->m_graphs->copyTo(interv, *a_finerLevel.m_graphs, interv);
    delete finerPtr;
  }
  exchange();
}
/*********************************************/
template <int order>
void
LevelGeometryService<order>::
coarsenGraph(LevelGeometryService&      a_finerLevel)
{
  //fills finer level ghost cells 
  DataIterator dit = m_grids.dataIterator();
  for (int ibox = 0; ibox < dit.size(); ibox++)
  { 
    EBGraph & fineEBGraph = (*a_finerLevel.m_graphs)[dit[ibox]];
    Box        vofRegion       = m_grids[dit[ibox]];
    vofRegion.grow(m_geomgrow);
    vofRegion &= m_domain;
    EBGraph& coarEBGraph = (*m_graphs)[dit[ibox]];
    //this deals with the fine to coarse thing
    Bx vofRegionBx = ProtoCh::getProtoBox(vofRegion);
    coarEBGraph.coarsenVoFs( fineEBGraph, vofRegionBx);
    coarEBGraph.coarsenFaces(fineEBGraph, vofRegionBx);
  }
}
/*********************************************/
template <int order>
void
LevelGeometryService<order>::
coarsenData(LevelGeometryService&      a_finerLevel)
{
  DataIterator dit = m_grids.dataIterator();
  for (int ibox = 0; ibox < dit.size(); ibox++)
  { 
    const EBGraph& ebgraph = (*m_graphs)[dit[ibox]];
    Box          vofRegion =     m_grids[dit[ibox]];
    vofRegion.grow(m_geomgrow);
    vofRegion &= m_domain;
    Bx vofBx = ProtoCh::getProtoBox(vofRegion);
    (*m_voludata)[dit[ibox]].define(vofBx, ebgraph);
    (*m_ebfadata)[dit[ibox]].define(vofBx, ebgraph);
    (*m_xfacdata)[dit[ibox]].define(vofBx, ebgraph);
    (*m_yfacdata)[dit[ibox]].define(vofBx, ebgraph);
#if DIM==3
    (*m_zfacdata)[dit[ibox]].define(vofBx, ebgraph);
#endif

    //EBGraph& coarEBGraph = (*m_graphs)[dit[ibox]];
    //this deals with the fine to coarse thing
    Bx vofRegionBx = ProtoCh::getProtoBox(vofRegion);
    for(auto bit = vofRegionBx.begin(); bit != vofRegionBx.end(); ++bit)
    {
      const Point& pt = *bit;
      if(ebgraph.isIrregular(pt))
      {
        vector<EBIndex<CELL> > vofs = ebgraph.getVoFs(pt);
        coarsenVoFData(vofs, a_finerLevel, dit[ibox]);
        for(int ivof = 0; ivof < vofs.size(); ivof++)
        {
          const EBIndex<CELL>& vof = vofs[ivof];
          for(int idir = 0; idir < DIM; idir++)
          {
            for(SideIterator sit; sit.ok(); ++sit)
            {
              Point otherpt = pt + sign(sit())*Point::Basis(idir);
              Bx domainbx = getProtoBox(m_domain);
              Bx regionbx = getProtoBox(vofRegion);
              ///boundary faces are coarsenable.
              ///faces one in from vof boundary are coarsenable.
              if( (!domainbx.contains(otherpt)) || (regionbx.contains(otherpt)))
              {
                if(idir == 0)
                {
                  vector<EBIndex<XFACE> > xfaces = ebgraph.getXFaces(vof, sit());
                  coarsenXFaceData(xfaces, a_finerLevel, dit[ibox]);
                }
                else if(idir == 1)
                {
                  vector<EBIndex<YFACE> > yfaces = ebgraph.getYFaces(vof, sit());
                  coarsenYFaceData(yfaces, a_finerLevel, dit[ibox]);
                }
#if DIM==3
                else if(idir == 2)
                {
                  vector<EBIndex<ZFACE> > zfaces = ebgraph.getZFaces(vof, sit());
                  coarsenZFaceData(zfaces, a_finerLevel, dit[ibox]);
                }
#endif
                else
                {
                  PROTO_ASSERT(false, "bad idir in coarsening");
                }
              } //end check to see if we are in coarsening range
            }  //end loop over sides
          }   //end loop over directions
        }    //end loop over vofs in a cell
      }     //end irregular check
    }      //end loop over cells
  }       //Attention, all planets of the solar federation:
}        //We have assumed control.

/*********************************************/
template <int order>
void 
LevelGeometryService<order>::
coarsenVoFData(const vector< EBIndex<CELL> >  & a_coarVoFs,
               const LevelGeometryService     & a_finerLevel,
               const DataIndex                & a_dit)
{
  const EBGraph& coarEBGraph = (*m_graphs)[a_dit];
  const EBGraph& fineEBGraph = (*a_finerLevel.m_graphs)[a_dit];
  for(int ivof = 0; ivof < a_coarVoFs.size(); ivof++)
  {
    const EBIndex<CELL> & coarVoF = a_coarVoFs[ivof];
    IndMomDIM& coarMom  = (*m_voludata)[a_dit](coarVoF, 0);
    vector< EBIndex<CELL> > fineVoFs = coarEBGraph.refine(coarVoF);
    coarMom.setToZero();
    for(int ifine = 0; ifine < fineVoFs.size(); ifine++)
    {
      const EBIndex<CELL> & fineVoF = fineVoFs[ifine];
      IndMomDIM finemom;
      if(fineEBGraph.isCovered(fineVoF.m_pt))
      {
        finemom.setToZero(); //not sure if this can happen
      }
      else if(fineEBGraph.isRegular(fineVoF.m_pt))
      {
        finemom.setToRegular(a_finerLevel.m_dx);
      }
      else
      {
        finemom = (*(a_finerLevel.m_voludata))[a_dit](fineVoF, 0);
      }
      coarMom += finemom;  //moments are not scaled so this is a straight addition
    }
  }
}

/*********************************************/
template <int order>
void 
LevelGeometryService<order>::
coarsenXFaceData(const vector< EBIndex<XFACE> >  & a_coarFaces,
                 const LevelGeometryService      & a_finerLevel,
                 const DataIndex                 & a_dit)
{
  const EBGraph& coarEBGraph =               (*m_graphs) [a_dit];
  const EBGraph& fineEBGraph = (*(a_finerLevel.m_graphs))[a_dit];
  for(int iface = 0; iface < a_coarFaces.size(); iface++)
  {
    const EBIndex<XFACE> & coarFace = a_coarFaces[iface];
    IndMomSDMinOne & coarMom  = (*m_xfacdata)[a_dit](coarFace, 0);
    vector< EBIndex<XFACE> > fineFaces = coarEBGraph.refine<XFACE>(coarFace, fineEBGraph);
    coarMom.setToZero();
    for(int ifine = 0; ifine < fineFaces.size(); ifine++)
    {
      const EBIndex<XFACE> & fineFace = fineFaces[ifine];
      IndMomSDMinOne fineMom; 
      //by construction should be regular or irregular
      if ((*(a_finerLevel.m_xfacdata))[a_dit].contains(fineFace))
      {
        fineMom  = (*(a_finerLevel.m_xfacdata))[a_dit](fineFace, 0);
      }
      else
      {
        fineMom.setToRegular(a_finerLevel.m_dx);
      }
      coarMom += fineMom;  //moments are not scaled so this is a straight addition
    }
  }
}


/*********************************************/
template <int order>
void 
LevelGeometryService<order>::
coarsenYFaceData(const vector< EBIndex<YFACE> >  & a_coarFaces,
                 const LevelGeometryService      & a_finerLevel,
                 const DataIndex                 & a_dit)

{
  const EBGraph& coarEBGraph =               (*m_graphs) [a_dit];
  const EBGraph& fineEBGraph = (*(a_finerLevel.m_graphs))[a_dit];
  for(int iface = 0; iface < a_coarFaces.size(); iface++)
  {
    const EBIndex<YFACE> & coarFace = a_coarFaces[iface];
    IndMomSDMinOne & coarMom  = (*m_yfacdata)[a_dit](coarFace, 0);
    vector< EBIndex<YFACE> > fineFaces = coarEBGraph.refine<YFACE>(coarFace, fineEBGraph);
    coarMom.setToZero();
    for(int ifine = 0; ifine < fineFaces.size(); ifine++)
    {
      const EBIndex<YFACE> & fineFace = fineFaces[ifine];
      IndMomSDMinOne  fineMom;
      //by construction should be regular or irregular
      if((*(a_finerLevel.m_yfacdata))[a_dit].contains(fineFace))
      {
        fineMom  = (*(a_finerLevel.m_yfacdata))[a_dit](fineFace, 0);
      }
      else
      {
        fineMom.setToRegular(a_finerLevel.m_dx);
      }

      coarMom += fineMom;  //moments are not scaled so this is a straight addition
    }
  }
}
/*********************************************/
#if DIM==3
template <int order>
void 
LevelGeometryService<order>::
coarsenZFaceData(const vector< EBIndex<ZFACE> >  & a_coarFaces,
                 const LevelGeometryService      & a_finerLevel,
                 const DataIndex                 & a_dit)
{
  const EBGraph& coarEBGraph =               (*m_graphs) [a_dit];
  const EBGraph& fineEBGraph = (*(a_finerLevel.m_graphs))[a_dit];
  for(int iface = 0; iface < a_coarFaces.size(); iface++)
  {
    const EBIndex<ZFACE> & coarFace = a_coarFaces[iface];
    IndMomSDMinOne & coarMom  = (*m_zfacdata)[a_dit](coarFace, 0);
    vector< EBIndex<ZFACE> > fineFaces = coarEBGraph.refine<ZFACE>(coarFace, fineEBGraph);
    coarMom.setToZero();
    for(int ifine = 0; ifine < fineFaces.size(); ifine++)
    {
      const EBIndex<ZFACE> & fineFace = fineFaces[ifine];
      IndMomSDMinOne  fineMom;
      //by construction should be regular or irregular
      if((*(a_finerLevel.m_zfacdata))[a_dit].contains(fineFace))
      {
       fineMom  = (*(a_finerLevel.m_zfacdata))[a_dit](fineFace, 0);
      }
      else
      {
        fineMom.setToRegular(a_finerLevel.m_dx);
      }

      coarMom += fineMom;  //moments are not scaled so this is a straight addition
    }
  }
}
#endif

               
/*********************************************/
template <int order>
void 
LevelGeometryService<order>::
fillGraph(HostBoxData<int>&    a_regIrregCovered,
          vector<inode>&       a_nodes,
          const Box&           a_validRegion,
          const Box&           a_ghostRegion,
          const Box&           a_domain) const
{
  PR_TIME("LevelGeometryService::fillGraph");

  PR_assert(m_domain.contains(a_ghostRegion));
  Bx ghostRegion = ProtoCh::getProtoBox(a_ghostRegion);
  Bx validRegion = ProtoCh::getProtoBox(a_validRegion);
  a_regIrregCovered.define(ghostRegion);

  PointSet ivsirreg;
  double thrshd = 1.0e-10;
  {
    PR_TIME("boxiterator loop");
    for (auto bit = ghostRegion.begin(); bit != ghostRegion.end(); ++bit)
    {
      const Point iv =*bit;
      RvgDim cellCenter;
      for (int idir = 0;idir < DIM; ++idir)
      {
        cellCenter[idir] = m_dx*(iv[idir] +0.5) + m_origin[idir];
      }

      //member data: sign(chosen from -1,0,1) of each vertex,
      //location of each edge intersection, cellCenter,normal and gradNormal

      IndexTM<double, DIM> vectdx;
      vectdx.setAll(m_dx);

      int degreeP = order + m_degreeP;
      IFData<DIM> edgeData(*m_baseIF, vectdx, cellCenter,  degreeP);


      //create a CutCellMoment object, in order to detect whether any face coincides with the interface
      CutCellMoments <DIM> cutCell(edgeData);
      if (cutCell.isCovered())
      {
        //set covered cells to -1
        a_regIrregCovered(iv, 0) = -1;
      }
      else if (cutCell.isRegular())
      {
        //set regular cells to 1
        a_regIrregCovered(iv, 0) =  1;
      }
      else
      {
        //set irregular cells to 0
        //irregular if any face coincides with interface and edgeData.m_allVerticesIn = true
        a_regIrregCovered(iv, 0) =  0;
        if (ghostRegion.contains(iv))
        {
          ivsirreg |= iv;
        }
      }
    }
  }
//if a regular is next to a  covered, change to irregular with correct arcs and so on.
  for (auto bit = ghostRegion.begin(); bit != ghostRegion.end() ; ++bit)
  {
    const Point iv =*bit;

    if(a_regIrregCovered(iv, 0) == -1)
    {
      fixRegularCellsNextToCovered(a_nodes, a_regIrregCovered, a_ghostRegion, iv);
    }
  }

//now loop through irregular cells and make nodes for each  one.
  PointSet ivsdrop; //volumes too small to keep

  for (PointSetIterator ivsit(ivsirreg); ivsit.ok(); ++ivsit)
  {
    Point iv = ivsit();
    PR_TIME("fillGraph::endOfirregularCellLoop");
    inode newNode;

    fillNewNode(newNode,
                ivsirreg,
                ivsit());

      
    double volFrac = newNode.m_volumeMoments[IvgDim::Zero];
    if (volFrac < thrshd)
    {
      ivsdrop |= iv;
      a_regIrregCovered(iv, 0) = -1;
    }//CP record these nodes to be removed
    else
    {
      //bool fixMoments = true;
      //checkNodeMoments(newNode, m_dx, fixMoments, s_relativeTol);      
      a_nodes.push_back(newNode);
    }
  } //end loop over cells in the box
// CP: fix sweep that removes cells with volFrac less than a certain threshold
  for(PointSetIterator ivsit(ivsdrop); ivsit.ok(); ++ivsit)
  {
    const Point& iv = ivsit();
  
    for (int faceDir = 0; faceDir < DIM; faceDir++)
    {
      for (SideIterator sit; sit.ok(); ++sit)
      {
        int isign = sign(sit());
        Point otherIV = iv + isign*Point::Basis(faceDir);
        if (ghostRegion.contains(otherIV))
        {
          if (a_regIrregCovered(otherIV,0) == 0)
          {
            // i am in the case where the other cell
            // is also irregular.   I just made a previously
            // irregular cell covered so I have to check to
            // see if it had any faces pointed this way.
            int inode = -1;
            bool found = false;
            for (int ivec = 0; ivec < a_nodes.size() && ! found; ivec++)
            {
              if (a_nodes[ivec].m_cell == otherIV)
              {
                inode = ivec;
                found = true;
              }
            }
            if (!found && validRegion.contains(otherIV))
            {
              PR_error("something wrong in our logic2");
            }
            if (found)
            {
              int arcindex = ebp_index(faceDir, flip(sit()));
              a_nodes[inode].m_arc[         arcindex].resize(0);
              a_nodes[inode].m_faceMoments [arcindex].setToZero();
            }
          }
        }//valid region
      }//sit
    }//facedir
        
    //also need to fix regular cells next to new covered cell
    fixRegularCellsNextToCovered(a_nodes, a_regIrregCovered, a_ghostRegion, iv);

  }//ivsdrop
}
/*************/
template <int order>
void
LevelGeometryService<order>::
fixRegularCellsNextToCovered(vector<inode>       & a_nodes, 
                             HostBoxData<int>    & a_regIrregCovered,
                             const Box           & a_validRegion,
                             const Point         & a_iv) const

{
  Bx domain = getProtoBox(m_domain);
  Bx ivbox(a_iv, a_iv);
  Bx grownBox = ivbox.grow(1) & domain;

  PointSet ivstocheck(grownBox);
  ivstocheck -= a_iv;
  Bx  ghostRegion   = a_regIrregCovered.box();
  Bx  validRegion = getProtoBox(a_validRegion);
  //first check neighbors in each direction.  
  //If any of these are regular, they are replaced 
  //by irregular cells with a boundary face facing the covered cell.
  for(int idir = 0; idir < DIM; idir++)
  {
    for(SideIterator sit; sit.ok(); ++sit)
    {
      int ishift = sign(sit());
      Point ivshift = a_iv + ishift*Point::Basis(idir);
      ivstocheck -= ivshift;
      int bfvalshift = -1;
      if(ghostRegion.contains(ivshift))
      {
        bfvalshift = a_regIrregCovered(ivshift, 0);
      }
      if(bfvalshift  == 1)
      {
        a_regIrregCovered(ivshift, 0) =  0;

        if(validRegion.contains(ivshift))
        {
          inode newNode;
          getFullNodeWithCoveredFace(newNode, 
                                     a_regIrregCovered,
                                     ivshift);

          a_nodes.push_back(newNode);
        }

      }
    }
  }
  
  //next we loop through the remaining cells (corner cells in 2d, corner and edge cells in 3D)
  //if any of these are regular, we change them to irregular 
  for(PointSetIterator ivsit(ivstocheck); ivsit.ok(); ++ivsit)
  {
    const Point& iv = ivsit();
    if(ghostRegion.contains(iv))
    {
      if(a_regIrregCovered(iv, 0) == 1)
      {
        a_regIrregCovered(iv, 0) = 0;
        inode newNode;
        newNode.makeRegular(iv, m_dx);
        a_nodes.push_back(newNode);
      }
    }
  }
}

/**********************************************/
template <int order>
void
LevelGeometryService<order>::
getFullNodeWithCoveredFace(inode                   & a_newNode, 
                           const HostBoxData<int>  & a_regIrregCovered,
                           const Point             & a_iv) const
{

  //set all moments to regular and then zero out the appropriate face mometns
  a_newNode.makeRegular(a_iv, m_dx);
//    int coveredDir;
//    Side::LoHiSide coveredSide;
//    bool found = false;

  Bx domain = getProtoBox(m_domain);
  for (int faceDir = 0; faceDir < DIM; faceDir++)
  {
    for(SideIterator sit; sit.ok(); ++sit)
    {
      int ishift = sign(sit());
      int nodeInd = ebp_index(faceDir, sit());
      Point ivshift = a_iv + ishift*Point::Basis(faceDir);
      vector<int> arc;
      if(!domain.contains(ivshift))
      {
        // boundary arcs always -1
        arc.resize(1,-1);
        a_newNode.m_faceMoments[nodeInd].setToRegular(m_dx);
      }
      else if (a_regIrregCovered(ivshift, 0) >= 0)
      {
        //irregular cell or regular cell
        //compute vof internals returns something special if 
        //connected to a regular cell but EBGraph treats both the  same.
        //it just  knows that the cell index of a regular cell is 0
        arc.resize(1,0);
        a_newNode.m_faceMoments[nodeInd].setToRegular(m_dx);
      }
      else if (a_regIrregCovered(ivshift, 0) < 0)
      {
//          found = true;
//          coveredDir= faceDir;
//          coveredSide = sit();
        // covered face!
        arc.resize(0);
        a_newNode.m_faceMoments[nodeInd].setToZero();
      }
      else
      {
        PR_error("logic error");
      }
          
      a_newNode.m_arc[nodeInd]          = arc;
    }
  }
}

/******/
template <int order>
void 
LevelGeometryService<order>::
fillNewNode(inode               &     a_node,
            const PointSet      &     a_ivsIrreg,
            const Point         &     a_iv) const
{
  PR_TIME("LevelGeometryService::ComputeVofInternals");

  //for each CutCellMoments<dim>, we record the cell Center
  //(in physical coordinates at the global dimension)
  RvgDim cellCenter;
  for (int idir = 0;idir < DIM; ++idir)
  {
    cellCenter[idir] = m_dx*(a_iv[idir] +0.5) + m_origin[idir];
  }

  // member data: sign (chosen from -1,0,1) of each vertex,
  // location of each edge intersection, cellCenter, normal and gradNormal
  //int degreeP = m_degreeP;
  int degreeP = order + m_degreeP;
  int orderP  = 0;
  IndexTM<double, DIM> vectdx;
  vectdx.setAll(m_dx);
  IFData<DIM> edgeData(*m_baseIF,vectdx,cellCenter,  degreeP);

  //construct data holders for all moments
  MinimalCCCM<DIM> computeThisVof(edgeData);

  //compute the moments and save answers in thisVof
  computeThisVof.computeMoments(orderP,degreeP);

  CutCellMoments<DIM> thisVof = computeThisVof.m_cutCellMoments;

  vector<int> loArc[DIM];
  vector<int> hiArc[DIM];

  a_node.m_cell        = a_iv;
  a_node.m_cellIndex = 0;

  //now the arcs (0 is lo, 1 is high)
  fillArc(loArc, thisVof, 0, a_ivsIrreg, a_iv);
  fillArc(hiArc, thisVof, 1, a_ivsIrreg, a_iv);
  for(int idir = 0; idir < DIM; idir++)
  {
    int indexlo = ebp_index(idir, Side::Lo);
    int indexhi = ebp_index(idir, Side::Hi);
    a_node.m_arc[indexlo] = loArc[idir];
    a_node.m_arc[indexhi] = hiArc[idir];
  }

  a_node.m_volumeMoments = thisVof.m_moments;          
  a_node.m_EBMoments     = thisVof.m_EBmoments;

  IndexTM<double,DIM> point;
  for (int idir = 0; idir < DIM; ++idir)
  {
    double cellCent = (a_iv[idir]+ 0.5)*m_dx;
    point[idir] = cellCent;
  }  

  NormalDerivativeNew<DIM> normalDerivative;

  int maxOrder = order;

  IndexTM<double, DIM> itmpoint;
    
  convertToITM(itmpoint, point);
  IFSlicer<DIM> ifSlicer(*m_baseIF);
  NormalDerivativeNew<DIM>::NormalDerivativeMap ndMap = 
    normalDerivative.calculateAll(maxOrder,
                                  itmpoint,
                                  &ifSlicer);

  for (int idir = 0; idir < DIM; ++idir)
  {
    MomItDIM momit;
    for(momit.reset(); momit.ok(); ++momit)
    {
      double derivVal = ndMap[momit()][idir];
      a_node.m_normalPartialDeriv[idir][momit()] = derivVal;
    }
  }

  
  Iv2 bdId;
  for(SideIterator sit; sit.ok(); ++sit)
  {
    int hilo = 0;
    if(sit() == Side::Hi) hilo = 1;
      
    for (int idir = 0; idir < DIM; ++idir)
    {
      int iindex = ebp_index(idir, sit());
      PR_assert((a_node.m_arc[iindex].size() == 1) || (a_node.m_arc[iindex].size() == 0));
      if(a_node.m_arc[iindex].size() == 1)
      {
        bdId[BDID_HILO] = hilo;
        bdId[BDID_DIR] = idir;
        const CutCellMoments<DIM-1>& bdccm = thisVof.getBdCutCellMoments(bdId);
        a_node.m_faceMoments[iindex] = bdccm.m_moments;
      }
      else
      {
        a_node.m_faceMoments[iindex].setToZero();
      }
    }
  }
}

//records connectivity between vofs
template <int order>
void 
LevelGeometryService<order>::
fillArc(vector<int>                     a_arc[DIM],
        CutCellMoments<DIM>       &     a_cutCellMoments,
        const int                 &     a_hilo,
        const PointSet            &     a_ivsIrreg,
        const Point               &     a_curriv) const
{
  Iv2 bdId;
  //a_hilo is 0 or 1
  bdId[BDID_HILO] = a_hilo;
  Bx domain = getProtoBox(m_domain);
  for (int idir = 0; idir < DIM; ++idir)
  {
    bdId[BDID_DIR] = idir;
    bool covered =  a_cutCellMoments.getBdCutCellMoments(bdId).isCovered();

    if (covered)
    {
      a_arc[idir].resize(0);
    }
    else
    {
      a_arc[idir].resize(1);

      //otherIV is the iv in the idir direction on the a_hilo side
      Point otherIV = a_curriv;
      otherIV[idir] += (a_hilo*2) - 1;

      if (domain.contains(otherIV))
      {
        int otherCellIndex;
        if (a_ivsIrreg.contains(otherIV))
        {
          otherCellIndex = 0;
        }
        else
        {
          //arc to regular cell
          otherCellIndex = -2;
        }
        a_arc[idir][0]=otherCellIndex;
      }
      else if (!domain.contains(otherIV))
      {
        //boundary arcs always -1
        a_arc[idir][0] = -1;
      }
    }
  }
}

#endif
