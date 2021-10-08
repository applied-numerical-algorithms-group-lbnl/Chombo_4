#ifndef __Proto_EBGraphIMPLEM__
#define __Proto_EBGraphIMPLEM__

namespace Proto
{

  /*****/
  void 
  EBGraph::
  copy(const EBGraph              & a_src,
       const Box                  & a_srcBox,
       unsigned int                 a_srcComp,
       const Box                  & a_destBox,
       unsigned int                 a_destComp,
       unsigned int                 a_numcomp)
  {
    PR_assert(a_srcBox == a_destBox);  //periodic and EB don't mix
    Box region = a_srcBox & m_region & a_src.m_region;
    if(isRegular(region) && a_src.isRegular(region))
    {
      return;
    }
    else if(isCovered(region) && a_src.isCovered(region))
    {
      return;
    }
    else if(a_src.isCovered(region) && region.contains(m_region))
    {
      defineAsAllCovered(m_validBox, m_region, m_domain);
      return;
    }
    else if(a_src.isRegular(region) && region.contains(m_region))
    {
      defineAsAllRegular(m_validBox, m_region, m_domain);
      return;
    }
    else if (isAllRegular() && a_src.isAllCovered())
    {
      defineAsFullIrregular(m_validBox, m_region, m_domain);
      for(BoxIterator bit = region.begin(); bit != region.end(); ++bit)
      {
        m_internals->m_graph(*bit, 0).defineAsCovered();
      }
      
    }
    else if (isAllCovered() && a_src.isAllRegular())
    {
      defineAsEmptyIrregular(m_validBox, m_region, m_domain);
      for(BoxIterator bit = region.begin(); bit != region.end(); ++bit)
      {
        m_internals->m_graph(*bit, 0).defineAsRegular();
      }
      
    }
  else
    {
      //one or both has irregular cells.
      //because i am sick of combinatorics,
      //use basefab copy to transfer the data.

      //see if we need to generate a source fab
      HostBoxData<CellNode, 1>* srcFabPtr = NULL;
      bool needToDelete = false;
      if(a_src.hasIrregular())
      {
        srcFabPtr = (HostBoxData<CellNode, 1>*)(&(a_src.m_internals->m_graph));
      }
      else
      {
        srcFabPtr = new HostBoxData<CellNode, 1>(region);
        needToDelete = true;

        CellNode val;
        if(a_src.isAllRegular())
        {
          val.defineAsRegular();
        }
        else if(a_src.isAllCovered())
        {
          val.defineAsCovered();
        }
        else
        {
          PROTO_ASSERT(false, "logic error");
        }
        srcFabPtr->setVal(val);
      }
      //define our graph if I need to, leave alone otherwise
      if(isAllRegular())
      {
        defineAsFullIrregular(m_validBox, m_region, m_domain);
      }
      else if(isAllCovered())
      {
        defineAsEmptyIrregular(m_validBox, m_region, m_domain);
      }
      PointSet ivsInter;
      if(a_src.hasIrregular())
      {      
        ivsInter = a_src.m_internals->m_irregIVS;
      }
      ivsInter &= region;
      m_internals->m_irregIVS |= ivsInter;
      m_internals->m_graph.copy(*srcFabPtr, region, 0, region, 0, 1);
      if(needToDelete)
      {
        delete srcFabPtr;
      }
    }
  }
  /***
  **/
  std::size_t 
  EBGraph::
  charsize(const Box&   a_bx,
           unsigned int a_startcomp,
           unsigned int a_numcomps) const
  {
    using ProtoSer::linearCharBufSize;
    size_t retval = 0;
    //m_tag
    retval += sizeof(TAG);
    // m_region
    retval += linearCharBufSize(m_region);
    // m_validBox
    retval += linearCharBufSize(m_validBox);
    // m_domain
    retval += linearCharBufSize(m_domain);
    //if we have a graph, include its components;
    if(m_tag == HasIrregular)
    {
      const auto& graph = m_internals->m_graph;
      const auto& irreg = m_internals->m_irregIVS;
      //first the irregular cells
      Box gbox = graph.box();
      retval += linearCharBufSize(irreg);
      retval += linearCharBufSize(gbox);
      for(auto bit = gbox.begin();  bit != gbox.end(); ++bit)
      {
        auto pt = *bit;
        retval += linearCharBufSize(graph(pt, 0));
      }
    }
//    //begin debug
//    using std::cout;
//    using std::endl;
//    cout << "ebgraph::charsize valid box = " << m_validBox  << ", charsiz = " << retval << endl;
//    //begin debug
    return retval;
  }
  /*
    create a version of me as a  stream of bytes.
  */
  void
  EBGraph::
  energize(void*  a_buf) const
  {
    using ProtoSer::linearOutToCharBuf;
    using ProtoSer::linearCharBufSize;

    size_t charsiz = 0;
    char* charbuf = (char *) a_buf;
    //m_tag
    TAG* tagbuf = (TAG *) charbuf;
    *tagbuf = m_tag;
    charbuf += sizeof(TAG);
    charsiz += sizeof(TAG);
    // m_region
    linearOutToCharBuf(charbuf, m_region);
    charbuf += linearCharBufSize(m_region);
    charsiz += linearCharBufSize(m_region);
    // m_validBox
    linearOutToCharBuf(charbuf, m_validBox);
    charbuf += linearCharBufSize(m_validBox);
    charsiz += linearCharBufSize(m_validBox);
    // m_domain
    linearOutToCharBuf(charbuf ,m_domain);
    charbuf += linearCharBufSize(m_domain);
    charsiz += linearCharBufSize(m_domain);
    //if we have a graph, include its components;
    if(m_tag == HasIrregular)
    {
      const auto& graph = m_internals->m_graph;
      const auto& irreg = m_internals->m_irregIVS;
      //first the irregular cells
      linearOutToCharBuf(charbuf, irreg);
      charbuf += linearCharBufSize(irreg);
      charsiz += linearCharBufSize(irreg);
      //now the graph box
      Box gbox = graph.box();
      linearOutToCharBuf(charbuf, gbox);
      charbuf += gbox.linearSize();
      charsiz += gbox.linearSize();
      //now all the graph stuff
      for(auto bit = gbox.begin();  bit != gbox.end(); ++bit)
      {
        const auto& pt = *bit;
        const auto& node = graph(pt, 0);
        linearOutToCharBuf(charbuf, node);
        charsiz += linearCharBufSize(node);
        charbuf += linearCharBufSize(node);
      }
    }
    //begin debug
//    using std::cout;
//    using std::endl;
//    cout << "ebgraph::energize: valid box = " << m_validBox  << ", charsiz = " << charsiz << endl;
    //end   debug
  }
  /***/
  EBGraph::
  EBGraph(void*  a_buf)
  {
    m_isDefined = true;
    using ProtoSer::linearInFromCharBuf;
    using ProtoSer::linearCharBufSize;

    size_t charsiz = 0;
    char* charbuf = (char *) a_buf;
    //m_tag
    TAG* tagbuf = (TAG *) charbuf;
    m_tag   = *tagbuf;
    charbuf += sizeof(TAG);
    charsiz += sizeof(TAG);
    // m_region
    linearInFromCharBuf(m_region, charbuf);
    charbuf += linearCharBufSize(m_region);
    charsiz += linearCharBufSize(m_region);
    // m_validBox
    linearInFromCharBuf(m_validBox, charbuf);
    charbuf += linearCharBufSize(m_validBox);
    charsiz += linearCharBufSize(m_validBox);
    // m_domain
    linearInFromCharBuf(m_domain, charbuf);
    charbuf += linearCharBufSize(m_domain);
    charsiz += linearCharBufSize(m_domain);
    //if we have a graph, include its components;
    m_internals = shared_ptr<localData>(new localData());
    if(m_tag == HasIrregular)
    {
      auto& graph = m_internals->m_graph;
      auto& irreg = m_internals->m_irregIVS;
      //first the irregular cells
      linearInFromCharBuf(irreg, charbuf);
      charbuf += linearCharBufSize(irreg);
      charsiz += linearCharBufSize(irreg);
      //now the graph box
      Box gbox;
      linearInFromCharBuf(gbox, charbuf);
      charbuf += linearCharBufSize(gbox);
      charsiz += linearCharBufSize(gbox);
      //define the data holder
      graph.define(gbox);
      //now all the graph stuff
      for(auto bit = gbox.begin();  bit != gbox.end(); ++bit)
      {
        auto pt = *bit;
        auto& node = graph(pt, 0);
        linearInFromCharBuf(node, charbuf);
        charsiz += linearCharBufSize(node);
        charbuf += linearCharBufSize(node);
      }
    }
//    //begin debug
//    using std::cout;
//    using std::endl;
//    cout << "ebgraph::ebgraph(buf) valid box = " << m_validBox  << ", charsiz = " << charsiz << endl;
//    //end   debug
  }
  /***
     The darkness drops again; but now I know   
     That twenty centuries of stony sleep
     Were vexed to nightmare by a rocking cradle,   
     And what rough beast, its hour come round at last,   
     Slouches towards Bethlehem to be born?
     --WB Yeats.
  **/
  void  
  EBGraph::
  linearIn(void*                  a_buf,
           const Box       &      a_bx,
           unsigned int           a_startcomp,
           unsigned int           a_numcomps)
  {
    EBGraph srcGraph(a_buf);
    copy(srcGraph, a_bx, 0, a_bx, 0, 1);
  }

  /**
  ***/    
  void 
  EBGraph::
  linearOut(void*               a_buf,
            const Box       &   a_bx,
            unsigned int        a_startcomp,
            unsigned int        a_numcomps) const
  {
    energize(a_buf);
  }


  /*****/
  template <CENTERING cent>
  vector< EBIndex<cent> > 
  EBGraph::
  refine(const EBIndex<cent> & a_ind,
         const EBGraph       & a_finer) const
  {
    vector< EBIndex<cent> >  retval;
    
    
    EBIndex<CELL> voflo = a_ind.getVoF(Side::Lo);
    EBIndex<CELL> vofhi = a_ind.getVoF(Side::Hi);
    if(!a_ind.m_isBoundary)
    {
      vector<EBIndex<CELL> > vofslo = refine(voflo);
      vector<EBIndex<CELL> > vofshi = refine(vofhi);
      for(int ilo = 0; ilo < vofslo.size(); ilo++)
      {
        for(int ihi = 0; ihi < vofshi.size(); ihi++)
        {
          if(a_finer.isConnected(vofslo[ilo], vofshi[ihi]))
          {
            EBIndex<cent> face;
            face.m_pt      = vofshi[ihi].m_pt;
            face.m_vofIDMe = vofshi[ihi].m_vofIDMe;
            face.m_vofIDLo = vofslo[ilo].m_vofIDMe;
            face.m_isBoundary = false;
            retval.push_back(face);
          }
        }
      }
    }
    else
    {
      //from here we know it it is a boundary face
      if(m_domain.contains(voflo.m_pt))
      {
        vector<EBIndex<CELL> > vofslo = refine(voflo);
        int idir = 0;
        if(cent == YFACE)
        {
          idir = 1;
        }
        else if(cent == ZFACE)
        {
          idir = 2;
        }
        for(int ilo = 0; ilo < vofslo.size(); ilo++)
        {
          EBIndex<cent> face;
          face.m_pt      = vofslo[ilo].m_pt + Point::Basis(idir);
          face.m_vofIDMe = -1;
          face.m_vofIDLo = vofslo[ilo].m_vofIDMe;
          face.m_isBoundary = true;
          retval.push_back(face);
        }
      }
      else
      {
        vector<EBIndex<CELL> > vofshi = refine(vofhi);
        for(int ihi = 0; ihi < vofshi.size(); ihi++)
        {
          EBIndex<cent> face;
          face.m_pt      = vofshi[ihi].m_pt;
          face.m_vofIDLo = -1;
          face.m_vofIDMe = vofshi[ihi].m_vofIDMe;
          face.m_isBoundary = true;
          retval.push_back(face);
        }
      }
    }
    return retval;
  }


  /*****/
  template <  >
  inline vector< EBIndex<BOUNDARY> > 
  EBGraph::
  refine(const EBIndex<BOUNDARY> & a_ind,
         const EBGraph           & a_finer) const
  {
    

    vector<EBIndex<BOUNDARY> > retval;
    EBIndex<CELL> vof;
    vof.m_pt = a_ind.m_pt;
    vof.m_vofIDMe = a_ind.m_vofIDMe;
    vof.m_vofIDLo = 0;
    vof.m_isBoundary = false;
    vector<EBIndex<CELL> > vofs = refine(vof);
    for(int ivof = 0; ivof < vofs.size(); ivof++)
    {
      const Point& pt  = vofs[ivof].m_pt;
      const int& vofid = vofs[ivof].m_vofIDMe;
      if(a_finer.isIrregular(vof.m_pt))
      {
        EBIndex<BOUNDARY> cutcell;
        cutcell.m_pt      = pt;
        cutcell.m_vofIDMe = vofid;
        cutcell.m_vofIDLo = 0;
        cutcell.m_isBoundary = false;
        retval.push_back(cutcell);
      }
    }
    return retval;
  }
  /*****/
  template < >
  inline vector< EBIndex<CELL> > 
  EBGraph::
  refine(const EBIndex<CELL> & a_ind,
         const EBGraph       & a_finer) const
  {
    
      return  refine(a_ind);
  }



  /****/
  vector< EBIndex<CELL> > 
  EBGraph::
  refine(const EBIndex<CELL> & a_vof) const
  {
    vector<EBIndex<CELL> > retval;
    if(!isCovered(a_vof.m_pt)) //covered returns an empty vector
    {
      if(isRegular(a_vof.m_pt))
      {
        Box finedom = m_domain.refine(2); //not really used for the cell case but required in the API
        retval = EBIndex<CELL>::getRefinedRegular(a_vof, finedom);
      }
      else
      {
        const CellNode               & cellnode =   m_internals->m_graph(a_vof.m_pt, 0);
        const CellNode::IrrGraphNode & irrnode  = (*cellnode.m_cellList)[a_vof.m_vofIDMe];
        retval = irrnode.m_finer;
      }
    }
    return retval;
  }

  /*****/
  template <CENTERING cent>
  vector< EBIndex<cent> > 
  EBIndex<cent>::getRefinedRegular(const EBIndex<cent>  & a_coar,
                                   const Box            & a_fineDomain) 
  {
    vector< EBIndex<cent> > retval;
    if(cent == BOUNDARY)
    {
      //boundary cells are not regular
    }
    else if(cent == CELL)
    {
      Box coarBox(a_coar.m_pt, a_coar.m_pt);
      Box fineBox = coarBox.refine(2);
      for(auto bit = fineBox.begin(); bit != fineBox.end(); ++bit)
      {
        EBIndex<cent> finevof;
        finevof.m_pt = *bit;
        finevof.m_vofIDMe = 0;
        finevof.m_vofIDLo = 0;
        finevof.m_isBoundary = false;
        retval.push_back(finevof);
      }
    }
    else
    {
      Box coarBox(a_coar.m_pt, a_coar.m_pt);
      Box fineBox = coarBox.refine(2);
      //Box faceBox = fineBox.grow(idir, -1);
      for(auto bit = fineBox.begin(); bit != fineBox.end(); ++bit)
      {
        Point pt = *bit;
        EBIndex<cent> fineind;
        fineind.m_pt = pt;
        fineind.m_isBoundary = a_coar.m_isBoundary;
        fineind.m_vofIDMe    =  a_coar.m_vofIDMe;
        fineind.m_vofIDLo    =  a_coar.m_vofIDLo;

        retval.push_back(fineind);
      }
    }
    return retval;
  }

  ///
  template <CENTERING cent>
  EBIndex<cent>  
  EBGraph::
  coarsen(const EBIndex<cent> & a_ind) const
  {
    PR_assert(!isCovered(a_ind.m_pt));
    EBIndex<cent> retval;
    if(cent == CELL)
    {
      retval.m_pt = a_ind.m_pt.coarsen(2);
      retval.m_isBoundary = false;
      retval.m_vofIDLo = 0;
      if(!m_domain.contains(a_ind.m_pt))
      {
        retval.m_vofIDMe = -1;
      }
      else if(isRegular(a_ind.m_pt))
      {
        retval.m_vofIDMe = 0;
      }
      else
      {
        const CellNode               & cellnode = m_internals->m_graph(a_ind.m_pt, 0);
        const CellNode::IrrGraphNode & irrnode  = (*cellnode.m_cellList)[a_ind.m_vofIDMe];
        retval.m_vofIDMe = irrnode.m_coarser;
      }
    }
    else
    {
      EBIndex<CELL> voflo = coarsen(a_ind.getVoF(Side::Lo));
      EBIndex<CELL> vofhi = coarsen(a_ind.getVoF(Side::Hi));

      retval.m_pt      = vofhi.m_pt;
      retval.m_vofIDMe = vofhi.m_vofIDMe;
      retval.m_vofIDLo = voflo.m_vofIDMe;
    }
    return retval;
  }

  ///
  template <CENTERING cent> 
  bool
  EBGraph::
  getAdjacentFace(EBIndex<cent>       & a_adjFace, 
                  const EBIndex<cent> & a_startFace,
                  int a_tanDir, Side::LoHiSide a_side) const
  {
    PR_assert(a_tanDir != a_startFace.direction());
    bool hasFace = true;
    Point hiiv = a_startFace.m_pt + sign(a_side)*Point::Basis(a_tanDir, 1);
    Point loiv = hiiv - Point::Basis(a_startFace.direction(), 1);
    if((!m_domain.contains(loiv))  || (!m_domain.contains(hiiv)))
    {
      hasFace = false;
      return hasFace;
    }
    EBIndex<CELL> loVoF = a_startFace.getVoF(Side::Lo);
    EBIndex<CELL> hiVoF = a_startFace.getVoF(Side::Hi);

    if(a_tanDir == 0)
    {
      vector<EBIndex<XFACE> > faceslo = getXFaces(loVoF,  a_side);
      vector<EBIndex<XFACE> > faceshi = getXFaces(hiVoF,  a_side);
      EBIndex<CELL> voflo, vofhi;
      bool found = false;
      for(int ilo = 0; ilo < faceslo.size(); ilo++)
      {
        voflo = faceslo[ilo].getVoF(a_side);
        for(int ihi = 0; ihi < faceshi.size(); ihi++)
        {
          vofhi = faceshi[ihi].getVoF(a_side);
          if(isConnected(voflo, vofhi))
          {
            found = true;
            break;
          }
        }
      }
      if(!found)
      {
        hasFace = false;
        return hasFace;
      }
      else
      {
        a_adjFace = EBIndex<cent>(voflo, vofhi);
      }
    }
    else if(a_tanDir == 1)
    {
      vector<EBIndex<YFACE> > faceslo = getYFaces(loVoF,    a_side);
      vector<EBIndex<YFACE> > faceshi = getYFaces(hiVoF,    a_side);
      EBIndex<CELL> voflo, vofhi;
      bool found = false;
      for(int ilo = 0; ilo < faceslo.size(); ilo++)
      {
        voflo = faceslo[ilo].getVoF(a_side);
        for(int ihi = 0; ihi < faceshi.size(); ihi++)
        {
          vofhi = faceshi[ihi].getVoF(a_side);
          if(isConnected(voflo, vofhi))
          {
            found = true;
            break;
          }
        }
      }
      if(!found)
      {
        hasFace = false;
        return hasFace;
      }
      else
      {
        a_adjFace = EBIndex<cent>(voflo, vofhi);
      }
    }
#if DIM==3
    else if(a_tanDir == 2)
    {
      vector<EBIndex<ZFACE> > faceslo = getZFaces(loVoF,  a_side);
      vector<EBIndex<ZFACE> > faceshi = getZFaces(hiVoF,  a_side);
      EBIndex<CELL> voflo, vofhi;
      bool found = false;
      for(int ilo = 0; ilo < faceslo.size(); ilo++)
      {
        voflo = faceslo[ilo].getVoF(a_side);
        for(int ihi = 0; ihi < faceshi.size(); ihi++)
        {
          vofhi = faceshi[ihi].getVoF(a_side);
          if(isConnected(voflo, vofhi))
          {
            found = true;
            break;
          }
        }
      }
      if(!found)
      {
        hasFace = false;
        return hasFace;
      }
      else
      {
        a_adjFace = EBIndex<cent>(voflo, vofhi);
      }
    }
#endif
    else
    {
      PR_error("bogus tandir");
    }
    return hasFace;
  }

  
  template <CENTERING cent>
  void
  CellNode::
  setFaceStuff(EBIndex<cent>             &   a_face,
               const int                 &   a_inputIndex,
               const int                 &   a_otherIndex,
               const Point               &   a_inputIV,
               const Point               &   a_otherIV,
               const Side::LoHiSide      &   a_sd,
               const Box                 &   a_domain) const
  {
    a_face.m_isBoundary = false;
    int otherIndex = a_otherIndex;
    int inputIndex = a_inputIndex;
    if(!a_domain.contains(a_otherIV)) 
    {
      otherIndex = -1;
      a_face.m_isBoundary = true;
    }
    if(a_sd == Side::Hi)
    {
      a_face.m_pt      = a_otherIV;
      a_face.m_vofIDMe = otherIndex;
      a_face.m_vofIDLo = inputIndex;
    }
    else 
    {
      a_face.m_pt      = a_inputIV;
      a_face.m_vofIDMe = inputIndex;
      a_face.m_vofIDLo = otherIndex;
    }
  }

    
  template <CENTERING cent>
  void
  CellNode::
  getFaceInfo(vector< EBIndex<cent> >  & a_faces,
              const int                & a_dir,
              const  EBIndex<CELL>     & a_vof,
              const Side::LoHiSide     & a_sd,
              const Box                & a_domain) const
  {
    PR_assert(a_domain.contains(a_vof.m_pt)) ;
    if(isCovered())
    {
      a_faces.resize(0);
    }
    else
    {

      Point inputIV    = a_vof.m_pt;
      Point otherIV    = a_vof.m_pt +sign(a_sd)*Point::Basis(a_dir, 1);
      int inputIndex   = a_vof.m_vofIDMe;
      if(isRegular())
      {
        a_faces.resize(1);

        //in regular land this ought to be zero
        //if it is not, then there is a mismatch between graphs and the cellnode 
        PR_assert(inputIndex == 0);
        int otherIndex = 0;
        setFaceStuff(a_faces[0], inputIndex, otherIndex, inputIV, otherIV, a_sd, a_domain);
      }
      else  //irregular case
      {
        vector<IrrGraphNode>& nodeVec = *m_cellList;
        IrrGraphNode vofnode = nodeVec[inputIndex];
        int arcindex = ebp_index(a_dir, a_sd);
        vector<int> arcvec = vofnode.m_arc[arcindex];
        a_faces.resize(arcvec.size());
        for(int iface = 0; iface < a_faces.size(); iface++)
        {
          int otherIndex = arcvec[iface];
          setFaceStuff(a_faces[iface], inputIndex, otherIndex, inputIV, otherIV, a_sd, a_domain);
        }

      }
    }
  }
  /*******************************/
  vector<EBIndex<XFACE> >
  CellNode::
  getXFaces(const  EBIndex<CELL> & a_vof,
            const Side::LoHiSide & a_sd,
            const Box            & a_domain) const
  {
    PR_assert(a_domain.contains(a_vof.m_pt));
    vector<EBIndex<XFACE> > retval;
    getFaceInfo(retval,  0, a_vof, a_sd, a_domain);
    
    return retval;
  }

  /*******************************/
  vector<EBIndex<YFACE> >
  CellNode::
  getYFaces(const  EBIndex<CELL>        & a_vof,
            const Side::LoHiSide & a_sd,
            const Box            & a_domain) const
  {
    PR_assert(a_domain.contains(a_vof.m_pt));

    vector<EBIndex<YFACE> > retval;
    getFaceInfo(retval,  1, a_vof, a_sd, a_domain);

    return retval;
  }

  /*******************************/
#if DIM==3
  vector<EBIndex<ZFACE> >
  CellNode::
  getZFaces(const  EBIndex<CELL> & a_vof,
            const Side::LoHiSide & a_sd,
            const Box            & a_domain) const
  {
    PR_assert(a_domain.contains(a_vof.m_pt));

    vector<EBIndex<ZFACE> > retval;
    getFaceInfo(retval, 2, a_vof, a_sd, a_domain);

    return retval;
  }
#endif
  vector< EBIndex<CELL> >
  CellNode::
  getVoFs(const Point& a_pt) const
  {
    vector<EBIndex<CELL> > retval(0);
    if(isCovered())
    {
      //nothing to return
    }
    else
    {
      if(isRegular())
      {
        retval.resize(1);
      }
      else
      {
        retval.resize(m_cellList->size());
      }
      for(int ivof = 0; ivof < retval.size(); ivof++)
      {
        retval[ivof] = EBIndex<CELL>(a_pt, ivof);
      }
    }
    return retval;
  }


  vector< EBIndex<BOUNDARY> >
  CellNode::
  getCutFaces(const Point& a_pt) const
  {
    vector<EBIndex<BOUNDARY> > retval;
    //only something to do if irregular
    if(isIrregular())
    {
      retval.resize(m_cellList->size());
      for(int ivof = 0; ivof < retval.size(); ivof++)
      {
        retval[0].m_pt = a_pt;
        retval[0].m_vofIDMe = ivof;
      }
    }
    return retval;
  }
 //this deals with the fine to coarse thing
  void
  EBGraph::
  coarsenVoFs(const Box     &  a_validBox,
              EBGraph       &  a_fineEBGraph,
              const Box     &  a_vofRegionCoar)
  {
    m_isDefined = true;
    m_validBox = a_validBox;
    m_domain = a_fineEBGraph.m_domain.coarsen(2);
    m_region = a_vofRegionCoar;
    Box regionFine = a_vofRegionCoar.refine(2) ;
    regionFine &= a_fineEBGraph.m_region;
    if(a_fineEBGraph.isCovered(regionFine))
    {
      m_tag = AllCovered;
    }
    else if(a_fineEBGraph.isRegular(regionFine))
    {
      m_tag = AllRegular;
    }
    else
    {
      m_tag = HasIrregular;
      m_internals = shared_ptr<localData>(new localData());
      m_internals->m_graph.define(m_region);
      Box cfregion = a_fineEBGraph.m_region.coarsen(2);
      Box region = m_region & cfregion;
      for(BoxIterator bit = region.begin(); bit != region.end(); ++bit)
      {
        Point pt = *bit;
        Box pointBox(pt, pt);
        Box fineBox = pointBox.refine(2);
        if(a_fineEBGraph.isRegular(fineBox))
        {
          m_internals->m_graph(pt, 0).defineAsRegular();
        }
        else if(a_fineEBGraph.isCovered(fineBox))
        {
          m_internals->m_graph(pt, 0).defineAsCovered();
        }
        else
        {
          CellNode& coarCellNode = m_internals->m_graph(pt, 0);
          coarCellNode.defineAsIrregular();
          m_internals->m_irregIVS |= pt;
          //this gets the collections of connected vofs
          vector< set<EBIndex<CELL> > > fineVoFSets
                = a_fineEBGraph.getVoFSets(fineBox);
          coarCellNode.m_cellList->resize(fineVoFSets.size());
          for (int iset = 0; iset < fineVoFSets.size(); iset++)
          {
            set<EBIndex<CELL> >& thisSet = fineVoFSets[iset];
            (*coarCellNode.m_cellList)[iset].m_finer.resize(thisSet.size());
            int icell = 0;
            for(auto iter = thisSet.begin(); iter != thisSet.end(); ++iter)
            {
              const EBIndex<CELL>& fineVoF = *iter;
              (*coarCellNode.m_cellList)[iset].m_finer[icell] = fineVoF;
              icell++;
              //if a fine cell is over a coarse multivalued cell,
              //redefine the fine cell as a full irregular cell and reset its coarser value
              //this stuff used to happen in fixFineToCoarse
              if( (iset > 0) && (a_fineEBGraph.isRegular(fineVoF.m_pt)))
              {
                a_fineEBGraph.m_internals->m_graph(fineVoF.m_pt,0).defineAsFullIrregular(fineVoF.m_pt, a_fineEBGraph.m_domain);
                a_fineEBGraph.m_internals->m_irregIVS |= fineVoF.m_pt;
              }
              if(a_fineEBGraph.isIrregular(fineVoF.m_pt))
              {
                vector<CellNode::IrrGraphNode>& vecirr = *a_fineEBGraph.m_internals->m_graph(fineVoF.m_pt, 0).m_cellList;
                vecirr[fineVoF.m_vofIDMe].m_coarser = iset;
              }
            }
          }
          
        }
      }
    }
  }

  /*******************************/
  ///used in coarsening--get collections of connected vofs
  /*******************************/
  vector< set<EBIndex<CELL> > >
  EBGraph::
  getVoFSets(const Box& a_region) const
  {
    //I do not know if this will be faster than the chombo 3 version
    //which just uses vectors but this is at least cleaner
    vector< set<EBIndex<CELL> > > retval;
    vector< EBIndex<CELL> > allVoFs;
    for(BoxIterator bit = a_region.begin(); bit != a_region.end(); ++bit)
    {
      Point pt = *bit;
      vector< EBIndex<CELL> > vofs = getVoFs(pt);
      allVoFs.insert(allVoFs.end(), vofs.begin(), vofs.end());
    }
    vector<bool> beenAdded(allVoFs.size(), false);
    for(int ivof = 0; ivof < allVoFs.size(); ivof++)
    {
      const EBIndex<CELL>& thisVoF = allVoFs[ivof];
      if(!beenAdded[ivof])
      {
        set<EBIndex<CELL > > thisVoFSet;
        thisVoFSet.insert(thisVoF);
        beenAdded[ivof]= true;
        bool doneAdding = false;
        while(!doneAdding)
        {
          //making temporary set because we need to iterate through thisvofset inside
          //and can be addding to it.
          set<EBIndex<CELL> > tempSet = thisVoFSet;
          ///i can start at ivof+1 here because 0 to ivof
          //has always been added.
          for (int jvof = ivof+1; jvof < allVoFs.size(); jvof++)
          {
            const EBIndex<CELL>& testVoF = allVoFs[jvof];
            for(auto iter = thisVoFSet.begin(); iter != thisVoFSet.end(); ++iter)
            {
              EBIndex<CELL> innerVoF = *iter;
              if(!beenAdded[jvof] && isConnected(innerVoF, testVoF))
              {
                tempSet.insert(testVoF);
                beenAdded[jvof] = true;
              }
            }
          }
          doneAdding = (tempSet.size() == thisVoFSet.size());
          if(!doneAdding)
          {
            thisVoFSet = tempSet;
          }
        }
        retval.push_back(thisVoFSet);
      }
    }
    return retval;
  }
  /*******************************/
  /// region coar is a cell centered region that represents
  ///how far out we know the cell centered data.
  /*******************************/
  void
  EBGraph::
  coarsenFaces(const EBGraph  &  a_fineEBGraph,
               const Box      &  a_regionCoar)  
  {
    for(BoxIterator bit = a_regionCoar.begin(); bit != a_regionCoar.end(); ++bit)
    {
      const Point& pt = *bit;
      if(isIrregular(pt))
      {
        vector< EBIndex<CELL> >vofsCoar = this->getVoFs(pt);
        for(int ivof = 0; ivof < vofsCoar.size(); ivof++)
        {
          CellNode::IrrGraphNode& node = (*m_internals->m_graph(pt, 0).m_cellList)[ivof];
          for(int idir = 0; idir < DIM; idir++)
          {
            for(SideIterator sit; sit.ok(); ++sit)
            {
              int iindex = ebp_index(idir, sit());
              node.m_arc[iindex] = coarsenFaces(a_fineEBGraph, a_regionCoar, vofsCoar[ivof], idir, sit());
            }
          }
        }
      }
    }
  }
  /*******************************/
  bool
  EBGraph::
  hasBoundaryFaces(const EBIndex<CELL> & a_vof, int a_idir, Side::LoHiSide a_sd) const
  {
    bool retval = false;
    if(a_idir == 0)
    {
      vector< EBIndex<XFACE> > faces = getXFaces(a_vof, a_sd);
      for(int ifac = 0; ifac < faces.size(); ifac++)
      {
        if(faces[ifac].m_isBoundary)
        {
          retval = true;
        }  
      }
    }
    else if(a_idir == 1)
    {
      vector< EBIndex<YFACE> > faces = getYFaces(a_vof, a_sd);
      for(int ifac = 0; ifac < faces.size(); ifac++)
      {
        if(faces[ifac].m_isBoundary)
        {
          retval = true;
        }  
      }
    }
#if DIM==3
    else if(a_idir == 2)
    {
      vector< EBIndex<ZFACE> > faces = getZFaces(a_vof, a_sd);
      for(int ifac = 0; ifac < faces.size(); ifac++)
      {
        if(faces[ifac].m_isBoundary)
        {
          retval = true;
        }  
      }
    }
#endif
    else
    {
      PROTO_ASSERT(false, "bogus idir");
    }
    return retval;
  }
  /*******************************/
  vector<int>
  EBGraph::
  coarsenFaces(const EBGraph        & a_fineEBGraph, 
               const Box            & a_regionCoar,
               const EBIndex<CELL>  & a_vofCoar,
               const int            & a_idir,
               const Side::LoHiSide & a_sd)
  {
    vector<int> retval;
    Point localPt  = a_vofCoar.m_pt;
    Point otherPt = localPt + sign(a_sd)*Point::Basis(a_idir);
    vector< EBIndex<CELL> > localFineVoFs = refine(a_vofCoar);
    if(!m_domain.contains(otherPt))
    {
      //boundary faces.
      //if there are any boundary faces on the fine level,
      //make one here too
      bool hasBoundaryFaceFine = false;
      for (int ithis = 0; ithis < localFineVoFs.size() && !hasBoundaryFaceFine ; ithis++)
      {
        const EBIndex<CELL> & fineVoF = localFineVoFs[ithis];
        if(a_fineEBGraph.hasBoundaryFaces(fineVoF, a_idir, a_sd))
        {
          hasBoundaryFaceFine = true;
        }   
      }
      if (hasBoundaryFaceFine)
      {
        //remember that -1 is the secret code for boundary face arcs
        retval.push_back(-1);
      }
    }
    else if(a_regionCoar.contains(otherPt))  
    {

      //if the other cell is in the domain, get all the vofs at that point
      //refine each one and see if it is connected to the vofs in the refine 
      //verison of this one.  if so, add the face.
      //this means we are walking the graph up to the finer level, seeing
      //what is connected there and making sure that there are coarse connections
      //where there are fine connections.

      //faces are coarsened one face in from boundary where vofs are coarsened
      //the follow line of code shows why this has to be.
      //if this causes failure, you need more ghost cells in your EBIS.
      vector< EBIndex<CELL> > otherCoarVoFs = getVoFs(otherPt);
      for(int iother = 0; iother < otherCoarVoFs.size(); iother++)
      {
        const EBIndex<CELL> otherCoarVoF = otherCoarVoFs[iother];
        bool addThisFace = false;
        vector< EBIndex<CELL> > otherFineVoFs = refine(otherCoarVoF);
        for(int iotherFine = 0; iotherFine < otherFineVoFs.size(); iotherFine++)
        {
          const EBIndex<CELL>& otherFineVoF = otherFineVoFs[iotherFine];
          for(int ilocalFine = 0; ilocalFine < localFineVoFs.size(); ilocalFine++)
          {
            const EBIndex<CELL> & localFineVoF = localFineVoFs[ilocalFine];
            if(a_fineEBGraph.isConnected(localFineVoF, otherFineVoF))
            {
              addThisFace = true;
            }
          }
        }
        if(addThisFace)
        {
          retval.push_back(otherCoarVoF.m_vofIDMe);
        }
      }
    }

    return retval;
  }
    
  /*******************************/
  void
  EBGraph::
  defineAsAllRegular(const Box& a_validBox, const Box& a_region, const Box & a_domain)
  {
    m_isDefined = true;
    m_validBox  = a_validBox;
    m_domain    = a_domain;
    m_region    = a_region;
    m_tag       = AllRegular;
    m_internals = shared_ptr<localData>();
  }

  /*******************************/
  void 
  EBGraph::
  defineAsAllCovered(const Box& a_validBox, const Box& a_region, const Box & a_domain)
  {

    m_isDefined = true;
    m_validBox  = a_validBox;
    m_domain    = a_domain;
    m_region    = a_region;
    m_tag       = AllCovered;
    m_internals = shared_ptr<localData>();
  }


/*******************************/
  void EBGraph::defineFromPrimitives(const Box                & a_validBox,
                                     const HostBoxData<int>   & a_regIrregCovered,
                                     const vector<IrregGraph> & a_irregGraph,
                                     const Box                & a_region,
                                     const Box                & a_domain,
                                     bool singleValuedCellsOnly)
  {
    PR_TIME("EBGraph::buildGraph");
    //this eliminates all the stuff in Chombo's version referring back to the original vector index and searching though them.
    //this simplifies the construction process enormously.
    //multivalued cells can still happen via coarsening.
    PROTO_ASSERT(singleValuedCellsOnly, "EBGraph::defineFromPrimitives is currently not general enough to handle multivalued inputs");
    m_isDefined = true;
    m_validBox  = a_validBox;
    m_domain    = a_domain;
    m_region    = a_region;
    m_tag       = HasIrregular;
    m_internals = shared_ptr<localData>(new localData());

    m_internals->m_graph.define(m_region);

    // set regular and covered cells
    int num_irreg = 0;
    for (BoxIterator bit = m_region.begin(); bit != m_region.end(); ++bit)
    {
      const Point& iv = *bit;
      if (a_regIrregCovered(iv, 0) == 1) //regular cell
      {
        m_internals->m_graph(iv, 0).defineAsRegular();
      }
      else if (a_regIrregCovered(iv, 0) == -1) // covered cell
      {
        m_internals->m_graph(iv, 0).defineAsCovered();
      }
      else if (a_regIrregCovered(iv, 0) == 0)
      {
        m_internals->m_graph(iv, 0).defineAsIrregular();
        num_irreg++;
      }
      else 
      {
        PR_error("invalid flag");
      }
    }
    //it can be greater if there are multivalued cells
    PROTO_ASSERT(a_irregGraph.size() >= num_irreg, "mismatch in ebgraph::define inputs");

    // now for irregular cells
    for (int ivecIrreg = 0; ivecIrreg < a_irregGraph.size(); ivecIrreg++)
    {
      const IrregGraph & inputNode = a_irregGraph[ivecIrreg];
      const Point     & iv = inputNode.m_cell;
      m_internals->m_irregIVS |= iv;

      CellNode::IrrGraphNode currNode;
      for (int idir = 0; idir < DIM; idir++)
      {
        for (SideIterator sit; sit.ok(); ++sit)
        {
          int irregIndex = ebp_index(idir, sit());
          const vector<int>& irregArcs = inputNode.m_arc[irregIndex];
          //comes from the singleValuedOnly thing
          PR_assert(irregArcs.size() <= 1);
          currNode.m_arc[irregIndex] = irregArcs;
//begin debugging checks          
          using std::cout;
          using std::endl;
          Point otheriv = iv + sign(sit())*Point::Basis(idir);
          Box graphbox = m_internals->m_graph.box();
          if(graphbox.contains(otheriv))
          {
            const CellNode& othernode = m_internals->m_graph(otheriv, 0);
            int ideb =0;
            if(othernode.isCovered() && (irregArcs.size() > 0))
            {
              ideb = 1;
            }
            if(othernode.isRegular() && (irregArcs.size() != 1))
            {
              ideb = 2;
            }
            if(ideb != 0)
            {
              cout << "found an inconsistency in incoming graph between " << iv << " and " << otheriv << " = " << ideb  << endl;
            }
          }
//end debug

        }
      }
      m_internals->m_graph(iv, 0).m_cellList->push_back(currNode);
    }
  }


/*******************************/
  void EBGraph::defineAsFullIrregular(const Box              & a_validBox,
                                      const Box              & a_region,
                                      const Box              & a_domain)
                                      
  {
    PR_TIME("EBGraph::buildGraph");

    m_isDefined = true;
    m_validBox = a_validBox;
    m_domain = a_domain;
    m_region = a_region;
    m_tag = HasIrregular;
    m_internals = shared_ptr<localData>(new localData());

    m_internals->m_graph.define(m_region);

    for (BoxIterator bit = m_region.begin(); bit != m_region.end(); ++bit)
    {
      m_internals->m_graph(*bit, 0).defineAsRegular();
    }
  }
/*******************************/
  void EBGraph::defineAsEmptyIrregular(const Box              & a_validBox,
                                       const Box              & a_region,
                                       const Box              & a_domain)
                                      
  {
    PR_TIME("EBGraph::buildGraph");

    m_isDefined = true;
    m_validBox = a_validBox;
    m_domain = a_domain;
    m_region = a_region;
    m_tag = HasIrregular;
    m_internals = shared_ptr<localData>(new localData());

    m_internals->m_graph.define(m_region);

    for (BoxIterator bit = m_region.begin(); bit != m_region.end(); ++bit)
    {
      m_internals->m_graph(*bit, 0).defineAsCovered();
    }
  }
/*******************************/
  vector< EBIndex<CELL> > EBGraph::getVoFs(const Point& a_iv) const
  {
    PR_assert(isDefined());
    PR_assert(m_region.contains(a_iv));

    vector< EBIndex<CELL> > retvec;
    if (m_tag == AllRegular)
    {
      retvec.push_back( EBIndex<CELL> (a_iv, 0));
    }
    else if (m_tag == AllCovered)
    {
      //return an empty vector
    }
    else if (m_tag == HasIrregular)
    {
      const CellNode& node = m_internals->m_graph(a_iv, 0);
      retvec = node.getVoFs(a_iv);
    }
    return retvec;
  }

/*******************************/
  vector< EBIndex<BOUNDARY> > EBGraph::getCutFaces(const Point& a_iv) const
  {
    PR_assert(isDefined());

    vector< EBIndex<BOUNDARY> > retvec;
    if (m_tag == HasIrregular) //otherwise return an empty vector
    {
      const CellNode& node = m_internals->m_graph(a_iv, 0);
      if(node.isIrregular())
      {
        retvec = node.getCutFaces(a_iv);
      }
    }
    return retvec;
  }

/*******************************/
  bool EBGraph::isRegular(const Point& a_iv) const
  {
    PR_assert(isDefined());
    bool retval;
    if (m_tag == AllRegular)
    {
      retval = true;
    }
    else if (m_tag == AllCovered)
    {
      retval = false;
    }
    else if (m_tag == HasIrregular)
    {
      const CellNode& node = m_internals->m_graph(a_iv, 0);
      retval = node.isRegular();
    }
    else
    {
      retval = false;
      PR_error("EBGraph::isRegular:Bogus Tag");
    }
    return retval;
  }


  template<CENTERING cent>
  bool
  EBGraph::
  isFaceRegular(const EBIndex<cent>& a_face) const
  {
    EBIndex<CELL> voflo = a_face.getVoF(Side::Lo);
    EBIndex<CELL> vofhi = a_face.getVoF(Side::Hi);
    bool retval = false;
    if(m_domain.contains(voflo.m_pt) && isRegular(voflo.m_pt))
    {
      retval = true;
    }
    if(m_domain.contains(vofhi.m_pt) && isRegular(vofhi.m_pt))
    {
      retval = true;
    }
    return retval;
  }
  /****/
  template<CENTERING cent>
  bool
  EBGraph::
  isFaceCovered(const EBIndex<cent>& a_face) const
  {
    EBIndex<CELL> voflo = a_face.getVoF(Side::Lo);
    EBIndex<CELL> vofhi = a_face.getVoF(Side::Hi);
    bool retval = false;
    if(m_domain.contains(voflo.m_pt) && isCovered(voflo.m_pt))
    {
      retval = true;
    }
    if(m_domain.contains(vofhi.m_pt) && isCovered(vofhi.m_pt))
    {
      retval = true;
    }
    return retval;
  }
/*******************************/
  bool EBGraph::isIrregular(const Point& a_iv) const
  {
    PR_assert(isDefined());

    bool retval;
    if (m_tag == AllRegular)
    {
      retval = false;
    }
    else if (m_tag == AllCovered)
    {
      retval = false;
    }
    else if (m_tag == HasIrregular)
    {
      const CellNode& node = m_internals->m_graph(a_iv, 0);
      retval = node.isIrregular();
    }
    else
    {
      retval = false;
      PR_error("EBGraph::isIrregular:Bogus Tag");
    }
    return retval;
  }


/*******************************/
  bool EBGraph::isCovered(const Point& a_iv) const
  {
    PR_assert(isDefined());

    bool retval;
    if (m_tag == AllRegular)
    {
      retval = false;
    }
    else if (m_tag == AllCovered)
    {
      retval = true;
    }
    else if (m_tag == HasIrregular)
    {
      const CellNode& node = m_internals->m_graph(a_iv, 0);
      retval = node.isCovered();
    }
    else
    {
      retval = false;
      PR_error("EBGraph::isIrregular:Bogus Tag");
    }
    return retval;
  }

/*******************************/
  bool EBGraph::isCovered(const Box& a_box) const
  {
    PR_assert(isDefined());

    bool retval;
    if (m_tag == AllRegular)
    {
      retval = false;
    }
    else if (m_tag == AllCovered)
    {
      retval = true;
    }
    else if (m_tag == HasIrregular)
    {

      retval = true;
      BoxIterator bit(a_box);
      for (BoxIterator bit = a_box.begin(); bit != a_box.end(); ++bit)
      {
        if (!isCovered(*bit))
        {
          retval = false;
          break;
        }
      }
    }
    else
    {
      retval = false;
      PR_error("EBGraph::isCovered:Bogus Tag");
    }
    return retval;
  }

/*******************************/
  bool EBGraph::isRegular(const Box& a_box) const
  {
    PR_assert(isDefined());

    bool retval;
    if (m_tag == AllRegular)
    {
      retval = true;
    }
    else if (m_tag == AllCovered)
    {
      retval = false;
    }
    else if (m_tag == HasIrregular)
    {
      retval = true;
      for (BoxIterator bit = a_box.begin(); bit != a_box.end(); ++bit)
      {
        if (!isRegular(*bit))
        {
          retval = false;
          break;
        }
      }
    }
    else
    {
      retval = false;
      PR_error("EBGraph::isRegular:Bogus Tag");
    }
    return retval;
  }
/*******************************/
  vector< EBIndex<XFACE> > 
  EBGraph::
  getXFaces(const  EBIndex<CELL>        & a_vof,
            const Side::LoHiSide & a_sd) const
  {
    PR_assert(isDefined());

    vector< EBIndex<XFACE> > retval;
    if (m_tag == AllRegular)
    {
      int idir = 0;
      Point otherIV = a_vof.m_pt
        + sign(a_sd)*Point::Basis(idir);
      int otherCellInd = 0;
      if (!m_domain.contains(otherIV))
      {
        otherCellInd = -1;
      }
      retval.resize(1);
      if(a_sd == Side::Lo)
      {
        retval[0].m_pt      = a_vof.m_pt;
        retval[0].m_vofIDMe = a_vof.m_vofIDMe;
        retval[0].m_vofIDLo = otherCellInd;
      }
      else
      {
        retval[0].m_pt      = otherIV;
        retval[0].m_vofIDMe = otherCellInd;
        retval[0].m_vofIDLo = a_vof.m_vofIDMe;
      }
      retval[0].m_isBoundary = (!m_domain.contains(otherIV));
    }
    else if (m_tag == AllCovered)
    {
      //return empty vector
    }
    else if (m_tag == HasIrregular)
    {
      const Point& iv = a_vof.m_pt;
      const CellNode& node = m_internals->m_graph(iv, 0);
      retval = node.getXFaces(a_vof, a_sd, m_domain);
    }

    return retval;
  }


/*******************************/
  vector< EBIndex<YFACE> > 
  EBGraph::
  getYFaces(const  EBIndex<CELL>        & a_vof,
            const Side::LoHiSide & a_sd) const
  {
    PR_assert(isDefined());

    vector< EBIndex<YFACE> > retval;
    if (m_tag == AllRegular)
    {
      int idir = 1;
      Point otherIV = a_vof.m_pt
        + sign(a_sd)*Point::Basis(idir);
      int otherCellInd = 0;
      if (!m_domain.contains(otherIV))
      {
        otherCellInd = -1;
      }
      retval.resize(1);
      if(a_sd == Side::Lo)
      {
        retval[0].m_pt      = a_vof.m_pt;
        retval[0].m_vofIDMe = a_vof.m_vofIDMe;
        retval[0].m_vofIDLo = otherCellInd;
      }
      else
      {
        retval[0].m_pt      = otherIV;
        retval[0].m_vofIDMe = otherCellInd;
        retval[0].m_vofIDLo = a_vof.m_vofIDMe;
      }
      retval[0].m_isBoundary = (!m_domain.contains(otherIV));
    }
    else if (m_tag == AllCovered)
    {
      //return empty vector
    }
    else if (m_tag == HasIrregular)
    {
      const Point& iv = a_vof.m_pt;
      const CellNode& node = m_internals->m_graph(iv, 0);
      retval = node.getYFaces(a_vof, a_sd, m_domain);
    }

    return retval;
  }
#if DIM==3
/*******************************/
  vector< EBIndex<ZFACE> > 
  EBGraph::
  getZFaces(const  EBIndex<CELL>        & a_vof,
            const Side::LoHiSide & a_sd) const
  {
    PR_assert(isDefined());

    vector< EBIndex<ZFACE> > retval;
    if (m_tag == AllRegular)
    {
      int idir = 2;
      Point otherIV = a_vof.m_pt
        + sign(a_sd)*Point::Basis(idir, 1);
      int otherCellInd = 0;
      if (!m_domain.contains(otherIV))
      {
        otherCellInd = -1;
      }
      retval.resize(1);
      if(a_sd == Side::Lo)
      {
        retval[0].m_pt      = a_vof.m_pt;
        retval[0].m_vofIDMe = a_vof.m_vofIDMe;
        retval[0].m_vofIDLo = otherCellInd;
      }
      else
      {
        retval[0].m_pt      = otherIV;
        retval[0].m_vofIDMe = otherCellInd;
        retval[0].m_vofIDLo = a_vof.m_vofIDMe;
      }
      retval[0].m_isBoundary = (!m_domain.contains(otherIV));
    }
    else if (m_tag == AllCovered)
    {
      //return empty vector
    }
    else if (m_tag == HasIrregular)
    {
      const Point& iv = a_vof.m_pt;
      const CellNode& node = m_internals->m_graph(iv, 0);
      retval = node.getZFaces(a_vof, a_sd, m_domain);
    }

    return retval;
  }
#endif


/*******************************/
  bool EBGraph::isConnected(const  EBIndex<CELL> & a_vof1,
                            const  EBIndex<CELL> & a_vof2) const
  {
    PR_TIME("EBGraph::isConnected");

    PR_assert(isDefined());


    const Point& iv1 = a_vof1.m_pt;
    const Point& iv2 = a_vof2.m_pt;
    if(iv1 == iv2) return false;

    PR_assert(m_region.contains(iv1));
    PR_assert(m_region.contains(iv2));
    EBIndex<CELL>  vofLo, vofHi;
    //use the lexical greaterThan because it will
    //always give the correct answer since only one component
    //is different.
    if (iv1 < (iv2))
    {
      vofLo = a_vof1;
      vofHi = a_vof2;
    }
    else
    {
      vofLo = a_vof2;
      vofHi = a_vof1;
    }
    int direction;
    bool dirfound;
    const Point& ivLo = vofLo.m_pt;
    const Point& ivHi = vofHi.m_pt;

    dirfound = false;
    for (int idir = 0; ((idir<DIM) && !dirfound); idir++)
    {
      if ((ivHi - ivLo) == Point::Basis(idir, 1))
      {
        direction = idir;
        dirfound = true;
      }
    }
    //if not neigboring intvects, no way it can be connected
    if (!dirfound) return false;

    bool voffound = false;
    if(direction == 0)
    {
      vector<EBIndex<XFACE> > faces = getXFaces(vofLo, Side::Hi);
      for (int iface = 0; iface < faces.size(); iface++)
      {
        const EBIndex<XFACE>& face = faces[iface];
        if (face.m_vofIDMe == vofHi.m_vofIDMe)
        {
          voffound = true;
        }
      }
    }
    else if(direction == 1)
    {
      vector<EBIndex<YFACE> > faces = getYFaces(vofLo, Side::Hi);
      for (int iface = 0; iface < faces.size(); iface++)
      {
        const EBIndex<YFACE>& face = faces[iface];
        if (face.m_vofIDMe == vofHi.m_vofIDMe)
        {
          voffound = true;
        }
      }
    }
#if DIM==3
    else if(direction == 2)
    {
      vector<EBIndex<ZFACE> > faces = getZFaces(vofLo, Side::Hi);
      for (int iface = 0; iface < faces.size(); iface++)
      {
        const EBIndex<ZFACE>& face = faces[iface];
        if (face.m_vofIDMe == vofHi.m_vofIDMe)
        {
          voffound = true;
        }
      }
    }
#endif
    else
    {
      PR_error("bogus direction");
    }
    return voffound;
  }

}

#endif
