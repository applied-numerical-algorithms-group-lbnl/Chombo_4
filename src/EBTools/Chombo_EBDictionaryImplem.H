
#ifndef _Chombo_EBDictionaryImplem_H_
#define _Chombo_EBDictionaryImplem_H_



///
/**
   a_ghost is the number of ghost cells your data will hold.
   It must be at least as big as the stencil span in each direction..
*/
template <int order, typename data_t, CENTERING srcCenter, CENTERING dstCenter>
EBDictionary<order, data_t, srcCenter, dstCenter>::
EBDictionary(const shared_ptr<GeometryService<order> >    & a_geoserv,
             const DisjointBoxLayout                      & a_grids,
             const Box                                   & a_domain,
             data_t a_dx, Point  a_ghost)
{
  m_geoserv  = a_geoserv;
  m_grids .resize(1, a_grids );
  m_domain.resize(1, a_domain);
  m_dx    .resize(1, a_dx    );
  
  m_srcGhost = a_ghost;
  m_dstGhost = a_ghost;
  m_isClosed = false;

  DataIterator dit = m_grids[0].dataIterator();
  m_registeredStencils.resize(1);
  m_registeredStencils[0].resize(dit.size());
}    


template <int order, typename data_t, CENTERING srcCenter, CENTERING dstCenter>
EBDictionary<order, data_t, srcCenter, dstCenter>::
EBDictionary(const shared_ptr<GeometryService<order> >   & a_geoserv,     
             const vector<DisjointBoxLayout>             & a_grids,       
             const vector<Box>                           & a_domain,      
             const vector<data_t>                        & a_dx,          
             Point  a_ghost)
{
  m_geoserv  =  a_geoserv;
  m_grids    =  a_grids  ;
  m_domain   =  a_domain ;
  m_dx       =  a_dx     ;
  
  m_srcGhost = a_ghost;
  m_dstGhost = a_ghost;
  m_isClosed = false;

  m_registeredStencils.resize(m_grids.size());
  for(int ilev = 0; ilev < m_grids.size(); ilev++)
  {
    DataIterator dit = m_grids[ilev].dataIterator();
    m_registeredStencils[ilev].resize(dit.size());
  }
}    

///
template <int order, typename data_t, CENTERING srcCenter, CENTERING dstCenter>
void
EBDictionary<order, data_t, srcCenter, dstCenter>::
registerStencil(string a_stencilName,
                string a_ebbcName,
                const vector<EBIndex<dstCenter> >               & a_dstVoFs,
                const vector<LocalStencil<srcCenter,  data_t> > & a_stencil,
                const Proto::Box                                & a_srcValid,
                const Proto::Box                                & a_dstValid,
                const Proto::Box                                & a_srcDomain,
                const Proto::Box                                & a_dstDomain,
                const Point                                     & a_srcGhost,
                const Point                                     & a_dstGhost,
                bool                                              a_needDiagonalWeights,
                unsigned int                                      a_ibox)
{
  PR_assert(!m_isClosed);

  stencilKey key;
  key.m_ebbcName    = a_ebbcName;
  key.m_stencilName = a_stencilName;
  Box ch_srcDom = ProtoCh::getBox(a_srcDomain);
  Box ch_dstDom = ProtoCh::getBox(a_dstDomain);
  int dstIndex = getIndex(ch_dstDom);

  
  auto ldgraphS = m_geoserv->getGraphs(ch_srcDom);
  auto ldgraphD = m_geoserv->getGraphs(ch_dstDom);

  DataIterator dit = m_grids[dstIndex].dataIterator();

  const auto& srcGraph = (*ldgraphS)[dit[a_ibox]];
  const auto& dstGraph = (*ldgraphD)[dit[a_ibox]];
  Proto::Box srcValid  = srcGraph.validBox();
  Proto::Box dstValid  = dstGraph.validBox();
  std::map<stencilKey, shared_ptr<ebstencil_t> >& boxmap = m_registeredStencils[dstIndex][a_ibox];
  shared_ptr<ebstencil_t>
    newsten(new ebstencil_t(a_dstVoFs, a_stencil, 
                            srcGraph, dstGraph, srcValid, dstValid, a_srcDomain, a_dstDomain,
                            a_srcGhost, a_dstGhost,  a_needDiagonalWeights));
  boxmap[key] = newsten;
  
}
///
template <int order, typename data_t, CENTERING srcCenter, CENTERING dstCenter>
void
EBDictionary<order, data_t, srcCenter, dstCenter>::
registerStencil(string a_stencilName,
                string a_domainBCName[2*DIM],
                string a_ebbcName,
                Box    a_srcDomain,
                Box    a_dstDomain,
                bool   a_needDiagonalWeights,
                Point  a_dstGrow)
{
  PR_assert(!m_isClosed);

  stencilKey key;
  key.m_ebbcName    = a_ebbcName;
  key.m_stencilName = a_stencilName;
  int dstIndex = getIndex(a_dstDomain);

  
  auto graphS = m_geoserv->getGraphs(    a_srcDomain);
  auto graphD = m_geoserv->getGraphs(    a_dstDomain);
  auto voldat = m_geoserv->getVoluData(  a_srcDomain);
  auto ebfdat = m_geoserv->getEBFaceData(a_srcDomain);
  auto xfadat = m_geoserv->getXFaceData( a_srcDomain);
  auto yfadat = m_geoserv->getYFaceData( a_srcDomain);
  auto zfadat = m_geoserv->getZFaceData( a_srcDomain);

  DataIterator dit = m_grids[dstIndex].dataIterator();
  for(int ibox = 0; ibox < dit.size(); ibox++)
  {
    Proto::Box gridbxS  = (*graphS)[dit[ibox]].validBox();
    Proto::Box gridbxD  = (*graphD)[dit[ibox]].validBox();
    Proto::Box domainbxS = ProtoCh::getProtoBox(a_srcDomain);
    Proto::Box domainbxD = ProtoCh::getProtoBox(a_dstDomain);
    std::map<stencilKey, shared_ptr<ebstencil_t> >& boxmap = m_registeredStencils[dstIndex][ibox];
    shared_ptr<ebstencil_t> 
      newsten
	 ( new ebstencil_t(a_stencilName, a_domainBCName, a_ebbcName, 
                              gridbxS, gridbxD, domainbxS, domainbxD,
                              m_srcGhost, m_dstGhost,
                              (*graphS)[dit[ibox]],(*graphD)[dit[ibox]],
                              (*voldat)[dit[ibox]], (*ebfdat)[dit[ibox]],
                              (*xfadat)[dit[ibox]], (*yfadat)[dit[ibox]], (*zfadat)[dit[ibox]],
                              m_dx[dstIndex], false, a_needDiagonalWeights, a_dstGrow));  //false is to turn off eriodic

 
    boxmap[key] = newsten;
  }
}



///
template <int order, typename data_t, CENTERING srcCenter, CENTERING dstCenter>
shared_ptr< EBStencil<order, data_t, srcCenter, dstCenter> >
EBDictionary<order, data_t, srcCenter, dstCenter>::
getEBStencil(const string & a_stencilName,
             const string & a_ebbcName, 
             const Box    & a_srcDomain,
             const Box    & a_dstDomain,               
             const int    & a_boxid)
{
  stencilKey key;
  key.m_ebbcName    = a_ebbcName;
  key.m_stencilName = a_stencilName;
  shared_ptr<ebstencil_t>   retval;
  int dstIndex = getIndex(a_dstDomain);
  
  if(m_registeredStencils[dstIndex][a_boxid].find(key) == m_registeredStencils[dstIndex][a_boxid].end())
  {
    PR_error("stencil not found");
  }
  else
  {
    retval = m_registeredStencils[dstIndex][a_boxid][key];
  }
  return retval;
}

#endif
