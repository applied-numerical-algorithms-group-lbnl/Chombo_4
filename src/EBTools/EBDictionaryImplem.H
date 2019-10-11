
#ifndef _Chomboo_EBDictionaryImplem_H_
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
             Point  a_srcGhost,  Point  a_dstGhost, data_t a_dx)
{
  m_geoserv  = a_geoserv;
  m_grids .resize(1, a_grids );
  m_domain.resize(1, a_domain);
  m_dx    .resize(1, a_dx    );
  
  m_srcGhost = a_srcGhost;
  m_dstGhost = a_dstGhost;
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
             Point  a_srcGhost,  Point  a_dstGhost)                        
{
  m_geoserv  =  a_geoserv;
  m_grids    =  a_grids  ;
  m_domain   =  a_domain ;
  m_dx       =  a_dx     ;
  
  m_srcGhost = a_srcGhost;
  m_dstGhost = a_dstGhost;
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
                string a_domainBCName,
                string a_ebbcName,
                Box    a_srcDomain,
                Box    a_dstDomain,
                bool   a_needDiagonalWeights)
{
  PR_assert(!m_isClosed);

  stencilKey key;
  key.m_ebbcName    = a_ebbcName;
  key.m_stencilName = a_stencilName;
  int dstIndex = getIndex(a_dstDomain);

  
  shared_ptr< LevelData< EBGraph  >  > graphS = m_geoserv->getGraphs(    a_srcDomain);
  shared_ptr< LevelData< EBGraph  >  > graphD = m_geoserv->getGraphs(    a_dstDomain);
  shared_ptr< LevelData< VoluData >  > voldat = m_geoserv->getVoluData(  a_srcDomain);
  shared_ptr< LevelData< EBFaData >  > ebfdat = m_geoserv->getEBFaceData(a_srcDomain);
  shared_ptr< LevelData< XFacData >  > xfadat = m_geoserv->getXFaceData( a_srcDomain);
  shared_ptr< LevelData< YFacData >  > yfadat = m_geoserv->getYFaceData( a_srcDomain);
  shared_ptr< LevelData< ZFacData >  > zfadat = m_geoserv->getZFaceData( a_srcDomain);

  DataIterator dit = m_grids[dstIndex].dataIterator();
  for(int ibox = 0; ibox < dit.size(); ibox++)
  {
    Proto::Box gridbxS  = (*graphS)[dit[ibox]].validBox();
    Proto::Box gridbxD  = (*graphD)[dit[ibox]].validBox();
    Proto::Box domainbxS = ProtoCh::getProtoBox(a_srcDomain);
    Proto::Box domainbxD = ProtoCh::getProtoBox(a_dstDomain);
    shared_ptr<ebstencil_t> 
      newsten(new ebstencil_t(a_stencilName, a_domainBCName, a_ebbcName, 
                              gridbxS, gridbxD, domainbxS, domainbxD,
                              m_srcGhost, m_dstGhost,
                              (*graphS)[dit[ibox]],(*graphD)[dit[ibox]],
                              (*voldat)[dit[ibox]], (*ebfdat)[dit[ibox]],
                              (*xfadat)[dit[ibox]], (*yfadat)[dit[ibox]], (*zfadat)[dit[ibox]],
                              m_dx[dstIndex], false, a_needDiagonalWeights));  //false is to turn off periodic
    std::map<stencilKey, shared_ptr<ebstencil_t> >& boxmap = m_registeredStencils[dstIndex][ibox];
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
