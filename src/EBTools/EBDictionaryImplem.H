
#ifndef _Proto_EBDictionaryImplem_H_
#define _Proto_EBDictionaryImplem_H_



namespace Proto
{
  ///
  /**
     a_ghost is the number of ghost cells your data will hold.
     It must be at least as big as the stencil span in each direction..
  */
  template <int order, typename data_t, CENTERING srcCenter, CENTERING dstCenter>
  EBDictionary<order, data_t, srcCenter, dstCenter>::
  EBDictionary(const shared_ptr<GeometryService<order> >    & a_geoserv,
               const DisjointBoxLayout                      & a_grids,
               Point  a_srcGhost,  Point  a_dstGhost,
               data_t a_dx,  bool a_periodiconly)
  {
    m_geoserv  = a_geoserv;
    m_grids    = a_grids;
    m_srcGhost = a_srcGhost;
    m_dstGhost = a_dstGhost;
    m_dx       = a_dx;
    m_periodiconly     = a_periodiconly;
    m_isClosed = false;
    m_registeredStencils.resize(m_grids.size());
  }    
    
  ///
  template <int order, typename data_t, CENTERING srcCenter, CENTERING dstCenter>
  void
  EBDictionary<order, data_t, srcCenter, dstCenter>::
  registerStencil(string a_stencilName,
                  string a_domainBCName,
                  string a_ebbcName)
  {
    PR_assert(!m_isClosed);

    stencilKey key;
    key.m_ebbcName    = a_ebbcName;
    key.m_stencilName = a_stencilName;
    Box domain = m_grids.getDomain();
    shared_ptr< LevelData< EBGraph  >  > graphs = m_geoserv->getGraphs(domain);
    shared_ptr< LevelData< VoluData >  > voldat = m_geoserv->getVoluData(domain);
    shared_ptr< LevelData< EBFaData >  > ebfdat = m_geoserv->getEBFaceData(domain);
    shared_ptr< LevelData< XFacData >  > xfadat = m_geoserv->getXFaceData(domain);
    shared_ptr< LevelData< YFacData >  > yfadat = m_geoserv->getYFaceData(domain);
    shared_ptr< LevelData< ZFacData >  > zfadat = m_geoserv->getZFaceData(domain);


    for(int ibox = 0; ibox < m_grids.size(); ibox++)
    {
      shared_ptr<ebstencil_t> 
        newsten(new ebstencil_t(a_stencilName, a_domainBCName, a_ebbcName, 
                                m_grids[ibox], domain, m_srcGhost, m_dstGhost,
                                (*graphs)[ibox], (*voldat)[ibox], (*ebfdat)[ibox],
                                (*xfadat)[ibox], (*yfadat)[ibox], (*zfadat)[ibox],
                                m_dx, m_periodiconly));
      std::map<stencilKey, shared_ptr<ebstencil_t> >& boxmap = m_registeredStencils[ibox];
      boxmap[key] = newsten;
    }
  }



  ///
  template <int order, typename data_t, CENTERING srcCenter, CENTERING dstCenter>
  shared_ptr< EBStencil<order, data_t, srcCenter, dstCenter> >
  EBDictionary<order, data_t, srcCenter, dstCenter>::
  getEBStencil(const string & a_stencilName,
               const string & a_ebbcName, 
               const int    & a_boxid)
  {
    stencilKey key;
    key.m_ebbcName    = a_ebbcName;
    key.m_stencilName = a_stencilName;
    shared_ptr<ebstencil_t>   retval;
    if(m_registeredStencils[a_boxid].find(key) == m_registeredStencils[a_boxid].end())
    {
      PR_error("stencil not found");
    }
    else
    {
      retval = m_registeredStencils[a_boxid][key];
    }
    return retval;
  }
}
#endif
