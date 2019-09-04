
#ifndef _Chombo_EBDictionary_H_
#define _Chombo_EBDictionary_H_


#include "LevelData.H"
#include "Proto.H"
#include "EBProto.H"
#include "EBChombo.H"
#include "GeometryService.H"

#include "NamespaceHeader.H"


///
/**
   makes geometric description from an implicit function.
*/
template <int order, typename data_t, CENTERING srcCenter, CENTERING dstCenter>
class EBDictionary
{

public:
  typedef EBStencil<order, data_t, srcCenter, dstCenter>  ebstencil_t;

  typedef IndexedMoments<DIM  , order> IndMomDIM;
  typedef IndexedMoments<DIM-1, order> IndMomSDMinOne;
  typedef MomentIterator<DIM  , order> MomItDIM;
  typedef MomentIterator<DIM-1, order> MomItSDMinOne;
  typedef IrregNode<order> inode;
  typedef HostIrregData<CELL,      IndMomDIM     , 1>  VoluData;
  typedef HostIrregData<BOUNDARY,  IndMomDIM     , 1>  EBFaData;
  typedef HostIrregData<XFACE, IndMomSDMinOne, 1>  XFacData;
  typedef HostIrregData<YFACE, IndMomSDMinOne, 1>  YFacData;
  typedef HostIrregData<ZFACE, IndMomSDMinOne, 1>  ZFacData;

  ///
  /**
     a_ghost is the number of ghost cells your data will hold.
     It must be at least as big as the stencil span in each direction..
  */
  EBDictionary(const shared_ptr<GeometryService<order> >  & a_geoserv,
               const DisjointBoxLayout                    & a_grids,
               const Box                                  & a_domain,
               Point  a_srcGhost,  Point  a_dstGhost, data_t a_dx);


  inline void registerStencil(string a_stencilName,
                              string a_domainBCName,
                              string a_ebbcName);


  ///get the stencil for a particular box at a particular level of refinement
  inline shared_ptr< ebstencil_t >
  getEBStencil(const string & a_stencilName, 
               const string & a_ebbcName, 
               const int    & a_boxid);

  ///
  inline int numLevels() const 
  {
    return m_registeredStencils.size();
  }

  ///after this is done, you cannot register any more stencils.  This will delete the moment information.
  void  close()
  {
    m_geoserv = shared_ptr<GeometryService<order> >();
  }
protected:


  struct stencilKey
  {
    string m_ebbcName;
    string m_stencilName;
    bool operator< (const stencilKey& a_key) const
    {
      bool retval = false;
      if(m_stencilName == a_key.m_stencilName)
      {
        retval= (m_ebbcName < a_key.m_ebbcName);
      }
      else
      {
        retval= (m_stencilName < a_key.m_stencilName);
      }
      return retval;
    }
  };
 
  /// one per box in the grid
  vector< map< stencilKey, shared_ptr< ebstencil_t> > >     m_registeredStencils;
    
  shared_ptr<GeometryService<order> > m_geoserv;
  Point                               m_srcGhost;
  Point                               m_dstGhost;
  DisjointBoxLayout                   m_grids;
  bool                                m_isClosed;
  Box                                 m_domain;
  data_t m_dx;
private:
  EBDictionary();
  EBDictionary(  const EBDictionary<order, data_t, srcCenter, dstCenter>& a_inp);
  void operator=(const EBDictionary<order, data_t, srcCenter, dstCenter>& a_inp);


};

#include "EBDictionaryImplem.H"
#include "NamespaceFooter.H"
#endif
