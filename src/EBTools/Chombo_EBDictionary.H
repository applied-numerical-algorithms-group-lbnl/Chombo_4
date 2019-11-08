
#ifndef _Chombo_EBDictionary_H_
#define _Chombo_EBDictionary_H_


#include "Proto.H"
#include "EBProto.H"
#include "Chombo_LevelData.H"
#include "Chombo_EBChombo.H"
#include "Chombo_GeometryService.H"
#include "Chombo_NamespaceHeader.H"


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
     This interface is for a single level
     a_ghost is the number of ghost cells your data will hold.
     It must be at least as big as the stencil span in each direction..
  */
  EBDictionary(const shared_ptr<GeometryService<order> >  & a_geoserv,
               const DisjointBoxLayout                    & a_grids,
               const Box                                  & a_domain,
               data_t a_dxPoint,   Point a_srcGhost,  Point  a_dstGhost);


  ///maybe we need an amr structure.
  EBDictionary(const shared_ptr<GeometryService<order> >   & a_geoserv,
               const vector<DisjointBoxLayout>             & a_grids,
               const vector<Box>                           & a_domain,
               const vector<data_t>                        & a_dx,
               Point  a_srcGhost,  Point  a_dstGhost);                      
  

  int getIndex(const Box& a_dstDomain)
  {
    bool found = false;
    int retval =-1;
    for(int ilev = 0; ilev < m_domain.size(); ilev++)
    {
      if(a_dstDomain == m_domain[ilev])
      {
        found = true;
        retval = ilev;
      }
    }
    PROTO_ASSERT(found, "dictionary did not find input domain");
    return retval;
  }
  
  ///
  /**
     diagonal weights are for relaxation--expensive to compute so leave them out if you do not need them
     dest grow is how many cells to grow valid box for destination.
  */
  inline void registerStencil(string a_stencilName,
                              string a_domainBCName,
                              string a_ebbcName,
                              Box    a_srcDomain,
                              Box    a_dstDomain,
                              bool needDiagonalWeights = false,
                              Point  a_destGrow = Point::Zeros());


  ///get the stencil for a particular box at a particular level of refinement
  inline shared_ptr< ebstencil_t >
  getEBStencil(const string & a_stencilName, 
               const string & a_ebbcName,
               const Box    & a_srcDomain,
               const Box    & a_dstDomain,               
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
  vector< vector< map< stencilKey, shared_ptr< ebstencil_t> > > >   m_registeredStencils;
    
  shared_ptr<GeometryService<order> > m_geoserv;
  Point                               m_srcGhost;
  Point                               m_dstGhost;
  bool                                m_isClosed;

  vector<DisjointBoxLayout>           m_grids;
  vector<Box>                         m_domain;
  vector<data_t>                      m_dx;
private:
  EBDictionary();
  EBDictionary(  const EBDictionary<order, data_t, srcCenter, dstCenter>& a_inp);
  void operator=(const EBDictionary<order, data_t, srcCenter, dstCenter>& a_inp);


};

#include "Chombo_EBDictionaryImplem.H"
#include "Chombo_NamespaceFooter.H"
#endif
