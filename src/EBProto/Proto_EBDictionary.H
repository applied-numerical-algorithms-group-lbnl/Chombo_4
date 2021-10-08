
#ifndef _Proto_EBDictionary_H_
#define _Proto_EBDictionary_H_


#include "Proto_RealVect.H"
#include "Proto_Box.H"
#include "Proto_Point.H"
#include "implem/Proto_LSquares.H"
#include "Proto_BaseIF.H"
#include "implem/Proto_IrregNode.H"
#include "implem/Proto_CutCellMoments.H"
#include "Proto_HostBoxData.H"
#include "implem/Proto_HostIrregData.H"
#include "implem/Proto_PointSet.H"
#include "implem/Proto_LevelData.H"
#include "Proto_EBGraph.H"
#include "Proto_EBStencil.H"
#include "EBProto.H"

using std::shared_ptr;
using std::vector;

namespace Proto
{
  /// class that creates and stores EBStencils
  /**
     class that creates and stores EBStencils
     You register a particular stencil at a particular refinement
     and after that you can always retrieve it.
     This is a simple implementation.  The Chombo version is better.
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
    typedef HostIrregData<CELL    ,  IndMomDIM , 1>  VoluData;
    typedef HostIrregData<BOUNDARY,  IndMomDIM , 1>  EBFaData;
    typedef HostIrregData<XFACE, IndMomSDMinOne, 1>  XFacData;
    typedef HostIrregData<YFACE, IndMomSDMinOne, 1>  YFacData;
    typedef HostIrregData<ZFACE, IndMomSDMinOne, 1>  ZFacData;

    ///
    /**
       geoserv is the holder for all momeent information \\
       a_grids is the layout over which that data is held. \\
       a_ghost is the number of ghost cells your data will hold.
       It must be at least as big as the stencil span in each direction..
    */
    EBDictionary(const shared_ptr<GeometryService<order> >  & a_geoserv,
                 const DisjointBoxLayout                    & a_grids,
                 Point  a_srcGhost,  Point  a_dstGhost, data_t a_dx,
                 bool a_periodiconly);


    ///
    /**
       This is how you store a stencil.   The string must match one of the appropriate centering in EBStencilArchive.
       stencilName is the name of the stencil \\
       domainBC name is the name of the domain boundary condition to be applied (this can be turned off) \\
       domainBC name is the name of the Cut cell boundary condition to be applied (this can be turned off)

     */
    inline void registerStencil(string a_stencilName,
                                string a_domainBCName,
                                string a_ebbcName);


    ///
    /** Retrieve the stencil for a particular box at a particular level of refinement.  The input strings must match the register call */
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
    data_t m_dx;
    bool   m_periodiconly;
  private:
    EBDictionary();
    EBDictionary(  const EBDictionary<order, data_t, srcCenter, dstCenter>& a_inp);
    void operator=(const EBDictionary<order, data_t, srcCenter, dstCenter>& a_inp);


  };
}
#include "implem/Proto_EBDictionaryImplem.H"
#endif
