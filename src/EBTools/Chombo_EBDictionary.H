
#ifndef _Chombo_EBDictionary_H_
#define _Chombo_EBDictionary_H_


#include "Proto.H"
#include "EBProto.H"
#include "Chombo_LevelData.H"
#include "Chombo_EBChombo.H"
#include "Chombo_GeometryService.H"
#include "Chombo_NamespaceHeader.H"


///names that Proto::EBStencilArchive will recognize
namespace StencilNames
{
  ///second order Poisson stencil a la Schwartz, et al.
  static const string Poisson2("Second_Order_Poisson");

  ///Signals Dirichlet BCs
  static const string Dirichlet("Dirichlet");

  ///Signals Neumann  BCs
  static const string  Neumann("Neumann");
  
  ///placeholder for when there are no bcs for an operator
  static const string     NoBC("no_bcs");

  ///root for the string used for the volume weighted avergaging operator 
  /**
     Used in Trebotich and Graves
     This is for the non-conservative divergence.
     Add the radius to get the stencil for ncdiv = sum_rad (kappa div)/sum_rad(kappa)
  */
  static const string NCDivergeRoot("Volume_Weighted_Averaging_rad_"); 
  
  ///for correcting tangential velocity at faces
  static const string TanVelCorrect("Tangential_Velocity_Correction");

  ///copy values in cell fab to fab over cut faces
  /**
     In Trebotich and Graves, when evaluating the conservative gradient,
     the value of the pressure is required at the cut face.  One  way this is done
     is to just use the pressure at the cell center for the cut face pressure. 
     This stencil does this.   
   */
  static const string CopyCellValueToCutFace("copy_cell_value_to_cut_face");

  ///extrapolate face values to domain faces
  /**
     extrapolates face values to domain boundaries
     domval = 2*nearval - farval
   */
  static const string ExtrapToDomainFace("extrap_to_domain_face");


  /// Stencil for addng the contribution of the cut face flux F to (kappa div F)
  /**
     stencil weight = -alpha*normal[dir]/dx 
     where alpha is the boundary fraction of the cut face and dx is the grid spacing.
   */
  static const string CutFaceIncrementToKappaDiv("cut_face_increment_to_kappa_div_dir_");
  
  /// for simple average from faces to cells
  /**
    This is to get the velocity to faces in projection.
  */
  static const string AveCellToFace("AverageCellToFace"); 

  ///for extrap bcs used to get one-sided differences in slopes
  static const string LinearExtrapolation("LinearExtrapolation"); 

  /// for simple average from  cells to faces
  /**
    This is to get the cc  gradient in the projection from mac grad
  */
  static const string AveFaceToCell("AverageFaceToCell"); 

  ///for volume-weighted redistribution
  /**
     add a radius to get the stencil for volume-weighted redistribution (originally Chern and Colella but everywhere).
   */
  static const string SmushRoot("Volume_Weighted_Redistribution_rad_"); 

  ///for interpolating from face centers to face centroids
  static const string InterpToFaceCentroid("InterpolateToFaceCentroid"); 

  ///for (cell centered) one-sided differences on the low side
  /**
     For (cell centered) one-sided differences on the low side.  Append the string for 
     the direction (0, 1, or 2) to complete it
  */
  static const string SlopeLoRoot("Slope_Low_");   

  ///for (cell centered) one-sided differences on the high side
  /**
     For (cell centered) one-sided differences on the high side.  Append the string for 
     the direction (0, 1, or 2) to complete it
  */
  static const string SlopeHiRoot("Slope_High_");   

  ///for taking the divergence of face centered stuff to cell centered result
  static const string DivergeFtoC("Divergence");   

  ///for getting stuff from low  side cells to faces
  /**
     For getting stuff from low  side cells to faces.  Commonly used in Godunov-style 
     extrapolation to get the extrapolated states to face holders so the Riemann problem
     can be solved.
   */
  static const string CellToFaceLo("Cell_To_Face_Low"); 

  ///for getting stuff from high side cells to faces
  /**
     For getting stuff from high side cells to faces.  Commonly used in Godunov-style 
     extrapolation to get the extrapolated states to face holders so the Riemann problem
     can be solved.
   */
  static const string CellToFaceHi("Cell_To_Face_High"); 

  ///2nd order face gradient of cell based data
  /**
     This is the mac gradient
   */
  static const string MACGradient("MAC_Gradient"); 

}

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
               data_t a_dxPoint,   Point a_ghost);


  ///maybe we need an amr structure.
  EBDictionary(const shared_ptr<GeometryService<order> >   & a_geoserv,
               const vector<DisjointBoxLayout>             & a_grids,
               const vector<Box>                           & a_domain,
               const vector<data_t>                        & a_dx,
               Point  a_ghost);
  

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
    if(!found)
      MayDay::Error("dictionary did not find input domain");

    return retval;
  }
  
  ///
  /**
     diagonal weights are for relaxation--expensive to compute so leave them out if you do not need them
     dest grow is how many cells to grow valid box for destination.
  */
  inline void registerStencil(string a_stencilName,
                              string a_domainBCName[2*DIM],
                              string a_ebbcName,
                              Box    a_srcDomain,
                              Box    a_dstDomain,
                              bool needDiagonalWeights = false,
                              Point  a_destGrow = Point::Zeros());


  inline void registerStencil(string a_stencilName,
                              string a_domainBCName,
                              string a_ebbcName,
                              Box    a_srcDomain,
                              Box    a_dstDomain,
                              bool   a_needDiagonalWeights = false,
                              Point  a_destGrow = Point::Zeros())
  {
    string domBCNames[2*DIM];
    for(int ivec = 0; ivec < 2*DIM; ivec++)
    {
      domBCNames[ivec] = a_domainBCName;
    }
    registerStencil(a_stencilName,
                    domBCNames,
                    a_ebbcName,
                    a_srcDomain,
                    a_dstDomain,
                    a_needDiagonalWeights,
                    a_destGrow);
  }


  
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
