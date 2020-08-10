#ifndef __proto_cutcellmoments_H__
#define __proto_cutcellmoments_H__
#include <map>
using std::map;


#include "Proto_IndexTM.H"
#include "Proto_IFData.H"
#include "Proto_LSProblem.H"
#include "Proto_Notation.H"

using std::map;

namespace Proto
{
  template <int dim> 
  class CutCellMoments
  {
  public:
    typedef IndexTM<int,dim>  IvDim;
    typedef IndexTM<double,dim> RvDim;

    typedef map<IvDim,double >  PthMoment;

    typedef map<IndexTM<int,dim-1>,double > PthMomentLesserDimension;

    typedef map<IndexTM<int,1>,double >     OneDMoments;
    typedef map<int,IvDim>                LocPthMoment;
    typedef map<IvDim,int > PthMomentLoc;

    typedef map<Iv2,CutCellMoments<dim-1> > BdCutCellMoments;

    // Constructors
    inline
    CutCellMoments();
    inline
    CutCellMoments(const CutCellMoments<dim>& a_cutCellMoments);

    // This is used to build bd CutCellMoments
    inline
    CutCellMoments(const IFData<dim>& a_info);

    // Destructor
    inline
    ~CutCellMoments();

    // Returns a boundary element of one dim less 
    inline const CutCellMoments<dim - 1> getBdCutCellMoments(const Iv2& a_bdId) const;

    inline double changeMomentCoordinates(PthMoment               & a_refinedMomentMap,
                                          const IndexTM<int,dim>  & a_monomial,
                                          const IndexTM<double,dim> & a_refinedCenterDelta);

    inline void changeMomentCoordinatesToCellCenter();

    inline void changeMomentCoordinatesToParentCenter();

    inline void initialize(CutCellMoments<dim> & a_refinedCutCell);

    inline void initializeMap(PthMoment & a_map1,
                              PthMoment & a_map2);

    inline void initializeMap(PthMomentLesserDimension & a_map1,
                              PthMomentLesserDimension & a_map2);


    // et the value of the moment depending on the values of the booleans
    // m_allVerticesIn/Out these functions are used under refinement, when
    // moments on the faces are needed even on a covered/regular cell
    inline double getBdMoment(const IvDim             & a_mono,
                              const IFData<dim+1>     & a_IFData,
                              const IndexTM<double,dim> & a_refinedCenterDelta,
                              PthMoment                 a_fullCellMap = PthMoment());

    inline double getBdEBMoment(const IvDim             & a_mono,
                                const IFData<dim+1>     & a_IFData,
                                const IndexTM<double,dim> & a_refinedCenterDelta);

    inline void addBdMoments(CutCellMoments<dim>     & a_coarseCutCell,
                             const IFData<dim+1>     & a_IFData,
                             const int               & a_degreePmax,
                             const bool              & a_useConstraints,
                             const IndexTM<double,dim> & a_refinedCenterDelta,
                             const IndexTM<int,dim>  & a_localHilo);

    // Integrates a monomial over a full cell
    inline double fullCellQuadrature(const IndexTM<int,dim>      & a_mono,
                                     const CoordinateSystem<dim> & a_coord);

    // Output methods that check a_mono is in the map
    inline double getMoment(const IvDim   & a_mono,
                            const EBorVol & a_EBorVOL) const;

    // Methods for reading geom data that do sanity checks on results
    inline double getVol(const EBorVol & a_EBorVol) const;

    inline RvDim getCentroid(const EBorVol & a_EBorVOL) const;

    inline double getResidual(const int & a_iDegree,
                              const int & a_normJ) const;

    inline void setResidual(const double& a_value,
                            const int & a_iDegree,
                            const int & a_normJ);

    inline std::vector<double> sliceResidual(const int & a_iDegree) const;

    inline bool isCovered() const;

    inline bool isRegular() const;

    // Operators
    inline void operator=(const CutCellMoments<dim>& a_cutCellMoments);

    //volume moments
    PthMoment             m_moments;

    //eb moments
    PthMoment             m_EBmoments;

    //lower dimensional cut cells
    BdCutCellMoments      m_bdCutCellMoments;

    //edge interesections,normals and derivatives of normals
    IFData<dim>           m_IFData;

    //indicates that a boundary CutCellMoment coincides with the interface
    bool                  m_bdCCOn;

    //residual from the least squares problem at the highest dimension
    std::vector<std::vector<double> > m_residual;

    //number of active constraints
    int                   m_numActiveBounds;

    // records whether this cutCellMoment or any elements of BdCutCellMoments are using the zero vector for a normal
    bool                  m_badNormal;
  };

// One dimensional cutCellMoments
  template <> class CutCellMoments<1>
  {
  public:
    typedef map<IndexTM<int,1>,double> OneDMoments;

    // Constructors
    inline
    CutCellMoments();
    inline
    CutCellMoments(const CutCellMoments<1> & a_cutCellMoments);

    inline
    CutCellMoments(const IFData<1>& a_info);

    // Destructor
    inline
    ~CutCellMoments();

    inline double changeMomentCoordinates(OneDMoments           & a_refinedMap,
                                          const IndexTM<int,1>  & a_monomial,
                                          const IndexTM<double,1> & a_refinedCenterDelta);

    inline void changeMomentCoordinatesToCellCenter();

    inline void changeMomentCoordinatesToParentCenter();

    inline void initialize(CutCellMoments<1> & a_refinedCutCell);

    inline void initializeMap(OneDMoments & a_map1,
                              OneDMoments & a_map2);

    inline double getBdMoment(const IndexTM<int,1>  & a_mono,
                              const IFData<2>       & a_IFData,
                              const IndexTM<double,1> & a_refinedCenterDelta,
                              OneDMoments             a_fullCellMap = OneDMoments());

    inline double getBdEBMoment(const IndexTM<int,1>  & a_mono,
                                const IFData<2>       & a_IFData,
                                const IndexTM<double,1> & a_refinedCenterDelta);

    inline void addBdMoments(CutCellMoments<1>     & a_coarseCutCell,
                             const IFData<2>       & a_IFData,
                             const int             & a_degreePmax,
                             const bool            & a_useConstraints,
                             const IndexTM<double,1> & a_refinedCenterDelta,
                             const IndexTM<int,1>  & a_localHilo);

    // Output method that check a_mono is in the map
    inline double getMoment(const IndexTM<int,1> & a_mono,
                            const EBorVol        & a_EBorVOL) const;

    inline double getMoment(const IndexTM<int,1>& a_mono) const;

    // Methods for reading geom data that do sanity checks on results
    inline double getVol(const EBorVol& a_EBorVol) const;

    inline IndexTM<double,1> getCentroid(const EBorVol& a_EBorVol) const;

    inline bool isCovered() const;

    inline bool isRegular() const;

    // Operators
    inline void operator=(const CutCellMoments<1>& a_cutCellMoments);

    // Member data
    OneDMoments  m_moments;
    IFData<1>    m_IFData;
    bool         m_bdCCOn;
    int          m_numActiveBounds;
    bool         m_badNormal;

    // This is a dummy quantity: no EB in 1D
    OneDMoments  m_EBmoments;
  };

}

#include "Proto_CutCellMomentsImplem.H"

#endif
