
#ifndef _Proto_LSPROBLEM_H_
#define _Proto_LSPROBLEM_H_

#include <map>

#include "Proto_LSquares.H"
#include "Proto_IFData.H"
#include "Proto_IndexTM.H"

using std::map;
using std::vector;

namespace Proto
{
  template <int dim> class CutCellMoments;

  template <int dim> class LSProblem
  {
    typedef IndexTM<int,dim>IvDim;
    typedef map<IvDim,int > PthMomentLoc;
    typedef map<int,IvDim> LocPthMoment;

  public:

    //destructor will free m_matrix and free matrixes in m_rhs
    inline
    ~LSProblem();

    //constructor for just accessing fill map function
    // useConstraints ignored
    inline
    LSProblem(const int  & a_degreeP,
              const bool & a_useConstraints);

    //constructor
    inline
    LSProblem(const int                & a_orderAccuracy,
              const int                & a_degreeP,
              const bool               & a_useConstraints,
              const IndexTM<double,dim>  & a_normal);

    //Solves over determined system matrix*unknown = rhs
    inline
    int invertNormalEq(const vector<double> & a_rhs,
                       vector<double>       & a_residual);


    //accessor functions
    inline
    int getDegree()
    {
      return m_degreeP;
    }

    inline
    int getOrderAccuracy()
    {
      return m_order;
    }

    inline void getMatrix(double ** a_matrix);

    inline void getRhs(vector<double>& a_rhs);

    //getUnknowns will return the geometry info at the current level
    inline void getUnknowns(vector<double>& a_unknown);

    // Compute max and min values of monomial over the cell in local coordinates
    inline void monoMaxMin(double                   & a_maxVal,
                           double                   & a_minVal,
                           const IndexTM<int,dim> & a_mono,
                           const IFData<dim>      & a_IFData);

    inline void computeBounds(const IndexTM<double,dim>   & a_dx,
                              const CutCellMoments<dim> & a_ccm);


    inline const LocPthMoment& getLocMonomialMapDegreeP() const
    {
      return (const LocPthMoment&)m_locMonoP;
    }

    inline const PthMomentLoc& getMonomialLocMapDegreeP() const
    {
      return (const PthMomentLoc&)m_monoLocP;
    }

    inline const PthMomentLoc& getMonomialLocMapDegreePLess1() const
    {
      return (const PthMomentLoc&)m_monoLocPLess1;
    }

    inline double getUnknown(int loc)
    {
      return m_unknowns[loc];
    }

    inline int getNumberDegP()
    {
      return m_numP;
    }

    inline int getNumberDegPLess1()
    {
      return m_numPLess1;
    }

    inline int numActiveBounds() const
    {
      return m_numActiveBounds;
    }




    inline void setRhs(const vector<double>& a_rhs); //deprecated

    //constructs matrix of overdetermined system
    inline void setMatrix();

    inline void momentBounds(double              & a_lobnd,
                             double              & a_hibnd,
                             const IvDim       & a_mono,
                             const IFData<dim> & a_ifData);

    //\partial/\partial xidir (mono) = coeff*Dmono
    inline void differentiate(int        & a_coeff,
                              IvDim      & a_Dmono,
                              int        & a_idir,
                              const IvDim& a_mono);

    inline int nChooseR(int a_n,
                        int a_r);

    //makes a list of monomials of degreeP
   inline  void fillMap(PthMomentLoc & a_monoLoc,
                        LocPthMoment & a_locMono,
                        const int    & a_degree);

    //computes number of monomials of degree a_monoDegree
    inline int numMonomials(const int & a_monoDegree);

    inline int factorial(const int & a_n,
                  const int & a_m=0);

    //array tools
    inline void allocArray(const int & a_rows,
                    const int & a_cols,
                    double**    & a_A);

    inline
    void freeArray(const int & a_rows,
                   const int & a_cols,
                   double**    & a_A);

    //order of accuracy
    int m_order;

    //degree of monomial
    int m_degreeP;

    //number of active constraints
    int m_numActiveBounds;

    //use constrained LS
    bool m_useConstraints;

    //normal
    IndexTM<double,dim> m_normal;

    //monomials of degree P
    PthMomentLoc m_monoLocP;
    LocPthMoment m_locMonoP;

    //number of monomials of degreeP
    int m_numP;

    //monomials of degree P-1
    PthMomentLoc m_monoLocPLess1;
    LocPthMoment m_locMonoPLess1;

    //number of monomials of degreeP
    int m_numPLess1;


    //matrix
    double**  m_matrix;
    vector<double> m_unknowns;
    vector<double> m_rhs;

    vector<double> m_lowerBound;
    vector<double> m_upperBound;
  };

  template<>class LSProblem<1>
  {
  public:

    //number of monomials of degree P
    int m_numP;

    //number of monomials of degree P-1
    int m_numPLess1;

    //monomials of degree P
    double m_pMoments;

    //used in filling the map
    int m_ithMoment;

    //Destructor
    inline
    ~LSProblem();

    //copy constructor
    inline
    LSProblem(const LSProblem<1> & a_lsProblem);

    //empty Constructor
    inline
    LSProblem();

    //count num monomials
    inline int recursiveCount(const int & a_degreeP);

    inline void setNumMonomials();

    //equals operator
    inline void operator=(const LSProblem & a_lSProblem);
  };
}

#include "Proto_LSProblemImplem.H"

#endif
