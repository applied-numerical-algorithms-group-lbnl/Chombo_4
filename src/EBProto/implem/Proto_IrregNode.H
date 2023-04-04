
#ifndef _Proto_IRREGNODE_H_
#define _Proto_IRREGNODE_H_


#include "base/Proto_Point.H"
#include "Proto_RealVect.H"
#include "base/Proto_LoHiSide.H"
#include "Proto_IndexedMoments.H"
#include "EBProto.H"
#include "Proto_IrregGraph.H"
using std::vector;
namespace Proto
{
///
/**
   Node for construction of geometric information.
*/
  template <int order>
  class IrregNode
  {
  public:
    typedef IndexedMoments<DIM  , order> IndMomDIM;
    typedef IndexedMoments<DIM-1, order> IndMomSDMinOne;

    ///
    IrregNode()
    {;    }

    ///
    ~IrregNode()
    {;    }


    ///
    Point     m_cell;
    int       m_cellIndex;

    ///
    /**
       Indicies into a_nodes to show connectivity.
       If the arc is to  an irregular
       cell, the index is the unique index of the vof
       in the cell.
       For arcs to regular cells, the arc = -2
       If the arc is to the domain boundary, arc  = -1.
    */
    vector<int>      m_arc[2*DIM];


    ///for covered cells
    inline void setMomentsToZero()
    {
      m_volumeMoments.setToZero();
      m_EBMoments.setToZero();
      for(int idir = 0; idir < DIM; idir++)
      {
        m_normalPartialDeriv[idir].setToZero();
        m_EBNormalMomentsDir[idir].setToZero();
      }

      for(int iarc = 0; iarc < 2*DIM; iarc++)
      {
        m_faceMoments[iarc].setToZero();
      }
    }

    //for regular cells
    inline void setMomentsToRegular(const double& a_dx)
    {
      m_volumeMoments.setToRegular(a_dx);
      m_EBMoments.setToZero();
      for(int idir = 0; idir < DIM; idir++)
      {
        m_normalPartialDeriv[idir].setToZero();
        m_EBNormalMomentsDir[idir].setToZero();
      }

      for(int iarc = 0; iarc < 2*DIM; iarc++)
      {
        m_faceMoments[iarc].setToRegular(a_dx);
      }
    }

    ///data for irregular nodes

    /// the moment at  the VoF associated with the  monomial with the input exponents
    /**
       Given VoF variables x, y, z,   p = mono(0), q = mono(1), r = mono(2),
       returns integral_over_VoF(x^p y^q z^r dV)
    **/
    IndMomDIM m_volumeMoments;


    void makeRegular(const Point& a_pt, const double& a_dx)
    {
      m_cell = a_pt;
      m_cellIndex = 0;
      setMomentsToRegular(a_dx);
      for(int iarc = 0; iarc < 2*DIM; iarc++)
      {
        m_arc[iarc] = vector<int>(1,0);
      }
    }
    
    void print()
    {
      using std::cout;
      using std::endl;
      cout << "vol moments:" << endl;
      m_EBMoments.print();
      cout << "eb moments:" << endl;
      m_volumeMoments.print();
      for(int idir = 0; idir < 2*DIM; idir++)
      {
        cout << "face moments for index :" << idir << endl;
        m_faceMoments[idir].print();
      }
    }
    /// the moment at  the irregular face associated with the  monomial with the input exponents 
    /**
       Given VoF variables x, y, z,   p = mono(0), q = mono(1), r = mono(2),
       returns integral_over_irregular_area((x^p y^q z^r) dA)
    **/
    IndMomDIM m_EBMoments;
  
    IndMomDIM m_normalPartialDeriv[DIM];
    IndMomDIM m_EBNormalMomentsDir[DIM];
    ///face centered moments
    /// the moment at  the face associated with the  monomial with the input exponents 
    /**
       Given face variables x, y,    p = mono(0), q = mono(1)
       returns integral_over_face_area((x^p y^q ) dA)
    **/
    IndMomSDMinOne m_faceMoments[2*DIM];


  };



  template <int order>
  inline void getIrregGraphFromNodes(vector<IrregGraph>              & a_graph,
                                     const vector<IrregNode<order> > & a_nodes)
  {
    a_graph.resize(a_nodes.size());
    for(int ivec = 0; ivec < a_nodes.size();ivec++)
    {
      a_graph[ivec].m_cell      = a_nodes[ivec].m_cell;
      a_graph[ivec].m_cellIndex = a_nodes[ivec].m_cellIndex;
      for(int iarc = 0; iarc < 2*DIM; iarc++)
      {
        a_graph[ivec].m_arc[iarc] = a_nodes[ivec].m_arc[iarc];
      }
    }
  }
}


#endif
