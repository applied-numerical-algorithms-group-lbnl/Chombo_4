#ifndef _MOMENTS_H_
#define _MOMENTS_H_
using std::vector;

#include "Proto_RealVect.H"

namespace Proto
{
/////////////////////////////////class edgeMo///////////////////////////////////
///
/**
   This class stores data about an edge of a control volume
*/
  class edgeMo
  {
  public:
    edgeMo()
    {
    }

    void define (const RealVect& a_loPt,
                 const RealVect& a_hiPt,
                 const bool& a_intersectLo,
                 const int& a_direction,
                 const bool& a_covered,
                 const bool& a_regular,
                 const bool& a_dontKnow);

    ///access member data
    RealVect getLo() const;
    RealVect getHi() const;
    bool getIntersectLo() const;

    double     getEdgeLength() const;
    RealVect getEdgeCentroid() const;

    bool isCovered() const;
    bool isRegular() const;
    bool dontKnow()  const;

    int direction() const;

    ///integrate monomial
    double moment(const Point& a_exponent) const;

  private:
    RealVect m_Lo;
    RealVect m_Hi;

    bool m_intersectLo;

    double     m_edgeLength;
    RealVect m_edgeCentroid;

    bool m_covered;
    bool m_regular;
    bool m_dontKnow;

    int m_varOfInt;
  };

/////////////////////////////////class faceMo///////////////////////////////////
  class faceMo
  {
  public:
    faceMo()
    {
    }

    void define(const edgeMo a_edges[4],
                const int& a_faceNormal,
                const bool& a_covered,
                const bool& a_regular,
                const bool& a_dontKnow);

    void setFaceArea(const double& a_area);
    double getFaceArea() const;

    void     setFaceCentroid(const RealVect& a_centroid);
    RealVect getFaceCentroid() const;

    bool isRegular() const;
    bool isCovered() const;
    bool dontKnow()  const;

    double   getBdLength() const;
    void   getNormal(double a_normal[2]) const;
    void   getEdges(edgeMo a_Edges[4]) const;
    edgeMo retrieveEdge(int& a_iEdge) const;

    int  getFaceNormal()const;
    void makeNormal();
    void normalize(double a_normalVec[2]);

  private:
    edgeMo m_edges[4]; //lexigraphic order:0Lo,0Hi,1Lo,1Hi,2Ho,2Hi

    double m_normalVec[2];

    bool m_covered;
    bool m_regular;
    bool m_dontKnow;

    RealVect m_centroid;
    double     m_areaFrac;
    int      m_faceNormal;
    double     m_bdLength;

    // disallowed for all the usual reasons
    void operator=(const faceMo& a_input);

    // disallowed for all the usual reasons
    faceMo(const faceMo& a_input);
  };

/////////////////////////////////class vofMo////////////////////////////////////
  class vofMo
  {
  public:
    vofMo()
    {
    }

    void define(const faceMo a_faces[6]);

    void makeNormal();
    void setNormal(double a_normalVec[3]);
    void getNormal(double a_normalVec[3]) const;
    void normalize(double a_normalVec[3]);

    void getFaces(faceMo a_faces[6]) const;

    double getBdArea() const;

  private:
    faceMo m_faces[6];

    double m_normalVec[3];
    double m_bdArea;
  };

//////////////////////classMoments//////////////////////////////////////
  class Moments
  {
  public:
    Moments();

    vector<double> momentCalc3D(const int& a_order,
                                vofMo&     a_vof);

    vector<double> momentCalc2D(const int&    a_order,
                                const faceMo& a_face);

    void  listOfMoments(const int&       a_order,
                        vector<Point>& a_exponents);

    void makeMatrix(const vector<Point>& a_list,
                    const vector<Point>& a_listPlus,
                    const vector<double>&    a_normalVec,
                    double**                 a_A,
                    const int&             a_faceNormal=999
      ); //use default faceNormal for 3D matrix or when SpaceDim==2
  };
}
#include "Proto_MomentsImplem.H"
#endif
