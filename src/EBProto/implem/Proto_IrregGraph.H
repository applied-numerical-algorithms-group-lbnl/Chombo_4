
#ifndef _Proto_IRREGGRAPH_H_
#define _Proto_IRREGGRAPH_H_


#include "base/Proto_Point.H"

using std::vector;
namespace Proto
{


  ///simple struct to allow me to keep EBGraph untemplated
  struct IrregGraph
  {
  public:

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
    
    size_t linearSize( ) const
    {
      size_t retval = sizeof(Point); //m_cell
      retval += sizeof(int); //m_cellindex
      for(int iside = 0; iside < 2*DIM; iside++)
      {
        retval += sizeof(unsigned int); //vector size;
        unsigned int vecsize = m_arc[iside].size();
        retval += vecsize*sizeof(int);
      }
      return retval;
    }

    void linearIn(const char* a_buf)
    {
      char* charbuf = (char *) a_buf;
      Point* pointbuf = (Point *) charbuf;
      m_cell = *pointbuf;
      charbuf += sizeof(Point);
      
      int* intbuf = (int *) charbuf;
      m_cellIndex = *intbuf;
      charbuf += sizeof(int);

      for(int iside = 0; iside < 2*DIM; iside++)
      {
        unsigned int* sizeptr = (unsigned int *) charbuf;
        unsigned int  vecsize = *sizeptr;
        charbuf += sizeof(unsigned int);

        m_arc[iside].resize(vecsize);
        for(int ivec = 0; ivec < vecsize; ivec++)
        {
          int* arcbuf = (int *)charbuf;
          m_arc[iside][ivec] = *arcbuf;
          charbuf += sizeof(int);
        }
      }

    }

    void linearOut(char* a_buf) const
    {
      char* charbuf = a_buf;
      Point* pointbuf = (Point *) charbuf;
      *pointbuf = m_cell;
      charbuf += sizeof(Point);
      
      int* intbuf = (int *) charbuf;
      *intbuf = m_cellIndex;
      charbuf += sizeof(int);

      for(int iside = 0; iside < 2*DIM; iside++)
      {
        unsigned int* sizeptr = (unsigned int *)charbuf;
        *sizeptr = m_arc[iside].size();
        charbuf += sizeof(unsigned int);

        for(int ivec = 0; ivec < m_arc[iside].size(); ivec++)
        {
          int* arcbuf = (int *)charbuf;
          *arcbuf = m_arc[iside][ivec];
          charbuf += sizeof(int);
        }
      }
    }
  };

}


#endif
