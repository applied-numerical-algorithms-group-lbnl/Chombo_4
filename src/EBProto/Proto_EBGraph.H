#ifndef __PROTO_EBGraph__
#define __PROTO_EBGraph__

#include "implem/Proto_IrregGraph.H"
#include "EBProto.H"
#include "Proto_RealVect.H"
using std::vector;
using std::set;

namespace Proto
{

  ///
  /// There is an important lesson to be learned here about template usage. On the surface this class makes a lot of sense as a template, it does something fairly different depending upon the centering. The problem with this becomes, now anything that uses this must also be dependent on the centering also. As a result, there is no sane way to loop over each direction or do things in a direction independent way, leading to lots of code reuse
  /// REVIEW Consider changing this to store the centering as a trait of the index, rather than as a template
  template<CENTERING cent>
  class EBIndex
  {
  public: //private: //FIXME this is bad design, especially for such a simple class

    Point m_pt;
    //for faces, the high cell. for cells, this is the cell
    int   m_vofIDMe; 
    //for faces, the low cell. for cells, undefined
    int   m_vofIDLo;

    ///
    bool m_isBoundary;

  public:
    EBIndex()
      {
        m_isBoundary = false;
        m_vofIDLo = 0;
        m_vofIDMe = 0;
      }

    ~EBIndex()
      {
      }

    ///for vof construction
    EBIndex(Point a_pt, int a_cellInd) 
    {
      m_pt = a_pt;
      m_vofIDMe = a_cellInd;
      m_vofIDLo = 0;
      m_isBoundary = false;
    }

    // not great, but at least we can stop calling the member data directly
    Point& getPoint()
      {
        return m_pt;
      }

    bool& isBoundary()
      {
        return m_isBoundary;
      }

    size_t linearSize( ) const
    {
      size_t retval = 0;
      retval += (3+DIM)*sizeof(int);
      return retval;
    }

    ///
    static vector< EBIndex<cent> > getRefinedRegular(const EBIndex<cent>& a_pt,
                                                     const Box&           a_fineDomain) ;

    ///for serialization (chombo uses this in parallel code)
    void linearIn(const char* a_buf)
    {
      int* intbuf = (int *)a_buf;
      for(int idir = 0; idir < DIM; idir++)
      {
        m_pt[idir] = *intbuf;
        intbuf++;
      }
      m_vofIDMe = *intbuf;
      intbuf++;
      m_vofIDLo = *intbuf;
      intbuf++;

      int ibound = *intbuf;
      m_isBoundary = (ibound != 0);
      intbuf++;
    }


    ///for serialization (chombo uses this in parallel code)
    void linearOut(char* a_buf) const
    {
      int* intbuf = (int *) a_buf;
      for(int idir = 0; idir < DIM; idir++)
      {
        *intbuf = m_pt[idir];
        intbuf++;
      }
      *intbuf = m_vofIDMe;
      intbuf++;
      *intbuf = m_vofIDLo;
      intbuf++;

      int ibound = 0;
      if(m_isBoundary) ibound = 1;

      *intbuf = ibound;
      intbuf++;
    }

    void print() const
    {
#if DIM==2
      printf("((%d, %d), %d, %d)", m_pt[0], m_pt[1], m_vofIDMe, m_vofIDLo);
#endif
#if DIM==3
      printf("((%d, %d, %d), %d, %d)", m_pt[0], m_pt[1], m_pt[2], m_vofIDMe, m_vofIDLo);
#endif
    }

    CUDA_DECORATION
    bool operator< (const EBIndex<cent>& a_ind) const
    {
      bool retval = false;
      if(m_pt == a_ind.m_pt)
      {
        retval = (m_vofIDMe < a_ind.m_vofIDMe);
      }
      else
      {
        retval = (m_pt < a_ind.m_pt);
      }
      return retval;
    }

    CUDA_DECORATION
    bool operator==(const EBIndex<cent>& a_ind) const
    {
  
      return ((m_pt == a_ind.m_pt) && (m_vofIDMe == a_ind.m_vofIDMe));
    }


    bool operator!=(const EBIndex<cent>& a_ind) const
    {
      return !(*this ==  a_ind);
    }

    //for sussing out a face's direction
    inline int direction() const
    {
      int retval = -1;
      if(cent == XFACE)
      {
        retval = 0;
      }
      else if(cent == YFACE)
      {
        retval = 1;
      }
      else if(cent == ZFACE)
      {
        retval = 2;
      }
      return retval;
    }
    //for getting the vof of a cut face 
    inline EBIndex<CELL> getVoF() const
    {
      EBIndex<CELL> retval;
      retval.m_pt = m_pt;
      retval.m_vofIDMe = m_vofIDMe;
      retval.m_vofIDLo = 0;
      retval.m_isBoundary = false;
      return retval;
    }


    //for getting the  boundary face  of a vof
    inline EBIndex<BOUNDARY> getCutFace() const
    {
      EBIndex<BOUNDARY> retval;
      retval.m_pt = m_pt;
      retval.m_vofIDMe = m_vofIDMe;
      retval.m_vofIDLo = 0;
      retval.m_isBoundary = false;
      return retval;
    }

    //for getting the vofs of a coordinate aligned face
    inline EBIndex<CELL> getVoF(Side::LoHiSide a_sd) const
    {
      EBIndex<CELL> retval;
      if(a_sd == Side::Lo)
      {
        int dir = direction();
        retval.m_pt = m_pt - Point::Basis(dir);
        retval.m_vofIDMe = m_vofIDLo;
        retval.m_vofIDLo = 0;
      }
      else
      {
        retval.m_pt = m_pt;
        retval.m_vofIDMe = m_vofIDMe;
        retval.m_vofIDLo = 0;
      }
      retval.m_isBoundary = 0;
      return retval;
    }
    

    void defineAsRegular(const Point& a_pt)
    {
      m_pt = a_pt;
      m_vofIDMe = 0;
      m_vofIDLo = 0;
      m_isBoundary = 0;
    }
    //for making a face
    EBIndex(const EBIndex<CELL>& a_vof1,const EBIndex<CELL>& a_vof2)
    {
      if(a_vof1.m_pt < a_vof2.m_pt)
      {
        m_pt      = a_vof2.m_pt;
        m_vofIDMe = a_vof2.m_vofIDMe;
        m_vofIDLo = a_vof1.m_vofIDMe;
      }
      else
      {
        m_pt      = a_vof1.m_pt;
        m_vofIDMe = a_vof1.m_vofIDMe;
        m_vofIDLo = a_vof2.m_vofIDMe;
      }
      m_isBoundary = ((m_vofIDMe < 0) || (m_vofIDLo < 0));
    }

  };

  template <CENTERING dstCent, CENTERING srcCent>
  EBIndex<dstCent> getEBIndex(const EBIndex<srcCent>& a_src)
    {
      EBIndex<dstCent> retval;
      retval.m_pt         = a_src.m_pt;
      retval.m_vofIDMe    = a_src.m_vofIDMe;
      retval.m_vofIDLo    = a_src.m_vofIDLo;
      retval.m_isBoundary = a_src.m_isBoundary;
      return retval;
    }


  template <CENTERING cent>
  inline RealVect indexLoc(const EBIndex<cent>& a_index,const RealVect& a_origin, double a_dx) 
  {
    //default implementation is  for faces
    int  facedir = a_index.direction();
    RealVect retval;
    for(int idir = 0; idir < DIM; idir++)
    {
      retval[idir] = a_origin[idir] + a_dx*(0.5 + double(a_index.m_pt[idir]));
    }
    retval[facedir] = a_origin[facedir] + a_dx*(double(a_index.m_pt[facedir]));
    return retval;
  }

  template < >
  inline RealVect indexLoc(const EBIndex<CELL>& a_index,const RealVect& a_origin, double a_dx) 
  {
    //default implementation is  for faces
    RealVect retval;
    for(int idir = 0; idir < DIM; idir++)
    {
      retval[idir] = a_origin[idir] + a_dx*(0.5 + double(a_index.m_pt[idir]));
    }
    return retval;
  }

  template <>
  inline RealVect indexLoc(const EBIndex<BOUNDARY>& a_index,const RealVect& a_origin, double a_dx) 
  {
    //default implementation is  for faces
    RealVect retval;
    for(int idir = 0; idir < DIM; idir++)
    {
      retval[idir] = a_origin[idir] + a_dx*(0.5 + double(a_index.m_pt[idir]));
    }
    return retval;
  }

  
  ///
  /**
     This is a list showing the connectivity of a given cell.  This
     has the special property of also being able to set itself to regular
     or covered by setting the current nodeimplem pointer to 1 or 0.
  */
  class CellNode
  {
   public:

    //each cell  contains a pointer to vector of these
    //the pointer has a special value for regular and covered cells.
    struct IrrGraphNode
    {
      //cell index is implied by place in the vector<IrrGraphNode>
      //grid index is implied by place in the HostData
      ///
      vector<int> m_arc[2*DIM];
      int         m_coarser;
      vector<EBIndex<CELL> > m_finer;

      IrrGraphNode& operator=(const IrrGraphNode& a_ndin)
      {
        if(&a_ndin != this)
        {
          for(int iarc = 0; iarc < 2*DIM; iarc++)
          {
            m_arc[iarc] = a_ndin.m_arc[iarc];
          }
          m_coarser = a_ndin.m_coarser;
          m_finer   = a_ndin.m_finer;
        }
        return *this;
      }
      ///for serialization (chombo uses this in parallel code)
      size_t linearSize( ) const
      {
        size_t charsiz = 0;
        //all the arcs
        for(int iside = 0; iside < 2*DIM; iside++)
        {
          charsiz += sizeof(int); //vector size;
          int vecsize = m_arc[iside].size();
          charsiz += vecsize*sizeof(int);
        }
        //coarser
        charsiz += sizeof(int);
        //finer
        int finesize = m_finer.size();
        charsiz += sizeof(int);
        for(int ivec = 0; ivec < finesize; ivec++)
        {
          charsiz += m_finer[ivec].linearSize();
        }
        return charsiz;
      }

      ///for serialization (chombo uses this in parallel code)
      void linearIn(const char* a_buf)
      {
        char* charbuf  = (char *) a_buf;
        size_t charsiz = 0;
        int* intbuf = (int *) charbuf;
        for(int iside = 0; iside < 2*DIM; iside++)
        {
          intbuf = (int *) charbuf;
          int vecsize = *intbuf;
          charbuf += sizeof(int);
          charsiz += sizeof(int);
          
          m_arc[iside].resize(vecsize);
          for(int ivec = 0; ivec < vecsize; ivec++)
          {
            intbuf = (int * )charbuf;
            m_arc[iside][ivec] = *intbuf;
            charbuf += sizeof(int);
            charsiz += sizeof(int);
          }
        }

        //coarser
        intbuf = (int *) charbuf;
        m_coarser = *intbuf;
        charbuf += sizeof(int);
        charsiz += sizeof(int);
        
        //finer
        intbuf = (int *) charbuf;
        int vecsize = *intbuf;
        charbuf += sizeof(int);
        charsiz += sizeof(int);
        
        m_finer.resize(vecsize);
        for(int ivec = 0; ivec < vecsize; ivec++)
        {
          m_finer[ivec].linearIn(charbuf);
          charbuf += m_finer[ivec].linearSize();
          charsiz += m_finer[ivec].linearSize();
        }
      }


      ///for serialization (chombo uses this in parallel code)
      void linearOut(char* a_buf) const
      {
        char* charbuf  = (char *) a_buf;
        size_t charsiz = 0;
        int* intbuf = (int *) charbuf;
        for(int iside = 0; iside < 2*DIM; iside++)
        {
          intbuf = (int *) charbuf;
          int vecsize = m_arc[iside].size();
          *intbuf = vecsize;
          charbuf += sizeof(int);
          charsiz += sizeof(int);
          
          for(int ivec = 0; ivec < vecsize; ivec++)
          {
            intbuf = (int * )charbuf;
            *intbuf = m_arc[iside][ivec];
            charbuf += sizeof(int);
            charsiz += sizeof(int);
          }
        }

        //coarser
        intbuf = (int *) charbuf;
        *intbuf = m_coarser;
        charbuf += sizeof(int);
        charsiz += sizeof(int);
        
        //finer
        intbuf = (int *) charbuf;
        int vecsize = m_finer.size();
        *intbuf = vecsize;
        charbuf += sizeof(int);
        charsiz += sizeof(int);

        for(int ivec = 0; ivec < vecsize; ivec++)
        {
          m_finer[ivec].linearOut(charbuf);
          charbuf += m_finer[ivec].linearSize();
          charsiz += m_finer[ivec].linearSize();
        }
      }
    };

    template <CENTERING cent>
    inline void
    getFaceInfo(vector<EBIndex<cent> > & a_faces,
                const int            & a_dir,
                const  EBIndex<CELL> & a_vof,
                const Side::LoHiSide & a_sd,
                const Box            & a_domain) const;


    template <CENTERING cent>
    inline void
    setFaceStuff(EBIndex<cent>             &   a_face,
                 const int                 &   a_inputIndex,
                 const int                 &   a_otherIndex,
                 const Point               &   a_inputIV,
                 const Point               &   a_otherIV,
                 const Side::LoHiSide      &   a_sd,
                 const Box                 &   a_domain) const;
    
    CellNode()
    {
      m_cellList = 0;
    }

    ~CellNode()
    {
      if(isIrregular())
      {
        delete m_cellList;
      }
    }
    ///
    /**
       Return true if the node is covered (m_cellList==0).
    */
    inline bool isCovered() const
    {
      return m_cellList==0;
    }

    ///
    /**
       Return true if the node is regular.
    */
    inline bool isRegular() const
    {
      unsigned long long int cellList = (unsigned long long int)(m_cellList);
      return cellList==1;
    }

    ///
    /**
       Return true if the node is neither regular or covered.
    */
    inline bool isIrregular() const
    {
      return ((!isRegular()) && (!isCovered()));
    }

    ///
    /**
       Set Node to regular.  If previously set to irregular,
       deletes their memory.
    */
    inline void defineAsRegular()
    {
      if(isIrregular()) delete m_cellList;
      m_cellList = (vector<IrrGraphNode>*)(1);
    }


    //resets the size of the irregular list to 1 and sets all the arcs to 0
    //or -1 depending on whether they are boundary arcs
    inline void defineAsFullIrregular(const Point& a_pt, const Box& a_domain)
    {
      if(!isIrregular())
      {
        m_cellList = new vector<IrrGraphNode>();
      }
      m_cellList->resize(1);
      (*m_cellList)[0].m_coarser = 0;

      EBIndex<CELL> vof;
      vof.defineAsRegular(a_pt);
      (*m_cellList)[0].m_finer = EBIndex<CELL>::getRefinedRegular(vof, a_domain);
      for(int idir = 0; idir < DIM; idir++)
      {
        for(SideIterator sit; sit.ok(); ++sit)
        {
          int index = ebp_index(idir, sit());
          Point otherPt = a_pt + sign(sit())*Point::Basis(idir);
          int iarc = 0;
          if(!a_domain.contains(otherPt))
          {
            iarc = -1;
          }
          (*m_cellList)[0].m_arc[index].resize(1, iarc);
        }
      }
    }

        //resets the size of the irregular list to 0
    inline void defineAsIrregular()
    {
      if(!isIrregular())
      {
        m_cellList = new vector<IrrGraphNode>();
      }
      m_cellList->resize(0);
    }

    ///
    /**
       Set Node to covered.  If previously set to irregular,
       deletes their memory.
    */
    inline void defineAsCovered()
    {
      if(isIrregular()) delete m_cellList;
      m_cellList = (vector<IrrGraphNode>*)(0);
    }

    ///
    /**
       Get the faces in the direction and side
       for the vof in the list.  if the vof's cell
       index is not found in the list, abort.
       Use the input vof's grid index for a the grid index
       of the list.
    */
    inline vector<EBIndex<XFACE> >
    getXFaces(const  EBIndex<CELL> & a_vof,
              const Side::LoHiSide & a_sd,
              const Box            & a_domain) const;


    inline vector<EBIndex<YFACE> >
    getYFaces(const  EBIndex<CELL> & a_vof,
              const Side::LoHiSide & a_sd,
              const Box            & a_domain) const;


#if DIM==3
    inline vector<EBIndex<ZFACE> >
    getZFaces(const  EBIndex<CELL>  & a_vof,
              const Side::LoHiSide  & a_sd,
              const Box             & a_domain) const;

#endif


    ///
    /**
       Return all the vofs in the list, using the input
       intvect for the gridIndex
    */
    inline vector< EBIndex<CELL> >
    getVoFs(const Point & a_iv) const;


    ///return vofs labelled as irregular (only places where an EB can happen)
    inline vector<EBIndex<BOUNDARY> >
    getCutFaces(const Point & a_iv) const;



    ///
    /**
       The connectivity data at this point.
       If m_cellList == 0, node is covered.
       If m_cellList == 1, node is regular and parent is single-valued.
       otherwise, it is a real list (node may still be regular).
    */
    vector<IrrGraphNode>* m_cellList;


    CellNode(const CellNode& a_ebin)
    {
      *this = a_ebin;
    }

    CellNode& operator=(const CellNode& a_ebin)
    {
      if(&a_ebin != this)
      {

        if(a_ebin.isRegular())
        {
          defineAsRegular();
        }
        else if(a_ebin.isCovered())
        {
          defineAsCovered();
        }
        else
        {
          defineAsIrregular();
          *m_cellList = *a_ebin.m_cellList;
        }
      }
      return *this;
    }

    ///for serialization (chombo uses this in parallel code)
    static bool constantMemorySize()
    {
      return false;
    }

    ///for serialization (chombo uses this in parallel code)
    size_t linearSize( ) const
    {
      size_t retval = sizeof(unsigned int);  //0 or 1 says whether this is a real vector 
      if( isRegular() || isCovered())
      {
        retval += sizeof(unsigned int);
      }
      else
      {
        retval += sizeof(unsigned int); //size of vector
        for(int ivec = 0; ivec < m_cellList->size(); ivec++)
        {
          retval += (*m_cellList)[ivec].linearSize();
        }
      }
      return retval;
    }

    ///for serialization (chombo uses this in parallel code)
    inline void linearIn(const char* a_buf)
    {
      char* charbuf = (char *) a_buf;
      unsigned int* intbuf = (unsigned int*) charbuf;
      unsigned int  hasIrreg = *intbuf;
      charbuf += sizeof(unsigned int);

      if(hasIrreg == 0)
      {
        unsigned int* isregbuf  = (unsigned int*) charbuf;
        unsigned int  isregular = *isregbuf;
        charbuf += sizeof(unsigned int);
        if(isregular == 0)
        {
          defineAsCovered();
        }
        else
        {
          defineAsRegular();
        }
      }
      else
      {
        defineAsIrregular();
        unsigned int* intbuf = (unsigned int*) charbuf;
        unsigned int  vecsize = *intbuf;
        charbuf += sizeof(unsigned int);
        m_cellList->resize(vecsize);
        for(int ivec = 0; ivec < m_cellList->size(); ivec++)
        {
          (*m_cellList)[ivec].linearIn(charbuf);
          charbuf += (*m_cellList)[ivec].linearSize();
        }
      }
    }


    ///for serialization (chombo uses this in parallel code)
    inline void linearOut(char* a_buf) const
    {
      unsigned int hasIrreg = 0;
      if(isIrregular()) 
      {
        hasIrreg = 1;
      }
      char* charbuf = a_buf;
      { //just for scoping
        unsigned int* intbuf = (unsigned int*) charbuf;
        *intbuf = hasIrreg;
        charbuf += sizeof(unsigned int);
      }
      if(hasIrreg == 0)
      {
        unsigned int* intbuf = (unsigned int*) charbuf;
        unsigned int  isregular = 0;
        if(isRegular())
        {
          isregular = 1;
        }
        *intbuf = isregular;
        charbuf += sizeof(unsigned int);
      }
      else
      {
        unsigned int* intbuf = (unsigned int*) charbuf;
        unsigned int  vecsize = m_cellList->size();
        *intbuf  = vecsize;
        charbuf += sizeof(unsigned int);

        for(int ivec = 0; ivec < m_cellList->size(); ivec++)
        {
          (*m_cellList)[ivec].linearOut(charbuf);
          charbuf += (*m_cellList)[ivec].linearSize();
        }
      }
    }



  };


/// Geometric description within a box
/**
   EBGraph represents the geometric information
   of the domain at a given refinement
   within the boundaries of a particular box.
*/
  class EBGraph
  {
  public:
    ///
    /**
       The EBGraph can be all regular or all
       covered or contain irregular cells.
    */
    enum TAG
    {
      AllRegular = 0,
      AllCovered,
      HasIrregular
    };

    ///
    /**
       Makes an EBGraph whose contents are undefined.
    */
    EBGraph()
    {
      m_isDefined = false;
    }

    ///
    template<CENTERING cent>
    inline 
    bool isFaceRegular(const EBIndex<cent>& a_face) const;

    ///
    template<CENTERING cent>
    inline 
    bool isFaceCovered(const EBIndex<cent>& a_face) const;

    ///needed for Proto::LevelData.  meaningless here becasue we need the graph.
    EBGraph(const Box& a_box, bool a_stackAlloc)
    {
    }

    
    ///given vof on a coarser level, get vofs on THIS level it was coarsened from
    /**
       This is useful when you only have the graph on the finer level but need to refine the
       coarser vof.   Trust me, it comes up.
     */
    inline vector<EBIndex<CELL> >
    getVoFsThatCoarsenTo(const EBIndex<CELL>&   a_coarVoF) const
    {
      vector<EBIndex<CELL> > retval;

      Box coarptbox(a_coarVoF.m_pt, a_coarVoF.m_pt);
      Box finebox = coarptbox.refine(2);
      for(BoxIterator bit = finebox.begin(); bit != finebox.end(); ++bit)
      {
        vector<EBIndex<CELL> > finevofs = getVoFs(*bit);
        for(int ifine = 0; ifine < finevofs.size(); ifine++)
        {
          EBIndex<CELL>& finevof = finevofs[ifine];
          EBIndex<CELL>  cofiVoF = coarsen(finevof);
          if(cofiVoF == a_coarVoF)
          {
            retval.push_back(finevof);
          }
        }
      }

      return retval;
    }

    static int preAllocatable()
    {
      return 2; // dyanmic allocatable.
    }

    static int memTypeAllocation()
    {
      return 0;
    }

    inline  
    void copy(const EBGraph              & a_dsrc,
              const Box                  & a_srcBox,
              unsigned int                 a_srcComp,
              const Box                  & a_destBox,
              unsigned int                 a_destComp,
              unsigned int                 a_numcomp);
    ///
    virtual ~EBGraph()
    {
    }

    ///returns false if no such unique face exists (multivalued, outside of the domain, whatever)
    /**
       tandir should be shomething different from face.direction
     */
    template <CENTERING cent> 
    inline bool
    getAdjacentFace(EBIndex<cent>       & a_adjFace, 
                    const EBIndex<cent> & a_startFace,
                    int a_tanDir, Side::LoHiSide a_side) const;
    ///
    /**
       Return the VoFs in the cell.
    */
    inline vector< EBIndex<CELL> > getVoFs(const Point& a_iv) const;


    inline vector< EBIndex<CELL> >
    getAllVoFsWithinRadius(const Point & a_iv, unsigned int a_radius) const
    {
      vector<EBIndex<CELL> > retval;
      Box cellBox(a_iv, a_iv);
      Box grownBox = cellBox.grow(a_radius);
      grownBox &= m_domain;

      for(BoxIterator bit = grownBox.begin(); bit != grownBox.end(); ++bit)
      {
        vector< EBIndex<CELL> > cellvofs = getVoFs(*bit);
        retval.insert(retval.end(), cellvofs.begin(), cellvofs.end());        
      }
      return retval;
    }



    ///
    /**
       Return true if every cell in the EBGraph is a regular cell.
    */
    inline bool isAllRegular() const
    {
      return (m_tag == AllRegular);
    }

    ///
    /**
       Return true neither all regular nor all covered
    */
    inline bool hasIrregular() const
    {
      return (m_tag == HasIrregular);
    }


    ///
    /**
       Return true if  every cell in the EBGraph is a covered cell.
    */
    inline bool isAllCovered() const
    {
      return (m_tag == AllCovered);
    }
    ///
    /**
       Return true if a_iv is an irregular cell.
    */
    inline bool isIrregular(const Point& a_iv) const;

    ///
    /**
       Return true if every cell in a_box is a
       covered cell.
    */
    inline bool isCovered(const Box& a_box) const;

    ///
    /**
       Return true if every cell in a_box is a
       regular cell.
    */
    inline
    bool isRegular(const Box  & a_box) const;

    ///
    inline
    bool isRegular(const Point& a_iv) const;

    ///
    inline
    bool isCovered(const Point& a_iv) const;


    ///
    inline
    bool isConnected(const  EBIndex<CELL> & a_vof1,
                     const  EBIndex<CELL> & a_vof2) const;

    ///
    /**
       Return the faces on the side and direction
       of the input VoF.
    */
    inline
    vector< EBIndex<XFACE> > getXFaces(const  EBIndex<CELL>      & a_vof,
                                       const Side::LoHiSide      & a_sd) const;


    ///
    /**
       Return the faces on the side and direction
       of the input VoF.
    */
    inline
    vector< EBIndex<YFACE> > getYFaces(const  EBIndex<CELL>       & a_vof,
                                       const Side::LoHiSide      & a_sd) const;


#if DIM==3
    ///
    /**
       Return the faces on the side and direction
       of the input VoF.
    */
    inline
    vector< EBIndex<ZFACE> > getZFaces(const  EBIndex<CELL>    & a_vof,
                                       const Side::LoHiSide    & a_sd) const;

#endif

    inline vector<EBIndex<CELL> > getNeighborVoFs(const  EBIndex<CELL> & a_vof,
                                                  const  unsigned int  & a_dir,
                                                  const Side::LoHiSide & a_sd) const
    {
      vector< EBIndex<CELL> > retval;
      if(a_dir == 0)
      {
        vector< EBIndex<XFACE> > faces = getXFaces(a_vof, a_sd);
        for(int iface = 0; iface < faces.size(); iface++) 
        {
          if(!faces[iface].m_isBoundary)
          {
            retval.push_back(faces[iface].getVoF(a_sd));
          }
        }
      }
      else if(a_dir == 1)
      {
        vector< EBIndex<YFACE> >faces = getYFaces(a_vof, a_sd);
        for(int iface = 0; iface < faces.size(); iface++) 
        {
          if(!faces[iface].m_isBoundary)
          {
            retval.push_back(faces[iface].getVoF(a_sd));
          }
        }
      }
#if DIM==3
      else if(a_dir == 2)
      {
        vector< EBIndex<ZFACE> >faces = getZFaces(a_vof, a_sd);
        for(int iface = 0; iface < faces.size(); iface++) 
        {
          if(!faces[iface].m_isBoundary)
          {
            retval.push_back(faces[iface].getVoF(a_sd));
          }
        }
      }
#endif
      else
      {
        PROTO_ASSERT(false, "bogus idir");
      }
      return retval;
    }

    ///
    /**
       Completely construct the graph from IrregGraphs.
       This makes a DENSE representation of the graph.  To make
       sparse representations of all regular or all covered graphs,
       use setToAllRegular or setToAllCovered.\\
       irregGraph contains a complete list of irregular vofs. \\
       validRegion is the the region of the graph \\
       domain is the domain of compuation \\
       regIrregCovered = 1 for regular, 0 for irregular, -1 for covered.
       regIrregCovered must contain the valid region grown by 1. \\
       In chombospeak, \\
       regIrregCovered.box().contains(grow(validRegion, 1) & domain); \\
       otherwise an error will result.
       SingleValuedOnly signals that the geometry generator can only generate
       single-valued cut cells (GeometryService has this quality).  
       Coarsening can still create multi-valued cells.
    */
    inline
    void
    defineFromPrimitives(const Box               &     a_validBox,
                         const HostBoxData<int>  &     a_regIrregCovered,
                         const vector<IrregGraph>&     a_irregGraph,
                         const Box               &     a_validRegion,
                         const Box               &     a_domain,
                         bool singleValuedOnly = true);

    ///
    /**
       Set the graph to all regular cells (sparse);
    */
    inline void defineAsAllRegular(const Box& a_validBox, const Box& a_region, const Box& a_domain);


    ///
    /**
       Set the graph to all regular cells (sparse);
    */
    inline
    void defineAsAllCovered(const Box& a_validBox, const Box& a_region, const Box& a_domain);

    ///
    inline
    const Box& getRegion() const
    {
      return m_region;
    }
    
    inline
    const Box& validBox() const
    {
      return m_validBox;
    }


    ///
    inline
    const Box& box() const
    {
      return m_region;
    }

    ///
    inline const Box& getDomain() const
    {
      return m_domain;
    }

    ///
    inline bool isDefined() const
    {
      return m_isDefined;
    }

    ////get the cut cells
    inline PointSet getIrregCells(const Box& a_box) const
    {
      PR_TIME("EBGraph::getIrregCells");
      PR_assert(m_isDefined);

      PointSet irreg;
      if(hasIrregular())
      {
        irreg = m_internals->m_irregIVS;
        irreg &= a_box;
      }
      return irreg;
    }

    //get the cut vofs
    inline vector< EBIndex<CELL> > getIrregVoFs(const Box& a_box) const
    {
      PR_TIME("EBGraph::getIrregVoFs");
      Box dombox = a_box & m_domain;

      PointSet irreg = getIrregCells(dombox);
      vector< EBIndex<CELL> > vofs;
      for(PointSetIterator ivsit(irreg);  ivsit.ok(); ++ivsit)
      {
        vector<EBIndex<CELL> > vofpt = getVoFs(ivsit());
        vofs.insert(vofs.end(), vofpt.begin(), vofpt.end());
      }
      return vofs;
    }

    ///used in coarsening--get collections of connected vofs
    inline
    vector< set<EBIndex<CELL> > >
    getVoFSets(const Box& a_region) const;

    inline 
    bool
    hasBoundaryFaces(const EBIndex<CELL> & a_vof, int a_idir, Side::LoHiSide a_sd);

    ///
    inline vector< EBIndex<BOUNDARY> > getCutFaces(const Point& a_pt) const;

    ///get the cut faces (basically the same as above
    inline vector< EBIndex<BOUNDARY> > getCutFaces(const Box& a_box) const
    {
      Box dombox = a_box & m_domain;
      PointSet irreg = getIrregCells(dombox);
      vector< EBIndex<BOUNDARY> > vofs;
      for(PointSetIterator ivsit(irreg);  ivsit.ok(); ++ivsit)
      {
        vector<EBIndex<BOUNDARY> > vofpt = getCutFaces(ivsit());
        vofs.insert(vofs.end(), vofpt.begin(), vofpt.end());
      }
      return vofs;
    }


    inline vector<EBIndex<CELL> > getAllVoFs(const Box& a_box) const
    {
      vector<EBIndex<CELL> > allvofs;
      for (auto bit = a_box.begin(); bit != a_box.end(); ++bit)
      {
        vector<EBIndex<CELL> > cellvofs = getVoFs(*bit);
        allvofs.insert(allvofs.end(), cellvofs.begin(), cellvofs.end());
      }
      return allvofs;
    }

    ///for faces
    template <CENTERING cent>
    inline
    vector< EBIndex<cent> > refine(const EBIndex<cent> & a_ind,
                                   const EBGraph       & a_finer) const;


    template <CENTERING cent>
    inline EBIndex<cent>   coarsen(const EBIndex<cent> & a_cent) const;

    inline vector< EBIndex<CELL> > 
    refine(const EBIndex<CELL> & a_vof) const;
    ///
    inline
    void
    coarsenVoFs(const Box     & a_validBox,
                EBGraph       &  a_fineEBGraph,
                const Box     &  a_vofRegionCoar);
    ///
    inline
    void
    coarsenFaces(const EBGraph  &  a_fineEBGraph,
                 const Box      &  a_faceRegionCoar);

    inline
    vector<int>
    coarsenFaces(const EBGraph        & a_fineEBGraph, 
                 const Box            & a_regionCoar,
                 const EBIndex<CELL>  & a_vofCoar,
                 const int            & a_idir,
                 const Side::LoHiSide & a_sd);

    inline bool
    hasBoundaryFaces(const EBIndex<CELL> & a_vof, int a_idir, Side::LoHiSide a_sd) const;

    ///needed for leveldata
    EBGraph(const Box& a_grid)
    {
      m_isDefined = false;
    }

    //for serialization --defines this object from a byte stream
    inline EBGraph(void * a_buf);

    //for serialization --puts this object into a byte stream
    inline void energize(void *a_buf) const;
    
    ///Chombo leveldata serialization routines
    //this one defines the size of the aformentioned byte stream
    inline
    virtual std::size_t charsize(const Box&   a_bx,
                                 unsigned int a_startcomp,
                                 unsigned int a_numcomps) const;

    ///
    inline
    virtual void  linearIn(void*                  a_buf,
                           const Box       &      a_bx,
                           unsigned int           a_startcomp,
                           unsigned int           a_numcomps);
    
    ///
    inline
    virtual void linearOut(void*               a_buf,
                           const Box       &   a_bx,
                           unsigned int        a_startcomp,
                           unsigned int        a_numcomps) const;



  private:

    //define everything as an irregular cell but with all the cells regular
    inline
    void defineAsFullIrregular(const Box& a_validBox, const Box& a_region, const Box& a_domain);
    //define everything as an irregular cell but with all the cells covered
    inline
    void defineAsEmptyIrregular(const Box& a_validBox, const Box& a_region, const Box& a_domain);
    ///
    /**
       If this is allregular or allcovered,
       the HostBoxData below is undefined.
    */
    TAG m_tag;

    ///box over which this graph is defined
    Box m_region;

    //valid box of the grid
    Box m_validBox;
    
    ///computational domain at this level of refinement
    Box m_domain;

    ///
    bool m_isDefined;

    struct localData
    {

      ///
      /**
         the graph description of the box.
         If the pointer CellNode.m_cellList == 0, the cell is covered.
         If the pointer CellNode.m_cellList == 1, the cell is regular.
         Otherwise, the pointer goes to the node of the graph.
      */
      HostBoxData<CellNode, 1> m_graph;

      ///
      PointSet m_irregIVS;

    };


    // this is null in the case of m_tag == AllRegular || AllCovered
    shared_ptr<localData> m_internals;


    };
  template <CENTERING cent> 
  inline vector<EBIndex<cent> > 
  getFaces(const  EBIndex<CELL>   & a_vof,
           const Side::LoHiSide   & a_sd,
           const EBGraph          & a_graph) 
  {
    vector<EBIndex<cent> > retval;
    PROTO_ASSERT(false, "should not get here");
    return retval;
  }

  template< > 
  inline vector<EBIndex<XFACE> > 
  getFaces(const  EBIndex<CELL>   & a_vof,
           const Side::LoHiSide   & a_sd,
           const EBGraph          & a_graph) 
  {
    return a_graph.getXFaces(a_vof, a_sd);
  }

  template< > 
  inline vector<EBIndex<YFACE> > 
  getFaces(const  EBIndex<CELL>    & a_vof,
           const Side::LoHiSide    & a_sd,
           const EBGraph           & a_graph)
  {
    return a_graph.getYFaces(a_vof, a_sd);
  }
#if DIM==3
  template< > 
  inline vector<EBIndex<ZFACE> >
  getFaces(const  EBIndex<CELL>     & a_vof,
           const Side::LoHiSide     & a_sd,
           const EBGraph            & a_graph)
  {
    return a_graph.getZFaces(a_vof, a_sd);
  }
#endif
  template <CENTERING cent> 
  inline vector<EBIndex<cent> >
  getIrregFaces(const Box& a_box, const EBGraph& a_graph)
  {
    PR_TIME("EBGraph::getIrregFaces");
    //need to do this is a set at first to avoid double counting faces.
    std::set< EBIndex<cent> > resultSet;  
    Box dombox = a_box & a_graph.getDomain();
    PointSet irreg = a_graph.getIrregCells(dombox);
    for(PointSetIterator ivsit(irreg);  ivsit.ok(); ++ivsit)
    {

      vector< EBIndex<CELL> > vofpt = a_graph.getVoFs(ivsit());
      for(int ivof = 0; ivof < vofpt.size(); ivof++)
      {
        for(SideIterator sit; sit.ok(); ++sit)
        {
          vector<EBIndex<cent> > faces=   getFaces<cent>(vofpt[ivof],sit(), a_graph);
          for(int ifac = 0; ifac < faces.size(); ifac++)
          {
            resultSet.insert(faces[ifac]);
          }
        }
      }
    }
    vector< EBIndex<cent> > faces(resultSet.size());
    int ivec = 0;
    for (auto i = resultSet.begin(); i != resultSet.end();   ++i )
    {
      faces[ivec] = ( *i );
      ivec++;
    }
    return faces;
  }


  template <CENTERING cent> 
  inline vector<EBIndex<cent> >
  getAllFaces(const Box& a_box, const EBGraph& a_graph)
  {
    PR_TIME("EBGraph::getAllFaces");
    //need to do this is a set at first to avoid double counting faces.
    std::set< EBIndex<cent> > resultSet;  
    Box dombox = a_box & a_graph.getDomain();
    PointSet irreg(dombox);
    for(PointSetIterator ivsit(irreg);  ivsit.ok(); ++ivsit)
    {
      vector< EBIndex<CELL> > vofpt = a_graph.getVoFs(ivsit());
      for(int ivof = 0; ivof < vofpt.size(); ivof++)
      {
        for(SideIterator sit; sit.ok(); ++sit)
        {
          vector<EBIndex<cent> > faces = getFaces<cent>(vofpt[ivof],sit(), a_graph);
          for(int ifac = 0; ifac < faces.size(); ifac++)
          {
            resultSet.insert(faces[ifac]);
          }
        }
      }
    }
      
    vector< EBIndex<cent> > faces(resultSet.size());
    int ivec = 0;
    for (auto i = resultSet.begin(); i != resultSet.end();   ++i )
    {
      faces[ivec] = ( *i );
      ivec++;
    }
      
    return faces;
  }

  template < >
  inline vector<EBIndex<BOUNDARY> >
  getIrregFaces(const Box& a_box, const EBGraph& a_graph)
  {
    return a_graph.getCutFaces(a_box);
  }
  template < >
  inline vector<EBIndex<BOUNDARY> >
  getAllFaces(const Box& a_box, const EBGraph& a_graph)
  {
    return a_graph.getCutFaces(a_box);
  }


  template <CENTERING cent> 
  inline vector<EBIndex<cent> >
  getAllIndicies(const Box& a_box, const EBGraph& a_graph)
  {
    return getAllFaces<cent>(a_box, a_graph);
  }

  template <CENTERING cent> 
  inline vector<EBIndex<cent> >
  getIrregIndicies(const Box& a_box, const EBGraph& a_graph)
  {
    return getIrregFaces<cent>(a_box, a_graph);
  }

  template < > inline
  vector<EBIndex<CELL> >
  getAllIndicies(const Box& a_box, const EBGraph& a_graph)
  {
    return a_graph.getAllVoFs(a_box);
  }


  template < > inline
  vector<EBIndex<CELL> >
  getIrregIndicies(const Box& a_box, const EBGraph& a_graph)
  {
    return a_graph.getIrregVoFs(a_box);
  }


  template <CENTERING cent> 
  inline vector<EBIndex<cent> >
  getSideFaces(const EBIndex<CELL> a_vof, Side::LoHiSide a_sd, const EBGraph& a_graph)
  {
    return getFaces<cent>(a_vof, a_sd, a_graph);
  }

}

#include "implem/Proto_EBGraphImplem.H"
#endif
