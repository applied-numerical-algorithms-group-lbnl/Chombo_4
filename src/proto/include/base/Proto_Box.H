#pragma once
#ifndef _PROTO_BOX_
#define _PROTO_BOX_

#include "Proto_Point.H"
#include "Proto_LoHiSide.H"
#include <cstdlib> //for size_t
#include <iostream>
#include <cmath>

namespace Proto 
{

    // forward declarations
    class BoxIterator;
#ifdef PR_HDF5
    class HDF5Handler;
#endif
    /// An interval in DIM dimensional space.
    /**
        A Box is a region in \f$ \mathbb{Z}^{DIM} \f$ specified by two corner Point objects, <code> high </code> and <code> low </code> INCLUSIVELY.
        Equivalently, a Box is a collection of \f$DIM\f$ linear intervals in \f$\mathbb{Z}^{DIM}\f$ 
        \f$[low[i] ,high[i]] : i \in [0,DIM)\f$. Boxes are written using the notation <code> [low, high] </code>.
     */
    class Box
    {
#ifdef PR_HDF5
        friend class Proto::HDF5Handler;
#endif
        public:
        typedef BoxIterator iterator; //Box::iterator aliases to BoxIterator

        ///////////////////////////////////////////////////////////////////////////////////////////  
        /** @name Constructors */
        ///@{

        /// Default Constructor
        /**
          Builds the empty Box <code>[(-1,-1,...,-1), (0,0,...,0)]</code>
          */
        inline Box();

      //useful for bc evaluation
      Box faceBox(int a_idir, const Side::LoHiSide a_side) const
      {
        Box retval = adjCellSide(a_idir, 1, a_side);
        if(a_side == Side::Lo)
        {
          retval =  retval.shift(Point::Basis(a_idir, 1));
        }
        return retval;
      }
      
        /// Two-Point Constructor
        /**
          Builds the non-trivial Box <code>[low, high]</code>

          \param a_low    low corner
          \param a_high   high corner
          */ 
        inline Box(const Point& a_low, const Point& a_high);

        /// Primitive Two-Point Constructor
        /**
          Used for building a Box on platforms where plain-old-data is more convenient
          
          \param a_low    C-Array representing this.low()
          \param a_high   C-Array representing this.high()
          */ 
        inline Box (const int* a_low, const int* a_high) { define(a_low, a_high); }

        /// Size Constructor
        /**
          Builds the non-trivial Box <code>[(0,...,0), a_sizes - (1,...,1)]</code>.
          This is the box with its low corner at the origin and <code>a_sizes[dir]</code>
          Points in direction <code>dir</code>.

          \param a_sizes   sizes
          */ 
        inline Box(const Point& a_sizes);

        /// Copy Constructor
        CUDA_DECORATION
        inline Box(const Box& a_box); 

        /// Primitive Two-Point Lazy Constructor
        /**
          Weak construct this Box if necessary

          \param a_lo   C-Array representing this.low()
          \param a_hi   C-Array representing this.high()
          */
        inline void define(const int* a_lo, const int* a_hi);
        ///@}

        ///////////////////////////////////////////////////////////////////////////////////////////////  
        /** @name Static Functions */
        ///@{
        /// Build Cube Box
        /**
          Creates the Box <code>[Point::Zeros(), Point::Ones(a_size-1)]</code>, 
          a cube in \f$\bf{Z}^{DIM}\f$ of side length a_size

          \param a_size   side length
          */
        inline static Box Cube(int a_size);

        /// Build Kernel Box
        /**
          Creates a Box of size <code> 2*a_radius + 1</code> on all sides with (0,...,0) at the center.
          Useful for iterating through boundary cells, defining symmetric Stencils, etc.

          \param a_radius the number of "rings" of cells around the center cell.

          Example:
          @code
          //DIM=2
          using namespace Proto;
          auto K0 = Box::Kernel(1);
          K0.print(); // [(-1,-1), (1, 1)] 
          auto K1 = Box::Kernel(2);
          K1.print(); // [(-2, -2), (2, 2)]
          @endcode
        */
        inline static Box Kernel(int a_radius);

        /// Stencil Index
        /**
          @private
          Returns the linear index of <code>a_pt</code> inside the <code>Box</code> defined by <code>[a_low, a_high].
          \nEquivalent to <code>Box::index</code> but does not require a <code>Box</code> instance.

          \param a_pt     A Point
          \param a_low    The low Point of a Box
          \param a_high   The high Point of a Box
          */
        CUDA_DECORATION
            inline static unsigned int sindex(Point a_pt, Point a_low, Point a_high); 

        /// Stencil Offset
        /**
          @private
          Used internally by <code>Stencil</code>, not recommended for public use.

          \param a_pt     A Point
          \param a_low    The low Point of a Box
          \param a_high   The high Point of a Box
          */
        // TODO: Doesn't appear to be used anywhere -CLG
        //CUDA_DECORATION
        //inline static unsigned int soffset(const Point& a_pt, const Point& a_low, const Point& a_high);
        ///@}

        ///////////////////////////////////////////////////////////////////////////////////////////////  
        /** @name Accessors And Queries */
        ///@{

        /// Access Low Corner
        /**
          This function returns by value; a Box cannot be altered through the output of <code>low()</code>.
          */
        CUDA_DECORATION
        inline Point low() const {return m_low;}

        /// Access High Corner
        /**
          This function returns by value; a Box cannot be altered through the output of <code>high()</code>.
          */
        CUDA_DECORATION
        inline Point high() const {return m_high;} 

        //TODO: needs renaming
        //useful for bc evaluation
        inline Point boundary(const Side::LoHiSide a_side) const;

        /// Edge Size
        /**
          Returns the "edge length" of this Box on a given axis

          /param a_dim    direction axis in <code>[0,DIM)</code>

          Example:
          @code
          //DIM=2
          using namespace Proto;
          Box B(Point::Zeros(), Point({1,3}));
          B.print(); //prints [(0,0),(1,3)]
          std::cout << B.size(0) << std::endl; //prints 2
          std::cout << B.size(1) << std::endl; //prints 4
          @endcode
        */
        CUDA_DECORATION
        inline std::size_t size(unsigned char a_dim) const;

        /// Volumetric Size
        /**
          Returns the "volume" (in number of points) of this box.

          Example:
          @code
          //DIM=2
          using namespace Proto;
          Box B(Point::Zeros(), Point({1,3}));
          B.print(); //prints [(0,0),(1,3)]
          std::cout << B.size() << std::endl; //prints 8 (=2*4)
          @endcode
        */
        CUDA_DECORATION
        inline std::size_t size() const; 
        
        /// All Sizes
        /**
          Returns all sizes at once in the form of a Point.  
        */
        CUDA_DECORATION
        inline Point sizes() const; 

        /// Contains Point Query
        /**
          Checks if a_pt is inside of *this.
          Note that <code>this->low() </code> and <code> this->high() </code> are both considered inside *this

          Example:
          @code
          //DIM=2
          using namespace Proto;
          Box B(Point::Zeros(), Point({1,3}));
          std::cout << B.contains(Point({0,0})) << std::endl; //prints true
          std::cout << B.contains(Point({1,1})) << std::endl; //prints true
          std::cout << B.contains(Point({1,3})) << std::endl; //prints true
          std::cout << B.contains(Point({1,4})) << std::endl; //prints false
          @endcode
        */
        CUDA_DECORATION
            inline bool contains(const Point& a_pt) const;

        /// Contains Box Query
        /**
          Check if Box a_rhs is a subset of *this. Returns true even if a_rhs is not a proper subset.
          Equivalent to <code> this->contains(a_rhs.low()) && this->contains(a_rhs.high() </code>

          Example:
          @code
          //DIM=2
          using namespace Proto;
          Box B0(Point({0,0}),Point({3,3}));
          Box B1(Point({1,1}),Point({2,2}));
          Box B2(Point({1,1}),Point({4,4}));
          Box B2(Point({-1,-1}),Point({2,2}));
          std::cout << B0.contains(B0) << std::endl; //prints true
          std::cout << B0.contains(B1) << std::endl; //prints true
          std::cout << B0.contains(B2) << std::endl; //prints false
          std::cout << B0.contains(B3) << std::endl; //prints false
          @endcode
        */
        CUDA_DECORATION
            inline bool contains(const Box& a_rhs) const;

        /// Point on Boundry Query
        /**
          Check of <code>a_p</code> is part of the boundary of this Box.

          \param a_p  A Point
          */
        inline bool onBoundary(const Point& a_p) const;

        /// Empty Query
        /**
          Check if *this contains no Points
          */
        CUDA_DECORATION
            inline bool empty() const {return (m_size < 1);}

        /// Point to Linear Index
        /**
          Returns a linear index in <code>[0,this->size())</code> associated with <code>a_pt</code>.
          Fails by assertion if <code>a_pt</code> is not inside <code>*this</code>. 
          Indices start at <code>this->low()</code> which has an index of 0. Indices increase
          fastest in dimension 0 and slowest in dimension <code>DIM-1</code>. <code>this->high()</code>
          has an index of <code>this->size()-1</code> as is consistent with the < operator of <code>Point</code>.

          \param a_pt     a Point inside *this
          */
        CUDA_DECORATION
            inline unsigned int index(const Point& a_pt) const;

        /// Access Point by Index
        /**
          Return a Point associated with <code> a_index </code>.
          This is the inverse of the <code>index(Point)</code> function. 

          \param a_index  an index in <code>[0,this->size())</code>
          */
        CUDA_DECORATION
            inline Point operator[](unsigned int a_index) const;


        CUDA_DECORATION
            inline Point operator()(unsigned int a_idx, unsigned int a_idy, unsigned a_idz) const;
        ///@}



        ///////////////////////////////////////////////////////////////////////////////////////////////  
        /** @name Operators */
        ///@{

        ///  Intersection Operator 
        /**
          Returns a new Box which is the intersection of <code>*this</code> and <code>*a_rightBox</code>.

          \param a_rightBox  Another Box

          Example:
          @code
          //DIM=2
          using namespace Proto;
          Box B0 = Box::Cube(3);
          Box B1 = Box::Cube(3).shift(Point::Ones());
          Box B2 = B0 & B1;
          B0.print(); //prints [(0,0),(2,2)]
          B1.print(); //prints [(1,1),(3,3)]
          B2.print(); //prints [(1,1),(2,2)]
          @endcode
        */
        inline Box operator&(const Box& a_rightBox) const;


        /// In Place Intersection
        /**
          Performs intersection on *this with <code>a_rhs</code>

          \param a_rhs    Another Box 
          */
        inline void operator&=(const Box& a_rhs); 

        /// Add Point Operator
        /**
          Returns a new Box which has been grown such that it contains a_pt. 
          if <code>this->contains(a_pt)</code> this function returns *this.

          \param a_pt   Point to be added

          Example:
          @code
          //DIM=2
          using namespace Proto;
          auto B0 = Box::Cube(3);
          B0.print(); //[(0,0), (2,2)]
          Point p(-1,3);
          auto B1 = B0 & p;
          B1.print(); //[(-1,0), (2,3)]
          @endcode
        */
        inline Box operator&(const Point& a_pt) const;

        /// In Place Add Point Operator
        /**
          Alters <code>*this</code> so that it now includes <code>a_pt</code>.
          If <code>this->contains(a_pt)</code>, this function does nothing.

          \param a_pt     A Point
          */
        inline void operator&=(const Point& a_pt);

        //inline void operator&=(const Box& a_rightBox);

        /// Equality Operator
        /**
          Two Boxes are considered equal if they have identical (==) low() and high()

          \param a_rhsBox     A Box to test equality with
          */
        inline bool operator==(const Box& a_rhsBox) const;

        /// Inequality Operator
        inline bool operator!=(const Box& a_rhsBox) const { return !(*this == a_rhsBox); }

        /// Less Than Operator
        /**
          Establishes an ordering for Boxes. This ordering is defined by comparing
          the low corner Point of the two Boxes. 

          \param a_rhsBox     A Box.
          */
        inline bool operator<(const Box& a_rhsBox) const { return m_low < a_rhsBox.low(); }

        /// Modulus Operator
        /**
          Convenience operator for <code>Box::mod</code>.
          Returns the periodic image of <code>a_pt</code> contained in <code>*this</code>

          \param a_pt   A Point to mod by    
          */
        CUDA_DECORATION
        inline Point operator%(const Point& a_pt) const;

        /// Modulus Function
        /**
          Returns the periodic image of a_pt that is inside *this.

          \param a_pt   Point divisor

          Example:
          @code
          //DIM = 2;
          using namespace Proto;
          Box B = Box::Cube(3); // [(0,0), (2,2)]
          std::cout << B.mod(Point({3,3})) << std::endl; //prints (0,0)
          std::cout << B.mod(Point({-1,-1})) << std::endl; //prints (2,2)
          @endcode
        */
        CUDA_DECORATION
        inline Point mod(const Point& a_pt) const;

        ///@}

        ///////////////////////////////////////////////////////////////////////////////////////////////  
        /** @name Transformations */
        ///@{

        /// Shift Transformation
        /**
          Creates a new Box shifted in <code>a_direction</code> by <code>a_offset</code> 
          units relative to <code>*this</code>.     

          \param a_direction  an int in /f$[0,DIM)/f$ specifying a direction
          \param a_offset     number of points to offset. Use a negative value for a negative offset.

          Example:
          @code
          //DIM=2
          using namespace Proto;
          Box B1 = Box::Cube(2);
          Box B2 = B1.shift(0,1);
          Box B3 = B1.shift(1,-1);
          B1.print(); //prints [(0,0), (1,1)]         
          B2.print(); //prints [(1,0), (2,1)]         
          B3.print(); //prints [(0,-1), (1,0)]
          @endcode         
        */
        inline Box shift(int a_direction, int a_offset) const;

        /// Point Shift Transformation
        /**
          Creates a new Box shifted by <code>a_pt</code> relative to <code>*this</code>.
          New Box is: <code>[low() + a_pt, high() + a_pt]</code>
          \param a_pt     shift offset

          Example:
          @code
          //DIM=2
          using namespace Proto;
          Box B1(Point::Zeros(),Point::Ones());
          Box B2 = B1.shift(Point({2,-3}));
          B1.print(); //prints [(0,0), (1,1)]
          B2.print(); //prints [(2,-3), (3,-2)]
          @endcode
        */
        inline Box shift(const Point& a_pt) const;

        /// Shift To Origin
        /**
          Creates a new Box shifted such that <code>this->low() == (0,...,0)</code>.
          Equivalent to <code>this->shift(this->low())</code>.
          */
        inline Box toOrigin() const;

        /// Isotropic Grow Operation
        /**
          Returns a new Box which is larger in every direction by a_numpoints units.
          If a_numpoints is negative, the box will shrink. 
          \param a_numpoints  number of points to grow by in each direction

          Example:
          @code
          //DIM=2
          using namespace Proto;
          Box B0(Point::Zeros(),Point::Ones());
          Box B1 = B0.grow(3);
          Box B2 = B1.grow(-2);
          B0.print(); //prints [(0, 0), (1, 1)]
          B1.print(); //prints [(-3,-3), (4, 4)]
          B2.print(); //prints [(-1,-1), (2, 2)]
          @endcode
        */
        inline Box grow(int a_numpoints) const;

        /// Anisotropic Grow Operation
        /**
          Returns the reshaped Box: [low() - a_pt, high() + a_pt]
          To grow a box only upwards / downwards, see <code>Box::extrude</code>.

          \param a_pt growth offset

          Example:
          @code
          //DIM=2
          using namespace Proto;
          Box B0(Point::Zeros(),Point::Ones()*3);
          Box B1 = B0.grow(Point({-1,1})); //shrink in dimension 0, grow in dimension 1
          B0.print(); //prints [(0, 0), (3, 3)]
          B1.print(); //prints [(1, -1), (2, 4)]
          @endcode
        */
        inline Box grow(const Point& a_pt) const;

        /// Anisotropic Grow (Int Version)
        /**
          Returns grow(Point::Basis(a_dir, a_dist))
          Useful in places where plain-old-data inputs are preferable.
          To grow a grow a box only upwards / downwards, see <code>Box::extrude</code>.

          \param a_dir     A direction in [0,DIM)
          \param a_dist  Distance to grow
          */
        inline Box grow(int a_dir, int a_dist) const;

        /// One-Sided Grow
        /**
          Grows the box along coordinate <code>a_dir</code> a distance <code> a_dist </code>
          on side <code>a_side</code> only. This function is equivalent to
          <code>extrude(Point::Basis(a_dir, sign(a_side)*a_dist)) </code>.

          \param a_dir    A direction in [0, DIM)
          \param a_side   Either Side::Lo or Side::Hi
          \param a_dist   How far to grow
          */ 
        inline Box grow(int a_dir, Side::LoHiSide a_side, int a_dist) const; 

        /// Extrude
        /**
          AKA directional grow. Returns a new Box with the <code> a_dir </code> direction extruded a distance <code> a_dist </code).

          \param a_dir    Direction to extrude in.
          \param a_dist   (Optional) Distance to extrude (default: 1)

          Example:
          @code
          //DIM=3
          using namespace Proto;
          Box B0 = Box::Cube(4).shift(Point::Ones());
          Box B1 = B0.extrude(Point::Ones(), 2);
          Box B2 = B0.extrude(Point::Basis(0,-1),3);
          Box B3 = B0.extrued(Point(-1, 1, 0);
          B0.print(); //prints [(1, 1, 1), (4, 4, 4)]
          B1.print(); //prints [(1, 1, 1), (6, 6, 6)]
          B2.print(); //prints [(-2, 1, 1), (4, 4, 4)]
          B3.print(); //prints [(0, 1, 1), (4, 5, 4)]
          @endcode
        */
        inline Box extrude(const Point& a_dir, int a_dist = 1) const;

        //TODO: Better to use grow(dir, side, dist) in my opinion -CLG
        /// Extrude (Primitive Version)
        /**
          AKA directional grow. Returns a new Box with the <code> a_dir </code> direction extruded a distance <code> a_dist </code).
          By default, the extrusion is upwards. If <code>a_upper</code> is false, extrudes downwards.

          \param a_dir    Axis to extrude along
          \param a_dist   (Optional) Distance to extrude (default: 1)
          \param a_upper  (Optional) Extrude upwards? (default: true)

          Example:
          @code
          //DIM=2
          using namespace Proto;
          Box B0 = Box::Cube(4).shift(Point::Ones());
          Box B1 = B0.extrude(2,3,false);
          Box B2 = B0.extrude(2,3,true);
          B0.print(); //prints [(1, 1, 1), (4, 4, 4)]
          B1.print(); //prints [(1, 1, -2), (4, 4, 4)]
          B2.print(); //prints [(1, 1, 1), (4, 4, 7)]
          @endcode
        */
        inline Box extrude(int a_dir, int a_dist = 1, bool a_upper = true) const;

        /// Grow High Side
        /**
          Equivalent to <code>grow(a_dir, Side::Hi, a_dist)</code>

          \param a_dir    Coordinate axis to grow along
          \param a_dist   Distance to grow
          */
        inline Box growHi(int a_dir, int a_dist) const;

        /// Grow Low Side
        /**
          Equivalent to <code>grow(a_dir, Side::Lo, a_dist)</code>

          \param a_dir    Coordinate axis to grow along
          \param a_dist   Distance to grow
          */
        inline Box growLo(int idir, int igrow) const;

        /// Isotropic Coarsen Operation
        /**
          Returns a new Box coarsened by a factor of <code>a_numpoints</code>. Fails if user tries to 
          coarsen using a non-positive ratio.
          If the limits of <code>*this</code> are not multiples of the coarsening ratio
          (e.g. if <code>this-> coarsenable() != true</code>) the resulting Box is not guaranteed to
          be a subset of <code>*this</code>. In this situation, use <code>Box::taperCoarsen</code>

          \param a_ratio     Coarsening ratio

          Example:
          @code
          //DIM=2
          using namespace Proto;
          Box B0 = Box::Cube(4);
          Box B1 = B0.coarsen(2);
          B0.print(); //prints [(0, 0), (3, 3)]
          B1.print(); //prints [(0, 0), (1, 1)]
          Box B2 = Box::Cube(3).shift(Point::Ones(2));
          Box B3 = B2.coarsen(2);
          B2.print(); //prints [(2, 2), (4, 4)]
          B3.print(); //prints [(1, 1), (2, 2)]
          @endcode
        */
        inline Box coarsen(unsigned int a_ratio) const;

        /// Anisotropic Coarsen Operation
        /**
          Returns a new Box coarsened in each direction according to <code>a_pt</code>.
          Fails if user tries to coarsen using a non-positive ratio.
          (*this).coarsen(a_pt).refine(a_pt) will always contain *this.

          \param a_pt    Coarsening ratios

          Example:
          @code
          //DIM=2
          using namespace Proto;
          Box B0 = Box::Cube(4);
          //coarsen in only the 1 direction
          Box B1 = B0.coarsen(Point({1,2}));
          B1.print(); //prints [(0,0),(3,1)]
          @endcode
        */
        inline Box coarsen(const Point& a_pt) const;

        /// (Isotropic) Coarsenable Query
        /**
          Returns true if *this will coarsen normally by a given ratio. 

          \param a_ratio    Coarsening ratio
          */
        inline bool coarsenable(const int& a_ratio) const;

        /// (Anisotropic) Coarsenable Query
        /**
          Returns true if *this will coarsen normally by a given ratio. 

          \param a_ratio    Coarsening ratio
          */
        inline bool coarsenable(const Point& a_ratio) const;

        /// Tapered Coarsen
        /**
          This function is identical to Box::coarsen when Box::coarsenable() is true.
          For non-coarsenable Boxes, the new Box's limits are rounded such that the result
          Is always a subset of *this. Specifically, the lower limits are always rounded UP.
          This function is used heavily in the automatic computation of domains and ranges
          for pointwise and Stencil operators, but generally has limited use otherwise.

          \param a_ref  Coarsening ratio
          */
        inline Box taperCoarsen(const Point& a_ref) const;

        /// Anisotropic Refine Operation
        /**
          Returns a new Box refined in each direction according to a_pt.
          \param a_pt    Refinement ratios

          Example:
          @code
          //DIM=2
          using namespace Proto;
          Box B0 = Box::Cube(2);
          Box B1 = B0.refine(Point({1,2}));
          B0.print(); //prints [(0, 0), (1, 1)]
          B1.print(); //prints [(0, 0), (1, 3)]
          @endcode
        */
        inline Box refine(const Point& a_pt) const;

        /// Isotropic Refine Operation
        /**
          Returns a new Box refined in all directions by a_numpoints.
          
          \param a_ratio     Refinement ratio

          Example:
          @code
          //DIM=2
          using namespace Proto;
          Box B0 = Box::Cube(2);
          Box B1 = B0.refine(2);
          B0.print(); //prints [(0, 0), (1, 1)]
          B1.print(); //prints [(0, 0), (3, 3)]

          Box B2 = Box::Cube(2).shift(Point::Ones());
          Box B3 = B2.refine(2);
          B2.print(); //prints [(1, 1), (2, 2)]
          B3.print(); //prints [(2, 2), (5, 5)]
          @endcode
        */
        inline Box refine(unsigned int a_ratio) const;

        /// Edge
        /**
          Returns the subset on the boundary of *this in a given direction with specified thickness.
          Very handy for adjacency and boundary based computations. This function is slightly more general
          than <code>face(...)</code> in that it can return regions of all codimensionalities (faces, edges, corners, etc.)
          The returned <code>Box</code> is always a subset of <code>*this</code>

          \param a_dir  Direction of desired edge with respect to the cell center of *this
          \param a_dist (Optional) Thickness of the output (default: 1)

          Examples:
          @code
          //DIM = 2
          using namespace Proto;
          auto B0 = Box::Cube(4).shift(Point::Ones()); //[(1,1), (4,4)]
          auto B1 = B0.edge(Point::Basis(0));          //[(4,1), (4,4)]
          auto B2 = B0.edge(Point::Ones(), 2);         //[(3,3), (4,4)]
          auto B3 = B0.edge(Point::Basis(1,-1), 2);    //[(1,1), (4,2)]
          @endcode
        */
        inline Box edge(const Point& a_dir, int a_dist = 1) const;

        /// Face
        /**
          Returns the face (e.g. region of codimension 1) of this box in a specified direction.
          The output box is always contained in *this. 
          For adjacent faces, use <code>adjacent(...)</code>. 
          If <code>a_dist >= size(a_dir) </code>, the rturned box = *this.

          \param a_dir     Coordinate direction in [0, DIM)
          \param a_side    Side::Lo or Side::Hi
          \param a_dist    (Optional) Thickness of the output. (default: 1)

          Examples:
          @code
          //DIM = 2
          using namespace Proto;
          auto B0 = Box::Cube(4).shift(Point::Ones()); //[(1,1), (4,4)]
          auto B1 = B0.face(0, Side::Hi);              //[(4,1), (4,4)]
          auto B2 = B0.face(1, Side::Lo, 2);           //[(1,1), (4,2)]
          auto B3 = B0.face(0, Side::Lo, 5);           //[(1,1), (4,4)]
          @endcode
        */
        inline Box face(int a_dir, Side::LoHiSide a_side, int a_dist = 1) const;

        //TODO: redundant with face, unused
        //Box faceBox(int a_idir, const Side::LoHiSide a_side) const;

        //TODO: Redundant with face
        /// Flatten
        /**
          Returns a copy of *this with dimension a_dir flattened to a thickness of 1.
          Useful for creating Boxes of dimensionality less than DIM (e.g. a plane in 3D).

          \param a_dir    Direction to flatten
          \param a_upper  Flatten upwards?

          Example:
          @code
          //DIM=3
          using namespace Proto;
          Box B0 = Box::Cube(4).shift(Point::Ones());
          Box B1 = B0.flatten(2);
          Box B2 = B0.flatten(2,true);
          B0.print(); //prints [(1, 1, 1), (4, 4, 4)]
          B1.print(); //prints [(1, 1, 1), (4, 4, 1)]
          B2.print(); //prints [(1, 1, 4), (4, 4, 4)]
          @endcode
        */
        inline Box flatten(const int a_dir, bool a_upper = false) const;

        /// Adjacent Cells
        /**
          Returns a box adjacent to *this in a given direction with a thickness
          a_dist in the normal direction. If there are multiple non-zero entries
          in a_dir, a_dist will be applied to the thickness in all of those directions

          \param a_dir  "Normal" direction of the desired adjacent cell
          \param a_dist (Optional) "thickness" of the desired adjacent cell (default: this->size(a_dir)) 

          Examples:
          @code
          //DIM = 2;
          Box B0 = Box::Cube(8);  // [(0,0) ,  (7,7)]
          B0.adjacent((1,0) , 2); // [(8,0) ,  (9,7)]
          B0.adjacent((0,-1), 2); // [(0,-2), (7,-1)]
          B0.adjacent((-1,1), 2); // [(-2,8), (-1,9)]
          B0.adjacent((1,0));     // [(8,0) , (15,7)], 8 x 8 cube directly adjacent in (1, 0) direction
          @endcode 
        */
        inline Box adjacent(const Point& a_dir, int a_dist) const;

        /// Adjacent Cells (Side)
        /**
          Overload of adjacent that uses a Side input.

          \param a_dir    Coordinate direction (0 for x direction, etc.)
          \param a_side   Low or high side   
          \param a_dist (optional) "thickness" of the desired adjacent cell (default: this->size(a_dir))

          Examples:
          @code
          //DIM = 2;
          Box B0 = Box::Cube(8);  // [(0,0) ,  (7,7)]
          B0.adjacent((0, Side::Hi, 2); // [(8, 0), (9, 7)]
          B0.adjacent((0, Side::Lo, 2); // [(0,-2), (7,-1)]
          B0.adjacent((0, Side::Hi);    // [(8, 0), (15,7)]
          @endcode 
        */
        inline Box adjacent(int a_dir, Side::LoHiSide a_side, int a_dist) const;

        //TODO: Redundant
        inline Box adjCellLo(int a_dir, int a_dist) const;

        //TODO: Redundant
        inline Box adjCellHi(int a_dir, int a_dist) const;

        //TODO: Redundant
        inline Box adjCellSide(int a_idir, int a_length, Side::LoHiSide a_sd) const
        {
            Box retval = adjCellHi(a_idir, a_length);
            if(a_sd == Side::Lo)
            {
                retval = adjCellLo(a_idir, a_length);
            }
            return retval;
        }

        ///@}

        ///////////////////////////////////////////////////////////////////////////////////////////////  
        /** @name Utility */
        ///@{

        /// Iterator Begin
        /**
          See documentation for Proto::BoxIterator for a basic usage example.
          */
        inline BoxIterator begin() const; 

        /// Iterator End
        /**
          See documentation for Proto::BoxIterator for a basic usage example.
          */
        inline BoxIterator end() const; 

        /// Iterator Reverse Begin
        /**
          See documentation for Proto::BoxIterator for a basic usage example.
          */
        inline BoxIterator rbegin() const; 

        /// Iterator Reverse End
        /**
          See documentation for Proto::BoxIterator for a basic usage example.
          */
        inline BoxIterator rend() const; 

        /// Linear Size
        /**
          Computes the number of bytes needed to store *this in linear buffer.
          */
        inline size_t linearSize() const;

        /// Read From Buffer
        /**
          Reads data into *this from a buffer
          */
        inline void linearIn(const char* a_buf);

        /// Write To Buffer
        /**
          Writes the data in *this to a buffer
          */
        inline void linearOut(char* a_buf) const;

        /// Print
        /**
          Prints *this using the format <code>[low, high]</code>
          */ 
        inline void print() const; 
        ///@}

        private:  

        inline void recomputeSize(); ///< Used to reevaluate the size of the box when it is changed.

        Point m_low; ///< Point object containing the lower bounds of the Box.
        Point m_high;  ///< Point object containing the upper bounds of the Box.
        int m_size; ///< "Volume" of the box.

    }; //end class Box

     /// OStream Operator 
    inline std::ostream& operator<<(std::ostream& a_os, const Box& a_box)
    {
        a_os << "[" << a_box.low() << "," << a_box.high() << "]";
        return a_os;
    }

     /// OStream Operator 
  inline std::istream& operator>>(std::istream& a_is, Box& a_box)
    {

      Point low, high;
      char dum1, dum2, dum3;
      a_is >> dum1 >> low >> dum2 >> high >> dum3;
      a_box = Box(low, high);
      return a_is;
    }
  
    /// Iterator for Boxes
    /**
       Iteration class which conforms to most of the syntax of std Iterator implementations.

       Example:
       @code
       Box box = Box::Cube(8);
       // STL style iteration
       for (auto biter = box.begin(); biter != box.end(); ++biter)
       {
           (*biter).print(); //prints all the points in box 
       }
       // Alternate style
       for (auto biter = box.begin(); biter.ok(); ++biter)
       {
           (*biter).print(); //prints all the points in box 
       }
       @endcode
    */
    class BoxIterator
    {
        public:
            /// Default Constructor
            inline BoxIterator(){};

            /// Index Constructor
            inline BoxIterator(const Box& a_box, int a_pos = 0);

            /// Point Constructor
            inline BoxIterator(const Box& a_box, const Point& a_pos);
            
            /// Reset Iterator
            inline void begin() { m_pos = 0; }

            /// Continue Iteration Query
            inline bool ok() const;

            /// Equality Operator
            /**
              Returns true for two iterators with identical iteration index and Box.
              The Box instances need not be the same. 
            */
            inline bool operator==(const BoxIterator& a_iter) const;

            /// Inequality Operator
            inline bool operator!=(const BoxIterator& a_iter) const { return !((*this) == a_iter); }

            /// Dereference Iterator
            /**
              Returns a Point
            */
            inline Point operator*() const { return m_box[m_pos]; };

            /// Get Index
            /**
              Returns the current index of the point *(*this).  
            */
            inline int operator()() const { return m_pos; };

            /// Prefix Increment Iterator
            inline BoxIterator& operator++();

            /// Prefix Decrement Iterator
            inline BoxIterator& operator--();

            /// Postfix Increment Iterator
            inline BoxIterator operator++(int);

            /// Postfix Decrement Iterator
            inline BoxIterator operator--(int);

        private:
            Box m_box;
            int m_pos;
    }; //end class BoxIterator

    #include "implem/Proto_BoxImplem.H"
} //end namespace Proto
#endif //end include guard
