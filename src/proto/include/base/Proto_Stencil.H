#ifndef _PROTO_STENCIL_H_
#define _PROTO_STENCIL_H_
#include "Proto_Timer.H"
#include "Proto_BoxData.H"
#include "Proto_DisjointBoxLayout.H"

#include <vector>
#include <tuple>
#include <iostream>
#include <iomanip> //for pretty printing
#include <set>

//biggest stencil size I have seen:
#define PR_MAX_COEFFS 343

namespace Proto {

    // Forward declarations
    template <typename T>
    class Stencil;
    
    template <typename T, unsigned int C, MemType MEMTYPE, unsigned char D, unsigned char E>
    class BoxData;

//=======================================================================================
// SHIFT ||
//=======++
    
    /** @defgroup stencil_operations Stencil Operations*/
    /*@{*/

    /// Stencil Shift
    /**
      \ingroup stencil_operations
      A shift is an alias for a Point which is used solely to provide
      a fluent syntax for creating Stencils. Refer to the documentation
      for Stencil for an illustrative example.
      */
    class Shift {
        public:

            ///////////////////////////////////////////////////////////////////////////////////////////////
            /** @name Constructors */
            ///@{

            /// Default Constructor
            Shift() : m_shift(Point::Zeros()){};

            /// Point Constructor
            /**
              \ingroup stencil_operations
              Builds a Shift from a Point

              \param a_pt   Point to base this Shift on
              */
            explicit Shift(const Point& a_pt){m_shift = a_pt;};

            /// Variadic Constructor
            /**
              \ingroup stencil_operations
              More or less identical to the variadic Point constructor.
              Builds a Shift using the first DIM arguments and ignores the rest.

              \param args   At least DIM int arguments, the first DIM of which define *this
              */
            template<typename... vals>
                inline explicit Shift(vals... args) : m_shift(args...) {}

            ///@}
            ///////////////////////////////////////////////////////////////////////////////////////////////
            /** @name Methods */
            ///@{

            /// Basis Shift
            /**
              \ingroup stencil_operations
              Shortcut for <code>Shift(Point::Basis(...))</code>

              \param a_dir    Direction of basis vector in [0,DIM)
              \param a_scale  (Optional) Amount to scale the basis vector by. (Default: 1).
              */
            inline static Shift Basis(int a_dir, int a_scale=1) { return Shift(Point::Basis(a_dir, a_scale)); }

            /// Zero Shift
            /**
              \ingroup stencil_operations
              Shortcut for <code>Shift(Point::Zeros())</code>
              */
            inline static Shift Zeros() { return Shift(Point::Zeros()); }

            /// Unit Shift
            /**
              \ingroup stencil_operations
              Shortcut for <code>Shift(Point::Ones(...))</code>

              \param a_scale    (Optional) Used to scale the output. (Default: 1).
            */
            inline static Shift Ones(int a_scale=1){return Shift(Point::Ones(a_scale));};

            /// Get Shift Point
            inline Point& shift(){ return m_shift; }

            /// Scalar Multiplication
            /**
              \ingroup stencil_operations
              Generates a Stencil<T> from the product of a scalar T coefficient and a Shift.
              This operator is what allows for a fluid Stencil construction syntax:

              \param a_coef     A coefficient

              Examples:
              @code
              // DIM == 1
              Stencil<double> S1 = 1.0*Shift(-1) - 2.0*Shift(0) + 1.0*Shift(1);

              // DIM - independant
              Stencil<double> S2 = (-2*DIM)*Shift::Zeros();
              for (int dir = 0; dir < DIM; dir++)
              {
                  S2 += 1.0*Shift::Basis(dir, -1);
                  S2 += 1.0*Shift::Basis(dir, +1);
              }
              @endcode
            */
            template<typename T>
                inline Stencil<T> operator*(T a_coef) const { return Stencil<T>(*this,a_coef); }

            /// Convolution
            /**
              \ingroup stencil_operations
              The product of two Shift objects is defined as their sum.
              
              \param a_shift    Another Shift.
            */
            inline Shift operator*(const Shift& a_shift) const { return Shift(m_shift + a_shift.m_shift); }

            /// Componentwise Access
            inline int& operator[](int a_dir)
            {
                PROTO_ASSERT(((a_dir >= 0) && (a_dir < DIM)),
                        "Shift::operator[] | Error: \
                        a_dir = %i is invalid. a_dir must be in [0,DIM=%i)",
                        a_dir, DIM);
                return m_shift[a_dir];
            }

            ///@}

        private:

            Point m_shift;
    
    }; // end class Shift

//=======================================================================================
// LAZYSTENCIL ||
//=============++

    // forward declaration

    /// An Unevaluated Stencil Operation
    /**
      \ingroup stencil_operations
      LazyStencil is an intermediate structure that holds the intermediate data for a
      Stencil operation.

      LazyStencil is not explicitly part of the user interface, and is only public by
      virtue of necessity.
    */
    template <typename T, unsigned int C, MemType MEMTYPE, unsigned char D, unsigned char E>
    struct LazyStencil {
        
        inline LazyStencil() {};
        inline LazyStencil(const Stencil<T>* a_stencil,
                const BoxData<T,C,MEMTYPE,D,E>* a_src,
                Box a_box, T a_scale);

        inline void apply(BoxData<T,C,MEMTYPE,D,E>& a_dest, bool a_overwrite);

        inline unsigned long long int size(){return m_stencil.size();}

        Box m_range;

        std::vector<const Stencil<T>*> m_stencil;
        std::vector<BoxData<T, C,MEMTYPE, D, E>*> m_src;
        std::vector<Box> m_box;
        std::vector<T> m_scale;
    };

//=======================================================================================
// STENCIL ||
//=========++

/// A Linear Stencil Operation
    /**
        \ingroup stencil_operations
        Encapsulates a linear stencil operation where coefficients are of type T.
        Stencil objects are built and used in a way that conforms to their nature as operators.
        For illustrative usage examples, refer to the following code snippets:

        Examples:
        Build a Stencil from Shifts and coefficients:
        \snippet Snippets.cpp proto_stencil_build
        Apply a Stencil with no Source / Dest Refinement to a BoxData:
        \snippet Snippets.cpp proto_stencil_apply

        The above examples illustrate the usage of Stencil to do computations in which the source and
        destination arrays are at the same refinement. Stencil is also capable of "prolong" and "restrict"
        type operations in which the destination and source respectively are refined. This functionality is
        useful for operations such as averaging and interpolation, or for algorithms like Red-Black Gauss Seidel
        Iteration in which it is necessary to evaluate a Stencil on "every other cell".

        To facilitate these more exotic operations, the Stencil API allows the user to designate a source and/or destination
        refinement ratio. If these values are different from (1,...,1), then input Box object will be interpreted as an
        "index range" instead of a physical Box domain. The following code snippets illustrate some examples of this functionality.

        Examples:
        Non-Trivial Source Refinement Ratio
        \snippet Snippets.cpp proto_stencil_average
        Non-Trivial Destination Refinement Ratio
        \snippet Snippets.cpp proto_stencil_dest_refine

        In the case of non-trivial destination refinement, an array of refRatio^DIM Stencils is often needed to fully populate the output.
        Proto provides a convenience structure called InterpStencil which is designed to mitigate this form of pedantry. See the associated
        documentation for additional information and example usage.
    */
    template <typename T>
    class Stencil {

        public:

        // Friend declarations
        template <typename TT, unsigned int C, MemType MEMTYPE , unsigned char D, unsigned char E>
        friend class LazyStencil;

        template <typename TT, unsigned int C , MemType MEMTYPE, unsigned char D, unsigned char E>
        friend class BoxData;

        template <typename TT, unsigned int C, MemType MEMTYPE, unsigned char D, unsigned char E>
        friend BoxData<TT,C,MEMTYPE,D,E>& operator|=(
                BoxData<TT,C,MEMTYPE,D,E>& a_dest,
                LazyStencil<TT,C,MEMTYPE,D,E>&& a_op);

        template <typename TT, unsigned int C, MemType MEMTYPE, unsigned char D, unsigned char E>
        friend BoxData<TT,C,MEMTYPE,D,E>& operator+=(
                BoxData<TT,C,MEMTYPE,D,E>& a_dest,
                LazyStencil<TT,C,MEMTYPE,D,E>&& a_op);
        
        ///////////////////////////////////////////////////////////////////////////////////////////
        /** @name Constructors */
        ///@{

        /// Default Constructor
        Stencil();

        /// General Constructor
        /**
          \ingroup stencil_operations
          Creates a Stencil with a single shift and coefficent.

          Not recommended for public use; see the Stencil class documentation for examples
          of how to build a Stencil with Shift - coefficient syntax.

          \param a_shift        Shift of this operation
          \param a_coef         Coefficient of this operation
          \param a_destRefratio (Optional) Destination refinement ratio. [Default: (1,...,1) ]
          \param a_destShift    (Optional) Destination shift.            [Default: (0,...,0) ]
          \param a_srcRefratio  (Optional) Source refinement ratio.      [Default: (1,...,1) ]
          */
        Stencil(Shift   a_shift,
                T       a_coef,
                Point   a_destRefratio  = Point::Ones(),
                Point   a_destShift     = Point::Zeros(),
                Point   a_srcRefratio   = Point::Ones());
        
        // Destructor
        ~Stencil();

        ///@}

        ///////////////////////////////////////////////////////////////////////////////////////////////
        /** @name Operators */
        ///@{

        /// Equality Operator
        /**
          \ingroup stencil_operations
          Equality between Stencils is determined by value

          \param a_stencil  Another Stencil
          */
        inline bool operator==(const Stencil<T>& a_stencil) const;

        /// Inquality Operator
        /**
          \ingroup stencil_operations
          \param a_stencil  Another Stencil
          */
        inline bool operator!=(const Stencil<T>& a_stencil) const {return !(*this == a_stencil);}
        
        /// Stencil Composition
        /**
          \ingroup stencil_operations
          The product of two Stencils is defined as their composition.

          \param a_stencil  Another Stencil
        */
        Stencil<T> operator*(const Stencil<T>& a_stencil) const;

        /// Scalar Multiplication
        /**
          \ingroup stencil_operations
          The product of a Stencil S and a coefficient v results in the scaling of all
          the coefficients of S by v.

          \param a_coef   Scaling coefficient
        */
        Stencil<T> operator*(const T a_coef) const;

        /// In Place Stencil Composition
        /**
          \ingroup stencil_operations
          \param a_stencil  Another Stencil
          */
        void operator*=(const Stencil<T>& a_stencil);

        /// In Place Scalar Multiplication
        void operator*=(const T a_coef);

        /// Stencil Addition
        /**
          \ingroup stencil_operations
          Adding two Stencils results in a new Stencil which is the union of the coefficent-Shift pairs of the inputs.
          If the two input Stencils share a common Shift, the associated coefficients will be added together.

          \param a_stencil  Another Stencil
          */
        Stencil<T> operator+(const Stencil<T>& a_stencil) const;

        /// Stencil Subtraction (Convenience)
        /**
          \ingroup stencil_operations
          Equivalent to adding <code> a_stencil*(-1) </code>

          \param a_stencil  Another Stencil
          */
        Stencil<T> operator-(const Stencil<T>& a_stencil) const;

        /// In Place Stencil Addition
        /**
          \ingroup stencil_operations
          \param a_stencil  Another Stencil
          */
        void operator+=(const Stencil<T>& a_stencil);

        /// In Place Stencil Subtraction
        /**
          \ingroup stencil_operations
          \param a_stencil  Another Stencil
          */
        void operator-=(const Stencil<T>& a_stencil);

        ///@}

        ///////////////////////////////////////////////////////////////////////////////////////////////
        /** @name Accessors and Queries*/
        ///@{

        /// Get Vector of Coefficients
        /**
          \ingroup stencil_operations
          Returned vector is read-only. Ordering corresponds to the output of Stencil::offsets()
          */
        inline const std::vector<T>& coefs() const {return m_coefs;};

        /// Get Vector of Coefficients (Primitive)
        /**
            @private
            \ingroup stencil_operations
            primitive version for testing purposes
            CURRENTLY UNUSED?
        */
        //inline const T* devCoefs() const {return d_coeff;};

        /// Get Vector of Offsets
        /**
          \ingroup stencil_operations
          Returned vector is read-only. Ordering corresponds to the output of Stencil::coefs()
          */
        inline const std::vector<Point>& offsets() const {return m_offsets;};

        /// Get Vector of Offsets
        /**
          @private
          \ingroup stencil_operations
          primitive version for testing purposes  
            CURRENTLY UNUSED?
        */
        //inline const Point* devOffsets() const {return d_offset;};

        /// Size
        /**
          \ingroup stencil_operations
          Defined as the number of coefficient-offset pairs in the stencil
          */
        inline unsigned long long int size() const {return m_coefs.size();}

        /// Span
        /**
          \ingroup stencil_operations
          Returns a Box which bounds all offsets in *this. Useful for automated ghost-cell checks.

          Example:
          @code
          //DIM = 2
          using namespace Proto;
          Stencil<T> S = 1.0*Shift(0,-1) +
          2.0*Shift(0,2) +
          3.0*Shift(1,3);
          std::cout << S.span() << std::endl; //prints [(-1,0), (2,3)]
          @endcode
        */
        inline Box span() const {return m_span;};

        /// Ghost
        /**
          Returns the span as a Point. Useful when interested only in isotropic ghost regions.
          This function will always allow for *at least* the correct number of ghost cells, and will
          overestimate in the case of non-isotropic spans.
        */
        inline Point ghost() const;

        /// Get Source Refinement Ratio
        /**
          \ingroup stencil_operations
        */
        inline Point& srcRatio(){return m_srcRefratio;};

        /// Get Source Refinement Ratio (Const)
        /**
          \ingroup stencil_operations
        */
        inline const Point& srcRatio() const {return m_srcRefratio;};

        /// Get Destination Refinement Ratio
        /**
          \ingroup stencil_operations
        */
        inline Point& destRatio(){return m_destRefratio;};

        /// Get Destination Refinement Ratio (Const)
        /**
          \ingroup stencil_operations
        */
        inline const Point& destRatio() const {return m_destRefratio;};

        /// Get Destination Shift
        /**
          \ingroup stencil_operations
        */
        inline Point& destShift(){return m_destShift;};

        /// Get Destination Shift (Const)
        /**
          \ingroup stencil_operations
        */
        inline const Point& destShift() const {return m_destShift;};

        /// Num Flops
        /**
          @private
          Compute the number of FLOPS needed to evaluate this Stencil on the box <code>a_box</code>
        */
        inline unsigned long long int numFlops(const Box& a_box) const;
        
        /// Stencil Closed Query
        /**
            Ask if *this has been properly closed.
        */
        inline bool closed() const {return m_isClosed;}

        ///@}

        ///////////////////////////////////////////////////////////////////////////////////////////////
        /** @name Stencil Binding and Application */
        ///@{

        /// Operate on BoxData
        /**
          \ingroup stencil_operations
          Operate *this on a BoxData. This function works in tandem with the namespace-defined
          operators += and |=. See the example in the description of Stencil.

          \param a_scr    Source BoxData
          \param a_scale  (Optional) Scale the output of the Stencil operation (Default: 1)
        */
        template<unsigned int C, MemType MEMTYPE, unsigned char D, unsigned char E>
        inline LazyStencil <T,C,MEMTYPE,D,E>
        operator()(const BoxData<T,C,MEMTYPE,D,E>&  a_src, T a_scale = 1) const;

        /// Operate on BoxData (Overload with Box Input)
        /**
          \ingroup stencil_operations
          Operate *this on a BoxData. This function works in tandem with the namespace-defined
          operators += and |=. See the example in the description of Stencil.

          \param a_scr    Source BoxData
          \param a_box    Confinement Box. Must be a subset of the allowable range of *this
          \param a_scale  (Optional) Scale the output of the Stencil operation (Default: 1)
        */
        template<unsigned int C, MemType MEMTYPE, unsigned char D, unsigned char E>
        inline LazyStencil <T,C,MEMTYPE,D,E>
        operator()(const BoxData<T,C,MEMTYPE,D,E>&  a_src, Box a_box, T a_scale = 1) const;
        ///@}

        ///////////////////////////////////////////////////////////////////////////////////////////////
        /** @name Utility */
        ///@{

        /// Invert Stencil
        /**
          \ingroup stencil_operations
          Inverts the coefficients of this stencil across a given dimension.
        */
        inline void invert(int a_dir);

        /// Transpose Stencil
        /**
          \ingroup stencil_operations
          Transposes *this across two directions in [0,DIM).
          After transpose, coefficients associated with
          the offset (a,...,b,...) will be associated instead with (b,...,a,...)
          */
        inline void transpose(unsigned char a, unsigned char b);

        /// Get Max Index Range
        /**
          @private
          \ingroup stencil_operations
          Given a domain, compute the largest possible iteration Box, taking source/ destination refinement into account.
          For Stencils without source or destination refinement, this function is identical to Stencil::range.
          This function is used primarily for Box inference and probably shouldn't be used publically. 

          \param a_domain     A computational domain
          */
        inline Box indexRange(Box a_domain) const;

        /// Get Max Index Domain
        /**
          @private
          \ingroup stencil_operations
          Given a domain, compute the largest possible iteration Box, taking source / destination refinement into account.
          The output of this function is always a valid input for a Stencil operation.
          For Stencils without source or destination refinement, this function is identical to Stencil::domain.
          This function is used primarily for Box inference and probably shouldn't be used publically. 

          \param a_range      A computational range
          */
        inline Box indexDomain(Box a_range) const;

        /// Get Max Range Box
        /**
          \ingroup stencil_operations
          Given a domain, compute the largest associated physical range, taking refinement AND destination shifting into account.
          The output of this function is NOT a valid input for a Stencil operation when refinement is present.
          This function is best used for defining output BoxData when the input domain is known.

          \param a_domain     A computational domain
          */
        inline Box range(Box a_domain) const;

        /// Get Min Domain Box
        /**
          \ingroup stencil_operations
          Given a range, compute the smallest associated physical domain, taking refinement AND destination shifting into account.
          The output of this function is NOT a valid input for a Stencil operation when refinement is present.
          This function is best used for defining input BoxData when the output range is known.

          \param a_range     A computational range
          */
        inline Box domain(Box a_range) const;

        /// Get Diagonal Value
        /**
          \ingroup stencil_operations
          Returns the coefficient multiplying the identity shift (0,0,...,0).
          If *this doesn't contain the identity shift, returns 0.
        */
        T diagonal() const;

        /// Print
        /**
          \ingroup stencil_operations
          Print *this to the command line. Useful for debugging.
        */
        inline void print() const;
        ///@}

        // End of Stencil Operations Doxygen Module
        /*@}*/

        /** @defgroup stencil_library Stencil Library*/
        /*@{*/

        ///////////////////////////////////////////////////////////////////////////////////////////////
        /** @name Stencil Library */
        ///@{

        /// Stencil Library: Derivative
        /**
          Built in implementation of compact differentiation stencils.
          Includes derivatives of order n >= 1 and accuracy m >= 2 where n + m <= 14.

          \param a_n      Degree of derivative (e.g. nth derivative.)
          \param a_dir    Coordinate of differentiation. Must be in [0,DIM)
          \param a_order  (Optional) Order of accuracy. An unsigned int >= 2. Maximum accuracy depends on a_n. (Default: 2)
          */
        static Stencil<T> Derivative(int a_n, int a_dir, int a_order = 2);

        /// Stencil Library: Laplacian
        /**
          Built in implementation of the 2nd order 2*DIM + 1 point Laplace operator.
          */
        static Stencil<T> Laplacian();
#if DIM == 2
        /// 9 Point Laplacian (Mehrstellen)
        static Stencil<T> Laplacian_9();
#elif DIM == 3
        /// 19 Point Laplacian (Mehrstellen)
        static Stencil<T> Laplacian_19();
        /// 27 Point Laplacian
        static Stencil<T> Laplacian_27();
#endif
        /// Stencil Library: Perpendicular Laplacian
        /**
          Built in implementation of Laplacian perpendicular to direction dir

          \param a_dir    Normal direction
          \param a_order  (Optional) Order of accuracy (default: 2 | supported: 2)
          */
        static Stencil<T> LaplacianFace(int a_dir, int a_order = 2);

        /// Stencil Library: Cell to Face Interpolation
        /**
          Interpolates cell averaged values to face averaged values.
          The <code>a_side</code> controls if the face averaged output is on the 
          upper or lower face of the cell at the center of the stencil. 
          For finer control of the Stencil's centering, see the functions
          <code>CellToFaceL</code> and <code>CelltoFaceH</code>. 

          \param a_dir      Coordinate direction in [0, DIM)
          \param a_side     Upper or lower side (Side::Hi or Side::Lo) (default: Side::Lo)
          \param a_order    (Optional) Order of accuracy. Valid values (default 4 | supported: 4, 5)
        */
        static Stencil<T> CellToFace(
                int             a_dir,
                Side::LoHiSide  a_side = Side::Lo,
                int             a_order = 4);

        /// Stencil Library: Cell to Face Differentiation
        /**
          Computes the normal derivative of a cell-averaged quantity as a face-averaged quantity.

          \param dir      Coordinate direction in [0, DIM)
          \param side     Upper or lower side (Side::Hi or Side::Lo) (default: Side::Lo)
          \param order    (Optional) Order of accuracy. Valid values (default 4 | supported: 4, 5)
        */
        static Stencil<T> DiffCellToFace(
                int             a_dir,
                Side::LoHiSide  a_side = Side::Lo,
                int             a_order = 4);

        /// Stencil Library: Upwind Cell to Face Interpolation
        /**
          Interpolates cell averaged values to face averaged values.
          The <code>a_side</code> controls if the face averaged output is on the 
          upper or lower face of the cell at the center of the stencil. 
          For even orders, this function is identical to CellToFace.
          For odd orders, the stencil is upwinded

          Example Stencil Footprints:
          @code
          // Legend | o: center stencil weight | ^: output location | x: other stencil weights
          auto S4L = CellToFaceL(0, Side::Lo, 4);
          // Footprint: 
          //    | x | x | o | x |
          //            ^
          auto S4H = CellToFaceL(0, Side::Hi, 4);
          // Footprint: 
          //    | x | o | x | x |
          //            ^
          auto S5L = CellToFaceL(0, Side::Lo, 5);
          // Footprint: 
          //    | x | x | x | o | x |
          //                ^
          auto S5H = CellToFaceL(0, Side::Hi, 5);
          // Footprint: 
          //    | x | x | o | x | x |
          //                ^
          @endcode

          \param dir      Coordinate direction in [0, DIM)
          \param side     Upper or lower side (Side::Hi or Side::Lo) (default: Side::Lo)
          \param order    (Optional) Order of accuracy. ( default: 4 | supported: 4, 5 )
        */
        static Stencil<T> CellToFaceL(
                int             a_dir,
                Side::LoHiSide  a_side = Side::Lo,
                int             a_order = 4);

        /// Stencil Library: Upwind Cell to Face Interpolation
        /**
          Interpolates cell averaged values to face averaged values.
          The <code>a_side</code> controls if the face averaged output is on the 
          upper or lower face of the cell at the center of the stencil. 
          For even orders, this function is identical to CellToFace.
          For odd orders, the Stencil is downwinded.
          
          Example Stencil Footprints:
          @code
          // Legend | o: center stencil weight | ^: output location | x: other stencil weights
          auto S4L = CellToFaceH(0, Side::Lo, 4);
          // Footprint: 
          //    | x | x | o | x |
          //            ^
          auto S4H = CellToFaceH(0, Side::Hi, 4);
          // Footprint: 
          //    | x | o | x | x |
          //            ^
          auto S5L = CellToFaceH(0, Side::Lo, 5);
          // Footprint: 
          //    | x | x | o | x | x |
          //            ^     
          auto S5H = CellToFaceH(0, Side::Hi, 5);
          // Footprint: 
          //    | x | o | x | x | x |
          //            ^    
          @endcode

          \param dir      Coordinate direction in [0, DIM)
          \param side     Upper or lower side (Side::Hi or Side::Lo ) (default: Side::Lo)
          \param order    (Optional) Order of accuracy. ( default: 4 | supported: 4, 5 )
        */
        static Stencil<T> CellToFaceH(
                int             a_dir,
                Side::LoHiSide  a_side = Side::Lo,
                int             a_order = 4);

        /// Stencil Library: Simple Average
        /**
          Averages data from a refined grid onto a coarsened grid.
          Refinement is assumed to be isotropic.
          Source data is refined relative to the destination data by <code>a_refRatio</code>

          \param a_refRatio   Refinement ratio.
          */
        static Stencil<T> AvgDown(int a_refRatio);

        /// Stencil Library: Anisotropic Average
        /**
          Anisotropic overload of <code>AvgDown</code>.

          \param a_refRatio   Refinement ratio
        */
        static Stencil<T> AvgDown(Point a_refRatio);

        /// Stencil Library: Sum
        /**
          Undivided average (e.g. the sum of all elements).
          Source data is refined relative to the destination data by <code>a_refRatio</code>

          \param a_refRatio   Refinement ratio 
        */
        static Stencil<T> Sum(int a_refRatio);

        /// Stencil Library: Sum
        /**
          Anisotropic overload of <code>Sum</code>.

          \param a_refRatio     Refinement ratio
          */
        static Stencil<T> Sum(Point a_refRatio);

        //TODO: I don't think this does the right thing for Side::Hi. Needs better description.
        /// Stencil Library: Simple Average over a Face
        /**
          Averaging operation between face-averaged data sets with normal <code>a_normDir</code>.

          \param a_normDir      Normal coordinate direction
          \param a_side         Cell side to average over
          \param a_refRatio     Non-normal coordinate refinement ratios
        */
        static Stencil<T> AvgDownFace(int a_normDir, Side::LoHiSide a_side, int a_refRatio);

        //TODO: I don't think this does the right thing for Side::Hi. Needs better description.
        /// Stencil Library: Simple Average over a Face
        /**
          Averaging operation between face-averaged data sets with normal <code>a_normDir</code>.

          \param a_normDir      Normal coordinate direction
          \param a_side         Cell side to average over
          \param a_refRatio     Non-normal coordinate refinement ratios
        */
        static Stencil<T> AvgDownFace(int a_dir, Side::LoHiSide a_side, Point a_refRatio);

        //TODO: Give this a Side::LoHiSide argument so it can handle high-sided fluxes
        /// Stencil Library: Flux Divergence
        /**
          Simple flux differencing stencil: OUT(i) = IN(i+1) - IN(i)
          Assumes the low-side fluxes are stored in cell i.

          \param a_dir    Coordinate axis in the flux normal direction.
        */
        static Stencil<T> FluxDivergence(int a_dir);
        ///@}

        // End of Stencil Operations Doxygen Module
        /*@}*/

        ///////////////////////////////////////////////////////////////////////////////////////////////
        //TODO: Functional syntax convension is (output, inputs...)
        /// Apply Stencil Helper function
        /**
          \ingroup stencil_operations
          Manually apply *this to a source and destination BoxData in a functional programming style.
          Not recommended for public use, but useful for debugging in some cases.

          \param a_src        Source data
          \param a_dest       Output data
          \param a_bx         Iteration box for this computation. Can be created with indexDomain or indexRange functions.
          \param a_replace    (Optional) If true, zero out the data in a_dest and replace with the computation result. If true, solution will be added to existing data. Default false.
          \param a_scale      (Optional) Scale the computation by some value.
          */
        template<unsigned int C, MemType MEMTYPE, unsigned char D, unsigned char E>
        void apply(
                const BoxData<T,C,MEMTYPE,D,E>& a_src,
                BoxData<T,C,MEMTYPE,D,E>&       a_dest,
                const Box&                      a_bx,
                bool                            a_initToZero = false,
                const T                         a_scale = 1) const;

        //TODO: Documentation 
        template<unsigned int C, unsigned char D, unsigned char E>
        void protoApply(
                const BoxData<T,C,HOST,D,E>& a_src,
                BoxData<T,C,HOST,D,E>&       a_dst,
                const Box&                      a_box,
                bool                            a_initToZero = false,
                T                               a_scale = 1);
        
        //TODO: Documentation 
        template<unsigned int C, unsigned char D, unsigned char E>
        void protoApply(
                const BoxData<T,C,DEVICE,D,E>& a_src,
                BoxData<T,C,DEVICE,D,E>&       a_dst,
                const Box&                      a_box,
                bool                            a_initToZero = false,
                T                               a_scale = 1);

        //TODO: Documentation
        template<unsigned int C, unsigned char D, unsigned char E>   
            void hostApply(  const BoxData<T,C,MemType::HOST,D,E>&  a_src,
                    BoxData<T,C,MemType::HOST,D,E>&  a_dst,
                    const Box&               a_box,
                    bool               a_initToZero = false,
                    T                  a_scale = 1) const;
        
    private:
        /// Add Coefficient-Offset Pair
        /**
          @private
          Helper function that encapsulates all of the proper checks needed
          to add a coefficient-offset pair to a Stencil. Any new function
          that needs to add data to an existing Stencil should call this.
          */
        void addCoef(T a_coef, Point a_offset);
       
        // Helper object for device Stencils
        /* Appears to be unused. -CLG 4/19/2022
        struct coeff_holder
        {
            T elements[64];
        };
        */
        //Must be called before computations can be performed on the device
        inline void closeForDevice();

        std::vector<T>      m_coefs;        ///< Coefficients of the Stencil.
        std::vector<Point>  m_offsets;      ///< Offsets associated with the Stencil.
        Point               m_srcRefratio;  ///< Refinement of source data
        Point               m_destRefratio; ///< Refinement of destination data
        Point               m_destShift;    ///< Output shift in (refined) destination data. Meaningless without destination refinement
        Box                  m_span;         ///< Bounding Box defining the largest offsets of *this. Defines Stencil's spatial footprint

        bool m_isClosed;
#ifdef PROTO_CUDA
        //this copies to the device
        T*    d_coeff;
        Point* d_offset;
#endif
        // Appears to be unused. -CLG 4/19/2022
        //Stencil<T>::coeff_holder c_coeff;
    }; // end class Stencil

///////////////////////////////////////////////////////////////////////////////////////////////
/** @name Non-Member Functions */
///@{

/// Coefficient Shift Product "Constructor"
/**
    \ingroup stencil_operations
    Syntactical sugar that allows a Stencil to be constructed with coef*Shift(Point) syntax

    Example:
    @code
    Stencil<double> S = 3.7*Shift(Point::Basis(0));
    @endcode
*/
template <typename T>
inline Stencil<T> operator*(T a_coef, Shift a_shift) { return Stencil<T>(a_shift, a_coef); }

/// Scalar Multiplication of Stencil Coefficients
/**
    \ingroup stencil_operations
    Allows for pre multiplication by a T scalar
*/
template <typename T>
inline Stencil<T> operator*(T a_coef, const Stencil<T> a_stencil) { return a_stencil*a_coef; }

/// Stencil Unary Negation
/**
    \ingroup stencil_operations
*/
template <typename T>
inline Stencil<T> operator-(Stencil<T> a_stencil) { return a_stencil*(-1); }

/// Application by Replacement
/**
    \ingroup stencil_operations
    Applies a Stencil and REPLACES a subset of a_dest with the result:

    Usage:
    @code
    Stencil<double> L = Stencil<double>::Laplacian(); //2*DIM+1 point Laplacian operator
    Box rangeBox = Box::Cube(64);                     //define the data range
    Box domainBox = L.domain(rangeBox);               //compute the domain from the range
    BoxData<double> Src(domainBox);
    BoxData<double> Dst(rangeBox);
    // ... initialize Src ...
    Dst |= L(Src);
    // NB:  This call computes the Laplacian of Src and replaces the values of Dst
    //      inside rangeBox with the result. Only the cells in rangeBox & Dst.box()
    //      are affected. If rangeBox & Dst.box() is empty, |= is a null-op.
    @endcode

    TODO: These examples are missing.
    See the main documentation for Stencil for additional examples

    \param a_dest   Destination array
    \param a_op     Uncomputed Stencil operation (output of Stencil::operator())
*/
template <typename T, unsigned int C, MemType MEMTYPE, unsigned char D, unsigned char E>
inline BoxData<T,C,MEMTYPE,D,E>& operator|=(
        BoxData<T,C,MEMTYPE,D,E>&       a_dest,
        LazyStencil<T,C,MEMTYPE,D,E>&&  a_op);

/// Application by Increment
/**
    \ingroup stencil_operations
    Applies a Stencil and ADDS a subset of a_dest with the result

    Usage:
    @code
    Stencil<double> L = Stencil<double>::Laplacian(); //2*DIM+1 point Laplacian operator
    Box rangeBox = Box::Cube(64);                     //define the data range
    Box domainBox = L.domain(rangeBox);               //compute the domain from the range
    BoxData<double> Src(domainBox);
    BoxData<double> Dst(rangeBox);
    // ... initialize Src ...
    Dst += L(Src);
    // NB:  This call computes the Laplacian of Src and addds the values of Dst
    //      inside rangeBox with the result. Only the cells in rangeBox & Dst.box()
    //      are affected. If rangeBox & Dst.box() is empty, |= is a null-op.
    @endcode

    TODO: These examples are missing.
    See the main documentation for Stencil for additional examples

    \param a_dest   Destination array
    \param a_op     Uncomputed Stencil operation (output of Stencil::operator())
*/
template <class T, unsigned int C, MemType MEMTYPE, unsigned char D, unsigned char E>
BoxData<T,C,MEMTYPE,D,E>& operator+=(BoxData<T,C,MEMTYPE,D,E>& a_dest, LazyStencil<T,C,MEMTYPE,D,E>&& a_op);

///@}
#include "implem/Proto_StencilImplem.H"
#include "implem/Proto_StencilDefs.H"
} //end Proto namespace
#endif
