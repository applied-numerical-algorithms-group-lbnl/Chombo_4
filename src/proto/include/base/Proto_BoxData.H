#pragma once
#ifndef _PROTO_RECTMDARRAY_H_
#define _PROTO_RECTMDARRAY_H_
#include "Proto_Timer.H"

#include "Proto_Memory.H"
#include "Proto_MemType.H"

#include "Proto_CInterval.H"
#include "Proto_Box.H"
#include "Proto_PowerItoI.H"
#include "Proto_Stencil.H"
#include "Proto_Reduction.H"
#include "Proto_Stack.H"
#include "Proto_DisjointBoxLayout.H"

#include <memory>   //for shared_ptr
#include <iostream>
#include <fstream>  // for parallel pretty printing
#include <iomanip>  // for pretty printing
#include <limits>   // for max / min functions
#include <cstdlib>  // for size_t
#include <cmath>    // floating point abs
#include <map>
#include <unordered_map>
 
// This is to facilitate automated testing.
#ifdef PROTO_MEM_CHECK
struct BoxDataMemCheck
{
    static BoxDataMemCheck& instance()
    {
        static BoxDataMemCheck inst;
        return inst;
    }
    
    BoxDataMemCheck()
    {
        m_numCopies = 0;
        m_numSlowCopies = 0;
        m_numFastCopies = 0;
        m_numMoveAssign = 0;
    }
    
    inline static int& numCopies()     { return instance().m_numCopies; }
    inline static int& numFastCopies() { return instance().m_numFastCopies; }
    inline static int& numSlowCopies() { return instance().m_numSlowCopies; }
    inline static int& numMoveAssign()  { return instance().m_numMoveAssign; }

    inline static void clear()
    {
        numCopies() = 0;
        numFastCopies() = 0;
        numSlowCopies() = 0;
        numMoveAssign() = 0;
    }
    
    private:

    int m_numCopies;
    int m_numSlowCopies;
    int m_numFastCopies;
    int m_numMoveAssign;
};
#endif

using std::shared_ptr;
namespace Proto 
{
    // Forward Declarations
    template<typename T> class Stencil;
    template<typename T, unsigned int C, MemType MEM, unsigned char D, unsigned char E>
    class LazyStencil;

    template<MemType MEM>
    inline void null_deleter_boxdata(void* ptr) { }

    enum BoxDataOp
    {
        Invalid = -1,
        Copy    = 0,
        Add     = 1,
        Assign,
        Subtract, 
        Multiply, 
        Divide, 
        //ThresholdMax,
        //ThresholdMin,
        NUM_OP_TYPES
    };

/** @defgroup pointwise_operations Pointwise Operations*/

/** @addtogroup pointwise_operations*/
/*@{*/

    /// Pointwise Variable
    /**
      Var is used for implementing pointwise operators on BoxDatas.
      See the documentation for forall for more information.
      Template parameters refer to the type and dimension of data held by the associated BoxData.
      See BoxData documentation.

      \tparam T Type of data in array (bool, int, double, etc.)
      \tparam C Number of components in first data index. Defaults to 1. 
      \tparam D Number of components in second data index. Defaults to 1. 
      \tparam E Number of components in third data index. Defaults to 1.
      */
    template<typename T, unsigned int C=1, MemType MEM=MEMTYPE_DEFAULT,
        unsigned char D=1, unsigned char E=1>
    class Var
    {
        public:


#ifdef PROTO_CUDA //__CUDA_ARCH__
        __device__
        inline T& getValDevice(unsigned int a_c, unsigned char a_d = 0, unsigned char a_e = 0)
        {
            int idx = threadIdx.x + blockIdx.x*blockDim.x;
            int idy = blockIdx.y;
    #if DIM == 3
            int idz = blockIdx.z;
            return (m_ptrs[a_c + C*a_d + C*D*a_e])[idx + boxDimX * (idy + idz * boxDimY)];
    #else
            return (m_ptrs[a_c + C*a_d + C*D*a_e])[idx + idy * boxDimX];
    #endif
        }

    #ifdef PROTO_HIP
        __host__
        inline T& getValDevice(unsigned int a_c, unsigned char a_d = 0, unsigned char a_e = 0)
        {
            return *(m_ptrs[a_c + C*a_d + C*D*a_e]);
        }
    #endif
#else
        inline T& getValDevice(unsigned int a_c, unsigned char a_d = 0, unsigned char a_e = 0)
        {
            return *(m_ptrs[a_c + C*a_d + C*D*a_e]);
        }
#endif
        /// Pointwise Accessor
        /**
          Access component (c,d,e) of the <code>BoxData<T,C,MEM,D,E></code> 
          associated with *this.

          \param a_c   First component index
          \param a_d   Second component index (default: 0)
          \param a_e   Third component index  (default: 0)
        */
        CUDA_DECORATION
        __attribute__((always_inline))
        T& operator()(unsigned int a_c, unsigned char a_d = 0, unsigned char a_e = 0)  
        {
#ifdef PROTO_CUDA
            if(MEM==MemType::DEVICE)
                return getValDevice(a_c,a_d,a_e);
            else
#endif
                return *(m_ptrs[a_c + C*a_d + C*D*a_e]);
        }

        /// Pointwise Accessor (Const)
        /**
          Access component (c,d,e) of the <code>const BoxData<T,C,D,E></code> 
          associated with <code>*this</code>.

          \param a_c   First component index
          \param a_d   Second component index (default: 0)
          \param a_e   Third component index  (default: 0)
          */
#ifdef PROTO_CUDA //__CUDA_ARCH__
        __device__
        inline const T& getValDevice(
                unsigned int a_c,
                unsigned char a_d = 0,
                unsigned char a_e = 0) const
        {
            int idx = threadIdx.x + blockIdx.x*blockDim.x;
            int idy = blockIdx.y;
#if DIM == 3
            int idz = blockIdx.z;
            return (m_ptrs[a_c + C*a_d + C*D*a_e])[idx + boxDimX * (idy + idz * boxDimY)];
#else
            return (m_ptrs[a_c + C*a_d + C*D*a_e])[idx + idy * boxDimX];
#endif
        }

#ifdef PROTO_HIP
        __host__
        inline const T& getValDevice(
                unsigned int a_c,
                unsigned char a_d = 0,
                unsigned char a_e = 0) const
        {
            return *(m_ptrs[a_c + C*a_d + C*D*a_e]);
        }
#endif
#else
        inline const T& getValDevice(
                unsigned int a_c,
                unsigned char a_d = 0,
                unsigned char a_e = 0) const
        {
            return *(m_ptrs[a_c + C*a_d + C*D*a_e]);
        }
#endif

        CUDA_DECORATION
        __attribute__((always_inline))
        const T& operator()(
                unsigned int a_c,
                unsigned char a_d = 0,
                unsigned char a_e = 0) const 
        {
#ifdef PROTO_CUDA
            if(MEM==MemType::DEVICE)
                return getValDevice(a_c,a_d,a_e);
            else
#endif
                return *(m_ptrs[a_c + C*a_d + C*D*a_e]);
        }

        CUDA_DECORATION
        inline Var&  operator+=(unsigned int& a_increment) 
        {
            for (int ii = 0; ii < C*D*E; ii++)
            {
                m_ptrs[ii] += a_increment;
            }
            return *this;
        }

        // test
        CUDA_DECORATION
        inline Var&  operator+=(Point& a_p)
        {
#if DIM == 3
            unsigned int shift = a_p[0] + (a_p[1] + a_p[2]*boxDimY)*boxDimX;
#else
            unsigned int shift = a_p[0] + a_p[1]*boxDimX;
#endif
            return *this+=shift;
        }

        CUDA_DECORATION
        inline Var&  operator+=(const Point& a_p)
        {
#if DIM == 3
            unsigned int shift = a_p[0] + (a_p[1] + a_p[2]*boxDimY)*boxDimX;
#else
            unsigned int shift = a_p[0] + a_p[1]*boxDimX;
#endif
            return *this+=shift;
        }

        CUDA_DECORATION
        inline Var& operator++() 
        {
            for (int ii = 0; ii < C*D*E; ii++)
            {
                ++m_ptrs[ii];
            }
            return *this;
        }

        CUDA_DECORATION
        inline Var& operator--() 
        {
            for (int ii = 0; ii < C*D*E; ii++)
            {
                --m_ptrs[ii];
            }
            return *this;
        }

        unsigned int boxDimX;
        unsigned int boxDimY;
        unsigned int subBoxDimX;
        unsigned int subBoxDimY;
        T* m_ptrs[C*D*E];
    }; // End Class Var

    ///////////////////////////////////////////////////////////////////////////////////////////////
    /// Multidimensional Rectangular Array
    /**
      BoxData is the main dataholder class of Proto. It contains an 
      array of data of type T on a domain defined by a Box.
      Each data point is associated with a Point in the Box.
      C, D, and E define the tensor structure of the data:
      C,D,E = 1:            Data is scalar (default)
      C = N, D,E = 1:       Data is vector valued with N components
      C = M, D = N, E = 1:  Data is MxN matrix valued

      \tparam T Type of data in array (bool, int, double, etc.)
      \tparam C (Optional) Size of first component axis. Defaults to 1. 
      \tparam D (Optional) Size of second component axis. Defaults to 1. 
      \tparam E (Optional) Size of third component axis. Defaults to 1.
      \tparam MEM (Optional) Where memory is allocated (device or host). Defaults to host.
      */
    template <class T=double, unsigned int C=1,
    MemType MEM=MEMTYPE_DEFAULT,
    unsigned char D=1, unsigned char E=1>
    class BoxData 
    {
        public:

        typedef Var<T,C,MEM,D,E> reference;
        
        /// Forbidden Functions
        BoxData(const BoxData<T,C,MEM,D,E>& a_src) = delete;
        BoxData& operator=(const BoxData<T,C,MEM,D,E>& a_src) = delete;
       

        //////////////////////////////////////////////////////////////////////////////////////
        /** @name Constructors and Define*/
        ///@{


        /// Default Constructor
        /**
          Creates an empty BoxData with an empty Box.
          <code>BoxData::define</code> must be called to allocate memory.
        */
        BoxData();

        /// Box Constructor
        /**
          Creates and allocates an uninitialized array according to the space 
          needed by <code>a_box</code> and the number of components.

          \param a_box      Box defining the domain of *this
        */
        explicit BoxData(const Box& a_box);

        /// Box Constructor With Initialization
        /**
          Creates, allocates, and initializes a BoxData using a constant value. 
          The Stack is automatically used if it is active when this function is called.

          \param a_box      Box defining the domain of *this
          \param a_init     T value used to initialize all data values
        */
        BoxData(const Box& a_box, T a_init);

        /// Move constructor
        /**
          Necessary in cases where we return BoxData by value, but don't want an  actual deep copy
          The Stack is automatically used if it is active when this function is called.

          \param a_src    Source data
        */
        BoxData(BoxData<T,C,MEM,D,E>&& a_src);

        /// Box Define
        /**
          Allocates or reallocates memory for a previously declared BoxData.
          The Stack is automatically used if it is active when this function is called.

          \param a_box              Box defining the domain of *this
        */  
        void define(const Box& a_box);


        /// Alias Define
        /**@private */
        template<unsigned int ncomp>
        void define(BoxData<T,ncomp,MEM,1,1>& a_input,
                    unsigned int&                 a_comp);

        /// Raw Pointer Alias Constructor
        /**
          Builds a vector valued BoxData by aliasing a raw pointer to T which should be of size
          <code> a_box.size()*a_ncomp </code>. The data in <code> a_ptr </code> is not copied.
          This constructor is not recommended for public use, but is provided to facilitate
          compatability with third party libraries.

          \param a_ptr   MEM Source buffer of type T to which *this will be aliased.
          \param a_box   Box defining the domain of *this
          \param a_ncomp Dummy input for the number of components. Completely unused.
        */
        BoxData(const T* a_ptr, const Box& a_box, int a_ncomp = C);

        /// Raw Pointer Alias Define
        /**
          Defines a vector valued BoxData by aliasing a raw pointer to T which should be of size
          <code> a_box.size()*a_ncomp </code>. The data in <code> a_ptr </code> is not copied.
          This constructor is not recommended for public use, but is provided to facilitate
          compatability with third party libraries.

          \param a_ptr   MEM Source buffer of type T to which *this will be aliased.
          \param a_box   Box defining the domain of *this
          \param a_ncomp Dummy input for the number of components. Completely unused.
        */
        void define(const T* a_ptr, const Box& a_box, int a_ncomp = C);


        /// LazyStencil Constructor
        /**
          @private
          Builds a BoxData from a lazily evaluated Stencil operation. Called implicitly.
          Allows for the syntax:

          @code
          Box domain = Box::Cube(64);
          BoxData<double> ones(domain, 1);
          Stencil<double> triple = 3*Shift::Zeros();
          BoxData<double> threes = triple(ones);
          //NB: Using "auto" here will result in an error:
          auto threes = triple(ones); // runtime error
          @endcode

          \param a_op   The output of a Stencil or InterpStencil apply operation.
          */
        BoxData(LazyStencil<T,C,MEM,D,E>&& a_op);

        /// Slice Constructor
        /**
          @private
          Used internally. Not part of the public interface.
        */
        BoxData(shared_ptr<T> a_data,const T* a_ptr,const Box& a_box);

        /// Destructor.
        ~BoxData(); 

        ///@}

        ///////////////////////////////////////////////////////////////////////////////////////
        /** @name Data Movement */
        ///@{

        /// Move Assignment Operator
        /**
          Moves data in a_src to *this without performing a deep copy.

          \param a_src   Source Data 
        */
        BoxData& operator=(BoxData<T,C,MEM,D,E>&& a_src);

        /// Copy on Intersection.
        /**
          Copy data from <code> *this </code> to <code> a_dest </code> 
          within the domain intersection. Does nothing if the domain intersection is empty.
          Note that <code>a_dest</code> must have the same MemType as *this.

          \param a_dest  Destination data holder
          */
        template<MemType MEM_DEST>
        void copyTo(BoxData<T,C,MEM_DEST,D,E>& a_dest) const;

        /// Copy Region
        /**
          Copy with a prescribed Box argument and optional shift.
          Explicitly, this function copies the subset of data from <code> *this </code> 
          contained in <code> a_srcBox </code> into the region of <code> a_dest </code> 
          defined by <code> a_srcBox.shift(a_destShift) </code>.

          This function will fail if <code> a_srcBox </code> is not contained in both
          <code> this->box() </code> and <code> a_dest.box().shift(a_destShift) </code>.

          Note that <code>a_dest</code> must have the same MemType as *this.

          \param a_dest       Destination data holder
          \param a_srcBox     Region of data to copy from *this
          \param a_destShift  (Optional) Determines region of a_dest to copy data to.
          */
        template<MemType MEM_DEST>
        void copyTo(
                BoxData<T,C,MEM_DEST,D,E>&   a_dest,
                const Box&                  a_srcBox,
                const Point&                a_destShift = Point::Zeros()) const;

        ///Copy with src and dest boxes.                                          
        /**
          Copies from a <code>a_srcBox</code> in <code>*this</code> into 
          a <code>a_destBox</code> in <code>a_dest</code>. Used in ::Proto::Copier.
          Note that <code>a_dest</code> must have the same MemType as *this.
          TODO: This needs debugging. It's unclear what this function is supposed to do.
        */
        template<MemType MEM_DEST>
        void copyTo(BoxData<T,C,MEM_DEST,D,E>& a_dest,
                const Box&    a_srcBox,
                const Box&    a_destBox) const;

        /// General Copy
        /**
          The most general form of copyTo.
          Copies a prescribed set of components from <code> *this </code>
          into <code> a_dest </code> within prescribed region with a possible shift.
          <code>a_srcComps </code> and <code> a_destComps </code> must have the same size.
          The range of the copy in <code>a_dest</code> is 
          <code>a_srcBox.shift(a_destShift) & a_dest.box()</code>. 
          Note that <code>a_dest</code> must have the same MemType as *this.

          Example Usage:

          \param a_dest       Destination data holder
          \param a_srcBox     Region of data to copy from *this
          \param a_srcComps   Components of *this to copy from    
          \param a_destShift  Determines region of a_dest to copy data to
          \param a_destComps  Components of a_dest to copy into.
        */
        template <unsigned int Cdest, MemType MEM_DEST,
                 unsigned char Ddest, unsigned char Edest>
        void copyTo(
                BoxData<T,Cdest,MEM_DEST,Ddest,Edest>& a_dest,
                const Box&      a_srcBox,
                CInterval       a_srcComps,
                const Point&    a_destShift,
                CInterval       a_destComps) const;

        /// General Copy From
        /**
          Nearly identical to the most general version of <code> copyTo </code> 
          except the copy direction is reversed.
          Only valid for vector or scalar valued BoxData.
          This function is provided to facilitate interoperability with third party libraries, 
          and is not recommended for public use unless necessary.

          \param a_dsrc     Source data
          \param a_srcBox   Source domain to copy from
          \param a_srcComp  First component to copy from source
          \param a_destBox  Destination domain to copy into
          \param a_destComp First component to copy into destiation
          \param a_numcomp  Number of components to copy
          */
        template<unsigned int Csrc>
        void copy(
                const BoxData<T,Csrc,MEM,D,E>& a_dsrc,
                const Box&      a_srcBox,
                unsigned int    a_srcComp,
                const Box&      a_destBox,
                unsigned int    a_destComp,
                unsigned int    a_numcomp);

        ///@}

        ///////////////////////////////////////////////////////////////////////////////////////
        /** @name Accessors */
        ///@{

        /// Compute Index
        /**
          Computes the index in <code>[0,this->size())</code> associated with a given
          Point and components

          \param a_pt   A Point in the domain of <code> *this </code>
          \param a_c  (Optional) First tensor index.  (defaults : 0)
          \param a_d  (Optional) Second tensor index. (defaults : 0)
          \param a_e  (Optional) Third tensor index.  (defaults : 0)
          */
        CUDA_DECORATION
        inline size_t index(const Point a_pt,
                unsigned int a_c = 0,
                unsigned int a_d = 0,
                unsigned int a_e = 0) const;
        
        //TODO: Why do we need so many accessors? Sury we can prune these down
        //TODO: Documentation claims this is read only, but it returns a writeable value?
        /// Point Accessor 
        /**
          Read-only access to data stored in <code> *this </code>.
          This function is read-only to facilitate cross-platform code (e.g. GPU compatability).
          This function is for debugging ONLY, and will not work on the device.
          Pointwise write operations should be done through forall.

          \param a_pt   A Point in the domain of <code> *this </code>
          \param a_c  (Optional) First tensor index.  (defaults : 0)
          \param a_d  (Optional) Second tensor index. (defaults : 0)
          \param a_e  (Optional) Third tensor index.  (defaults : 0)
          */
        CUDA_DECORATION
        inline T& operator()(
                const Point& a_pt,
                unsigned int  a_c = 0,
                unsigned char a_d = 0,
                unsigned char a_e = 0) const;

        /// Index Accessor (Const)
        /**
          Access *this using a linear index. Useful for interacting with
          BoxData like a regular buffer, but not recommended.

          \param a_index   Index in [0,this->size() = m_box.size()*C*D*E)
          */
        inline const T* operator[](unsigned int a_index) const; 

        /// Index Accessor (Non-Const)
        /**
          Access *this using a linear index. Useful for interacting with
          BoxData like a regular buffer, but not recommended.

          \param a_index   Index in [0,this->size() = m_box.size()*C*D*E)
          */
        inline T* operator[](unsigned int a_index); 

        /// Pointer Accessor (Const)
        /**
          Return a const pointer to a data point in <code> *this </code>
          given a Point and tensor indices.

          \param a_p  Point in the domain of this 
          \param a_c  (Optional) First tensor index.  (defaults : 0)
          \param a_d  (Optional) Second tensor index. (defaults : 0)
          \param a_e  (Optional) Third tensor index.  (defaults : 0)
          */
        inline const T* data(
                const Point&  a_p,
                unsigned int  a_c = 0,
                unsigned char a_d = 0,
                unsigned char a_e = 0) const;
            
        /// Pointer Accessor (Non-Const)
        /**
          Return a pointer to a data point in <code> *this </code>
          given a Point and tensor indices.

          \param a_p  Point in the domain of this 
          \param a_c  (Optional) First tensor index.  (defaults : 0)
          \param a_d  (Optional) Second tensor index. (defaults : 0)
          \param a_e  (Optional) Third tensor index.  (defaults : 0)
          */
        inline T* data(
                const Point&  a_p,
                unsigned int  a_c = 0,
                unsigned char a_d = 0,
                unsigned char a_e = 0);

        /// Buffer Accessor (Non-Const)
        /** 
          Return the contiguous data buffer where the data in <code>*this</code> is stored,
          starting at a specified component. Not recommended for public use.
          
          \param a_c  (Optional) First tensor index.  (defaults : 0)
          \param a_d  (Optional) Second tensor index. (defaults : 0)
          \param a_e  (Optional) Third tensor index.  (defaults : 0)
        */
        inline T* data(
                unsigned int  a_c = 0,
                unsigned char a_d = 0,
                unsigned char a_e = 0);

        /// Buffer Accessor (Const)
        /** 
          Return the contiguous data buffer where the data in <code>*this</code> is stored.
          Not recommended for public use. 
          */
        inline const T* data(
                unsigned int  a_c = 0,
                unsigned char a_d = 0,
                unsigned char a_e = 0) const;

        /// Shared Buffer Accessor
        /**
            Returns the shared_ptr object that is common to all aliases. If this not 
            an alias, this function is identical to data().
        */
        inline ::std::shared_ptr<T>  aliasData() { return m_data; }
        
        /// Shared Buffer Accessor (Const)
        /**
            Returns the shared_ptr object that is common to all aliases. If this not 
            an alias, this function is identical to data().
        */
        inline ::std::shared_ptr<T>  aliasData() const { return m_data; }
        
        // Get Domain Box
        /**
          Return's the domain Box of <code> *this </code>. 
          Includes ghost cells if it was built it that way.
          */
        inline Box box() const {return m_box;};

        // Get Size
        /**
          Returns the number of data points in <code> *this </code>.
          Return value is equal to <code> m_box.size()*C*D*E </code>
          */

        inline std::size_t size() const {return m_box.size()*C*D*E;};

        /// Defined Query
        /**
          Returns true if this is defined (e.g., it has memory allocated to it).
          This will return false for default constructed objects and true otherwise
        */
        inline bool defined() const {return bool(m_data);};

        //TODO: should be private
        /// Create Pointwise Access Variable (Non-Const)
        /**
          @private
          */
        inline Var<T,C,MEM,D,E> var(const Point& a_pt);

        //TODO: should be private
        /// Create Pointwise Access Variable (Const)
        /**
          @private
          */
        inline Var<T,C,MEM,D,E> var(const Point& a_pt) const;

        ///@}

        ///////////////////////////////////////////////////////////////////////////////////////
        /** @name Algebraic Operations */
        ///@{

        /// Algebraic Op Helper
        /** @private */
        template<BoxDataOp>
        inline void operatorT(const BoxData<T,C,MEM,D,E>& a_src);

        /// Algebraic Op Helper
        /** @private */
        template<BoxDataOp>
        inline void operatorT(const T a_scale);

        /// Pointwise Addition on Intersection
        /**
          Not recommended for performance unless necessary.
          Try to fuse arithmetic operations into forall or Stencil operations if possible.

          \param a_rhs    Another BoxData
          */
        inline BoxData<T,C,MEM,D,E>& operator+=(const BoxData<T,C,MEM,D,E>& a_rhs);

        /// Pointwise Subtraction on Intersection
        /**
          Not recommended for performance unless necessary.
          Try to fuse arithmetic operations into forall or Stencil operations if possible.

          \param a_rhs    Another BoxData
          */
        inline BoxData<T,C,MEM,D,E>& operator-=(const BoxData<T,C,MEM,D,E>& a_rhs);

        /// Pointwise Multiplication on Intersection
        /**
          Not recommended for performance unless necessary.
          Try to fuse arithmetic operations into forall or Stencil operations if possible.

          \param a_rhs    Another BoxData
          */
        inline BoxData<T,C,MEM,D,E>& operator*=(const BoxData<T,C,MEM,D,E>& a_rhs);

        /// Pointwise Division on Intersection
        /**
          Not recommended for performance unless necessary.
          Try to fuse arithmetic operations into forall or Stencil operations if possible.

          \param a_rhs    Another BoxData
          */
        inline BoxData<T,C,MEM,D,E>& operator/=(const BoxData<T,C,MEM,D,E>& a_rhs);

        /// Pointwise Addition by Scalar
        /**
          Not recommended for performance unless necessary.
          Try to fuse arithmetic operations into forall or Stencil operations if possible.

          \param a_scale  A scalare of type T
          */
        inline BoxData<T,C,MEM,D,E>& operator+=(T a_scale);

        /// Pointwise Subtraction by Scalar
        /**
          Not recommended for performance unless necessary.
          Try to fuse arithmetic operations into forall or Stencil operations if possible.

          \param a_scale  A scalare of type T
          */
        inline BoxData<T,C,MEM,D,E>& operator-=(T scale);

        /// Pointwise Multiplication by Scalar
        /**
          Not recommended for performance unless necessary.
          Try to fuse arithmetic operations into forall or Stencil operations if possible.

          \param a_scale  A scalare of type T
          */
        inline BoxData<T,C,MEM,D,E>& operator*=(T scale);

        /// Pointwise Division by Scalar
        /**
          Not recommended for performance unless necessary.
          Try to fuse arithmetic operations into forall or Stencil operations if possible.

          \param a_scale  A scalare of type T
          */
        inline BoxData<T,C,MEM,D,E>& operator/=(T scale);

        ///@}

        ///////////////////////////////////////////////////////////////////////////////////////
        /** @name Utility */
        ///@{

        /// Initialize All Values
        /**
          Set all values equal to a constant. 

          \param a_val    A constant value. 
        */
        inline void setVal(const T& a_val);
        
        // TODO: Deprecated because Threshold ops don't work with std::complex
        /// Clamp Values 
        /**
          Threshold data between an input max and min  

          \param a_min      Minimum value
          \param a_max      Maximum value
        */
        //inline void clampVal(const T& a_min, const T& a_max);

        /// Initialize All Values as Zero
        /**
          Set all values equal to a zero.

        */
        inline void setToZero();
 
        /// Set All Values in Box
        /**
          Selectively initializes values for a specified Box. If <code>a_box</code>
          is not contained in <code>this->box()</code>, the initialization is done
          on the intersection.
          Avoid calling this function as much as possible in performance critical code

          \param a_val    A constant value.
          \param a_box    Domain to set to a_val
          */
        inline void setVal(const T& a_val, const Box& a_box);

        /// Set All Values of Component in Box
        /**
          Avoid calling this function as much as possible in performance critical code

          \param a_val    A constant value. 
          \param a_box    Domain to set to a_val
          \param a_c      First index
          \param a_d      (Optional) Second index. (default: 0)
          \param a_e      (Optional) Third index. (default: 0)
        */
        inline void setVal(const T& a_val,
                const Box& a_box,
                int a_c,
                int a_d = 0,
                int a_e = 0);

        /// Generic Reduction (Global)
        /**
          Computes a reduction operation
          Stores the result in a Reduction operator.

          \param a_rxn    A reduction operator 
          */
        template<Proto::Operation OP>
        inline void reduce(Reduction<T,OP,MEM>& a_Rxn) const;

        /// Generic Reduction (Componentwise)
        /**
          Computes a reduction operation on a given component
          Stores the result in a Reduction operator.

          \param a_rxn    A reduction operator 
          \param a_c      First tensor index.
          \param a_d      (Optional) Second index. (default: 0)
          \param a_e      (Optional) Third index. (default: 0)
          */
        template<Proto::Operation OP>
        inline void reduce( Reduction<T,OP,MEM>& a_Rxn, int a_c, int a_d = 0, int a_e = 0) const;

        /// Absolute Maximum Value (Global)
        /**
          Returns the maximum absolute value over the entire data set
          */
        inline T absMax() const;

        /// Absolute Maximum Value (Componentwise)
        /**
          Returns the maximum absolute value of a given component

          \param a_c      First tensor index.
          \param a_d      Second tensor index.
          \param a_e      Third tensor index.
          */
        inline T absMax(int a_c, int a_d = 0, int a_e = 0) const;

        /// Absolute Maximum Value (Global, Reduction)
        /**
          Computes the global maximum absolute value
          Stores the result in a Reduction operator.

          \param a_rxn    A reduction operator 
          */
        inline void absMax(Reduction<T,Abs,MEM>& a_Rxn) const;

        /// Absolute Maximum Value (Componentwise, Reduction)
        /**
          Computes the maximum absolute value of a given component
          Stores the result in a Reduction operator.

          \param a_rxn    A reduction operator 
          \param a_c      First tensor index.
          \param a_d      Second tensor index.
          \param a_e      Third tensor index.
          */
        inline void absMax(Reduction<T,Abs,MEM>& a_Rxn, int a_c, int a_d = 0, int a_e = 0) const;

        /// Minimum Value (Global)
        /**
          Returns the minimum value over the entire data set
          */
        inline T min() const;

        /// Minimum Value (Componentwise)
        /**
          Returns the minimum value of a given component

          \param a_c      First tensor index.
          \param a_d      Second tensor index.
          \param a_e      Third tensor index.
          */
        inline T min(int a_c, int a_d = 0, int a_e = 0) const;

        /// Minimum Value (Global, Reduction)
        /**
          Computes the minimum value over the entire data set
          Stores the result in a reduction operator

          \param a_rxn    A reduction operator
          */
        inline void min(Reduction<T, Min, MEM>& a_rxn) const;

        /// Minimum Value (Componentwise, Reduction)
        /**
          Computes the minimum value of a given component
          Stores the result in a reduction operator

          \param a_rxn    A reduction operator
          \param a_c      First tensor index.
          \param a_d      Second tensor index.
          \param a_e      Third tensor index.
          */
        inline void min(Reduction<T, Min, MEM>& a_rxn, int a_c, int a_d = 0, int a_e = 0) const;

        /// Maximum Value (Global)
        /**
          Returns the maximum value over the entire data set
          */
        inline T max() const;

        /// Maximum Value (Componentwise)
        /**
          Returns the maximum value of a given component

          \param a_c      First tensor index.
          \param a_d      Second tensor index.
          \param a_e      Third tensor index.
          */
        inline T max(int a_c, int a_d = 0, int a_e = 0) const;

        /// Max (Global, Reduction)
        /**
          Computes the maximum of values over the entire data set
          Stores the result in a reduction operator

          \param a_rxn    A reduction operator.
          */
        inline void max(Reduction<T, Max, MEM>& a_rxn) const;

        /// Max (Componentwise, Reduction)
        /**
          Computes the maximum value of a given component
          Stores the result in a reduction operator

          \param a_rxn    A reduction operator.
          \param a_c      First tensor index.
          \param a_d      Second tensor index.
          \param a_e      Third tensor index.
          */
        inline void max(Reduction<T, Max, MEM>& a_rxn, int a_c, int a_d = 0, int a_e = 0) const;

        /// Sum (Global)
        /**
          Return the global sum
          */
        inline T sum() const;

        /// Sum (Componentwise)
        /**
          Returns the sum of a particular component

          \param a_c      First tensor index.
          \param a_d      Second tensor index.
          \param a_e      Third tensor index.
          */
        inline T sum(int a_c, int a_d = 0, int a_e = 0) const;

        /// Sum (Global, Reduction)
        /**
          Computes the global sum
          Stores the result in a reduction operator

          \param a_rxn    A reduction operator.
          */
        inline void sum(Reduction<T, Sum, MEM>& a_rxn) const;

        /// Sum (Componentwise)
        /**
          Computes the sum of a particular component
          Stores the result in a reduction operator

          \param a_rxn    A reduction operator.
          \param a_c      First tensor index.
          \param a_d      Second tensor index.
          \param a_e      Third tensor index.
          */
        inline void sum(Reduction<T, Sum, MEM>& a_rxn, int a_c, int a_d = 0, int a_e = 0) const;

        /// Integrate (Isotropic Grid)
        /**
          Approximate the definite integral using an isotropic grid spacing as
          <code>sum(a_c, a_d, a_e) * pow(a_dx, DIM) </code>.
          
          \param a_dx   Isotropic grid spacing.
          \param a_c    First tensor index.
          \param a_d    Second tensor index.
          \param a_e    Third tensor index.
        */
        inline T integrate(T a_dx, int a_c = 0, int a_d = 0, int a_e = 0) const;
        
        /// Integrate (Anisotropic Grid)
        /**
          Approximate the definite integral using an isotropic grid spacing as
          <code>sum(a_c, a_d, a_e) * pow(a_dx, DIM) </code>.
          
          \param a_dx   Anisotropic grid spacings.
          \param a_c    First tensor index.
          \param a_d    Second tensor index.
          \param a_e    Third tensor index.
        */
        inline T integrate(
                std::array<T, DIM> a_dx,
                int a_c = 0, int a_d = 0, int a_e = 0) const;

        /// Shift Domain
        /**
          Shifts the domain Box of this, moving all data along with it.
          e.g. data associated with Point p will now be associated with Point p + a_shift.

          \param a_shift  A Point interpreted as a shift vector.
          */
        inline void shift(const Point& a_shift);

        /// Buffer Size
        inline size_t linearSize() const;

        /// Buffer Write
        /**
          Write a subset of data from <code> *this </code> into a C-array buffer.

          \param a_buffer Destination buffer
          \param a_box    Domain to write from
          \param a_comps  Components to write from
          */
        inline void linearOut(
                void*       a_buffer,
                const Box&  a_box,
                CInterval   a_comps) const;

        //FIXME: Only works for vector/scalar valued BoxData
        /// Buffer Write (Primitive)
        /**
          Write a subset of data from <code> *this </code> into a C-array buffer.
          Overload using plain-old-data inputs, useful on GPUs

          \param a_buffer       Destination buffer
          \param a_box          Domain to write from
          \param a_startcomp    First component to write from
          \param a_numcomps     Number of components to write
          */
        inline void linearOut(
                void*           a_buffer,
                const Box&      a_box,
                unsigned int    a_startcomp,
                unsigned int    a_numcomps) const;

        /// Buffer Read
        /**
          Read data into *this from a C-array buffer populated by
          <code>BoxData<T,C,MEM,D,E>::linearOut(...)</code>.
          <code>this->box()</code> may be shifted with respect to the BoxData 
          that wrote to the buffer, but the dimensions of the Box as well as the 
          tensor structure must be the same.

          \param a_buffer Source buffer
          \param a_box    Domain to read from
          \param a_comps  Components to read from
          */
        inline void linearIn(
                void*       a_buffer,
                const Box&  a_box,
                CInterval   a_comps);


      inline T* rawPtr() const
      {
        return m_rawPtr;
      }
      
      inline ::std::shared_ptr<T>  getData() const
      {
        return m_data; ///< Data array
      }

      /// Buffer Read
        /**
          Read data into *this from a C-array buffer populated by
          <code>BoxData<T,C,MEM,D,E>::linearOut(...)</code>.
          Primitive data overload for use with GPUs.

          \param a_buffer Source buffer
          \param a_box    Domain to read from
          \param a_comps  Components to read from
          */
        inline void linearIn(
                void*               a_buffer,
                const ::Proto::Box& a_bx,
                unsigned int        a_startcomp,
                unsigned int        a_numcomps);

        /// Contains CInterval
        /**
          Checks if an CInterval object is a subset of the component space
          of *this. This function is mostly used internally, but it can
          be useful for checking inputs to CopyTo and the slicing functions
          */
        inline bool contains(CInterval a_interval) const;

        /// Check Aliasing
        /**
          Returns true if <code> *this </code> and <code> a_src </code> are aliased
          to the same data buffer. This function also returns true if one array
          is a slice of another.
          */
        template<unsigned int CC, unsigned char DD, unsigned char EE>
        inline bool isAlias(const BoxData<T, CC, MEM, DD, EE>& a_src) const;

        /// Print
        /**
          Default Print method. Outputs Domain Box and extrema of *this
        */
        void print() const;

        /// Print Data
        /**
          Pretty prints *all* of the data in *this. Prettiness may vary depending
          on the domain size. This function is purely for debugging purposes and should not 
          be used in performance code.

          \param a_prec   (Optional) Desired precision for fixed decimal output. Defaults to 2.
          */    
        void printData(int a_prec = 2) const;

        /// Print Data in Box
        /**
          Pretty prints the data in *this within a given Box for all components.
          Input Boxes of size 16^DIM or smaller are recommended to maintain prettiness.
          This function is purely for debugging purposes and should not be used in 
          performance code.

          \param a_box    Desired subset for printing
          \param a_prec   (Optional) Desired precision for fixed decimal output. Defaults to 2.
        */    
        void printData(const Box& a_box, int a_prec = 2) const;

        /// Print Component Data in Box
        /**
          Pretty prints the data in *this within a given Box for a single component.
          Input Boxes of size 16^DIM or smaller are recommended to maintain prettiness.
          This function is purely for debugging purposes and should not be used in 
          performance code.


          \param a_box    Desired subset for printing
          \param a_c      First tensor index. Defaults to 0.
          \param a_d      Second tensor index. Defaults to 0.
          \param a_e      Third tensor index. Defaults to 0.
          \param a_prec   (Optional) Desired precision for fixed decimal output. Defaults to 2
        */    
        void printData(const Box& a_box, int a_c, int a_d, int a_e, int a_prec = 2) const;
        ///@}
        
        ///////////////////////////////////////////////////////////////////////////////////////
        // Friends

        template<class TT>
            friend class Stencil;

        template<class TT, unsigned int CC, unsigned char DD, unsigned char EE,MemType MM>
            friend BoxData<TT,CC,MM,DD,EE> alias(BoxData<TT,CC,MM,DD,EE>& a_original,
                    const Point&              a_shift);

        template<class TT, unsigned int CC, unsigned char DD, unsigned char EE,MemType MM>
            friend const BoxData<TT,CC,MM,DD,EE> alias(const BoxData<TT,CC,MM,DD,EE>& a_original,
                    const Point&                a_shift);
        
        template<typename _T, unsigned int _C, MemType _MEM, unsigned char _D, unsigned char _E>
            friend BoxData<_T,1,_MEM,1,1> slice(
                    const BoxData<_T,_C,_MEM,_D,_E>& a_src,
                    unsigned int a_c,
                    unsigned int a_d,
                    unsigned int a_e);

        template<typename _T, unsigned int _C, unsigned char _CC, MemType _MEM>
            friend BoxData<_T,_CC,_MEM,1,1> slice(
                    const BoxData<_T,_C,_MEM,1,1>& a_src,
                    unsigned int a_nstart);

        // TODO: Why?
        size_t size(const Box&      a_box,
                const size_t a_comps) const
        {
            return a_box.size() * (sizeof(T) * a_comps);
        }


        // TODO: Why?                
        size_t charsize(const ::Proto::Box&      a_bx,
                unsigned int a_startcomp,
                unsigned int a_numcomps) const
        {
            return a_bx.size()*sizeof(T)*a_numcomps;
        }
        
        // TODO: Chombo Compatibility?
        static int preAllocatable()
        {
            return 0; // static preAllocatable
        }
        
        // TODO: What is this for?
        static int memTypeAllocation()
        {
            return 1;
        }

        private:

        bool                    m_stackAlloc;
        Box                     m_box; ///< Box defining the domain of *this
        ::std::shared_ptr<T>    m_data; ///< Data array
        T*                      m_rawPtr; ///< Raw pointer to the data
    }; // end class BoxData

///////////////////////////////////////////////////////////////////////////////////////////////
/** @name Alias and Slice Operators 
    The alias and slice operations facilitate BoxData operations while avoiding unnecessary copies.
    See the sample code below for an explanation of the syntax.

    Example <code> alias </code> usage:
    \snippet Snippets.cpp proto_alias
    
    Example <code> slice </code> usage:
    \snippet Snippets.cpp proto_slice
*/

///@{

    /// Alias (Non-Const)
    /**
      Creates a read-write alias to a mutable source BoxData with an optional shift.
      
      \param a_original Source data
      \param a_shift    (Optional) Domain shift for the aliased BoxData.
    */
    template<class T, unsigned int C=1, unsigned char D=1, unsigned char E=1,
        MemType MEM=MEMTYPE_DEFAULT>
    BoxData<T,C,MEM,D,E> alias(
            BoxData<T,C,MEM,D,E>&   a_original,
            const Point&                a_shift = Point::Zeros());

    /// Alias (Const)
    /**
      Creates a read-only alias to a const source BoxData with an optional shift.
      
      \param a_original Source data
      \param a_shift    (Optional) Domain shift for the aliased BoxData.
    */
    template<class T, unsigned int C = 1, unsigned char D = 1, unsigned char E = 1,
        MemType MEM=MEMTYPE_DEFAULT>
    const BoxData<T,C,MEM,D,E> alias(
            const BoxData<T,C,MEM,D,E>& a_original,
            const Point&                    a_shift=Point::Zeros());

    /// Slice Arbitrary Component (Non-Const)
    /**
      Creates a scalar BoxData read-write alias to a prescribed component. 
    
      \param a_src  Source data
      \param a_c    First component index
      \param a_d    (Optional) Second component index (default: 0)
      \param a_e    (Optional) Third component index (default: 0)
    */
    template<typename T, unsigned int C, MemType MEM=MEMTYPE_DEFAULT,
        unsigned char D=1, unsigned char E=1>
    BoxData<T,1,MEM,1,1> slice(
            const BoxData<T,C,MEM,D,E>& a_src,
            unsigned int a_c,
            unsigned int a_d = 0,
            unsigned int a_e = 0);

    /// Vector Slice (Non-Const)
    /**
      Creates a vector BoxData read-write alias to a prescribed component range
      
      \param a_src          Source data
      \param a_startcomp    First component of the slice
    */
    template<typename T, unsigned int C, unsigned char CC, MemType MEM=MEMTYPE_DEFAULT>
    BoxData<T,CC,MEM,1,1> slice(
            const BoxData<T,C,MEM,1,1>& a_src,
            unsigned int                    a_startcomp);

///@}

#include "implem/Proto_BoxDataImplem.H"
} //end Proto namespace
#endif //end include guard
