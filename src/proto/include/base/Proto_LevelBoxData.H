#pragma once
#ifndef _PROTO_LEVELBOXDATA_H_
#define _PROTO_LEVELBOXDATA_H_

#include "Proto_MemType.H"
#include "Proto_Memory.H"
#include "Proto_Point.H"
#include "Proto_Box.H"
#include "Proto_Copier.H"
#include "Proto_DisjointBoxLayout.H"
#include "Proto_DataIterator.H"
#include "Proto_NeighborIterator.H"
#include "Proto_WriteBoxData.H"
#include "Proto_MayDay.H"
#include "Proto_SPMD.H"
#include "Proto_Operator.H"
#include <cstdlib> //for size_t
#include <iostream>
#include <cstring> // for Writing data to .vtk files.

namespace Proto 
{
    // Forward Declarations
    template<typename T, unsigned int C, MemType MEM, Centering CTR>
    class LevelBoxData;

// =======================================================================
// LEVEL BOX DATA COPIER OP
    
    /// Level Copier Operator
    /**
        Copier operator which provides the necessary interface to copy to and from
        LevelBoxData using the Copier interface. See Proto_Copier.H.
    */
    template
    <typename T, unsigned int C, MemType SRC_MEM, MemType DST_MEM, Centering CTR>
    class LevelCopierOp
    {
        public:
        LevelBoxData<T,C,SRC_MEM,CTR>*  m_src;
        LevelBoxData<T,C,DST_MEM,CTR>*  m_dst;

        inline LevelCopierOp(){};
        
        inline LevelCopierOp(
            LevelBoxData<T, C, SRC_MEM, CTR>& a_src,
            LevelBoxData<T, C, DST_MEM, CTR>& a_dst);
        
        /// Get Serialized Size
        inline int linearSize(const Box& a_bx,
                const DataIndex<BoxPartition>& a_index) const;
        
        /// Serial Write To Buffer
        inline void linearOut(void* a_buf, 
                const Box& a_bx,
                const DataIndex<BoxPartition>& a_index) const;
        
        /// Serial Read From Buffer
        inline void linearIn(void* a_buf, 
                const Box& a_bx,
                const DataIndex<BoxPartition>& a_index) const;
        
        /// Local Copy
        inline void localCopy(
            const Box& a_domain, const DataIndex<BoxPartition>& a_domainIndex,
            const Box& a_range,  const DataIndex<BoxPartition>& a_rangeIndex) const;
    };
// =======================================================================
// LEVEL BOX DATA COPIERS
    
    /// Level Copier
    /**
        A Copier used to implement copying of data between LevelBoxData objects.
    */
    template<typename T, unsigned int C, MemType SRC_MEM, MemType DST_MEM, Centering CTR>
    class LevelCopier : public Copier<LevelCopierOp<T, C, SRC_MEM, DST_MEM, CTR>, BoxPartition, SRC_MEM, DST_MEM>
    {
        public:
        
        inline LevelCopier(){};
        
        /// Build Copier Motion Plan
        inline void buildMotionPlans(LevelCopierOp<T, C, SRC_MEM, DST_MEM, CTR>& a_op);
    };
    
    /// Exchange Copier
    /**
        A Copier used to execute the data copies necessary for filling ghost regions
        within a LevelBoxData.
    */
    template<typename T, unsigned int C, MemType MEM, Centering CTR>
    class LevelExchangeCopier : public Copier<LevelCopierOp<T, C, MEM, MEM, CTR>, BoxPartition, MEM, MEM>
    {
        public:

        inline LevelExchangeCopier(){};
        
        /// Build Copier Motion Plan
        inline void buildMotionPlans(LevelCopierOp<T, C, MEM, MEM, CTR>& a_op);
    };

// =======================================================================
// LEVEL BOX DATA
    
    /// Level Box Data
    /**
        A distributed data holder implementing the abstraction of a union of logically rectangular patches.
    */
    template
    <typename T, unsigned int C=1, MemType MEM=MEMTYPE_DEFAULT, Centering CTR = PR_CELL>
    class LevelBoxData
    {
        public:

        /// Default constructor
        /**
            Constructs an undefined LevelBoxData. Can be initialized through define(...)
        */
        inline LevelBoxData();

        /// Layout Constructor
        /**
            Constructs a non-trivial LevelBoxData by calling define.
        */
        inline LevelBoxData(const DisjointBoxLayout& a_dbl, const Point& a_ghost);

        /// Slice Constructor
        /**
            Constructs a LevelBoxData by aliasing C compontents of <code>a_src</code> starting with
            component <code>a_startComp</code>. The operation fails by assertion if 
            <code>a_startComp + C > CC </code>

            \param a_src        Source Data
            \param a_startComp  Starting component
        */
        template<unsigned int CC>
        inline LevelBoxData(const LevelBoxData<T, CC, MEM, CTR>& a_src, unsigned int a_startComp);

        /// Define
        /**
            Lazily constructs a non-trivial LevelBoxData
        */
        inline void define(const DisjointBoxLayout& a_dbl, const Point& a_ghost);

        /// Size
        /**
            Returns the number of boxes in the union across all processes.
        */
        inline unsigned int size() const;

        /// Initialize To Zero
        /**
            Sets all values in all patches equal to zero. This function is mostly for 
            debugging purposes and should be used sparingly if at all in production code.
        */
        inline void setToZero(int a_comp = -1);
        
        /// Initialize as Constant
        /**
            Set's all values in all patches equal to a constant value (including ghost cells.)
        
            \param a_value  A constant value
        */
        inline void setVal(T a_value, int a_comp = -1);
        
        /// Set the Level Boundary Ghost Cells
        /**
            Set's all values in non-periodic boundary ghost cells equal to a constant value.
        
            \param a_value  A constant value
        */
        inline void setBoundary(T a_value, int a_comp = -1);
        
        // Deprecated because this functionality doesn't work with std::complex
        /// Clamp Values
        /**
            Clamps all values of a component in *this between a min and max (including ghost cells).
            
            \param a_min  Min value
            \param a_max  Max value
        */
        //inline void clampVal(T a_min, T a_max, int a_comp = -1);

        /// Function Initialization
        /**
            Initializes *this using a function with a signature similar to that used by
            Proto::forallInPlace_p. This function is purely for convenience. 

            Func signature:
            void a_func(Point& a_point, Var<T, C, MEM>& a_thisData, ... )

            Unlike Proto::forall and its variants, this function cannot manage additional
            data holders (e.g. Srcs should not include a LevelBoxData or similar complex object).
            Elements of a_srcs should be plain-old-data such as single values, Points, or Boxes. 
        */
        template<typename Func, typename... Srcs>
        inline void initialize(Func& a_func, Srcs... a_srcs);
        
        /// Function Initialization (With Convolution)
        /**
            Initialize *this using a user defined function for the cell-centered quantity
            and then convert the result to a cell-averaged quantity through 4th order
            convolution.
        */
        template<typename Func, typename... Srcs>
        inline void initConvolve(Func& a_func, Srcs... a_srcs);

        /// Exchange
        /**
            Copies data from the valid regions of *this into ghost regions. When MPI is
            enabled, this function also takes care of any necessary communication between
            patches on different processes.
        */
        inline void exchange();

        /// Copy To
        /**
            Copies data from the valid regions of this into the valid regions of another
            LevelBoxData. This function copies both valid data and ghost regions if applicable.
            If <code>a_dest</code> has a larger ghost region then <code>*this</code>, the user
            is expected to call <code>exchange</code> to correctly fill the extra ghost cells.

            \param a_dest   The destination to which data in *this will be copied.
        */
        template<MemType DST_MEM>
        inline void copyTo(LevelBoxData<T, C, DST_MEM, CTR>& a_dest) const;
        
        template<MemType DST_MEM>
        inline void copyToSimple(LevelBoxData<T, C, DST_MEM, CTR>& a_dest) const;
       
        /// Copy With Coarsen
        /**
            Coarsen's the data in *this by applying a geometric average before copying 
            into the destination.
        */
        inline void coarsenTo(LevelBoxData& a_dest, Point a_refRatio) const;

        /// Linear Offset
        /**
          Returns the serial index of the first data element of this stored on
          this a_proc.
         */
        inline unsigned int offset(int a_proc) const;
        inline unsigned int offset() const { return offset(Proto::procID()); }

        /// Patch Size
        /**
          Computes the size (in data cells) of a single patch in *this.
          The resulting patch ACCOUNTS FOR GHOST CELLS
         */
        inline unsigned int patchSize() const;

        /// Linear In (Patch, Box)
        /**
          Reads data from a serial buffer, populating the patch associated with 
          a_index within the subset specified by a_box. 

          \param a_buffer     The read buffer
          \param a_box        A box defining the subset to read into
          \param a_index      The index defining the patch to read into
         */ 
        inline void linearIn(void* a_buffer, 
                const Box& a_box,
                const DataIndex<BoxPartition>& a_index);

        /// Linear Out (Patch, Box)
        /**
          Writes data to a serial buffer, from the patch associated with 
          a_index within the subset specified by a_box. 

          \param a_buffer     The write buffer
          \param a_box        A box defining the subset to write from
          \param a_index      The index defining the patch to write from
         */ 
        inline void linearOut(void* a_buffer, 
                const Box& a_box,
                const DataIndex<BoxPartition>& a_index) const;

        /// Serial Size (Patch, Box)
        /**
          Computes the size in bytes of the serial data buffer needed to hold
          the patch associated with a_index within the subset specified by a_box.

          \param a_box        A box defining the subset of a patch
          \param a_index      The index defining the patch
         */ 
        inline size_t linearSize(const Box& a_box,
                const DataIndex<BoxPartition>& a_index) const;
        /// Linear In (Patch)
        /**
          Overload of LinearIn that always reads in a full patch
         */
        inline void linearIn(void* a_buffer, 
                const DataIndex<BoxPartition>& a_index);
        /// Linear Out (Patch)
        /**
          Overload of LinearOut that always writes out a full patch
         */
        inline void linearOut(void* a_buffer, 
                const DataIndex<BoxPartition>& a_index) const;
        /// Size (Patch)
        /**
          Overload of Size that always computes the size of a full patch
         */
        inline size_t linearSize( const DataIndex<BoxPartition>& a_index) const;

        /// Linear In (All Local Data)
        /**
          Overload of LinearIn that reads all local data
         */
        inline void linearIn(void* a_buffer);
        /// Linear Out (All Local Data)
        /**
          Overload of LinearOut that writes out all local data
         */
        inline void linearOut(void* a_buffer) const;
        /// Size (All Local Data)
        /**
          Overload of Size that computes the full local size of the data in *this
         */
        inline size_t linearSize() const;

        /// Get Box Layout
        inline DisjointBoxLayout layout() const { return m_layout; }

        /// Get Ghost Size
        inline Point ghost() const { return m_ghost; }
       
        /// Increment
        /**
            Add a scaled multiple of another LevelBoxData.
        */
        inline void increment(LevelBoxData<T, C, MEM, CTR>& a_data, T a_scale = 1.0);
       
        /// Multiply By Scalar
        inline void operator*=(T a_scale);
        
        /// Reduction
        /**
            Computes a reduction for a component over all valid cells (ghost region is ignored)
        */
        template<Proto::Operation OP>
        inline double reduce(unsigned int a_comp = 0) const;
       
        /// Maximum Absolute Value
        /**
            Computes the maximum absolute value for a component over all valid cells (ghost region is ignored)
        */
        inline double absMax(unsigned int a_comp = 0) const;
        
        /// Sum
        /**
            Computes the sum of a component over all valid cells (ghost region is ignored)
        */
        inline double sum(unsigned int a_comp = 0) const;
        
        /// Max
        /**
            Computes the max value of a component over all valid cells (ghost region is ignored)
        */
        inline double max(unsigned int a_comp = 0) const;
        
        /// Min
        /**
            Computes the min value of a component over all valid cells (ghost region is ignored)
        */
        inline double min(unsigned int a_comp = 0) const;
        
        /// Integral
        /**
            Compute the integral over this using an isotropic grid spacing.

            \param a_dx     Isotropic grid spacing
            \param a_comp   A component to integrate over
        */
        inline double integrate(double a_dx, unsigned int a_comp = 0) const;
        
        /// Integral (Anisotropic)
        /**
            Compute the integral over this using an isotropic grid spacing.

            \param a_dx     Isotropic grid spacing
            \param a_comp   A component to integrate over
        */
        inline double integrate(std::array<double, DIM> a_dx, unsigned int a_comp = 0) const;
        
        /// Integral of Absolute Value
        /**
            Compute the integral of the absolute value over this using an isotropic grid spacing.

            \param a_dx     Isotropic grid spacing
            \param a_comp   A component to integrate over
        */
        inline double integrateAbs(double a_dx, unsigned int a_comp = 0) const;
        
        /// Integral of Absolute Value (Anisotropic)
        /**
            Compute the integral of the absolute value over this using an isotropic grid spacing.

            \param a_dx     Isotropic grid spacing
            \param a_comp   A component to integrate over
        */
        inline double integrateAbs(std::array<double, DIM> a_dx, unsigned int a_comp = 0) const;
      
        /// Query Layout Compatibility
        template<typename _T, unsigned int _C, MemType _MEM, Centering _CTR>
        inline bool compatible(const LevelBoxData<_T, _C, _MEM, _CTR>& a_data) const;

        /// Query Layout Compatibility
        inline bool compatible(const DisjointBoxLayout& a_layout) const;
        /// Query Layout Compatibility
        inline bool compatible(const DataIterator<BoxPartition>& a_iter) const;
        /// Query Layout Compatibility
        inline bool compatible(const DataIndex<BoxPartition>& a_index) const;

        //=======================================
        // DATA ITERATOR INTERFACE

        /// Get Patch
        inline BoxData<T,C,MEM>& operator[](const DataIndex<BoxPartition>& a_index);
        /// Get Const Patch
        inline const BoxData<T,C,MEM>& operator[](const DataIndex<BoxPartition>& a_index) const;
        
        /// Get Iterator
        inline DataIterator<BoxPartition> begin() const;
        inline DataIterator<BoxPartition> end() const;

        /// For debugging purposes.
        int s_verbosity = 0;
        
        private: 

        std::vector<shared_ptr<BoxData<T, C, MEM> > > m_data;
        Point                                   m_ghost;
        DisjointBoxLayout                       m_layout;
        bool                                    m_isDefined;
        LevelExchangeCopier<T, C, MEM, CTR>     m_exchangeCopier;
        
        // disallow copy constructors and assignment operators
        // to avoid very hard-to-find performance problems
        LevelBoxData<T, C, MEM, CTR>& operator= (const LevelBoxData<T, C, MEM, CTR>& a_rhs);
        LevelBoxData<T, C, MEM, CTR>(const LevelBoxData<T, C, MEM, CTR>& a_rhs);
    }; // end class LevelBoxData
    
    /// Interpolate Boundaries
    /**
        General utility function for interpolating data into coarse-fine boundaries.
        This version creates a temporary LevelBoxData representing the coarsened fine region.
    */
    template<typename T, unsigned int C, MemType MEM, Centering CTR>
    void interpBoundaries(
        LevelBoxData<T, C, MEM, CTR>& a_crse,
        LevelBoxData<T, C, MEM, CTR>& a_fine,
        InterpStencil<T>&       a_interp);
    
    /// Interpolate Boundaries
    /**
        General utility function for interpolating data into coarse-fine boundaries.
        This version does not create any temporary LevelBoxData but requires an additional
        input dataholder representing the coarsened fine region.
    */
    template<typename T, unsigned int C, MemType MEM, Centering CTR>
    void interpBoundaries(
        LevelBoxData<T, C, MEM, CTR>& a_crse,
        LevelBoxData<T, C, MEM, CTR>& a_fine,
        LevelBoxData<T, C, MEM, CTR>& a_crseFine,
        InterpStencil<T>&       a_interp);

    /// Average Down
    /**
        General utility function for averaging fine data onto coarse data
    */
    template<typename T, unsigned int C, MemType MEM, Centering CTR>
    void averageDown(
        LevelBoxData<T, C, MEM, CTR>& a_crse,
        LevelBoxData<T, C, MEM, CTR>& a_fine,
        Point                         a_refRatio);
    
    /// Average Down (Scalar Refinement Ratio)
    /**
        General utility function for averaging fine data onto coarse data
    */
    template<typename T, unsigned int C, MemType MEM, Centering CTR>
    void averageDown(
        LevelBoxData<T, C, MEM, CTR>& a_crse,
        LevelBoxData<T, C, MEM, CTR>& a_fine,
        int                           a_refRatio);
    
    /// Average Down
    /**
        General utility function for averaging fine data onto coarse data
    */
    template<typename T, unsigned int C, MemType MEM, Centering CTR>
    void averageDown(
        LevelBoxData<T, C, MEM, CTR>& a_crse,
        LevelBoxData<T, C, MEM, CTR>& a_fine,
        LevelBoxData<T, C, MEM, CTR>& a_crseFine,
        Point                         a_refRatio);
    
    /// Average Down (Scalar Refinement Ratio)
    /**
        General utility function for averaging fine data onto coarse data
    */
    template<typename T, unsigned int C, MemType MEM, Centering CTR>
    void averageDown(
        LevelBoxData<T, C, MEM, CTR>& a_crse,
        LevelBoxData<T, C, MEM, CTR>& a_fine,
        LevelBoxData<T, C, MEM, CTR>& a_crseFine,
        int                           a_refRatio);
    

#include "implem/Proto_LevelBoxDataImplem.H"
} //end namespace Proto
#endif
