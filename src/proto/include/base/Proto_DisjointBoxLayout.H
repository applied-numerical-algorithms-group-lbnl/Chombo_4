#pragma once
#ifndef _PROTO_DISJOINT_BOX_LAYOUT_
#define _PROTO_DISJOINT_BOX_LAYOUT_

#include "Proto_MayDay.H"
#include "Proto_Point.H"
#include "Proto_Box.H"
#include "Proto_ProblemDomain.H"
#include "Proto_BoxPartition.H"
#include "Proto_DataIterator.H"
#include "Proto_SPMD.H"
#include <cstdlib> //for size_t
#include <iostream>
#include <iomanip>
#include <stack>
#include <memory>
#include <array>
#include <unordered_map>

namespace Proto 
{ 
    /// Disjoint Box Layout
    /**
      Layout of disjoint, logically rectangular regions of DIM-dimensional space.
      This object is a simplified generalization of the DisjointBoxLayout class
      used in Chombo, now specialized to use fixed box sizes. 
    */
    class DisjointBoxLayout
    { 
        template<typename P> friend class DataIndex;
        template<typename P> friend class DataIterator;
        friend class NeighborIterator;
        friend std::ostream& operator<< (std::ostream& os, const DisjointBoxLayout& a_dbl);

        public:

        /// Default Constructor 
        inline DisjointBoxLayout();

        /// Constructor (Full Domain)
        /**
          Creates a DisjointBoxLayout by tiling an input ProblemDomain completely with boxes of
          a predetermined size. If MPI is used, boxes will be distributed over available
          processors using Morton ordering.
          If <code>a_problemDomain.coarsenable(a_boxSize) != true</code> this results in an error.

          \param a_problemDomain    A ProblemDomain
          \param a_boxSize          A set of box sizes
          */
        inline DisjointBoxLayout(
                const ProblemDomain   & a_problemDomain, 
                const Point           & a_boxSize);

        //TODO: Either improve the interface or the documentation
        /// Constructor (Partial Domain)
        /**
          Creates a DisjointBoxLayout by partially tiling an input ProblemDomain with boxes
          according to a set of Points. Each element of a_coarsenedPatches represents a box 
          in the layout which has been coarsened by boxSize (resulting in a Box with a
          single element). In parallel applications, <code>a_coarsenedPatches </code> should contain
          the patches in the layout on ALL ranks, not just locally.

          \param a_problemDomain    ProblemDomain containing the DBL.
          \param a_coarsenedPatches An array of points corresponding to patches in the tiled layout
          \param a_boxSize          Possibly anisotropic size of each Box in the layout
        */
        inline DisjointBoxLayout(
                const ProblemDomain   & a_problemDomain,
                const vector<Point>   & a_coarsenedPatches,
                const Point           & a_boxSize);

        /// Constructor (Partial Domain - Simple)
        /**
          Creats a DisjointBoxLayout by partially covering an input ProblemDomain with all patches
          in an input region. The sub-region must also be coarsenable by the box size, and must be
          a contained within the problem domain.
          Useful for quick debugging, but not practically useful for most real problems.
          
          \param a_problemDomain    ProblemDomain containing the DBL.
          \param a_region           Sub-region of the domain which will be tiled in the layout.
          \param a_boxSize          Possibly anisotropic size of each Box in the layout.
        
        */ 
        inline DisjointBoxLayout(
                const ProblemDomain&    a_problemDomain,
                const Box&              a_region, 
                const Point&            a_boxSize);

        /// Direct Define
        /**
            Creates a DisjointBoxLayout which is linked to an existing BoxPartition object. 
            This constructor is used internally and is not recommended for public use.
        */
        inline void define(
                std::shared_ptr<BoxPartition>   a_partition,
                const Point&                    a_boxSize);
        
        /// Copy Constructor
        /**
          Copies share pointers to the internal structure hence this is relatively cheap.

          \param a_layout   Another DisjointBoxLayout
        */
        inline DisjointBoxLayout(const DisjointBoxLayout& a_layout);

        /// Define (Full Domain)
        /**
          Construct a DisjointBoxLayout lazily or redefine an existing one. 
          */
        inline void define(
                const ProblemDomain   & a_problemDomain, 
                const Point           & a_boxSize);

        ///  Define (Partial Domain)
        /**
          Construct a DisjointBoxLayout lazily or redefine an existing one. 
          */
        inline void define(
                const ProblemDomain   & a_problemDomain,
                const vector<Point>   & a_coarsenedPatches,
                const Point           & a_boxSize);

        ///  Define (Partial Domaini - Simple)
        /**
          Construct a DisjointBoxLayout lazily or redefine an existing one. 
          */
        inline void define(
                const ProblemDomain   & a_problemDomain,
                const Box             & a_region,
                const Point           & a_boxSize);

        /// Load Balance
        /**
            Distribute the load of this layout over the range of processors whose
            indices span <code>[a_startProc, a_endProc-1]</code>.
        */
        inline void loadBalance(unsigned int a_startProc, unsigned int a_endProc);
       
        /// Load Assign
        /**
            Manually assign the load of this layout using the syntax
            <code> loadAssign(P0, N0, P1, N1, ..., Pk, Nk) </code>
            where process Pi is assigned Ni boxes. The sum over all Ni should be
            equal to the size of this layout.
        */
        template<typename... Args>
        inline void loadAssign(Args... a_args);

        /// Load Assign
        /**
            Manually assign the load of this layout. The input is a vector of
            <code>pair<int, unsigned int></code> where the first entry of each pair is a 
            process and the second is the number of boxes to be assigned to that process.
        */
        inline void loadAssign(std::vector<std::pair<int, unsigned int>>& a_assignment);

        /// Assignment
        /**
          Copies share pointers to the same internal structure
        */
        inline DisjointBoxLayout& operator=(const DisjointBoxLayout& a_input);

        /// Equality Operator
        inline bool operator==(const DisjointBoxLayout& a_input) const;

        /** @name Access operators. 
        */
        
        /// Get All Boxes
        /**
            Get a copy of the patches in the layout. The output is a vector
            of <code>pair<Point, int></code> where for each element
            <code>item</code>, <code>item.first</code> contains a Point corresponding
            to the Box <code>Box(item.first, item.first).refine(boxSize())</code>. 
            <code>item.second</code> contains the integer rank of the process which
            handles the patch.

            Note that this function returns data for all patches on all parallel processes.
        */
        inline const std::vector<pair<Point, unsigned int>> boxes() const;

        /// Index Access
        /**
          Return the Box associated with a DataIndex.
          Fails if the index is not compatible with *this.

          \param a_index    A DataIndex
        */
        inline Box operator[](const DataIndex<BoxPartition>& a_index) const;

        /// Box Indexing
        /**
          Return the Box associated with a DataIndex. Identical to <code>operator[]</code>
          Fails if the index is not compatible with *this.

          \param a_index    A DataIndex
        */
        inline Box box(const DataIndex<BoxPartition>& a_index) const;

        /// Point Indexing
        /**
          Return the Point version of the Box associated with a DataIndex.
          Fails if the index is not compatible with *this.

          \param a_index    A DataIndex
          */
        inline Point point(const DataIndex<BoxPartition>& a_index) const;

        /// Get Process ID
        /**
          Get the process ID associated with a DataIndex.
          Fails if the index is not compatible with *this.
        
          \param a_index    A DataIndex
        */
        inline int procID(const DataIndex<BoxPartition>& a_index) const; //consistent syntax
        //inline int procid(const DataIndex& a_index) const;

        /// Patch Offset
        /**
          Returns the linear offset (in number of boxes) of the first patch on proc a_proc.
          Used for serialization.
        
          \param a_proc     An MPI rank
        */
        inline unsigned int offset(int a_proc) const;

        /// This Patch Offset
        /**
          Returns the linear offset (in number of boxes) of the first patch on this proc.
          Used for linearization.
        */
        inline unsigned int offset() const { return offset(Proto::procID()); }

        /// Get Index From Point
        /**
          Find the DataIndex associated with an input patch point.
          This function is the inverse of DisjointBoxLayout::point(DataIndex).

          Input Points must be in <code>this->patchDomain()</code> or be a valid periodic 
          image thereof, otherwise this function will fail by assertion.

          \param a_patchPoint   A point representing a patch in the layout or a periodic image
          */
        inline DataIndex<BoxPartition> index(Point& a_patchPoint) const;

        /// Size
        /**
          Return the number of Boxes contained in *this including those on other processors.
        */
        inline unsigned int size() const;

        /// Local Size
        /**
          Return the number of Boxes on this processor.
          */
        inline unsigned int localSize() const;

        /// Global Indexing
        /**
          Get the DataIndex associated with a global integer index.
          Inputs are integers in <code>[0, this->size())</code>

          \param a_intIndex     A valid integer index
        */
        inline DataIndex<BoxPartition> index(unsigned int a_intIndex) const;

        /// Local Indexing
        /**
          Get the DataIndex associated with a local integer index.
          Inputs are integers in <code>[0, this->localSize())</code>
          
          \param a_intIndex     A valid integer index
          */
        inline DataIndex<BoxPartition> localIndex(unsigned int a_myIndexInt) const;

        /// Get Box Size
        /**
          Return the constant size of all boxes in the layout.
        */
        inline Point boxSize() const{return m_boxSize;};

        /// Get Problem Domain
        inline ProblemDomain domain() const {return m_problemDomain;};

        /// Get Problem Domain in Patch Space
        /**
          Returns the problem domain box in patch-space.
          */
        inline ProblemDomain patchDomain() const {return domain().coarsen(boxSize()); }

        /// Find Index
        /**
          Given a patch Point, return the DataIndex of the associated Box if
          it is in the layout including periodic images. If the Point does not
          represent a patch in this layout, this function returns <code>this->end()</code>
        
          \param a_pt   A patch Point.
        */
        inline DataIndex<BoxPartition> find(const Point& a_pt) const;

        /// Find Tile Point
        /**
            Returns true if the tile a_pt (or a periodic image thereof) is contained in the layout
        */
        inline bool contains(const Point& a_pt) const;

        /// On Domain Boundary
        /**
            Returns true if <code>a_pt</code> corresponds to a patch that is adjacent to a
            non-periodic domain boundary. If <code>a_pt</code> does not represent a valid patch
            in this layout, this function fails by assertion.

            \param a_pt A valid patch Point.
        */
        inline bool onDomainBoundary(const Point& a_pt) const;

        /// On Level Boundary
        /**
            Returns true if <code>a_pt</code> corresponds to a patch that is adjacent to a
            non-periodic layout boundary (including domain boundaries). In other words, this function
            returns false if and only if <code>a_pt</code> represents a patch that is completely
            surrounded by other patches in the layout.
            If <code>a_pt</code> does not represent a valid patch in this layout, this function
            fails by assertion.

            \param a_pt A valid patch Point.
        */
        inline bool onLevelBoundary(const Point& a_pt) const;

        /// Iterator End
        inline DataIterator<BoxPartition> end() const;

        /// Get Iterator
        inline DataIterator<BoxPartition> begin() const;// {return DataIterator(*this).begin(); }

        /// Coarsenable Query
        /**
          Checks if the layout is coarsenable by a (possibly anisotropic) refinement ratio.
          This property requires the ProblemDomain to be coarsenable as well as
          boxSize() % a_refRatio == (0,0,....)
        
          \param a_refRatio Refinement ratios
        */
        inline bool coarsenable(const Point& a_refRatio) const;

        /// Compatibility Query
        /**
          Checks if *this and another layout have the same patch configuration.
          This is to say, both layouts have boxes in the same tile positions and on 
          the same processors, but the boxes themselves may be different sizes.
          If true, objects built from these layouts can be processed together in the same
          DataIterator loop.
        
          \param a_layout   A DisjointBoxLayout
        */
        inline bool compatible(const DisjointBoxLayout& a_layout) const;
        
        /// Compatibility Query
        /**
          Checks if this layout can be iterated through using a DataIterator instance.

          \param a_layout   A DataIterator
        */
        inline bool compatible(const DataIterator<BoxPartition>& a_iter) const;

        /// Compatibility Query
        /**
          Checks if a DataIndex is associated with a compatible DataIterator.

          \param a_layout   A DataIterator
        */
        inline bool compatible(const DataIndex<BoxPartition>& a_iter) const;

        //TODO: Improve documentation
        /// Simple Coarsen
        /**
          Coarsens the layout by an input (possibly anisotropic) refinement ratio.
          coarsenable(a_refRatio) must be true or this function will result in an error. 
          The output of this function will always have the same size as the input and 
          will always be compatible in the DataIterator sense.

          The coarsened layout has:
          <code>crseLayout.domain().box() == fineLayout.problemDomain.box().coarsen(refRatio)</code>
          <code>crseLayout[ii] == fineLayout[ii].coarsen(refRatio) </code> for each ii in [0, size())
        
          \param a_refRatio     A refinement ratio
        */
        inline DisjointBoxLayout coarsen(const Point& a_refRatio) const;

        /// Coarsen
        /**
          Creates a new DisjointBoxLayout with a specified box size that has been coarsened by 
          a valid ratio. coarsenable(refRatio) must be true or this function will result in an error.
          The output of this function will in general have a different size than the input and
          is not guaranteed to be compatible in the DataIterator sense.

          \param a_refRatio A refinement ratio
          \param a_boxSize  Size of the output layout's boxes

        */
        inline DisjointBoxLayout coarsen(const Point& a_refRatio, const Point& a_boxSize) const;

        //TODO: Do we implement a maximum box size?
        /// Refine
        /**
          Refines the layout by an input (possibly anisotropic) refinement ratio.
          Unline coarsen(...), this function should always be possible.
          The output of this function will always have the same size as the input and 
          will always be compatible in the DataIterator sense.

          The refined layout has:
          <code>fineLayout.domain().box() == crseLayout.problemDomain.box().refine(refRatio)</code>
          <code>fineLayout[ii] == crseLayout[ii].refine(refRatio)</code> for each ii in [0, size())
          
          \param a_refRatio A refinement ratio
        */
        inline DisjointBoxLayout refine(const Point& a_refRatio) const;

        /// Get Bounding Box
        /**
          Returns the smallest Box containing all of the boxes in the layout.
          Output is guaranteed to be coarsenable by boxSize().
          */
        inline Box boundingBox() const;

        /// Check Radial Symmetry
        /**
          @private
          For debugging mesh refinement procedures
        */
        inline bool radialSymmetry(Point a_origin) const;

        /// Check Mirror Symmetry
        /**
          @private
          For debugging mesh refinement procedures
          */
        inline bool mirrorSymmetry(Point a_origin, int a_coord) const;

        /// Print
        inline void print(std::string a_name = "") const;
        
#ifdef PROTO_CUDA
        // TODO: Used in Chombo4::LevelBoxData. Unclear if we still need this.
        static protoStream_t getCurrentStream()
        {
            static protoStream_t currStream;
            static bool init=false;
            if(!init)
            {
                protoStreamCreate(&currStream);
                init=true;
            }
            return currStream;
        }
        
#endif
        private:
        
        ProblemDomain             m_problemDomain;
        Point                     m_boxSize;  
        //DataIndex<BoxPartition>   m_end;     
        shared_ptr<BoxPartition>  m_partition;
};

/// Stream output for DBL.
inline std::ostream& operator<< (std::ostream& os, const DisjointBoxLayout& a_dbl);
#include "implem/Proto_DisjointBoxLayoutImplem.H"
} //end namespace Proto

#endif
