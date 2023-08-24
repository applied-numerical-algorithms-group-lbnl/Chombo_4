#pragma once
#ifndef _PROTO_BOX_PARTITION_
#define _PROTO_BOX_PARTITION_

#include "Proto_Morton.H"
#include <unordered_map>
namespace Proto
{
    /// Box Partition (Name Pending?)
    /**
     * This data structure defines the "layoutness" of a collection of patches. 
     * More concretely, Box Partition encapsulates the information of which tile locations
     * in a domain have a patch, which tiles are empty, and which MPI processes are
     * responsible for which patch. 
     *
     * Two objects with the same (or equivalent) Box Partitions are "compatible" with one
     * another in the sense of DataIterator because the two objects in question are guaranteed
     * to have the same number of patches for a given MPI process.
    */
    class BoxPartition
    {
        public: 
       
        /// Constructor
        /**
         * Create a BoxPartition from a "patch domain" and an array of patches
         * in that domain. Here patches are assumed to have a uniform size so as
         * to evenly tile a rectangular layout. By default, patches will be distributed
         * evenly across all available processes.
         *
         * \param a_patchDomain     A ProblemDomain where each Point corresponds to a Box
         * \param a_patches         A vector of Point each representing a patch
        */
        inline BoxPartition(
                const ProblemDomain& a_patchDomain,
                const std::vector<Point>& a_patches);
        
        /// Constructor
        /**
         * Create a BoxPartition from a "patch domain" and an array of patches
         * in that domain. Here patches are assumed to have a uniform size so as
         * to evenly tile a rectangular layout. Patches will be distributed evenly
         * across the range of processes in <code>[a_startProc, a_endProc)</code>
         *
         * \param a_patchDomain     A ProblemDomain where each Point corresponds to a Box
         * \param a_patches         A vector of Point each representing a patch
        */
        inline BoxPartition(
                const ProblemDomain& a_patchDomain,
                const std::vector<Point>& a_patches,
                unsigned int a_startProc,
                unsigned int a_endProc);
        
        /// Define
        /**
         * Create a BoxPartition from a "patch domain" and an array of patches
         * in that domain. Here patches are assumed to have a uniform size so as
         * to evenly tile a rectangular layout. Patches will be distributed evenly
         * across the range of processes in <code>[a_startProc, a_endProc)</code>
         *
         * \param a_patchDomain     A ProblemDomain where each Point corresponds to a Box
         * \param a_patches         A vector of Point each representing a patch
        */
        inline void define(
                const ProblemDomain& a_patchDomain,
                const std::vector<Point>& a_patches,
                unsigned int a_startProc,
                unsigned int a_endProc);
        
        /// Load Balance
        /**
         * Rebalance the load of this by distributing the input patches evenly 
         * across the range of processes in <code>[a_startProc, a_endProc)</code>
         * The input patches must be contained in the existing ProblemDomain of this
         *
         * \param a_patches         A vector of patches represented as Points
         * \param a_patchDomain     A ProblemDomain where each Point corresponds to a Box
         * \param a_patches         A vector of Point each representing a patch
        */
        inline void loadBalance(
                std::vector<Point>& a_patches,
                unsigned int a_startProc,
                unsigned int a_endProc);

        /// Load Assign
        /**
         * Manually assign the load of this using the syntax
         * <code> loadAssign(patches, P0, N0, P1, N1, ..., Pk, Nk) </code>
         * where process Pi is assigned Ni boxes. The sum over all Ni should be
         * equal to the length of the patch vector.
         * 
         * \param a_patches         A vector of patches represented as Points
         * \param a_patchDomain     A ProblemDomain where each Point corresponds to a Box
         * \param a_patches         A vector of Point each representing a patch
        */
        template<typename... vals>
        inline void loadAssign(
                std::vector<Point>& a_patches,
                vals... args);

        /// Load Assign
        /**
          * Manually assign the load of this layout. The first input is a vector of
          * <code>pair<int, unsigned int></code> where the first entry of each pair is a 
          * process and the second is the number of boxes to be assigned to that process.
          * The second input is the set of patches themselves represented as Points.
          *
          * \param a_assignment     A vector of proc-assignment pairs
          * \param a_patches        A vector of patches represented as Points
        */
        inline void loadAssign(
                std::vector<std::pair<int, unsigned int>>& a_assignment,
                std::vector<Point>& a_patches);
       
        /// Load Assign
        /**
          * Manually (re)assign the load of this layout. The input is a vector of
          * <code>pair<int, unsigned int></code> where the first entry of each pair is a 
          * process and the second is the number of boxes to be assigned to that process.
          * The patches themselves are unchanged by this function.
          *
          * \param a_assignment     A vector of proc-assignment pairs
          * \param a_patches        A vector of patches represented as Points
        */
        inline void loadAssign(
                std::vector<std::pair<int, unsigned int>>& a_assignment);
        
        /// Compatibility Query
        /**
         * Determine if this has the same layout structure as another BoxPartition.
        */
        inline bool compatible(const BoxPartition& a_rhs);
   
        /// Number of Processes
        /**
         * Query the number of processes assigned to this BoxPartition
        */
        inline unsigned int numProcs() const;

        /// Number of Boxes (Global)
        /**
         * Query the number of patches in this BoxPartition across all processes
        */
        inline unsigned int numBoxes() const;
    
        /// Number of Boxes (Local)
        /**
         * Query the number of patches in this BoxPartition on a specified process
        */
        inline unsigned int numBoxes(unsigned int a_proc) const;

        /// Access Partition
        /**
         * Access the vector of pair<Point, unsigned int> which maps each patch
         * (represented as a Point) to a process number
        */
        inline const std::vector<std::pair<Point, unsigned int>>& partition() const
        { return m_partition; }

        /// Get Processor Starting Index
        /**
         * Get the index in global array at which a specified process's data begins.
         * Returns this->size() if the processor is not assigned to this.
        */
        inline unsigned int procStartIndex(unsigned int a_proc) const;
        
        /// Get Processor Ending Index
        /**
         * Get the index in global array at which a specified process's data ends
         * Returns this->size() if the processor is not assigned to this.
        */
        inline unsigned int procEndIndex(unsigned int a_proc) const;
        
        /// Get Patch Index
        /**
         * Get the index in global array corresponding to a specified patch.
         * Returns this->size() if the patch is not in this partition's domain.
        */
        inline unsigned int find(Point a_pt);
       
        inline const ProblemDomain& domain() const {return m_patchDomain; }

        /// Print
        inline void print() const;
        
        private:
        
        template<typename... vals>
        inline void unpack(
                std::vector<Point>& a_patches,
                unsigned int        a_globalIndex,
                int                 a_proc,
                unsigned int        a_num,
                vals...             a_args);

        inline void assign(
                std::vector<Point>& a_patches,
                unsigned int        a_globalIndex,
                int                 a_proc,
                unsigned int        a_num);

        ProblemDomain m_patchDomain; ///< Domain in patch space
        std::unordered_map<uint64_t, int> m_indexMap; ///< Maps Morton index to global index
        std::unordered_map<unsigned int, std::pair<unsigned int, unsigned int>>      m_procMap; ///< Maps processor number to global index
        //std::vector<unsigned int> m_procIndices; ///< Maps processor number to global index space
        //unsigned int m_startProc; ///< The first proc on which this partition has data
        std::vector<std::pair<Point, unsigned int>> m_partition; ///< Maps each patch to a proc
    };

#include "implem/Proto_BoxPartitionImplem.H"
} // end namespace Proto
#endif //end include guard
