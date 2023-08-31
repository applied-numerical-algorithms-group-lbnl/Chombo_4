#pragma once
#ifndef _PROTO_MB_BOX_PARTITION_
#define _PROTO_MB_BOX_PARTITION_

#include "Proto.H"
#include "Proto_MBProblemDomain.H"
#include <unordered_map>
namespace Proto
{
    typedef std::pair<Point, unsigned int> MBPatchID_t;
    
    class MBBoxPartition
    {
        public:

        inline MBBoxPartition(
                const MBProblemDomain& a_domain,
                const std::vector<MBPatchID_t>& a_patches);
        inline void define(
                const MBProblemDomain& a_domain,
                const std::vector<MBPatchID_t>& a_patches);
        inline void loadBalance();
        inline bool compatible(const MBBoxPartition& a_rhs);
        inline unsigned int numProcs() const;
        inline unsigned int numBoxes() const;
        inline unsigned int numBoxes(unsigned int a_proc);
        inline unsigned int numBlocks() const;
        inline uint64_t procStartIndex(unsigned int a_proc) const;    
        inline uint64_t procEndIndex(unsigned int a_proc) const;    
        inline uint64_t find(Point a_patch, unsigned int a_block);
        inline const std::vector<std::pair<MBPatchID_t, unsigned int>>& partition() const
        { return m_partition; }
        inline std::shared_ptr<BoxPartition> blockPartition(unsigned int a_block) const;
        inline void print();
    
        private:

        inline void assign(
                std::vector<Point>& a_patches,
                unsigned int        a_globalIndex,
                int                 a_proc,
                unsigned int        a_num);

        inline void buildLocalMaps();
        inline void buildGlobalMaps();

        MBProblemDomain m_patchDomain;
        std::map<MBPatchID_t, uint64_t> m_indexMap; ///< Maps local to global
        std::unordered_map<unsigned int, std::pair<uint64_t, uint64_t>> m_procMap; ///< Maps proc to global indices
        std::vector<std::shared_ptr<BoxPartition>> m_blockPartitions;
        std::vector<std::pair<MBPatchID_t, unsigned int>> m_partition;
    };
#include "implem/Proto_MBBoxPartitionImplem.H"
} // end namespace Proto

#endif //end include guard
