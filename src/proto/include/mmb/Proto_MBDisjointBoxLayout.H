#pragma once
#ifndef _PROTO_MB_BOX_LAYOUT_
#define _PROTO_MB_BOX_LAYOUT_

#include "Proto.H"
#include "Proto_MBProblemDomain.H"
#include "Proto_MBBoxPartition.H"

namespace Proto
{
    class MBDisjointBoxLayout
    {
        public:
        typedef DataIterator<MBBoxPartition> iterator_t;
        typedef DataIndex<MBBoxPartition> index_t;
        
        inline MBDisjointBoxLayout(
                const MBProblemDomain& a_domain,
                Point a_boxSize);

        inline MBDisjointBoxLayout(
                const MBProblemDomain& a_domain,
                const std::vector<Point>& a_boxSizes);
        
        inline MBDisjointBoxLayout(
                const MBProblemDomain& a_domain,
                const std::vector<MBPatchID_t>& a_patches,
                const std::vector<Point>& a_boxSizes);
        
        inline void define(
                const MBProblemDomain& a_domain,
                const std::vector<Point>& a_boxSizes);
        
        inline void define(
                const MBProblemDomain& a_domain,
                const std::vector<MBPatchID_t>& a_patches,
                const std::vector<Point>& a_boxSizes);
        
        inline MBProblemDomain domain() const; 
        inline unsigned int numBlocks() const;
        inline unsigned int numBoxes() const;
        inline unsigned int numBoxes(unsigned int a_proc) const;
        inline unsigned int numProcs() const;
        inline DataIndex<BoxPartition> blockIndex(const index_t& a_index) const;
        inline Box operator[](const index_t& a_index);
        inline iterator_t begin();
        inline iterator_t end();
        inline Point point(const index_t& a_index) const; 
        inline Box box(const index_t& a_index) const; 
        inline unsigned int block(const index_t& a_index) const; 
        inline DisjointBoxLayout& layout(unsigned int a_block);

        private:
        MBProblemDomain                 m_domain;
        std::vector<DisjointBoxLayout>  m_layouts;
        std::shared_ptr<MBBoxPartition> m_partition;
    };
#include "implem/Proto_MBDisjointBoxLayoutImplem.H"
} // end namespace Proto
#endif //end include guard
