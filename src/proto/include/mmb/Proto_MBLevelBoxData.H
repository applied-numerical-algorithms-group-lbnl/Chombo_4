#pragma once
#ifndef _PROTO_MB_LEVEL_DATA_
#define _PROTO_MB_LEVEL_DATA_

#include "Proto.H"
#include "Proto_MBDisjointBoxLayout.H"
#include "Proto_MBBoundaryData.H"

namespace Proto
{
    template<typename T, unsigned int C, MemType MEM, Centering CTR=PR_CELL>
    class MBLevelBoxData
    {
        public:
        typedef DataIterator<MBBoxPartition> iterator_t;
        typedef DataIndex<MBBoxPartition> index_t;
        //typedef std::pair<unsigned int, Point> boundsKey_t;
        // PatchID, BlockID, Dir
        typedef std::tuple<Point, unsigned int, Point> boundsKey_t;

        inline MBLevelBoxData(
                MBDisjointBoxLayout& a_layout,
                Point         a_ghost);

        inline MBLevelBoxData(
                MBDisjointBoxLayout&     a_layout,
                std::array<Point, DIM+1> a_ghost);
        
        inline void define(
                MBDisjointBoxLayout&     a_layout,
                std::array<Point, DIM+1> a_ghost);
        
        template<typename Func, typename... Srcs>
        inline void initialize(Func& a_func, Srcs... a_srcs);
        
        template<typename Func, typename... Srcs>
        inline void initConvolve(Func& a_func, Srcs... a_srcs);

        inline BoxData<T, C, MEM>& operator[](index_t a_index);
        inline const BoxData<T, C, MEM>& operator[](index_t a_index) const;
        inline LevelBoxData<T, C, MEM, CTR>& blockData(unsigned int a_block);
        inline const LevelBoxData<T, C, MEM, CTR>& blockData(unsigned int a_block) const;
        inline std::vector<MBBoundaryData<T, C, MEM>>& bounds(
                index_t a_index, Point a_dir);
        inline const std::vector<MBBoundaryData<T, C, MEM>>& bounds(
                index_t a_index, Point a_dir) const;
        inline unsigned int numBlocks() const { return m_layout->numBlocks(); }
        inline iterator_t begin();
        inline iterator_t end();
        private:
    
        std::array<Point, DIM+1> m_ghost;
        MBDisjointBoxLayout*     m_layout;
        std::vector<std::shared_ptr<LevelBoxData<T, C, MEM, CTR>>> m_data;
        std::map<boundsKey_t, std::vector<MBBoundaryData<T, C, MEM>>> m_bounds;
    };
#include "implem/Proto_MBLevelBoxDataImplem.H"
} // end namespace Proto
#endif //end include guard
