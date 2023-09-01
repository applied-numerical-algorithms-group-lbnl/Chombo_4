#pragma once
#ifndef _PROTO_MB_MAP_
#define _PROTO_MB_MAP_

#include "Proto_MBDisjointBoxLayout.H"
#include "Proto_MBLevelBoxData.H"

namespace Proto
{
    template<typename T, MemType MEM>
    class MBMap
    {
        public:
      
        /// Constructor
        /**
            Create a MBMap and allocate space for data members

            \param a_layout  The MMB domain
            \param a_order   Order of accuracy for differentiation (default: 4)
        */
        MBMap(const MBDisjointBoxLayout& a_layout, unsigned int a_ghost, unsigned int a_order = 4);
        
        /// Compute Metric Terms
        /**
            Compute the map, cofactor matrices, and Jacobian from a user defined function
            a_func must have the signature:
            
            template<typename T, MemType MEM>
            void a_func(Point& a_pt, Var<T, DIM, MEM>& a_X, unsigned int a_block, ...)

            Where <code>a_block</code> is the block index and ... are any optional arguments
            needed and which conform to the input requirements of Proto::forall. 
        */
        template<typename Func, typename... Srcs> 
        inline void compute(const Func& a_func, Srcs&&... a_srcs);

        /// Map Access
        inline const BoxData<T, DIM, MEM>& map(unsigned int a_blockIndex, DataIndex a_patchIndex) const
        inline const BoxData<T, DIM, MEM>& map(const MBDataIndex& a_index) const
        inline const LevelBoxData<T, DIM, MEM>& map(unsigned int a_blockIndex) const
        
        /// Cofactors Access
        inline const BoxData<T, DIM, MEM>& cofactors(unsigned int a_blockIndex, DataIndex a_patchIndex, int a_dir) const
        inline const BoxData<T, DIM, MEM>& cofactors(const MBDataIndex& a_index, int a_dir) const
        inline const LayoutData<std::array<std::array<BoxData<T, DIM, MEM>, DIM>>>& cofactors(unsigned int a_blockIndex) const

        /// Jacobian Access
        inline const BoxData<T, DIM, MEM>& jacobian(unsigned int a_blockIndex, DataIndex a_patchIndex) const
        inline const BoxData<T, DIM, MEM>& jacobian(const MBDataIndex& a_index) const
        inline const LevelBoxData<T, DIM, MEM>& jacobian(unsigned int a_blockIndex) const

        private:
        
        MBGraph* m_graph;
        unsigned int m_order;

        MBLevelBoxData<T, DIM, MEM>     m_X;
        //MBLevelBoxData<T, DIM*DIM, MEM> m_NT;
        MBLayoutData<std::array<std::array<BoxData<T, DIM, MEM>, DIM>>> m_NT;
        MBLevelBoxData<T, DIM, MEM>     m_J;


        bool m_defined;  ///< Has data been allocated
        bool m_computed; ///< Has data been computed
    }
} // end namespace Proto

#include "implem/Proto_MBMapImplem.H"
#endif // end include guard
