#pragma once
#ifndef _PROTO_LEVEL_BOUNDARY_REGISTER_
#define _PROTO_LEVEL_BOUNDARY_REGISTER_
#include "Proto_DisjointBoxLayout.H"

namespace Proto
{
    template<typename T, unsigned int C = 1, MemType MEM = MEMTYPE_DEFAULT>
    class BoundaryRegister
    {
        public:
        
        BoxData<T, C, MEM>& operator()(int a_dir, Side::LoHiSide a_side)
        {
            return m_boxes[a_dir][(int)a_side];
        }

        private:
        std::array<std::array<BoxData<T, C, MEM>, 2>, DIM> m_boxes; 
    }; // end class BoundaryRegister

    template<typename T, unsigned int C = 1, MemType MEM = MEMTYPE_DEFAULT>
    class LevelBoundaryRegister
    {
    
        public:
        
        LevelBoundaryRegister(){}
        
        LevelBoundaryRegister(DisjointBoxLayout a_layout, Point a_ghost)
        {
            define(a_layout, a_ghost);
        }
        
        inline void define(DisjointBoxLayout& a_layout, Point a_ghost)
        {
            m_layout = a_layout;
            m_registers.clear();
            m_registers.resize(a_layout.localSize());
            
            Box domainBox = a_layout.domain().box();
            for (auto iter = a_layout.begin(); iter.ok(); ++iter)
            {
                Box b0 = iter.box();
                SideIterator siter;
                for (siter.begin(); siter.ok(); ++siter)
                {
                    for (int dir = 0; dir < DIM; dir++)
                    {   
                        // Doesn't handle periodic boundaries (dealt with via exchange)
                        if (a_layout.domain().isPeriodic(dir)) { continue; }
                        Box boundBox = b0.adjacent(dir, *siter, a_ghost[dir]);
                        if (!domainBox.contains(boundBox))
                        {
                            m_registers[*iter](dir, *siter).define(boundBox);
                        }
                    }
                }
            }
        }

        inline void setRegisterValue(DataIndex a_index, int a_dir, Side::LoHiSide a_side, int a_comp, T a_value)
        {
            auto& data = getRegister(a_index, a_dir, a_side);
            if (!data.box().empty())
            {
                data.setVal(a_value, data.box(), a_comp);
            }

        }
       
        template<typename Func, typename... Srcs>
        inline void setRegisterFunc(DataIndex a_index, int a_dir, Side::LoHiSide a_side, Func& a_func, Srcs... a_srcs)
        {
            auto& data = m_registers[a_index](a_dir, a_side);
            if (!data.box().empty())
            {
                forallInPlace_p(a_func, data, a_srcs...);
            }
        }
        
        inline void setRegisterData(DataIndex a_index, int a_dir, Side::LoHiSide a_side, BoxData<T, C, MEM>& a_src)
        {
            auto& data = m_registers[a_index](a_dir, a_side);
            Box b = data.box() & a_src.box();
            if (!b.empty())
            {
                a_src.copyTo(data);
            }
        }

        inline BoxData<T, C, MEM>& getRegister(DataIndex a_index, int a_dir, Side::LoHiSide a_side)
        {
            return m_registers[a_index](a_dir, a_side);
        }

        inline DisjointBoxLayout& layout() { return m_layout; }

        inline void print()
        {
            for (int ii = 0; ii < m_registers.size(); ii++)
            {
                SideIterator siter;
                for (siter.begin(); siter.ok(); ++siter)
                {
                    for (int dir = 0; dir < DIM; dir++)
                    {   
                        auto& data = m_registers[ii](dir, *siter);
                        if (!data.box().empty())
                        {
                            Proto::pout() << "data address: " << data[0] << std::endl;
                            Proto::pout() << "index (local): " << ii << " | dir: " << dir << " | side: " << *siter << " | box: " << data.box() << std::endl;
                        }
                    }
                }
                
            }
        }

        private:

        std::vector<BoundaryRegister<T, C, MEM>> m_registers; 
        DisjointBoxLayout m_layout;

    }; //end class LevelBoundaryRegister
} // end namespace Proto
#endif //end include guard
