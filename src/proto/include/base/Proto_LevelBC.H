#pragma once

#ifndef _PROTO_LEVEL_BC_
#define _PROTO_LEVEL_BC_

#include "Proto_LevelBoundaryRegister.H"

namespace Proto
{
    enum BCType
    {
        DIRICHLET,
        NEUMANN,
    };

    template<BCType TYPE, typename T, unsigned int C = 1, MemType MEM = MEMTYPE_DEFAULT>
    class LevelBC
    {
        public:
            
        LevelBC(DisjointBoxLayout& a_layout, Point a_ghost)
        {
            m_register.define(a_layout, a_ghost);
            for (auto iter = a_layout.begin(); iter.ok(); ++iter)
            {
                SideIterator siter;
                for (siter.begin(); siter.ok(); ++siter)
                {
                    for (int dir = 0; dir < DIM; dir++)
                    {
                        for (int comp = 0; comp < C; comp++)
                        {
                            m_register.setRegisterValue(*iter, dir, *siter, comp, 0);
                        }
                    }
                }
            }
        }

        inline void set(int a_dir, Side::LoHiSide a_side, int a_comp, T a_value)
        {   
            for (auto iter = m_register.layout().begin(); iter.ok(); ++iter)
            {
                m_register.setRegisterValue(*iter, a_dir, a_side, a_comp, a_value);
            }
        }
        
        template<typename Func, typename... Srcs>
        inline void set(int a_dir, Side::LoHiSide a_side, Func& a_func, Srcs... a_srcs)
        {   
            for (auto iter = m_register.layout().begin(); iter.ok(); ++iter)
            {
                m_register.setRegisterFunc(*iter, a_dir, a_side, a_func, a_srcs...);
            }
        }
       
        //TODO: Need to think more about how/when/where the BCs are actually applied
        //      We must apply them at the BoxData scope, but we also must communicate
        //      The location of the BoxData within the domain (through a DataIndex)
        //      We may also want an overload that applies the BC on a whole level
        //      through fixing up LPhi. More thought is needed. -CLG
        inline void apply(BoxData<T, C, MEM>& a_data, DataIndex a_index)
        {
            SideIterator siter;
            for (siter.begin(); siter.ok(); ++siter)
            {
                for (int dir = 0; dir < DIM; dir++)
                {
                    auto& boundData = m_register.getRegister(a_index, dir, *siter);
                    if (boundData.box().empty()) { continue; }
                    //std::cout << "Applying BC | index: " << (int)a_index << " | dir: " << dir << " | side: " << *siter << std::endl;
                    switch (TYPE)
                    {
                        case NEUMANN: 
                        {
                            Point shift = Point::Zeros();
                            if (*siter == Side::Lo)
                            {
                                shift = Point::Basis(dir, Side::Hi);
                            }
                            boundData.copyTo(a_data, boundData.box(), shift);
                            break;
                        }
                        case DIRICHLET:
                        {
                            // need to extrapolate data to ghost region such that BC at face is satisfied.
                            // interpolate using known dbc at boundary plus values in the interior
                            MayDay<void>::Abort("LevelBC::apply | Error: Dirichlet BCs are not yet implemented.");
                        }
                    }
                }
            }
        }
            

        private:

        LevelBoundaryRegister<T, C, MEM> m_register;

    }; // end class LevelBC
} // end namespace Proto
#endif // end include guard
