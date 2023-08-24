#pragma once
#ifndef _PROTO_INTERP_STENCIL_TP_
#define _PROTO_INTERP_STENCIL_TP_
#include "Proto.H"

namespace Proto
{
    template<typename T>
    class InterpStencilTP
    {
        public:

        InterpStencilTP(){};
        InterpStencilTP(int a_order, int a_refRatio);
        void define(int a_order, int a_refRatio);
        
        template<unsigned int C, MemType MEM=MEMTYPE_DEFAULT>
        inline void apply(BoxData<T, C, MEM>& a_output, BoxData<T, C, MEM>& a_input, T a_scale = 1.0);
        inline Point ghost() const;        
        private:

        Stencil<T> build1D5pt(array<double,5>& a_coefs, int a_dir);

        std::array<InterpStencil<T>, DIM> m_interp;

    }; //end class InterpStencilTP
#include "implem/Proto_InterpStencilTPImplem.H"
} // end Proto namespace
#endif //end include guard
