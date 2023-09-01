#ifndef _PROTO_POWERITOI_H_
#define _PROTO_POWERITOI_H_

#include <cassert>
#include <cmath>

namespace Proto {

    /// Template Based Integer Exponentiation
    template <unsigned int P>
        inline int ipow(int M){return M*ipow<P-1>(M);}

    template <>
        inline int ipow<0>(int M){return 1;}

    // Dynamic version
    inline int ipow(int a_base, unsigned int a_exp)
    {
        int result = 1;
        for(;;)
        {
            if (a_exp & 1){result *= a_base;}
            a_exp >>= 1;
            if (!a_exp){break;}
            a_base *= a_base;
        }
        return result;
    }
} //end Proto namespace

#endif //end include guard
