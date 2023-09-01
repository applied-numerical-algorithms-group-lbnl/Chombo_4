#pragma once
#ifndef _PROTO_INTERVAL_
#define _PROTO_INTERVAL_

namespace Proto {

    /// Component-Space Interval
    /**
      Defines a subspace of tensor components.
      Used chiefly for copying and and linearized buffer IO.
      Usually the user will not need to explicitly create a CInterval.
      See the documentation for BoxData::copyTo for an example. 
      */
    class CInterval 
    {
        public:
            /// Bounds Constructor
            /**
              Builds the interval:
              {{a_c0, a_c1},{a_d0, a_d1},{a_e0, a_e1}}
            */
            inline CInterval(
                    unsigned int a_c0,      unsigned int a_c1,
                    unsigned int a_d0 = 0,  unsigned int a_d1 = 0,
                    unsigned int a_e0 = 0,  unsigned int a_e1 = 0);
            
            /// List Constructor
            /**
              Build a CInterval using initializer list syntax.

              Example:
              @code
              CInterval I0{{c0, c1},{d0, d1},{e0, e1}};
              // OR
              CInterval I1({{c0, c1},{d0, d1},{e0, e1}});
              // Or, for a single component index:
              CInterval I2{c0, c1};
              // OR
              CInterval I3({c0, c1});
            @endcode
            */
            inline CInterval(std::initializer_list<std::initializer_list<unsigned int>> a_lst);
   
            /// Equality
            /**
            */
            inline bool operator==(CInterval& a_rhs);
            
            /// Lower Bound
            /**
              Retrieve the lower bound along the given component axis

              \param a_comp   A component axis in [0,3)
              */
            inline unsigned int low(unsigned int a_comp) const;

            /// Upper Bound
            /**
              Retrieve the upper bound component
              for component axis <code> a_comp </code>

              \param a_comp   A component axis in [0,3)
              */
            inline unsigned int high(unsigned int a_comp) const;

            /// Contains Query
            /**
              Query if <code>*this</code> contains <code>index</code>
              of component <code>comp</code>

              \param a_index  An index
              \param a_comp   A component axis in [0,3)
              */
            inline bool contains(unsigned int a_index, unsigned int a_comp) const;

            /// Size Query
            /**
              Returns the number of components in <code> *this </code>
              along the component axis <code> a_comp </code>

              \param a_comp   A component axis in [0,3)
              */
            inline unsigned int size(unsigned int a_comp) const;

            /// Print
            inline void print() const;

        private:
            unsigned int m_comps[3][2]; ///< bounds of the interval
    };

    /// CInterval IOStream Operator
    inline std::ostream& operator<<(std::ostream& a_os, const CInterval& a_int);


#include "implem/Proto_CIntervalImplem.H"
} // end namespace Proto
#endif // end include guard
