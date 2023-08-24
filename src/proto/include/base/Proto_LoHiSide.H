#ifndef _Proto_LOHISIDE_H_
#define _Proto_LOHISIDE_H_

namespace Proto
{
    /*
        TODO: I have no clue what the "Side" class is doing aside from obfuscating syntax.
              This file should contain a Side enum, a couple namespace functions, and SideIterator.
              I'm not changing this because of Chombo4 dependence -CLG
    */
   
    // forward declarations
    class SideIterator;
    
    /// Side
    /**
        Abstraction for high / low side of a grid cell. 
    */
    class Side
    {
    public:
        /// Side Enum
        enum LoHiSide
        {
            Invalid=-1, ///< Invalid 
            Lo=0, ///< Low side; normal is in negative coordinate direction
            Hi, ///< High side; normal is in the positive coordinate direction
            NUMSIDES ///< Number of sides
        };

        /// Constructor
        Side() {}

        /// Destructor
        ~Side() {}
        
        /// Get Begin Iterator
        inline static SideIterator begin();
        
        /// Get End Iterator
        inline static SideIterator end(); 
    };
    
    /// Flip Side
    /**
      Utility function for flipping between Side::Lo and Side::Hi

      \param a_side   Side::Lo or Side::Hi
      */
    inline Side::LoHiSide flip(const Side::LoHiSide& a_side) { return (a_side == Side::Lo) ? Side::Hi : Side::Lo; }

    /// Sign
    /**
        Converts a Side to +1 (Hi) or -1 (Lo);
    */
    inline int sign(const Side::LoHiSide& a_side) { return (a_side == Side::Lo) ? -1 : 1; }

    /// Side Iterator
    /**
        Convenience iterator to allow iteration through Lo, Hi.
    */
    class SideIterator
    {
        public:
        /// Constructor
        SideIterator(int a_start = 0) { m_current = a_start; }

        /// Destructor
        ~SideIterator() {}

        /// Begin
        /**
            Set the iterator to the first side.
        */
        inline void begin() { m_current = 0; }

        /// Increment
        /**
            Advance the iterator to the next face.
        */
        inline void operator++() { ++m_current; }
        
        /// Equality Operator
        inline bool operator==(const SideIterator& a_iter) { return m_current == a_iter.m_current; }
        
        /// Inequality Operator
        inline bool operator!=(const SideIterator& a_iter) { return m_current != a_iter.m_current; }

        /// Dereference
        inline Side::LoHiSide operator*() { return this->operator()(); }

        /// Get Value
        /**
            Returns the current Side
        */
        inline Side::LoHiSide operator () () const
        {
            switch (m_current)
            {
                case 0:
                    return Side::Lo;
                case 1:
                    return Side::Hi;
                default:
                    return Side::Invalid;
            }
        }

        /// Continue Query
        /**
           Returns <tt>true</tt> if there is a next side, <tt>false</tt> if not.
        */
        inline bool ok() const { return ((m_current > -1) && (m_current < Side::NUMSIDES)); }
        
        private:
        int m_current;
    };

    SideIterator Side::begin() { return SideIterator(0); }
    SideIterator Side::end() { return SideIterator(2); }

    }
#endif
