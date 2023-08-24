#pragma once
#ifndef _PROTO_DATA_INDEX_
#define _PROTO_DATA_INDEX_
#include <vector>
#include <iostream>
#include "Proto.H"
//#include "Proto_BoxPartition.H"
namespace Proto
{
    class DisjointBoxLayout;
    template<typename PP> class DataIterator;

    /// Data Index
    /**
     * Data Index is a special object used in tandem with DataIterator.
     * See DataIterator for further documentation.
     *
     * \tparam P    A "Partition" object which defines "layoutness" of the iterable object (see BoxPartition)
    */
    template<typename P>
    class DataIndex
    {
        friend class DisjointBoxLayout;
        template<typename PP> friend class DataIterator;
        template<typename PP>
        friend std::ostream& operator<<(std::ostream& a_os, const DataIndex<PP>& a_index);

        public:
        
        /// Default Constructor
        DataIndex(){};

        /// Explicit Constructor
        inline DataIndex(const std::shared_ptr<P> a_partition, unsigned int a_index)
        {
            PROTO_ASSERT(a_partition != nullptr,
                    "DataIndex::Constructor | Error: Partition object pointer is null.");
            m_partition =  a_partition;
            m_currentInt = a_index;     
        };
        
        /// Destructor
        ~DataIndex(){};
        
        /// Equality Operator
        inline bool operator==(const DataIndex<P>& a_di) const {return (m_partition == a_di.m_partition) && (m_currentInt == a_di.m_currentInt);};
       
        /// Equality Operator (Overload)
        /**
         * Provided to increase expressiveness and backwards compatibility
        */
        inline bool operator==(const DataIterator<P>& a_iter) const;
        
        /// Inequality Operator
        inline bool operator!=(const DataIndex<P>& a_di) const {return !(*this==a_di);};
        
        /// Inequality Operator (Overload)
        /**
         * Provided to increase expressiveness and backwards compatibility
        */
        inline bool operator!=(const DataIterator<P>& a_iter) const;

        /// Integer Cast
        /**
            Data Index implicitly casts to the LOCAL index of the iterable.
            This behavior is used internally to increase code clarity.
        */
        inline operator int() const { return local(); };
    
        /// Local Value
        /**
         *  Returns the integer value associated with the LOCAL (in the sense of MPI process) index of the iterate.
        */
        inline int local() const { return m_currentInt - m_partition->procStartIndex(procID()); } 
        
        /// Local Value
        /**
         *  Returns the integer value associated with the GLOBAL (in the sense of MPI process) index of the iterate.
        */
        inline int global() const { return m_currentInt; } 

        /// Compatibility Query
        /**
         * Returns true if the partition associated with <code>*this</code> is compatible 
         * with that associated with another index.
        */
        bool compatible(const DataIndex<P>& a_index) const
        { return m_partition->compatible(*a_index.m_partition); }
        
        /// Compatibility Query
        /**
         * Returns true if the partition associated with <code>*this</code> is compatible 
         * with another partition.
        */
        bool compatible(const P& a_partition) const
        { return m_partition->compatible(a_partition); }

        protected:

        shared_ptr<P> m_partition;
        int m_currentInt;

    }; // end class DataIndex;
    
    /// IO Stream Operator
    template<typename P>
    inline std::ostream& operator<<(std::ostream& os, const DataIndex<P>& a_di)
    {
        os << "DataIndex: partition pointer: " << a_di.m_partition;
        os <<  " | global(local) index: " << a_di.global();
        os << "(" << a_di.local() << ")" << std::endl;
        return os;
    }
}
#endif
