#pragma once
#ifndef _PROTO_MB_DATA_INDEX_
#define _PROTO_MB_DATA_INDEX_

namespace Proto
{
    class MBDataIterator;
    class MBDataIndex
    {
        public:
        
        MBDataIndex( const shared_ptr<MBPartition> a_partition, unsigned int a_index)
            : m_partition(a_partition), m_index(a_index);
        
        inline bool operator==(const MBDataIndex& a_rhs);
        inline bool operator!=(const MBDataIndex& a_rhs);
        inline bool operator==(const MBDataIterator& a_rhs);
        inline bool operator!=(const MBDataIterator& a_rhs);
        inline unsigned int local() const; 
        inline unsigned int global() const;
        inline operator int() const { return local(); }

        private:

        std::shared_ptr<MBPartition> m_partition;
        unsigned int m_index;

    }
#include "implem/Proto_MBDataIndex.H"
} // end namespace Proto
#endif // end include guard
