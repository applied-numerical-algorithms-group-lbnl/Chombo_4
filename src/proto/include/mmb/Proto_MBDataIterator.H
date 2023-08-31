#pragma once
#ifndef _PROTO_MB_DATA_ITERATOR_
#define _PROTO_MB_DATA_ITERATOR_

#include "Proto_MBDataIndex.H"
#include "Proto_MBPartition.H"

namespace Proto
{
    class MBDataIterator
    {
        public:
        
        MBDataIterator(const shared_ptr<MBPartition> a_partition);
        inline MBDataIterator& begin();
        inline MBDataIterator& end();
        inline MBDataIterator& operator++();
        inline const MBDataIndex& operator*() const;
        inline bool operator==(const MBDataIndex& a_rhs) const;
        inline bool operator!=(const MBDataIndex& a_rhs) const;
        inline bool operator==(const MBDataIterator& a_rhs) const;
        inline bool operator!=(const MBDataIterator& a_rhs) const;

        private:

        std::shared_ptr<MBPartition> m_partition;
        MBDataIndex m_current;
        MBDataIndex m_end;
    };
#include "implem/Proto_MBDataIterator.H"
} // end namespace Proto
#endif // end include guard
