MBDataIterator::MBDataIterator(const shared_ptr<MBPartition> a_partition)
{
    m_partition = a_partition;
    m_current = MBDataIndex(a_partition, a_partition->procStartIndex(procID()));
    m_end = MBDataIndex(a_partition, a_partition->numBoxes());
}

inline MBDataIterator& MBDataIterator::begin()
{
    m_current = MBDataIndex(m_partition, m_partition->procStartIndex(procID()));
    return *this;
}
inline MBDataIterator& end();
inline MBDataIterator& operator++();
inline const MBDataIndex& operator*() const;
inline bool operator==(const MBDataIndex& a_rhs) const;
inline bool operator!=(const MBDataIndex& a_rhs) const;
inline bool operator==(const MBDataIterator& a_rhs) const;
inline bool operator!=(const MBDataIterator& a_rhs) const;
