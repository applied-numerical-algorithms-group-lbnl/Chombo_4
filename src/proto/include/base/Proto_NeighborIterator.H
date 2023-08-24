#pragma once
#ifndef _PROTO_NEIGHBOR_ITERATOR__
#define _PROTO_NEIGHBOR_ITERATOR__
#include "Proto_DataIndex.H"
#include "Proto_DisjointBoxLayout.H"
namespace Proto
{
    class NeighborIterator
    {
        protected:
        
        DisjointBoxLayout m_layout;
        Box m_bx;
        BoxIterator m_bxit;
        Point m_shift;
        DataIndex<BoxPartition> m_current;
        DataIndex<BoxPartition> m_end;
        ProblemDomain m_dProbC;
        
        public:
        
        /// Constructor takes DBL, and Box that we are intersecting with the DBL.
        inline NeighborIterator(const DisjointBoxLayout& a_dbl, const Box& a_bx);
        //~NeighborIterator();

        /// Iterator functions begin, end, ++.
        inline NeighborIterator& begin();
        inline const DataIndex<BoxPartition>& end() const{return m_end;};
        inline bool ok() const {return (m_current != m_end);}
        inline NeighborIterator& operator++();

        /// Dereference operator returns the DataIndex corresponding to the current iterate.
        inline const DataIndex<BoxPartition>& operator*() const{return m_current;};
        /// The Point corresponding to the Box for that index. Undefined if invalid.
        inline Point point() const;
        /// The region to be copied from corresponding to that index. Undefined if invalid.
        inline Box srcBox() const;
        /// The region to be copied from corresponding to that index. Undefined if invalid.
        inline Box destBox() const;
        /// The amount by which the srcBox is shifted to obtain destBox a_bx.
        inline Point shift();
        // Processor ID for current iterate.
        inline int procID(){return m_layout.procID(m_current);}; 
    };
} // end namespace Proto
#include "implem/Proto_NeighborIteratorImplem.H"
#endif // end include guard
