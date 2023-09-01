#pragma once
#ifndef _PROTO_DATA_ITERATOR_
#define _PROTO_DATA_ITERATOR_
#include <vector>
#include "Proto_DataIndex.H"
//#include "Proto_BoxPartition.H"

namespace Proto
{
    // Implementation Details
    // DataIterator is templated on an object which defines "layout-ness".
    // In Proto, this was historically encoded in the deprecated DBLInternals object
    // (which has been refactored into BoxPartition).
    // 
    // Valid implementations of the class P must implement the following interface:
    //
    // class P {
    //      public:
    //      unsigned int numBoxes();                    // total number of patches in the layout
    //      unsigned int procStartIndex(int a_proc);    // return the global index of the first patch on a specified processor
    //      unsigned int procEndIndex(int a_proc);      // return the global index of the first patch on a specified processor
    //      bool compatible(const P& a_p);              // determine if two instances of P are compatible (e.g. if they have the layout-ness)
    // };
    
    /// Distributed Data Iterator
    /**
        Iterate through the local contents of a data layout specified by the template class P.
        Unlike normal iterators, DataIterator<P> dereferences into a DataIndex<P> which can
        be used to index into various data holders that are specifically designed to work
        in tandem with it. 

        Supported objects:
        DisjointBoxLayout   (P = BoxPartition)
        LevelBoxData        (P = BoxPartition)

        Proto objects are specifically designed to use STL-like iteration syntax which 
        makes the specific implementation of DataIterator (and it's template argument) 
        invisible to the end user. See the code below for sample usage:

        @code
        DisjointBoxLayout layout(...);
        LevelBoxData<double> data(layout, ...);

        for (auto iter : layout)
        {
            // iter is a DataIndex object which can be used to index into data holders
            Box b = layout[iter];
            BoxData<double>& patch = data[iter];
            // etc
        }
        @endcode
    */
    template<typename P>
    class DataIterator
    {
        public:
       
        /// Manual Constructor
        /**
            Construct a DataIterator from an existing DisjointBoxLayout.
            This constructor is mostly used internally. Users can rely on 
            the STL style iteration syntax of container objects
            (e.g. LevelBoxData or DisjointBoxLayout)
        */ 
        inline DataIterator(std::shared_ptr<P> a_partition);
        
        /// Set To Start
        /**
            Modifies the *this and returns it
        */
        inline DataIterator<P>& begin();

        /// Set To End 
        /**
            modifies *this and returns it
        */
        inline DataIterator<P>& end();
        
        /// Continue Query
        /**
            Shortcut to check if *this has reached its end.
        */
        inline bool ok() const;

        /// Increment 
        inline DataIterator<P>& operator++();
       
        /// Equality
        inline bool operator==(const DataIterator<P>& a_rhs) const { return m_current == a_rhs.m_current; }
        
        /// Equality
        inline bool operator==(const DataIndex<P>& a_rhs) const { return m_current == a_rhs; }
        
        /// Inquality
        inline bool operator!=(const DataIterator<P>& a_rhs) const { return m_current != a_rhs.m_current; }
        
        /// Inquality
        inline bool operator!=(const DataIndex<P>& a_rhs) const { return m_current != a_rhs; }

        /// Dereference
        /**
            Unlike normal iterators, dereferencing a DataIterator returns a DataIndex.
            DataIndex is accepted by all "level scope" data holders (e.g. LevelBoxData or DisjointBoxLayout).
        */
        inline const DataIndex<P>& operator*() const;

        inline bool compatible(const P& a_partition) const;

        inline bool compatible(const DataIndex<P>& a_index) const;
        
        inline bool compatible(const DataIterator<P>& a_iter) const;
        
        protected:
        
        std::shared_ptr<P> m_partition;
        DataIndex<P> m_current;
        DataIndex<P> m_end;
    };


#include "implem/Proto_DataIteratorImplem.H"
} // end Proto namespace.
#endif

