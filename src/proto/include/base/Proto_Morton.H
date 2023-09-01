#pragma once
#ifndef _PROTO_MORTON_H_
#define _PROTO_MORTON_H_
#include <array>
#include <cstdint>
#include <algorithm>
#include "Proto_Point.H"

#if DIM==2
#define LOGSIZE 16 // 21 // DIM*LOGSIZE must be less than the size (in bits) of (uint64_t).
#define MORTONSIZE 65536 // 2097152 // = 2^LOGSIZE.
#endif
#if DIM==3
#define LOGSIZE 10 // DIM*LOGSIZE must be less than the size (in bits) of (uint64_t).
#define MORTONSIZE 1024 // = 2^LOGSIZE.
#endif
#if DIM==4
#define LOGSIZE 16 // DIM*LOGSIZE must be less than the size (in bits) of (uint64_t).
#define MORTONSIZE 65536 // = 2^LOGSIZE.
#endif
#if DIM==5
#define LOGSIZE 12 // DIM*LOGSIZE must be less than the size (in bits) of (uint64_t).
#define MORTONSIZE 4096 // = 2^LOGSIZE.
#endif
#if DIM==6
#define LOGSIZE 10 // DIM*LOGSIZE must be less than the size (in bits) of (uint64_t).
#define MORTONSIZE 1024 // = 2^LOGSIZE.
#endif

namespace Proto
{
    /// Morton Indexer
    /** 
       Utility class for computing the Morton Index of a DIM-tuple corresponding to the 
       bits of each element of the tuple.
    */
    class Morton
    {
        public:
        
        /// Compute Morton Index
        /**
          Compute the Morton index of a Point.

          \param a_pt     A Point
          \param a_morton A Morton indexer
          */
        inline static uint64_t index(const Point& a_pt);
        
        /// Morton Sort
        /**
            Sort a vector of Point in place by Morton index.
            \param a_pts    A vector of Points to be sorted
        */
        inline static void sort(vector<Point>& a_pts);
        
        private:

        std::array<vector<uint64_t>,DIM> m_morton1D;
        
        // Public construction is not allowed. Morton is a singleton.
        inline Morton();
       
        // Access the singleton Morton object
        inline static Morton& instance();
        
        // Sorting condition
        inline static bool compareSecond(pair<Point,uint64_t> a_a,pair<Point,uint64_t> a_b)
        {
            return a_a.second < a_b.second;
        }
    };
#include "implem/Proto_MortonImplem.H"
}// end Proto namespace.
#endif
