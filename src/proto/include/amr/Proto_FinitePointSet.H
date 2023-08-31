#pragma once
#ifndef _PROTO_FINITE_POINT_SET_
#define _PROTO_FINITE_POINT_SET_

#include "Proto_Box.H"
#include "Proto_ProblemDomain.H"

namespace Proto
{
    
    class FinitePointSet
    {
        public:
        
        // no reason to be doing these things
        FinitePointSet(const FinitePointSet& a_pointSet) = delete;
        FinitePointSet& operator=(const FinitePointSet& a_pointSet) = delete;

        inline FinitePointSet(Box a_domain, std::array<bool, DIM> a_periodic);
        inline FinitePointSet(Box a_domain, bool a_periodic);
        inline FinitePointSet(ProblemDomain a_domain);
        inline ~FinitePointSet();
        
        inline bool add(Point a_point);
        inline void exchange();
        inline void exchange_old();
        inline void clear();
        inline std::vector<Point> points();
        inline const std::vector<Point>& localPoints() const {return m_pointBuffer;}
        private:

        ProblemDomain m_domain;
        std::vector<Point> m_pointBuffer;
        std::set<Point> m_points;
        int* m_recvBuffer;
        int* m_sendBuffer;
        int  m_recvBufferSize;
        int  m_sendBufferSize;

        int m_localBufferSize;
        int m_globalBufferSize;
        int* m_bufferSizes;
        int* m_bufferOffsets;
        int* m_localBuffer;
        int* m_globalBuffer;
    };

#include "implem/Proto_FinitePointSetImplem.H"
} // end namespace Proto

#endif //end include guard
