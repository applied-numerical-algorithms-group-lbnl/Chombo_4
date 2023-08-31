#pragma once

#ifndef _PROTO_COORD_PERMUTATION_
#define _PROTO_COORD_PERMUTATION_

#include "Proto.H"
#include <array>

namespace Proto
{
    class CoordPermutation
    {
        public:

        /// Trivial Permutation
        static CoordPermutation identity();
        /// Clockwise (-90* rotation) in the XY plane
        static CoordPermutation cw();
        /// Counter-Clockwise (+90* rotation) in the XY plane
        static CoordPermutation ccw();
        /// Invert all coordinates
        static CoordPermutation reverse();
        /// Rotate the plane containing two coordinates such that c0 -> c1
        /**
         * Note that swapping the order of the inputs reverses the direction or rotation.
        */
        static CoordPermutation rotatePlane(unsigned int a_c0, unsigned int a_c1);

        /// Default / Identity constructor
        CoordPermutation();
        /// Constructor
        /** Syntax:
         *  <code>CoordPermutation p {{c0, c1, s}, {c0, c1, s}};</code>
         *  Which is meant to be read as "map coordinate c0 to coordinate c1 with
         *  an optional reflection s".
         *  Any coordinates omitted from the constructor are understood to be mapped
         *  through a positive identity
         *
         *  Example: 2D rotation in 3D:
         *  <code>CoordPermutation rot{{0,1,1}, {1,0,-1}};</code>
         *  coord 0 -> +coord 1
         *  coord 1 -> -coord 0
         *  coord 2 -> +coord 2
        */
        CoordPermutation(std::initializer_list<std::array<int, 3>> a_args);
        CoordPermutation(std::vector<std::array<int, 3>>& a_args);
        inline void define(std::vector<std::array<int, 3>>& a_args);
        
        /// Direct Matrix Construction
        CoordPermutation(std::array<std::array<int, DIM>, DIM> a_matrix);

        /// Permute the coordinates of a Point
        Point operator()(Point a_pt) const;
        
        /// Compute and return the inverse Permutation
        CoordPermutation inverse() const; 

        Point rotateCell(Point a_point, Box a_box) const; 
        Point rotatePoint(Point a_point, Box a_box) const; 
        
        std::array<std::array<int, DIM>, DIM> matrix() const {return m_matrix; }
        
        bool operator==(const CoordPermutation& a_rhs) const;

        void print() const;

        void printMatrix() const;


        private:
        bool isValid() const; // checks if the permutation is valid

        std::array<std::array<int, DIM>, DIM> m_matrix; // discrete rotation matrix
    };

    CoordPermutation operator*(const CoordPermutation& a_A, const CoordPermutation& a_B);
#include "implem/Proto_CoordPermutationImplem.H"
} // end namespace Proto
#endif //end include guard
