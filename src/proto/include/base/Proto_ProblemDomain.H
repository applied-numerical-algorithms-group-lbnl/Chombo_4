#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif
#ifndef _PROTO_PROBLEMDOMAIN_H_
#define _PROTO_PROBLEMDOMAIN_H_
#include <array>
#include "Proto.H"

using namespace std;
namespace Proto
{
    /// Represents a rectangular domain over which a problem can be defined, including periodic images.
    /**
        A ProblemDomain is specified by a Box defining the rectangular set of Points over which a 
        problem may be defined, plus an array<bool,DIM> specifying whether or not each of the
        coordinate directions is periodic. The Box is assumed to have the low corner to be all Zeros,
        and specifies the overall range of indices that contain all the boxes in a DisjointBoxLayout. 
    */
    class ProblemDomain
    {
        public:
        
        ProblemDomain(){};
        
        /// Anisotropic Constructor
        /**
            Constructs a ProblemDomain with Box a_bx, and periodic directions given by a_periodic.
        
            \param a_box        Domain Box
            \param a_isPeriodic Periodic flags
        */
        inline ProblemDomain(const Box& a_box,const array<bool,DIM> a_isPeriodic);
        
        /// Isotropic Constructor
        /**
            Isotropic overload of the main constructor. Creates a ProblemDomain that is periodic
            or non-periodic in every direction

            \param a_box        Domain Box
            \param a_isPeriodic Periodic flag
        */
        inline ProblemDomain(const Box& a_box, bool a_isPeriodic);
        
        /// Define
        /**
            Lazily constructs a ProblemDomain with Box a_bx, and periodic directions given 
            by a_periodic.
        
            \param a_box        Domain Box
            \param a_isPeriodic Periodic flags
        */
        inline void define(const Box& a_box,const array<bool,DIM> a_isPeriodic);

        /// Box Intersection
        /**
            Returns the Box that is the intersection of a_bx with this->box() 
            and all of its periodic images. Points in a_box which have no periodic
            image in the ProblemDomain are truncated.

            This function does no periodic shifting. See <code>image</code> for that functionality.

            \param a_box    A Box
        */
        inline Box operator&(Box a_box) const;
      
        /// Contains Point
        /**
            Checks if a_pt or any of it's periodic images are contained in this ProblemDomain.
            
            \param a_pt    A Point
        */
        inline bool contains(const Point& a_pt) const;

        /// Periodic Image
        /**
            Computes the periodic image of a Point. This function results in an error if
            the input does not have a periodic image, so one should call 
            <code>contains(Point)</code> first.

            \param a_pt    A Point with a valid periodic image
        */ 
        inline Point image(const Point& a_pt) const;
        
        /// Periodic Image
        /**
            Computes the periodic image of a Box. Points with no periodic image are truncated.
            If <code>a_box</code> contains a periodic boundary in its interior, this function
            fails by assertion.

            \param a_pt    A Point with a valid periodic image
        */ 
        inline Box image(const Box& a_box) const;
        
        /// Check If Coarsenable
        /**
            Returns true if *this can be tiled by boxes of size <code>a_boxSize</code>.
            This is a requirement for using *this to construct a DisjointBoxLayout with
            Boxes of this size.

            \param a_boxSize    A possibly anisotropic Box size
        */
        inline bool coarsenable(Point a_boxSize) const;

        /// Query Size
        /**
            Size in number of points in each direction.
        */
        inline Point sizes() const { return m_box.sizes(); }
        
        /// Equality
        /**
            Two ProblemDomain objects are "==" if they have the same
            Box domain and periodicity.

            \param a_domain  A ProblemDomain
        */
        inline bool operator==(const ProblemDomain& a_domain) const;
        
        /// Inequality
        /**
            Two ProblemDomain objects are "==" if they have the same
            Box domain and periodicity.

            \param a_domain  A ProblemDomain
        */
        inline bool operator!=(const ProblemDomain& a_domain) const;

        /// Get Box
        inline Box box() const {return m_box;};
        
        /// Get Periodicity
        inline array<bool,DIM> periodicity() const { return m_isPeriodic; }
      
        /// Check Periodic Direction
        /**
            Check if this domain is periodic in direction a_dir.

            \param a_dir    A coordinate direction in [0,DIM)
        */
        inline bool isPeriodic(int a_dir) {return m_isPeriodic[a_dir]; }

        /// Coarsen (Anisotropic)
        /**
            Returns a ProblemDomain with a box given by this->box.coarsen(a_refRatio).
            Fails if this->coarsenable(a_refRatio) is false.

            \param a_refRatio   A vector of refinement ratios
         */
        inline ProblemDomain coarsen(Point a_refRatio) const;
        
        /// Coarsen (Isotropic)
        /**
            Returns a ProblemDomain with a box given by this->box.coarsen(a_refRatio).
            Fails if this->coarsenable(a_refRatio) is false
            
            \param a_refRatio   An isotropic refinement ratio
         */
        inline ProblemDomain coarsen(int a_refRatio) const;

        /// Refine (Anisotropic)
        /**
            Returns a ProblemDomain with a box given by this->box.refine(a_refRatio).
            
            \param a_refRatio   A vector of refinement ratios
        */
        inline ProblemDomain refine(Point a_refRatio) const;
        
        /// Refine (Isotropic)
        /**
            Returns a ProblemDomain with a box given by this->box.refine(a_refRatio).
            
            \param a_refRatio   An isotropic refinement ratio
        */
        inline ProblemDomain refine(int a_refRatio) const;
        
        protected:
        
        Box             m_box;
        array<bool,DIM> m_isPeriodic;
        bool            m_isDefined;
        
    }; // end class ProblemDomain

    /// Stream output for ProblemDomain.
    inline std::ostream& operator<< (std::ostream& os, const ProblemDomain& a_pd);
#include "implem/Proto_ProblemDomainImplem.H"
} // end namespace Proto  
#endif //end include guard
