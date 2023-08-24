#pragma once
#ifndef _PROTO_AMR_GRID_
#define _PROTO_AMR_GRID_

#include "Proto_FinitePointSet.H"
#include "Proto_DisjointBoxLayout.H"
#include "Proto_LevelBoxData.H"
#include "Proto_InterpStencilTP.H"

//#define PR_AMR_REFRATIO 2

namespace Proto
{
    // NB: Tags are encoded using short instead of bool so that Stencils can be used to
    // manipulate them (e.g. for buffering). -CLG
    typedef LevelBoxData<short, 1, MEMTYPE_DEFAULT, PR_CELL> LevelTagData;
    typedef BoxData<short, 1, MEMTYPE_DEFAULT> TagData;
    
    /// AMR Grid
    /**
        AMRGrid is a collection of DisjointBoxLayout objects that define
        properly nested unions of logically rectangular space.
    */
    class AMRGrid
    {
        public:

        /// Compute Tags
        /**
            Populate tags using the magnitude of the undivided difference of the a specified component of a dataset.
           
            \param tags         The tag data to be populated. Must have compatible layout with the level being checked.
            \param data         The level of data to compute tags on
            \param comp         Component of data to use
            \param bufferSize   Distance to buffer the computed tags by
            \param threshold    Threshold above whcih a tag is generated
        */
        template<typename T, unsigned int C=1, MemType MEM = MEMTYPE_DEFAULT, Centering CTR = PR_CELL>
        inline static void computeTags(
            LevelTagData&                 a_tags, 
            LevelBoxData<T, C, MEM, CTR>& a_data,
            Point                         a_bufferSize,
            T                             a_threshold,
            unsigned int                  a_comp = 0);
        
        /// Compute Tags (Isotropic)
        /**
            Computes tags with an isotropic buffer radius. 
        */
        template<typename T, unsigned int C=1, MemType MEM = MEMTYPE_DEFAULT, Centering CTR = PR_CELL>
        inline static void computeTags(
            LevelTagData&                 a_tags, 
            LevelBoxData<T, C, MEM, CTR>& a_data,
            int                           a_bufferSize,
            T                             a_threshold,
            unsigned int                  a_comp = 0);

        /// Buffer Tags
        /**
            Grows a level of tag data by bufferSize.
        */
      inline static void buffer(LevelTagData& a_tags, Point a_bufferSize);
        
        /// Buffer Tags (Isotropic)
        /**
            Grows a level of tag data by bufferSize.
        */
        inline static void buffer(LevelTagData& a_tags, int a_bufferSize);
        
        /// Trivial Constructor
        AMRGrid() {m_defined = false;}

        /// Multi-Level Constructor
        /**
            Create an AMRGrid using a vector of DisjointBoxLayouts.
            TODO: What to do if the DBLs are not properly nested?

            \param layouts      A vector of DisjointBoxLayouts. 
            \param refRatios    Refinement ratios for each *coarse* level. 
            \param maxLevels    The maximum number of levels allowable in this grid.
        */
        inline AMRGrid(
            std::vector<DisjointBoxLayout>& a_layouts,
            std::vector<Point>& a_refRatios,
            int a_maxLevels);
        
        /// Single-Level Constructor
        /**
            Create an AMRGrid by initializing only the coarsest level.

            \param layout       Layout of the coarsest level.
            \param refRatios    Refinement ratios for each *coarse* level. 
            \param maxLevels    The maximum number of levels allowable in this grid.
        */
        inline AMRGrid(
            DisjointBoxLayout& a_layout,
            std::vector<Point>& a_refRatios,
            int a_maxLevels);
        
        /// Multi-Level Define
        /**
            Lazy multi-level construction
            
            \param layouts      A vector of DisjointBoxLayouts. 
            \param refRatios    Refinement ratios for each *coarse* level. 
            \param maxLevels    The maximum number of levels allowable in this grid.
        */
        inline void
        define(
            std::vector<DisjointBoxLayout>& a_layouts,
            std::vector<Point>& a_refRatios,
            int a_maxLevels);        
        
        /// Single-Level Define
        /**
            Lazy single-level construction
            
            \param layout       Layout of the coarsest level.
            \param refRatios    Refinement ratios for each *coarse* level. 
            \param maxLevels    The maximum number of levels allowable in this grid.
        */
        inline void
        define(
            DisjointBoxLayout& a_layout,
            std::vector<Point>& a_refRatios,
            int a_maxLevels);        
        
        /// Layout Access (Const)
        inline const DisjointBoxLayout&
        operator[](unsigned int a_level) const;
        
        /// Layout Access (Non-Const)
        inline DisjointBoxLayout&
        operator[](unsigned int a_level);

        /// Regrid
        /**
            Given tag data on a level, modify or generate a grid on the next finer level.
            Input level must correspond to the grid on which input tags are defined. 
            Input level must be coarser than the finest permissible level in the grid.
            For grids with fewer than their predefined max number of levels, this
            function can add at most one level of grid.
            Proper nesting is NOT enforced by this function (see enforceNesting).
        */
        inline void regrid(LevelTagData& a_tags, unsigned int a_level, Point a_boxSize);
        inline void regrid(LevelTagData& a_tags, unsigned int a_level);
        
        /// Enforce Nesting
        /**
            Enforce the proper nesting criteria on a given level.
            Input level must be an "internal" level. That is to say, greater than 0 
            and finer than the finest level in the grid.
        */
        inline void enforceNesting(unsigned int a_level, int a_nestingDistance = 1);

        /// Trim a_level grids so that they properly nest in a_level-1 grids.  
        inline void enforceNesting2(unsigned int a_level, int a_nestingDistance = 1);
      
        /// Query Number of Levels
        inline int numLevels() const { return m_layouts.size(); }

        /// Query Max Number of Levels
        inline int maxLevels() const { return m_maxLevels; }

        /// Query Refinement Ratio
        inline Point refRatio(int a_level) const;

        /// Query Compatibility
        /**
            Determines if the layouts are compatible in a data iterator sence on 
            each level of this and another grid. If the two grids have different
            numbers of levels, this check is only performed for levels that 
            are defined for both grids.

            \param grid     Another AMRGrid
        */
        inline bool compatible(const AMRGrid& a_grid) const;
      
        /// Augments set of a_level tags by set intersecting a_level+2 grids.
        inline void addFinerGrids(LevelTagData& a_tags,int a_level);
        inline void addFinerTags(LevelTagData& a_tags, unsigned int a_level);
     
        inline void print() const;

        private:

        bool m_defined;
        std::vector<DisjointBoxLayout> m_layouts;
        std::vector<Point> m_refRatios;
        int m_nestingDistance;
        int m_maxLevels;
    };

#include "implem/Proto_AMRGridImplem.H"
} // end namespace Proto
#endif //end include guard
