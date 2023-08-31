#pragma once
#ifndef _PROTO_AMR_DATA_
#define _PROTO_AMR_DATA_

#include "Proto_AMRGrid.H"
#include "Proto_LevelBoxData.H"

namespace Proto
{
    
    /// AMR Data Hierarchy
    /**
        A nested hierarchy of data defined on an AMRGrid.
        Each level is stored as a LevelBoxData<T, C, MEM, CTR>.

        \tparam T       Type of data in array (int, double, etc.)
        \tparam C       Number of components
        \tparam MEM     Proto::MemType. HOST or DEVICE. 
        \tparam CTR     Centering of underlying LevelBoxData  
    */
    template<typename T, unsigned int C = 1, MemType MEM = MEMTYPE_DEFAULT, Centering CTR = PR_CELL>
    class AMRData
    {
        public:
        
        /// Trivial Constructor
        AMRData() {m_defined = false;}

        /// Grid Constructor
        /**
            Creates data holders for each layout level in the input grid with
            specified ghost data space.

            \param grid     An AMR Grid. Data is allocated only for existing levels
            \param ghost    Specifies the number of ghost cells to allocate
        */
        AMRData(AMRGrid& a_grid, Point a_ghost) { define (a_grid, a_ghost); }

        /// Single Level Constructor
        /**
            Creates an AMRGrid with a single level.

            \param layout       Coarsest level layout
            \param maxLevels    Maximum number of levels allowable in the grid
            \param ghost    Specifies the number of ghost cells to allocate
        */
        /*
        AMRData(DisjointBoxLayout& a_layout, int a_maxLevels, Point a_ghost)
        {
            AMRGrid grid(a_layout, a_maxLevels);
            define(grid, a_ghost);
        }
        */
        /// Define
        /**
            Lazy construction. If *this is not yet defined, data holders are created based on the
            inputs.
            TODO: If *this was previously defined, redefine it based on the new grid
            TODO: Determine how exiting data should be moved / reinitialized
        */
        inline void
        define(AMRGrid& a_grid, Point a_ghost);
      
       /// Reinitialize grids and data due to regridding.
      /**
         For all levels starting at level a_lbase + 1, initialize new data holders.
          (1) Define new DisjointBoxLayout, LevelBoxData.
          (2) Interpolate / initialize from next coarser level.
          (3) Copy on intersection from ald data at this level.
       */
        inline void
        regrid(AMRGrid& a_newgrid, int a_lbase, int a_order);
      
        /// Grid Access (Const)
        inline const AMRGrid& grid() const {return m_grid; }
        
        /// Grid Access (Non-Const)
        inline AMRGrid& grid() { return m_grid;}
        
        /// Level Data Access (Const)
        inline const LevelBoxData<T, C, MEM, CTR>&
        operator[](unsigned int a_level) const;
        
        /// Level Data Access (Non-Const)
        inline LevelBoxData<T, C, MEM, CTR>&
        operator[](unsigned int a_level);

        /// Zero Initialize
        inline void setToZero();

        /// Function Initialization
        /**
            Initializes *this using a function with a signature similar to that used by
            Proto::forallInPlace_p. The input function must have the signature:

            <code>
            void a_func(Point& a_point, Var<T, C, MEM>& a_thisData, std::array<T, DIM> a_dx, ... )
            </code>

            Unlike Proto::forall and its variants, this function cannot manage additional
            data holders (e.g. Srcs should not include a LevelBoxData or similar complex object).
            Unlike LevelBoxData::initialize(...), valid input functions must also include a 
            double argument representing the grid spacing, even if it is not used by the function
            itself. Elements of a_srcs should be plain-old-data such as single values,
            Points, or Boxes. 
        
            \param dx0      Isotropic grid spacing on the coarsest level
            \param func     Initialization function
            \param srcs     Variadic arguments after the grid spacing in func
        */
        template<typename Func, typename... Srcs>
        inline void initialize(T a_dx0, Func& a_func, Srcs... a_srcs);
        
        template<typename Func, typename... Srcs>
        inline void initialize(std::array<T, DIM>& a_dx0, Func& a_func, Srcs... a_srcs);
        
        /// Function Initialization (With Convolution)
        /**
            Initialize *this as in <code>initialize</code> and then convolve the result using
            a 4th-order convolution stencil to yield a cell-averaged field useful for finite
            volume operations. 
            TODO: Do the correct thing here if CTR != PR_CELL

            \param dx0      Isotropic grid spacing on the coarsest level
            \param func     Initialization function
            \param srcs     Variadic arguments after the grid spacing in func
        */
        template<typename Func, typename... Srcs>
        inline void initConvolve(double a_dx0, Func& a_func, Srcs... a_srcs);

        /// Average Down All Levels
        /**
            Synchronizes data between levels by recursively replacing coarse level data
            with the geometric average of overlying fine data where it exists.

            TODO: Do the correct thing here if CTR != PR_CELL
            TODO: Remove this from the AMRData API and make standalone(?)
        */
        inline void averageDown();
        
        /// Average Down Single Level
        /**
            Replace the data on a specified level with the geometric average of the overlying
            fine data should any exist. If the specified level is the finest in the 
            hierarchy, this function is a null-op. There are more than one level finer than
            the specified level, the average down is NOT computed recursively as in 
            <code>averageDown()</code>. 
            TODO: Do the correct thing here if CTR != PR_CELL
            TODO: Remove this from the AMRData API and make standalone(?)

            \param crseLevel    The level onto which finer level data will be averaged. 
        */
        //inline void averageDown(int a_crseLevel);

        /// Interpolate Fine Boundaries
        /**
            Populate coarse-fine boundary ghost cells with data interpolated from the next 
            coarser level for all levels > 0.
             
            TODO: figure out how this should be defined for CTR != 0 (if at all)
            TODO: Remove this from the AMRData API and make standalone
        
            \param a_order  Order of interpolation. Valid values: 3, 5 (default 5)
        */
        //inline void interp(int a_order = 5);

        /// Interpolate Fine Boundaries
        /**
            Populate coarse-fine boundary ghost cells with data interpolated from the next 
            coarser level for the specified fine level ( > 0 )
             
            TODO: figure out how this should be defined for CTR != 0 (if at all)
            TODO: Remove this from the AMRData API and make standalone
        
            \param a_order  Order of interpolation. Valid values: 3, 5 (default 5)
        */
        //inline void interpLevel(int a_level, int a_order = 5);

        /// Increment
        /**
            Add a scaled multiple of another AMRData to this. This function will fail
            if the two AMRData do not share an AMRGrid configuration.
            TODO: Implement += and -= operators for completeness

            \param data     Another AMRData with the same AMRGrid configuration as *this
            \param scale    Optional scaling for the right-hand side data 
        */
        inline void increment(AMRData<T, C, MEM, CTR>& a_data, T a_scale = 1.0);

        /// Compute Integral
        /**
            Compute the discrete integral over all valid data.
            TODO: Right now this is just the integral over level 0. This is the correct
            thing to do so long as averageDown has been called prior. 

            \param cdx  Coarsest grid spacing (Level 0)
        */
        inline T integrate(T a_cdx, unsigned int a_c = 0);
        inline T integrate(const std::array<T, DIM>& a_cdx, unsigned int a_c = 0);
        
        /// Compute Integral of Absolute Value
        /**
            Compute the discrete integral over all valid data.
            TODO: Right now this is just the integral over level 0. This is the correct
            thing to do so long as averageDown has been called prior. 

            \param cdx  Coarsest grid spacing (Level 0)
        */
        inline T integrateAbs(T a_cdx, unsigned int a_c = 0);
        inline T integrateAbs(const std::array<T, DIM>& a_cdx, unsigned int a_c = 0);
        
        /// Compute Max Norm
        /**
            Compute the max norm of a specified component.
            
            \param c    A component in [0, C)
        */
        inline T absMax(unsigned int a_c = 0);
       
        /// Exchange
        /**
            Exchanges ghost cell data on all levels.
        */
        inline void exchange();
        
        /// Copy To
        /**
            Copies data from this into <code>a_rhs</code> on all levels shared
            by the two AMRData objects.

            \param a_rhs    Another AMRData.
        */
        inline void copyTo(AMRData<T, C, MEM, CTR>& a_rhs);
            
        /// Get Number of Levels
        /**
            Returns the number of levels in *this with a defined LevelBoxData.
            This value will in general be different from the maximum number of levels.
            TODO:
            if numLevels() and grid().numLevels() are different, this is probably a bug. -CLG
        */
        inline int numLevels() const { return m_data.size(); }

        /// Get Max Levels
        inline int maxLevels() const { return m_grid.maxLevels(); }

        private:

        bool m_defined;
        AMRGrid m_grid;
        Point m_ghost;
        std::vector<std::shared_ptr<LevelBoxData<T, C, MEM, CTR>>> m_data;
    };
    
    typedef AMRData<short, 1, MEMTYPE_DEFAULT, PR_CELL> AMRTagData;

#include "implem/Proto_AMRDataImplem.H"
}
#endif //end include guard
