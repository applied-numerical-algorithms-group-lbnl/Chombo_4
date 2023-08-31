#ifndef _PROTO_INTERP_STENCIL_H_
#define _PROTO_INTERP_STENCIL_H_
#include "Proto_Lapack.H"
#include "Proto_Stencil.H"

namespace Proto {

/// Interpolation Stencil
/**
    \ingroup stencil_operations
    InterpStencil is a convenient abstraction which encapsulates the functionality of an operation in which
    the destination refinement is non-trivial. This kind of operation can be done with normal Stencil objects,
    however for a destination refinement of R, an array of R^DIM Stencils will be needed in general.
*/
template <class T>
class InterpStencil
{
public:
    /// Default Constructor
    inline InterpStencil() { m_closed = false; }

    /// Isotropic Constructor
    /**
        \ingroup stencil_operations
        Builds an InterpStencil assuming a destination refinement of <code> a_ratio </code> in all directions.

        \param a_ratio  destination refinement
    */
    inline InterpStencil(int a_ratio) { define(Point::Ones(a_ratio)); }

    /// Anisotropic Constructor
    /**
        \ingroup stencil_operations
        Builds an InterpStencil with a destination refinement of a_ratio;

        \param a_ratio  destination refinement
    */
    inline InterpStencil(Point a_ratio) { define(a_ratio); }

    /// Define
    /**
        \ingroup stencil_operations
        Builds an InterpStencil with a destination refinement a_ratio.

        \param a_ratio  destination refinement
    */
    inline void define(Point a_ratio);
    
    /// Close
    /**
        @private
        Computes all of the destRatio and shift metadata for each Stencil in the
        InterpStencil array. This function is called automatically when the InterpStencil
        is used the first time, so there is no need for a user to mess with it.
    */
    inline void close();

    /// Get Read-Write Stencil Subcomponent
    /**
        \ingroup stencil_operations
        Returns a reference to the Stencil associated with the destShift a_p.

        \param a_p  A destination shift in [(0,...,0), (destRatio-1,...,destRatio-1)]
    */
    inline Stencil<T>& operator()(Point a_p);

    /// Get Read-Only Stencil Subcomponent
    /**
        \ingroup stencil_operations
        Returns a const reference to the Stencil associated with the destShift a_p.

        \param a_p  A destination shift in [(0,...,0), (destRatio-1,...,destRatio-1)]
    */
    inline const Stencil<T>& get(Point a_p) const;

    /// Span
    /** \ingroup stencil_operations
        Computes the composite span of all Stencils in the InterpStencil
    */
    inline Box span() const;

    /// Ghost Size
    /**
        \ingroup stencil_operations
        Conservatively computes the amount of ghost data needed to apply this InterpStencil.
        May overestimate in the case of assymmetric InterpStencil footprints.
    */
    inline Point ghost() const;

    /// Lazy Apply
    /**
        \ingroup stencil_operations
        Gathers inputs for the composite Stencil operation and stores them in an intermediate
        LazyStencil to be consumed by a left hand side operator ("=", "|=", or "+=").
        The syntax is identical to that used by Stencil.

        \param a_src    Source data
        \param a_scale  (Optional) value to scale the output of *this by.
    */
    template<unsigned int C, MemType MEMTYPE, unsigned char D, unsigned char E>
    inline LazyStencil<T,C,MEMTYPE,D,E> operator()(
            const BoxData<T,C,MEMTYPE,D,E>& a_src,
            T                               a_scale = 1);
    
    /// Lazy Apply
    /**
        \ingroup stencil_operations
        Gathers inputs for the composite Stencil operation and stores them in an intermediate
        LazyStencil to be consumed by a left hand side operator ("=", "|=", or "+=").
        The syntax is identical to that used by Stencil.
        Explicitly defining a Box input is not recommended unless absolutely necessary.

        \param a_src    Source data
        \param a_box    Iteration box
        \param a_scale  (Optional) value to scale the output of *this by.
    */
    template<unsigned int C, MemType MEMTYPE, unsigned char D, unsigned char E>
    inline LazyStencil<T,C,MEMTYPE,D,E> operator()(
            const BoxData<T,C,MEMTYPE,D,E>& a_src,
            Box                             a_box,
            T                               a_scale = 1);
    
//=======================================================================================
// InterpStencil Library ||
//=======================++

    /// Piecewise Constant Interpolation
    /**
        \ingroup stencil_library
        \param a_ratio  Anisotropic destination refinement
    */
    static inline InterpStencil<T> PiecewiseConstant(Point a_ratio);
    /// Piecewise Linear Interpolation
    /**
        \ingroup stencil_library
        \param a_ratio  Anisotropic destination refinement
    */
    static inline InterpStencil<T> PiecewiseLinear(Point a_ratio);
    /// Piecewise Linear Interpolation
    /**
        Computed based on tensor products of 1 dimensional interpolations.
        \ingroup stencil_library
        \param a_ratio  Isotropic destination refinement
    */
    static inline InterpStencil<T> Quadratic(int a_ratio);

    /// Build High-Order Interpolation Stencil (Simple)
    /**
        Builds a symmetric, <code>a_order</code> order accurate interpolation operator
        with an isotropic refinement ratio.
        
        \param a_order      Order of accuracy
        \param a_refRatio   Isotropic refinement ratio
    */
    static inline InterpStencil<T> Build(int a_order, int a_refRatio);
    
    /// Build High-Order Interpolation Stencil (Simple, Anisotropic)
    /**
        Builds a symmetric, <code>a_order</code> order accurate interpolation operator
        with an isotropic refinement ratio.
        
        \param a_order      Order of accuracy
        \param a_refRatio   Isotropic refinement ratio
    */
    static inline InterpStencil<T> Build(int a_order, Point a_refRatio);

    /// Build General High-Order Interpolation Stencil (General)
    /**
        Builds an InterpStencil with a custom footprint using a Least-Squares solve.

        \ingroup stencil_library
        \param a_shiftMax       Maximum Manhattan distance from the Stencil center
        \param a_shiftKernel    Maximum span of the InterpStencil
        \param a_order          Order of polynomial fit (e.g. 2 for quadratic)
        \param a_refRatio       Desired ratio of interpolation. Must be a power of 2
    */
    static inline InterpStencil<T> Build(
            int a_shiftMax,
            Box a_shiftKernel,
            int a_order,
            int a_refRatio);

    /// Build General High-Order Interpolation Stencil (General, Anisotropic)
    /**
        Builds an InterpStencil with a custom footprint using a Least-Squares solve.
        
        \ingroup stencil_library
        \param a_shiftMax       Maximum Manhattan distance from the Stencil center
        \param a_shiftKernel    Maximum span of the InterpStencil
        \param a_order          Order of polynomial fit (e.g. 2 for quadratic)
        \param a_refRatio       Desired ratio of interpolation. Must be a power of 2
    */
    static inline InterpStencil<T> Build(
            int a_shiftMax,
            Box a_shiftKernel,
            int a_order,
            Point a_refRatio);

    /// Build General High-Order Interpolation Stencil
    /**
        Builds an InterpStencil with a custom footprint using a Least-Squares solve.
        
        \ingroup stencil_library
        \param a_shifts         Desired shifts of the stencil.
        \param a_maxOrder       Highest desired order of Polynomial interpolant.
        \param a_refRatio       Desired ratio of interpolation. Must be a power of 2.
    */
    static inline InterpStencil<T> Build(
            std::vector<Point>& a_shifts,
            int                 a_maxOrder,
            Point               a_refRatio);
    
    /// Build General High-Order Interpolation Stencil
    /**
        Builds an InterpStencil with a custom footprint using a Least-Squares solve.
        Orders are specified as Point such that (1,2,3) is equivalent to x*y^2*z^3.

        \ingroup stencil_library
        \param a_shifts         Desired shifts of the stencil.
        \param a_orders         Polynomial exponents in vector form.
        \param a_refRatio       Desired ratio of interpolation. Must be a power of 2
    */
    static inline InterpStencil<T> Build(
            std::vector<Point>&       a_shifts,
            const std::vector<Point>& a_orders,
            Point                     a_refRatio);

    /// Return Shift Kernel
    /**
        \ingroup stencil_operations
    */
    inline Box kernel() const {return m_k;}
    /// Return Destination Refinement Ratio
    /**
        \ingroup stencil_operations
    */
    inline Point ratio() const {return m_r;}
    /// Iterate Over Stencil Components
    /**
        \ingroup stencil_operations
    */
    inline typename std::vector<Stencil<T>>::iterator begin(){return m_s.begin();}
    /// Iterate Over Stencil Components
    /**
        \ingroup stencil_operations
    */
    inline typename std::vector<Stencil<T>>::iterator end(){return m_s.end();}

    /// Empty Query
    /**
        \ingroup stencil_operations
    */
    inline bool empty(){return (m_s.size() <= 0);}

    /// Size Query
    /**
        Returns the number of Stencil components in *this.
        \ingroup stencil_operations
    */
    inline unsigned long long int size() const {return m_s.size();}

    inline void print() const;
private:
    Point m_r;                      ///< Destination refinement ratio
    Box m_k;                        ///< Kernel of destination shifts
    std::vector<Stencil<T>> m_s;    ///< Array of Stencils
    bool m_closed = false;          ///< Has this been closed properly
};

#include "implem/Proto_InterpStencilImplem.H"
#include "implem/Proto_InterpStencilDefs.H"
} // end namespace Proto
#endif // End Include Guard
